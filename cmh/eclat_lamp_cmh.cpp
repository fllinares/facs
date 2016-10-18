// AUTHOR: FELIPE LLINARES


// -----------------------------INCLUDES-----------------------------------------


#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <algorithm>
#include <numeric>
#include <vector>
#include <unordered_set>
#include <set>
#include <unistd.h>

#include <math.h>

#include "shared_vars.h"
#include "tarone_cmh.cpp"
#include "time_keeping.cpp"

using namespace std;

// ------------------------DATA STRUCTURES--------------------------------------


// Structure holding a vertical representation of a (conditional) transaction database
struct Database{
	// List of lists holding the transactions in which each item appears
	vector< vector<int> > transactions;
	// List of boolean vectors, such that transactions_bool[i][j]==true iff j is contained in transactions[i]. This
	// redundancy will in fact help make the intersection of transaction lists more efficiently
	vector< vector<bool> > transactions_bool;
	// List holding the corresponding item IDs
	vector<int> items;
	// List holding the corresponding supports
	vector< vector<int> > supports;
	// List holding the minimum attainable p-values
	vector<double> minpvals;
	// List holding the corresponding values of the minimum attainable p-value lower envelope
	vector<double> lower_envelope_minpvals;
};


// -----------------------------GLOBAL VARIABLES----------------------------------


int max_itemset_size; // Default: no limit in itemset size


// Number of distinct items in the transaction database
int n_items;
// Mapping of internal item labeling (consecutive ints) to original item label
vector<int> item_label_map;
// Total number of samples
int n_samples;
// Number of samples in minority class
int n1;


// Number of enumerated itemsets
long long n_enumerated;
long long n_enumerated_closed;

// Array of booleans containing class labels
vector<bool> labels;

// Output stream to write results summary into. Set to stdout by default,
// unless optional arguments override that behaviour
ostream *summary_file;

// ----------------------------------DEBUG FUNCTIONS ---------------------------------


void print_transaction_database(Database &db){
	for(int i=0; i < db.items.size(); ++i){
		cout << i << "/" << item_label_map[i] << " (";
		for(int j=0; j < n_cat-1; ++j) cout << db.supports[i][j] << ",";
		cout << db.supports[i][n_cat-1] << "): ";
		for(int j=0; j < db.transactions[i].size()-1; ++j) cout << db.transactions[i][j] << ",";
		cout << db.transactions[i][db.transactions[i].size()-1] << endl;
	}
}


void print_vector(vector<int> v){
	for(int i=0; i < v.size()-1; ++i) cout << v[i] << ",";
	cout << v[v.size()-1] << endl;
}


// ----------------------------------I/O-FUNCTIONS--------------------------------


// Read transactions file, returning a Database struct containing a vertical representation
// of the transaction database.
// INPUTS:
// 		filename: Path to text file containing transactions database
//		db: Reference to Database object in which the final transaction database will be stored.
//			This will be a vertical representation, with perfect extensions of the root removed
//		pexs: Reference to vector of integers where the internal item labels of all items which are perfect
//			  extensions of the root will be stored
int read_transactions(const string filename, Database &db, vector<int> &pexs){
	// List of transactions
	vector< vector<int> > T_list;
	// Temporal transaction database (in vertical representation)
	Database tmp_db;


    // String to store each line read from the file
    string line;
    // And a Stringstream to parse it
    stringstream ss_line;
    // Integer to store each index read from each field in a line of the file
    int idx;
    // Number of transactions read
    int n = 0;
    // Set of all items in original transaction database
    set<int> items;

    // Open file
    ifstream f(filename);
    // Check if file was correctly opened
    if(f.is_open()){
    	// Process file line by line
    	while(getline(f,line)){
    		// Vector to store current transaction
    		vector<int> T;

    		// Move line to stringstream for parsing
			ss_line << line;
			// Parse each field (assume that they are separated by some sort
			// of whitespaces or tabs), appending them to current transaction
			while(ss_line >> idx) T.push_back(idx);
			// Increase counter of number of transactions read and reset
			// stringstream flags for the next iteration
			n++; ss_line.clear();

			T_list.push_back(T);
			items.insert(T.begin(), T.end());
    	}
    	// Close file
    	f.close();

    	// Total number of items in transaction database
    	n_items = items.size();
		tmp_db.transactions.resize(n_items);
    	// And relabel them as consecutive integers, keeping track of the mapping back to
    	// their original labels.
		for(auto item : items) item_label_map.push_back(item);

    	// Build temporary transaction database
    	int n_trans = 0, tmp_idx;
		vector<int> sum_of_transaction_sizes(n_items); //sum_of_transaction_sizes[i]=sum of sizes of all transactions containing item i
    	for(auto trans : T_list){
    		// Add transaction index to transaction lists of all items contained in transaction. For that, we need to map
    		// the original label index to the consecutive indexing using internally, using find() and distance() for that
    		// purpose
			for(auto item : trans){
				tmp_idx = distance(item_label_map.begin(), find(item_label_map.begin(), item_label_map.end(), item));
				tmp_db.transactions[tmp_idx].push_back(n_trans);
				sum_of_transaction_sizes[tmp_idx] += trans.size();
			}
    		// Increment transaction index
    		n_trans++;
    	}

		// Sort items in decreasing order of cumulative size of all transactions they belong to (rationale, items
		// contained in very large transactions will lead to denser conditional transaction databases, and should be
		// processed in the "shallow" part of the enumeration tree for extra performance)
		vector<int> idx_sort(n_items);
		iota(idx_sort.begin(), idx_sort.end(), 0);
		sort(idx_sort.begin(), idx_sort.end(), [&](size_t i, size_t j){ return sum_of_transaction_sizes[i] > sum_of_transaction_sizes[j]; });

		// Create final transaction database, ordered according to idx_sort above and storing perfect extensions
		// of the root separately
		int current_support; // temporal variable to store support of each item
    	for(int i=0; i < n_items; ++i) {
    		current_support = tmp_db.transactions[idx_sort[i]].size();
    		if(current_support == n_trans) pexs.push_back(idx_sort[i]);
    		else{
    			db.transactions.push_back(tmp_db.transactions[idx_sort[i]]);
				db.transactions_bool.push_back(vector<bool>(n_trans));
				for(auto trans_idx: db.transactions.back()) db.transactions_bool.back()[trans_idx] = true;
    			db.items.push_back(idx_sort[i]);
    			//Compute supports for each category
    			db.supports.push_back(vector<int>(n_cat));
    			for(auto trans_idx : db.transactions.back()) db.supports.back()[cats[trans_idx]]++;
    			db.minpvals.push_back(compute_minpval(db.supports.back()));
    			db.lower_envelope_minpvals.push_back(compute_lower_envelope_minpval(db.supports.back()));
    		}
    	}

    	// Check that number of lines read, and total number of transactions as measured by sum of transaction
    	// cardinalities match
    	if (n != n_trans){
    		cerr << "Error @ read_transactions: Number of transactions read from file: " << n << " doesn't match sum of transaction cardinalities: " << n_trans << ".!!!" << endl;
    		exit(-1);
    	}

    	// Return transaction database
        return n_trans;

    }
    else{
        cerr << "Error @ read_transactions: Unable to open transactions file " << filename << endl;
        exit(-1);
    }
}


// Reads the file containing the labels, whose format is assumed to be either a long
// binary column, with one row per sample, or one long binary row, containing one
// column per sample.
// This function will read the file, compute the number of positive examples (stored in the
// global variable n1) and store all the actual labels in the global variable labels (which
// is also allocated within this function).
int read_labels_file(const string filename){
    // String to store each line read from the file
    string line;
    // And a Stringstream to parse it
    stringstream ss_line;
    // Integer to store each index read from each field in a line of the file
    bool idx;
    // Number of elements processed = number of labels read
    int n = 0;
    // Initialize n1 to 0
    n1 = 0;

    // Open file
    ifstream f(filename);
    // Check if file was correctly opened
    if(f.is_open()){
        while(getline(f,line)){
            // Move line to stringstream for parsing
            ss_line << line;
            // Read next phenotype
            while(ss_line >> idx) {
            	labels.push_back(idx);
            	if(labels[n++]) n1++;
            }
            ss_line.clear();
        }
        // Close file
        f.close();
        // Exit call, returning number of labels read
        return n;
    }
    else{
        cerr << "Error @ read_labels_file: Unable to open labels file " << filename << endl;
        exit(-1);
    }
}


// Reads the file containing the covariate categories, whose format is assumed to be either a long
// column, with one row per sample, or one long binary row, containing one column per sample.
// The categories must be encoded as consecutive integers 0, 1, 2, ..., n_cat-1
// The function will also set the global variable n_cat, containing the number of distinct categories.
int read_covariates_file(const string filename){
	// Set object to find number of categories
	unordered_set<int> unique_cats;
    // String to store each line read from the file
    string line;
    // And a Stringstream to parse it
    stringstream ss_line;
    // Integer to store each index read from each field in a line of the file
    int idx;
    // Number of elements processed = number of covariate values read
    int n = 0;

    // Open file
    ifstream f(filename);
    // Check if file was correctly opened
    if(f.is_open()){
        while(getline(f,line)){
            // Move line to stringstream for parsing
            ss_line << line;
            // Read next phenotype
            while(ss_line >> idx) {
            	cats.push_back(idx);
            	unique_cats.insert(idx);
            	n++;
            }
            ss_line.clear();
        }
        // Close file
        f.close();
        // Set number of unique categories
        n_cat = unique_cats.size();
        // Exit call, returning number of covariate values read
        return n;
    }
    else{
        cerr << "Error @ read_covariates_file: Unable to open covariates file " << filename << endl;
        exit(-1);
    }
}


// Output itemset
void output_itemset(double pval, int a, vector<int> &x_t, vector<int> &iset, vector<int> &pexs){
	sig_itemsets_file << pval << '\t' << a << '\t';
	for(auto x : x_t) sig_itemsets_file << x << '\t';
	for(int i=0; i < iset.size()-1; ++i) sig_itemsets_file << item_label_map[iset[i]] << '\t';
	if(pexs.size() > 0){
		sig_itemsets_file << item_label_map[iset[iset.size()-1]] << '\t';
		for(int i=0; i < pexs.size()-1; ++i) sig_itemsets_file << item_label_map[pexs[i]] << '\t';
		sig_itemsets_file << item_label_map[pexs[pexs.size()-1]] << endl;
	}
	else sig_itemsets_file << item_label_map[iset[iset.size()-1]] << endl;
}


// Output itemset
void output_testable_itemset(double pval){

}


// Return results of the execution
void output_results_summary(){
    // Compute corrected significance threshold
    delta_opt = target_fwer/n_testable;
    *summary_file << "DATASET CHARACTERISTICS:" << endl;
    *summary_file << "\t" << "n_samples = " << n_samples << ", n1 = " << n1 << endl;
    for(int c=0; c<n_cat; ++c) *summary_file << "\t\tCat " << c << ": " << "n_samples = " << n_samples_t[c] << ", n1 = " << n1_t[c] << endl;
    *summary_file << "RESULTS:" << endl;
    if(max_itemset_size==0) *summary_file << "Maximum itemset size to be processed: unlimited" << endl;
    else *summary_file << "Maximum itemset size to be processed: " << max_itemset_size << endl;
    *summary_file << "Associated testability threshold: " << pth << endl;
    *summary_file << "Number of itemsets enumerated: " << n_enumerated << endl;
	*summary_file << "Number of closed itemsets enumerated: " << n_enumerated_closed << endl;
    *summary_file << "Number of testable itemsets: " << n_testable << endl;
    *summary_file << "Corrected significance threshold at target FWER " << target_fwer << ": " << delta_opt << endl;
}

// -----------------------RECURSIVE ITEMSET ENUMERATION-------------------------------


int depth(Database &db, vector<int> &iset, vector<int> &pexs, vector< vector<bool> *> &elim){
	// Declare variables for loops
	int i, k, kk;
	// Declare variables to hold supports
	int current_support, max_support = 0;
	// Variable to check if itemsets are closed
	bool closed;
	// Variable to hold the value of the minimum attainable p-value and its corresponding lower envelope
	double minpval, lower_envelope_minpval;
	// Number of items in the current (conditional) transaction database
	int db_size = db.items.size();
	// Number of perfect extensions found during the processing of a new itemset
	int n_proj_pexs;
	// Output of recursion
	int r_out;
	// Variables to compute itemset pvalue (only used if compute_pals==true)
	int a;
	double pval;

	// Insert all transaction lists which must be used during closure check
	for(kk=(db_size-1); kk > 0; --kk) elim.push_back(&(db.transactions_bool[kk]));

	// Start processing all possible child nodes
	for(k=0; k < db_size; ++k){
		current_support = db.transactions[k].size();
		// Between the time the item {db.items[k]} was added to the current conditional database during the processing
		// of the parent itemset, and now, the prunability threshold might have changed, due to the earlier depth-first
		// processing of its siblings. Therefore, it's worth checking again if the itemset can be pruned.
		if(not IS_TESTABLE(db.minpvals[k])){
			// If the itemset is not testable, it might be prunable. In this case, we must compute the minpval lower envelope
			// if it was not yet computed
			if(db.lower_envelope_minpvals[k] == -2) db.lower_envelope_minpvals[k] = compute_lower_envelope_minpval(db.supports[k]);
			// Based on the value of the minpval lower envelope, check if itemset is prunable
			if(not IS_NOT_PRUNABLE(db.lower_envelope_minpvals[k])){
				if(k < (db_size-1)) elim.pop_back();
				continue;
			}
		}
		n_enumerated++;  // Increase count of enumerated nodes
		max_support = (current_support > max_support) ? current_support : max_support; //max_support=max(current_support,max_support)

		// Check closure
		closed = true;
		for(kk=(elim.size()-1); kk >=0; --kk) {
			closed = false;
			for(auto trans_idx : db.transactions[k]) {
				if (not (*elim[kk])[trans_idx]) {
					closed = true;
					break;
				}
			}
			if (not closed) break;
		}
		if(not closed) {
			if(k < (db_size-1)) elim.pop_back();
			continue;
		}

		// Build new conditional transaction database
		Database proj_db; n_proj_pexs = 0;
		for(kk=0; kk < k; ++kk){
			proj_db.transactions.push_back(vector<int>()); vector<int> &U = proj_db.transactions.back();
			proj_db.transactions_bool.push_back(vector<bool>(db.transactions_bool[k])); vector<bool> &U_hash = proj_db.transactions_bool.back();
			// Compute intersection
			for(auto trans_idx : db.transactions[k]){
				if (db.transactions_bool[kk][trans_idx]) U.push_back(trans_idx);
				else U_hash[trans_idx] = false;
			}
			// Check if item is a perfect extension
			if(U.size() == current_support){
				pexs.push_back(db.items[kk]);
				n_proj_pexs++;
				proj_db.transactions.pop_back();
				proj_db.transactions_bool.pop_back();
			}
			// If not, we might want to add the item to the next conditional transaction database
			else{
				proj_db.supports.push_back(vector<int>(n_cat)); vector<int> &proj_supports = proj_db.supports.back();
				//Compute supports for each category
				for(auto trans_idx : U) proj_supports[cats[trans_idx]]++;
				minpval = compute_minpval(proj_supports);
				// If the itemset is testable, it cannot be prunable
				if(IS_TESTABLE(minpval)){
					proj_db.items.push_back(db.items[kk]);
					proj_db.minpvals.push_back(minpval);
					// Avoid computing the minpval lower envelope. Set its value to -2 (an arbitrary value outside its normal range),
					// to indicate that the computation was skipped
					proj_db.lower_envelope_minpvals.push_back(-2);
				}
				// If it is not testable, it could be prunable or not. We need to compute the minpval lower envelope to answer that
				else{
					lower_envelope_minpval = compute_lower_envelope_minpval(proj_supports);
					if(IS_NOT_PRUNABLE(lower_envelope_minpval)){
						proj_db.items.push_back(db.items[kk]);
						proj_db.minpvals.push_back(minpval);
						proj_db.lower_envelope_minpvals.push_back(lower_envelope_minpval);
					}
					else{
						proj_db.transactions.pop_back();
						proj_db.transactions_bool.pop_back();
						proj_db.supports.pop_back();
					}
				}
			}
		}
		// Invoke next step of depth-first search recursively
		iset.push_back(db.items[k]);
		r_out = (proj_db.items.size() > 0) ? depth(proj_db, iset, pexs, elim) : 0;
		// If the itemset has no perfect extension and is testable, process it
		if((r_out < current_support) and IS_TESTABLE(db.minpvals[k])) {
			// Increase counter of testable itemsets, and upgrade histogram of supports
			freq_cnt[bucket_idx(db.minpvals[k])]++; n_testable++;
			// Check FWER condition, adjusting threshold if violated until the condition is restored
			while((n_testable*pth) > target_fwer) decrease_threshold();

			// If we have to compute pvalues for testable, and the pattern is still testable, we assume
			// that it will remain testable. This won't be true for some patterns, leading to a few unnecessary
			// pvalues being computed. However, the extra runtime will in practice be still smaller than that
			// of executing the entire enumeration process a second time, once the final testability threshold
			// is known
			if(compute_pvals and IS_TESTABLE(db.minpvals[k])){
				// Compute cell-count
				a = 0;
				for(auto trans_idx : db.transactions[k]) if(labels[trans_idx]) a++;
				pval = compute_pval(a, db.supports[k]);
				//TODO: This step is greedy. In pathological cases, some significant itemsets could be lost...
				if(pval <= tentative_sig_ths[tentative_sig_ths.size()-1]) output_itemset(pval, a, db.supports[k], iset, pexs);
				if(output_testable) testable_itemsets_file << pval << endl;
			}
		}
		if(r_out < current_support) n_enumerated_closed++; // Increase count of enumerated closed itemsets

		// Undo changes to variables shared across recursion iterations (iset, pexs and elim)
		for(i=0; i<n_proj_pexs; ++i) pexs.pop_back();
		iset.pop_back();
		if(k < (db_size-1)) elim.pop_back();
	}
	return max_support;
}


// ----------------------------------ENTRY POINT --------------------------------------


void show_usage_instructions(char *program_name){
    cout << "Usage:" << "\t" << program_name << " [-option] [option_argument]" << endl;
    cout << "Mandatory options:" << endl;
    cout << "\t" << "option:" << "\t" << "-t path_to_transactions_file" << endl; 
    cout << "\t" << "option:" << "\t" << "-l path_to_labels_file" << endl; 
    cout << "\t" << "option:" << "\t" << "-c path_to_covariates_file" << endl; 
    cout << "\t" << "option:" << "\t" << "-f target_FWER" << endl; 
    cout << "Optional options:" << endl;
    cout << "\t" << "option:" << "\t" << "-o path_to_results_summary_file" << endl; 
    cout << "\t\t" << "Default: results summary displayed in stdout" << endl;
    cout << "\t" << "option:" << "\t" << "-p path_to_profiling_summary_file" << endl; 
    cout << "\t\t" << "Default: profiling results displayed in stdout" << endl;
    cout << "\t" << "option:" << "\t" << "-s path_to_significant_itemsets_file" << endl; 
    cout << "\t\t" << "Default: significance of itemsets is not evaluated at all, only testability" << endl;
	cout << "\t" << "option:" << "\t" << "-d path_to_testable_itemsets_file" << endl;
	cout << "\t\t" << "Note: requires usage of option -s as well!" << endl;
	cout << "\t\t" << "Default: testable itemsets not shown in output" << endl;
}


int main(int argc, char *argv[]){
	// Declare some variables which cannot be initialized inside the switch statement
	char opt;

	ofstream summary_file_obj;
	ofstream profiling_file_obj;

	// Some variables to hold input arguments
	string covariates_filename;
	string transactions_filename;
	string labels_filename;
	string sig_itemsets_filename;
	double fwer;

	// Main variables used during depth-first traversal of itemset lattice
	Database db;
	vector<int> iset;
	vector<int> pexs;
	vector< vector<bool> *> elim;
	int r_out;

	// Record starting timestamp
	t_init = measureTime();

	// Process arguments (boolean flags to check that mandatory arguments are inputed)
	max_itemset_size = 0; // Default: no maximum itemset size
	summary_file = &cout; // Set stream to output results to stdout by default
	compute_pvals = false; // Default: do not compute p-values for testable itemsets
	output_testable = false; // Default: do not output testable itemsets
	bool flag_transactions_file = false, flag_labels_file = false, flag_covariates_file = false, flag_fwer = false;
	while( (opt = getopt(argc, argv, "t:l:c:f:M:o:p:s:d:")) != -1){
		switch(opt){
			// Option t: Path to transactions file
			case 't':
				transactions_filename = string(optarg); flag_transactions_file = true;
				break;
			// Option l: Path to labels file
			case 'l':
				labels_filename = string(optarg); flag_labels_file = true;
				break;
			// Option c: Path to covariates file
			case 'c':
				covariates_filename = string(optarg); flag_covariates_file = true;
				break;
			// Option f: Target FWER
			case 'f':
				fwer = atof(optarg); flag_fwer = true;
				break;
			// Option M: Maximum itemset size
			case 'M':
				max_itemset_size = atoi(optarg);
				break;
			// Option o: Write results summary to file instead than stdout
			case 'o':
				summary_file_obj.open(optarg);
				if(summary_file_obj.is_open()) summary_file = &summary_file_obj;
				else{
					cerr << "Error @ main: Unable to create output file " << optarg << endl;
					exit(-1);
				}
				break;
			// Option p: Write code profiling measurements to file
			case 'p':
				profiling_file_obj.open(optarg);
				if(profiling_file_obj.is_open()) timing_file = &profiling_file_obj;
			    else{
			        cerr << "Error @ main: Unable to create output profiling file " << optarg << endl;
			        exit(-1);
			    }
				break;
			// Option s: Evaluate significance of testable itemsets
			case 's':
				sig_itemsets_filename = string(optarg);
				sig_itemsets_file.open(sig_itemsets_filename + string(".tmp"));
				if(sig_itemsets_file.is_open()) compute_pvals = true;
				else{
					cerr << "Error @ main: Unable to create temporal significant itemsets output file " << optarg << ".tmp" << endl;
					exit(-1);
				}
				break;
				// Option d: Output all testable itemsets
			case 'd':
				testable_itemsets_file.open(optarg);
				if(testable_itemsets_file.is_open()) output_testable = true;
				else{
					cerr << "Error @ main: Unable to create testable itemsets output file " << optarg << ".tmp" << endl;
					exit(-1);
				}
				break;
			// Not supported option
			default:
				cerr << "Error @ main: Incorrect arguments!!!" << endl;
				show_usage_instructions(argv[0]);
				exit(-1);
		}
	}
	if( (not flag_transactions_file) or (not flag_labels_file) or (not flag_covariates_file) or (not flag_fwer)){
		cerr << "Error @ main: Incorrect arguments!!!" << endl;
		show_usage_instructions(argv[0]);
		exit(-1);
	}

	// Read values for covariate for each transaction
	tic = measureTime();
	int n_cov_vals = read_covariates_file(covariates_filename);
	// Read labels for each transaction
	int n_labels = read_labels_file(labels_filename);
	toc = measureTime();
	time_IO += (toc-tic);

	// Check if number of covariate values and labels read match
	if(n_labels != n_cov_vals){
		cerr << "Error @ main: Number of labels read: " << n_labels << " does not match number of covariate values read: " << n_cov_vals << "!!!" << endl;
		show_usage_instructions(argv[0]);
		exit(-1);
	}

	// Initialize Tarone-related code (must be done before reading transactions, as read_transactions() already
	// computes the minpvals and lower envelope minpvals of singleton itemsets
	tic = measureTime();
	init_tarone(fwer, n_cov_vals); // Here we hope that n_cov_vals==n_labels==n_trans, we'll check if that holds later and assign the unique value to global var n_samples
	toc = measureTime();
	time_initialisation += (toc-tic);

	// Read transactions
	tic = measureTime();
	int n_trans = read_transactions(transactions_filename, db, pexs);

	// Make sure that there is a label for each transaction
	if(n_trans == n_labels) n_samples = n_trans;
	else{
		cerr << "Error @ main: Number of transactions read: " << n_trans << " does not match number of labels read: " << n_labels << "!!!" << endl;
		show_usage_instructions(argv[0]);
		exit(-1);
	}
	toc = measureTime();
	time_IO += (toc-tic);

	// Start depth-first travesal of itemset lattice
	tic = measureTime();
	n_enumerated = 0;  // Initialize number of enumerated itemsets to zero
	n_enumerated_closed = 0;  // Same for number of enumerated closed itemsets
	r_out = depth(db, iset, pexs, elim);
	// Itemsets which are perfect extensions of the root cannot be testable, hence they are ignored
	toc = measureTime();
	time_enumeration += (toc-tic);

	// Output results summary
	tic = measureTime();
	output_results_summary();

	// Process significant itemsets, if the option was chosen
	if(compute_pvals){
		 long long n_sig = process_significant_itemsets(sig_itemsets_filename);
		 *summary_file << "Number of significant itemsets: " << n_sig << endl;
	}
	toc = measureTime();
	time_postprocessing += (toc-tic);

	// Record final timestamp
	t_end = measureTime();

	// Report code profiling measurements
	profileCode();

	// Close files (if they are open)
	if(summary_file_obj.is_open()) summary_file_obj.close();
	if(profiling_file_obj.is_open()) profiling_file_obj.close();
	if(compute_pvals) sig_itemsets_file.close();
	if(output_testable) testable_itemsets_file.close();

	return 0;
}
