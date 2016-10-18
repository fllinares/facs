// AUTHOR: FELIPE LLINARES

#ifndef _tarone_cmh_2k_cpp_
#define _tarone_cmh_2k_cpp_


// -----------------------------INCLUDES-----------------------------------------


#include <stdio.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <algorithm>
#include <vector>
#include <unistd.h>
#include <math.h>

#include "chi2.h"
#include "time_keeping.cpp"
#include "shared_vars.h"

using namespace std;


// ---------------------------CONSTANT DEFINES -------------------------------------


#define NGRID 500 //Number of samples in the grid of tentative corrected significance thresholds, not counting 1 as it is a trivial threshold
#define LOG10_MIN_PVAL -30.0 //Minimum tentative corrected significance threshold = 10^{LOG10_MIN_PVAL}


// -----------------------------GLOBAL VARIABLES----------------------------------

// Number of different categories for the covariate
int n_cat;
// Category (integer between 0 and n_cat-1) of each sample
vector<int> cats;

// Number of testable itemsets
long long n_testable;

// Target FWER
double target_fwer;

// Current P-value threshold
double delta;

// Final corrected significance threshold
double delta_opt;

// Number of observations per table
vector<int> n_samples_t;
// Number of observations in positive class per table
vector<int> n1_t;
// Now some precomputed quantities to save time when computing the CMH test statistic
vector<int> n2_t; // Number of observations in negative class per table
vector<int> hypercorner_bnd; // $min(n_{1,j},n_{2,j})$ for each of the n_cat tables
vector<double> gamma_t; // n_{1,j}/n_j for each of the n_cat tables
vector<double> gammabin_t; // (n_{1,j}/n_j)*(1 - n_{1,j}/n_j) for each of the n_cat tables
// And some others to save time when computing the maximum CMH test statistic in the "bottom left" hypercorner
vector<double> f_vals, g_vals, betas;
double f_sum, g_sum, Tcmh_max_corner_l, Tcmh_max_corner_r, Tcmh_aux_corner;
vector<int> idx_betas_sorted;

// Auxiliary vector to evaluate the min p-value lower envelope via brute force
vector<unsigned long long> bitmasks; //bitmasks[k] = 2^k


// Grid of logarithmically spaced corrected significance thresholds. The sequence is ordered
// from larger thresholds to smaller thresholds.
vector<double> pgrid;
// Current tentative corrected significance threshold and index of it in the grid
double pth;
int idx_th;
// Step size in the grid
double log10_p_step;


// A (n_samples+1)-dimensional vector such that
// freq_cnt[j] = #itemsets with support=j processed so far
vector<long long> freq_cnt;

// If compute_pvals==true, pvalues for testable itemsets are compute on the fly, using a greedy choice
// for the significance threshold.
bool compute_pvals;
vector<double> tentative_sig_ths;
ofstream sig_itemsets_file;


// ----------------------------INITIALIZATION FUNCTIONS--------------------------------------------


/* Initialize Tarone related global variables and constants */
void init_tarone(double fwer, int n_samples){
	int j;
	double log10_p;
	int n_samples_over_2 = (n_samples % 2) ? (n_samples-1)/2 : n_samples/2;  //floor(n_samples/2)
    // Initialize some constants
    target_fwer = fwer;
    // Set number of testable itemsets initially to 0
    n_testable = 0;

    // Initialize grid of candidate corrected significance thresholds
    //pgrid.resize(NGRID+1);
	for(log10_p=0,log10_p_step=-LOG10_MIN_PVAL/NGRID,j=0; j<=NGRID; log10_p-=log10_p_step, j++) pgrid.push_back(pow(10,log10_p));
	// Initialize threshold values
	idx_th = 1; pth = pgrid[idx_th];


	// Compute number of observations, and number of observations in each class per category
	n_samples_t.resize(n_cat); n1_t.resize(n_cat); n2_t.resize(n_cat);
	for(int i=0; i<n_samples; ++i){
		n_samples_t[cats[i]]++;
		if(labels[i]) n1_t[cats[i]]++;
		else n2_t[cats[i]]++;
	}

	// Compute auxiliary quantities for fast evaluation of CMH test and its pruning criterion
	hypercorner_bnd.resize(n_cat); gamma_t.resize(n_cat); gammabin_t.resize(n_cat);
	for(int c=0; c<n_cat; ++c){
		hypercorner_bnd[c] = (n1_t[c] <= n2_t[c]) ? n1_t[c] : n2_t[c]; //min(n1_t[c],n2_t[c])
		gamma_t[c] = ((double)n1_t[c])/n_samples_t[c];
		gammabin_t[c] = gamma_t[c]*(1-gamma_t[c]);
	}

	// Allocate memory for p-value and minimum attainable p-value computation "registers"
	f_vals.resize(n_cat); g_vals.resize(n_cat); betas.resize(n_cat); idx_betas_sorted.resize(n_cat);

    // Allocate memory for support histogram
    freq_cnt.resize(NGRID+1);

    // Initialize the bitmasks for brute-force evaluation of the minimum attainable p-value lower envelope
	bitmasks.push_back(1);
	for(int c=1; c<n_cat; ++c) bitmasks.push_back(2*bitmasks[c-1]);

    // Push the first tentative significance threshold for greedy p-value evaluation procedure
    if(compute_pvals) tentative_sig_ths.push_back(fwer);
}


// ----------------------------MAIN TARONE FUNCTIONALITY------------------------------------------


/* Decrease the minimum p-value threshold one level
 */
void decrease_threshold(){
	// Remove the intervals which become untestable after the change
	n_testable -= freq_cnt[idx_th];
	// Change threshold
	idx_th++; pth = pgrid[idx_th];

	if(compute_pvals) tentative_sig_ths.push_back(target_fwer/n_testable);
}


// Return true if itemset is testable
#define IS_TESTABLE(minpval) (minpval <= pth)


// Return true if itemset x is not prunable
#define IS_NOT_PRUNABLE(lower_envelope_minpval) (lower_envelope_minpval <= pth)


/* Computes the CMH p-value as a function of the margins x, n1 and n and the cell counts a for the n_cat tables */
double compute_pval(int a, vector<int> &x_t){
	double num = a, den = 0;
	for(int c=0; c<n_cat; ++c){
		num -= x_t[c]*gamma_t[c];
		den += x_t[c]*(1-((double)x_t[c])/n_samples_t[c])*gammabin_t[c];
	}
	num *= num; //num = num^2
	if(den==0) return 1;
	else return Chi2_sf(num/den,1);
}


/* Computes the minimum attainable CMH p-value depending on the margins x, n1 and n for the n_cat tables */
double compute_minpval(vector<int> &x_t){
	double left_tail_num = 0, right_tail_num = 0, den = 0;
	double aux1, aux2;
	for(int c=0; c<n_cat; ++c){
		aux1 = x_t[c]-n2_t[c]; aux2 = x_t[c]*gamma_t[c];
		left_tail_num += ((aux1 > 0) ? aux1 : 0) - aux2;
		right_tail_num += ((x_t[c] > n1_t[c]) ? n1_t[c] : x_t[c]) - aux2;
		den += x_t[c]*(1-((double)x_t[c])/n_samples_t[c])*gammabin_t[c];
	}
	left_tail_num *= left_tail_num; right_tail_num *= right_tail_num;
	if(den==0) return 1;
	else return Chi2_sf(((left_tail_num > right_tail_num) ? left_tail_num : right_tail_num)/den,1);
}

double compute_lower_envelope_minpval(vector<int> &x_t){
	double lower_envelope_minpval = 1;
	double minpval;
	// Variables for looping across 2^K cases
	unsigned long long idx_mask, max_mask;
	vector<int> z_in(n_cat);

	// If for any of the n_cat tables, its margin x is larger than the minimum of n1 and n_samples-n1, then
	// we cannot prune the itemset (we are not in the "bottom-left" hypercorner)
	for(int c=0; c<n_cat; ++c) if(x_t[c] > hypercorner_bnd[c]) return -1;

	//tic2 = measureTime();

	max_mask = 1 << n_cat; Tcmh_max_corner_l = 0;
	for(idx_mask=1;idx_mask<max_mask;idx_mask++){
		for(int c=0; c<n_cat; ++c){
			// Skip if the bitmask contains a zero
			if(idx_mask & bitmasks[c]) z_in[c] = x_t[c];
			else z_in[c] = 0;
		}
		minpval = compute_minpval(z_in);
		lower_envelope_minpval = (lower_envelope_minpval <= minpval) ? lower_envelope_minpval : minpval;
	}

	//toc2 = measureTime();
	//time_pruning_check += (toc2-tic2);

	return lower_envelope_minpval;
}



//double compute_lower_envelope_minpval(vector<int> &x_t){
//	// Variables for looping across 2^K cases
//	unsigned long long idx_mask, max_mask;
//	// Variables for computing minimum attainable p-value in each case
//	double left_tail_num, right_tail_num, den;
//	double aux;
//
//	// If for any of the n_cat tables, its margin x is larger than the minimum of n1 and n_samples-n1, then
//	// we cannot prune the itemset (we are not in the "bottom-left" hypercorner)
//	for(int c=0; c<n_cat; ++c) if(x_t[c] > hypercorner_bnd[c]) return -1;
//
//	max_mask = 1 << n_cat; Tcmh_max_corner_l = 0;
//	for(idx_mask=1;idx_mask<max_mask;idx_mask++){
//		left_tail_num = 0, right_tail_num = 0, den = 0;
//		for(int c=0; c<n_cat; ++c){
//			// Skip if the bitmask contains a zero
//			if(idx_mask & bitmasks[c]){
//				aux = x_t[c]*gamma_t[c];
//				left_tail_num -= aux;
//				right_tail_num += x_t[c] - aux;
//				den += x_t[c]*(1-((double)x_t[c])/n_samples_t[c])*gammabin_t[c];
//			}
//		}
//		left_tail_num *= left_tail_num; right_tail_num *= right_tail_num;
//		if(den==0) Tcmh_aux_corner = 0;
//		else Tcmh_aux_corner = ((left_tail_num > right_tail_num) ? left_tail_num : right_tail_num)/den;
//		Tcmh_max_corner_l = (Tcmh_max_corner_l >= Tcmh_aux_corner) ? Tcmh_max_corner_l : Tcmh_aux_corner;
//	}
//	aux = Chi2_sf(Tcmh_max_corner_l,1);
//	cout << aux << endl;
//	return Chi2_sf(Tcmh_max_corner_l,1);
//}


// Map p-value or minimum attainable p-value to its corresponding bucket in the grid of threshold candidates
inline int bucket_idx(double pval){
	int idx;
	idx = (int)floor(-log10(pval)/log10_p_step);
	if(idx<0) idx = 0;
	if(idx>NGRID) idx = NGRID;
	return idx;
}


// Process the greedy p-value evaluation structure
long long process_significant_itemsets(string filename){
	// String to store each line read from the file
	string line;
	// And a Stringstream to parse it
	stringstream ss_line;
	// Double to store p-value in each line of the file
	double pval;
	// Integer to store each value read from each field in a line of the file
	int val;
	// Number of significant itemsets
	long long n_sig = 0;

	// Close temporary file containing testable itemsets which might be significant
	sig_itemsets_file.close();
	// And open same file again in read mode
	ifstream sig_itemsets_file_tmp_read(filename + string(".tmp"));
	if(sig_itemsets_file_tmp_read.is_open()){
		// Open final file containing significant itemsets which are actually significant
		sig_itemsets_file.open(filename);
		if(sig_itemsets_file.is_open()){
			// Process file line by line
			while(getline(sig_itemsets_file_tmp_read,line)){
				ss_line << line;
				// Retrieve p-value (first element in line)
				ss_line >> pval;
				// Retrieve cell count a
				//ss_line >> val;
				// Retrive margin x
				//ss_line >> val;
				while(ss_line >> val); // Do stuff with itemset
				// Check if itemset is truly significant. If it is, write it to final file
				if(pval <= delta_opt){
					sig_itemsets_file << line << "\n";
					n_sig++;
				}
				// stringstream flags for the next iteration
				ss_line.clear();
			}
		}
		else{
			cerr << "Error @ process_significant_itemsets: Unable to open output significant itemsets file " << filename << endl;
			exit(-1);
		}
	}
	else{
		cerr << "Error @ process_significant_itemsets: Unable to open temporal significant itemsets file " << filename << ".tmp" << endl;
		exit(-1);
	}

	// Output list of tentative significance thresholds for debugging purposes
	sig_itemsets_file << "\n" << "TENTATIVE SIGNIFICANCE THRESHOLDS:" << "\n";
	for(int i=0; i < tentative_sig_ths.size()-1; ++i) sig_itemsets_file << tentative_sig_ths[i] << "\t";
	sig_itemsets_file << tentative_sig_ths[tentative_sig_ths.size()-1] << endl;
	// Finally, check if some of the tentative significance thresholds were too strict, and thus patterns might have been missed
	double min_tentative_sig_th = *min_element(tentative_sig_ths.begin(), tentative_sig_ths.end());
	if(min_tentative_sig_th < delta_opt) sig_itemsets_file << "WARNING!: The greedy p-value evaluation approach used threshold " << min_tentative_sig_th << " during the procedure. Some significant patterns might be lost." << endl;
	else sig_itemsets_file << "Greedy p-value evaluation approach was successful. Have a nice day :)" << endl;
	// Close output file
	sig_itemsets_file.close();
	// Close temporary file and delete it
	sig_itemsets_file_tmp_read.close();
	remove((filename + string(".tmp")).c_str());

	return n_sig;
}

#endif
