// AUTHOR: FELIPE LLINARES

#ifndef _tarone_chi_cpp_
#define _tarone_chi_cpp_


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
#include "shared_vars.h"

using namespace std;

// -----------------------------GLOBAL VARIABLES----------------------------------

// Number of testable itemsets
long long n_testable;

// Region thresholds: Sigma_k = [sl1,sl2] U [su1,su2]
// We always have su1=n_samples-sl2 and su2=n_samples-sl1, but keeping each variable separate
// reduces the number of computations at the expense of a tiny amount of memory
int sl1, sl2, su1, su2;

// Flag variable to keep track of the last change done to region Sigma_k
// If flag==1, the last change was sl1++ (shrink on extremes of the W)
// If flag==0, the last change was sl2-- (shrink on center of the W)
bool W_flag;

// Target FWER
double target_fwer;

// Current P-value threshold
double delta;

// Final corrected significance threshold
double delta_opt;

// Array with all values of minimum attainable P-value in [0,n_samples] pre-computed
vector<double> psi;
// And some constants which are useful to precompute
double class_ratio, class_ratio_bin;

// A (n_samples+1)-dimensional vector such that
// freq_cnt[j] = #itemsets with support=j processed so far
vector<long long> freq_cnt;

// If compute_pvals==true, pvalues for testable itemsets are compute on the fly, using a greedy choice
// for the significance threshold.
bool compute_pvals;
vector<double> tentative_sig_ths;
ofstream sig_itemsets_file;


bool output_testable;
ofstream testable_itemsets_file;

// ----------------------------INITIALIZATION FUNCTIONS--------------------------------------------


/* Precompute minimum attainable P-values $\psi(x)$ for all x in [0,n_samples] and store them in array psi */
void psi_init(){
	double num, den;
	int x;
	int n_samples_over_2 = (n_samples % 2) ? (n_samples-1)/2 : n_samples/2;  //floor(n_samples/2)

	// Allocate memory for psi
	psi.resize(n_samples+1);

	/* Initialize caches with appropriate values */

	// Precompute some useful constants
	class_ratio = ((double)n1)/n_samples; class_ratio_bin = class_ratio*(1-class_ratio);

	// First compute the left side of "the W", i.e. the range [0,n]
	psi[0] = 1;
	for(x=1; x<=n1; x++) {
		num = x*(1-class_ratio); num *= num;
		den = x*(1-((double)x)/n_samples)*class_ratio_bin;
		psi[x] = Chi2_sf(num/den,1);
	}

	// Now, compute the minimum attainable p-values in the range [N-N_over_2,N]
	for(x=(n1+1); x<=n_samples_over_2; x++) {
		num = n1*(1-((double)x)/n_samples); num *= num;
		den = x*(1-((double)x)/n_samples)*class_ratio_bin;
		psi[x] = Chi2_sf(num/den,1);
	}

	// Finally, since $\psi(x)$ is symmetric around N_over_2, complete the right half by symmetry
	for(x=(n_samples_over_2+1); x<=n_samples; x++) psi[x] = psi[n_samples-x];
}


/* Initialize Tarone related global variables and constants */
void init_tarone(double fwer){
	int n_samples_over_2 = (n_samples % 2) ? (n_samples-1)/2 : n_samples/2;  //floor(n_samples/2)
    // Initialize some constants
    target_fwer = fwer;
    sl1 = 1; sl2 = n_samples_over_2; su1 = n_samples-sl2; su2 = n_samples-sl1;
    W_flag = true;
    // Set number of testable itemsets initially to 0
    n_testable = 0;

    // Initialize cache for psi(x)
    psi_init();
	// Initialize corrected significance threshold
    delta = psi[sl1]; //$\psi(1)$

    // Allocate memory for support histogram
    freq_cnt.resize(n_samples+1);

    // Push the first tentative significance threshold for greedy p-value evaluation procedure
    if(compute_pvals) tentative_sig_ths.push_back(fwer);
}


// ----------------------------MAIN TARONE FUNCTIONALITY------------------------------------------

/* Decrease the minimum p-value threshold one level
 * Main operations that need to be performed are:
 * 1) Figure out whether we have to shrink "the W" on the left side or the right side, that is, if the current region
 *    is Sigma_{k} = [sl1,sl2] U [n_individuals-sl2,n_individuals-sl1], we need to figure out if Sigma_{k+1} is of
 *    the form Sigma_{k+1} = [sl1+1,sl2] U [n_individuals-sl2,n_individuals-sl1-1] (shrink left side) or
 *    Sigma_{k+1} = [sl1,sl2-1] U [n_individuals-sl2+1,n_individuals-sl1-1] (shrink right side).
 *    This is done with help of a binary flag that remembers which of the two types of region change happened the
 *    last time the threshold was decreased.
 * 2) Update variables sl1, sl2 and delta accordingly
 * 3) Recompute the number of testable items by removing those being excluded from the testable region due to the
 *    threshold change
 * */
void decrease_threshold(){
	if(W_flag){ // W_flag==true means the last call to decrease_threshold() shrunk "the W" on the left side
		// Update number of testable intervals
		n_testable -= freq_cnt[sl1]; n_testable -= freq_cnt[su2];
		sl1++; su2--; // Shrink Sigma_k on extremes of the W
		// Check what the new case will be
		if (psi[sl1] >= psi[sl2]) delta = psi[sl1];
		else{ delta = psi[sl2]; W_flag = 0; }
	}
	else{ // W_flag==0 means the last call to decrease_threshold() shrunk "the W" on the right side
		// Update number of testable intervals
		if(sl2==su1) n_testable -= freq_cnt[sl2];//(beware of case sl2==su1 since it could lead to discounting the same thing twice!)
		else {n_testable -= freq_cnt[sl2]; n_testable -= freq_cnt[su1];}
		sl2--; su1++; // Shrink Sigma_k on center of the W
		// Check what the new case will be
		if (psi[sl1] >= psi[sl2]){ delta = psi[sl1]; W_flag = 1; }
		else delta = psi[sl2];
		//No need to update minimum support in this case, since sl1 remains the same
	}
	if(compute_pvals) tentative_sig_ths.push_back(target_fwer/n_testable);
}


// Return true if support x is testable
#define IS_TESTABLE(x) (psi[(x)] <= delta)


// Evaluate $\chi^2$ p-value
double compute_pval(int a, int x){
	double aux, num, den;
	aux = ((double)x)/n_samples;
	num = a-n1*aux; num = pow(num,2);
	den = x*(1-aux)*class_ratio_bin;
	return Chi2_sf(num/den,1);
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
