// AUTHOR: FELIPE LLINARES

#ifndef _bonf_cmh_cpp_
#define _bonf_cmh_cpp_


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


// ---------------------------CONSTANT DEFINES -------------------------------------


#define NGRID 500 //Number of samples in the grid of tentative corrected significance thresholds, not counting 1 as it is a trivial threshold
#define LOG10_MIN_PVAL -30.0 //Minimum tentative corrected significance threshold = 10^{LOG10_MIN_PVAL}


// -----------------------------GLOBAL VARIABLES----------------------------------

// Number of different categories for the covariate
int n_cat;
// Category (integer between 0 and n_cat-1) of each sample
vector<int> cats;

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


// If compute_pvals==true, pvalues for testable itemsets are compute on the fly, using a greedy choice
// for the significance threshold.
bool compute_pvals;
ofstream sig_itemsets_file;


// ----------------------------INITIALIZATION FUNCTIONS--------------------------------------------


/* Initialize Tarone related global variables and constants */
void init_tarone(double fwer, int n_samples){
	int n_samples_over_2 = (n_samples % 2) ? (n_samples-1)/2 : n_samples/2;  //floor(n_samples/2)
    // Initialize some constants
    target_fwer = fwer;


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

}

void set_significance_threshold(){
	delta_opt = target_fwer/n_enumerated_closed;
}


// ----------------------------MAIN TARONE FUNCTIONALITY------------------------------------------




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

	// Close output file
	sig_itemsets_file.close();
	// Close temporary file and delete it
	sig_itemsets_file_tmp_read.close();
	remove((filename + string(".tmp")).c_str());

	return n_sig;
}

#endif
