# FACS

**F**ast **A**utomatic **C**onditional **S**earch

This is a C implementation of the algorithm described in the following paper:

+ L. Papaxanthos\*, F. Llinares-López\*, D. Bodenham and K. Borgwardt (2016)
**Finding significant combinations of features in the presence of categorical covariates**, *NIPS 2016*

The paper can be found [here](http://papers.nips.cc/paper/6345-finding-significant-combinations-of-features-in-the-presence-of-categorical-covariates). 

A short spotlight video summarising the main aspects of the paper can be found [here](https://www.youtube.com/watch?v=_E6VsjvxVdA).

## Compilation
It is possible to compile the C source code by using the provided makefile:
```
$ cd ./facs
$ make
```
By default, the code will be compiled using ```gcc``` and the resulting executable files will be stored in the folder compiled. If available, using Intel's ```icc``` instead of ```gcc``` might lead to improved performance.

## Usage

After compilation, the folder compiled will contain all executable files. The executable for FACS will be named **facs** by default. By executing it with no input arguments, a list of all mandatory and optional input arguments will be displayed:

```
$ ./compiled/facs
> Error @ main: Incorrect arguments!!!
> Usage:	./compiled/facs [-option] [option_argument]
> Mandatory options:
> 	option:	-t path_to_transactions_file
> 	option:	-l path_to_labels_file
> 	option:	-c path_to_covariates_file
> 	option:	-f target_FWER
> Optional options:
> 	option:	-o path_to_results_summary_file
> 		Default: results summary displayed in stdout
> 	option:	-p path_to_profiling_summary_file
> 		Default: profiling results displayed in stdout
> 	option:	-s path_to_significant_itemsets_file
> 		Default: significance of itemsets is not evaluated at all, only testability
> 	option:	-d path_to_testable_itemsets_file
> 		Note: requires usage of option -s as well!
> 		Default: testable itemsets not shown in output
```

The input data format used by FACS represents each dataset using three separate files: 

+ The *transactions file* contains the (binary) feature matrix. Each row corresponds to a different sample in the dataset, represented as a whitespace-separated list of the indices corresponding to non-zero features in the sample.
+ The *labels file* contains the (binary) class labels. Each row corresponds to a different sample in the dataset and must be either 0 or 1.
+ The *covariates file* contains the (categorical) covariates. Each row corresponds to a different sample in the dataset and must be a value in 0, 1, ..., k-1, where k is the number of distinct categories for the covariate.

The folder example_data contains two toy datasets, tictactoe and mushroom, in order to illustrate the input data format used by FACS.

As an illustration, FACS can be run on the tictactoe data as follows:

```
$ mkdir output_tictactoe
$ ./compiled/facs -t example_data/tictactoe_transactions.dat -l example_data/tictactoe_labels.dat -c example_data/tictactoe_covariates.dat -f 0.05 -o output_tictactoe/summary.txt -p output_tictactoe/profiling.txt -s output_tictactoe/significant_itemsets.txt -d output_tictactoe/testable_itemsets.txt

```

Upon execution, the folder output_tictactoe will now contain four output files:

+ The *summary file*, specified with option ```-o```, displays general information, such as dataset size, resulting Tarone's testability threshold, the corrected significance threshold and the total number of significantly associated feature combinations found in the data, among others.
+ The *profiling file*, specified with option ```-p```, displays the total runtime and memory usage of the execution of FACS, as well as a breakdown of the runtime for different key steps of the algorithm.
+ The *significant itemsets file*, specified with option ```-s```, contains all feature combinations found by FACS to be significantly associated with the class labels given the categorical covariate. Each row of the file corresponds to a different significant feature combination. The first entry of the row is the p-value of the CMH test for the feature combination; all remaining entries are the indices of the features in the combination. Note that feature indices will not be necessarily ordered.
+ The *testable itemsets file*, specified with option ```-d```, contains the CMH test p-values of all testable feature combinations, include those not deemed significant. While this is useful for certain purposes (e.g. to evaluate inflation due to confounding) the resulting file can be rather large. Therefore, we recommend not using the option unless necessary.

## Baselines

By default, the provided makefile will also compile the different baseline algorithms described in our paper. The corresponding executables are:

+ **facs_2k** and **facs_Mk**: The executables for two baseline versions of FACS that do **not** use our novel approach to compute the minimum attainable p-value lower envelope. While the computational complexity of FACS scales as O(k log(k)), the naive baseline implementations in eclat_lamp_cmh_2k and eclat_lamp_cmh_Mk scale as O(2^k) and O(m^k), m >> 2, respectively.
+ **bonf_cmh**: The executable for a baseline approach that uses the CMH test along a naive Bonferroni correction. As FACS, the use of the CMH test allows to correct for a confounding categorical covariate. However, this baseline does not use Tarone's testability criterion, leading to a decrease in statistical power and computational efficiency.
+ **lamp_chi**: An implementation of the state-of-the-art significant pattern mining algorithm [LAMP](http://a-terada.github.io/lamp/) that uses Pearson's Chi-squared test. As FACS, LAMP uses Tarone's testability criterion, leading to improved statistical power and drastically reduced runtime compared to Bonferroni-based approaches. However, it is unable to correct for confounding covariates.

## Contact 

Any questions can be directed to Laetitia Papaxanthos or Felipe Llinares-López: {laetitia.papaxanthos, felipe.llinares} [at] bsse.ethz.ch  



