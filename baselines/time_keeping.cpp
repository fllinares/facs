#ifndef _time_keeping_cpp_
#define _time_keeping_cpp_

/* LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */
#include <fstream>
#include <time.h>
#include <sys/resource.h>
/* END OF LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */

using namespace std;


/* GLOBAL VARIABLES (TIME SPENT) */


double t_init = 0, t_end = 0;
double time_IO = 0;
double time_initialisation = 0;
double time_enumeration = 0;
double time_postprocessing = 0;
//double time_pruning_check = 0;
double tic,toc;
//double tic2,toc2;


// Output stream to write profiling measurements into. Set to stdout by default,
// unless optional inputs override the behaviour
ostream *timing_file = &cout;


// -------------------------------PROFILING-FUNCTIONS-----------------------------


// Measure running time
//double measureTime(){
//  struct rusage t;
//  struct timeval tv, ts;
//  getrusage(RUSAGE_SELF, &t);
//  tv = t.ru_utime;
//  ts = t.ru_stime;
//  return tv.tv_sec + ts.tv_sec + ((double)tv.tv_usec + (double)ts.tv_usec) * 1e-6;
//}

double measureTime(){
	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	return ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}


// Measure peak memory usage
size_t measurePeakMemory(){
  struct rusage t;
  getrusage(RUSAGE_SELF, &t);
  return (size_t)t.ru_maxrss;
}


// ----------------------------------I/O-FUNCTIONS--------------------------------


// Display execution time and memory consumption
void profileCode(){
	size_t peak_memory = measurePeakMemory();

	*timing_file << "CODE PROFILING" << endl;
	*timing_file << "Total execution time: " << (t_end-t_init) << " (s)." << endl;
	*timing_file << "\t" << "File I/O time: " << time_IO << " (s)." << endl;
	*timing_file << "\t" << "Initialisation time: " << time_initialisation << " (s)." << endl;
	*timing_file << "\t" << "Enumeration time: " << time_enumeration << " (s)." << endl;
	//*timing_file << "\t\t" << "Pruning check time: " << time_pruning_check << " (s)." << endl;
	*timing_file << "\t" << "Postprocessing time: " << time_postprocessing << " (s)." << endl;
	*timing_file << "Peak memory usage: " << peak_memory << " (bytes)." << endl;
}
#endif
