#ifndef CONFIG_H
#define CONFIG_H

#include <cfloat>
#include <string>

enum calculation_method { LINKFIND, COMPUTEPAIRS, ALEXANDER};
enum output_location { LOC_NONE, LOC_BIRTH, LOC_DEATH};
enum file_format { DIPHA, PERSEUS, NUMPY };


struct Config {
	std::string filename = "";
	std::string output_filename = "output.csv"; //default output filename
	file_format format;
//    calculation_method method = ALEXANDER;
	calculation_method method = LINKFIND;
	double threshold = DBL_MAX;
	int maxdim=2;;  // compute PH for these dimensions
	bool print = false; // flag for printing to std
	bool embedded = false; // embed image in the sphere (for alexander duality)
	output_location location = LOC_BIRTH; // flag for saving location
	int min_cache_size = 0; // num of minimum non-zero entries of a reduced column to be cached
};

#endif
