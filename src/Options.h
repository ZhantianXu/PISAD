/*
 * Options.h
 *
 *  Created on: Oct 14, 2020
 *      Author: cjustin
 */

#ifndef OPTIONS_H
#define OPTIONS_H 1

#include <limits>
#include <stdint.h>
#include <string>
#include <vector>

using namespace std;

/**
 * Global variables that are mostly constant for the duration of the
 * execution of the program.
 */
namespace opt {

int verbose = 0;
unsigned threads = 1;
unsigned k = 21;
unsigned coverage = 30;
double scoreThresh = 0.63;

std::vector<std::string> snp;
string summary = "";
double covThresh = std::numeric_limits<double>::max();

unsigned minCov = 1;
bool information = false;
uint64_t genomeSize = 6200000000;

string ref;
unsigned window = 31;
unsigned multi = 20;

string debug = "";
string name = "";
} // namespace opt
#endif
