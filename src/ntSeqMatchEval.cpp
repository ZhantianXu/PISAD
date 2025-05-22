/*
 * Compares count lists using Fisher's Exact Test per pair and produces a table
 * with the combinations of all pairs
 */

#include "CompareCounts.hpp"
#include "src/Options.h"
#include "src/Util.h"
#include <assert.h>
#include <getopt.h>
#include <iostream>
#include <limits.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

#include <omp.h>

using namespace std;

#define PROGRAM "pisadEval"

void printHelpDialog() {
  const string dialog =
      "Usage: " PROGRAM "[OPTION]... [FILES...]\n"
      "Optional options:\n"
      "  -t, --threads = INT    Number of threads to run.[1]\n"
      "  -h, --help             Display this dialog.\n";
  cerr << dialog << endl;
  exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[]) {
  // switch statement variable
  int c;

  // control variables
  bool die = false;
  int OPT_VERSION = 0;

  // long form arguments
  static struct option long_options[] = {
      {"threads", required_argument, NULL, 't'}, {NULL, 0, NULL, 0}};

  int option_index = 0;
  while ((c = getopt_long(argc, argv, "t:h", long_options, &option_index)) !=-1) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 't': {
      stringstream convert(optarg);
      if (!(convert >> opt::threads)) {
        cerr << "Error - Invalid parameter t: " << optarg << endl;
        return 0;
      }
      break;
    }
    case 'h': {
      printHelpDialog();
      break;
    }
    case '?': {
      die = true;
      break;
    }
    }
  }

  vector<string> inputFiles;
  while (optind < argc) {
    inputFiles.emplace_back(argv[optind]);
    optind++;
  }

  int existingFileCount = 0;
  for (const auto &file : inputFiles) {
    if (Util::fexists(file)) {
      existingFileCount++;
    }
  }

  if (existingFileCount == 0) {
    cerr << "Error: Input File is wrong" << endl;
    die = true;
  }

  // Check needed options
  if (inputFiles.size() == 0) {
    cerr << "Error: Need Input File" << endl;
    die = true;
  }
  if (die) {
    cerr << "Try '--help' for more information.\n";
    exit(EXIT_FAILURE);
  }

  double time = omp_get_wtime();

  CompareCounts comp(inputFiles);

  comp.computeScore();

  cerr << "Time: " << omp_get_wtime() - time << " s Memory: " << Util::getRSS()
       << " kbytes" << endl;
  return 0;
}
