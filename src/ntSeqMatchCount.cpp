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

#include "FingerPrint.hpp"

#if _OPENMP
#include <omp.h>
#endif

using namespace std;

#define PROGRAM "pisadCount"

void printHelpDialog() {
  const char dialog[] =
      "Usage: " PROGRAM " -s [FASTA] [OPTION]... [FILES...]\n"
      "Required options:\n"
      "  -s, --snp = STR        variant sketch (one or more) [required]\n"
      "Optional options:\n"
      "  -t, --threads = INT    Number of threads to run.[1]\n"
      "  -m, --maxCov = INT     k-mer coverage threshold for early "
      "termination. [inf]\n"
      "  -i, --information      extra debug information.\n"
      "  -k, --kmer = INT       k-mer size used. [21]\n"
      "  -h, --help             Display this dialog.\n"
      "  -o, --output           Evaluation file path\n";

      cerr << dialog << endl;
  exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[]) {
  // switch statement variable
  int cov;

  // control variables
  bool die = false;
  int OPT_VERSION = 0;
  bool valid = false;
  // long form arguments
  static struct option long_options[] = {
      {"threads", required_argument, NULL, 't'},
      {"maxCov", required_argument, NULL, 'm'},
      {"snp", required_argument, NULL, 's'},
      {"kmer", required_argument, NULL, 'k'},
      {"help", no_argument, NULL, 'h'},
      {"version", no_argument, &OPT_VERSION, 1},
      {"information", no_argument, NULL, 'i'},
      {"output", required_argument, NULL, 'n'},
      {NULL, 0, NULL, 0}};

  int option_index = 0;
  while ((cov = getopt_long(argc, argv, "s:t:o:ihk:m:", long_options,
                            &option_index)) != -1) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (cov) {
    case 'h': {
      printHelpDialog();
      break;
    }
    case 'i': {
      opt::information = true;
      break;
    }
    case 'o': {
      stringstream convert(optarg);
      if (!(convert >> opt::name)) {
        cerr << "Error - Invalid parameter n: " << optarg << endl;
        return 0;
      }
      break;
    }
    case 's': {
      stringstream ss(optarg);
      string file;
      while (getline(ss, file, ',')) {
        if (!file.empty()) {
          opt::snp.push_back(file);
          valid = true;
        }
      }
      break;
    }
    case 'm': {
      stringstream convert(optarg);
      if (!(convert >> opt::covThresh)) {
        cerr << "Error - Invalid parameter m: " << optarg << endl;
        return 0;
      }
      break;
    }
    case 'k': {
      stringstream convert(optarg);
      if (!(convert >> opt::k)) {
        cerr << "Error - Invalid parameter k: " << optarg << endl;
        return 0;
      }
      break;
    }
    case 't': {
      stringstream convert(optarg);
      if (!(convert >> opt::threads)) {
        cerr << "Error - Invalid parameter t: " << optarg << endl;
        return 0;
      }
      break;
    }
    case '?': {
      die = true;
      break;
    }
    }
  }

#if defined(_OPENMP)
  if (opt::threads > 0)
    omp_set_num_threads(opt::threads);
#endif

  if (OPT_VERSION) {
    // printVersion();
  }

  if (opt::k > 32) {
    die = true;
    cerr << "k cannot be greater than 32" << endl;
  }

  vector<string> inputFiles;
  while (optind < argc) {
    inputFiles.emplace_back(argv[optind]);
    assert(Util::fexists(inputFiles.back()));
    optind++;
  }

  string extracted = "";
  // Check needed options
  if (!inputFiles.empty()) {
    string filename = inputFiles[0];
    cout << filename << endl;

    size_t lastSlash = filename.find_last_of("/");
    size_t lastDot = std::string::npos;
    const std::string exts[] = {".fastq.gz", ".fq.gz", ".fastq", ".fq"};
    for (const auto &ext : exts) {
      size_t pos = filename.rfind(ext);
      if (pos != std::string::npos && (lastDot == std::string::npos || pos > lastDot)) {
        lastDot = pos;
      }
      }
    // size_t lastDot = filename.find_last_of(".");

    if (lastSlash != string::npos && lastDot != string::npos &&
        lastSlash < lastDot) {
      extracted = filename.substr(lastSlash + 1, lastDot - lastSlash - 1);
      cout << "Extracted string: " << extracted << endl;
    } else {
      cout << "Invalid filename format!" << endl;
      die = true;
    }
  } else {
    cout << "Error: No input files!" << endl;
    die = true;
  }

  if (!valid) {
    std::cerr
        << "Error: -s option requires a valid comma-separated list of files."
        << std::endl;
    die = true;
  }

  if (die) {
    cerr << "Try '--help' for more information.\n";
    exit(EXIT_FAILURE);
  }

  double time = omp_get_wtime();

  FingerPrint fp;
  cerr << "read fa" << endl;
  cerr << "Time: " << omp_get_wtime() - time
       << " s Memory: " << (Util::getRSS() / (1024 * 1024)) << " G" << endl;
  fp.computeCounts(inputFiles);
  // fp.computeCountsProducerConsumer(inputFiles);
  // fp.printOptionalHeader();
  // fp.printCountsMax();
  // fp.printCountsAllCounts();
  cerr << "begin print" << endl;
  cerr << "Time: " << omp_get_wtime() - time
       << " s Memory: " << (Util::getRSS() / (1024 * 1024)) << " G" << endl;
  fp.printCountsMode(extracted);
  cerr << fp.printInfoSummary() << endl;
  cerr << "Time: " << omp_get_wtime() - time
       << " s Memory: " << (Util::getRSS() / (1024 * 1024)) << " G" << endl;
  return 0;
}
