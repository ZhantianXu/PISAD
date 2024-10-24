/*
 * Compares count lists using Fisher's Exact Test per pair and produces a table
 * with the combinations of all pairs
 */

#include <sstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include "src/Options.h"
#include "src/Util.h"
#include "CompareCounts.hpp"

#include <omp.h>

using namespace std;

#define PROGRAM "ntsmEval"

// void printVersion()
// {
// 	const char VERSION_MESSAGE[] = PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
// 										   "Written by Justin Chu <cjustin@ds.dfci.harvard.edu>\n"
// 										   "\n"
// 										   "Copyright 2021 Dana-Farber Cancer Institute\n";
// 	cerr << VERSION_MESSAGE << endl;
// 	exit(EXIT_SUCCESS);
// }

void printHelpDialog()
{
	const string dialog =
		"Usage: " PROGRAM " [FILES...]\n"
		"Processes sets of counts files and compares their similarity.\n"
		"The first file is the statistics of the target file, followed by the statistics of one or more other files.\n"
		"  -t, --threads              Number of threads to run.\n";

	cerr << dialog << endl;
	exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
	// switch statement variable
	int c;

	// control variables
	bool die = false;
	int OPT_VERSION = 0;

	// long form arguments
	static struct option long_options[] = {{"threads", required_argument, NULL, 't'}, {"debug", required_argument, NULL, 'b'}, {"version", no_argument, &OPT_VERSION, 1}, {"verbose", no_argument, NULL, 'v'}, {NULL, 0, NULL, 0}};

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "t:v", long_options,
							&option_index)) != -1)
	{
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c)
		{
		case 't':
		{
			stringstream convert(optarg);
			if (!(convert >> opt::threads))
			{
				cerr << "Error - Invalid parameter t: "
					 << optarg << endl;
				return 0;
			}
			break;
		}
		case 'v':
		{
			opt::verbose++;
			break;
		}
		case '?':
		{
			die = true;
			break;
		}
		}
	}

#if defined(_OPENMP)
	if (opt::threads > 0)
		omp_set_num_threads(opt::threads);
#endif

	if (OPT_VERSION)
	{
		// printVersion();
	}

	vector<string> inputFiles;
	while (optind < argc)
	{
		inputFiles.emplace_back(argv[optind]);
		assert(Util::fexists(inputFiles.back()));
		optind++;
	}

	// Check needed options
	if (inputFiles.size() == 0)
	{
		cerr << "Error: Need Input File" << endl;
		die = true;
	}
	if (die)
	{
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	double time = omp_get_wtime();

	CompareCounts comp(inputFiles); // 初始化读文件

	if (inputFiles.size() == 1)
	{
		if (opt::verbose > 1)
		{
			cerr << "error! need more file" << endl;
		}
	}
	else
	{
		comp.computeScore();
	}

	cerr << "Time: " << omp_get_wtime() - time << " s Memory: " << Util::getRSS() << " kbytes" << endl;
	return 0;
}
