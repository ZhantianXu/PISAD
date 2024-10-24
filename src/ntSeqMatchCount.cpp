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

#include "FingerPrint.hpp"

#if _OPENMP
#include <omp.h>
#endif

using namespace std;

#define PROGRAM "ntsmCount"

void printHelpDialog()
{
	const char dialog[] =
		"Usage: " PROGRAM " -s [FASTA] [OPTION]... [FILES...]\n"
		"  -t, --threads = INT    Number of threads to run.[1]\n"
		"  -m, --maxCov = INT     k-mer coverage threshold for early termination. [inf]\n"
		"  -d, --dupes            Allow shared k-mers between sites to be counted.\n"
		"  -i, --information      extra debug information.\n"
		"  -s, --snp = STR        variant sketch (one or more) [required]\n"
		"  -k, --kmer = INT       k-mer size used. [19]\n"
		"  -h, --help             Display this dialog.\n"
		"  -v, --verbose          Display verbose output.\n"
		"  -n, --name             name.\n"
		"      --version          Print version information.\n";

	cerr << dialog << endl;
	exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
	// switch statement variable
	int cov;

	// control variables
	bool die = false;
	int OPT_VERSION = 0;

	// long form arguments
	static struct option long_options[] = {
		{"threads", required_argument, NULL, 't'},
		{"maxCov", required_argument, NULL, 'm'},
		{"dupes", required_argument, NULL, 'd'},
		{"snp", required_argument, NULL, 's'},
		{"kmer", required_argument, NULL, 'k'},
		{"help", no_argument, NULL, 'h'},
		{"version", no_argument, &OPT_VERSION, 1},
		{"information", no_argument, NULL, 'i'},
		{"name", required_argument, NULL, 'n'},
		{NULL, 0, NULL, 0}};

	int option_index = 0;
	while ((cov = getopt_long(argc, argv, "s:t:n:ivhk:m:d", long_options,
							  &option_index)) != -1)
	{
		istringstream arg(optarg != NULL ? optarg : "");
		switch (cov)
		{
		case 'h':
		{
			printHelpDialog();
			break;
		}
		case 'd':
		{
			opt::dupes = true;
			break;
		}
		case 'i':
		{
			opt::information = true;
			break;
		}
		case 'n':
		{
			stringstream convert(optarg);
			if (!(convert >> opt::name))
			{
				cerr << "Error - Invalid parameter n: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 's':
		{
			stringstream ss(optarg);
			string file;
			while (getline(ss, file, ','))
			{
				if (!file.empty())
				{
					opt::snp.push_back(file); // 将每个解析出来的文件路径添加到 opt::snp
				}
			}
			break;
		}
		case 'm':
		{
			stringstream convert(optarg);
			if (!(convert >> opt::covThresh))
			{
				cerr << "Error - Invalid parameter m: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 'k':
		{
			stringstream convert(optarg);
			if (!(convert >> opt::k))
			{
				cerr << "Error - Invalid parameter k: " << optarg << endl;
				return 0;
			}
			break;
		}
		case 't':
		{
			stringstream convert(optarg);
			if (!(convert >> opt::threads))
			{
				cerr << "Error - Invalid parameter t: " << optarg << endl;
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

	if (opt::k > 32)
	{
		die = true;
		cerr << "k cannot be greater than 32" << endl;
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

	FingerPrint fp;
	cerr << "read fa" << endl;
	cerr << "Time: " << omp_get_wtime() - time << " s Memory: " << (Util::getRSS() / (1024 * 1024)) << " G" << endl;
	fp.computeCounts(inputFiles);
	// fp.computeCountsProducerConsumer(inputFiles);
	// fp.printOptionalHeader();
	// fp.printCountsMax();
	// fp.printCountsAllCounts();
	cerr << "begin print" << endl;
	cerr << "Time: " << omp_get_wtime() - time << " s Memory: " << (Util::getRSS() / (1024 * 1024)) << " G" << endl;
	fp.printCountsMode();
	cerr << fp.printInfoSummary() << endl;
	cerr << "Time: " << omp_get_wtime() - time << " s Memory: " << (Util::getRSS() / (1024 * 1024)) << " G" << endl;
	return 0;
}
