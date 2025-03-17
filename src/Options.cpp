

#include <string>
#include <vector>
// #include "src/Options.h"

using namespace std;

namespace opt {
extern int verbose;
extern unsigned threads;
extern unsigned k;
extern string summary;
extern double covThresh;

extern unsigned minCov;
extern bool dupes;
// extern uint64_t minSites;

extern string ref;
extern unsigned window;
extern unsigned multi;
extern double scoreThresh;

extern string debug;
extern std::vector<std::string> snp;
} // namespace opt
