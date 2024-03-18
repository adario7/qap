#include <ilcplex/cplex.h>

// seems faster on bigger problems
extern bool PARAM_DYNAMIC_SEARCH;
// greatly improves the LP
extern bool PARAM_LOCAL_L;
// makes nodes enumeration slower, slightly improves the LP
extern bool PARAM_LOCAL_L_PAIRS;
// makes nodes enumeration slower, slightly improves the LP
extern bool PARAM_LOCAL_L_ALL;
// makes nodes enumeration slighlty slower, improves the LP
extern bool PARAM_LOCAL_M;
// whether cuts are calculated only once per node, even after a new relaxation
extern bool PARAM_CUT_ONCE;
// don't bother adding cuts if less then a minimum are found
extern int PARAM_CUTS_MIN;
// -1 to disable
extern double PARAM_TIME_LIMIT;
// single vs multi thread
extern bool PARAM_SINGLE_THREAD;
// deterministic vs opportunistic
extern bool PARAM_OPPORTUNISTIC;
// memory limit
extern int PARAM_MEMLIMIT;
// display nodes as they are explored
extern bool PARAM_LIVE_SOL;

void read_parameters(int argc, char** argv);

void apply_parameters(CPXENVptr env, CPXLPptr lp);
