#include <ilcplex/cplex.h>

// seems faster on bigger problems
extern bool PARAM_DYNAMIC_SEARCH;
// local L cuts
extern bool PARAM_LOCAL_L;
// F cuts between fixed-free pairs
extern bool PARAM_LOCAL_L_PAIRS;
// F cuts between free variables and all fixed ones
extern bool PARAM_LOCAL_L_ALL;
// cuts for a better big M at the local node
extern bool PARAM_LOCAL_M;
// cuts between free pairs using local Ls
extern bool PARAM_LOCAL_FREE;
// when true, only appy cuts once per node
extern bool PARAM_CUT_ONCE;
// when != 0, don't add cuts that are violated by less then a fraction of sqrt(obj)
extern double PARAM_REL_DELTA;
// calculate expensive cuts only when there are no more simple cuts
extern bool PARAM_EXP_LATER;
// don't bother adding cuts if less then a minimum are found
extern int PARAM_CUTS_MIN;
// -1 = disable all cut types, 0 = usa a cut factor of zero, 1 = default cplex cuts
extern int PARAM_CPLEX_CUTS;
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
// the random seed for the run
extern int PARAM_SEED;

void read_parameters(int argc, char** argv);

void apply_parameters(CPXENVptr env, CPXLPptr lp);
