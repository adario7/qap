#include <iostream>
#include <string>

#include <params.hh>

using namespace std;

// seems faster on bigger problems
bool PARAM_DYNAMIC_SEARCH = true;
// greatly improves the LP
bool PARAM_LOCAL_L = true;
// makes nodes enumeration slower, slightly improves the LP
bool PARAM_LOCAL_L_PAIRS = false;
// makes nodes enumeration slower, slightly improves the LP
bool PARAM_LOCAL_L_ALL = false;
// makes nodes enumeration slighlty slower, improves the LP
bool PARAM_LOCAL_M = true;
// whether cuts are calculated only once per node: 0 = never, 1 = only expesive cuts, 2 = always
int PARAM_CUT_ONCE = 0;
// don't bother adding cuts if less then a minimum are found
int PARAM_CUTS_MIN = 4;
// -1 = disable all cut types, 0 = usa a cut factor of zero, 1 = default cplex cuts
int PARAM_CPLEX_CUTS = 0;
// -1 to disable
double PARAM_TIME_LIMIT = 12 * 3600;
// single vs multi thread
bool PARAM_SINGLE_THREAD = false;
// deterministic vs opportunistic
bool PARAM_OPPORTUNISTIC = false;
// memory limit
int PARAM_MEMLIMIT = 14*1024 - 256;
// display nodes as they are explored
bool PARAM_LIVE_SOL = false;

void read_parameters(int argc, char** argv) {
	int argi = 1;

	auto next_arg = [&]() {
		if (argi == argc) {
			cout << "missing parameter value" << endl;
			abort();
		}
		return string(argv[argi++]);
	};

	while (argi < argc) {
		string arg = next_arg();
		if (arg == "ds") {
			PARAM_DYNAMIC_SEARCH = stoi(next_arg()) != 0;
		} else if (arg == "l") {
			PARAM_LOCAL_L = stoi(next_arg()) != 0;
		} else if (arg == "p") {
			PARAM_LOCAL_L_PAIRS = stoi(next_arg()) != 0;
		} else if (arg == "a") {
			PARAM_LOCAL_L_ALL = stoi(next_arg()) != 0;
		} else if (arg == "m") {
			PARAM_LOCAL_M = stoi(next_arg()) != 0;
		} else if (arg == "co") {
			PARAM_CUT_ONCE = stoi(next_arg()) != 0;
		} else if (arg == "live") {
			PARAM_LIVE_SOL = stoi(next_arg()) != 0;
		} else if (arg == "cm") {
			PARAM_CUTS_MIN = stoi(next_arg());
		} else if (arg == "tl") {
			PARAM_TIME_LIMIT = stod(next_arg());
		} else if (arg == "opp") {
			PARAM_OPPORTUNISTIC = stoi(next_arg()) != 0;
		} else if (arg == "ml") {
			PARAM_MEMLIMIT = stoi(next_arg());
		} else if (arg == "st") {
			PARAM_SINGLE_THREAD = stoi(next_arg()) != 0;
		} else if (arg == "cc") {
			PARAM_CPLEX_CUTS = stoi(next_arg());
		} else {
			cout << "unknown parameter: " << arg << endl;
			abort();
		}
	}
}

#define _c(what) if (int _error = what) { \
	cout << "CPX error: " #what << endl; cout << "CPX error: " << _error << endl; abort(); }

void apply_parameters(CPXENVptr env, CPXLPptr lp) {
	_c(CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON)); // enable solver logging
	_c(CPXsetintparam(env, CPXPARAM_Preprocessing_Symmetry, 5));
	if (PARAM_TIME_LIMIT != -1)
		_c(CPXsetdblparam(env, CPXPARAM_TimeLimit, PARAM_TIME_LIMIT));
	if (!PARAM_DYNAMIC_SEARCH)
		_c(CPXsetintparam(env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL));
	if (PARAM_SINGLE_THREAD) {
		_c(CPXsetintparam(env, CPXPARAM_Threads, 1));
	} else if (PARAM_OPPORTUNISTIC) {
		_c(CPXsetintparam(env, CPXPARAM_Parallel, CPX_PARALLEL_OPPORTUNISTIC));
	}
	if (PARAM_MEMLIMIT > 0) {
		_c(CPXsetintparam(env, CPXPARAM_MIP_Strategy_File, 2));
		_c(CPXsetintparam(env, CPXPARAM_Emphasis_Memory, 1));
		_c(CPXsetdblparam(env, CPXPARAM_WorkMem, PARAM_MEMLIMIT));
	}
	if (PARAM_CPLEX_CUTS == -1) {
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_BQP, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_Cliques, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_Covers, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_Disjunctive, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_FlowCovers, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_PathCut, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_Gomory, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_GUBCovers, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_Implied, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_LocalImplied, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_LiftProj, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_MIRCut, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_MCFCut, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_RLT, -1));
		_c(CPXsetintparam(env, CPXPARAM_MIP_Cuts_ZeroHalfCut, -1));
	}
	if (PARAM_CPLEX_CUTS == 0) {
		_c(CPXsetdblparam(env, CPXPARAM_MIP_Limits_CutsFactor, 0));
	}
}
