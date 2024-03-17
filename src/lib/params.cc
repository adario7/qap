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
// whether cuts are calculated only once per node, even after a new relaxation
bool PARAM_CUT_ONCE = false;
// don't bother adding cuts if less then a minimum are found
int PARAM_CUTS_MIN = 4;
// -1 to disable
double PARAM_TIME_LIMIT = -1;
// single vs multi thread
bool PARAM_SINGLE_THREAD = false;
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
        } else if (arg == "st") {
            PARAM_SINGLE_THREAD = stoi(next_arg()) != 0;
        } else {
			cout << "unknown parameter: " << arg << endl;
			abort();
		}
	}
}
