#include "myHeader.h"

using namespace std;

string intToString(int t){
	std::string ch;
	ostringstream outs;
	outs << t; // Convert value into a string.
	ch = outs.str();
return ch;
} 
