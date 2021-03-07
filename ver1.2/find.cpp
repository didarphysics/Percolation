#include"myHeader.h"

extern std::vector<int> state;

int findRoot(int pos){
	if(state[pos]<0){
		return pos;
	}
	else{
		return state[pos]=findRoot(state[pos]);
	}
}
