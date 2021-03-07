#include "myHeader.h"

using namespace std;

extern std::vector<int> state;
///*
extern std::vector<int> stackOfCluster;

void push(int pos){
	bool flag=true;
	for(int i=0;i<stackOfCluster.size();i++){
		if(stackOfCluster[i]==pos){
			flag=false;
		}
	}
	if(flag){
		stackOfCluster.push_back(pos);
	}
}

void destroy(int pos){
	int index=-1;
	for(int i=0;i<stackOfCluster.size();i++){
		if(stackOfCluster[i]==pos){
			index=i;
		}
	}
	if(index>=0){
		stackOfCluster.erase(stackOfCluster.begin()+index);
	}
}

void update(int pos){
//	if(checkPercolation(pos)==1.0){
	push(pos);
//	}
}

//*/

double calculateEntropy(){
	double result=0.0;
	double p=0.0;
	int totalSite=lsize*lsize;
	for(int i=0;i<stackOfCluster.size();i++){
		p=-1.0*state[stackOfCluster[i]]/(double)totalSite;
		result -= p*log(p);
	}
	return result;
}
