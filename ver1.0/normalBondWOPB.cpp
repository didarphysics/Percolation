#include<stdio.h>
#include<iostream>
#include<cmath>
#include<time.h>
#include<cstdlib>
#include<stdlib.h>
#include<vector>
#include<algorithm>

#include "myHeader.h"

using namespace std;

struct bondPair {

int position1;		//the format of position is x*lsize+y
int position2;

};

int ensemble_size=1;
int connectedSiteNumber=0;
int currentRootSite=-1;


vector<bondPair> bondIndex;
vector<int> state;
vector<int> orderOfBond;
vector<double> empty;
vector< vector<double> > percolationState;
vector<double> avgSP;
vector<double> canonicalSp;

void coldStart(){
	bondPair temp;
	for(int i=0;i<lsize;i++){
		for(int j=0;j<lsize-1;j++){
			temp.position1=i*lsize+j;
			temp.position2=i*lsize+j+1;
			bondIndex.push_back(temp);		
		}	
	}
	for(int i=0;i<lsize-1;i++){
		for(int j=0;j<lsize;j++){
			temp.position1=i*lsize+j;
			temp.position2=(i+1)*lsize+j;
			bondIndex.push_back(temp);		
		}	
	}
	for(int i=0;i<lsize*lsize;i++){
		state.push_back(-1);
	}
	for(int i=0;i<ensemble_size;i++){
		percolationState.push_back(empty);
		for(int j=0;j<numberOfBond;j++){
			percolationState[i].push_back(0.0);
		}
	}
	for(int i=0;i<numberOfBond;i++){
		avgSP.push_back(0.0);
		canonicalSP.push_back(0.0);
	}
}

void initialize(){
	connectedSiteNumber=0;
	currentRootSite=-1;
	orderOfBond.clear();
	for(int i=0;i<numberOfBond;i++){
			orderOfBond.push_back(i);
	}
	random_shuffle(orderOfBond.begin(),orderOfBond.end());
	for(int i=0;i<lsize*lsize;i++){
		state[i]=-1;
	}
}

int findRoot(int pos){
	if(state[pos]<0){
		return pos;
	}
	else{
		return state[pos]=findRoot(state[pos]);
	}
}

void connectBond(int bI){
	int pos1=bondIndex[bI].position1;
	int pos2=bondIndex[bI].position2;
	int rootSite1;
	int rootSite2;

	if(state[pos1]==-1){
		connectedSiteNumber +=1;
	}
	if(state[pos2]==-1){
		connectedSiteNumber +=1;
	}

	rootSite1=findRoot(pos1);
	rootSite2=findRoot(pos2);
	currentRootSite=rootSite1;
	if(rootSite1!=rootSite2){
		if(state[rootSite1]<state[rootSite2]){
			state[rootSite1] += state[rootSite2];
			state[rootSite2]=rootSite1;
		}
		else{
			state[rootSite2] += state[rootSite1];
			state[rootSite1]=rootSite2;
			currentRootSite=rootSite2;
		}
	}
	


}

void calculateParameters(int step,int bondNumber){
	percolationState[step][bondNumber]=checkPercolation();
}

int main(){
	coldStart();

	for(int step=0;step<ensemble_size;step++){
		initialize();
		for(int i=0;i<numberOfBond;i++){
			connectBond(orderOfBond[i]);
			calculateParameters(step,i);
		}
	}
	canonicalTech(percolationState,numberOfBond,avgSP,canonicalSP);
	/*
	for(int i=0;i<bondIndex.size();i++){
		cout<<i<<"th bond, position1 x= "<<bondIndex[i].position1/lsize<<" and y= "<<bondIndex[i].position1%lsize<<endl;	
		cout<<i<<"th bond, position2 x= "<<bondIndex[i].position2/lsize<<" and y= "<<bondIndex[i].position2%lsize<<endl;
		cout<<"********************************************************"<<endl;
	}
	*/
	
	return 0;
}
