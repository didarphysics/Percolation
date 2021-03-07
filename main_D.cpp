//code for finding D

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <time.h>
#include "stdlib.h"
#include <stdio.h>

FILE *fp;

using namespace std;
int lsize=10;
int ensemble=10;
int currentLabel=0;
double recipocal_lattice=1.0/(lsize*lsize);
double S=0.0;
double D=0.0;

#define up(x) ((x+1)%lsize)            // boundary condition
#define down(x) ((x+lsize-1)%lsize)

vector< vector<int> > label;
vector<int> empty1;
vector<int> cluster_size;//(lsize*lsize+1,0);
vector <int> stack;
vector <float> empty2;
vector < vector<float> > entropy_ensemble;
vector < vector <float> > D_ensemble;


void initialize(){
	currentLabel=0;      // for each ensemble the initial condition
	cluster_size.clear();
	cluster_size.resize(lsize*lsize+1,0);
	//cluster_size[0] = lsize*lsize;
	stack.clear();
	label.clear();

	for(int i=0;i<lsize;i++){    // initialize the 2d vector
		label.push_back(empty1);
		for(int j=0;j<lsize;j++){
			label[i].push_back(-1);
		}
	}

}

void size_workshop (int x,int y){

    if(label[up(x)][y]!=-1 && label[up(x)][y]!=currentLabel){
		stack.push_back(label[up(x)][y]);
	}
	if(label[down(x)][y]!=-1 && label[down(x)][y]!=currentLabel){
		stack.push_back(label[down(x)][y]);
	}
	if(label[x][up(y)]!=-1 && label[x][up(y)]!=currentLabel){
		stack.push_back(label[x][up(y)]);
	}
	if(label[x][down(y)]!=-1 && label[x][down(y)]!=currentLabel){
		stack.push_back(label[(x)][down(y)]);
	}

	if (stack.size()>0){
		cluster_size[currentLabel]=1;
		//cout<<" okay 1"<<endl;
		for(int i=0;i<stack.size();i++){
				cluster_size[currentLabel]+= cluster_size[stack[i]];
		//		cout<<"cluster size a  "<<cluster_size[currentLabel]<<endl;
				cluster_size[stack[i]]=0;
		}
	}
	else{
		cluster_size[currentLabel]	=1;
	//	cout<<"cluster size b  "<<cluster_size[currentLabel]<<endl;
	}

 }

void clustering(int x,int y){
	label[x][y]=currentLabel;
	if(label[up(x)][y]!=-1 && label[up(x)][y]!=currentLabel){
		clustering(up(x),y);
	}
	if(label[down(x)][y]!=-1 && label[down(x)][y]!=currentLabel){
		clustering(down(x),y);
	}
	if(label[x][up(y)]!=-1 && label[x][up(y)]!=currentLabel){
		clustering(x,up(y));
	}
	if(label[x][down(y)]!=-1 && label[x][down(y)]!=currentLabel){
		clustering(x,down(y));
	}
}

void entropy_workshop(){
	register int deposition=0;
	S=0.0;            // before each depostion  entropy must ber zero
	D=0.0;
	register double probability=0.0;

	//cout<<" recipocal_lattice"<<recipocal_lattice<<end;

	for (int c=1;c<=lsize*lsize;c++){  // culster_suze starts from 1
		deposition+= cluster_size[c];   //it will give us the total deposition at this certain time
	}
	for (int c=1;c<=deposition;c++){
		probability= cluster_size[c]/(double)deposition;  //it has some specific name, which i forget at this moment

		if(probability!=0.0){                        //to make it faster
    		S+= probability*log10(probability);
    		D+=(probability - recipocal_lattice)*(probability - recipocal_lattice);
    		//cout<<" D " <<D<<endl;
    	}
	}
	//entropy_ensemble[en].push_back(S);
//	cout<<"deposition  "<<(double)i/(lsize*lsize)<<"  entropy "<<S <<endl;
//	fprintf(fp,"%lf		%lf\n",(double)i/(lsize*lsize),S);
	//cout<<" entropy "<<S<<endl;
   // return D;
}

int main(){
	srand(time(NULL));
	register int r,x,y;
//	fp = fopen ("D_entropy_100_500.txt","w");
	register double sum_D=0.0;
	register double sum_entropy=0.0;
//	cout<<"okay 1"<<endl;
	///*
	for (register int l=0;l<=lsize*lsize;l++){
		entropy_ensemble.push_back(empty2);
		for(register int en=0;en<ensemble;en++){
			entropy_ensemble[l].push_back(0.0);
		}
	}
	//*/
	for (register int l=0;l<=lsize*lsize;l++){
		D_ensemble.push_back(empty2);
		for(register int en=0;en<ensemble;en++){
			D_ensemble[l].push_back(0.0);
		}
	}

	for (register int en=0;en<ensemble;en++){
		cout<<"ensemble "<<en+1<<endl;
		initialize();

		for(int i=0;i<lsize*lsize;i++){
			r=rand()%(lsize*lsize);        // generating a random number between 0 to lsize^2
			x=r%lsize;
			y=r/lsize;
			if(label[x][y]==-1){
				++currentLabel;

			//	cout<<"I am okay1"<<endl;
				size_workshop(x,y);
				clustering(x,y);
				entropy_workshop();
				D_ensemble[i][en]=D;
				entropy_ensemble[i][en]=S;
			//	cout<<"okay 2"<<endl;
				//cout<<"entropy "<<entropy_ensemble[en]<<endl;
			//	cout<<"I am okay2"<<endl;
			}
			else{
				i--;
			}


		}

		cout<<"okay 2"<<endl;
		for (int d=0;d<D_ensemble.size();d++){
			sum_D=0.0;
			sum_entropy=0.0;
			for(register int en=0;en<D_ensemble[d].size();en++){
				sum_D+=D_ensemble[d][en];
				//cout<<" "<<entropy_ensemble[l][en];
			}

	//cout<<"sum  "<<sum_entropy/4<<endl;
			cout<<"d "<<sum_D/D_ensemble[d].size()<<endl;
		//	fprintf(fp,"%lf		%lf\n",(double)l/(lsize*lsize),sum_D/D_ensemble[l].size());
		}
	fclose(fp);
	}
}
