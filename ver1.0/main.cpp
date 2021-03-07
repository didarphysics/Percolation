
#include "myHeader.h"

using namespace std;

struct bondPair {

int position1;		//the format of position is x*lsize+y
int position2;

};

int alreadyPercolated=0;

int connectedSiteNumber=0;
int currentRootSite=-1;

vector<bondPair> bondIndex;
vector<int> state;
vector<int> orderOfBond;
vector<int> empty;
vector< vector<int> > percolationState;
vector<double> avgSP;
vector<double> canonicalSP;

vector<double> p_c;
double avgp_c=0.0;
double errp_c=0.0;

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
	for(int i=0;i<numberOfBond;i++){
		percolationState.push_back(empty);
		for(int j=0;j<ensemble_size;j++){
			percolationState[i].push_back(0);
		}
	}
	for(int i=0;i<numberOfBond;i++){
		avgSP.push_back(0.0);
		canonicalSP.push_back(0.0);
	}
	for(int i=0;i<ensemble_size;i++){
		p_c.push_back(0.0);
	}
}

void initialize(){
	alreadyPercolated=0;
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

void avgTech(){
	for(int i=0;i<p_c.size();i++){
		avgp_c +=p_c[i];
	}
	avgp_c /=(double)p_c.size();

	for(int i=0;i<p_c.size();i++){
		errp_c +=(p_c[i]-avgp_c)*(p_c[i]-avgp_c);
	}
	errp_c /=(double)p_c.size();
	errp_c = sqrt(errp_c);
}

void calculateParameters(int step,int bondNumber){
	if(alreadyPercolated==0){
		if(checkPercolation(currentRootSite)==1.0){
			alreadyPercolated=1;
			p_c[step]=(double)bondNumber/(double)(numberOfBond);
		}
	}
	if(alreadyPercolated==1){
		percolationState[bondNumber][step]=1.0;
	}
}

int main(){
	coldStart();
	//ofstream	sp("res/spaningProbability_L100_E1000.data");
	ofstream	exp("resExp/exponentNu", std::ofstream::out | std::ofstream::app);

	std::string name = "res/spanningProbability_L" + intToString(lsize)+"_E"+intToString(ensemble_size) + ".data";
	ofstream	sp(name.c_str());


	for(int step=0;step<ensemble_size;step++){
		initialize();
		for(int i=0;i<numberOfBond;i++){
			connectBond(orderOfBond[i]);
			calculateParameters(step,i+1);
		}
		cout<<step<<endl;
	}
	canonicalTech(percolationState,numberOfBond,avgSP,canonicalSP);
	avgTech();
	/*
	for(int i=0;i<bondIndex.size();i++){
		cout<<i<<"th bond, position1 x= "<<bondIndex[i].position1/lsize<<" and y= "<<bondIndex[i].position1%lsize<<endl;
		cout<<i<<"th bond, position2 x= "<<bondIndex[i].position2/lsize<<" and y= "<<bondIndex[i].position2%lsize<<endl;
		cout<<"********************************************************"<<endl;
	}
	*/
	for(int i=0;i<canonicalSP.size();i++){
		sp<<(double)(i+1)/(double)(numberOfBond)<<" "<<canonicalSP[i]<<endl;
	}
	exp<<lsize<<" "<<errp_c<<" "<<avgp_c<<endl;
	sp.close();
	exp.close();
	return 0;
}
