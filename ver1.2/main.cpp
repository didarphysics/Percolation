#include "myHeader.h"

using namespace std;

struct bondPair {
	int position1;		//the format of position is x*lsize+y
	int position2;
};

int connectedSiteNumber=0;
int currentRootSite=-1;

vector<bondPair> bondIndex;
vector<int> state;
vector<int> orderOfBond;
vector<double> empty;
vector< vector<double> > Entropy;
vector<double> avgEntropy;
vector<double> canonicalEntropy;

vector<int> stackOfCluster;

void coldStart(){
	bondPair temp;
	bondIndex.clear();
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
	state.clear();
	for(int i=0;i<lsize*lsize;i++){
		state.push_back(-1);
	}
	Entropy.clear();
	for(int i=0;i<numberOfBond;i++){
		Entropy.push_back(empty);
		for(int j=0;j<ensemble_size;j++){
			Entropy[i].push_back(0.0);
		}
	}
	avgEntropy.clear();
	canonicalEntropy.clear();
	for(int i=0;i<numberOfBond;i++){
		avgEntropy.push_back(0.0);
		canonicalEntropy.push_back(0.0);
	}
}

void initialize(){
	connectedSiteNumber=0;
	currentRootSite=-1;
	orderOfBond.clear();
	stackOfCluster.clear();
	/*
	for(int i=0;i<lsize*lsize;i++){
		stackOfCluster.push_back(i);
	}
//*/
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
			destroy(rootSite2);
			state[rootSite1] += state[rootSite2];
			state[rootSite2]=rootSite1;
		}
		else{
			destroy(rootSite1);
			state[rootSite2] += state[rootSite1];
			state[rootSite1]=rootSite2;
			currentRootSite=rootSite2;
		}
	}
	update(currentRootSite);
}

int numberOfOccupiedSite(){
	int result=0;

	for(int i=0;i<lsize*lsize;i++){
		if(state[i]!=-1){
			result++;	
		}
}

return result;
}

void calculateParameters(int step,int bondNumber){
	Entropy[bondNumber][step]=calculateEntropy();
}

int main(){
	for(int fileNumber=0;fileNumber<numberOfFile;fileNumber++){
		coldStart();
		std::string name = "res/collection/entropy_L" + intToString(lsize)+"_E"+intToString(ensemble_size)+"_"+intToString(fileNumber)+".data";
		ofstream	myFile(name.c_str());
		//ofstream	exp("resExp/exponentNu", std::ofstream::out | std::ofstream::app);
		for(int step=0;step<ensemble_size;step++){
			initialize();
			for(int i=0;i<numberOfBond;i++){
				connectBond(orderOfBond[i]);
				calculateParameters(step,i);
			}
			cout<<step<<"_"<<fileNumber<<endl;
		}
		getAvg(Entropy,numberOfBond,avgEntropy);
		/*
		for(int i=0;i<bondIndex.size();i++){
			cout<<i<<"th bond, position1 x= "<<bondIndex[i].position1/lsize<<" and y= "<<bondIndex[i].position1%lsize<<endl;	
			cout<<i<<"th bond, position2 x= "<<bondIndex[i].position2/lsize<<" and y= "<<bondIndex[i].position2%lsize<<endl;
			cout<<"********************************************************"<<endl;
		}
		*/
		for(int i=0;i<canonicalEntropy.size();i++){
			myFile<<(double)(i+1)/(double)(numberOfBond)<<" "<<avgEntropy[i]<<endl;
		}
		myFile.close();
	}
	
	combine();
	return 0;
}
