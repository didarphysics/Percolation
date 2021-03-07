#include "myHeader.h"

using namespace std;

void combine(){
	vector< vector<double> > Data;
	vector< vector<double> > newData;
	vector<double> p;
	vector<double> avgData;
	vector<double> canonicalData;
	vector<double> empty;
	double k;
	for(int fileNumber=0;fileNumber<numberOfFile;fileNumber++){
		std::string name = "res/collection/entropy_L" + intToString(lsize)+"_E"+intToString(ensemble_size) +"_"+intToString(fileNumber)+ ".data";
		ifstream	myFile(name.c_str());
		Data.push_back(empty);
		
		do{
			myFile>>k;
			if(fileNumber==0){
				p.push_back(k);
				avgData.push_back(0.0);		
			}
			myFile>>k;
			Data[fileNumber].push_back(k);
		}while(!myFile.eof());
		myFile.close();
	}
	//cout<<avgData.size()<<"	"<<numberOfBond<<endl;
	for(int i=0;i<numberOfBond;i++){
		newData.push_back(empty);
		for(int j=0;j<numberOfFile;j++){
			newData[i].push_back(0.0);
		}
	}
	//cout<<"hello"<<endl;
	for(int i=0;i<numberOfBond;i++){
		for(int j=0;j<numberOfFile;j++){
			newData[i][j]=Data[j][i];
			//std::cout<<i<<" "<<j<<endl;
		}
	}

	//cout<<"hey"<<endl;
	
/*
	for(int i=0; i<avgData.size();i++){
		for(int j=0;j<numberOfFile;j++){
			avgData[i] += Data[j][i];
		}
		avgData[i] /=(double) numberOfFile;
	}
*/
	for(int i=0;i<numberOfBond;i++){
		avgData.push_back(0.0);
		canonicalData.push_back(0.0);
	}
	
	canonicalTech(newData,numberOfBond,numberOfFile,avgData,canonicalData);
	std::string name = "res/entropy_L" + intToString(lsize)+"_E"+intToString(ensemble_size*numberOfFile) +".data";
	ofstream	myFile(name.c_str());
	for(int i=0;i<canonicalData.size();i++){
		myFile<<p[i]<<" "<<canonicalData[i]<<endl;
	}
	
}
