#include "myHeader.h"

void canonicalTech(std::vector< std::vector<double> > &sample,int sampleSize,int sampleSize2, std::vector<double> &avgSample,std::vector<double> &canonicalSample){
	for(int i=0;i<sampleSize;i++){
		for(int j=0;j<sampleSize2;j++){
			avgSample[i] +=sample[i][j];
		}
		avgSample[i] /=(double)(sampleSize2);
	}
	///*
	std::vector<double> bionom;
	
	bionom.push_back(1);
	for(int i=0;i<=sampleSize;i++){
		bionom.push_back(0.0);
	}
	
	for(int j=1;j<=sampleSize;j++){
		int n;
		double prob=j*1.0/sampleSize;
		bionom[j]=1.0;
		
		for(n=j+1;n<=sampleSize;++n){
			bionom[n] = bionom[n-1]*(sampleSize-n+1)*1.0/n*prob/(1-prob);
		}
		for(n=j-1;n>=0;--n){
			bionom[n] =bionom[n+1]*(n+1)*1.0/(sampleSize-n)*(1-prob)/prob;
		}
		double sum=0.0;
		for(int i=0;i<=sampleSize;++i){
			sum +=bionom[i];
		}
		for(int i=0;i<=sampleSize;++i){
			bionom[i] /=sum;
		}
		
		sum=0;
		
		for(n=1;n<=sampleSize;++n){
			sum +=avgSample[n-1]*bionom[n];
		}
		canonicalSample[j-1]=sum;
	}
	//*/

}
