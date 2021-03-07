#include "myHeader.h"

void getAvg(std::vector< std::vector<double> > &sample,int sampleSize, std::vector<double> &avgSample){
	for(int i=0;i<sampleSize;i++){
		for(int j=0;j<ensemble_size;j++){
			avgSample[i] +=sample[i][j];
		}
		avgSample[i] /=(double)(ensemble_size);
	}
	
}
