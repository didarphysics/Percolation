#include "myHeader.h"

extern std::vector<int> state;

double checkPercolation(int pos){
	for(int i=0;i<lsize;i++){
		if(findRoot(i)==pos){	//0*lsize+i; x=0;y=i
			for(int j=0;j<lsize;j++){
				if(findRoot((lsize-1)*lsize+j)==pos){		//(lsize-1)*lsize+j;x=lsize-1;y=j
					return 1.0;
				}
			}
			break;
		}
	}
	for(int i=0;i<lsize;i++){
		if(findRoot(i*lsize)==pos){	//i*lsize+0; y=0;x=i
			for(int j=0;j<lsize;j++){
				if(findRoot(j*lsize+(lsize-1))==pos){		//i*lsize+(lsize-1);y=lsize-1;x=i
					return 1.0;
				}
			}
			break;
		}
	}
	return 0.0;
}
