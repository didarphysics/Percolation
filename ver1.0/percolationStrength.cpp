#include"myHeader.h"

extern std::vector<int> state;
std::vector<int> stackOfPercolatedCluster;

void push(int pos){
	bool flag=true;
	for(int i=0;i<stackOfPercolatedCluster.size;i++){
		if(stackOfPercolatedCluster[i]==pos){
			flag=false;
		}
	}
	if(flag){
		stackOfPercolatedCluster.push_back(pos);
	}
}

void update(){
	for(int i=0;i<lsize;i++){
		if(state[i]>0){	//0*lsize+i; x=0;y=i
			for(int j=0;j<lsize;j++){
				if(state[(lsize-1)*lsize+j]==state[i]){		//(lsize-1)*lsize+i;x=lsize-1;y=i
					push(state[i]);
				}
			}
		}
	}

	for(int i=0;i<lsize;i++){
		if(state[i*lsize]>0){	//i*lsize+0; y=0;x=i
			for(int j=0;j<lsize;j++){
				if(state[i*lsize+(lsize-1)]==state[i*lsize]){		//i*lsize+(lsize-1);y=lsize-1;x=i
					push(state[i*lsize]);
				}
			}
		}
	}
	return 0.0;
}

double calculatePS(){
	update();
	double result=0.0;
	for(int i=0;i<stackOfPercolatedCluster.size;i++){
		result -=state[pos];
	}
	result /= (double)(lsize*lsize);
	return result;
}
