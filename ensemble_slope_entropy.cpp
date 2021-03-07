#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <time.h>
#include "stdlib.h"
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <sstream>
#include <limits>
#include <numeric>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
using namespace std;


#define L 250                   // lattice size 
int N =(L*L);                   // we actually need l*l so writing l*l many time we simply give it N total lattice point
#define EMPTY (N+1)            // to make all the lattice parent at first
#define recipocal_lattice 1.0/N // for disequilibrium we need this
#define version 14
extern int errno;               // for creating folder
int ensemble =500;
vector <double> avg_S;
vector <double> avg_D;
vector<double> p;

string conversion_int_to_string(int integar){
   ostringstream ch;
    ch << integar;
    return ch.str();  
 }

void combine(){
  ///*
  vector< vector<double> > Data_S;
  vector< vector<double> > Data_D;
  vector< vector<double> > newData_S;
  vector< vector<double> > newData_D;
  
  vector<double> avgData_S;
  vector<double> avgData_D;
  vector<double> canonicalData_S;
  vector<double> canonicalData_D;
  vector<double> empty;
  //*/
  double k1,k2,k3,k4,k;
  for(int fileNumber=0;fileNumber<=version;fileNumber++){
    string name = "L_"+conversion_int_to_string(L)+"/entropy_disequilibrium_L" + conversion_int_to_string(L) +"_E_"+ conversion_int_to_string(ensemble) +"_V_" + conversion_int_to_string(fileNumber)+".txt";
    ifstream  myFile(name.c_str());
    Data_S.push_back(empty);
    Data_D.push_back(empty);

    do{
      myFile>>k1>>k2>>k3>>k4;
      Data_S[fileNumber].push_back(k2);
      Data_D[fileNumber].push_back(k3);
    }while(!myFile.eof());
    myFile.close();
  }

  //cout<<avgData.size()<<" "<<numberOfBond<<endl;
 // cout<<"i am okay 1"<<endl;
  for(int i=0;i<N;i++){
    newData_S.push_back(empty);
    newData_D.push_back(empty);
    for(int j=0;j<version;j++){
      newData_S[i].push_back(0.0);
      newData_D[i].push_back(0.0);
    }
  }
  //cout<<"hello"<<endl;
 // cout<<"i am okay 2"<<endl;
  for(int i=0;i<N;i++){
    for(int j=0;j<version;j++){
      newData_S[i][j]=Data_S[j][i];
      newData_D[i][j]=Data_D[j][i];
      //std::cout<<i<<" "<<j<<endl;
    }
  }

 // cout<<"i am okay 3"<<endl;
  for(int i=0;i<N;i++){
    avg_S.push_back(0.0);
    avg_D.push_back(0.0);
    //canonicalData.push_back(0.0);
  }
  for(int i=0;i<newData_S.size();i++){
    for(int j=0;j<newData_S[i].size();j++){
     // cout<<"i am okay 4"<<endl;
      avg_S[i] +=newData_S[i][j];
      avg_D[i] +=newData_D[i][j];
    }
    //cout<<"i am okay 4"<<endl;
    avg_S[i] /=(double)(version);
    avg_D[i] /=(double)(version);
  }
  //cout<<"hey"<<endl;
  //cout<<" "<<newData_S[1].size()<<endl;
  //for(int i=0;i<N;i++)
  //cout<<"  "<<(i+1)/(double)N<<"   "<<avg_S[i]<<"   "<<avg_D[i] <<endl;
}


void canonical_ensemble(){
  std::string name1 = "res/slope_entropy_L_" + conversion_int_to_string(L)+"_E"+conversion_int_to_string(ensemble*version) +".txt";
  ofstream  myFile1(name1.c_str());
   
  double binom[(L*L)+1];
    for (int j = 1; j <=avg_S.size(); j++)
    {  
      // cout<<"i am okay!"<<endl;
        int n, i;
        double prob = j*1.0/(L*L);
        binom[j] = 1;
       
        for (n = j+1; n <= L*L; ++n)
            binom[n] = binom[n-1]*(L*L-n+1)*1.0/n*prob/(1-prob);
       // printf("binomila_1 %lf \n",binom[N] );
        
        for (n = j - 1; n >= 0; --n)
            binom[n] = binom[n+1]*(n+1)*1.0/(L*L-n)*(1-prob)/prob;
        
        //printf("binomila_2 %lf \n",binom[N] );
        
        double sum = 0;
        for (i = 0; i <= L*L; ++i) sum += binom[i];
        for (i = 0; i <= L*L; ++i) binom[i] /= sum;          // upto this is same for all
        
        double sum_entropy = 0.0, sum_D=0.0;
        for (n = 1; n <= L*L; ++n) sum_entropy += avg_S[n]*binom[n];
        //for (n = 1; n <= L*L; ++n) sum_D += avg_D[n]*binom[n];
        //double complexcity = sum_D*sum_entropy;
        myFile1.precision(10);
      if(j>(N*0.2))
        myFile1<<fixed<<prob<<"\t\t"<<sum_entropy<<endl;
        //fprintf(fp1,"%lf   %lf    %lf    %lf    %lf    %lf\n",prob,complexcity, S[j],sum_entropy,D[j],sum_D);
    }
}


int main(){


  combine();
  canonical_ensemble();
  return 0;
}
