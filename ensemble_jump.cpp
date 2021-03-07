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
using namespace std;


#define L 300                   // lattice size
int N =(L*L);                   // we actually need l*l so writing l*l many time we simply give it N total lattice point
#define EMPTY (N+1)            // to make all the lattice parent at first
#define recipocal_lattice 1.0/N // for disequilibrium we need this
#define version 1 //never write version 0

int ensemble =449000;
vector <double> avg_C;
vector<double> p;

string conversion_int_to_string(int integar){
   ostringstream ch;
    ch << integar;
    return ch.str();
 }

void combine(){
  ///*
  vector< vector<double> > Data_C;

  vector< vector<double> > newData_C;


  vector<double> avg_data_C;
  vector<double> canonicalData_C;

  vector<double> empty;
  //*/
  double k1,k2;
  for(int fileNumber=0;fileNumber<=version;fileNumber++){
    //string name = "L_"+conversion_int_to_string(L)+"/largest_cluster_L" + conversion_int_to_string(L) +"_E_"+ conversion_int_to_string(ensemble) +"_V_" + conversion_int_to_string(fileNumber) +".txt";
    //string name = "res/largest_cluster_L150_E_100000_1.txt";
    string name ="res/largest_cluster_L300_E_449000_1.txt";
    ifstream  myFile(name.c_str());
    Data_C.push_back(empty);


    do{
      myFile>>k1>>k2;
      Data_C[fileNumber].push_back(k2);
      //Data_D[fileNumber].push_back(k3);
    }while(!myFile.eof());
    myFile.close();
  }

  //cout<<avg_data_Cize()<<" "<<numberOfBond<<endl;
 // cout<<"i am okay 1"<<endl;
  for(int i=0;i<N;i++){
    newData_C.push_back(empty);

    for(int j=0;j<1;j++){       //change is before run another code
      newData_C[i].push_back(0.0);

    }
  }
  //cout<<"hello"<<endl;
  for(int i=0;i<N;i++){
    for(int j=0;j<1;j++){
      newData_C[i][j]=Data_C[j][i];

      //std::cout<<i<<" "<<j<<endl;
    }
  }

 // cout<<"i am okay 3"<<endl;
  for(int i=0;i<N;i++){
    avg_C.push_back(0.0);

    //canonicalData.push_back(0.0);
  }
  for(int i=0;i<newData_C.size();i++){
    //cout<<"size "<<newData_C[i].size()<<endl;
    for(int j=0;j<newData_C[i].size();j++){
     // cout<<"i am okay 4"<<endl;
      avg_C[i] +=newData_C[i][j];
      //cout<<"newData_C "<<newData_C[i][j]<<endl;

    }
   // cout<<"i am okay 4"<<endl;
    //cout<<" avg_C "<<avg_C[i]<<endl;

  }
}


void canonical_ensemble(){
  std::string name1 = "res/order_parameter_L_" + conversion_int_to_string(L)+"_E"+conversion_int_to_string(ensemble*(version+1)) +".txt";
  ofstream  myFile1(name1.c_str());

  double binom[(L*L)+1];
    for (int j = 1; j <=avg_C.size(); j++)
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

        double sum_order_parameter=0.0;
        for (n = 0; n < L*L; ++n) {
          sum_order_parameter+= avg_C[n]*binom[n];
          //cout<<binom[n]<<endl;

        }
        //for (n = 1; n <= L*L; ++n) sum_D += avg_D[n]*binom[n];

        myFile1.precision(10);
        myFile1<<fixed<<prob<<"\t\t"<<sum_order_parameter<<"\t\t"<<endl;
        //fprintf(fp1,"%lf   %lf    %lf    %lf    %lf    %lf\n",prob,complexcity, S[j],sum_entropy,D[j],sum_D);
    }
}


int main(){
  combine();
  canonical_ensemble();
  return 0;
}
