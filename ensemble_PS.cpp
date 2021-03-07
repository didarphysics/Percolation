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


#define L 300                  // lattice size 
int N =(L*L);                   // we actually need l*l so writing l*l many time we simply give it N total lattice point
#define EMPTY -(N+1)            // to make all the lattice parent at first
#define recipocal_lattice 1.0/N // for disequilibrium we need this
#define version  3            //never write version 0
#define version_to 3


int ensemble =10000;
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
  double k1,k2,k3;
  for(int fileNumber=0;fileNumber<version;fileNumber++){
    //string name="res/percolation_strength_L_300_E_40000.txt";
    //string name="res/percolation_strength_L_50_E_80000.txt";
    string name= "L_"+conversion_int_to_string(L)+"/per_strength_L_" + conversion_int_to_string(L) +"_E_"+ conversion_int_to_string(ensemble)+"_v_"+conversion_int_to_string(fileNumber)+".txt";
    ifstream  myFile(name.c_str());
    Data_C.push_back(empty);
    
	cout<<"okay "<<fileNumber<< endl;
    do{
      myFile>>k1>>k2>>k3;                //change it before run
    //  myFile>>k1>>k2;
      Data_C[fileNumber].push_back(k2);
     //Data_C[fileNumber].push_back(k2);
      //Data_D[fileNumber].push_back(k3);
    }while(!myFile.eof());
    myFile.close();
  }

  //cout<<avg_data_Cize()<<" "<<numberOfBond<<endl;
 // cout<<"i am okay 1"<<endl;
  for(int i=0;i<N;i++){
    newData_C.push_back(empty);
    
    for(int j=0;j<version;j++){       //change is before run another code
      newData_C[i].push_back(0.0);
      
    }
  }
  //cout<<"hello"<<endl;
  cout<<"i am okay 2"<<endl;
  for(int i=0;i<N;i++){
    for(int j=0;j<version;j++){
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
    avg_C[i] /=(double)version;
    //cout<<" avg_C "<<avg_C[i]<<endl;
    
  }
  cout<<"hey"<<endl;
  //cout<<" "<<newData_S[1].size()<<endl;

  /*
  canonicalTech(newData_S,newData_D,version,avg_data_C,avg_data_C,canonicalData_S,canonicalData_D);
  std::string name1 = "res/entropy_L" + conversion_int_to_string(L)+"_E"+conversion_int_to_string(ensemble*version) +".txt";
  ofstream  myFile1(name1.c_str());
  for(int i=0;i<canonicalData_S.size();i++){
    myFile1.precision(10);
    double complexcity=canonicalData_S[i]*canonicalData_D[i];
    myFile1<<fixed<<p[i]<<"\t\t "<<canonicalData_S[i]<<"\t\t"<<canonicalData_D[i]<<"\t\t"<<complexcity<< endl;
  }
  //*/
}

void Jump_workshop(){
  string name2= "res/per_strength_L" + conversion_int_to_string(L) +"_E_"+ conversion_int_to_string(ensemble*version_to) +".txt";
  ofstream  sp1(name2.c_str());
  

  double difference=0.0;;
  for(int i=0;i<N;i++){
    //cout<<"size"<<avg_S.size()<<endl;
    
    //x_difference = ((i+1)/(double)N ) - (i/(double)N);
    difference = (avg_C[i+1]-avg_C[i]);
    
    //cout<<difference<<endl;
    sp1.precision(10);
    sp1<<fixed<<(i+1)/(double)N<<"\t\t"<<difference<<endl;
   //fprintf(fp, "%lf    %lf      %lf\n",i/(double)N, S[i] ,Slope_entropy );
   }
}

void slope_percolation_strength(){
  string name2= "res/percolation_strength_slope_L_"+conversion_int_to_string(L)+"_E_"+conversion_int_to_string(ensemble)+"_v_"+conversion_int_to_string(version_to)+".txt";
  ofstream  sp1(name2.c_str());
  double x_difference=0.0;  
  double Slope_PS =0.0,difference=0.0;;
  for(int i=0;i<N-1;i++){
    //cout<<"size"<<avg_S.size()<<endl;
    
    x_difference = ((i+1)/(double)N ) - (i/(double)N);
    Slope_PS = (avg_C[i+1]-avg_C[i]);
    //difference    = avg_D[i+1]-avg_D[i];
    //cout<<difference<<endl;
    sp1.precision(10);
    sp1<<fixed<<i/(double)N<<"\t\t"<<Slope_PS/(x_difference)<<endl;
   //fprintf(fp, "%lf    %lf      %lf\n",i/(double)N, S[i] ,Slope_entropy );
   }
}


void canonical_ensemble(){
  std::string name1 = "res/percolation_strength_L_" + conversion_int_to_string(L)+"_E_"+conversion_int_to_string(ensemble*(version+1)) +".txt";
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
  //Jump_workshop();
 //slope_percolation_strength();

  return 0;
}
