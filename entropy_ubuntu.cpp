// this is the main code for simulating entropy of L*L latice using ziff-newman algorithm
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <time.h>
#include "stdlib.h"
#include <stdio.h>
using namespace std;

#define L 200                   // lattice size
int N =(L*L);                   // we actually need l*l so writing l*l many time we simply give it N total lattice point
#define EMPTY (-N-1)            // to make all the lattice parent at first
#define recipocal_lattice 1.0/N // for disequilibrium we need this
FILE *fp;

int ensemble =1;
 vector <vector <int> > nn;               // Nearest neighbors
 vector <int> empty;                      // for making 2d vector for neighbors
 vector<int>  parent;                     // Array of pointers
 vector <int> Occupationu_order;                       //In which Occupationu_order lattice point will occupy
 vector < vector <double> > entropy;
 vector <vector < double> > disequilibrium;
 vector <double> empty1;                    // for initializing 2d vector for ensemble

 extern int errno;

void initialize(){
    parent.clear();
    Occupationu_order.clear();
    nn.clear();

    for (int l=0; l<N; l++) parent.push_back(EMPTY);  // all latice point are initialize by some unique value

    for(int l=0;l<N;l++){
      nn.push_back(empty);
      for(int a=0;a<4;a++){
        nn[l].push_back(0);
      }
    }

 }

void boundaries(){            // boundary condition

    int i;

    for (i=0; i<N; i++) {
        nn[i][0] = (i+1)%N;   // first portion of entry is neighbor to the right
        nn[i][1] = (i+N-1)%N;  // second portion of entry in nn is the neighbor to the left
        nn[i][2] = (i+L)%N;    //third portion of entry in nn is the neighbor below
        nn[i][3] = (i+N-L)%N;  //fourth portion of entry in nn is the neighbor above

        if (i%L==0) nn[i][1] = i+L-1;    // maps left side to right side
        if ((i+1)%L==0) nn[i][0] = i-L+1; // Maps right side to left side
    }
 }

 // We must generate a random Occupationu_order for which the sites will be occupied
void permutation(){
      int i,j;
      int temp;
      srand(time(NULL));
      double r;

      for (i=0; i<N; i++) Occupationu_order.push_back(i);  // first they are all giving their correspoding lattice number
      // here we will decide the occupation Occupationu_order of the lattice
        //random_shuffle(Occupationu_order.begin(),Occupationu_order.end());


      for (i=0; i<N; i++) {                  // here we will decide the occupation Occupationu_order of the lattice
          r = ((double) rand() / (RAND_MAX));
          j = i + (N-i)*r;
          temp = Occupationu_order[i];
          Occupationu_order[i] = Occupationu_order[j];
          Occupationu_order[j] = temp;
       }


}


 /* We are going to need to have a method of finding the root of a cluster ... */
 /* a Recussive method is chosen to do this */

int findroot(int i){           // here the root of the cluster is founded or we search the parent of the cluster

    if (parent[i]<0) return i;   // roots are given - sign so if it is - than root is founded of that cluster
    return parent[i] = findroot(parent[i]); // we make the respective child to parent to make code fast, this is called path compresion
 }

void entropy_disequilibrium(int en,int l){
    double S=0.0,D=0.0,w=0.0;
    int total_clusterSize=0;


    for(int size=0;size<N;size++){

       if(parent[size] < 0){

          w=(-1)*parent[size]/(double)(l+1);        // probability of taking a cluster of size parent[k]
          // both are - so w becomes +
          //if(w!=0.0) {     // to avoid nan value cz if w is zero entropy cannt be measured
            S +=(-1)*w*log10(w);
            D += (w - recipocal_lattice)*(w - recipocal_lattice);
           // }
        }
    }
       // cout<<" "<<i/(double)N<<" "<<(slope+1)<<endl;
    entropy[l][en] = S;
    disequilibrium[l][en] = D;
}

void percolate()
 {
  for(int en=0;en<ensemble;en++){

    if(en%1==0){ cout<<"ensemble "<<en+1<<endl; }

    initialize();
    boundaries();
    permutation();

    int l,n;
    int cell1,cell2;
    int root1,root2;
    int big=0;            // to find the biggest cluster

    //for (i=0; i<N; i++) parent.push_back(EMPTY);  // all latice point are initialize by some unique value

    for (l=0; l<N; l++) {
      root1 = cell1 = Occupationu_order[l];
      parent[cell1] = -1;           // -1 means occupy

      for (n=0; n<4; n++) {
        cell2 = nn[cell1][n];

        if (parent[cell2]!=EMPTY) {    // if occupy than check the neighbour
          root2 = findroot(cell2);  // check the root of neighbour

         if (root2!=root1) {       // if the dont represent same root than have to marge

          if (parent[root1]>parent[root2]) {  // marge the smaller cluster to the bigger cluster
            parent[root2] += parent[root1];   // remeber roots are given - value.so -2 > -4
            parent[root1] = root2;  //ptr r1 then becomes value of its root , r2
            root1 = root2;          // still un clear about its jobs
          }
          else {                         //   But if ptr r1 contains more elements that ptr the complete the opposite action
            parent[root1] += parent[root2];
            parent[root2] = root1;
           }
             // if (-parent[root1]>big) big = -parent[root1];
          }
       }
       }
        entropy_disequilibrium(en,l);
    }
   }
 }

string conversion_int_to_string(int integar){
   ostringstream ch;
    ch << integar;
    return ch.str();
 }

void avg(int v){
    string name = "L_"+conversion_int_to_string(L)+"/entropy_disequilibrium_L" + conversion_int_to_string(L) +"_E_"+ conversion_int_to_string(ensemble) +"_V_" + conversion_int_to_string(v) +".txt";
    ofstream  sp(name.c_str());
    double C = 0.0;

    for (int l=0;l<N;l++) {
     double  sum_ensemble_entropy=0.0;
     double  sum_ensemble_disequilibrium =0.0;

     for(int en=0;en<ensemble;en++) {
       // sum_ensemble_entropy += entropy[l][en];
        sum_ensemble_entropy=accumulate(entropy[l].begin(), entropy[l].end(), 0.0);  //using STL to sum
        sum_ensemble_disequilibrium=accumulate(disequilibrium[l].begin(),disequilibrium[l].end(),0.0);
        //sum_ensemble_disequilibrium += disequilibrium[l][en];
       }

     C =sum_ensemble_entropy*sum_ensemble_disequilibrium/(ensemble*ensemble);
     sp.precision(10);      // it will print 10 number
     //fprintf(fp, "%lf     %lf     %lf     %lf\n",(double)(l+1)/N,sum_ensemble_entropy/ensemble,sum_ensemble_disequilibrium/ensemble, C);
     sp<<fixed<<(double)(l+1)/N<<"\t\t"<<sum_ensemble_entropy/ensemble<<"\t\t"<<sum_ensemble_disequilibrium/ensemble<<"\t\t"<<C<<endl;
    }
    //fclose(fp);
  sp.close();
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
  double k1,k2,k3,k4;
  for(int fileNumber=0;fileNumber<=version_to;fileNumber++){
    string name = "L_"+conversion_int_to_string(L)+"/entropy_disequilibrium_L" + conversion_int_to_string(L) +"_E_"+ conversion_int_to_string(ensemble) +"_V_" + conversion_int_to_string(fileNumber)+".txt";
    ifstream  myFile(name.c_str());
    Data_S.push_back(empty);
    Data_D.push_back(empty);

    do{
      myFile>>k1>>k2>>k3>>k4;
      Data_S[fileNumber].push_back(k2);
      Data_D[fileNumber].push_back(k3);
    }
    while(!myFile.eof());
    myFile.close();
  }

  for(int a=0;a<Data_S.size();a++){
    for(int k=0;k<Data_S[a].size();k++){
     }
  }
  for(int i=0;i<N;i++){
    newData_S.push_back(empty);
    newData_D.push_back(empty);
    for(int j=0;j<version_to;j++){
      newData_S[i].push_back(0.0);
      newData_D[i].push_back(0.0);
    }
  }

  for(int i=0;i<N;i++){
    for(int j=0;j<version_to;j++){
      newData_S[i][j]=Data_S[j][i];
      newData_D[i][j]=Data_D[j][i];
    }
  }

  for(int i=0;i<N;i++){
    avg_S.push_back(0.0);
    avg_D.push_back(0.0);

  }

  for(int i=0;i<newData_S.size();i++){
    for(int j=0;j<newData_S[i].size();j++){
      avg_S[i] +=newData_S[i][j];
      avg_D[i] +=newData_D[i][j];
    }
    avg_S[i] /=(double)(version_to);
    avg_D[i] /=(double)(version_to);
  }
}


void canonical_ensemble(){
  std::string name1 = "res/entropy_L_" + conversion_int_to_string(L)+"_E"+conversion_int_to_string(ensemble*(version_to+1)) +".txt";
  ofstream  myFile1(name1.c_str());

  double binom[(L*L)+1];
    for (int j = 1; j <=avg_S.size(); j++)
    {

        int n, i;
        double prob = j*1.0/(L*L);
        binom[j] = 1;

        for (n = j+1; n <= L*L; ++n)
            binom[n] = binom[n-1]*(L*L-n+1)*1.0/n*prob/(1-prob);


        for (n = j - 1; n >= 0; --n)
            binom[n] = binom[n+1]*(n+1)*1.0/(L*L-n)*(1-prob)/prob;

        double sum = 0;
        for (i = 0; i <= L*L; ++i) sum += binom[i];
        for (i = 0; i <= L*L; ++i) binom[i] /= sum;          // upto this is same for all

        double sum_entropy = 0.0, sum_D=0.0;
        for (n = 1; n <= L*L; ++n) sum_entropy += avg_S[n]*binom[n];
        for (n = 1; n <= L*L; ++n) sum_D += avg_D[n]*binom[n];
        double complexcity = sum_D*sum_entropy;
        myFile1.precision(10);
        myFile1<<fixed<<prob<<"\t\t"<<sum_entropy<<"\t\t"<<sum_D<<"\t\t"<<complexcity<<endl;
    }
}


int main(){

  int e1,e2;
  struct stat sb;
  string nf = "L_"+conversion_int_to_string (L);
  char const *name_of_folder1 = nf.c_str();
  char const *name_of_folder2 = "res";

  e1 = stat(name_of_folder1, &sb);
  e2 = stat(name_of_folder2, &sb);
  //printf("e=%d errno=%d\n",e,errno);
    if(e1!=0)
    e1 = mkdir(name_of_folder1, S_IRWXU);
    if(e2!=0)
    e1 = mkdir(name_of_folder2, S_IRWXU);

   for(int v=version_from ;v<=version_to;v++){
      cout<<"version "<<v<<endl;
      entropy.clear();
      disequilibrium.clear();
      for (int in=0;in<N;in++){

         entropy.push_back(empty1);
         disequilibrium.push_back(empty1);

         for(int en=0;en<ensemble;en++) {

           entropy[in].push_back(0.0);
           disequilibrium[in].push_back(0.0);
          }
        }
    percolate();
    avg(v);
    }
  combine();
  canonical_ensemble();
  return 0;
}
