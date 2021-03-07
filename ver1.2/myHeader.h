#include<stdio.h>
#include<iostream>
#include<cmath>
#include<time.h>
#include<cstdlib>
#include<stdlib.h>
#include<vector>
#include<algorithm>
#include<fstream>
#include<string>
#include<sstream>

#define lsize 150
#define numberOfBond 2*lsize*(lsize-1)
#define ensemble_size 100
#define numberOfFile 1

void canonicalTech(std::vector< std::vector<double> > &sample,int sampleSize,int sampleSize2, std::vector<double> &avgSample,std::vector<double> &canonicalSample);
double calculateEntropy();
int findRoot(int pos);
void update(int pos);
double checkPercolation(int pos);
void destroy(int pos);
std::string intToString(int t);
void combine();
void getAvg(std::vector< std::vector<double> > &sample,int sampleSize, std::vector<double> &avgSample);
