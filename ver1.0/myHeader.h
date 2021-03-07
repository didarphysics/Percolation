#include <stdio.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>

#define lsize 50
#define numberOfBond 2*lsize*(lsize-1)
#define ensemble_size 2

void canonicalTech(std::vector< std::vector<int> > &sample,int sampleSize, std::vector<double> &avgSample,std::vector<double> &canonicalSample);
double checkPercolation(int pos);
int findRoot(int pos);
std::string intToString(int t);
