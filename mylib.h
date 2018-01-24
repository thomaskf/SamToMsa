#ifndef _MYLIB_
#define _MYLIB_

#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <cstring>

#define MAXDIM 100000

using namespace std;

static string i2s[10] = {"0","1","2","3","4","5","6","7","8","9"};

static long double cumlogI[MAXDIM];

static bool initialCumLogI = false;

string int2str(int i);

void tokenizer(string seq, string separators, vector<string>* result);

void randomizePos(vector<int>& arr, int n);

void randomizePos(vector<int>& arr);

void statistics(long double* arr, int n, long double& min, long double& max, long double& avg, long double& median, long double& stdev);

void statistics(vector<long double>& arr, long double& min, long double& max, long double& avg, long double& median, long double& stdev);

// return the value of topK
int getTopK(vector<int>& arr, int topK); 

long combin(int n, int k);

void initializeCumLogI();

long double logCombin(int n, int k);

long double logBinomial(int x, int n, long double p);

// remove all spaces inside the string
void removeSpaces(string& s);

// compute the BIC value
long double BIC(long double logL, int degFree, int dataSize);

// compute the AIC value
long double AIC(long double logL, int degFree, int dataSize);

// compute the Adjusted BIC value
long double ABIC(long double logL, int degFree, int dataSize);

// compute the CAIC value
long double CAIC(long double logL, int degFree, int dataSize);

// compute the AICc value
long double AICc(long double logL, int degFree, int dataSize);

// double to string
string doublToStr(double d, double decimalPlace);

// round the double to a certain number of decimal place
double roundNumber(double d, int decimalPlace);

// compute log(x+y), given log(x) and log(y)
long double log_x_plus_y(long double logx, long double logy);

#endif
