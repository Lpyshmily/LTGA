#ifndef _PSO_H_
#define _PSO_H_
#include<stdlib.h>
#include<time.h>
#include<iostream>
#include<iomanip>
#include<math.h>
#include<fstream>
#include<process.h>
using namespace std;

void PSO(double (*ObjFun)(const double* X, const double* para), double* xbest, double& fbest, const double* para,
		 int D, int Np, int ItMax=500, int ItOut=100, double OmegaMin=0.4, double OmegaMax=0.9, double C1Min=0.5,
		 double C1Max=2.5, double C2Min=0.5, double C2Max=2.5, double Vmax=0.8);
#endif