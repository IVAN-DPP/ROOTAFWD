//Author: Edwin Munevar

/* Miscelaneous.h
   This code contains the variable and function declarations for the Miscelaneous.C routine.*/

#ifndef Miscelaneous_h
#define Miscelaneous_h

#include "Libraries.h"
#include "Miscelaneous.C"

using namespace std;

// DirPathName without the final /
// for example, ./ -> DirPathName = "." or ./Tables/ -> DirPathName = "./Tables"
// Extention without point, for example file.root -> Extention = "root"

void ListFilesAtDir(string,vector<string>&);

void ListFilesAtDir(string,string,vector<string>&);

//int LoadPolTable(int plane, char *PolTableList);

//double GetPol(int plane, double edge, int eid, int poltype, double lowThresh, double highThresh);

int LoadPolTable(int, char *);

double GetPol(int, double, int, int, double, double);

double GetPol(int, double, double, int, double, double);

bool badSCpad(int, int, int, int);

bool Cut_dcustBadTaggerCounters(int, int, int);

TVector3 DocaCalculation(int, int, double&);

//void Good_photon(int&,double,int,double,double,const TLorentzVector &,vector<double> &,vector<int> &);

void Good_photon(int&,double,double,double,const TLorentzVector &,vector<double> &,vector<int> &);

TLorentzVector define4Vector(double,double,double,double,double);

TLorentzVector eloss_func(const TLorentzVector,TVector3,float);

TLorentzVector lab_to_cm(TLorentzVector, TLorentzVector, TLorentzVector);

//void Setup_dcustTripFile_PartialName(string, char *);
void Setup_dcustTripFile_PartialName(string locPartialTripFileName, char *opt, int &Flag);
void Setup_dcustTripFile(string locTripFileName, int &TripFlag);
bool Cut_dcustTripEvents(int locEventNumber);

//bool badSCpadpos(int, int);
//bool badSCpadneg(int, int);

bool Cut_dcustJunkRuns(string);
bool Cut_dcustJunkFiles(string);

//void GFluxFunction(int, char *, double [], double [], ofstream &);
//string Construct_dcustGFluxFileName(string locPartialGFluxFileName);


#endif
