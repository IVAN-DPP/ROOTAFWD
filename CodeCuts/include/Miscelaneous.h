//Author: Edwin Munevar

/* Miscelaneous.h
   This code contains the variable and function declarations for the Miscelaneous.C routine.*/

#ifndef Miscelaneous_h
#define Miscelaneous_h

#include "Libraries.h"
#include "../src/Miscelaneous.C"

using namespace std;

// DirPathName without the final /
// for example, ./ -> DirPathName = "." or ./Tables/ -> DirPathName = "./Tables"
// Extention without point, for example file.root -> Extention = "root"

void ListFilesAtDir(string,vector<string>&);

void ListFilesAtDir(string,string,vector<string>&);

int LoadPolTable(int, char *,map<vector<float>,int>&);

double GetPol(int, double, int, int, double, double);

double GetPol(int, double, double, int, double, double);

void GetPolAv(vector<float>,vector<vector<int>>&,vector<vector<double>>&,double);

void GetPolAvTable();
#endif
