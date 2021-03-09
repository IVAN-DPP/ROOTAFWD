/*******************************************/
// Author: Ivan Dario Piernagorda Peña     
// Date:   22/02/2021                  
// Title:  Divide histogram with same events                              
/*******************************************/

#ifndef HISTOBINNING_H
#define HISTOBINNING_H

#include "Libraries.h"

using namespace std;

class HistoBinning{
private:

  TH1 *PrincipalHisto	= NULL;
  UInt_t Parts 		= 0;
  UInt_t Nb		= 0;		// #Binnins
  Double_t Xmin		= 0;		// Xmin of Histo
  Double_t Xmax		= 0;		// Xmax of Histo
  Double_t DX		= 0;		// Xmax-Xmin of Histo begining int Xmin
  Double_t A		= 0;		// Numeric Area
  Double_t At		= 0;		// Theoric Area
  Double_t Ev		= 0;		// Events times DX
  Double_t DXX		= 0;		// Position of X value in histogram

  vector<Double_t> Point;
  vector<Double_t> Area;
  vector<Double_t> Areat;
  vector<Double_t> AreaError;
  vector<Double_t> Events;
  
public:
  HistoBinning(){}
  HistoBinning(TH1 * = NULL,UInt_t = 0);
  void DoHistoBinning();
  void PrintLevel(UInt_t);
  void GetLatexTable(string,string,string);
  void GetPoints(vector<double> &,vector<double> &);
  vector<double> GetPartitions();
};


#endif
