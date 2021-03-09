/*******************************************/
// Author: Ivan Dario Piernagorda Peña     
// Date:   22/02/2021                  
// Title:  Divide histogram with same events                              
/*******************************************/

#ifndef HISTOBINNING_C
#define HISTOBINNING_C

#include "../include/Libraries.h"
#include "../include/HistoBinning.h"


HistoBinning::HistoBinning(TH1 *h,UInt_t part){
  PrincipalHisto= h;
  Parts 	= part;
  Nb		= PrincipalHisto->GetNbinsX();
  Xmin		= PrincipalHisto->GetXaxis()->GetXmin();
  Xmax		= PrincipalHisto->GetXaxis()->GetXmax();
  DX		= abs(Double_t(Xmax-Xmin)/Double_t(Nb));
  DXX		= Xmin;
  At		= PrincipalHisto->Integral("width");
}

void HistoBinning::DoHistoBinning(){

  UInt_t ii 	= 0;
  Double_t AtT	= 0;
  
  for(UInt_t j=1;j<Parts;j++)

    for (UInt_t i =ii; i < Nb; i++) {

      Double_t f  =PrincipalHisto->GetBinContent(i);
      Double_t f1 =PrincipalHisto->GetBinContent(i+1);
      A+=f*DX;
      if(i == 0) Point.push_back(Xmin);
      if(A + f1*DX > j*(At/Parts)){ 
	ii=i;
	Area.push_back(abs(A-Ev));
	AreaError.push_back(abs(A-Ev-At/Parts));
	Areat.push_back((At/Parts));
	Point.push_back(DXX);
	Events.push_back(abs(A/DX-Ev/DX));
	Ev=A;
	AtT=At;
	if(j == Parts-1) {
	  Point.push_back(Xmax);
	  Events.push_back(PrincipalHisto->Integral("width")/DX-A/DX);
	}
	break;
      }
      DXX+=DX;
    }
  
}

void HistoBinning::PrintLevel(UInt_t Lev){
  
  switch(Lev){
  case 0:
    for(UInt_t i=0;i<Area.size();i++)
      cout << setprecision(2) 	<< fixed
	   << "Numeric Area: " 	<< Area[i] 	<< "\t"
	   << "Theoric Area: " 	<< Areat[i] 	<< "\t"
	   << "to: " << i 	<< " partition"	<< endl;
    break;
  case 1:
    for(UInt_t i=0;i<Area.size();i++)
      cout << setprecision(2) 	<< fixed
	   << "Numeric Area: " 	<< Area[i] 	<< "\t"
	   << "Theoric Area: " 	<< Areat[i] 	<< "\t"
	   << "Delta Error : "	<< AreaError[i]	<< "\t"
	   << "to: " << i 	<< " partition"	<< endl;
    break;
  case 2:
    for(UInt_t i=0;i<Area.size();i++)
      cout << setprecision(2) 	<< fixed
	   << "Numeric Area: " 	<< Area[i] 	<< "\t"
	   << "Theoric Area: " 	<< Areat[i] 	<< "\t"
	   << "Delta Error : "	<< AreaError[i]	<< "\t"
	   << "Point : "	<< Point[i+1]	<< "\t"
	   << "to: " << i 	<< " partition"	<< endl;
    break;
  default :
    cout << "\nThis is not an option\n";
    break;
  }
}

void HistoBinning::GetLatexTable(string Path,string Caption,string Label){

  cout << setprecision(2);
  fstream LTXT;
  LTXT.open(Path,ios::out);
  
  LTXT << "\\begin{table}[H]" << endl
       << "\t\\centering"       << endl
       << "\t\\begin{tabular}{|c|c|c|}" << endl
       << "\t\\hline" << endl
       << "\t\tPartition & Angular Range & Events \\\\ \\hline\n";

  for (UInt_t i = 0; i < Events.size(); i++) {
    LTXT << "\t\t" << i+1 << "\t&"
	 << "$" <<Point[i] << "\\leq\\cos{\\theta^{cm}_{K^+}}\\leq" << Point[i+1] << "$" << "\t&"
	 << "$" <<Events[i] << "$"<< "\t\\\\ \\hline\n";
  }


  LTXT << "\t\\end{tabular}" << endl
       << "\t\\caption{"     << Caption << "}" << endl
       << "\t\\label{tab:"  << Label   << "}" << endl
       << "\\end{table}"     << endl;
  
  cout << "We Create: " << Path.substr(Path.find_last_of("/")+1,Path.size()-Path.find_last_of("/")) << endl;
  LTXT.close();

}

void HistoBinning::GetPoints(vector<double>& AvValue,vector<double>& Error){
  double Err = 0;
  double AvV = 0;
  for (UInt_t i = 0; i < Point.size()-1; i++) {
    if(Point[i] <= 0){
      Err = (Point[i]-Point[i+1])/2;
      AvV = Err+Point[i+1];
    }
    else {
      Err = (Point[i+1]-Point[i])/2;
      AvV = Err+Point[i];
    }

    AvValue.push_back(AvV);
    Error.push_back(abs(Err));

  }

}

vector<double> HistoBinning::GetPartitions(){
  return Point;
}

#endif
