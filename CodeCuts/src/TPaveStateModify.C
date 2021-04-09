/*******************************************/
// Author: Ivan Dario Piernagorda Peña     
// Date:   15/02/2020                  
// Title:  It Modify the TPave in a Graph
/*******************************************/

#ifndef TPAVESTATEMODIFY_C
#define TPAVESTATEMODIFY_C

#include "../include/TPaveStateModify.h"
#include "../include/Libraries.h"

using namespace std;

TPaveStateModify::TPaveStateModify(TH1 *PH,UInt_t L){
  PrincipalHisto=PH;
  Histos.push_back(PH);
  gPad->Update();
  Pave = (TPaveStats*)PrincipalHisto->FindObject("stats");
  Pave->SetName("mystats");
  ListText = Pave->GetListOfLines();

  switch(L){
  case 1:
    Const = Pave->GetLineWith("Entries");
    ListText->Remove(Const);
    break;
  case 2:
    Const = Pave->GetLineWith("Entries");
    ListText->Remove(Const);
    Const = Pave->GetLineWith("Mean");
    ListText->Remove(Const);
    break;
  }
 
  
}


TPaveStateModify::TPaveStateModify(TH1 *PH,TH1 *SH){
  PrincipalHisto=PH;
  Histos.push_back(PH);
  Histos.push_back(SH);
  gPad->Update();
  Pave = (TPaveStats*)PrincipalHisto->FindObject("stats");
  Pave->SetName("mystats");
  ListText = Pave->GetListOfLines();
  Const = Pave->GetLineWith("Entries");
  ListText->Remove(Const);
}


TPaveStateModify::TPaveStateModify(TH1 *PH,vector<TH1*> SH){
  PrincipalHisto=PH;
  Histos.push_back(PH);
  for (UInt_t i = 0; i < SH.size(); i++) Histos.push_back(SH[i]);
  gPad->Update();
  Pave = (TPaveStats*)PrincipalHisto->FindObject("stats");
  Pave->SetName("mystats");
  ListText = Pave->GetListOfLines();
  Const = Pave->GetLineWith("Entries");
  ListText->Remove(Const);
}

void TPaveStateModify::BoxOptStatActive(int it){
  k = Histos.at(it)->GetKurtosis();
  K = Histos.at(it)->GetKurtosis(11);
  s = Histos.at(it)->GetSkewness();
  S = Histos.at(it)->GetSkewness(11);
  i = Histos.at(it)->Integral();
  o = abs(Histos.at(it)->Integral(Histos.at(it)->GetNbinsX(),Histos.at(it)->GetNbinsX()+1));
  u = abs(Histos.at(it)->Integral(0,Histos.at(it)->GetNbinsX())-Histos.at(it)->Integral());
  r = Histos.at(it)->GetRMS();
  R = Histos.at(it)->GetRMSError();;
  m = Histos.at(it)->GetMean();
  M = Histos.at(it)->GetMeanError();
  e = Histos.at(it)->GetEntries();
  n = Histos.at(it)->GetName();
}

void TPaveStateModify::BoxOptStat(string TS,int Precision){

  for(UInt_t it = 0; it < Histos.size(); it++){
    BoxOptStatActive(it);
    string TSTemp = TS;

    while(true){
      ostringstream Obj;
      string StrObj;
      string STR;
      TLatex *myt = NULL;

  
      Obj.clear();
      switch(TSTemp.front()){

      case 'k':
	Obj << std::fixed;
	Obj << std::setprecision(Precision);
	Obj << k;
	StrObj = Obj.str();

	STR = "KRT = " + StrObj;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
	break;

    
      case 'K':

	Obj << std::fixed;
	Obj << std::setprecision(Precision);
	Obj << K;
	StrObj = Obj.str();

	STR = "KRT Error = " + StrObj;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
	break;

      case 's':

	Obj << std::fixed;
	Obj << std::setprecision(Precision);
	Obj << s;
	StrObj = Obj.str();

	STR = "SKN = " + StrObj;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
	break;
    
      case 'S':

	Obj << std::fixed;
	Obj << std::setprecision(Precision);
	Obj << S;
	StrObj = Obj.str();

	STR = "SKN Error = " + StrObj;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
    

      case 'i':
	
	Obj << int(i);
	StrObj = Obj.str();
	
	STR = "Integral = " + StrObj;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
	
	break;
	
      case 'o':
	
	Obj << int(o);
	StrObj = Obj.str();

	STR = "Overflow = " + StrObj;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
	
	break;
	
      case 'u':
	
	Obj << int(u);
	StrObj = Obj.str();

	STR = "Underflow = " + StrObj;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
	
	break;

      case 'r':
    
	Obj << std::fixed;
	Obj << std::setprecision(Precision);
	Obj << r;
	StrObj = Obj.str();

	STR = "RMS = " + StrObj;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
	break;

      case 'R':

	Obj << std::fixed;
	Obj << std::setprecision(Precision);
	Obj << R;
	StrObj = Obj.str();

	STR = "RMS Error = " + StrObj;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
	break;

      case 'm':

	Obj << std::fixed;
	Obj << std::setprecision(Precision);
	Obj << m;
	StrObj = Obj.str();

	STR = "Valor medio = " + StrObj;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
	break;
    
      case 'M':
    
	Obj << std::fixed;
	Obj << std::setprecision(Precision);
	Obj << M;
	StrObj = Obj.str();
    
	STR = "Error medio = " + StrObj;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
	break;

      case 'e':

	Obj << std::fixed;
	Obj << std::setprecision(0);
	Obj << e;
	StrObj = Obj.str();
    
	STR = "Entradas = " + StrObj;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
	break;

      case 'n':
	STR = "Nombre = " + n;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
	break;
      }

      TSTemp.erase(TSTemp.begin());
      if(TSTemp.empty()) break;
    }
  }
}

void TPaveStateModify::BoxSize(double X, double Y){
  Pave->SetX1NDC(X);
  Pave->SetY1NDC(Y);
}

void TPaveStateModify::BoxPosition(double X1, double Y1, double X2, double Y2){
  Pave->SetX1(X1);
  Pave->SetY1(Y1);
  Pave->SetX2(X2);
  Pave->SetY2(Y2);
}

void TPaveStateModify::BoxTextSize(double TS){
  Pave->SetTextSize(TS);
}

void TPaveStateModify::SaveChanges(){
  PrincipalHisto->SetStats(0);
  gPad->Modified();
}

#endif
