#ifndef TPAVESTATEMODIFY_C
#define TPAVESTATEMODIFY_C

#include "../include/Libraries.h"
#include "../include/TPaveStateModify.h"

TPaveStateModify::TPaveStateModify(TH1 *PH,TH1 *SH){
  PrincipalHisto=PH;
  Histos.push_back(SH);
  gPad->Update();
  Pave = (TPaveStats*)PrincipalHisto->FindObject("stats");
  Pave->SetName("mystats");
  ListText = Pave->GetListOfLines();
}


TPaveStateModify::TPaveStateModify(TH1 *PH,vector<TH1*> SH){
  PrincipalHisto=PH;
  Histos = SH;
  gPad->Update();
  Pave = (TPaveStats*)PrincipalHisto->FindObject("stats");
  Pave->SetName("mystats");
  ListText = Pave->GetListOfLines();
}

void TPaveStateModify::BoxOptStatActive(int it){
  k = Histos.at(it)->GetKurtosis();
  K = Histos.at(it)->GetKurtosis(11);
  s = Histos.at(it)->GetSkewness();
  S = Histos.at(it)->GetSkewness(11);
  i = Histos.at(it)->GetIntegral();
  o;
  u;
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
    
      case 'i': break;
      case 'o': break;
      case 'u': break;
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

	STR = "Mean = " + StrObj;
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
    
	STR = "Mean Error = " + StrObj;
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
    
	STR = "Entries = " + StrObj;
	myt = new TLatex(0,0,STR.c_str());
	myt ->SetTextFont(42);
	myt ->SetTextSize(TextSize);
	myt ->SetTextColor(Histos.at(it)->GetFillColor());
	ListText->Add(myt);
	break;

      case 'n':
	STR = "Name = " + n;
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
