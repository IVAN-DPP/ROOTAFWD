#ifndef TPAVESTATEMODIFY_H
#define TPAVESTATEMODIFY_H

#include "Libraries.h"


class TPaveStateModify{
protected:
  TH1 *PrincipalHisto = NULL;
  vector<TH1*> Histos;
  TPaveStats *Pave = NULL;
  TList *ListText  = NULL;
  double XSize;
  double YSize;
  double XStat;
  double YStat;
  double TextSize;

  double k;
  double K;
  double s;
  double S;
  double *i         = NULL;
  double o;
  double u;
  double r;
  double R;
  double m;
  double M;
  double e;
  string n;
  string Modes;
  
public:
  TPaveStateModify(){}
  TPaveStateModify(TH1 *,TH1 *);
  TPaveStateModify(TH1 *,vector<TH1*>);
  void BoxOptStatActive(int);
  void BoxOptStat(string,int = 3);
  void BoxSize(double,double);
  void BoxPosition(double,double,double,double);
  void BoxTextSize(double);
  void SaveChanges();
};


#endif
