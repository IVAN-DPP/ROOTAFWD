#include "Libraries.h"

#ifndef HISTOGRAMS_C
#define HISTOGRAMS_C

using namespace std;

class Histograms{

 protected:

  TH2F *h_DeltaBe[3];
  TH2F *h_BeVSp[3];
  TH1F *h_DeltaTall[2];
  TH2F *h_DeltaTallvsp[2];
  TH1F *h_DeltaT[1];
  TH1F *h_eloss[3];

  //---- Missing mass ----//

  TH1F *h_MissingMass;
  TH1F *h_MissingMass_kaonpion;

  TH2F *h_MissingMass_vsMissingMasskaonpion;


  TH1F *h_MissingP;
  TH2F *h_MissingPvsMass;
  TH2F *h_MissingMassvsSigmaMass;
  TH2F *h_MissingPvsSigmaMass;
 
  TH1F *h_InvariantMass;
  TH1F *h_LambdaMass;

 public:
  Histograms(){}
  void DoHistograms();
  void DoCanvas();
};


void Histograms::DoHistograms(){

  // ----------------------------//
  
  h_DeltaBe[0]=new TH2F("h_DeltaBe_0","Proton ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  h_DeltaBe[1]=new TH2F("h_DeltaBe_1","Kaon ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.15, 0.15);
  h_DeltaBe[2]=new TH2F("h_DeltaBe_2","Pion ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);

  //--------------------------- // 

  h_BeVSp[0]=new TH2F("h_BeVSp_0","Proton ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp[1]=new TH2F("h_BeVSp_1","Kaon ;p [GeV/c];#beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp[2]=new TH2F("h_BeVSp_2","Pion ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  //----------------------------//
  
  h_DeltaTall[0]=new TH1F("h_DeltaTall_0","Kaon ;#Delta t [ns];counts;", 200, -10, 10);
  h_DeltaTall[1]=new TH1F("h_DeltaTall_1","Pion ;#Delta t [ns];counts;", 200, -10, 10);


  h_DeltaTallvsp[0]=new TH2F("h_DeltaTallvsp_0","Kaon ;#Delta t [ns];p [GeV/c];", 200, 0, 3, 200, -10, 10);
  h_DeltaTallvsp[1]=new TH2F("h_DeltaTallvsp_1","Pion ;#Delta t [ns];p [GeV/c];", 200, 0, 3, 200, -10, 10);
  h_DeltaT[0]=new TH1F("h_DeltaT_0","Kaon ;#Delta t [ns];counts;", 200, -10, 10);

  //-----------------------------//

  h_eloss[0]=new TH1F("h_eloss_0","Proton; eloss",50, 0, 100);
  h_eloss[1]=new TH1F("h_eloss_1","Kaon; eloss",50, 0, 20);
  h_eloss[2]=new TH1F("h_eloss_2","Pion; eloss",50, 0, 20);

  // ------ Missing mass ------- //

  h_MissingMass = new TH1F("h_missingmass","Missing mass;", 100, 0.7, 1.2);
  h_MissingMass_kaonpion = new TH1F("h_missingmass_kaonpion","Missing mass Kaon_pion;", 100, 0.7, 1.2);
  h_MissingMass_vsMissingMasskaonpion = new TH2F("MissingMass_correlation", "",100, 0.7, 1.2, 100, 0.7, 1.2);
  h_MissingP = new TH1F("h_missingp","Missing momentum;", 100, 0.0, 1.5);
  h_MissingPvsMass = new TH2F("h_missingpvsm","Missing;", 100, 0.0, 1.5, 100, 0.7, 1.2);
  h_MissingMassvsSigmaMass = new TH2F("h_missingmvsSigma","Missing;",100,1.0,1.5, 100, 0.7, 1.2);
  h_MissingPvsSigmaMass = new TH2F("h_missingPvsSigma","Missing;",100,1.0,1.5, 100, 0.0, 1.5);
  h_InvariantMass = new TH1F("h_InvariantMass","Invariant mass;", 100, 1.0, 1.5);
  h_LambdaMass = new TH1F("h_LambdaMass","Invariant mass;", 100, 1.08, 1.16);
}


void Histograms::DoCanvas(){
  
  TCanvas *c0=new TCanvas("c0","Delta Beta y Beta", 900, 500);
  c0->Divide(2,3);
  c0->cd(1);
  h_DeltaBe[0]->Draw("colz");
  c0->cd(3);
  h_DeltaBe[1]->Draw("colz");
  c0->cd(5);
  h_DeltaBe[2]->Draw("colz");

  c0->cd(2);
  h_BeVSp[0]->Draw("colz");
  c0->cd(4);
  h_BeVSp[1]->Draw("colz");
  c0->cd(6);
  h_BeVSp[2]->Draw("colz");
  

  TCanvas *c1=new TCanvas("c1","Delta T", 900, 500);
  c1->Divide(2,3);
  c1->cd(1);
  h_DeltaTall[0]->Draw(); 
  c1->cd(2);
  h_DeltaTall[1]->Draw();
  c1->cd(3);
  h_DeltaTallvsp[0]->Draw("colz"); 
  c1->cd(4);
  h_DeltaTallvsp[1]->Draw("colz");
  c1->cd(5);
  h_DeltaT[0]->Draw();

  TCanvas *c2=new TCanvas("c2","Energy loss", 900, 500);
  c2->Divide(3,1);
  c2->cd(1);
  h_eloss[0]->Draw();
  c2->cd(2);
  h_eloss[1]->Draw();
  c2->cd(3);
  h_eloss[2]->Draw();

  TCanvas *c3=new TCanvas("c3","Missing mass", 900, 500);
  c3->Divide(2,3);
  c3->cd(1);
  h_MissingMass->Draw();
  c3->cd(2);
  h_MissingP->Draw();
  c3->cd(3);
  h_InvariantMass->Draw();
  c3->cd(4);
  h_LambdaMass->Draw();
  c3->cd(5);
  h_MissingMass_kaonpion->Draw();
  c3->cd(6);
  h_MissingMass_vsMissingMasskaonpion->Draw();

}

#endif
