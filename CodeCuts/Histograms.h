#include "Libraries.h"

#ifndef HISTOGRAMS_C
#define HISTOGRAMS_C

using namespace std;

class Histograms{

 protected:
  
  TH2F *h_DeltaBe[3];
  TH2F *h_BeVSp[3];
  TH2F *h_DeltaBecut[3];
  TH2F *h_BeVSpcut[3];

  TF1 *FFits[3];               //Fits for do cuts
  TF1 *FFitsminus[3];
  string FFname[3];
  
  TH1F *h_DeltaTall[2];
  TH2F *h_DeltaTallvsp[2];
  TH1F *h_DeltaT[2];
  TH1F *h_eloss[3];

  //--> Beta vs P with PDG MASSES  
  TF1 *BeVSpProton;
  TF1 *BeVSpKaon;
  TF1 *BeVSpPion;
  
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

  TH2F *h_DeltaBVSInvariantMass;
  TH2F *h_DeltaBVSMissingMass;
  TH2F *h_DeltaBVSMissingMomentum;

  //--- Ellipse Cuts ---- //
  
  TEllipse *myEllipse;

   //-----Correlation Theta-Phi, ----------//
 
  TH2F *h_ThePhi_proton;
  TH2F *h_ThePhi_kaon;
  TH2F *h_ThePhi_pion;

   //--- Fiduciary cuts ---//
  TH2F *h_ThePhi_protoncut;
  TH2F *h_ThePhi_kaoncut;
  TH2F *h_ThePhi_pioncut;
  
  
 public:
  Histograms(){}
  void DoHistograms();
  void DoCanvas();
};


void Histograms::DoHistograms(){

  //------------------ Delta Beta ---------------//
  
  h_DeltaBe[0]=new TH2F("h_DeltaBe_0","Proton ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  h_DeltaBe[1]=new TH2F("h_DeltaBe_1","Kaon ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.15, 0.15);
  h_DeltaBe[2]=new TH2F("h_DeltaBe_2","Pion ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);

  
  h_BeVSp[0]=new TH2F("h_BeVSp_0","Proton ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp[1]=new TH2F("h_BeVSp_1","Kaon ;p [GeV/c];#beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp[2]=new TH2F("h_BeVSp_2","Pion ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);

  //--> Beta vs P with PDG MASSES

  BeVSpProton = new TF1("BeVSpProton","x/std::sqrt(std::pow(x,2)+std::pow(0.938,2))",0,2);
  BeVSpProton->SetLineColor(kBlack);
  BeVSpProton->SetLineStyle(2);
  
  BeVSpKaon = new TF1("BeVSpKaon","x/std::sqrt(std::pow(x,2)+std::pow(0.493,2))",0,2);
  BeVSpKaon->SetLineColor(kBlack);
  BeVSpKaon->SetLineStyle(2);

  
  BeVSpPion = new TF1("BeVSpPion","x/std::sqrt(std::pow(x,2)+std::pow(0.139,2))",0,2);
  BeVSpPion->SetLineColor(kBlack);
  BeVSpPion->SetLineStyle(2);
  
  //------- Fits for do cuts in Delta B ----------- //

  FFits[0] = new TF1("DBProtonFit","gaus",0,3);
  FFits[1] = new TF1("DBKaonFit","0.025", 0, 3);
  FFits[2] = new TF1("DBPionFit","0.05",0,3);


  
  
  //-------- Delta Beta Cuts ------- //

  
  h_DeltaBecut[0]=new TH2F("h_DeltaBe_0","Proton ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  h_DeltaBecut[1]=new TH2F("h_DeltaBe_1","Kaon ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.15, 0.15);
  h_DeltaBecut[2]=new TH2F("h_DeltaBe_2","Pion ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  
  h_BeVSpcut[0]=new TH2F("h_BeVSp_0","Proton ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSpcut[1]=new TH2F("h_BeVSp_1","Kaon ;p [GeV/c];#beta;",200, 0, 3, 200, 0, 1);
  h_BeVSpcut[2]=new TH2F("h_BeVSp_2","Pion ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);

   //-----Correlation Theta-Phi, ----------//

  h_ThePhi_proton=new TH2F("h_ThePhi_proton","Proton ;#theta #circ; #phi #circ;",200, -180, 180, 200, 0, 150);
  h_ThePhi_kaon=new TH2F("h_ThePhi_kaon","Kaon ;#theta #circ; #phi #circ;",200,  -180, 180, 200, 0, 150);
  h_ThePhi_pion=new TH2F("h_ThePhi_pion","Pion ;#theta #circ; #phi #circ;", 200,  -180, 180, 200, 0, 150);

 //-------------Fiduciary cuts-------------------------//
  h_ThePhi_protoncut=new TH2F("h_ThePhi_proton","Proton ;#theta #circ; #phi #circ;",200, -180, 180, 200, 0, 150);
  h_ThePhi_kaoncut=new TH2F("h_ThePhi_kaon","Kaon ; #theta #circ; #phi #circ;",200,  -180, 180, 200, 0, 150);
  h_ThePhi_pioncut=new TH2F("h_ThePhi_pion","Pion ; #theta #circ; #phi #circ", 200,  -180, 180, 200, 0, 150);


  
  //------------------ Photons, Delta T  ------------------ // 
  
  h_DeltaTall[0]=new TH1F("h_DeltaTall_0","Kaon ;#Delta t [ns];counts;", 200, -10, 10);
  h_DeltaTall[1]=new TH1F("h_DeltaTall_1","Pion ;#Delta t [ns];counts;", 200, -10, 10);


  h_DeltaTallvsp[0]=new TH2F("h_DeltaTallvsp_0","Kaon ;p [GeV/c];#Delta t [ns];", 200, 0, 3, 200, -10, 10);
  h_DeltaTallvsp[1]=new TH2F("h_DeltaTallvsp_1","Pion ;p [GeV/c];#Delta t [ns];", 200, 0, 3, 200, -10, 10);

  //------------Delta T with Cuts ----------- //
  
  h_DeltaT[0]=new TH1F("h_DeltaT_0","Kaon ;#Delta t [ns];counts;", 200, -10, 10);
  h_DeltaT[1]=new TH1F("h_DeltaT_1","Pion ;#Delta t [ns];counts;", 200, -10, 10);

  //--------------- Energy loss ----------- //

  h_eloss[0]=new TH1F("h_eloss_0","Proton; eloss",50, 0, 100);
  h_eloss[1]=new TH1F("h_eloss_1","Kaon; eloss",50, 0, 20);
  h_eloss[2]=new TH1F("h_eloss_2","Pion; eloss",50, 0, 20);

  //-------------- Reconstruction --------- //

  h_MissingMass = new TH1F("h_missingmass",
			   "Missing mass (Neutron); Missing Mass [GeV/c^{2}]; counts",
			   100, 0.7, 1.2);
  
  h_MissingMass_kaonpion = new TH1F("h_missingmass_kaonpion",
				    "Missing mass Kaon_pion ; Missing Mass [GeV/c^{2}]; counts",
				    100, 0.7, 1.2);
  
  h_MissingMass_vsMissingMasskaonpion = new TH2F("MissingMass_correlation",
						 "Missing Mass Correlation; W-Kaon [GeV/c^{2}];  W-Pion [GeV/c^{2}]",
						 100, 0.7, 1.2, 100, 0.7, 1.2);
  
  h_MissingP = new TH1F("h_missingp",
			"Missing momentum (Neutron); Missing momentum [GeV/c]; counts",
			100, 0.0, 1.5);
  
  h_MissingPvsMass = new TH2F("h_missingpvsm",
			      "Missing W vs P; Missing momentum [GeV/c]; Missing Mass [GeV/c^{2}]",
			      100, 0.0, 1.5, 100, 0.7, 1.2);
  
  h_MissingMassvsSigmaMass = new TH2F("h_missingmvsSigma",
				      "Invariant Mass Sigma vs Missing Mass Neutron (Kaon); Mass [GeV/c^{2}]; Mass [GeV/c^{2}]",
				      100,1.0,1.5, 100, 0.7, 1.2);
  
  h_MissingPvsSigmaMass = new TH2F("h_missingPvsSigma",
				   "Invariant Mass Sigma vs Missing Momentum Neutron (Kaon); Mass [GeV/c^{2}]; p [GeV/c] ",
				   100,1.0,1.5, 100, 0.0, 1.5);
  
  h_InvariantMass = new TH1F("h_InvariantMass",
			     "Invariant mass Sigma With Cuts; Mass [GeV/c^{2}]; counts ",
			     100, 1.0, 1.5);
  
  h_LambdaMass = new TH1F("h_LambdaMass",
			  "Invariant mass Lambda With Cuts; Mass [GeV/c^{2}]; counts ",
			  100, 1.08, 1.16);


  h_DeltaBVSInvariantMass = new TH2F("h_DeltaBVSInvariantMass",
				     "",
				     100,0, 2,100,-0.17, 0.17);
  
  h_DeltaBVSMissingMass = new TH2F("h_DeltaBVSMissingMass",
				   "",
				   100,0, 2,100,-0.17, 0.17);
  
  h_DeltaBVSMissingMomentum = new TH2F("h_DeltaBVSMissingMomentum",
				       "",
				       100,0, 1.5,100,-0.17, 0.17);


  //0.9395601,1.052416,0.08003898,0.03282563,0,360,45
  myEllipse = new TEllipse(0.9395601,1.052416,0.04,0.02,0,360,45);
}


void Histograms::DoCanvas(){

  
  //---------------- Delta B Without Cuts ----------------- //


  TCanvas *c0=new TCanvas("c0","Delta Beta y Beta", 1450, 500);
  c0->Divide(2,1);
  c0->cd(1);
  h_DeltaBe[0]->Draw("colz");
  h_DeltaBe[0]->Fit(FFits[0]);
  //Fit functions
  FFname[0]=FFits[0]->GetName();
  FFname[0].insert(0,"-1.0*");
  FFitsminus[0] = new TF1("DBProtonFitminus",FFname[0].c_str(),0,3);
  FFits[0]->Draw("same");
  FFitsminus[0]->Draw("same");
  

  c0->cd(2);
  h_BeVSp[0]->Draw("colz");
  BeVSpProton->Draw("same");
  c0->SaveAs("imagenes/ProtonDB_VS_P.eps");
  
  TCanvas *c01=new TCanvas("c01","Delta Beta y Beta", 1450, 500);
  c01->Divide(2,1);
  c01->cd(1);
  h_DeltaBe[1]->Draw("colz");
  h_DeltaBe[1]->Fit(FFits[1]);
  //Fit functions
  FFname[1]=FFits[1]->GetName();
  FFname[1].insert(0,"-1.0*");
  FFitsminus[1] = new TF1("DBProtonFitminus",FFname[1].c_str(),0,3);
  FFits[1]->Draw("same");
  FFitsminus[1]->Draw("same");

  
  c01->cd(2);
  h_BeVSp[1]->Draw("colz");
  BeVSpKaon->Draw("same");
  c01->SaveAs("imagenes/KaonDB_VS_P.eps");
  
  TCanvas *c02=new TCanvas("c02","Delta Beta y Beta", 1450, 500);
  c02->Divide(2,1);
  c02->cd(1);
  h_DeltaBe[2]->Draw("colz");
  h_DeltaBe[2]->Fit(FFits[2]);
  //Fit functions
  FFname[2]=FFits[2]->GetName();
  FFname[2].insert(0,"-1.0*");
  FFitsminus[2] = new TF1("DBProtonFitminus",FFname[2].c_str(),0,3);
  FFits[2]->Draw("same");
  FFitsminus[2]->Draw("same");
  
  c02->cd(2);
  h_BeVSp[2]->Draw("colz");
  BeVSpPion->Draw("same");
  c02->SaveAs("imagenes/PionDB_VS_P.eps");
  
  /*
  //---------------- Delta B With Cuts ----------------- //

  TCanvas *c0cut=new TCanvas("c0cut","Delta Beta and Beta with Cuts", 1450, 500);
  c0cut->Divide(2,1);
  c0cut->cd(1);
  h_DeltaBecut[0]->Draw("colz");
  c0cut->cd(2);
  h_BeVSpcut[0]->Draw("colz");
  BeVSpProton->Draw("same");
  
  TCanvas *c01cut=new TCanvas("c01cut","Delta Beta and Beta with Cuts", 1450, 500);
  c01cut->Divide(2,1);
  c01cut->cd(1);
  h_DeltaBecut[1]->Draw("colz");
  c01cut->cd(2);
  h_BeVSpcut[1]->Draw("colz");
  BeVSpKaon->Draw("same");
  
  TCanvas *c02cut=new TCanvas("c02cut","Delta Beta and Beta with Cuts", 1450, 500);
  c02cut->Divide(2,1);
  c02cut->cd(1);
  h_DeltaBecut[2]->Draw("colz");
  c02cut->cd(2);
  h_BeVSpcut[2]->Draw("colz");
  BeVSpPion->Draw("same");
  

  //------------------ Delta de T without Cuts ---------------- //
  
  TCanvas *c1=new TCanvas("c1","Delta T", 900, 500);
  c1->Divide(2,1);
  c1->cd(1);
  h_DeltaTall[0]->Draw(); 
  c1->cd(2);
  h_DeltaTall[1]->Draw();

  TCanvas *c11=new TCanvas("c11","Delta T", 900, 500);
  c11->Divide(2,1);
  c11->cd(1);
  h_DeltaTallvsp[0]->Draw("colz"); 
  c11->cd(2);
  h_DeltaTallvsp[1]->Draw("colz");


  TCanvas *c12=new TCanvas("c12","Delta T", 900, 500);
  c12->Divide(2,1);
  c12->cd(1);
  h_DeltaT[0]->Draw();
  c12->cd(2);
  h_DeltaT[1]->Draw();

  // ------------- Energy Loss ---------------- //
  
  TCanvas *c2=new TCanvas("c2","Delta Energy loss", 900, 500);
  c2->Divide(3,1);
  c2->cd(1);
  h_eloss[0]->Draw();
  c2->cd(2);
  h_eloss[1]->Draw();
  c2->cd(3);
  h_eloss[2]->Draw();

  */
  
  //-------------- Reconstruction --------- //
  
  TCanvas *c3=new TCanvas("c3","Missing mass", 900, 500);
  c3->Divide(2,1);
  c3->cd(1);
  h_MissingMass->Draw();
  c3->cd(2);
  h_MissingMass_kaonpion->Draw();
  c3->SaveAs("imagenes/MissingMass.eps");

  //h_MissingP->Draw();
  
  TCanvas *c31=new TCanvas("c31","Missing mass", 900, 500);
  c31->Divide(2,1);
  c31->cd(1);
  h_InvariantMass->Draw();
  c31->cd(2);
  h_LambdaMass->Draw();
  c31->SaveAs("imagenes/InvariantMass.eps");

  double x=-145,b=-155,y1=150,y2=0; //Coordenadas de las líneas

 
  TCanvas *c4=new TCanvas("c4","Theta-Phi correlation", 900, 500);
  c4->Divide(1,3);
  c4->cd(1);
  h_ThePhi_proton->Draw("colz");
   vector<TLine*> lineasproton(12);
 
  
  for (int i=0; i<11; i+=2) {
    
    lineasproton.at(i)= new TLine(x, y1, x, y2);    
    if (x > 0)
      lineasproton.at(i+1)= new TLine(x+10, y1, x+10, y2);
    
    else
      lineasproton.at(i+1)= new TLine(x-10, y1, x-10, y2);

    if (x == -25)
      x+=50;
    else
      x+=60;

    lineasproton.at(i)->SetLineWidth(2);
    lineasproton.at(i+1)->SetLineWidth(2);
    lineasproton.at(i)->SetLineColor(2);
    lineasproton.at(i+1)->SetLineColor(2);
    
    lineasproton.at(i)->Draw("same");
    lineasproton.at(i+1)->Draw("same");
  }
  c4->cd(2);
  h_ThePhi_kaon->Draw("colz");
  x=-145;
  vector<TLine*> lineaskaon(12);
  
  for (int i=0; i<11; i+=2) {
    
    lineaskaon.at(i)= new TLine(x, y1, x, y2);    
    if (x > 0)
      lineaskaon.at(i+1)= new TLine(x+10, y1, x+10, y2);
    
    else
      lineaskaon.at(i+1)= new TLine(x-10, y1, x-10, y2);

    if (x == -25)
      x+=50;
    else
      x+=60;

    lineaskaon.at(i)->SetLineWidth(2);
    lineaskaon.at(i+1)->SetLineWidth(2);
    lineaskaon.at(i)->SetLineColor(2);
    lineaskaon.at(i+1)->SetLineColor(2);
    
    lineaskaon.at(i)->Draw("same");
    lineaskaon.at(i+1)->Draw("same");
  }
  c4->cd(3);
  h_ThePhi_pion->Draw("colz");
  x=-145;
  vector<TLine*> lineaspion(12);
   
  for (int i=0; i<11; i+=2) {
    
    lineaspion.at(i)= new TLine(x, y1, x, y2);    
    if (x > 0)
      lineaspion.at(i+1)= new TLine(x+10, y1, x+10, y2);
    
    else
      lineaspion.at(i+1)= new TLine(x-10, y1, x-10, y2);

    if (x == -25)
      x+=50;
    else
      x+=60;

    lineaspion.at(i)->SetLineWidth(2);
    lineaspion.at(i+1)->SetLineWidth(2);
    lineaspion.at(i)->SetLineColor(2);
    lineaspion.at(i+1)->SetLineColor(2);
    
    lineaspion.at(i)->Draw("same");
    lineaspion.at(i+1)->Draw("same");
  }

   
  
  
  c4->SaveAs("imagenes/Fiduciarycuts.eps");

  

  TCanvas *c5=new TCanvas("c5","Fiduciary cuts", 900, 500);
  c5->Divide(1,3);
  c5->cd(1);
  h_ThePhi_protoncut->Draw("colz");
  c5->cd(2);
  h_ThePhi_kaoncut->Draw("colz");
  c5->cd(3);
  h_ThePhi_pioncut->Draw("colz");  
  
  TCanvas *c32=new TCanvas("c32","Missing mass", 900, 500);
  c32->cd(1);
  h_MissingMass_vsMissingMasskaonpion->Draw("colz");
  myEllipse->SetFillStyle(0);
  myEllipse->SetLineColor(kRed);
  myEllipse->Draw("same");
  c32->SaveAs("imagenes/Ellipse.eps");
  
  TCanvas *c40 = new TCanvas("c40","Delta Beta Vs Missingmass and Invariantmass", 900, 500);;
  c40->Divide(3,1);
  c40->cd(1);
  h_DeltaBVSInvariantMass->Draw("colz");
  c40->cd(2);
  h_DeltaBVSMissingMass->Draw("colz");
  c40->cd(3);
  h_DeltaBVSMissingMomentum->Draw("colz");
}



#endif
