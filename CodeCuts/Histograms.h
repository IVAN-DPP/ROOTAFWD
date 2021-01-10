#include "include/Libraries.h"

#ifndef HISTOGRAMS_C
#define HISTOGRAMS_C

using namespace std;

//Friends Functions
void LinesPTCuts();
void NameLinesInv(double, double, int, int);
void TPaveStateModify(TH1*,TH1*);
vector<double> HistoBinning(TH1 *, int, double, double, double,double);
Double_t fitf(Double_t *,Double_t *);
class Histograms{
  //Friends Functions
  friend void LinesPTCuts();       //Functions for do the lines in Theta-Phi Correlations
  friend void NameLinesInv(double, double, int, int);
  friend void TPaveStateModify(TH1*,TH1*);
  friend vector<double> HistoBinning(TH1 *, int, double, double, double,double); //Function for binning the graphic Cos(Kaon.Theta)
  friend Double_t fitf(Double_t *,Double_t *); //Function to asymmetry compute
  
protected:

  TH1F *h_Vertex                            = NULL;
   
  TH2F *h_DeltaBe[3]                        = {};
  TH2F *h_BeVSp[3]                          = {};
  TH2F *h_BeVSpT                            = NULL;
  TH2F *h_DeltaBecut[3]                     = {};
  TH2F *h_BeVSpcut[3]                       = {};
  
  
  TF1 *FFits[3]                             = {};               //Fits for do cuts
  TF1 *FFitsminus[3]                        = {};
  string FFname[3]                          = {};
  
  TH1F *h_DeltaTall[2]                      = {};
  TH2F *h_DeltaTallvsp[2]                   = {};
  TH1F *h_DeltaT[2]                         = {};
  TH1F *h_eloss[3]                          = {};
  
  //--> Beta vs P with PDG MASSES  

  TF1 *BeVSp[3]                             = {};


  //---- Get Coherent Edge ---- //

  TH1F *h_TagrEpho[3]                       = {};
  
  //---- Missing mass ----//
  
  TH1F *h_MissingMass                     = NULL;
  TH1F *h_MissingMass_kaonpion            = NULL;
  TH1F *h_MissingMasscut                  = NULL;
  TH1F *h_MissingMass_kaonpioncut         = NULL;
  
  TH2F *h_MissingMass_vsMissingMasskaonpion[2] = {};


  TH1F *h_MissingP[2]                      = {};
  TH1F *h_MissingPcut[2]                   = {};
  TH2F *h_MissingPvsMass[2]                = {};
  TH2F *h_MissingMassvsSigmaMass          = NULL;
  TH2F *h_MissingPvsSigmaMass             = NULL;

  TH1F *h_InvariantMass                   = NULL;
  TH1F *h_InvariantMasscut[4]                = {};

  //---- Lambda and Lambda Fit ---- //
  
  TH1F *h_LambdaMass                      = NULL;
  TF1 *lamdaMassFit                       = NULL; 
  
  
  TH2F *h_DeltaBVSInvariantMass           = NULL;
  TH2F *h_DeltaBVSMissingMass             = NULL;
  TH2F *h_DeltaBVSMissingMomentum         = NULL;

  //--- Ellipse Cuts ---- //
  
  TEllipse *myEllipse                     = NULL;
  double radx=0.034, rady=0.02, offsetx=0.937, offsety=1.047, angle=70*TMath::DegToRad();
  
  //-----Correlation Theta-Phi, ----------//
 
  TH2F *h_ThePhi[3]                         = {};

  //--- Fiduciary cuts ---//
  TH2F *h_ThePhicut[3]                      = {};
  
  //----Costheta-Kaon Boost-------------------//
  TH1F *h_KCosThetaCM                     = NULL;

  //-----Kaon phi-------------------//
  TH1F *h_kaonPhiPA[2]                      = {};
  TH1F *h_kaonPhiPE[2]                      = {};
    
  //---------Function to do asymmetry fit----//
  vector<vector<double> > MEASPhi{2}; //2 is the number of binning
  vector<vector<double> > MEASGammaP{2};
    
  TF1 *FuncAsym                            =NULL;
  TH1F *h_Asym[2]                         ={};
    
public:
  Histograms(){}
  void DoHistograms();
  void DoCanvas();
};


void Histograms::DoHistograms(){

  h_Vertex = new TH1F("h_Vertex","Vertex; distance [cm]; counts",200,0,-40);

  
  
  //------------------ Delta Beta ---------------//
  
  h_DeltaBe[0]=new TH2F("h_DeltaBe_0","Proton ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  h_DeltaBe[1]=new TH2F("h_DeltaBe_1","Kaon ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.15, 0.15);
  h_DeltaBe[2]=new TH2F("h_DeltaBe_2","Pion ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);

  
  h_BeVSp[0]=new TH2F("h_BeVSp_0","Proton ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp[1]=new TH2F("h_BeVSp_1","Kaon ;p [GeV/c];#beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp[2]=new TH2F("h_BeVSp_2","Pion ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);

  h_BeVSpT = new TH2F("h_BeVSpT","Particles ;Momentum (p) [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);

  //--> Beta vs P with PDG MASSES

  BeVSp[0] = new TF1("BeVSpProton","x/std::sqrt(std::pow(x,2)+std::pow(0.938,2))",0,3);
  BeVSp[0]->SetLineColor(kBlack);
  BeVSp[0]->SetLineStyle(2);
  
  BeVSp[1] = new TF1("BeVSpKaon","x/std::sqrt(std::pow(x,2)+std::pow(0.493,2))",0,3);
  BeVSp[1]->SetLineColor(kBlack);
  BeVSp[1]->SetLineStyle(2);

  
  BeVSp[2] = new TF1("BeVSpPion","x/std::sqrt(std::pow(x,2)+std::pow(0.139,2))",0,3);
  BeVSp[2]->SetLineColor(kBlack);
  BeVSp[2]->SetLineStyle(2);
  
  //------- Fits for do cuts in Delta B ----------- //

  FFits[0] = new TF1("DBProtonFit","0.02",0,3);
  FFits[1] = new TF1("DBKaonFit","0.025", 0, 3);
  FFits[2] = new TF1("DBPionFit","0.05",0,3);


  
  
  //-------- Delta Beta Cuts ------- //

  
  h_DeltaBecut[0] = new TH2F("h_DeltaBe_0cut","Proton ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  h_DeltaBecut[1] = new TH2F("h_DeltaBe_1cut","Kaon ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.15, 0.15);
  h_DeltaBecut[2] = new TH2F("h_DeltaBe_2cut","Pion ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  
  h_BeVSpcut[0] = new TH2F("h_BeVSp_0cut","Proton ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSpcut[1] = new TH2F("h_BeVSp_1cut","Kaon ;p [GeV/c];#beta;",200, 0, 3, 200, 0, 1);
  h_BeVSpcut[2] = new TH2F("h_BeVSp_2cut","Pion ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);

  //-----Correlation Theta-Phi, ----------//

  h_ThePhi[0] = new TH2F("h_ThePhi_proton","Proton ;#theta #circ; #phi #circ;",200, -180, 180, 200, 0, 150);
  h_ThePhi[1] = new TH2F("h_ThePhi_kaon","Kaon ;#theta #circ; #phi #circ;",200,  -180, 180, 200, 0, 150);
  h_ThePhi[2] = new TH2F("h_ThePhi_pion","Pion ;#theta #circ; #phi #circ;", 200,  -180, 180, 200, 0, 150);

  //-------------Fiducial cuts-------------------------//

  
  h_ThePhicut[0] = new TH2F("h_ThePhi_protoncut","Proton ;#theta #circ; #phi #circ;",200, -180, 180, 200, 0, 150);
  h_ThePhicut[1] = new TH2F("h_ThePhi_kaoncut","Kaon ; #theta #circ; #phi #circ;",200,  -180, 180, 200, 0, 150);
  h_ThePhicut[2] = new TH2F("h_ThePhi_pioncut","Pion ; #theta #circ; #phi #circ", 200,  -180, 180, 200, 0, 150);


  
  //------------------ Photons, Delta T  ------------------ // 
  
  h_DeltaTall[0]=new TH1F("h_DeltaTall_0","Kaon ;#Delta t [ns];counts;", 200, -10, 10);
  h_DeltaTall[1]=new TH1F("h_DeltaTall_1","Pion ;#Delta t [ns];counts;", 200, -10, 10);


  h_DeltaTallvsp[0]=new TH2F("h_DeltaTallvsp_0","Kaon ;p [GeV/c];#Delta t [ns];", 200, 0, 3, 200, -10, 10);
  h_DeltaTallvsp[1]=new TH2F("h_DeltaTallvsp_1","Pion ;p [GeV/c];#Delta t [ns];", 200, 0, 3, 200, -10, 10);

  //------------Delta T with Cuts ----------- //
  
  h_DeltaT[0]=new TH1F("h_DeltaT_0","Kaon ;#Delta t [ns];counts;", 200, -10, 10);
  h_DeltaT[1]=new TH1F("h_DeltaT_1","Pion ;#Delta t [ns];counts;", 200, -10, 10);

  //--------------- Energy loss ----------- //

  h_eloss[0]=new TH1F("h_eloss_0","Proton; eloss [MeV];counts;",50, 0, 100);
  h_eloss[1]=new TH1F("h_eloss_1","Kaon; eloss [MeV];counts;",50, 0, 20);
  h_eloss[2]=new TH1F("h_eloss_2","Pion; eloss [MeV];counts;",50, 0, 20);

  //---- Get Coherent Edge ---- //

  h_TagrEpho[0] = new TH1F("h_Tagrepho1","",100,0,2000);
  h_TagrEpho[1] = new TH1F("h_Tagrepho2","",100,0,2000);
  h_TagrEpho[2] = new TH1F("h_Tagrepho3","",100,0,2000);
  
  //-------------- Reconstruction --------- //

  h_MissingMass = new TH1F("h_missingmass",
			   "Missing mass Neutron; #gamma d #rightarrow K^{+} #pi^{-} X p [GeV/c^{2}]; counts",
			   100, 0.7, 1.2);
  
  h_MissingMass_kaonpion = new TH1F("h_missingmass_kaonpion",
				    "Missing mass Neutron ; #gamma d #rightarrow #pi^{+} #pi^{-} X p [GeV/c^{2}]; counts",
				    100, 0.7, 1.2);

  h_MissingMasscut = new TH1F("h_missingmasscut",
			      "Missing mass Neutron; #gamma d #rightarrow K^{+} #pi^{-} X p [GeV/c^{2}]; counts",
			      100, 0.7, 1.2);
  
  h_MissingMass_kaonpioncut = new TH1F("h_missingmass_kaonpioncut",
				       "Missing mass Neutron ; #gamma d #rightarrow #pi^{+} #pi^{-} X p [GeV/c^{2}]; counts",
				       100, 0.7, 1.2);

  
  h_MissingMass_vsMissingMasskaonpion[0] = new TH2F("MissingMass_correlation",
						    "Missing Mass Correlation; #gamma d #rightarrow K^{+} #pi^{-} X p [GeV/c^{2}];  #gamma d #rightarrow #pi^{+} #pi^{-} X p [GeV/c^{2}]",
						    100, 0.7, 1.2, 100, 0.7, 1.2);

  h_MissingMass_vsMissingMasskaonpion[1] = new TH2F("MissingMass_correlation_Cut",
						    "Missing Mass Correlation; #gamma d #rightarrow K^{+} #pi^{-} X p [GeV/c^{2}];  #gamma d #rightarrow #pi^{+} #pi^{-} X p [GeV/c^{2}]",
						    100, 0.7, 1.2, 100, 0.7, 1.2);
  
  h_MissingP[0] = new TH1F("h_missingpSigma",
			   "Missing momentum Neutron; Missing momentum [GeV/c]; counts",
			   100, 0.0, 1.);
  h_MissingP[1] = new TH1F("h_missingpLambda",
			   "Missing momentum Neutron; Missing momentum [GeV/c]; counts",
			   100, 0.0, 1.);

  h_MissingPcut[0] = new TH1F("h_missingpcutSigma",
			      "Missing momentum Neutron; Missing momentum [GeV/c]; counts",
			      100, 0.0, 1.5);
  h_MissingPcut[1] = new TH1F("h_missingpcutLambda",
			      "Missing momentum Neutron; Missing momentum [GeV/c]; counts",
			      100, 0.0, 1.5);
  



  
  h_MissingPvsMass[0] = new TH2F("h_missingpvsmSigma",
				 "Missing Momentum Vs Missing Mass (Lambda cut); #gamma d #rightarrow K^{+} #pi^{-} X p [GeV/c^{2}]; Missing Momentum (p) [GeV/c]",
				 100, 0.7, 1.2, 100, 0.0, 1.5);

  h_MissingPvsMass[1] = new TH2F("h_missingpvsmLambda",
				 "Missing Momentum Vs Missing Mass (Sigma cut); #gamma d #rightarrow K^{+} #pi^{-} X p [GeV/c^{2}]; Missing Momentum (p) [GeV/c]",
				 100, 0.7, 1.2, 100, 0.0, 1.5);


  // ------- This aren't important -------- //
  h_MissingMassvsSigmaMass = new TH2F("h_missingmvsSigma",
				      "Invariant Mass Sigma vs Missing Mass Neutron (Kaon); Invariant Mass (#pi^{-} n) [GeV/c^{2}]; #gamma d #rightarrow K^{+} #pi^{-} X p [GeV/c^{2}]",
				      100,1.0,1.5, 100, 0.7, 1.2);
  
  h_MissingPvsSigmaMass = new TH2F("h_missingPvsSigma",
				   "Invariant Mass Sigma vs Missing Momentum Neutron (Kaon); Mass [GeV/c^{2}]; p [GeV/c] ",
				   100,1.0,1.5, 100, 0.0, 1.5);
  
  // ----------------------------------- //

  h_InvariantMass = new TH1F("h_InvariantMass",
			     "Invariant mass ; Mass [GeV/c^{2}]; counts ",
			     100, 1.0, 1.5);
  
  
  h_InvariantMasscut[0] = new TH1F("h_InvariantMasscut3Sig",
				   "Invariant mass Sigma With Cuts; Mass [GeV/c^{2}]; counts ",
				   100, 1.0, 1.5);

  h_InvariantMasscut[1] = new TH1F("h_InvariantMasscut4Sig",
				   "Invariant mass Sigma With Cuts; Mass [GeV/c^{2}]; counts ",
				   100, 1.0, 1.5);

  h_InvariantMasscut[2] = new TH1F("h_InvariantMasscut5Sig",
				   "Invariant mass Sigma With Cuts; Mass [GeV/c^{2}]; counts ",
				   100, 1.0, 1.5);
  
  h_InvariantMasscut[3] = new TH1F("h_InvariantMasscut_ACP",
				   "Invariant mass Sigma With Cuts; Mass [GeV/c^{2}]; counts ",
				   100, 1.0, 1.5);
  
  //-------- Lambda and Lambda Fit ------ //
  
  h_LambdaMass = new TH1F("h_LambdaMass",
			  "Invariant mass Lambda With Cuts; Mass [GeV/c^{2}]; counts ",
			  100, 1.08, 1.16);

  lamdaMassFit = new TF1("lamdaMassFit","gaus",1.08,1.16);
  
  //----------------------------------//
  h_DeltaBVSInvariantMass = new TH2F("h_DeltaBVSInvariantMass",
				     "",
				     100,0, 2,100,-0.17, 0.17);
  
  h_DeltaBVSMissingMass = new TH2F("h_DeltaBVSMissingMass",
				   "",
				   100,0, 2,100,-0.17, 0.17);
  
  h_DeltaBVSMissingMomentum = new TH2F("h_DeltaBVSMissingMomentum",
				       "",
				       100,0, 1.5,100,-0.17, 0.17);


  myEllipse = new TEllipse(offsetx,offsety,radx,rady,0,360,70);


  //--------------KaonCosTheta Boost-------------//
  h_KCosThetaCM = new TH1F("h_KCosThetaCM","Kaon Cos_Theta Boost", 100, -0.4, 2.);

  
  //-----------------Kaon phi-------------------//
  
  h_kaonPhiPA[0] = new TH1F("h_kaonPhiPA[0]","First partition Kaon_Phi (PARA)",200, -180, 180);
  h_kaonPhiPA[1] = new TH1F("h_kaonPhiPA[1]","Second partition Kaon_Phi (PARA)",200, -180, 180);
  h_kaonPhiPE[0] = new TH1F("h_kaonPhiPE[0]","First partition Kaon_Phi (PERP)",200, -180, 180);
  h_kaonPhiPE[1] = new TH1F("h_kaonPhiPE[1]","Second partition Kaon_Phi (PERP)", 200, -180, 180);
  
  //---------------  Fit to Asymmetry ------------//

  FuncAsym = new TF1("FuncAsym",fitf,-180,180,4);
  FuncAsym->SetParameters(1.2,1,0.7,0.5);
  FuncAsym->SetParNames("FR","PR","Pave","Sigma");

  
  //----------------- Histograms Asymmetry -------//

  
  
  h_Asym[0] = new TH1F("h_Asym[0]","Asymmetry first partition", 100, -180, 180);
  h_Asym[1] = new TH1F("h_Asym[1]","Asymmetry Second partition",100, -180, 180);
}


void Histograms::DoCanvas(){

  gStyle->SetOptStat("e");
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  
  TCanvas *c00 = new TCanvas("c00","Vertex",950,500);
  c00->cd(0);
  h_Vertex->SetTitleSize(0.047, "xy");
  h_Vertex->Draw();
  TLine *VertexLines[2]={};
  VertexLines[0] = new TLine(-1.0,h_Vertex->GetMaximum(),-1.0,0);
  VertexLines[1] = new TLine(-39.0,h_Vertex->GetMaximum(),-39.0,0);

  VertexLines[0]->SetLineWidth(2);
  VertexLines[1]->SetLineWidth(2);
  VertexLines[0]->SetLineColor(2);
  VertexLines[1]->SetLineColor(2);

  VertexLines[0]->Draw("same");
  VertexLines[1]->Draw("same");

  c00->SaveAs("imagenes/Vertex.eps");

  //---------------- Delta B With/out Cuts ----------------- //

  //------ Proton ------//
  
  TCanvas *c0DB=new TCanvas("c0DB","Delta Beta Without Cuts", 1000, 500);
  c0DB->cd(1);
  //h_DeltaBe[0]->SetLabelSize(0.056, "XY");
  h_DeltaBe[0]->SetTitleSize(0.047, "XY");
  h_DeltaBe[0]->Draw("colz");
  h_DeltaBe[0]->Fit(FFits[0]);
  //Fit functions
  FFname[0]=FFits[0]->GetName();
  FFname[0].insert(0,"-1.0*");
  FFitsminus[0] = new TF1("DBProtonFitminus",FFname[0].c_str(),0,3);
  FFits[0]->Draw("same");
  FFitsminus[0]->Draw("same");
  c0DB->SaveAs("imagenes/ProtonDB_VS_P.eps");
  
  TCanvas *c0DBC=new TCanvas("c0DBC","Delta Beta With Cuts", 1000, 500);
  c0DBC->cd(1);
  h_DeltaBecut[0]->SetTitleSize(0.047, "XY");
  h_DeltaBecut[0]->Draw("colz");
  c0DBC->SaveAs("imagenes/ProtonDB_VS_P_C.eps");

  //------ Kaon ------//
  
  TCanvas *c1DB=new TCanvas("c1DB","Delta Beta Without Cuts", 1000, 500);
  c1DB->cd(1);
  h_DeltaBe[1]->SetTitleSize(0.047, "XY");
  h_DeltaBe[1]->Draw("colz");
  h_DeltaBe[1]->Fit(FFits[1]);
  //Fit functions
  FFname[1]=FFits[1]->GetName();
  FFname[1].insert(0,"-1.0*");
  FFitsminus[1] = new TF1("DBKaonFitminus",FFname[0].c_str(),0,3);
  FFits[1]->Draw("same");
  FFitsminus[1]->Draw("same");
  c1DB->SaveAs("imagenes/KaonDB_VS_P.eps");

  TCanvas *c1DBC=new TCanvas("c1DBC","Delta Beta With Cuts", 1000, 500);
  c1DBC->cd(1);
  h_DeltaBecut[1]->SetTitleSize(0.047, "XY");
  h_DeltaBecut[1]->Draw("colz");
  c0DBC->SaveAs("imagenes/KaonDB_VS_P_C.eps");

  
  //------ Pion ------//
  
  TCanvas *c2DB=new TCanvas("c2DB","Delta Beta Without Cuts", 1000, 500);
  c2DB->cd(1);
  h_DeltaBe[2]->SetTitleSize(0.047, "XY");
  h_DeltaBe[2]->Draw("colz");
  h_DeltaBe[2]->Fit(FFits[2]);
  //Fit functions
  FFname[2]=FFits[1]->GetName();
  FFname[2].insert(0,"-1.0*");
  FFitsminus[2] = new TF1("DBPionFitminus",FFname[0].c_str(),0,3);
  FFits[2]->Draw("same");
  FFitsminus[2]->Draw("same");
  c2DB->SaveAs("imagenes/PionDB_VS_P.eps");

  TCanvas *c2DBC=new TCanvas("c2DBC","Delta Beta With Cuts", 1000, 500);
  c2DBC->cd(1);
  h_DeltaBecut[2]->SetTitleSize(0.047, "XY");
  h_DeltaBecut[2]->Draw("colz");
  c2DBC->SaveAs("imagenes/PionDB_VS_P_C.eps");


  //---------------- B With/out Cuts ----------------- //

  gStyle->SetStatY(0.94);
  //------ Proton ------//
  
  TCanvas *c0B=new TCanvas("c0B","Beta Without Cuts", 1000, 500);
  c0B->cd(1);
  h_BeVSp[0]->SetTitleSize(0.047, "XY");
  h_BeVSp[0]->Draw("colz");  
  BeVSp[0]->Draw("same");
  c0B->SaveAs("imagenes/ProtonB_VS_P.eps");

  TCanvas *c0BC=new TCanvas("c0BC","Beta With Cuts", 1000, 500);
  c0BC->cd(1);
  h_BeVSpcut[0]->SetTitleSize(0.047, "XY");
  h_BeVSpcut[0]->Draw("colz");  
  BeVSp[0]->Draw("same");
  
  c0BC->SaveAs("imagenes/ProtonB_VS_P_C.eps");

  //------ Kaon ------//
  
  TCanvas *c1B=new TCanvas("c1B","Beta Without Cuts", 1000, 500);
  c1B->cd(1);
  h_BeVSp[1]->SetTitleSize(0.047, "XY");
  h_BeVSp[1]->Draw("colz");  
  BeVSp[1]->Draw("same");
  c0B->SaveAs("imagenes/KaonB_VS_P.eps");

  TCanvas *c1BC=new TCanvas("c1BC","Beta With Cuts", 1000, 500);
  c1BC->cd(1);
  h_BeVSpcut[1]->SetTitleSize(0.047, "XY");
  h_BeVSpcut[1]->Draw("colz");  
  BeVSp[1]->Draw("same");
  
  c1BC->SaveAs("imagenes/KaonB_VS_P_C.eps");

  
  //------ Pion ------//
  
  TCanvas *c2B=new TCanvas("c2B","Beta Without Cuts", 1000, 500);
  c2B->cd(1);
  h_BeVSp[2]->SetTitleSize(0.047, "XY");
  h_BeVSp[2]->Draw("colz");  
  BeVSp[2]->Draw("same");
  c2B->SaveAs("imagenes/PionB_VS_P.eps");
  TCanvas *c2BC=new TCanvas("c2BC","Beta With Cuts", 1000, 500);
  c2BC->cd(1);
  h_BeVSpcut[2]->SetTitleSize(0.047, "XY");
  h_BeVSpcut[2]->Draw("colz");  
  BeVSp[2]->Draw("same");
  
  c2BC->SaveAs("imagenes/PionB_VS_P_C.eps");

  // B Total

  TCanvas *cB=new TCanvas("cB","Beta without Cuts", 1000, 500);
  cB->cd(1);
  h_BeVSpT->SetTitleSize(0.043, "XY");
  h_BeVSpT->Draw("colz");
  BeVSp[0]->Draw("same");
  BeVSp[1]->Draw("same");
  BeVSp[2]->Draw("same");
  
  cB->SaveAs("imagenes/B_VS_P.eps");
  

  
  //------------------ Delta de T without Cuts ---------------- //

  /*
    TCanvas *c1=new TCanvas("c1","Delta T", 900, 500);
    c1->Divide(2,1);
    c1->cd(1);
    h_DeltaT[0]->Draw(); 
    c1->cd(2);
    h_DeltaT[1]->Draw();

  
    TCanvas *c11=new TCanvas("c11","Delta T", 900, 500);
    c11->Divide(2,1);
    c11->cd(1);
    h_DeltaTallvsp[0]->Draw("colz"); 
    c11->cd(2);
    h_DeltaTallvsp[1]->Draw("colz");
  */

  gStyle->SetStatY(0.9);
  //--- DT Kaon ---//

  TCanvas *c0T=new TCanvas("c0T","Delta T", 1200, 500);
  TLine *DTL1= new TLine( -1.0, 0., -1.0 ,h_DeltaTall[0]->GetMaximum());
  TLine *DTL2= new TLine( 1.0, 0., 1.0 ,h_DeltaTall[0]->GetMaximum());
  DTL1->SetLineWidth(2);
  DTL2->SetLineWidth(2);
  DTL1->SetLineColor(2);
  DTL2->SetLineColor(2);
  c0T->cd(1);
  h_DeltaTall[0]->SetLabelSize(0.045, "XY");
  h_DeltaTall[0]->SetTitleSize(0.043, "XY");
  h_DeltaTall[0]->Draw();
  DTL1->Draw("same");
  DTL2->Draw("same");
  c0T->SaveAs("imagenes/DeltaTcut_Kaon.eps");

  //--- DT Pion ---//
  
  TCanvas *c1T=new TCanvas("c1T","Delta T", 1200, 500);
  c1T->cd(1);
  h_DeltaTall[1]->SetLabelSize(0.045, "XY");
  h_DeltaTall[1]->SetTitleSize(0.043, "XY");
  h_DeltaTall[1]->Draw();
  TLine *DTL1Pion= new TLine( -1.0, 0., -1.0 ,h_DeltaTall[1]->GetMaximum());
  TLine *DTL2Pion= new TLine( 1.0, 0., 1.0 ,h_DeltaTall[1]->GetMaximum());
  DTL1Pion->SetLineWidth(2);
  DTL2Pion->SetLineWidth(2);
  DTL1Pion->SetLineColor(2);
  DTL2Pion->SetLineColor(2);
  DTL1Pion->Draw("same");
  DTL2Pion->Draw("same");
  c1T->SaveAs("imagenes/DeltaTcut_Pion.eps");


  //----- Detector geometry (Fiduciary cuts)---------//

  gStyle->SetStatY(0.94);
    
  //------ Proton ------//
  
  TCanvas *c0TP=new TCanvas("c0TP","Theta-Phi correlation", 1200, 500);
  c0TP->cd(1);
  h_ThePhi[0]->SetTitleSize(0.05, "XY");
  h_ThePhi[0]->Draw("colz");
  LinesPTCuts();
  c0TP->SaveAs("imagenes/Fiduciary_Proton.eps");
  
  //------ Kaon ------//
  
  TCanvas *c1TP=new TCanvas("c1TP","Theta-Phi correlation", 1200, 500);
  c1TP->cd(1);
  h_ThePhi[1]->SetTitleSize(0.05, "XY");
  h_ThePhi[1]->Draw("colz");
  LinesPTCuts();
  c1TP->SaveAs("imagenes/Fiduciary_Kaon.eps");

  //------ Pion ------//
  
  TCanvas *c2TP=new TCanvas("c2TP","Theta-Phi correlation", 1200, 500);
  c2TP->cd(1);
  h_ThePhi[2]->SetTitleSize(0.05, "XY");
  h_ThePhi[2]->Draw("colz");
  LinesPTCuts();
  c2TP->SaveAs("imagenes/Fiduciary_Pion.eps");

  // 3 in 1

  TCanvas *c3TP=new TCanvas("c3TP","Theta-Phi correlation", 1200, 500);
  c3TP->Divide(1,3);
  c3TP->cd(1);
  h_ThePhi[0]->SetLabelSize(0.1, "XY");
  h_ThePhi[0]->SetTitleSize(0.05, "XY");
  h_ThePhi[0]->Draw("colz");
  LinesPTCuts();
  
  c3TP->cd(2);
  h_ThePhi[1]->SetLabelSize(0.1, "XY");
  h_ThePhi[1]->SetTitleSize(0.05, "XY");
  h_ThePhi[1]->Draw("colz");
  
  LinesPTCuts();

  c3TP->cd(3);
  h_ThePhi[2]->SetLabelSize(0.1, "XY");
  h_ThePhi[2]->SetTitleSize(0.05, "XY");
  h_ThePhi[2]->Draw("colz");
  LinesPTCuts();
  
  c3TP->SaveAs("imagenes/Fiduciary.eps");

  
  //-----Cuts due to detector geometry (Fiduciary cuts)---------//

  
  //------ Proton ------//
  
  TCanvas *c0TPC=new TCanvas("c0TPC","Fiduciary cuts", 1200, 500);
  c0TPC->cd(1);
  h_ThePhicut[0]->SetTitleSize(0.05, "XY");
  h_ThePhicut[0]->Draw("colz");
  LinesPTCuts();
  c0TPC->SaveAs("imagenes/Fiduciarycuts_Proton.eps");
  
  //------ Kaon ------//
  
  TCanvas *c1TPC=new TCanvas("c1TPC","Fiduciary cuts", 1200, 500);
  c1TPC->cd(1);
  h_ThePhicut[1]->SetTitleSize(0.05, "XY");
  h_ThePhicut[1]->Draw("colz");
  LinesPTCuts();
  c1TPC->SaveAs("imagenes/Fiduciarycuts_Kaon.eps");

  //------ Pion ------//
  
  TCanvas *c2TPC=new TCanvas("c2TPC","Fiduciary cuts", 1200, 500);
  c2TPC->cd(1);
  h_ThePhicut[2]->SetTitleSize(0.05, "XY");
  h_ThePhicut[2]->Draw("colz");
  LinesPTCuts();
  c2TPC->SaveAs("imagenes/Fiduciarycuts_Pion.eps");

  // 3 in 1
  
  TCanvas *c3TPC=new TCanvas("c3TPC","Fiduciary cuts", 1200, 500);
  c3TPC->Divide(1,3);
  c3TPC->cd(1);
  h_ThePhicut[0]->SetLabelSize(0.1, "XY");
  h_ThePhicut[0]->SetTitleSize(0.05, "XY");
  h_ThePhicut[0]->Draw("colz");
  c3TPC->cd(2);
  h_ThePhicut[1]->SetLabelSize(0.1, "XY");
  h_ThePhicut[1]->SetTitleSize(0.05, "XY");
  h_ThePhicut[1]->Draw("colz");
  c3TPC->cd(3);
  h_ThePhicut[2]->SetLabelSize(0.1, "XY");
  h_ThePhicut[2]->SetTitleSize(0.05, "XY");
  h_ThePhicut[2]->Draw("colz");  
  
  c3TPC->SaveAs("imagenes/Fiduciarycuts.eps");

  
  //-------------------- Energy Loss -------------------- //

  gStyle->SetStatY(0.9);

  //------ Proton ------//
  
  TCanvas *c0EL=new TCanvas("c0EL","Delta Energy loss", 1200, 500);
  c0EL->cd(1);
  h_eloss[0]->SetTitleSize(0.04, "XY");
  h_eloss[0]->Draw();
  c0EL->SaveAs("imagenes/Energyloss_Proton.eps");

  //------ Kaon ------//
  
  TCanvas *c1EL=new TCanvas("c1EL","Delta Energy loss", 1200, 500);
  c1EL->cd(1);
  h_eloss[1]->SetTitleSize(0.04, "XY");
  h_eloss[1]->Draw();
  c1EL->SaveAs("imagenes/Energyloss_Kaon.eps");

  //------ Pion ------//
  
  TCanvas *c2EL=new TCanvas("c2EL","Delta Energy loss", 1200, 500);
  h_eloss[2]->SetTitleSize(0.04, "XY");
  h_eloss[2]->Draw();
  c2EL->SaveAs("imagenes/Energyloss_Pion.eps");

  //---- Get Coherent Edge ---- //

  TCanvas *c20 = new TCanvas("c20","Coh Edge", 900,500);
  c20->Divide(3);
  c20->cd(1);
  h_TagrEpho[0]->Draw();
  c20->cd(2);
  h_TagrEpho[1]->Draw();
  c20->cd(3);
  h_TagrEpho[2]->Draw();
  
  //******************* Reconstruction **********************//

  //----------- MM before Cuts for MM correlation ------- //

  //------ Kaon ------//
  
  TCanvas *c0MM=new TCanvas("c0MM","Missing mass", 1200, 500);
  c0MM->cd(1);
  h_MissingMass->SetTitleSize(0.043, "XY");
  h_MissingMass->Draw();
  c0MM->SaveAs("imagenes/MissingMass.eps");

  //---------- MM Correlation without cut ---------------//
  
  TCanvas *c0ELL=new TCanvas("c0ELL","Correlation of MM", 900, 500);
  c0ELL->cd(1);
  h_MissingMass_vsMissingMasskaonpion[0]->SetTitleSize(0.045, "XY");  
  h_MissingMass_vsMissingMasskaonpion[0]->Draw("colz");
  myEllipse->SetFillStyle(0);
  myEllipse->SetLineColor(kRed);
  myEllipse->Draw("same");
  c0ELL->SaveAs("imagenes/Ellipse.eps");

  //---------- MM Correlation with cut ---------------//
  
  TCanvas *c0ELLC=new TCanvas("c0ELLC","Correlation of MM", 900, 500);
  c0ELLC->cd(1);
  h_MissingMass_vsMissingMasskaonpion[1]->SetTitleSize(0.045, "XY");  
  h_MissingMass_vsMissingMasskaonpion[1]->Draw("colz");
  myEllipse->SetFillStyle(0);
  myEllipse->SetLineColor(kRed);
  myEllipse->Draw("same");
  c0ELLC->SaveAs("imagenes/Ellipse_C.eps");

  //----------- MM After Cuts for MM correlation ------- //


  //------ Kaon ------//
 
  
  TCanvas *c0MMC=new TCanvas("c0MMC","Missing mass", 1200, 500);
  c0MMC->cd(1);
  
  h_MissingMass->SetLabelSize(0.045, "XY");
  h_MissingMass->SetTitleSize(0.043, "XY");
  h_MissingMass->Draw();
  TPaveStateModify(h_MissingMass,h_MissingMasscut);
  // gPad->Update();

  // TPaveStats *ps = (TPaveStats*)h_MissingMass->FindObject("stats");
  // ps->SetName("mystats");
  // TList *listOfLines = ps->GetListOfLines();

  // ostringstream streamObj;
  // streamObj << std::fixed;
  // streamObj << std::setprecision(0);
  // streamObj << h_MissingMasscut->GetEntries();
  // string StrObj = streamObj.str();

  // string EntriesCut = "Entries = " + StrObj;
  // TLatex *myt = new TLatex(0,0,EntriesCut.c_str());
  // myt ->SetTextFont(42);
  // //myt ->SetTextSize(0.02);
  // myt ->SetTextColor(kRed);
  // listOfLines->Add(myt);

  // ps->SetTextSize(0.1);
  // h_MissingMass->SetStats(0);
  
  h_MissingMasscut->SetFillColor(kRed-7);
  h_MissingMasscut->Draw("same");
  
  c0MMC->SaveAs("imagenes/MissingMass_Kaon.eps");

  //------ Pion + ------//
  TCanvas *c1MMC=new TCanvas("c1MMC","Missing mass", 1200, 500);
  c1MMC->cd(1);
  h_MissingMass_kaonpion->SetLabelSize(0.045, "XY");
  h_MissingMass_kaonpion->SetTitleSize(0.043, "XY");
  h_MissingMass_kaonpion->Draw();
  TPaveStateModify(h_MissingMass_kaonpion,h_MissingMass_kaonpioncut);
  h_MissingMass_kaonpioncut->SetFillColor(kRed-7);
  h_MissingMass_kaonpioncut->Draw("same");

  c1MMC->SaveAs("imagenes/MissingMass_Pion.eps");

  
  TCanvas *c310=new TCanvas("c31","Missing Momentum", 1450, 500);
  c310->cd(1);
  TLine *MMPL = new TLine(0.2,0,0.2,4300);
  MMPL->SetLineWidth(2);
  MMPL->SetLineColor(kBlue);
  h_MissingP[0]->SetLabelSize(0.053, "XY");
  h_MissingP[0]->SetTitleSize(0.047, "XY");
  TLine *LineM= new TLine( 0.2, .0, 0.2 ,72.0);
  LineM->SetLineWidth(2);
  LineM->SetLineColor(2);
  h_MissingP[0]->Draw();
  h_MissingP[1]->SetLabelSize(0.053, "XY");
  h_MissingP[1]->SetTitleSize(0.047, "XY");
  h_MissingP[1]->SetFillColor(kRed-7);
  h_MissingP[1]->Draw("same");
  LineM->Draw("same");
  /*
    h_MissingPcut->SetFillColor(kBlue-7);
    h_MissingPcut->Draw("same");
    MMPL->Draw("same");
  */
  c310->SaveAs("imagenes/MissingMomentum.eps");
  
  
  
  TCanvas *c312=new TCanvas("c312","Missing Momentum vs Missing Mass", 1450, 500);
  c312->Divide(2,1);
  c312->cd(1);
  h_MissingPvsMass[0]->SetLabelSize(0.05, "XY");
  h_MissingPvsMass[0]->SetTitleSize(0.045, "XY");
  h_MissingPvsMass[0]->Draw("colz");
  

  c312->cd(2);
  h_MissingPvsMass[1]->SetLabelSize(0.05, "XY");
  h_MissingPvsMass[1]->SetTitleSize(0.045, "XY");
  h_MissingPvsMass[1]->Draw("colz");
  
  c312->SaveAs("imagenes/MissingMomentumCorrelation.eps");


  //-----------------Comparación de sigmas------------------//
  TCanvas *c32=new TCanvas("c32","Invariant mass", 900, 500);
  c32->Divide(2,1);
  c32->cd(1);
  /*
    TLine *LAML1= new TLine( 1.108, 0., 1.108, 4100);
    TLine *LAML2= new TLine( 1.124, 0., 1.124 ,4100);
    TLine *LAML3= new TLine( 1.11, 0., 1.11, 4100);
    TLine *LAML4= new TLine( 1.132, 0., 1.132 ,4100);
    TLine *LAML5= new TLine( 1.092, 0., 1.092, 4100);
    TLine *LAML6= new TLine( 1.14, 0., 1.14 ,4100);
    LAML1->SetLineWidth(2);
    LAML2->SetLineWidth(2);
    LAML3->SetLineWidth(2);
    LAML4->SetLineWidth(2);
    LAML5->SetLineWidth(2);
    LAML6->SetLineWidth(2);
    LAML1->SetLineColor(kGreen);
    LAML2->SetLineColor(kGreen);
    LAML3->SetLineColor(kBlue);
    LAML4->SetLineColor(kBlue);
    LAML5->SetLineColor(kMagenta);
    LAML6->SetLineColor(kMagenta);
  */
  h_LambdaMass->SetLabelSize(0.045, "XY");
  h_LambdaMass->SetTitleSize(0.043, "XY");
  h_LambdaMass->Draw();
  h_LambdaMass->Fit(lamdaMassFit);
  lamdaMassFit->Draw("same");
  // LAML1->Draw("same");
  // LAML2->Draw("same");
  // LAML3->Draw("same");
  // LAML4->Draw("same");
  // LAML5->Draw("same");
  // LAML6->Draw("same");
  NameLinesInv(1.116, 0.002, 12, 4);
  
  c32->cd(2);
  h_InvariantMass->SetLabelSize(0.045, "XY");
  h_InvariantMass->SetTitleSize(0.043, "XY");
  h_InvariantMass->Draw();
  h_InvariantMasscut[0]->SetFillColor(kGreen-7);
  h_InvariantMasscut[1]->SetFillColor(kBlue-7);
  h_InvariantMasscut[2]->SetFillColor(kMagenta-7);
  h_InvariantMasscut[0]->Draw("same");
  h_InvariantMasscut[1]->Draw("same");
  h_InvariantMasscut[2]->Draw("same");
 
  
  c32->SaveAs("imagenes/InvariantMass.eps");


  TCanvas *c40 = new TCanvas("c40","Delta Beta Vs Missing mass and Invariantmass", 900, 500);
  c40->Divide(3,1);
  c40->cd(1);
  h_DeltaBVSInvariantMass->SetLabelSize(0.03, "XY");
  h_DeltaBVSInvariantMass->SetTitleSize(0.02, "XY");
  h_DeltaBVSInvariantMass->Draw("colz");
  c40->cd(2);
  h_DeltaBVSMissingMass->SetLabelSize(0.03, "XY");
  h_DeltaBVSMissingMass->SetTitleSize(0.02, "XY");
  h_DeltaBVSMissingMass->Draw("colz");
  c40->cd(3);
  h_DeltaBVSMissingMomentum->SetLabelSize(0.03, "XY");
  h_DeltaBVSMissingMomentum->SetTitleSize(0.2, "XY");
  h_DeltaBVSMissingMomentum->Draw("colz");

  TCanvas *c41 = new TCanvas("c41","Por si las moscas", 1450, 500);
  c41->Divide(2,1);
  c41->cd(1);
  h_MissingMassvsSigmaMass->Draw("colz");
  c41->cd(2);
  h_MissingPvsSigmaMass->Draw("colz");

  c41->SaveAs("imagenes/CorrelacionesXSIAL.eps");

  TCanvas *c42 = new TCanvas("c42","Theta Kaon Boost", 1450, 500);
  //Principalhisto,Parts,Firstpoint,step,error,Timeofbreak
  vector<double> Points(3);
  HistoBinning(h_KCosThetaCM,2,0.5,0.001,0.5,10);
  c42->cd(1);
  h_KCosThetaCM->Draw();
  c42->SaveAs("imagenes/ThetaKaonBoost.eps");


  TCanvas *c28 = new TCanvas("","Phi distribution to Kaon", 1450, 500);
  c28->Divide(2,2);
  c28->cd(1);
  h_kaonPhiPA[0]->Draw();
  c28->cd(2);
  h_kaonPhiPA[1]->Draw();
  c28->cd(3);
  h_kaonPhiPE[0]->Draw();
  c28->cd(4);
  h_kaonPhiPE[1]->Draw();
  
  c28->SaveAs("imagenes/PhiDistribution.eps");
  
  /*
    TCanvas *c43 = new TCanvas("c43","Phi distribution to Kaon", 1450, 500);
    c43->Divide(2,2);
    c43->cd(1);
    h_KaonPhiCM_PARA[0]->Draw();
    c43->cd(2);
    h_KaonPhiCM_PARA[1]->Draw();
    c43->cd(3);
    h_KaonPhiCM_PERP[0]->Draw();
    c43->cd(4);
    h_KaonPhiCM_PERP[1]->Draw();

    c43->SaveAs("imagenes/Phi1Distribution.eps");
  */
   
  double PPara=0, PPerp=0;
  int iPara=0, iPerp=0;

  TCanvas *c44 = new TCanvas("","Asymmetry", 1450, 500);
  c44->Divide(1,2);

  for(UInt_t i=0; i<MEASGammaP.size(); i++){

    //---------Binning method---------------------
    for(UInt_t j=0; j<MEASGammaP.at(i).size(); j++){
      if(MEASGammaP[i][j]>0){
	PPara+=MEASGammaP[i][j];
	iPara++;
      }
      else {
	PPerp+=abs(MEASGammaP[i][j]);
	iPerp++;
      }
    }
    PPara=PPara/iPara;
    PPerp=PPerp/iPerp;
    FuncAsym->FixParameter(1,PPara/PPerp);
    FuncAsym->FixParameter(2,(PPara+PPerp)/2.0);
    FuncAsym->SetParLimits(3,-1.2,1.2);
    FuncAsym->SetParLimits(0,0.4,2.4);
    c44->cd(i+1);
    h_Asym[i]->Fit("FuncAsym","R");
    h_Asym[i]->Draw();

    //------------MaxLike Method------------
    ROOT::Math::Minimizer* minim = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    MaxLike Min(MEASPhi.at(i), MEASGammaP.at(i));
    ROOT::Math::Functor f(Min,2);
    minim->SetFunction(f);
    minim->SetVariable(0, "Sigma", 0, 0.01);
    minim->SetVariable(1, "Phi", 0, 0.01);
    minim->SetPrintLevel(1);
    minim->Minimize();

    const double *xs = minim->X();
    std::cout << "minim: f(" << xs[0] << "," << xs[1] <<"): "
	      << minim->MinValue()  << std::endl;
  }

  c44->SaveAs("imagenes/Asymmetry.eps");
   
}


//***************** Friend Functions ******************** //

void LinesPTCuts(){

  double x=-145,y1=150,y2=0; //Coordenadas de las líneas

  vector<TLine*> lines(12);
  
  
  for (int i=0; i<11; i+=2) {
    
    lines.at(i)= new TLine(x, y1, x, y2);    
    if (x > 0)
      lines.at(i+1)= new TLine(x+10, y1, x+10, y2);
    
    else
      lines.at(i+1)= new TLine(x-10, y1, x-10, y2);
    
    if (x == -25)
      x+=50;
    else
      x+=60;
    
    lines.at(i)->SetLineWidth(2);
    lines.at(i+1)->SetLineWidth(2);
    lines.at(i)->SetLineColor(2);
    lines.at(i+1)->SetLineColor(2);
    
    lines.at(i)->Draw("same");
    lines.at(i+1)->Draw("same");
  }
}
void NameLinesInv(double Average, double  Sigma, int NSig, int Binning){
  
  string NumSig;
  double Val=2*NSig;
  double y1=50, y2=200;
  TPaveText* ListTPave=NULL;
  while(true) {
    
    NumSig.clear();
    NumSig=std::to_string(-NSig);
    NumSig.insert(NumSig.size(), "#sigma");


    ListTPave = new TPaveText((Average-NSig*Sigma)-0.003,y1,(Average-NSig*Sigma)+0.003,y2,"br");
    ListTPave->SetBorderSize(0);
    ListTPave->SetFillColor(0);
    ListTPave->SetFillStyle(0);
    ListTPave->SetTextFont(42);
    ListTPave->AddText(NumSig.c_str());
    ListTPave->Draw();

    if(-NSig==(Val/2)) break;
    NSig-=Binning;

  }
}

// TPaveState Modified
void TPaveStateModify(TH1 *h1,TH1 *h2){
  gStyle->SetOptStat("me");
  gPad->Update();

  TPaveStats *ps = (TPaveStats*)h1->FindObject("stats");
  ps->SetName("mystats");
  TList *listOfLines = ps->GetListOfLines();
  ostringstream streamObj1;
  streamObj1 << std::fixed;
  streamObj1 << std::setprecision(0);
  streamObj1 << h2->GetEntries();
  string StrObj1 = streamObj1.str();


  ostringstream streamObj2;
  streamObj2 << std::fixed;
  streamObj2 << std::setprecision(3);
  streamObj2 << h2->GetMean();
  string StrObj2 = streamObj2.str();

  string EntriesCut = "Entries = " + StrObj1;
  string MeanCut = "Mean = " + StrObj2;
  TLatex *myt = new TLatex(0,0,EntriesCut.c_str());
  TLatex *MeanCutTex = new TLatex(0,0,MeanCut.c_str());
  myt->SetTextFont(42);
  myt->SetTextSize(0.04);
  myt->SetTextColor(kRed);
  listOfLines->Add(myt);

  MeanCutTex->SetTextFont(42);
  MeanCutTex->SetTextSize(0.04);
  MeanCutTex->SetTextColor(kRed);
  listOfLines->Add(MeanCutTex);

  
  ps->SetTextSize(0.04);
  ps->SetX1NDC(0.6);
  ps->SetY1NDC(0.75);
  h1->SetStats(0);
  gPad->Modified();

  gStyle->SetOptStat("e");
  
}

vector<double> HistoBinning(TH1 *PrincipalHisto,
			    int Parts,
			    double FirstPoint,
			    double step,
			    double error = 0.5,
			    double TimeOfBreak = 5){

  //Define Timers if we need stop the program
  unsigned Time0, Time1;
  Time0 = clock();
  double time = 0;

  //Define histograms and their areas
  TH1F *TempHisto = (TH1F*)PrincipalHisto->Clone("");
  double TheoricArea = PrincipalHisto->Integral("width");   //We compute the area of the principal histogram
  TempHisto->GetXaxis()->SetRangeUser(0,FirstPoint);        //We change the range of the TempHisto
  double NumericArea = TempHisto->Integral("width");        //We compute the area of the temp histogram
  vector<double> CoordinatesX(Parts-1);                           //To save the points the histograms
  double Point = 0;
  
  for(int j=1;j<=Parts-1;j++){
    
    Point=FirstPoint;
    TempHisto->GetXaxis()->SetRangeUser(0,Point);
    NumericArea = TempHisto->Integral("width");
    
    while(true){
      Time1 = clock();
      time = (double(Time1- Time0)/CLOCKS_PER_SEC);

      if(int(time)==TimeOfBreak){
	cout << "ALERT: ¡Stop the code for Time!\n";
	break;
      }
      
      if(NumericArea>TheoricArea*j/Parts){
	Point-=step;
	TempHisto->GetXaxis()->SetRangeUser(0,Point);
	NumericArea = TempHisto->Integral("width");
      }
      else if(NumericArea < TheoricArea*j/Parts){
	Point+=step;
	TempHisto->GetXaxis()->SetRangeUser(0,Point);
	NumericArea = TempHisto->Integral("width");
      }
      if(abs(NumericArea-(TheoricArea*j/Parts)) <= error)
	break;
    }
    cout << "********************************************************************\n"
	 << "\t\tDelta Error: " << abs(NumericArea-(TheoricArea*j/Parts)) << endl
	 << "Point for " << j << " : " << Point << "\t"
	 << "Numeric Area: " << NumericArea << "\t"
	 << "Theoric Area: " << TheoricArea*j/Parts << endl
	 << "Execution Time: " << time << endl
	 << "********************************************************************\n";

    CoordinatesX[j-1]=Point;

   
  }

  return CoordinatesX;
  
}

// Function to asymmetry

Double_t fitf(Double_t *x, Double_t *par){
  double num=0, denom=0;
  num=par[0]-1.0-(par[0]*par[1]+1.0)/(par[1]+1.0)*2*par[2]*par[3]*TMath::Cos(2*x[0]*TMath::DegToRad());
  denom=par[0]+1.0-   (par[0]*par[1]-1.0)/(par[1]+1.0)*2*par[2]*par[3]*TMath::Cos(2*x[0]*TMath::DegToRad());
  Double_t fitval =num/denom;
  return fitval;
}

#endif
