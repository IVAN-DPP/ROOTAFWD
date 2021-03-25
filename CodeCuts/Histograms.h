/*******************************************/
// Author: Ivan Dario Piernagorda Peña     
// Author: Karen Milena Rodriguez      
// Date:   15/02/2020                  
// Title:  Histograms                              
/*******************************************/


#ifndef HISTOGRAMS_C
#define HISTOGRAMS_C

#include "include/Libraries.h"
#include "include/MaxLike.h"
#include "src/TPaveStateModify.C"
#include "src/HistoBinning.C"

using namespace std;

//Friends Functions
void LinesPTCuts();
void NameLinesInv(double, double, int, int);
Double_t fitf(Double_t *,Double_t *);
class Histograms{
  //Friends Functions
  friend void LinesPTCuts();       //Functions for do the lines in Theta-Phi Correlations
  friend void NameLinesInv(double, double, int, int);
  friend Double_t fitf(Double_t *,Double_t *); //Function to asymmetry compute
  
protected:

  TH1F *h_Vertex                            		= NULL;

  TH1F *h_Mass[2]                           		= {};
  
  TH2F *h_DeltaBe[3]                        		= {};
  TH2F *h_BeVSp[3]                          		= {};
  TH2F *h_BeVSpT                            		= NULL;
  TH2F *h_DeltaBecut[3]                     		= {};
  TH2F *h_BeVSpcut[3]                       		= {};
  
  
  TF1 *FFits[3]                             		= {};               //Fits for do cuts
  TF1 *FFitsminus[3]                        		= {};
  string FFname[3]                          		= {};
  
  TH1F *h_DeltaTall[2]                      		= {};
  TH2F *h_DeltaTallvsp[2]                  		= {};
  TH1F *h_DeltaT[2]                         		= {};
  TH1F *h_eloss[3]                          		= {};
  TH2F *h_Celoss[3]                         		= {};
  
  //--> Beta vs P with PDG MASSES  

  TF1 *BeVSp[3]                             		= {};


  //---- Get Coherent Edge ---- //

  TH1F *h_TagrEpho[3]                       		= {};
  
  //---- Missing mass ----//
  
  TH1F *h_MissingMass              		       	= NULL;
  TH1F *h_MissingMasscut                  		= NULL;
  TH1F *h_MissingMass_kaonpioncut         		= NULL;
  
  TH2F *h_MissingMass_vsMissingMasskaonpion[2] 		= {};
  TH2F *h_MissingMass_vsMissingMassPi0[2]       	= {};

  TH1F *h_MissingMass_kaonpion            		= NULL;
  TH1F *h_MissingMass_pi0             			= NULL;
  TH1F *h_MissingMass_Pion				= NULL;
  TH1F *h_MissingMass_Lambda				= NULL;
  TH1F *h_MissingP[2]                      		= {};
  TH1F *h_MissingPcut[2]                   		= {};
  TH2F *h_MissingMvsIMMass[2]                		= {};
  TH2F *h_MissingMassvsSigmaMass          		= NULL;
  TH2F *h_MissingPvsSigmaMass             		= NULL;

  TH1F *h_InvariantMass                   		= NULL;
  TH1F *h_InvariantMasscut[4]                		= {};
  TH1F *h_IMSigmaComparation[3]				= {};
  TH1F *h_MissingMassFinal_Neutron			= NULL;

  //---- Lambda and Lambda Fit ---- //
  
  TH1F *h_LambdaMass                      		= NULL;
  TF1 *lamdaMassFit                       		= NULL; 
  
  
  TH2F *h_DeltaBVSInvariantMass           		= NULL;
  TH2F *h_DeltaBVSMissingMass             		= NULL;
  TH2F *h_DeltaBVSMissingMomentum         		= NULL;

  //------Missing mass Sigma--------//
  TH1F *h_MMassSigma                       		= NULL;
  TH1F *h_MMassSigmaCut                    		= NULL;
    
  //--Correlations between Invariant and missing mass (lambda and sigma)--//
  TH2F *h_InvMassLambda_vsInvMassSigma    		= NULL;
  TH2F *h_MMNeutron_vsMMassSigma[2]       		= {};
  TH2F *h_MMSigmaVS_IMSigma                    		= NULL;
  
  //-----Correlation Momentums (lambda and Sigma)----//
  TH2F *h_CorrelationMMomentum            		= NULL;

  //-----Correlation Beta VS momentum (Neutron)----//
  TH2F *h_BetaVsMomNeu                         		= NULL;

  

  //--- Ellipse Cuts ---- //
  
  TEllipse *myEllipse                     		= NULL;
  double radx=0.01626367, rady=0.02722365, offsetx=1.197115, offsety=1.203825, angle=360*TMath::DegToRad();
  double radx1=0.12744, rady1=0.01959062, offsetx1=1.152367, offsety1=1.199867, angle1=360*TMath::DegToRad();
  //-----Correlation Theta-Phi, ----------//
 
  TH2F *h_ThePhi[3]                         		= {};
  TF1 *F_ThePhiProt[7]                      		= {};
  
  //--- Fiduciary cuts ---//
  TH2F *h_ThePhicut[3]                      		= {};

  
  //-------Momentum proton------------------//

  TH1F *h_MomentumProton                    		= NULL;

    
  //----Costheta-Kaon Boost-------------------//
  TH1F *h_CosThetaCM[3]                      		= {};
  TH1F *h_CosThetaCMT[6]				= {};
  double PARTCOSK[9]					= {};
  vector<vector<double>> PART				= {};
  vector<vector<double>> AvValueG			= {};
  vector<vector<double>> ErrorG				= {};
  double PARTCOSK17[9]					= {};
  TH1F *h_Theta[3] 					= {};
  TH2F *h_CosThetaCorr[3]				= {};
  
  //-----Kaon phi-------------------//
  TH1F *h_kaonPhiPA1[10]                      		= {};
  TH1F *h_kaonPhiPE1[10]                      		= {};
  TH1F *h_kaonPhiPA2[10]                      		= {};
  TH1F *h_kaonPhiPE2[10]                      		= {};
  TH1F *h_kaonPhiPA3[10]                      		= {};
  TH1F *h_kaonPhiPE3[10]                      		= {};
    
  //---------Function to do asymmetry fit----//
  
  vector<vector<double> > MEASPhi{10}; 			//10 is the number of binning
  vector<vector<double> > MEASGammaP{10};

  map<float,vector<vector<double>>> MEASGamma;
  map<float,vector<vector<double>>> MEASPhip;


  TF1 *FuncAsym                             		= NULL;
  TH1F *h_Asym1[10]                           		= {};
  TH1F *h_Asym2[10]                           		= {};
  TH1F *h_Asym3[10]                           		= {};
   
public:
  Histograms(){}
  void DoHistograms();
  void DoCanvas();
  void DoCanvasVertexToDT();
  void DoCanvasFiduciaryCut();
  void DoCanvasEnergyLossAndCoh();
  void DoCanvasReconstruction();
  void DoCanvasReconstructionMM();
  void DoCanvasReconstructionKPi();
  void DoCanvasReconstructionMMCuts();
  void DoCanvasReconstructionBgMissingMomentum();
  void DoCanvasReconstructionIMCorrelations();
  void DoCanvasProtonMomentum();
  void DoCanvasMissingMomentumCorrelations();
  void DoCanvasIMSigma();
  void DoCanvasOthers();
  void DoCanvasCosPart();
  void DoCanvasAsym();
};


void Histograms::DoHistograms(){

  h_Vertex = new TH1F("h_Vertex","; Distance [cm]; Frequency",200,0,-40);

  h_Mass[0] = new TH1F("","; Mass [GeV/c^{2}]; Frequency",200,0,-1.5);
  h_Mass[1] = new TH1F("","; Mass [GeV/c^{2}]; Frequency",200,0,-1.5);
  
  //------------------ Delta Beta ---------------//
  
  h_DeltaBe[0]=new TH2F("h_DeltaBe_0"," ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  h_DeltaBe[1]=new TH2F("h_DeltaBe_1"," ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.15, 0.15);
  h_DeltaBe[2]=new TH2F("h_DeltaBe_2"," ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);

  
  h_BeVSp[0]=new TH2F("h_BeVSp_0"," ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp[1]=new TH2F("h_BeVSp_1"," ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp[2]=new TH2F("h_BeVSp_2"," ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);

  h_BeVSpT = new TH2F("h_BeVSpT"," ; p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);

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

  
  h_DeltaBecut[0] = new TH2F("h_DeltaBe_0cut"," ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  h_DeltaBecut[1] = new TH2F("h_DeltaBe_1cut"," ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.15, 0.15);
  h_DeltaBecut[2] = new TH2F("h_DeltaBe_2cut"," ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  
  h_BeVSpcut[0] = new TH2F("h_BeVSp_0cut"," ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSpcut[1] = new TH2F("h_BeVSp_1cut"," ;p [GeV/c];#beta;",200, 0, 3, 200, 0, 1);
  h_BeVSpcut[2] = new TH2F("h_BeVSp_2cut"," ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);

  //-----Correlation Theta-Phi, ----------//

  h_ThePhi[0] = new TH2F("h_ThePhi_proton"," ; Azimuthal Angle  #phi #circ; Polar Angle  #theta #circ;",200, -180, 180, 200, 0, 150);
  h_ThePhi[1] = new TH2F("h_ThePhi_kaon"," ; Azimuthal Angle  #phi #circ; Polar Angle  #theta #circ;",200,  -180, 180, 200, 0, 150);
  h_ThePhi[2] = new TH2F("h_ThePhi_pion"," ; Azimuthal Angle  #phi #circ; Polar Angle  #theta #circ;", 200,  -180, 180, 200, 0, 150);

  F_ThePhiProt[0] = new TF1("F_ThePhiProt[0]","-1/(0.01*(TMath::Abs(x-180)-28))+10",150,180);
  F_ThePhiProt[1] = new TF1("F_ThePhiProt[0]","-1/(0.01*(TMath::Abs(x-120)-28))+10",90,150);
  F_ThePhiProt[2] = new TF1("F_ThePhiProt[0]","-1/(0.01*(TMath::Abs(x-60)-28))+10",30,90);
  F_ThePhiProt[3] = new TF1("F_ThePhiProt[0]","-1/(0.01*(TMath::Abs(x-(0))-28))+10",-30,30);
  F_ThePhiProt[4] = new TF1("F_ThePhiProt[0]","-1/(0.01*(TMath::Abs(x-(-60))-28))+10",-90,-30);
  F_ThePhiProt[5] = new TF1("F_ThePhiProt[0]","-1/(0.01*(TMath::Abs(x-(-120))-28))+10",-157,-90);
  F_ThePhiProt[6] = new TF1("F_ThePhiProt[0]","-1/(0.01*(TMath::Abs(x-(-180))-28))+10",-180,-143);
  
  //-------------Fiducial cuts-------------------------//

  
  h_ThePhicut[0] = new TH2F("h_ThePhi_protoncut"," ;  Azimuthal Angle  #phi #circ; Polar Angle  #theta #circ;",200, -180, 180, 200, 0, 150);
  h_ThePhicut[1] = new TH2F("h_ThePhi_kaoncut"," ;  Azimuthal Angle  #phi #circ; Polar Angle  #theta #circ;",200,  -180, 180, 200, 0, 150);
  h_ThePhicut[2] = new TH2F("h_ThePhi_pioncut"," ;   Azimuthal Angle  #phi #circ; Polar Angle  #theta #circ;", 200,  -180, 180, 200, 0, 150);


  
  //------------------ Photons, Delta T  ------------------ // 
  
  h_DeltaTall[0]=new TH1F("h_DeltaTall_0"," ;#Delta t [ns];Frequency;", 200, -10, 10);
  h_DeltaTall[1]=new TH1F("h_DeltaTall_1"," ;#Delta t [ns];Frequency;", 200, -10, 10);


  h_DeltaTallvsp[0]=new TH2F("h_DeltaTallvsp_0"," ;p [GeV/c];#Delta t [ns];", 200, 0, 3, 200, -10, 10);
  h_DeltaTallvsp[1]=new TH2F("h_DeltaTallvsp_1"," ;p [GeV/c];#Delta t [ns];", 200, 0, 3, 200, -10, 10);

  //------------Delta T with Cuts ----------- //
  
  h_DeltaT[0]=new TH1F("h_DeltaT_0"," ;#Delta t [ns];Frequency;", 200, -10, 10);
  h_DeltaT[1]=new TH1F("h_DeltaT_1"," ;#Delta t [ns];Frequency;", 200, -10, 10);

  //--------------- Energy loss ----------- //

  h_eloss[0]=new TH1F("h_eloss_0","; (P_{eloss}-P_{meas})/P_{meas};Frequency;",50, 0, 0.05);
  h_eloss[1]=new TH1F("h_eloss_1","; (P_{eloss}-P_{meas})/P_{meas};Frequency;",50, 0, .05);
  h_eloss[2]=new TH1F("h_eloss_2","; (P_{eloss}-P_{meas})/P_{meas};Frequency",50, 0, .05);
  h_Celoss[0]=new TH2F("h_Celoss_0"," ;P_{meas} [GeV/c];(P_{eloss}-P_{meas})/P_{meas};", 200, 0, 2.5, 200, 0, 0.25);
  h_Celoss[1]=new TH2F("h_Celoss_1"," ;P_{meas} [GeV/c];(P_{eloss}-P_{meas})/P_{meas};", 200, 0, 3, 200, 0, 0.15);
  h_Celoss[2]=new TH2F("h_Celoss_2"," ;P_{meas} [GeV/c];(P_{eloss}-P_{meas})/P_{meas};", 200, 0, 2.5, 200, 0, 0.25);

  //---- Get Coherent Edge ---- //

  h_TagrEpho[0] = new TH1F("h_Tagrepho1","",100,0,2000);
  h_TagrEpho[1] = new TH1F("h_Tagrepho2","",100,0,2000);
  h_TagrEpho[2] = new TH1F("h_Tagrepho3","",100,0,2000);
  
  //-------------- Reconstruction --------- //

  h_MissingMass = new TH1F("h_missingmass",
			   "; MM(#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c^{2}]; Frequency",
			   100, 0.7, 1.2);

  // ---- MIS-identification ----/
  
  h_MissingMass_kaonpion = new TH1F("h_missingmass_kaonpion",
				    "; MM(#gamma d #rightarrow #pi^{+} #pi^{-} X p) [GeV/c^{2}]; Frequency",
				    100, 0.7, 1.2);

  h_MissingMass_pi0 = new TH1F("h_missingmass_pi0",
				    "; MM(#gamma d #rightarrow K^{+} #pi^{-}npX) [GeV/c^{2}]; Frequency",
				    100, -0.05, 0.4);
  
  h_MissingMasscut = new TH1F("h_missingmasscut",
			      "; MM(#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c^{2}]; Frequency",
			      100, 0.7, 1.2);
  
  h_MissingMass_kaonpioncut = new TH1F("h_missingmass_kaonpioncut",
				       " ; MM(#gamma d #rightarrow #pi^{+} #pi^{-} X p) [GeV/c^{2}]; Frequency",
				       100, 0.7, 1.2);

   h_MissingMass_Pion = new TH1F(" h_MissingMass_Pion",
				       " ; MM(#gamma d #rightarrow #pi^{+} #pi^{-} X p) [GeV/c^{2}]; Frequency",
				       100, 0.7, 1.2);
 
   h_MissingMass_Lambda  = new TH1F("h_MissingMass_Lambda",
				       " ; MM(#gamma d #rightarrow #pi^{+} #pi^{-} X p) [GeV/c^{2}]; Frequency",
				       100, 0.7, 1.2);

   h_MissingMassFinal_Neutron = new TH1F("h_MissingMassFinal_Neutron",
				       " ; MM(#gamma d #rightarrow #pi^{+} #pi^{-} X p) [GeV/c^{2}]; Frequency",
				       100, 0.7, 1.2);
  

  h_MissingMass_vsMissingMasskaonpion[0] = new TH2F("MissingMass_correlationKaonPion",
						    "; MM(#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c^{2}];  MM(#gamma d #rightarrow #pi^{+} #pi^{-} X p) [GeV/c^{2}]",
						    300, 0.7, 1.2, 300, 0.7, 1.2);

  h_MissingMass_vsMissingMasskaonpion[1] = new TH2F("MissingMass_correlationKaonPion_Cut",
						    "; MM(#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c^{2}];  MM(#gamma d #rightarrow #pi^{+} #pi^{-} X p) [GeV/c^{2}]",
						    300, 0.7, 1.2, 300, 0.7, 1.2);

 
  h_MissingMass_vsMissingMassPi0[0] = new TH2F("MissingMass_correlationPionKaon",
						    "; MM(#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c^{2}];  MM(#gamma d #rightarrow K^{+} #pi^{-}npX) [GeV/c^{2}]",
						      300, 0.7, 1.2, 300, 0, 1.2);
  h_MissingMass_vsMissingMassPi0[1] = new TH2F("h_MissingMass_vsMissingMassPi0_Cut",
						    "; MM(#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c^{2}];  MM(#gamma d #rightarrow K^{+} K^{-} X p) [GeV/c^{2}]",
						      300, 0.7, 1.2, 300, 0, 1.2);

  h_MissingP[0] = new TH1F("h_missingpSigma",
			   "; Missing momentum (#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c]; Frequency",
			   100, 0.0, 1.);
  h_MissingP[1] = new TH1F("h_missingpLambda",
			   "; Missing momentum (#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c]; Frequency",
			   100, 0.0, 1.);

  h_MissingPcut[0] = new TH1F("h_missingpcutSigma",
			      "; Missing momentum (#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c]; Frequency",
			      100, 0.0, 1.5);
  h_MissingPcut[1] = new TH1F("h_missingpcutLambda",
			      "; Missing momentum (#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c]; Frequency",
			      100, 0.0, 1.5);
  
  
  h_MissingMvsIMMass[0] = new TH2F("h_missingMvsmSigma",
				 "; IM(#pi^{-} n) [GeV/c^{2}]; MM(#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c^{2}]",
				 200, 1.0, 1.5, 200, 0.7, 1.2);

  h_MissingMvsIMMass[1] = new TH2F("h_missingMvsmLambda",
				 "; IM(#pi^{-} p) [GeV/c^{2}]; MM(#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c^{2}]",
				 200, 1.08, 1.16, 200, 0.85, 1.1);


  // ------- This aren't important -------- //
  
  h_MissingMassvsSigmaMass = new TH2F("h_missingmvsSigma",
				      "; Masa Invariante (#pi^{-} n) [GeV/c^{2}]; #gamma d #rightarrow K^{+} #pi^{-} X p [GeV/c^{2}]",
				      100,1.0,1.5, 100, 0.7, 1.2);
  
  h_MissingPvsSigmaMass = new TH2F("h_missingPvsSigma",
				   "; Masa [GeV/c^{2}]; p [GeV/c] ",
				   100,1.0,1.5, 100, 0.0, 1.5);
  
  // ----------------------------------- //

  h_InvariantMass = new TH1F("h_InvariantMass",
			     "; IM(#pi^{-} n) [GeV/c^{2}]; Frequency ",
			     100, 1.0, 1.5);
  
  
  h_InvariantMasscut[0] = new TH1F("h_InvariantMasscut3Sig",
				   "; IM(#pi^{-} n) [GeV/c^{2}]; Frequency ",
				   100, 1.0, 1.5);

  h_InvariantMasscut[1] = new TH1F("h_InvariantMasscut4Sig",
				   "; IM(#pi^{-} n) [GeV/c^{2}]; Frequency ",
				   100, 1.0, 1.5);

  h_InvariantMasscut[2] = new TH1F("h_InvariantMasscut5Sig",
				   "; IM(#pi^{-} n) [GeV/c^{2}]; Frequency ",
				   100, 1.0, 1.5);
  
  h_InvariantMasscut[3] = new TH1F("h_InvariantMasscut_ACP",
				   "; IM(#pi^{-} n) [GeV/c^{2}]; Frequency ",
				   100, 1.0, 1.5);

  h_IMSigmaComparation[0]  = new TH1F("h_IMSigmaComparation[0]",
				   "; IM(#pi^{-} n) [GeV/c^{2}]; Frequency ",
				   100, 1.0, 1.5);

  h_IMSigmaComparation[1] = new TH1F("h_IMSigmaComparation[1]",
				   "; IM(#pi^{-} n) [GeV/c^{2}]; Frequency ",
				   100, 1.0, 1.5);

  h_IMSigmaComparation[2] = new TH1F("h_IMSigmaComparation[2]",
				   "; IM(#pi^{-} n) [GeV/c^{2}]; Frequency ",
				   100, 1.0, 1.5);

  //----------Momentum proton---------------//
  
   h_MomentumProton = new TH1F("h_MomentumProton",
			  "; Proton momentum [GeV/c]; Frequency ",
			  100, 0, 2.0);

   //-------- Lambda and Lambda Fit ------ //
  
  h_LambdaMass = new TH1F("h_LambdaMass",
			  "; IM(#pi^{-} p) [GeV/c^{2}]; Frequency ",
			  100, 1.08, 1.16);

  lamdaMassFit = new TF1("lamdaMassFit","gaus",1.08,1.16);
  
  //--------- Others --------------//
  
  h_DeltaBVSInvariantMass = new TH2F("h_DeltaBVSInvariantMass",
				     "",
				     100,0, 2,100,-0.17, 0.17);
  
  h_DeltaBVSMissingMass = new TH2F("h_DeltaBVSMissingMass",
				   "",
				   100,0, 2,100,-0.17, 0.17);
  
  h_DeltaBVSMissingMomentum = new TH2F("h_DeltaBVSMissingMomentum",
				       "",
				       100,0, 1.5,100,-0.17, 0.17);


  //--------------Correlation Invariant masses (Lambda vs Sigma)----//

  h_InvMassLambda_vsInvMassSigma = new TH2F("h_InvMassLambda_vsInvMassSigma",
				      "; IM(#pi^{-} n) [GeV/c^{2}]; IM(#pi^{-} p) [GeV/c^{2}]",
				      200,1.06,1.4, 200, 1.0, 1.5);

  //--------------Correlation MM Sigma and IM Sigma-----------------//
  h_MMSigmaVS_IMSigma = new TH2F("h_MMSigmaVS_IMSigma",
				      "; IM(#pi^{-} p) [GeV/c^{2}]; MM(Sigma) [GeV/c^{2}]",
				      300,1.06,1.4, 300, 1.0, 1.35);
  

  
  //--------------Correlation momentums (Lambda vs Sigma)----//

  h_CorrelationMMomentum  = new TH2F("h_CorrelationMMomentum",
				     "; IM(#pi^{-} n) momentum [GeV/c^{2}]; IM(#pi^{-} p) momentum [GeV/c^{2}] [GeV/c^{2}]",
				     200,0, 1.6 , 200, 0, 1.6);

  // Ensayo momentum vs beta
  h_BetaVsMomNeu = new TH2F(" h_BetaVsMomNeu",
				     "; Missing Momentum (#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c];#Delta #beta_{Neutron}",
			    300, 0, 1.5, 300, -0.15, 0.15);

  

  
  //---------------Missing mass Sigma-----------------------------//
  h_MMassSigma = new TH1F("h_MMassSigma",
				   "; MM(#gamma d #rightarrow K^{+} X p) [GeV/c^{2}]; Frequency ",
				   100, 0.95, 1.45);

  //--------Missing mass Sigma cut--------------//
  h_MMassSigmaCut = new TH1F("h_MMassSigmaCut",
				   "; MM(#gamma d #rightarrow K^{+} X p) [GeV/c^{2}]; Frequency ",
				   100, 0.95, 1.45);
  
  
  //--------------Correlation MM neutron and Missing mass Sigma----//
  
  h_MMNeutron_vsMMassSigma[0] = new TH2F("h_MMNeutron_vsMMassSigma",
					 "; MM(#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c^{2}]; MM(#gamma d #rightarrow K^{+} X p) [GeV/c^{2}]",
					 200, 0.6, 1.5,200,1,1.4);
 


  h_MMNeutron_vsMMassSigma[1] = new TH2F("h_MMNeutron_vsMMassSigmaC",
					 "; MM(#gamma d #rightarrow K^{+} #pi^{-} X p) [GeV/c^{2}]; MM(#gamma d #rightarrow K^{+} X p) [GeV/c^{2}]",
					 200, 0.6, 1.5,200,1,1.4);
 


  //myEllipse = new TEllipse(offsetx,offsety,radx,rady,0,360,70);

  
  //--------------KaonCosTheta Boost-------------//

  h_CosThetaCM[0] = new TH1F("h_CosThetaCM[0]",";cos(#theta^{cm});Frequency", 100, -1, 1);
  h_CosThetaCM[1] = new TH1F("h_CosThetaCM[1]",";cos(#theta^{cm});Frequency", 900, -1, 1);
  h_CosThetaCM[2] = new TH1F("h_CosThetaCM[2]",";cos(#theta^{cm});Frequency", 100, -1, 1);

  h_CosThetaCMT[0] = new TH1F("h_CosThetaCMT[0]",";cos(#theta^{cm});Frequency", 900, -1, 1);
  h_CosThetaCMT[1] = new TH1F("h_CosThetaCMT[1]",";cos(#theta^{cm});Frequency", 900, -1, 1);
  h_CosThetaCMT[2] = new TH1F("h_CosThetaCMT[2]",";cos(#theta^{cm});Frequency", 900, -1, 1);
  h_CosThetaCMT[3] = new TH1F("h_CosThetaCMT[3]",";cos(#theta^{cm});Frequency", 900, -1, 1);
  h_CosThetaCMT[4] = new TH1F("h_CosThetaCMT[4]",";cos(#theta^{cm});Frequency", 900, -1, 1);
  h_CosThetaCMT[5] = new TH1F("h_CosThetaCMT[5]",";cos(#theta^{cm});Frequency", 900, -1, 1);

  PARTCOSK17[0] = -0.7;		  PARTCOSK17[5] = 0.36;
  PARTCOSK17[1] = -0.36;	  PARTCOSK17[6] = 0.44;
  PARTCOSK17[2] = -0.06;	  PARTCOSK17[7] = 0.52;
  PARTCOSK17[3] = 0.12;		  PARTCOSK17[8] = 0.6;
  PARTCOSK17[4] = 0.26;

  
  PARTCOSK[0] = -0.31;		PARTCOSK[5] = 0.53;
  PARTCOSK[1] = 0.009;	  	PARTCOSK[6] = 0.61;
  PARTCOSK[2] = 0.21;	  	PARTCOSK[7] = 0.70;
  PARTCOSK[3] = 0.34;	  	PARTCOSK[8] = 0.76;
  PARTCOSK[4] = 0.44;	
  
  h_Theta[0] = new TH1F("h_ThetaCM[0]","Proton Cos_Theta Boost", 100 ,0 ,3.14159);
  h_Theta[1] = new TH1F("h_ThetaCM[1]","Kaon Cos_Theta Boost", 100 ,0 ,3.14159);
  h_Theta[2] = new TH1F("h_ThetaCM[2]","Sigma Cos_Theta Boost", 100 ,0 ,3.14159);


  h_CosThetaCorr[0] = new TH2F("h_CosThetaCorr[0]","Proton Kaon", 100 ,-1 ,1, 100 ,-1 ,1);
  h_CosThetaCorr[1] = new TH2F("h_CosThetaCorr[1]","Proton Sigma", 100 ,-1 ,1, 100 ,-1 ,1);
  h_CosThetaCorr[2] = new TH2F("h_CosThetaCorr[2]","Sigma Kaon", 100 ,-1 ,1, 100 ,-1 ,1);
  
  //-----------------Kaon phi-------------------//

  vector<vector<float>> XLOW;
  
  for (UInt_t BinsFid = 2; BinsFid <= 6; BinsFid+=2){   	//#of bins/fiducial region. This number has to be edited!! 
    int NSector=6; 						//# of CLAS sectors
    int TotBins=(BinsFid*NSector)+NSector+1; 			//Size of the array xlow{}
    float Width=float(50)/BinsFid; 				//Width of each bin within fiducial region
    int Tot=TotBins;
    vector<float> xlow(Tot); xlow[0]=-180;
    int it=0,i=1,itP = 0;
    if(BinsFid%2 == 0)
      it=0.5*BinsFid+1;
    itP = it;
    while(it < Tot+itP-1){
      if((it%(BinsFid+1)) == 0)					//Coils regions
	xlow[i]=xlow[i-1]+10;
      else
	xlow[i]=xlow[i-1]+Width; 				//Fiducial regions
      it++;
      i++;
    }
    XLOW.push_back(xlow);
  }

  
  h_kaonPhiPA1[0] = new TH1F("h_kaonPhiPA1[0]","First partition Kaon_Phi (PARA)",(XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPA1[1] = new TH1F("h_kaonPhiPA1[1]","Second partition Kaon_Phi (PARA)",(XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPA1[2] = new TH1F("h_kaonPhiPA1[2]","Second partition Kaon_Phi (PARA)",(XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPA1[3] = new TH1F("h_kaonPhiPA1[3]","Second partition Kaon_Phi (PARA)",(XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPA1[4] = new TH1F("h_kaonPhiPA1[4]","Second partition Kaon_Phi (PARA)",(XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPA1[5] = new TH1F("h_kaonPhiPA1[5]","Second partition Kaon_Phi (PARA)",(XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPA1[6] = new TH1F("h_kaonPhiPA1[6]","Second partition Kaon_Phi (PARA)",(XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPA1[7] = new TH1F("h_kaonPhiPA1[7]","Second partition Kaon_Phi (PARA)",(XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPA1[8] = new TH1F("h_kaonPhiPA1[8]","Second partition Kaon_Phi (PARA)",(XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPA1[9] = new TH1F("h_kaonPhiPA1[9]","Second partition Kaon_Phi (PARA)",(XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPE1[0] = new TH1F("h_kaonPhiPE1[0]","First partition Kaon_Phi (PERP)",(XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPE1[1] = new TH1F("h_kaonPhiPE1[1]","Second partition Kaon_Phi (PERP)", (XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPE1[2] = new TH1F("h_kaonPhiPE1[2]","Second partition Kaon_Phi (PERP)", (XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPE1[3] = new TH1F("h_kaonPhiPE1[3]","Second partition Kaon_Phi (PERP)", (XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPE1[4] = new TH1F("h_kaonPhiPE1[4]","Second partition Kaon_Phi (PERP)", (XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPE1[5] = new TH1F("h_kaonPhiPE1[5]","Second partition Kaon_Phi (PERP)", (XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPE1[6] = new TH1F("h_kaonPhiPE1[6]","Second partition Kaon_Phi (PERP)", (XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPE1[7] = new TH1F("h_kaonPhiPE1[7]","Second partition Kaon_Phi (PERP)", (XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPE1[8] = new TH1F("h_kaonPhiPE1[8]","Second partition Kaon_Phi (PERP)", (XLOW[0].size()-1), XLOW[0].data());
  h_kaonPhiPE1[9] = new TH1F("h_kaonPhiPE1[9]","Second partition Kaon_Phi (PERP)", (XLOW[0].size()-1), XLOW[0].data());
  
  
  h_kaonPhiPA2[0] = new TH1F("h_kaonPhiPA2[0]","First partition Kaon_Phi (PARA)",(XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPA2[1] = new TH1F("h_kaonPhiPA2[1]","Second partition Kaon_Phi (PARA)",(XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPA2[2] = new TH1F("h_kaonPhiPA2[2]","Second partition Kaon_Phi (PARA)",(XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPA2[3] = new TH1F("h_kaonPhiPA2[3]","Second partition Kaon_Phi (PARA)",(XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPA2[4] = new TH1F("h_kaonPhiPA2[4]","Second partition Kaon_Phi (PARA)",(XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPA2[5] = new TH1F("h_kaonPhiPA2[5]","Second partition Kaon_Phi (PARA)",(XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPA2[6] = new TH1F("h_kaonPhiPA2[6]","Second partition Kaon_Phi (PARA)",(XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPA2[7] = new TH1F("h_kaonPhiPA2[7]","Second partition Kaon_Phi (PARA)",(XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPA2[8] = new TH1F("h_kaonPhiPA2[8]","Second partition Kaon_Phi (PARA)",(XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPA2[9] = new TH1F("h_kaonPhiPA2[9]","Second partition Kaon_Phi (PARA)",(XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPE2[0] = new TH1F("h_kaonPhiPE2[0]","First partition Kaon_Phi (PERP)",(XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPE2[1] = new TH1F("h_kaonPhiPE2[1]","Second partition Kaon_Phi (PERP)", (XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPE2[2] = new TH1F("h_kaonPhiPE2[2]","Second partition Kaon_Phi (PERP)", (XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPE2[3] = new TH1F("h_kaonPhiPE2[3]","Second partition Kaon_Phi (PERP)", (XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPE2[4] = new TH1F("h_kaonPhiPE2[4]","Second partition Kaon_Phi (PERP)", (XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPE2[5] = new TH1F("h_kaonPhiPE2[5]","Second partition Kaon_Phi (PERP)", (XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPE2[6] = new TH1F("h_kaonPhiPE2[6]","Second partition Kaon_Phi (PERP)", (XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPE2[7] = new TH1F("h_kaonPhiPE2[7]","Second partition Kaon_Phi (PERP)", (XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPE2[8] = new TH1F("h_kaonPhiPE2[8]","Second partition Kaon_Phi (PERP)", (XLOW[1].size()-1), XLOW[1].data());
  h_kaonPhiPE2[9] = new TH1F("h_kaonPhiPE2[9]","Second partition Kaon_Phi (PERP)", (XLOW[1].size()-1), XLOW[1].data());
  
  
  h_kaonPhiPA3[0] = new TH1F("h_kaonPhiPA3[0]","First partition Kaon_Phi (PARA)",(XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPA3[1] = new TH1F("h_kaonPhiPA3[1]","Second partition Kaon_Phi (PARA)",(XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPA3[2] = new TH1F("h_kaonPhiPA3[2]","Second partition Kaon_Phi (PARA)",(XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPA3[3] = new TH1F("h_kaonPhiPA3[3]","Second partition Kaon_Phi (PARA)",(XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPA3[4] = new TH1F("h_kaonPhiPA3[4]","Second partition Kaon_Phi (PARA)",(XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPA3[5] = new TH1F("h_kaonPhiPA3[5]","Second partition Kaon_Phi (PARA)",(XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPA3[6] = new TH1F("h_kaonPhiPA3[6]","Second partition Kaon_Phi (PARA)",(XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPA3[7] = new TH1F("h_kaonPhiPA3[7]","Second partition Kaon_Phi (PARA)",(XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPA3[8] = new TH1F("h_kaonPhiPA3[8]","Second partition Kaon_Phi (PARA)",(XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPA3[9] = new TH1F("h_kaonPhiPA3[9]","Second partition Kaon_Phi (PARA)",(XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPE3[0] = new TH1F("h_kaonPhiPE3[0]","First partition Kaon_Phi (PERP)",(XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPE3[1] = new TH1F("h_kaonPhiPE3[1]","Second partition Kaon_Phi (PERP)", (XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPE3[2] = new TH1F("h_kaonPhiPE3[2]","Second partition Kaon_Phi (PERP)", (XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPE3[3] = new TH1F("h_kaonPhiPE3[3]","Second partition Kaon_Phi (PERP)", (XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPE3[4] = new TH1F("h_kaonPhiPE3[4]","Second partition Kaon_Phi (PERP)", (XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPE3[5] = new TH1F("h_kaonPhiPE3[5]","Second partition Kaon_Phi (PERP)", (XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPE3[6] = new TH1F("h_kaonPhiPE3[6]","Second partition Kaon_Phi (PERP)", (XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPE3[7] = new TH1F("h_kaonPhiPE3[7]","Second partition Kaon_Phi (PERP)", (XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPE3[8] = new TH1F("h_kaonPhiPE3[8]","Second partition Kaon_Phi (PERP)", (XLOW[2].size()-1), XLOW[2].data());
  h_kaonPhiPE3[9] = new TH1F("h_kaonPhiPE3[9]","Second partition Kaon_Phi (PERP)", (XLOW[2].size()-1), XLOW[2].data());
  
  
  //---------------  Fit to Asymmetry ------------//

  FuncAsym = new TF1("FuncAsym",fitf,-180,180,4);
  FuncAsym->SetParameters(1.2,1,0.7,0.5);
  FuncAsym->SetParNames("FR","PR","Pave","Sigma");

  
  //----------------- Histograms Asymmetry -------//

  
  //Asymmetry
  MEASGamma[1.3]=MEASGammaP;		MEASPhip[1.3]=MEASPhi;
  MEASGamma[1.5]=MEASGammaP;		MEASPhip[1.5]=MEASPhi;
  MEASGamma[1.7]=MEASGammaP;		MEASPhip[1.7]=MEASPhi;	
  MEASGamma[1.9]=MEASGammaP;		MEASPhip[1.9]=MEASPhi;
  MEASGamma[2.1]=MEASGammaP;		MEASPhip[2.1]=MEASPhi;
  MEASGamma[2.3]=MEASGammaP;		MEASPhip[2.3]=MEASPhi;
  
}

void Histograms::DoCanvasVertexToDT(){

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

  
  TCanvas *cMt0 = new TCanvas("cMt0","Positive Particles",950,500);
  cMt0->cd(1);
  h_Mass[0]->Draw();
  cMt0->SaveAs("imagenes/MassPos.eps");
  TCanvas *cMt1 = new TCanvas("cMt1","Negative Paricles",950,500);
  cMt1->cd(1);
  h_Mass[1]->Draw();
  cMt1->SaveAs("imagenes/MassNeg.eps");
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

}

void Histograms::DoCanvasFiduciaryCut(){


  //----- Detector geometry (Fiduciary cuts)---------//

  gStyle->SetStatY(0.94);
    
  //------ Proton ------//
  
  TCanvas *c0TP=new TCanvas("c0TP","Theta-Phi correlation", 1200, 500);
  c0TP->cd(1);
  h_ThePhi[0]->SetTitleSize(0.05, "XY");
  h_ThePhi[0]->Draw("colz");
  F_ThePhiProt[0]->Draw("same");
  F_ThePhiProt[1]->Draw("same");
  F_ThePhiProt[2]->Draw("same");
  F_ThePhiProt[3]->Draw("same");
  F_ThePhiProt[4]->Draw("same");
  F_ThePhiProt[5]->Draw("same");
  F_ThePhiProt[6]->Draw("same");

  LinesPTCuts();
  c0TP->SaveAs("imagenes/Fiduciary_Proton.eps");
  
  //------ Kaon ------//
  
  TCanvas *c1TP=new TCanvas("c1TP","Theta-Phi correlation", 1200, 500);
  c1TP->cd(1);
  h_ThePhi[1]->SetTitleSize(0.05, "XY");
  h_ThePhi[1]->Draw("colz");
  F_ThePhiProt[0]->Draw("same");
  F_ThePhiProt[1]->Draw("same");
  F_ThePhiProt[2]->Draw("same");
  F_ThePhiProt[3]->Draw("same");
  F_ThePhiProt[4]->Draw("same");
  F_ThePhiProt[5]->Draw("same");
  F_ThePhiProt[6]->Draw("same");

  LinesPTCuts();
  c1TP->SaveAs("imagenes/Fiduciary_Kaon.eps");

  //------ Pion ------//
  
  TCanvas *c2TP=new TCanvas("c2TP","Theta-Phi correlation", 1200, 500);
  c2TP->cd(1);
  h_ThePhi[2]->SetTitleSize(0.05, "XY");
  h_ThePhi[2]->Draw("colz");
  F_ThePhiProt[0]->Draw("same");
  F_ThePhiProt[1]->Draw("same");
  F_ThePhiProt[2]->Draw("same");
  F_ThePhiProt[3]->Draw("same");
  F_ThePhiProt[4]->Draw("same");
  F_ThePhiProt[5]->Draw("same");
  F_ThePhiProt[6]->Draw("same");


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
  c0TPC->SaveAs("imagenes/Fiduciarycuts_Proton.eps");
  
  //------ Kaon ------//
  
  TCanvas *c1TPC=new TCanvas("c1TPC","Fiduciary cuts", 1200, 500);
  c1TPC->cd(1);
  h_ThePhicut[1]->SetTitleSize(0.05, "XY");
  h_ThePhicut[1]->Draw("colz");
  c1TPC->SaveAs("imagenes/Fiduciarycuts_Kaon.eps");

  //------ Pion ------//
  
  TCanvas *c2TPC=new TCanvas("c2TPC","Fiduciary cuts", 1200, 500);
  c2TPC->cd(1);
  h_ThePhicut[2]->SetTitleSize(0.05, "XY");
  h_ThePhicut[2]->Draw("colz");
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

}

void Histograms::DoCanvasEnergyLossAndCoh(){

  //-------------------- Energy Loss -------------------- //

  gStyle->SetStatY(0.9);

  //------ Proton ------//
  
  TCanvas *c0EL=new TCanvas("c0EL","Delta Energy loss", 1200, 500);
  c0EL->cd(1);
  h_eloss[0]->SetTitleSize(0.04, "XY");
  h_eloss[0]->Draw();
  c0EL->SaveAs("imagenes/Energyloss_Proton.eps");

  TCanvas *c0CEL=new TCanvas("c0CEL","Correlation Energy loss", 1200, 500);
  c0CEL->cd(1);
  h_Celoss[0]->SetTitleSize(0.04, "XY");
  h_Celoss[0]->Draw("colz");
  c0CEL->SaveAs("imagenes/CorrelationEnergyloss_Proton.eps");

 
  //------ Kaon ------//
  
  TCanvas *c1EL=new TCanvas("c1EL","Delta Energy loss", 1200, 500);
  c1EL->cd(1);
  h_eloss[1]->SetTitleSize(0.04, "XY");
  h_eloss[1]->Draw();
  c1EL->SaveAs("imagenes/Energyloss_Kaon.eps");

  TCanvas *c1CEL=new TCanvas("c1CEL","Correlation Energy loss", 1200, 500);
  c1CEL->cd(1);
  h_Celoss[1]->SetTitleSize(0.04, "XY");
  h_Celoss[1]->Draw("colz");
  c1CEL->SaveAs("imagenes/CorrelationEnergyloss_Kaon.eps");

  //------ Pion ------//
  
  TCanvas *c2EL=new TCanvas("c2EL","Delta Energy loss", 1200, 500);
  h_eloss[2]->SetTitleSize(0.04, "XY");
  h_eloss[2]->Draw();
  c2EL->SaveAs("imagenes/Energyloss_Pion.eps");

  TCanvas *c2CEL=new TCanvas("c2CEL","Correlation Energy loss", 1200, 500);
  c2CEL->cd(1);
  h_Celoss[2]->SetTitleSize(0.04, "XY");
  h_Celoss[2]->Draw("colz");
  c2CEL->SaveAs("imagenes/CorrelationEnergyloss_Pion.eps");

  //---- Get Coherent Edge ---- //

  TCanvas *c20 = new TCanvas("c20","Coh Edge", 900,500);
  c20->Divide(3);
  c20->cd(1);
  h_TagrEpho[0]->Draw();
  c20->cd(2);
  h_TagrEpho[1]->Draw();
  c20->cd(3);
  h_TagrEpho[2]->Draw();

}

void Histograms::DoCanvasReconstruction(){
  //******************* Reconstruction **********************//
}

void Histograms::DoCanvasReconstructionMM(){

  //----------- MM before Cuts for MM correlation ------- //

  //------ Kaon ------//
  
  TCanvas *c0MM=new TCanvas("c0MM","Missing mass", 1200, 500);
  c0MM->cd(1);
  h_MissingMass->SetTitleSize(0.043, "XY");
  h_MissingMass->Draw();
  c0MM->SaveAs("imagenes/MissingMass.eps");
}

void Histograms::DoCanvasReconstructionKPi(){

  //---------- MM Kaon-Pion Correlation without cut ---------------//
  
  TCanvas *c0ELL=new TCanvas("c0ELL","Correlation of MM", 900, 500);
  c0ELL->cd(1);
  h_MissingMass_vsMissingMasskaonpion[0]->SetTitleSize(0.045, "XY");  
  h_MissingMass_vsMissingMasskaonpion[0]->Draw("colz");
  TLine *MMCorrLine;
  MMCorrLine = new TLine(0.7, 0.98,1.2, 0.98);
  MMCorrLine->SetLineWidth(2);
  MMCorrLine->SetLineColor(2);
  MMCorrLine->Draw("same");
  
  //myEllipse->SetFillStyle(0);
  //myEllipse->SetLineColor(kRed);
  //myEllipse->Draw("same");
  c0ELL->SaveAs("imagenes/MIS_Identification_KPi.eps");

  //---------- MM Kaon-Pion Correlation with cut ---------------//
  
  TCanvas *c0ELLC=new TCanvas("c0ELLC","Correlation of MM", 900, 500);
  c0ELLC->cd(1);
  h_MissingMass_vsMissingMasskaonpion[1]->SetTitleSize(0.045, "XY");  
  h_MissingMass_vsMissingMasskaonpion[1]->Draw("colz");
  //myEllipse->SetFillStyle(0);
  //myEllipse->SetLineColor(kRed);
  //myEllipse->Draw("same");
  c0ELLC->SaveAs("imagenes/MIS_Identification_KPi_C.eps");
  
  //---------- MM Pion_0 Correlation without cut ---------------//
  
  TCanvas *c2ELL=new TCanvas("c2ELL","Correlation of MM", 900, 500);
  c2ELL->cd(1);
  h_MissingMass_vsMissingMassPi0[0]->SetTitleSize(0.045, "XY");  
  h_MissingMass_vsMissingMassPi0[0]->Draw("colz");
  c2ELL->SaveAs("imagenes/MIS_Identification_Pi0.eps");

  //---------- MM Pion_0 Correlation with cut ---------------//
  
  TCanvas *c2ELLC=new TCanvas("c2ELLC","Correlation of MM", 900, 500);
  c2ELLC->cd(1);
  h_MissingMass_vsMissingMassPi0[1]->SetTitleSize(0.045, "XY");  
  h_MissingMass_vsMissingMassPi0[1]->Draw("colz");
  c2ELLC->SaveAs("imagenes/MIS_Identification_Pi0_C.eps");

  
  gStyle->SetOptStat("me");
   
  TCanvas *c1MMC=new TCanvas("c1MMC","Missing mass", 1200, 500);
  c1MMC->cd(1);
  h_MissingMass_kaonpion->SetLabelSize(0.045, "XY");
  h_MissingMass_kaonpion->SetTitleSize(0.043, "XY");
  h_MissingMass_kaonpion->Draw();
  h_MissingMass_kaonpioncut->SetFillColor(kRed-7);
  TPaveStateModify MMC1(h_MissingMass_kaonpion,h_MissingMass_kaonpioncut);
  MMC1.BoxOptStat("em");
  MMC1.BoxSize(0.7,0.7);
  MMC1.BoxPosition(0.75,(h_MissingMass_kaonpion->GetMaximum()/2),0.85,h_MissingMass_kaonpion->GetMaximum()+100);
  MMC1.BoxTextSize(0.05);
  MMC1.SaveChanges();
  h_MissingMass_kaonpioncut->Draw("same");

  c1MMC->SaveAs("imagenes/MissingMass_Pion.eps");

  TCanvas *c2MM=new TCanvas("c2MM","Missing mass", 1200, 500);
  c2MM->cd(1);
  h_MissingMass_pi0->Draw();
  TLine *MMLinePi0;
  MMLinePi0 = new TLine(0.03, 0.0,0.03, 17500);
  MMLinePi0->SetLineWidth(2);
  MMLinePi0->SetLineColor(2);
  MMLinePi0->Draw("same");
  c2MM->SaveAs("imagenes/MissingMass_Pi0.eps");

}

void Histograms::DoCanvasReconstructionMMCuts(){

  //----------- MM After Cuts for MM correlation ------- //
  

  //------ Kaon ------//  
  

  TCanvas *c0MMC=new TCanvas("c0MMC","Missing mass", 1200, 500);
  c0MMC->cd(1);
  h_MissingMass->SetLabelSize(0.045, "XY");
  h_MissingMass->SetTitleSize(0.043, "XY");
  h_MissingMass->Draw();
  h_MissingMasscut->SetFillColor(kRed-7);
  h_MissingMass_Lambda->SetFillColor(kGreen-7);
  
  vector<TH1*> LHMM;
  LHMM.push_back(h_MissingMasscut);
  LHMM.push_back(h_MissingMass_Lambda);
  gStyle->SetOptStat(110);
  TPaveStateModify MMC0(h_MissingMass, LHMM);
  MMC0.BoxPosition(1.05, h_MissingMass->GetMaximum()/2, 1.2, h_MissingMass->GetMaximum()+100);
  MMC0.BoxTextSize(0.045);
  MMC0.BoxOptStat("em",2);
  MMC0.SaveChanges();
  h_MissingMasscut->Draw("same");
  h_MissingMass_Lambda->Draw("same");
 
  c0MMC->SaveAs("imagenes/MissingMass_Kaon.eps");

  //------Missing mass Neutron after cuts---------------//

  TCanvas *cMMF=new TCanvas("cMMF","Missing mass neutron final", 1200, 500);
  cMMF->cd(1);
  h_MissingMassFinal_Neutron->Draw();
  cMMF->SaveAs("imagenes/MissingMassFinalNeutron.eps");


}

void Histograms::DoCanvasReconstructionBgMissingMomentum(){
  
  //------------ Missing Momentum -----------------//
  
  TCanvas *MMP=new TCanvas("MMP","Missing Momentum", 1450, 500);
  MMP->cd(1);
  TLine *MMPL = new TLine(0.2,0,0.2,4300);
  MMPL->SetLineWidth(2);
  MMPL->SetLineColor(kBlue);
  h_MissingP[0]->SetLabelSize(0.053, "XY");
  h_MissingP[0]->SetTitleSize(0.047, "XY");
  TLine *LineM= new TLine( 0.2, .0, 0.2 ,h_MissingP[0]->GetMaximum());
  LineM->SetLineWidth(2);
  LineM->SetLineColor(2);
  h_MissingP[0]->Draw();
  gStyle->SetOptStat("e");
  h_MissingP[1]->SetLabelSize(0.053, "XY");
  h_MissingP[1]->SetTitleSize(0.047, "XY");
  h_MissingP[1]->SetFillColor(kRed-7);
  TPaveStateModify MMPStat(h_MissingP[0],h_MissingP[1]);
  MMPStat.BoxOptStat("e");
  MMPStat.BoxPosition(0.75,0.85*h_MissingP[0]->GetMaximum(),1,h_MissingP[0]->GetMaximum());
  MMPStat.BoxTextSize(0.04);
  MMPStat.SaveChanges();
  h_MissingP[1]->Draw("same");
  LineM->Draw("same");

  MMP->SaveAs("imagenes/MissingMomentum.eps");

}

void Histograms::DoCanvasReconstructionIMCorrelations(){
  TCanvas *c312=new TCanvas("c312","MM Neutron vs IM Sigma", 1450, 500);
  c312->cd(1);
  h_MissingMvsIMMass[0]->SetLabelSize(0.05, "XY");
  h_MissingMvsIMMass[0]->SetTitleSize(0.045, "XY");
  h_MissingMvsIMMass[0]->Draw("colz");
  c312->SaveAs("imagenes/MMvsInvMassCorrelationSigma.eps");

  TCanvas *c314=new TCanvas("c314","MM Neutron vs IM Lambda", 1450, 500);
  c314->cd(1);
  h_MissingMvsIMMass[1]->SetLabelSize(0.05, "XY");
  h_MissingMvsIMMass[1]->SetTitleSize(0.045, "XY");
  h_MissingMvsIMMass[1]->Draw("colz");
   
  c314->SaveAs("imagenes/MMvsInvMassCorrelationLambda.eps");


  //-----------------Std Desviation Comparation------------------//

  //------ Lambda -------//
  
  TCanvas *STDC = new TCanvas("STDC","Invariant mass", 1450, 500);
  STDC->cd(1);
  h_LambdaMass->SetLabelSize(0.045, "XY");
  h_LambdaMass->SetTitleSize(0.043, "XY");
  h_LambdaMass->Draw();
  TLine *LambdaLines[2]={};
  LambdaLines[0] = new TLine(1.1, h_LambdaMass->GetMaximum(),1.1,0);
  LambdaLines[1] = new TLine(1.132, h_LambdaMass->GetMaximum(),1.132,0);

  LambdaLines[0]->SetLineWidth(2);
  LambdaLines[1]->SetLineWidth(2);
  LambdaLines[0]->SetLineColor(3);
  LambdaLines[1]->SetLineColor(3);

  LambdaLines[0]->Draw("same");
  LambdaLines[1]->Draw("same");
  
  h_LambdaMass->Fit(lamdaMassFit);
  gStyle->SetOptFit(111);
  lamdaMassFit->Draw("same");
  NameLinesInv(1.116, 0.002, 12, 4);

  STDC->SaveAs("imagenes/InvariantMassComparation_Lambda.eps");

  //----- Sigma ------//

  TCanvas *IVM = new TCanvas("IVM","Invariant mass", 900, 500);
  IVM->cd(1);
  h_InvariantMass->SetLabelSize(0.045, "XY");
  h_InvariantMass->SetTitleSize(0.043, "XY");
  h_InvariantMass->Draw();
  gStyle->SetOptStat("me");
  h_InvariantMasscut[0]->SetFillColor(kGreen-7);
  h_InvariantMasscut[1]->SetFillColor(kBlue-7);
  h_InvariantMasscut[2]->SetFillColor(kMagenta-7);
  vector<TH1*> IVMHistos;
  IVMHistos.push_back(h_InvariantMasscut[0]);
  IVMHistos.push_back(h_InvariantMasscut[1]);
  IVMHistos.push_back(h_InvariantMasscut[2]);
  TPaveStateModify IVMStat(h_InvariantMass,IVMHistos);
  IVMStat.BoxOptStat("em", 2);
  IVMStat.BoxTextSize(0.04);
  IVMStat.BoxPosition(1.35,0.4*h_InvariantMass->GetMaximum(),1.5,h_InvariantMass->GetMaximum());
  IVMStat.SaveChanges();
  h_InvariantMasscut[0]->Draw("same");
  h_InvariantMasscut[1]->Draw("same");
  h_InvariantMasscut[2]->Draw("same");
 
  
  IVM->SaveAs("imagenes/InvariantMassComparation_Sigma.eps");

  //-------Sigma Comparation------------------------------//
  
  TCanvas *IVMC = new TCanvas("IVMC","Invariant mass comparation", 900, 500);
  IVMC->cd(1);
  h_IMSigmaComparation[0]->SetLabelSize(0.045, "XY");
  h_IMSigmaComparation[0]->SetTitleSize(0.043, "XY");
  h_IMSigmaComparation[0]->Draw();
  gStyle->SetOptStat("me");
  h_IMSigmaComparation[1]->SetFillColor(kRed-7);
  h_IMSigmaComparation[2]->SetFillColor(kGreen-7);
  vector<TH1*> IVMCHistos;
  IVMCHistos.push_back(h_IMSigmaComparation[1]);
  IVMCHistos.push_back(h_IMSigmaComparation[2]);
  TPaveStateModify IVMCStat(h_IMSigmaComparation[0], IVMCHistos);
  IVMCStat.BoxOptStat("em", 2);
  IVMCStat.BoxTextSize(0.04);
  IVMCStat.BoxPosition(1.35,0.4*h_IMSigmaComparation[0]->GetMaximum(),1.5,h_IMSigmaComparation[0]->GetMaximum());
  IVMCStat.SaveChanges();
  h_IMSigmaComparation[1]->Draw("same");
  h_IMSigmaComparation[2]->Draw("same");
 
  IVMC->SaveAs("imagenes/IMSigmaComparation.eps");

  //--------------Correlation Invariant masses (Lambda vs Sigma)----//
  TCanvas *cIMLS=new TCanvas("cIMLS","Correlation of Invariant masses", 900, 500);
  cIMLS->cd(1);
  h_InvMassLambda_vsInvMassSigma->SetTitleSize(0.045, "XY");  
  h_InvMassLambda_vsInvMassSigma->Draw("colz");
  TLine *IMCorrLine1;
  IMCorrLine1 = new TLine(1,1.1,1.4, 1.1);
  IMCorrLine1->SetLineWidth(2);
  IMCorrLine1->SetLineColor(2);
  IMCorrLine1->Draw("same");
  TLine *IMCorrLine2;
  IMCorrLine2 = new TLine(1,1.132,1.4, 1.132);
  IMCorrLine2->SetLineWidth(2);
  IMCorrLine2->SetLineColor(2);
  IMCorrLine2->Draw("same");
  
  cIMLS->SaveAs("imagenes/InvariantMassCorrelation.eps");

  //-------------------Correlation MM Sigma and IM Sigma------------//
  TCanvas *cMMIMS=new TCanvas("MMIMS","Correlation MM Sigma and IM Sigma", 900, 500);
  cMMIMS->cd(1);
  h_MMSigmaVS_IMSigma->SetTitleSize(0.045, "XY");  
  h_MMSigmaVS_IMSigma->Draw("colz");
  
  cMMIMS->SaveAs("imagenes/MissingMassSigmaVSInvariantMassSigma.eps");

  //--------------Correlation Invariant mass (Lambda) vs Missing Mass Sigma----//
  TCanvas *cIMMM=new TCanvas("cIMMM","Correlation Invariant mass (lambda) vs Missing Mass Sigma", 900, 500);
  cIMMM->cd(1);
  h_MMNeutron_vsMMassSigma[0]->SetTitleSize(0.045, "XY");  
  h_MMNeutron_vsMMassSigma[0]->Draw("colz");
  
  cIMMM->SaveAs("imagenes/InvariantMassMMSigmaCorrelation.eps");


  TCanvas *cIMMMC=new TCanvas("cIMMMC","Correlation Invariant mass (lambda) vs Missing Mass Sigma", 900, 500);
  cIMMMC->cd(1);
  h_MMNeutron_vsMMassSigma[1]->SetTitleSize(0.045, "XY");  
  h_MMNeutron_vsMMassSigma[1]->Draw("colz");
  
  cIMMMC->SaveAs("imagenes/InvariantMassMMSigmaCorrelation_C.eps");

}

void Histograms::DoCanvasProtonMomentum(){

  //-------------Momentum proton--------------//

  TCanvas *cMP=new TCanvas("cMP","Momentum proton", 900, 500);
  cMP->cd(1);
  h_MomentumProton->SetTitleSize(0.045, "XY");  
  h_MomentumProton->Draw();
  
  cMP->SaveAs("imagenes/MomentumProton.eps");    

}

void Histograms::DoCanvasMissingMomentumCorrelations(){

  //--------------Correlation Missing momentums--------------------//

  TCanvas *cMMom=new TCanvas("cMMom","Correlation Missing momentums", 900, 500);
  cMMom->cd(1);
  h_CorrelationMMomentum->SetTitleSize(0.045, "XY");  
  h_CorrelationMMomentum->Draw("colz");
  
  cMMom->SaveAs("imagenes/CorrelationMissingMomentums.eps");

  // ensayo momentum vs beta
  TCanvas *cMMomN=new TCanvas("cMMomN","Correlation Momentum", 900, 500);
  cMMomN->cd(1);
  h_BetaVsMomNeu->SetTitleSize(0.045, "XY");  
  h_BetaVsMomNeu->Draw("colz");
  TLine *CutDBLines[2]={};
  CutDBLines[0] = new TLine(h_BetaVsMomNeu->GetMinimum(), 0.005,1.5, 0.005);
  CutDBLines[1] = new TLine(h_BetaVsMomNeu->GetMinimum(), -0.005,1.5, -0.005);

  CutDBLines[0]->SetLineWidth(2);
  CutDBLines[1]->SetLineWidth(2);
  CutDBLines[0]->SetLineColor(2);
  CutDBLines[1]->SetLineColor(2);

  CutDBLines[0]->Draw("same");
  CutDBLines[1]->Draw("same");
  cMMomN->SaveAs("imagenes/MomentumVsBetaNeutron.eps");

}

void Histograms::DoCanvasIMSigma(){

  
  //---------------Missing mass Sigma-----------------------//
  
  TCanvas *cMMS=new TCanvas("cMMS","Missing mass Sigma", 900, 500);
  cMMS->cd(1);
  h_MMassSigma->SetTitleSize(0.045, "XY");  
  h_MMassSigma->Draw();
  
  cMMS->SaveAs("imagenes/MissingMassSigma.eps");

  //---------------Missing mass Sigma Cut-----------------------//
  
  TCanvas *cMMSC=new TCanvas("cMMSC","Missing mass Sigma Cut", 900, 500);
  cMMSC->cd(1);
  h_MMassSigmaCut->SetTitleSize(0.045, "XY");  
  h_MMassSigmaCut->Draw();
  
  cMMSC->SaveAs("imagenes/MissingMassSigmaCut.eps");

  //------------ Final Invarian Mass ---------------//
  
  TCanvas *IVMF = new TCanvas("IVMF","Invariant mass", 900, 500);
  IVMF->cd(1);

  gStyle->SetOptStat("me");
  h_InvariantMasscut[3]->Draw("same");
 
  
  IVMF->SaveAs("imagenes/InvariantMassFinall_Sigma.eps");


}
void Histograms::DoCanvasOthers(){

  //------------- Others ---------------- //
  
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

}

void Histograms::DoCanvas(){  
  DoCanvasVertexToDT();
  DoCanvasFiduciaryCut();
  DoCanvasEnergyLossAndCoh();
  DoCanvasReconstruction();
  DoCanvasReconstructionMM();
  DoCanvasReconstructionKPi();
  DoCanvasReconstructionMMCuts();
  DoCanvasReconstructionBgMissingMomentum();
  DoCanvasReconstructionIMCorrelations();
  DoCanvasProtonMomentum();
  DoCanvasMissingMomentumCorrelations();
  DoCanvasIMSigma();
  DoCanvasOthers();
}


void Histograms::DoCanvasCosPart(){

  for (UInt_t i = 0; i < 6; i++) {
    HistoBinning *BT = new HistoBinning(h_CosThetaCMT[i],10);
    vector<double> AvValue, Error;
    BT->DoHistoBinning();
    BT->GetPoints(AvValue, Error);
    AvValueG.push_back(AvValue);
    ErrorG.push_back(Error);
    PART.push_back(BT->GetPartitions());
    PART[i].pop_back();
    PART[i].erase(PART[i].begin(),PART[i].begin()+1);
  }

}

void Histograms::DoCanvasAsym(){
  
  TCanvas *cCM0 = new TCanvas("cCM0","cos Theta proton Boost", 1450, 500);
  cCM0->cd(1);
  h_CosThetaCM[0]->Draw();
  cCM0->SaveAs("imagenes/ThetaProtonBoost.eps");

  gStyle->SetOptStat("me");
  
  TCanvas *cCM1 = new TCanvas("cCM1","cos Theta Kaon Boost", 1450, 500);
  HistoBinning *BB = new HistoBinning(h_CosThetaCM[1],10);
  BB->DoHistoBinning();
  BB->PrintLevel(2);
  BB->GetLatexTable("./CosBin.tex", "Cos Theta Binning", "COST");
  vector<double> AvValue, Error;
  BB->GetPoints(AvValue, Error);
  cCM1->cd(1);
  h_CosThetaCM[1]->Draw();

  TLine *COSL[9];	
  double MaxGraphCosK = h_CosThetaCM[1]->GetMaximum();
  
  COSL[0] = new TLine(PARTCOSK[0],0, PARTCOSK[0], MaxGraphCosK);
  COSL[1] = new TLine(PARTCOSK[1],0, PARTCOSK[1], MaxGraphCosK);
  COSL[2] = new TLine(PARTCOSK[2],0, PARTCOSK[2], MaxGraphCosK);
  COSL[3] = new TLine(PARTCOSK[3],0, PARTCOSK[3], MaxGraphCosK);
  COSL[4] = new TLine(PARTCOSK[4],0, PARTCOSK[4], MaxGraphCosK);
  COSL[5] = new TLine(PARTCOSK[5],0, PARTCOSK[5], MaxGraphCosK);
  COSL[6] = new TLine(PARTCOSK[6],0, PARTCOSK[6], MaxGraphCosK);
  COSL[7] = new TLine(PARTCOSK[7],0, PARTCOSK[7], MaxGraphCosK);
  COSL[8] = new TLine(PARTCOSK[8],0, PARTCOSK[8], MaxGraphCosK);
  
  COSL[0]->Draw("same");
  COSL[1]->Draw("same");
  COSL[2]->Draw("same");
  COSL[3]->Draw("same");
  COSL[4]->Draw("same");
  COSL[5]->Draw("same");
  COSL[6]->Draw("same");
  COSL[7]->Draw("same");
  COSL[8]->Draw("same");

  TText *TextCosK[10];

  TextCosK[0] = new TText(AvValue[0],MaxGraphCosK/2,"Bin 1");
  TextCosK[1] = new TText(AvValue[1],MaxGraphCosK/2,"Bin 2");
  TextCosK[2] = new TText(AvValue[2],MaxGraphCosK/2,"Bin 3");
  TextCosK[3] = new TText(AvValue[3],MaxGraphCosK/2,"Bin 4");
  TextCosK[4] = new TText(AvValue[4],MaxGraphCosK/2,"Bin 5");
  TextCosK[5] = new TText(AvValue[5],MaxGraphCosK/2,"Bin 6");
  TextCosK[6] = new TText(AvValue[6],MaxGraphCosK/2,"Bin 7");
  TextCosK[7] = new TText(AvValue[7],MaxGraphCosK/2,"Bin 8");
  TextCosK[8] = new TText(AvValue[8],MaxGraphCosK/2,"Bin 9");
  TextCosK[9] = new TText(AvValue[9],MaxGraphCosK/2,"Bin 10");

  TextCosK[0]->SetTextAngle(90);	  TextCosK[5]->SetTextAngle(90);
  TextCosK[1]->SetTextAngle(90);	  TextCosK[6]->SetTextAngle(90);
  TextCosK[2]->SetTextAngle(90);	  TextCosK[7]->SetTextAngle(90);
  TextCosK[3]->SetTextAngle(90);	  TextCosK[8]->SetTextAngle(90);
  TextCosK[4]->SetTextAngle(90);	  TextCosK[9]->SetTextAngle(90);
	  

  TextCosK[0]->Draw("same");	  TextCosK[5]->Draw("same");
  TextCosK[1]->Draw("same");	  TextCosK[6]->Draw("same");
  TextCosK[2]->Draw("same");	  TextCosK[7]->Draw("same");
  TextCosK[3]->Draw("same");	  TextCosK[8]->Draw("same");
  TextCosK[4]->Draw("same");	  TextCosK[9]->Draw("same");

  cCM1->SaveAs("imagenes/ThetaKaonBoost.eps");
  
  TCanvas *cCM2 = new TCanvas("cCM2","cos Theta Sigma Boost", 1450, 500);
  cCM2->cd(1);
  h_CosThetaCM[2]->Draw();
  cCM2->SaveAs("imagenes/ThetaSigmaBoost.eps");


  TCanvas *cCMT = new TCanvas("cCMT","cos Theta Kaon Boost", 1450, 500);
  cCMT->Divide(2,3);
  for (UInt_t i = 0; i < 6; i++) {
    cCMT->cd(i+1);
    h_CosThetaCMT[i]->Draw();
    
    TLine *COSLT[9];	
    double MaxGraphCosKT = h_CosThetaCMT[i]->GetMaximum();
  
    COSLT[0] = new TLine(PART[i][0],0, PART[i][0], MaxGraphCosKT);
    COSLT[1] = new TLine(PART[i][1],0, PART[i][1], MaxGraphCosKT);
    COSLT[2] = new TLine(PART[i][2],0, PART[i][2], MaxGraphCosKT);
    COSLT[3] = new TLine(PART[i][3],0, PART[i][3], MaxGraphCosKT);
    COSLT[4] = new TLine(PART[i][4],0, PART[i][4], MaxGraphCosKT);
    COSLT[5] = new TLine(PART[i][5],0, PART[i][5], MaxGraphCosKT);
    COSLT[6] = new TLine(PART[i][6],0, PART[i][6], MaxGraphCosKT);
    COSLT[7] = new TLine(PART[i][7],0, PART[i][7], MaxGraphCosKT);
    COSLT[8] = new TLine(PART[i][8],0, PART[i][8], MaxGraphCosKT);
  
    COSLT[0]->Draw("same");
    COSLT[1]->Draw("same");
    COSLT[2]->Draw("same");
    COSLT[3]->Draw("same");
    COSLT[4]->Draw("same");
    COSLT[5]->Draw("same");
    COSLT[6]->Draw("same");
    COSLT[7]->Draw("same");
    COSLT[8]->Draw("same");

    TText *TextCosKT[10];

    TextCosKT[0] = new TText(AvValueG[i][0],MaxGraphCosKT/2,"Bin 1");
    TextCosKT[1] = new TText(AvValueG[i][1],MaxGraphCosKT/2,"Bin 2");
    TextCosKT[2] = new TText(AvValueG[i][2],MaxGraphCosKT/2,"Bin 3");
    TextCosKT[3] = new TText(AvValueG[i][3],MaxGraphCosKT/2,"Bin 4");
    TextCosKT[4] = new TText(AvValueG[i][4],MaxGraphCosKT/2,"Bin 5");
    TextCosKT[5] = new TText(AvValueG[i][5],MaxGraphCosKT/2,"Bin 6");
    TextCosKT[6] = new TText(AvValueG[i][6],MaxGraphCosKT/2,"Bin 7");
    TextCosKT[7] = new TText(AvValueG[i][7],MaxGraphCosKT/2,"Bin 8");
    TextCosKT[8] = new TText(AvValueG[i][8],MaxGraphCosKT/2,"Bin 9");
    TextCosKT[9] = new TText(AvValueG[i][9],MaxGraphCosKT/2,"Bin 10");

    TextCosKT[0]->SetTextAngle(90);	  TextCosKT[5]->SetTextAngle(90);
    TextCosKT[1]->SetTextAngle(90);	  TextCosKT[6]->SetTextAngle(90);
    TextCosKT[2]->SetTextAngle(90);	  TextCosKT[7]->SetTextAngle(90);
    TextCosKT[3]->SetTextAngle(90);	  TextCosKT[8]->SetTextAngle(90);
    TextCosKT[4]->SetTextAngle(90);	  TextCosKT[9]->SetTextAngle(90);
	  

    TextCosKT[0]->Draw("same");	  TextCosKT[5]->Draw("same");
    TextCosKT[1]->Draw("same");	  TextCosKT[6]->Draw("same");
    TextCosKT[2]->Draw("same");	  TextCosKT[7]->Draw("same");
    TextCosKT[3]->Draw("same");	  TextCosKT[8]->Draw("same");
    TextCosKT[4]->Draw("same");	  TextCosKT[9]->Draw("same");
  }
  cCMT->SaveAs("imagenes/ThetaKaonBinBoost.eps");
  
  TCanvas *cTH0 = new TCanvas("CTH0","Theta proton",1450,500);
  cTH0->Divide(1,3);
  cTH0->cd(1);
  h_Theta[0]->Draw();
  cTH0->cd(2);
  h_Theta[1]->Draw();
  cTH0->cd(3);
  h_Theta[2]->Draw();

  cTH0->SaveAs("imagenes/Theta.eps");

  TCanvas *cTHC0 = new TCanvas("CTHC0","Theta proton",1450,500);
  cTHC0->Divide(3,1);
  cTHC0->cd(1);
  h_CosThetaCorr[0]->Draw("colz");
  cTHC0->cd(2);
  h_CosThetaCorr[1]->Draw("colz");
  cTHC0->cd(3);
  h_CosThetaCorr[2]->Draw("colz");

  cTHC0->SaveAs("imagenes/ThetaCorr.eps");

  TCanvas *cPH0 = new TCanvas("","Phi PARA distribution to Kaon", 1450, 500);
  cPH0->Divide(5,2);
  for(int i = 0;i<10;i++){
    cPH0->cd(i+1);
    h_kaonPhiPA1[i]->Draw();
  }
  cPH0->SaveAs("imagenes/PhiDistributionPARA.eps");

  TCanvas *cPH1 = new TCanvas("","Phi PERP distribution to Kaon", 1450, 500);
  cPH1->Divide(5,2);
  for (int i = 0; i < 10; i++) {
    cPH1->cd(i+1);
    h_kaonPhiPE1[i]->Draw();
  }
  
  cPH1->SaveAs("imagenes/PhiDistributionPERP.eps");
  
  //--------- Binning method ----------------- //
  
  vector<double> Asym1, AsymE1;
  vector<double> Asym2, AsymE2;
  vector<double> Asym3, AsymE3;
      
  TCanvas *CanvasAsym[3]	= {};
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat("me");
  
  for (UInt_t k = 0; k < 3; k++) {				// # Of binning 2,4 and 6
    CanvasAsym[k] = new TCanvas("","Asymmetry", 1450, 500);
    CanvasAsym[k]->Divide(3,4);
    for(UInt_t i=0; i<10; i++){
      
      double PPara=0, PPerp=0;
      int iPara=0, iPerp=0;
      
      for(UInt_t j=0; j<MEASGamma[2.3].at(i).size(); j++){
	if(MEASGamma[2.3][i][j]>0){
	  PPara+=MEASGamma[2.3][i][j];
	  iPara++;
	}
	else {
	  PPerp+=abs(MEASGamma[2.3][i][j]);
	  iPerp++;
	}
      }
      PPara=PPara/iPara;
      PPerp=PPerp/iPerp;
      FuncAsym->FixParameter(1,PPara/PPerp);
      FuncAsym->FixParameter(2,(PPara+PPerp)/2.0);
      FuncAsym->SetParLimits(3,-1.2,1.2);
      FuncAsym->SetParLimits(0,0.4,2.4);
      if(i == 9) CanvasAsym[k]->cd(i+2);
      else	 CanvasAsym[k]->cd(i+1);
      if(k == 0){
	h_Asym1[i]->Draw();
	h_Asym1[i]->Fit(FuncAsym);
	Asym1.push_back(FuncAsym->GetParameter(3)*25.0*TMath::DegToRad()/sin(25.0*TMath::DegToRad()));
	AsymE1.push_back(FuncAsym->GetParError(3));
      }
      else if(k == 1){
	h_Asym2[i]->Draw();
	h_Asym2[i]->Fit(FuncAsym);
	Asym2.push_back(FuncAsym->GetParameter(3)*12.5*TMath::DegToRad()/sin(12.5*TMath::DegToRad()));
	AsymE2.push_back(FuncAsym->GetParError(3));
      }
      else {
	h_Asym3[i]->Draw();
	h_Asym3[i]->Fit(FuncAsym);
	Asym3.push_back(FuncAsym->GetParameter(3)*8.33*TMath::DegToRad()/sin(8.33*TMath::DegToRad()));
	AsymE3.push_back(FuncAsym->GetParError(3));
      }

      h_Asym1[i]->GetXaxis()->SetTitle("Azimuthal Angle  #phi #circ (Lab)");
      h_Asym1[i]->GetYaxis()->SetTitle("(N_{PARA}-N_{PERP})/(N_{PARA}+N_{PERP})");
      h_Asym1[i]->SetTitle(" ");
      h_Asym2[i]->GetXaxis()->SetTitle("Azimuthal Angle  #phi #circ (Lab)");
      h_Asym2[i]->GetYaxis()->SetTitle("(N_{PARA}-N_{PERP})/(N_{PARA}+N_{PERP})");
      h_Asym2[i]->SetTitle(" ");
      h_Asym3[i]->GetXaxis()->SetTitle("Azimuthal Angle  #phi #circ (Lab)");
      h_Asym3[i]->GetYaxis()->SetTitle("(N_{PARA}-N_{PERP})/(N_{PARA}+N_{PERP})");
      h_Asym3[i]->SetTitle(" ");
    }
    string AsymS;
    AsymS = "imagenes/Asymmetry" + std::to_string(k) + ".eps";
    CanvasAsym[k]->SaveAs(AsymS.c_str());
  } 

  //------------MaxLike Method------------//

  vector<vector<double>> MaxL(6), MaxLE(6);
  UInt_t it = 0;
  for(double k=1.3;k<=2.5;k+=0.2){
    cout << "//------Energy------// : " << k <<"\n";
    for (UInt_t i = 0; i < MEASPhip[float(k)].size(); i++) {
      ROOT::Math::Minimizer* minim = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
      MaxLike Min(MEASPhip[float(k)].at(i),MEASGamma[float(k)].at(i));
      ROOT::Math::Functor f(Min,1);
      minim->SetFunction(f);
      minim->SetVariable(0, "Sigma", 0, 0.01);
      minim->SetPrintLevel(1);
      minim->Minimize();
      const double *xs = minim->X();
      const double *ee = minim->Errors();
      MaxL[it].push_back(xs[0]);
      MaxLE[it].push_back(ee[0]);
    }
    it++;
  }

  // -------- Delta Phi Binning Comparation ------ //
  
  TCanvas *cASMB = new TCanvas("","",900,450);
  cASMB->cd(1);
  TMultiGraph *AsymBB = new TMultiGraph();
  TGraphErrors *AsymG[3];
  AsymG[0] = new TGraphErrors(Asym1.size(),&(AvValueG[5][0]),&(Asym1[0]),&(ErrorG[5][0]),&(AsymE1[0]));
  AsymG[1] = new TGraphErrors(Asym2.size(),&(AvValueG[5][0]),&(Asym2[0]),&(ErrorG[5][0]),&(AsymE2[0]));
  AsymG[2] = new TGraphErrors(Asym3.size(),&(AvValueG[5][0]),&(Asym3[0]),&(ErrorG[5][0]),&(AsymE3[0]));
  AsymG[0]->SetTitle("#Delta #phi = 25.0^{#circ}"); AsymG[0]->SetLineColor(kBlue);
  AsymG[1]->SetTitle("#Delta #phi = 12.5^{#circ}"); AsymG[1]->SetLineColor(kRed);
  AsymG[2]->SetTitle("#Delta #phi = 8.3^{#circ}"); AsymG[2]->SetLineColor(kGreen);
  AsymBB->Add(AsymG[0]);
  AsymBB->Add(AsymG[1]);
  AsymBB->Add(AsymG[2]);
  AsymBB->GetHistogram()->GetYaxis()->SetRangeUser(-1,1);
  AsymBB->Draw("AP");
  AsymBB->GetXaxis()->SetTitle("cos(#theta^{cm}_{K^{+}})");
  AsymBB->GetYaxis()->SetTitle("#Sigma");
  cASMB->BuildLegend();
  cASMB->SaveAs("imagenes/SigmaAsymBin.eps");


  // -------- Delta Phi Binning Systematic Errors ------ //
  
  TCanvas *cASMBSyE[3];
  TGraphErrors *BinSysError[3];
  vector<double> SysBin1, SysBinErrors1;
  vector<double> SysBin2, SysBinErrors2;
  vector<double> SysBin3, SysBinErrors3;
      

  for (UInt_t i = 0; i < Asym1.size(); i++) {
    SysBin1.push_back(Asym1[i]-Asym2[i]);	SysBinErrors1.push_back(sqrt(pow(AsymE1[i],2)+pow(AsymE2[i],2)));
    SysBin2.push_back(Asym1[i]-Asym3[i]);	SysBinErrors2.push_back(sqrt(pow(AsymE1[i],2)+pow(AsymE3[i],2)));
    SysBin3.push_back(Asym2[i]-Asym3[i]);	SysBinErrors3.push_back(sqrt(pow(AsymE2[i],2)+pow(AsymE3[i],2)));
  }


  TF1 *FitBinnError1[3];
  FitBinnError1[0]	= new TF1("E1","[0]",-1,1);
  FitBinnError1[0]->SetParNames("#Delta #Sigma _{syst}");
  
  cASMBSyE[0] = new TCanvas("","",900,450);  
  cASMBSyE[0]->cd(1);
  BinSysError[0] = new TGraphErrors(SysBin1.size(),AvValueG[5].data(),SysBin1.data(),ErrorG[5].data(),SysBinErrors1.data());
  BinSysError[0]->GetHistogram()->GetYaxis()->SetRangeUser(-0.5,0.5);
  BinSysError[0]->SetTitle("");
  BinSysError[0]->GetXaxis()->SetTitle("cos(#theta^{cm})");
  BinSysError[0]->GetYaxis()->SetTitle("#Delta #Sigma");
  BinSysError[0]->Draw("AP");
  BinSysError[0]->Fit(FitBinnError1[0]);
  string FitErrorsBin1 = FitBinnError1[0]->GetName();
  FitErrorsBin1 = FitErrorsBin1 + "+" + std::to_string(FitBinnError1[0]->GetParError(0));
  FitBinnError1[1] = new TF1("ER1",FitErrorsBin1.c_str(),-1,1);
  FitErrorsBin1 = FitBinnError1[0]->GetName();
  FitErrorsBin1 = FitErrorsBin1 + "-" + std::to_string(FitBinnError1[0]->GetParError(0));
  FitBinnError1[2] = new TF1("ER1",FitErrorsBin1.c_str(),-1,1);

  FitBinnError1[1]->SetLineStyle(2);
  FitBinnError1[2]->SetLineStyle(2);
  FitBinnError1[1]->Draw("same");
  FitBinnError1[2]->Draw("same");
  FitBinnError1[0]->Draw("same");
  cASMBSyE[0]->SaveAs("imagenes/SigmaAsymBinSysErrors1_2.eps");

  TF1 *FitBinnError2[3];
  FitBinnError2[0]	= new TF1("E1","[0]",-1,1);
  FitBinnError2[0]->SetParNames("#Delta #Sigma _{syst}");
  

  cASMBSyE[1] = new TCanvas("","",900,450);  
  cASMBSyE[1]->cd(1);
  BinSysError[1] = new TGraphErrors(SysBin2.size(),AvValueG[5].data(),SysBin2.data(),ErrorG[5].data(),SysBinErrors2.data());
  BinSysError[1]->GetHistogram()->GetYaxis()->SetRangeUser(-0.5,0.5);
  BinSysError[1]->SetTitle("");
  BinSysError[1]->GetXaxis()->SetTitle("cos(#theta^{cm})");
  BinSysError[1]->GetYaxis()->SetTitle("#Delta #Sigma");
  BinSysError[1]->Draw("AP");
  BinSysError[1]->Fit(FitBinnError2[0]);
  string FitErrorsBin2 = FitBinnError2[0]->GetName();
  FitErrorsBin2 = FitErrorsBin2 + "+" + std::to_string(FitBinnError2[0]->GetParError(0));
  FitBinnError2[1] = new TF1("ER1",FitErrorsBin2.c_str(),-1,1);
  FitErrorsBin2 = FitBinnError2[0]->GetName();
  FitErrorsBin2 = FitErrorsBin2 + "-" + std::to_string(FitBinnError2[0]->GetParError(0));
  FitBinnError2[2] = new TF1("ER1",FitErrorsBin2.c_str(),-1,1);

  FitBinnError2[1]->SetLineStyle(2);
  FitBinnError2[2]->SetLineStyle(2);
  FitBinnError2[1]->Draw("same");
  FitBinnError2[2]->Draw("same");
  FitBinnError2[0]->Draw("same");
  
  cASMBSyE[1]->SaveAs("imagenes/SigmaAsymBinSysErrors1_3.eps");

  TF1 *FitBinnError3[3];
  FitBinnError3[0]	= new TF1("E1","[0]",-1,1);
  FitBinnError3[0]->SetParNames("#Delta #Sigma _{syst}");
  
  cASMBSyE[2] = new TCanvas("","",900,450);  
  cASMBSyE[2]->cd(1);
  BinSysError[2] = new TGraphErrors(SysBin3.size(),AvValueG[5].data(),SysBin3.data(),ErrorG[5].data(),SysBinErrors3.data());
  BinSysError[2]->GetHistogram()->GetYaxis()->SetRangeUser(-0.5,0.5);
  BinSysError[2]->SetTitle("");
  BinSysError[2]->GetXaxis()->SetTitle("cos(#theta^{cm})");
  BinSysError[2]->GetYaxis()->SetTitle("#Delta #Sigma");
  BinSysError[2]->Draw("AP");
  BinSysError[2]->Fit(FitBinnError3[0]);
  string FitErrorsBin3 = FitBinnError3[0]->GetName();
  FitErrorsBin3 = FitErrorsBin3 + "+" + std::to_string(FitBinnError3[0]->GetParError(0));
  FitBinnError3[1] = new TF1("ER1",FitErrorsBin3.c_str(),-1,1);
  FitErrorsBin3 = FitBinnError3[0]->GetName();
  FitErrorsBin3 = FitErrorsBin3 + "-" + std::to_string(FitBinnError3[0]->GetParError(0));
  FitBinnError3[2] = new TF1("ER1",FitErrorsBin3.c_str(),-1,1);

  FitBinnError3[1]->SetLineStyle(2);
  FitBinnError3[2]->SetLineStyle(2);
  FitBinnError3[1]->Draw("same");
  FitBinnError3[2]->Draw("same");
  FitBinnError3[0]->Draw("same");
  
  cASMBSyE[2]->SaveAs("imagenes/SigmaAsymBinSysErrors2_3.eps");
  
  // ------ Max Like Binning method comparation ---- //
  
  TCanvas *cASM = new TCanvas("","",900,450);
  TGraphErrors *AsymM = new TGraphErrors(Asym1.size(),&(AvValueG[5][0]),MaxL[5].data(),&(ErrorG[5][0]),MaxLE[5].data());
  TMultiGraph *AsymT = new TMultiGraph();
  AsymM->SetLineColor(kBlack);
  AsymM->SetTitle("Max Like");
  AsymT->Add(AsymG[0]);
  AsymT->Add(AsymM);
  AsymT->Draw("AP");
  AsymT->GetXaxis()->SetTitle("cos(#theta^{cm}_{K^{+}})");
  AsymT->GetYaxis()->SetTitle("#Sigma");
  AsymT->GetHistogram()->GetYaxis()->SetRangeUser(-1,1);
  cASM->BuildLegend();
  cASM->SaveAs("imagenes/SigmaAsymBinMaxLike.eps");

  // ------ Max Like Binning method comparation ---- //

  vector<double> SysMaxL1, SysMaxLErrors1;      

  for (UInt_t i = 0; i < Asym1.size(); i++) {
    SysMaxL1.push_back(Asym1[i]-MaxL[5][i]);	SysMaxLErrors1.push_back(sqrt(pow(AsymE1[i],2)+pow(MaxLE[5][i],2)));
  }

  TF1 *FitMaxLError[3];
  FitMaxLError[0]	= new TF1("E1","[0]",-1,1);
  FitMaxLError[0]->SetParNames("#Delta #Sigma _{syst}");
  
  TCanvas *cMAXSys= new TCanvas("","",900,450);  
  cMAXSys->cd(1);
  TGraphErrors *MaxLSysError = new TGraphErrors(SysMaxL1.size(),AvValueG[5].data(),SysMaxL1.data(),ErrorG[5].data(),SysMaxLErrors1.data());
  MaxLSysError->GetHistogram()->GetYaxis()->SetRangeUser(-0.5,0.5);
  MaxLSysError->SetTitle("");
  MaxLSysError->GetXaxis()->SetTitle("cos(#theta^{cm})");
  MaxLSysError->GetYaxis()->SetTitle("#Delta #Sigma");
  MaxLSysError->Draw("AP");
  MaxLSysError->Fit(FitMaxLError[0]);
  string FitErrorsMaxL = FitMaxLError[0]->GetName();
  FitErrorsMaxL = FitErrorsMaxL + "+" + std::to_string(FitMaxLError[0]->GetParError(0));
  FitMaxLError[1] = new TF1("ER1",FitErrorsMaxL.c_str(),-1,1);
  FitErrorsMaxL = FitMaxLError[0]->GetName();
  FitErrorsMaxL = FitErrorsMaxL + "-" + std::to_string(FitMaxLError[0]->GetParError(0));
  FitMaxLError[2] = new TF1("ER1",FitErrorsMaxL.c_str(),-1,1);

  FitMaxLError[1]->SetLineStyle(2);
  FitMaxLError[2]->SetLineStyle(2);
  FitMaxLError[1]->Draw("same");
  FitMaxLError[2]->Draw("same");
  FitMaxLError[0]->Draw("same");
  cMAXSys->SaveAs("imagenes/SigmaAsymBinMaxLikeSys.eps");

  
  TCanvas *cASME[6];
  float BimEn = 1.3;

  for (UInt_t i = 0; i < 6; i++) {
    cASME [i] = new TCanvas("","",900,450);
    TGraphErrors *AsymME = new TGraphErrors(Asym1.size(),&(AvValueG[i][0]),MaxL[i].data(),&(ErrorG[i][0]),MaxLE[i].data());
    string Name;
    Name = "Max Like " + std::to_string(BimEn);
    AsymME->SetTitle(Name.c_str());
    AsymME->GetXaxis()->SetTitle("cos(#theta^{cm}_{K^{+}})");
    AsymME->GetYaxis()->SetTitle("#Sigma");
    AsymME->GetHistogram()->GetYaxis()->SetRangeUser(-1,1);
    AsymME->Draw("AP");
    string AsymS;
    AsymS = "imagenes/SigmaAsymMaxLike" + std::to_string(i) + ".eps";
    cASME[i]->SaveAs(AsymS.c_str());
    BimEn+=0.2;
  }

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

// Function to asymmetry

Double_t fitf(Double_t *x, Double_t *par){
  double num=0, denom=0;
  num=par[0]-1.0-(par[0]*par[1]+1.0)/(par[1]+1.0)*2*par[2]*par[3]*TMath::Cos(2*x[0]*TMath::DegToRad());
  denom=par[0]+1.0-(par[0]*par[1]-1.0)/(par[1]+1.0)*2*par[2]*par[3]*TMath::Cos(2*x[0]*TMath::DegToRad());
  Double_t fitval =num/denom;
  return fitval;
}

#endif
