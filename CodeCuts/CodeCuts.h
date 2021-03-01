/*******************************************/
// Author: Ivan Dario Piernagorda Peña     
// Author: Karen Milena Rodriguez      
// Date:   15/02/2020                  
// Title:  Cuts to Data                              
/*******************************************/

#include "./Histograms.h"
#include "./include/Libraries.h"
#include "./include/Miscelaneous.h"
#include "./src/DataEvent.cpp"

#ifndef CODECUTS_H
#define CODECUTS_H

class Codecuts : public Histograms {

private:

public:
  Codecuts(){DoHistograms();}
  void CodeCuts();
  void CodeCutsAsym();

};

void Codecuts::CodeCuts(){

  
  string treeName="g13b";  
  string fileNamePERP="./SKIMS/ListFiles.txt";
  DataEvent *myDataList=new DataEvent(fileNamePERP,treeName, 35);
  
  vector<string> PolTableName;
  ListFilesAtDir("./TABLES/ListTables", PolTableName);
  
  map<vector<float>,int> keysPlane;

  //Average Table Polarization
  vector<double> vecInitDob(2);   vector<vector<double>> AvP(10,vecInitDob);
  vector<int> vecInitInt(2);      vector<vector<int>> ItP(10,vecInitInt);

  //Table Events X 14 cuts
  vector<int> Events(19);

  const int NumbOfPolFiles=PolTableName.size();
  for (int i=0;i<NumbOfPolFiles;i++){
    char *cstr = const_cast<char *>(PolTableName[i].c_str());
    LoadPolTable(i,cstr,keysPlane);
  }

  //----- Save Vars To MaxLike and Binning Method ---//
  
  TFile *Cuts = new TFile("Cuts.root","RECREATE");
  TTree *FinalCut = new TTree("FinalCut","Final cuts for MaxLike and Binning method");
  float CohE, CohEN, PhotoP;
  int NumEv, CohP;
  TLorentzVector Proton, Kaon, Sig, WB;


  FinalCut->Branch("CohE",&CohE);
  FinalCut->Branch("CohEN",&CohEN);
  FinalCut->Branch("CohP",&CohP);
  FinalCut->Branch("NumEv",&NumEv);
  FinalCut->Branch("PhotoPol",&PhotoP);
  FinalCut->Branch("VProton",&Proton);
  FinalCut->Branch("VKaon",&Kaon);
  FinalCut->Branch("VSig",&Sig);
  FinalCut->Branch("VWBoost",&WB);
  
  
  
  while (myDataList->getEntry()<myDataList->getEntries()){

    myDataList->getNextEntry();
    if (myDataList->getEntry() % 1000 == 0){
      if(myDataList->getCoh_plan() == 0){
	fprintf (stderr, "Looped %s : %.2f percent, CohEdge %f BeamE %f \r", "PARA" ,myDataList->getEntry()*100.0/myDataList->getEntries(),myDataList->getCoh_edge_nom(),myDataList->getBeam_en());
	fflush (stderr);
      }
      else {
	fprintf (stderr, "Looped %s : %.2f percent, CohEdge %f BeamE %f \r", "PERP" ,myDataList->getEntry()*100.0/myDataList->getEntries(),myDataList->getCoh_edge_nom(),myDataList->getBeam_en());
	fflush (stderr);
      }
    }
      
      
    Events[0]++;         //Events Without Cuts

    //------------------ Delta Beta ---------------//
    double deltbeta[3]     = {};
    double deltbetacut[3]  = {};
      
    //---------- Vertex ------------ //
      
      
    h_Vertex->Fill(myDataList->getEVNT_vertex(1).Z());
    if(myDataList->getEVNT_vertex(1).Z()<-39.0  || myDataList->getEVNT_vertex(1).Z()>-1.0) continue;

    Events[1]++;         //Events With Vertex Cut
    
    for (int i=0;i<myDataList->getNum_chargedtracks();i++){
      deltbeta[i]=myDataList->getEVNT_track(i).Beta()-myDataList->getEVNT_bem(i);
      h_DeltaBe[i]->Fill(myDataList->getEVNT_track(i).Rho(),deltbeta[i]);
      h_BeVSp[i]->Fill(myDataList->getEVNT_track(i).Rho(),myDataList->getEVNT_bem(i));
      h_BeVSpT->Fill(myDataList->getEVNT_track(i).Rho(),myDataList->getEVNT_bem(i));

      if(i == 0 || i == 1)
	h_Mass[0]->Fill(myDataList->getEVNT_track(i).Mt());
      else
	h_Mass[1]->Fill(myDataList->getEVNT_track(i).Mt());
    }

    
    
    //------------------ Delta Beta with Cuts ---------------//
    
    for (int i=0;i<myDataList->getNum_chargedtracks();i++){
      deltbetacut[i]=myDataList->getEVNT_track(i).Beta()-myDataList->getEVNT_bem(i);
	
      if(deltbetacut[2] > 0.05  || deltbetacut[2] < -0.05) continue; 
      if(deltbetacut[0] > 0.02  || deltbetacut[0] < -0.02) continue;
      if(deltbetacut[1] > 0.025 || deltbetacut[1] < -0.025) continue;


      h_DeltaBecut[i]->Fill(myDataList->getEVNT_track(i).Rho(),deltbeta[i]);
      h_BeVSpcut[i]->Fill(myDataList->getEVNT_track(i).Rho(),myDataList->getEVNT_bem(i));

    }

    if(deltbetacut[2] > 0.05  || deltbetacut[2] < -0.05) continue;
    Events[2]++;         //Events With DeltaB cut
    if(deltbetacut[0] > 0.02  || deltbetacut[0] < -0.02) continue;
    Events[3]++;         //Events With DeltaB cut
    if(deltbetacut[1] > 0.025 || deltbetacut[1] < -0.025) continue;
    Events[4]++;         //Events With DeltaB cut
       
    //------------------Correlation Theta-Phi, -----------------------/
    
    h_ThePhi[0]->Fill(myDataList->getEVNT_track(0).Phi()*TMath::RadToDeg(), myDataList->getEVNT_track(0).Theta()*TMath::RadToDeg());  
    h_ThePhi[1]->Fill(myDataList->getEVNT_track(1).Phi()*TMath::RadToDeg(), myDataList->getEVNT_track(1).Theta()*TMath::RadToDeg());
    h_ThePhi[2]->Fill(myDataList->getEVNT_track(2).Phi()*TMath::RadToDeg(), myDataList->getEVNT_track(2).Theta()*TMath::RadToDeg());
    
    
    //-----Cuts due to detector geometry (Fiduciary cuts)---------//
    //------------------------------------------------------------//
      
    Double_t phiproton_cut;
    phiproton_cut= myDataList->getEVNT_track(0).Phi()*TMath::RadToDeg();
    
    if( (phiproton_cut <= -25)   && (phiproton_cut >= -35)  ) continue;
    if( (phiproton_cut <= -85)   && (phiproton_cut >= -95)  ) continue;
    if( (phiproton_cut <= -145)  && (phiproton_cut >= -155) ) continue;
    if( (phiproton_cut >= 25)    && (phiproton_cut <= 35)   ) continue;
    if( (phiproton_cut >= 85)    && (phiproton_cut <= 95)   ) continue;
    if( (phiproton_cut >= 145)   && (phiproton_cut <= 155)  ) continue;

    if(F_ThePhiProt[0]->Eval(phiproton_cut) > myDataList->getEVNT_track(0).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[1]->Eval(phiproton_cut) > myDataList->getEVNT_track(0).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[2]->Eval(phiproton_cut) > myDataList->getEVNT_track(0).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[3]->Eval(phiproton_cut) > myDataList->getEVNT_track(0).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[4]->Eval(phiproton_cut) > myDataList->getEVNT_track(0).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[5]->Eval(phiproton_cut) > myDataList->getEVNT_track(0).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[6]->Eval(phiproton_cut) > myDataList->getEVNT_track(0).Theta()*TMath::RadToDeg()) continue;
       
    h_ThePhicut[0]->Fill(phiproton_cut, myDataList->getEVNT_track(0).Theta()*TMath::RadToDeg());

    Events[5]++;         //Events With Phi-Theta Cuts
    //kaon cuts
    Double_t phikaon_cut;
    phikaon_cut= myDataList->getEVNT_track(1).Phi()*TMath::RadToDeg();
        
    if( (phikaon_cut <= -25)   && (phikaon_cut >= -35)  ) continue;
    if( (phikaon_cut <= -85)   && (phikaon_cut >= -95)  ) continue;
    if( (phikaon_cut <= -145)  && (phikaon_cut >= -155) ) continue;
    if( (phikaon_cut >= 25)    && (phikaon_cut <= 35)   ) continue;
    if( (phikaon_cut >= 85)    && (phikaon_cut <= 95)   ) continue;
    if( (phikaon_cut >= 145)   && (phikaon_cut <= 155)  ) continue;

    if(F_ThePhiProt[0]->Eval(phikaon_cut) > myDataList->getEVNT_track(1).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[1]->Eval(phikaon_cut) > myDataList->getEVNT_track(1).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[2]->Eval(phikaon_cut) > myDataList->getEVNT_track(1).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[3]->Eval(phikaon_cut) > myDataList->getEVNT_track(1).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[4]->Eval(phikaon_cut) > myDataList->getEVNT_track(1).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[5]->Eval(phikaon_cut) > myDataList->getEVNT_track(1).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[6]->Eval(phikaon_cut) > myDataList->getEVNT_track(1).Theta()*TMath::RadToDeg()) continue;

    h_ThePhicut[1]->Fill(phikaon_cut, myDataList->getEVNT_track(1).Theta()*TMath::RadToDeg());

    Events[6]++;         //Events With Phi-Theta Cuts   
      
    //Pion cuts
    
    Double_t phiPion_cut;
    phiPion_cut= myDataList->getEVNT_track(2).Phi()*TMath::RadToDeg();
      
    if( (phiPion_cut <= -25)   && (phiPion_cut >= -35)  ) continue;
    if( (phiPion_cut <= -85)   && (phiPion_cut >= -95)  ) continue;
    if( (phiPion_cut <= -145)  && (phiPion_cut >= -155) ) continue;
    if( (phiPion_cut >= 25)    && (phiPion_cut <= 35)   ) continue;
    if( (phiPion_cut >= 85)    && (phiPion_cut <= 95)   ) continue;
    if( (phiPion_cut >= 145)   && (phiPion_cut <= 155)  ) continue;

    if(F_ThePhiProt[0]->Eval(phiPion_cut) > myDataList->getEVNT_track(2).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[1]->Eval(phiPion_cut) > myDataList->getEVNT_track(2).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[2]->Eval(phiPion_cut) > myDataList->getEVNT_track(2).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[3]->Eval(phiPion_cut) > myDataList->getEVNT_track(2).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[4]->Eval(phiPion_cut) > myDataList->getEVNT_track(2).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[5]->Eval(phiPion_cut) > myDataList->getEVNT_track(2).Theta()*TMath::RadToDeg() ||
       F_ThePhiProt[6]->Eval(phiPion_cut) > myDataList->getEVNT_track(2).Theta()*TMath::RadToDeg()) continue;

    
    h_ThePhicut[2]->Fill(phiPion_cut, myDataList->getEVNT_track(2).Theta()*TMath::RadToDeg());

    Events[7]++;         //Events With Phi-Theta Cuts   
    
    //------------------ Photons, Delta T  ------------------ // 
      
    for (int i=0;i<myDataList->getNum_photons();i++){
      //LOOP OVER ALL PHOTONS
      if (fabs(deltbeta[1])<0.025){
	h_DeltaTall[0]->Fill(myDataList->getDelt_t_k(i));
	h_DeltaTallvsp[0]->Fill(myDataList->getEVNT_track(1).Rho(),myDataList->getDelt_t_k(i));
      }
      if (fabs(deltbeta[2])<0.05){
	h_DeltaTall[1]->Fill(myDataList->getDelt_t_pi(i));
	h_DeltaTallvsp[1]->Fill(myDataList->getEVNT_track(2).Rho(),myDataList->getDelt_t_pi(i));
      }
    }
      
    //------------Delta T with Cuts ----------- //
      
    if (myDataList->getNumph_k()==1)
      h_DeltaT[0]->Fill(myDataList->getDelt_t_k(myDataList->getIndex_k(0)));
    if(myDataList->getNumph_pi()==1)
      h_DeltaT[1]->Fill(myDataList->getDelt_t_pi(myDataList->getIndex_pi(0)));
      
    if (myDataList->getNumph_k()!=1) continue;
    Events[8]++;         //Events With Delta T Photon cut
    //if (myDataList->getNumph_pi()!=1) continue;
      
    //--------------- Energy loss ----------- //
    double CorrElossKa, CorrElossPr, CorrElossPi;

    CorrElossKa=(myDataList->geteloss_track(0).Rho() - myDataList->getEVNT_track(0).Rho() )/myDataList->getEVNT_track(0).Rho();
    CorrElossPr=(myDataList->geteloss_track(1).Rho() - myDataList->getEVNT_track(1).Rho() )/myDataList->getEVNT_track(1).Rho();
    CorrElossPi=(myDataList->geteloss_track(2).Rho() - myDataList->getEVNT_track(2).Rho() )/myDataList->getEVNT_track(2).Rho();
      
    h_eloss[0]->Fill(CorrElossKa);
    h_eloss[1]->Fill(CorrElossPr);
    h_eloss[2]->Fill(CorrElossPi);

    h_Celoss[0]->Fill(myDataList->getEVNT_track(0).Rho(), CorrElossKa);
    h_Celoss[1]->Fill(myDataList->getEVNT_track(1).Rho(), CorrElossPr);
    h_Celoss[2]->Fill(myDataList->getEVNT_track(2).Rho(), CorrElossPi);
    
    


    

    //--------------- Coh Edge -------------- //
    
      
      
    h_TagrEpho[0]->Fill(myDataList->getTAGR_epho(myDataList->getIndex_k(0))*1000.0);
    if (myDataList->getTAGR_epho(myDataList->getIndex_k(0))*1000.0 > myDataList->getCoh_edge()) continue;
    h_TagrEpho[1]->Fill(myDataList->getTAGR_epho(myDataList->getIndex_k(0))*1000.0);        Events[9]++;         //Events With Tager Epho
    if (myDataList->getTAGR_epho(myDataList->getIndex_k(0))*1000< myDataList->getCoh_edge()-200.0) continue;
    h_TagrEpho[2]->Fill(myDataList->getTAGR_epho(myDataList->getIndex_k(0))*1000.0);        Events[10]++;         //Events With Tager Epho
    if (fabs(myDataList->getCoh_edge()-myDataList->getCoh_edge_nom()*1000)>15)continue;     Events[11]++;         //Events With Tager Epho
    if (myDataList->getTrip_flag()!=0)continue;                                             Events[12]++;         //Events With Tager Epho
    if (myDataList->getCoh_plan()!=0 && myDataList->getCoh_plan()!=1)continue;              Events[13]++;         //Events With Tager Epho

    vector<float> Keys(3);
    if(myDataList->getCoh_edge_nom() == float(1.3)){
      Keys[0]=float(4.2);
      Keys[1]=myDataList->getCoh_edge_nom();
      Keys[2]=float(myDataList->getCoh_plan());
  }
    else{
      Keys[0]=float(myDataList->getBeam_en());
      Keys[1]=myDataList->getCoh_edge_nom();
      Keys[2]=float(myDataList->getCoh_plan());
  }
        
    double PhotoPol=0;
    PhotoPol=GetPol(keysPlane[Keys], myDataList->getCoh_edge(), myDataList->getTAGR_epho(myDataList->getIndex_k(0))*1000.0, 8, 0.2,0.3);
    
    if (PhotoPol<0.5) continue;
    
    GetPolAv(Keys,ItP,AvP,PhotoPol);
    Events[14]++;         //Events With PhotoPol Tables

    //-------------- Reconstruction --------- //
      
    TLorentzVector photon, deuteron, kaon, proton ,pion, Neutron, WBoost;			//Principal Reaction
    TLorentzVector kaonpion, kaonproton, pionkaon;						//MIS-identification particles
    TLorentzVector MMNeut_kaon, MMNeut_KPi, MMNeut_KP, MMNeut_PiK ,MMSigma;			//Missing mass
    TLorentzVector Sigma, Lambda;								//Invariant Mass
    
    photon.SetXYZM(0,0,myDataList->getTAGR_epho(myDataList->getIndex_k(0)),0);
    deuteron.SetXYZM(0,0,0,1.8756);
    
    double Px_kaonpion = myDataList->getEVNT_track(1).Rho()* sin(myDataList->getEVNT_track(1).Theta())* cos(myDataList->getEVNT_track(1).Phi());
    double Py_kaonpion = myDataList->getEVNT_track(1).Rho()* sin(myDataList->getEVNT_track(1).Theta())* sin(myDataList->getEVNT_track(1).Phi());
    double Pz_kaonpion = myDataList->getEVNT_track(1).Rho()*(myDataList->getEVNT_track(1).CosTheta());

    double Px_pionkaon = myDataList->getEVNT_track(2).Rho()* sin(myDataList->getEVNT_track(2).Theta())* cos(myDataList->getEVNT_track(2).Phi());
    double Py_pionkaon = myDataList->getEVNT_track(2).Rho()* sin(myDataList->getEVNT_track(2).Theta())* sin(myDataList->getEVNT_track(2).Phi());
    double Pz_pionkaon = myDataList->getEVNT_track(2).Rho()*(myDataList->getEVNT_track(2).CosTheta());

    // This mass is of Pion-, because we need remove the background of Pion-
    kaonpion.SetXYZM(Px_kaonpion, Py_kaonpion, Pz_kaonpion, 0.139);       			
    kaonproton.SetXYZM(Px_kaonpion, Py_kaonpion, Pz_kaonpion,0.9383);
    pionkaon.SetXYZM(Px_pionkaon, Py_pionkaon, Pz_pionkaon,0.4937);
    
    proton 	= myDataList->geteloss_track(0);
    kaon 	= myDataList->geteloss_track(1);
    pion 	= myDataList->geteloss_track(2);
    MMNeut_kaon = photon + deuteron - proton - kaon - pion;
    Neutron.SetXYZM(MMNeut_kaon.Px(), MMNeut_kaon.Py(), MMNeut_kaon.Pz(), 0.939);
    Sigma 	= pion + Neutron;
    Lambda 	= pion + proton;
    MMSigma  	= photon + deuteron - proton - kaon;           			// Correlation with invariant mass (lambda)
    WBoost   	= photon + deuteron; 						// to make Boost
      
    MMNeut_KPi 	= photon + deuteron - proton - kaonpion - pion;     	   	// This missing mass is with the Pion-
    MMNeut_KP	= photon + deuteron - proton - kaonproton - pion;     	   	// This missing mass is with the Proton
    MMNeut_PiK 	= photon + deuteron - proton - pionkaon - kaon;     	   	// This missing mass is with the Kaon-      
      
    h_MissingMass->Fill(MMNeut_kaon.M());

    h_MissingMass_kaonpion->Fill(MMNeut_KPi.M());
    h_MissingMass_kaonproton->Fill(MMNeut_KP.M());
    h_MissingMass_pionkaon->Fill(MMNeut_PiK.M());
    //h_MissingPvsIMMass->Fill(MMNeut_kaon.M(),MMNeut_kaon.P());
    h_MissingMass_vsMissingMasskaonpion[0]->Fill(MMNeut_kaon.M(), MMNeut_KPi.M());
    h_MissingMass_vsMissingMasskaonproton[0]->Fill(MMNeut_kaon.M(),MMNeut_KP.M());
    h_MissingMass_vsMissingMasspionkaon[0]->Fill(MMNeut_kaon.M(),MMNeut_PiK.M());
    
            
    
    if(MMNeut_KPi.M() < 0.98) continue;       //Cut from correlation MM
    h_MissingMass_vsMissingMasskaonpion[1]->Fill(MMNeut_kaon.M(), MMNeut_KPi.M());
    if(MMNeut_KP.M()  > 0.75) continue; 
    h_MissingMass_vsMissingMasskaonproton[1]->Fill(MMNeut_kaon.M(),MMNeut_KP.M());
    if(MMNeut_PiK.M() > 0.70) continue; 
    h_MissingMass_vsMissingMasspionkaon[1]->Fill(MMNeut_kaon.M(),MMNeut_PiK.M());
    Events[15]++;         //Events With NOT PION, YES Kaon
    // Double_t El = TMath::Power((proton.P()-offsetx)*cos(angle)+(MMSigma.M()-offsety)*sin(angle),2)/TMath::Power(radx,2)
    //   +TMath::Power((proton.P()-offsetx)*sin(angle)-(MMSigma.M()-offsety)*cos(angle),2)/TMath::Power(rady,2);
    // if(El <= 1) continue;
    



    
    //----------Correlación momentums vs missing mass------------------------//
      
    if( Lambda.M()<1.11 || Lambda.M()>1.132) 
      h_MissingPvsIMMass[0]->Fill(Sigma.M(),MMNeut_kaon.P());
    
    if( Sigma.M()<1.08 || Sigma.M()>1.3)
      h_MissingPvsIMMass[1]->Fill(Lambda.M(),MMNeut_kaon.P());

     h_MissingMassvsSigmaMass->Fill(Sigma.M(), MMNeut_kaon.P());
    
    h_MissingMasscut->Fill(MMNeut_kaon.M());
    h_MissingMass_kaonpioncut->Fill(MMNeut_KPi.M());
      
    // h_MissingP->Fill(MMNeut_kaon.P());
      
    h_MissingPvsSigmaMass->Fill(Sigma.M(), MMNeut_kaon.P());
          
    h_InvariantMass->Fill(Sigma.M());
    h_LambdaMass->Fill(Lambda.M());
    h_InvMassLambda_vsInvMassSigma->Fill(Sigma.M(), Lambda.M());
    h_MMassSigma->Fill(MMSigma.M());
    h_MMNeutron_vsMMassSigma[0]->Fill(MMNeut_kaon.M(),MMSigma.M());

    
    //------------------------Comparación de sigmas-------------//
    if( Lambda.M()<1.108 || Lambda.M()>1.124)
      h_InvariantMasscut[0]->Fill(Sigma.M());
      
    if( Lambda.M()<1.11 || Lambda.M()>1.132)
      h_InvariantMasscut[1]->Fill(Sigma.M());
      
    if( Lambda.M()<1.092 || Lambda.M()>1.14)
      h_InvariantMasscut[2]->Fill(Sigma.M());
      
      
    //------------- Comparación de missing momentums-----------//
      
    if( Lambda.M()<1.096 || Lambda.M()>1.136) 
      h_MissingP[0]->Fill(MMNeut_kaon.P());   
    if( Sigma.M()<1.08 || Sigma.M()>1.3)
      h_MissingP[1]->Fill(MMNeut_kaon.P());


    if ( Lambda.M()>=1.1 && Lambda.M()<=1.132) continue;                   //Cut for LamdaMass in +/- 8sigma
    Events[16]++;         //Events With Lambda cuts


    //--------Correlation Momentums--------------//
    h_CorrelationMMomentum->Fill(Sigma.P(), Lambda.P());
    
    // Events[17]++;         //Events With Mass neutron range
    if(MMNeut_kaon.P()<=0.2) continue;                                    //Cut for rescattering
    Events[18]++;         //Events With Neutron Rescattering Lambada cuts


    h_MMNeutron_vsMMassSigma[1]->Fill(MMNeut_kaon.M(),MMSigma.M());       //Cut MM neutron and MM sigma 
    h_InvariantMasscut[3]->Fill(Sigma.M());
   
    h_MMassSigmaCut->Fill(MMSigma.M());
      
    // Pion -  (Por si las moscas)
    h_DeltaBVSInvariantMass->Fill(Sigma.M(),deltbeta[2]);
    h_DeltaBVSMissingMass->Fill(MMNeut_kaon.M(),deltbeta[2]);
    h_DeltaBVSMissingMomentum->Fill(MMNeut_kaon.P(),deltbeta[2]);

    //------------Momentum proton----------------.//
    h_MomentumProton->Fill(proton.P());
   
    //-----------------BOOST------------------------------//

    CohE 	= myDataList->getCoh_edge();
    CohEN	= myDataList->getCoh_edge_nom();
    CohP	= myDataList->getCoh_plan();
    PhotoP	= PhotoPol;
    Proton 	= proton;
    Kaon	= kaon;
    Sig 	= Sigma;
    WB		= WBoost;
    NumEv	= myDataList->getHEAD_eventnum();
    
    FinalCut->Fill();
         
  }
  cout<<endl;

  FinalCut->Write();
  Cuts->Write();
  Cuts->Close();
  GetPolAvTable(ItP,AvP);
  GetPolAvTableLatex(ItP, AvP, "./PolTable.tex","Polarization Tables","poltab");
  GetEventPercentLatex(Events, "./EventCuts.tex", "Event Cut Tables", "eventtab");
  GetEventPercent(Events);
  //DoCanvas();
  
}


void Codecuts::CodeCutsAsym(){
 
  TChain *Cuts = new TChain("FinalCut");
  TLorentzVector *Proton = NULL , *Kaon = NULL, *Sigma = NULL, *WBoost = NULL;
  float CohE = 0, CohEN = 0, PhotoPol = 0;
  int NumEv = 0, Event = 0, CohP = 0;
  
  Cuts->Add("Cuts.root");
  Cuts->SetBranchAddress("VProton",&Proton);
  Cuts->SetBranchAddress("VKaon",&Kaon);
  Cuts->SetBranchAddress("VSig",&Sigma);
  Cuts->SetBranchAddress("VWBoost",&WBoost);
  Cuts->SetBranchAddress("CohE",&CohE);
  Cuts->SetBranchAddress("CohEN",&CohEN);
  Cuts->SetBranchAddress("CohP",&CohP);
  Cuts->SetBranchAddress("NumEv",&NumEv);
  Cuts->SetBranchAddress("PhotoPol",&PhotoPol);
  
  while(Event < Cuts->GetEntries()){

    if(Event < Cuts->GetEntries()){
      Cuts->GetEvent(Event);
      Event++;
    }
    
    if (Event % 1000 == 0){
      if(CohP == 0){
	fprintf (stderr, "Looped %s : %.2f percent, CohEdge %f \r", "PARA" ,Event*100.0/Cuts->GetEntries(),CohEN);
	fflush (stderr);
      }
      else {
	fprintf (stderr, "Looped %s : %.2f percent, CohEdge %f \r", "PERP" ,Event*100.0/Cuts->GetEntries(),CohEN);
	fflush (stderr);
      }
    }

  
    TVector3 b=WBoost->BoostVector();
    Proton->Boost(-b);
    Kaon->Boost(-b);
    Sigma->Boost(-b);
    double ProtonCosThetaCM=TMath::Cos(Proton->Theta());
    double KaonCosThetaCM=TMath::Cos(Kaon->Theta());
    double SigmaCosThetaCM=TMath::Cos(Sigma->Theta());
    double PhiCM=Kaon->Phi()*TMath::RadToDeg();

    
    if(CohEN == float(1.7)){

      
      h_CosThetaCM17[0]->Fill(ProtonCosThetaCM);
      h_CosThetaCM17[1]->Fill(KaonCosThetaCM);
      h_CosThetaCM17[2]->Fill(SigmaCosThetaCM);

      if (CohP == 0){
	if(KaonCosThetaCM < -0.70) 
	  { h_kaonPhiPA1[0]->Fill(PhiCM); h_kaonPhiPA2[0]->Fill(PhiCM); h_kaonPhiPA3[0]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= -0.70 && KaonCosThetaCM <= -0.36)
	  { h_kaonPhiPA1[1]->Fill(PhiCM); h_kaonPhiPA2[1]->Fill(PhiCM); h_kaonPhiPA3[1]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= -0.36 && KaonCosThetaCM <= -0.06)
	  { h_kaonPhiPA1[2]->Fill(PhiCM); h_kaonPhiPA2[2]->Fill(PhiCM); h_kaonPhiPA3[2]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= -0.06 && KaonCosThetaCM <= 0.12)
	  { h_kaonPhiPA1[3]->Fill(PhiCM); h_kaonPhiPA2[3]->Fill(PhiCM); h_kaonPhiPA3[3]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= 0.12 && KaonCosThetaCM <= 0.26)
	  { h_kaonPhiPA1[4]->Fill(PhiCM); h_kaonPhiPA2[4]->Fill(PhiCM); h_kaonPhiPA3[4]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= 0.26 && KaonCosThetaCM <= 0.36)
	  { h_kaonPhiPA1[5]->Fill(PhiCM); h_kaonPhiPA2[5]->Fill(PhiCM); h_kaonPhiPA3[5]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= 0.36 && KaonCosThetaCM <= 0.44)
	  { h_kaonPhiPA1[6]->Fill(PhiCM); h_kaonPhiPA2[6]->Fill(PhiCM); h_kaonPhiPA3[6]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= 0.44 && KaonCosThetaCM <= 0.52)
	  { h_kaonPhiPA1[7]->Fill(PhiCM); h_kaonPhiPA2[7]->Fill(PhiCM); h_kaonPhiPA3[7]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= 0.52 && KaonCosThetaCM <= 0.6)
	  { h_kaonPhiPA1[8]->Fill(PhiCM); h_kaonPhiPA2[8]->Fill(PhiCM); h_kaonPhiPA3[8]->Fill(PhiCM); }
	else if (KaonCosThetaCM > 0.6)
	  { h_kaonPhiPA1[9]->Fill(PhiCM); h_kaonPhiPA2[9]->Fill(PhiCM); h_kaonPhiPA3[9]->Fill(PhiCM); }
      }
      
      else if (CohP == 1){	
	if(KaonCosThetaCM < -0.70) 
	  { h_kaonPhiPE1[0]->Fill(PhiCM); h_kaonPhiPE2[0]->Fill(PhiCM); h_kaonPhiPE3[0]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= -0.70 && KaonCosThetaCM <= -0.36)
	  { h_kaonPhiPE1[1]->Fill(PhiCM); h_kaonPhiPE2[1]->Fill(PhiCM); h_kaonPhiPE3[1]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= -0.36 && KaonCosThetaCM <= -0.06)
	  { h_kaonPhiPE1[2]->Fill(PhiCM); h_kaonPhiPE2[2]->Fill(PhiCM); h_kaonPhiPE3[2]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= -0.06 && KaonCosThetaCM <= 0.12)
	  { h_kaonPhiPE1[3]->Fill(PhiCM); h_kaonPhiPE2[3]->Fill(PhiCM); h_kaonPhiPE3[3]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= 0.12 && KaonCosThetaCM <= 0.26)
	  { h_kaonPhiPE1[4]->Fill(PhiCM); h_kaonPhiPE2[4]->Fill(PhiCM); h_kaonPhiPE3[4]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= 0.26 && KaonCosThetaCM <= 0.36)
	  { h_kaonPhiPE1[5]->Fill(PhiCM); h_kaonPhiPE2[5]->Fill(PhiCM); h_kaonPhiPE3[5]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= 0.36 && KaonCosThetaCM <= 0.44)
	  { h_kaonPhiPE1[6]->Fill(PhiCM); h_kaonPhiPE2[6]->Fill(PhiCM); h_kaonPhiPE3[6]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= 0.44 && KaonCosThetaCM <= 0.52)
	  { h_kaonPhiPE1[7]->Fill(PhiCM); h_kaonPhiPE2[7]->Fill(PhiCM); h_kaonPhiPE3[7]->Fill(PhiCM); }
	else if (KaonCosThetaCM >= 0.52 && KaonCosThetaCM <= 0.6)
	  { h_kaonPhiPE1[8]->Fill(PhiCM); h_kaonPhiPE2[8]->Fill(PhiCM); h_kaonPhiPE3[8]->Fill(PhiCM); }
	else if (KaonCosThetaCM > 0.6)
	  { h_kaonPhiPE1[9]->Fill(PhiCM); h_kaonPhiPE2[9]->Fill(PhiCM); h_kaonPhiPE3[9]->Fill(PhiCM); }
      }
    }

    
    h_CosThetaCM[0]->Fill(ProtonCosThetaCM);
    h_CosThetaCM[1]->Fill(KaonCosThetaCM);
    h_CosThetaCM[2]->Fill(SigmaCosThetaCM);

    h_Theta[0]->Fill(Proton->Theta());
    h_Theta[1]->Fill(Kaon->Theta());
    h_Theta[2]->Fill(Sigma->Theta());
    
    h_ThetaCorr[0]->Fill(ProtonCosThetaCM,KaonCosThetaCM);
    h_ThetaCorr[1]->Fill(ProtonCosThetaCM,SigmaCosThetaCM);
    h_ThetaCorr[2]->Fill(SigmaCosThetaCM,KaonCosThetaCM);


      
    //--------------- Bins Cos Theta Kaon----------------//

    //0 is for PARA
    //1 is for PERP
  
    if (CohP == 0){
      if(KaonCosThetaCM < -0.72){
	MEASGamma[CohEN].at(0).push_back(PhotoPol);
	MEASPhip[CohEN].at(0).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= -0.72 && KaonCosThetaCM <= -0.48){
	MEASGamma[CohEN].at(1).push_back(PhotoPol);
	MEASPhip[CohEN].at(1).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= -0.48 && KaonCosThetaCM <= -0.16){
	MEASGamma[CohEN].at(2).push_back(PhotoPol);
	MEASPhip[CohEN].at(2).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= -0.16 && KaonCosThetaCM <= 0.04){
	MEASGamma[CohEN].at(3).push_back(PhotoPol);
	MEASPhip[CohEN].at(3).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= 0.04 && KaonCosThetaCM <= 0.18){
	MEASGamma[CohEN].at(4).push_back(PhotoPol);
	MEASPhip[CohEN].at(4).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= 0.18 && KaonCosThetaCM <= 0.28){
	MEASGamma[CohEN].at(5).push_back(PhotoPol);
	MEASPhip[CohEN].at(5).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= 0.28 && KaonCosThetaCM <= 0.38){
	MEASGamma[CohEN].at(6).push_back(PhotoPol);
	MEASPhip[CohEN].at(6).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= 0.38 && KaonCosThetaCM <= 0.46){
	MEASGamma[CohEN].at(7).push_back(PhotoPol);
	MEASPhip[CohEN].at(7).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= 0.46 && KaonCosThetaCM <= 0.54){
	MEASGamma[CohEN].at(8).push_back(PhotoPol);
	MEASPhip[CohEN].at(8).push_back(PhiCM);
      }
      else if (KaonCosThetaCM > 0.54){
	MEASGamma[CohEN].at(9).push_back(PhotoPol);
	MEASPhip[CohEN].at(9).push_back(PhiCM);
      }
    }
      
    else if (CohP == 1){
      if(KaonCosThetaCM < -0.72){
	MEASGamma[CohEN].at(0).push_back(-PhotoPol);
	MEASPhip[CohEN].at(0).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= -0.72 && KaonCosThetaCM <= -0.48){
	MEASGamma[CohEN].at(1).push_back(-PhotoPol);
	MEASPhip[CohEN].at(1).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= -0.48 && KaonCosThetaCM <= -0.16){
	MEASGamma[CohEN].at(2).push_back(-PhotoPol);
	MEASPhip[CohEN].at(2).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= -0.16 && KaonCosThetaCM <= 0.04){
	MEASGamma[CohEN].at(3).push_back(-PhotoPol);
	MEASPhip[CohEN].at(3).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= 0.04 && KaonCosThetaCM <= 0.18){
	MEASGamma[CohEN].at(4).push_back(-PhotoPol);
	MEASPhip[CohEN].at(4).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= 0.18 && KaonCosThetaCM <= 0.28){
	MEASGamma[CohEN].at(5).push_back(-PhotoPol);
	MEASPhip[CohEN].at(5).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= 0.28 && KaonCosThetaCM <= 0.38){
	MEASGamma[CohEN].at(6).push_back(-PhotoPol);
	MEASPhip[CohEN].at(6).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= 0.38 && KaonCosThetaCM <= 0.46){
	MEASGamma[CohEN].at(7).push_back(-PhotoPol);
	MEASPhip[CohEN].at(7).push_back(PhiCM);
      }
      else if (KaonCosThetaCM >= 0.46 && KaonCosThetaCM <= 0.54){
	MEASGamma[CohEN].at(8).push_back(-PhotoPol);
	MEASPhip[CohEN].at(8).push_back(PhiCM);
      }
      else if (KaonCosThetaCM > 0.54){
	MEASGamma[CohEN].at(9).push_back(-PhotoPol);
	MEASPhip[CohEN].at(9).push_back(PhiCM);
      }
    }


  }
  

  h_Asym1[0]=(TH1F*)h_kaonPhiPA1[0]->GetAsymmetry(h_kaonPhiPE1[0]);
  h_Asym1[1]=(TH1F*)h_kaonPhiPA1[1]->GetAsymmetry(h_kaonPhiPE1[1]);
  h_Asym1[2]=(TH1F*)h_kaonPhiPA1[2]->GetAsymmetry(h_kaonPhiPE1[2]);
  h_Asym1[3]=(TH1F*)h_kaonPhiPA1[3]->GetAsymmetry(h_kaonPhiPE1[3]);
  h_Asym1[4]=(TH1F*)h_kaonPhiPA1[4]->GetAsymmetry(h_kaonPhiPE1[4]);
  h_Asym1[5]=(TH1F*)h_kaonPhiPA1[5]->GetAsymmetry(h_kaonPhiPE1[5]);
  h_Asym1[6]=(TH1F*)h_kaonPhiPA1[6]->GetAsymmetry(h_kaonPhiPE1[6]);
  h_Asym1[7]=(TH1F*)h_kaonPhiPA1[7]->GetAsymmetry(h_kaonPhiPE1[7]);
  h_Asym1[8]=(TH1F*)h_kaonPhiPA1[8]->GetAsymmetry(h_kaonPhiPE1[8]);
  h_Asym1[9]=(TH1F*)h_kaonPhiPA1[9]->GetAsymmetry(h_kaonPhiPE1[9]);

  h_Asym2[0]=(TH1F*)h_kaonPhiPA2[0]->GetAsymmetry(h_kaonPhiPE2[0]);
  h_Asym2[1]=(TH1F*)h_kaonPhiPA2[1]->GetAsymmetry(h_kaonPhiPE2[1]);
  h_Asym2[2]=(TH1F*)h_kaonPhiPA2[2]->GetAsymmetry(h_kaonPhiPE2[2]);
  h_Asym2[3]=(TH1F*)h_kaonPhiPA2[3]->GetAsymmetry(h_kaonPhiPE2[3]);
  h_Asym2[4]=(TH1F*)h_kaonPhiPA2[4]->GetAsymmetry(h_kaonPhiPE2[4]);
  h_Asym2[5]=(TH1F*)h_kaonPhiPA2[5]->GetAsymmetry(h_kaonPhiPE2[5]);
  h_Asym2[6]=(TH1F*)h_kaonPhiPA2[6]->GetAsymmetry(h_kaonPhiPE2[6]);
  h_Asym2[7]=(TH1F*)h_kaonPhiPA2[7]->GetAsymmetry(h_kaonPhiPE2[7]);
  h_Asym2[8]=(TH1F*)h_kaonPhiPA2[8]->GetAsymmetry(h_kaonPhiPE2[8]);
  h_Asym2[9]=(TH1F*)h_kaonPhiPA2[9]->GetAsymmetry(h_kaonPhiPE2[9]);

  h_Asym3[0]=(TH1F*)h_kaonPhiPA3[0]->GetAsymmetry(h_kaonPhiPE3[0]);
  h_Asym3[1]=(TH1F*)h_kaonPhiPA3[1]->GetAsymmetry(h_kaonPhiPE3[1]);
  h_Asym3[2]=(TH1F*)h_kaonPhiPA3[2]->GetAsymmetry(h_kaonPhiPE3[2]);
  h_Asym3[3]=(TH1F*)h_kaonPhiPA3[3]->GetAsymmetry(h_kaonPhiPE3[3]);
  h_Asym3[4]=(TH1F*)h_kaonPhiPA3[4]->GetAsymmetry(h_kaonPhiPE3[4]);
  h_Asym3[5]=(TH1F*)h_kaonPhiPA3[5]->GetAsymmetry(h_kaonPhiPE3[5]);
  h_Asym3[6]=(TH1F*)h_kaonPhiPA3[6]->GetAsymmetry(h_kaonPhiPE3[6]);
  h_Asym3[7]=(TH1F*)h_kaonPhiPA3[7]->GetAsymmetry(h_kaonPhiPE3[7]);
  h_Asym3[8]=(TH1F*)h_kaonPhiPA3[8]->GetAsymmetry(h_kaonPhiPE3[8]);
  h_Asym3[9]=(TH1F*)h_kaonPhiPA3[9]->GetAsymmetry(h_kaonPhiPE3[9]);


  DoCanvasAsym();
}

#endif
