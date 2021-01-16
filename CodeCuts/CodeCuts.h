#include "./Histograms.h"
#include "./include/Libraries.h"
#include "./include/Miscelaneous.h"

#ifndef CODECUTS_H
#define CODECUTS_H

class Codecuts : public Histograms {

private:

public:
  Codecuts(){DoHistograms();}
  void CodeCuts();

};

void Codecuts::CodeCuts(){

  
  string treeName="g13b";  
  string fileNamePERP="ListFiles.txt";
  DataEvent *myDataList=new DataEvent(fileNamePERP,treeName, 35);
  
  vector<string> PolTableName;
  ListFilesAtDir("./ListTables", PolTableName);

  
  
  const int NumbOfPolFiles=PolTableName.size();
  for (int i=0;i<NumbOfPolFiles;i++){
    char *cstr = const_cast<char *>(PolTableName[i].c_str());
    LoadPolTable(i,cstr);
  }

  while (myDataList->getEntry()<myDataList->getEntries()){

    myDataList->getNextEntry();
    if (myDataList->getEntry() % 1000 == 0){
      if(myDataList->getCoh_plan() == 0){
	fprintf (stderr, "Looped %s : %.2f percent \r", "PARA" ,myDataList->getEntry()*100.0/myDataList->getEntries());
	fflush (stderr);
      }
      else {
	fprintf (stderr, "Looped %s : %.2f percent \r", "PERP" ,myDataList->getEntry()*100.0/myDataList->getEntries());
	fflush (stderr);
      }
    }
      
      
      
    //------------------ Delta Beta ---------------//
    double deltbeta[3]     = {};
    double deltbetacut[3]  = {};
      
    //---------- Vertex ------------ //
      
      
    h_Vertex->Fill(myDataList->getEVNT_vertex(1).Z());
    if(myDataList->getEVNT_vertex(1).Z()<-39.0  || myDataList->getEVNT_vertex(1).Z()>-1.0) continue;
      
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
    if(deltbetacut[0] > 0.02  || deltbetacut[0] < -0.02) continue;
    if(deltbetacut[1] > 0.025 || deltbetacut[1] < -0.025) continue;
    //Cut from Delta Beta vs Missingmass,Missing momentum, Invariantmass
    
    
    
    
   
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
      
    h_ThePhicut[0]->Fill(phiproton_cut, myDataList->getEVNT_track(0).Theta()*TMath::RadToDeg());
        
    //kaon cuts
    Double_t phikaon_cut;
    phikaon_cut= myDataList->getEVNT_track(1).Phi()*TMath::RadToDeg();
        
    if( (phikaon_cut <= -25)   && (phikaon_cut >= -35)  ) continue;
    if( (phikaon_cut <= -85)   && (phikaon_cut >= -95)  ) continue;
    if( (phikaon_cut <= -145)  && (phikaon_cut >= -155) ) continue;
    if( (phikaon_cut >= 25)    && (phikaon_cut <= 35)   ) continue;
    if( (phikaon_cut >= 85)    && (phikaon_cut <= 95)   ) continue;
    if( (phikaon_cut >= 145)   && (phikaon_cut <= 155)  ) continue;

    h_ThePhicut[1]->Fill(phikaon_cut, myDataList->getEVNT_track(1).Theta()*TMath::RadToDeg());

   
      
    //Pion cuts
    
    Double_t phiPion_cut;
    phiPion_cut= myDataList->getEVNT_track(2).Phi()*TMath::RadToDeg();
      
    if( (phiPion_cut <= -25)   && (phiPion_cut >= -35)  ) continue;
    if( (phiPion_cut <= -85)   && (phiPion_cut >= -95)  ) continue;
    if( (phiPion_cut <= -145)  && (phiPion_cut >= -155) ) continue;
    if( (phiPion_cut >= 25)    && (phiPion_cut <= 35)   ) continue;
    if( (phiPion_cut >= 85)    && (phiPion_cut <= 95)   ) continue;
    if( (phiPion_cut >= 145)   && (phiPion_cut <= 155)  ) continue;
      
    h_ThePhicut[2]->Fill(phiPion_cut, myDataList->getEVNT_track(2).Theta()*TMath::RadToDeg());
      
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
    //if (myDataList->getNumph_pi()!=1) continue;
      
    //--------------- Energy loss ----------- //
      
    h_eloss[0]->Fill(1000.0*(myDataList->geteloss_track(0).Rho() - myDataList->getEVNT_track(0).Rho() ));
    h_eloss[1]->Fill(1000.0*(myDataList->geteloss_track(1).Rho() - myDataList->getEVNT_track(1).Rho() ));
    h_eloss[2]->Fill(1000.0*(myDataList->geteloss_track(2).Rho() - myDataList->getEVNT_track(2).Rho() ));

    //--------------- Coh Edge -------------- //
    
      
      
    h_TagrEpho[0]->Fill(myDataList->getTAGR_epho(myDataList->getIndex_k(0))*1000.0);
    if (myDataList->getTAGR_epho(myDataList->getIndex_k(0))*1000.0 > myDataList->getCoh_edge()) continue;
    h_TagrEpho[1]->Fill(myDataList->getTAGR_epho(myDataList->getIndex_k(0))*1000.0);
    if (myDataList->getTAGR_epho(myDataList->getIndex_k(0))*1000< myDataList->getCoh_edge()-200.0) continue;
    h_TagrEpho[2]->Fill(myDataList->getTAGR_epho(myDataList->getIndex_k(0))*1000.0);
    if (fabs(myDataList->getCoh_edge()-myDataList->getCoh_edge_nom()*1000)>15)continue;
    if (myDataList->getTrip_flag()!=0)continue;
    if (myDataList->getCoh_plan()!=0 && myDataList->getCoh_plan()!=1)continue;
      


       
    double PhotoPol=GetPol(myDataList->getCoh_plan(), myDataList->getCoh_edge(), myDataList->getTAGR_epho(myDataList->getIndex_k(0))*1000.0, 8, 0.2,0.3); 
    if (PhotoPol<0.5) continue;

      

      
      
    //-------------- Reconstruction --------- //
      
    TLorentzVector photon, deuteron, kaon, kaonpion, proton, pion, Wneutron_kaon, Wneutron_pion, Sigma, Lambda, Neutron, WBoost;
    photon.SetXYZM(0,0,myDataList->getTAGR_epho(myDataList->getIndex_k(0)),0);
    deuteron.SetXYZM(0,0,0,1.8756);
    double Px_kaonpion = myDataList->getEVNT_track(1).Rho()* sin(myDataList->getEVNT_track(1).Theta())* cos(myDataList->getEVNT_track(1).Phi());
    double Py_kaonpion = myDataList->getEVNT_track(1).Rho()* sin(myDataList->getEVNT_track(1).Theta())* sin(myDataList->getEVNT_track(1).Phi());
    double Pz_kaonpion = myDataList->getEVNT_track(1).Rho()*(myDataList->getEVNT_track(1).CosTheta());
    kaonpion.SetXYZM(Px_kaonpion, Py_kaonpion, Pz_kaonpion, 0.139);       // This mass is of Pion-, because we need remove the background of Pion-
      
      
    proton = myDataList->geteloss_track(0);
    kaon = myDataList->geteloss_track(1);
    pion = myDataList->geteloss_track(2);
    Wneutron_kaon = photon + deuteron - proton - kaon - pion;
    Neutron.SetXYZM(Wneutron_kaon.Px(), Wneutron_kaon.Py(), Wneutron_kaon.Pz(), 0.939);
    Sigma = pion + Neutron;
    Lambda = pion + proton;
    WBoost = photon + deuteron; // to make Boost
      
    Wneutron_pion = photon + deuteron - proton - kaonpion - pion;         // This missing mass is with the Pion-
      
      
    h_MissingMass->Fill(Wneutron_kaon.M());
    h_MissingMass_kaonpion->Fill(Wneutron_pion.M());
    //h_MissingPvsMass->Fill(Wneutron_kaon.M(),Wneutron_kaon.P());
    h_MissingMass_vsMissingMasskaonpion[0]->Fill(Wneutron_kaon.M(), Wneutron_pion.M());

      
    //----------Correlación momentums vs missing mass------------------------//
      
    if( Lambda.M()<1.11 || Lambda.M()>1.132) 
      h_MissingPvsMass[0]->Fill(Wneutron_kaon.M(),Wneutron_kaon.P());
      
    if( Sigma.M()<1.08 || Sigma.M()>1.3)
      h_MissingPvsMass[1]->Fill(Wneutron_kaon.M(),Wneutron_kaon.P());
      
      
      
    Double_t El = TMath::Power((Wneutron_kaon.M()-offsetx)*cos(angle)+(Wneutron_pion.M()-offsety)*sin(angle),2)/TMath::Power(radx,2)
      +TMath::Power((Wneutron_kaon.M()-offsetx)*sin(angle)-(Wneutron_pion.M()-offsety)*cos(angle),2)/TMath::Power(rady,2);
      
    if(El > 1) continue;
      
      
    h_MissingMasscut->Fill(Wneutron_kaon.M());
    h_MissingMass_kaonpioncut->Fill(Wneutron_pion.M());
      
    // h_MissingP->Fill(Wneutron_kaon.P());
      
    h_MissingPvsSigmaMass->Fill(Sigma.M(), Wneutron_kaon.P());
    h_MissingMass_vsMissingMasskaonpion[1]->Fill(Wneutron_kaon.M(), Wneutron_pion.M());
          
    h_InvariantMass->Fill(Sigma.M());
    h_LambdaMass->Fill(Lambda.M());
      
      
    //------------------------Comparación de sigmas-------------//
    if( Lambda.M()<1.108 || Lambda.M()>1.124)
      h_InvariantMasscut[0]->Fill(Sigma.M());
      
    if( Lambda.M()<1.11 || Lambda.M()>1.132)
      h_InvariantMasscut[1]->Fill(Sigma.M());
      
    if( Lambda.M()<1.092 || Lambda.M()>1.14)
      h_InvariantMasscut[2]->Fill(Sigma.M());
      
      
    //------------- Comparación de missing momentums-----------//
      
    if( Lambda.M()<1.096 || Lambda.M()>1.136) 
      h_MissingP[0]->Fill(Wneutron_kaon.P());
      
      
      
    if( Sigma.M()<1.08 || Sigma.M()>1.3)
      h_MissingP[1]->Fill(Wneutron_kaon.P());
      
      
    if(Wneutron_kaon.P()<=0.2) continue;                                    //Cut for rescattering
      
      
    // h_MissingPcut->Fill(Wneutron_kaon.P());
      
      
    if ( Lambda.M()>=1.11 && Lambda.M()<=1.132) continue;                   //Cut for LamdaMass in +/- 3sigma
    if ( Wneutron_kaon.M()<=0.9 || Wneutron_kaon.M()>=0.96 ) continue;       //Cut from correlation MM
      
    h_InvariantMasscut[3]->Fill(Sigma.M());
    h_MissingMassvsSigmaMass->Fill(Sigma.M(), Wneutron_kaon.M());
      
    // Pion -  (Por si las moscas)
    h_DeltaBVSInvariantMass->Fill(Sigma.M(),deltbeta[2]);
    h_DeltaBVSMissingMass->Fill(Wneutron_kaon.M(),deltbeta[2]);
    h_DeltaBVSMissingMomentum->Fill(Wneutron_kaon.P(),deltbeta[2]);


   
    //-----------------BOOST------------------------------//

    TVector3 b=WBoost.BoostVector();
    kaon.Boost(-b);
    double KaonCosThetaCM=TMath::Cos(kaon.Theta());
    double KaonPhiCM=kaon.Phi()*TMath::RadToDeg();
    h_KCosThetaCM->Fill(KaonCosThetaCM);

    //---------------Bins Cos Theta Kaon----------------//

    //0 is for PARA
    //1 is for PERP
    if (myDataList->getCoh_plan()==0){
      if(KaonCosThetaCM < 0.668){
	MEASPhi.at(0).push_back(KaonPhiCM);
	MEASGammaP.at(0).push_back(PhotoPol);
	h_kaonPhiPA[0]->Fill(KaonPhiCM);
      }
      else if (KaonCosThetaCM > 0.668){
	MEASPhi.at(1).push_back(KaonPhiCM);
	MEASGammaP.at(1).push_back(PhotoPol);
	h_kaonPhiPA[1]->Fill(KaonPhiCM);
      }
    }
      
    else if (myDataList->getCoh_plan()==1){
      if(KaonCosThetaCM < 0.668){
	MEASPhi.at(0).push_back(KaonPhiCM);
	MEASGammaP.at(0).push_back(-PhotoPol);
	h_kaonPhiPE[0]->Fill(KaonPhiCM);
      }
      else if (KaonCosThetaCM > 0.668){
	MEASPhi.at(1).push_back(KaonPhiCM);
	MEASGammaP.at(1).push_back(-PhotoPol);
	h_kaonPhiPE[1]->Fill(KaonPhiCM);
      }
    }
     
    //---------------Asymmetry Analysis----------------//
      
    h_Asym[0]=(TH1F*)h_kaonPhiPA[0]->GetAsymmetry(h_kaonPhiPE[0]);
    h_Asym[1]=(TH1F*)h_kaonPhiPA[1]->GetAsymmetry(h_kaonPhiPE[1]);
      
      

      
  }
  cout<<endl;
    
  DoCanvas();
  
}

#endif
