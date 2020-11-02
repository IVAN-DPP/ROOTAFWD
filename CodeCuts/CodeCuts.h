#include "Histograms.h"
#include "Libraries.h"

#ifndef CODECUTS_H
#define CODECUTS_H

class Codecuts : public Histograms {

private:

public:
  Codecuts(){DoHistograms();}
  void CodeCuts();

};

void Codecuts::CodeCuts(){


  //RooRealVar *x = new RooRealVar("x","x",-10,-10); // Pendiente --p-p-p-p-p-p-p-p-p-para
  
  string treeName="g13b";
  //For 1 root file use constructor below
  //string fileNamePERP="/home/emuneva/Analysis/PARA/skim2_55065.root";
  //DataEvent *myDataPERP=new DataEvent(fileNamePERP,treeName);
  string fileNamePERP="../List1.txt";
  DataEvent *myDataPERP=new DataEvent(fileNamePERP,treeName, 35);


  while (myDataPERP->getEntry()<myDataPERP->getEntries()){
    myDataPERP->getNextEntry();
    if (myDataPERP->getEntry() % 1000 == 0){
      fprintf (stderr, "Looped over PERP %.2f percent \r", myDataPERP->getEntry()*100.0/myDataPERP->getEntries());
      fflush (stderr);
    }

    //------------------ Delta Beta ---------------//
    double deltbeta[3];
    double deltbetacut[3];
    
    if(myDataPERP->getEVNT_vertex(1).Z()<-39.0  || myDataPERP->getEVNT_vertex(1).Z()>-1.0) continue;
    for (int i=0;i<myDataPERP->getNum_chargedtracks();i++){
      deltbeta[i]=myDataPERP->getEVNT_track(i).Beta()-myDataPERP->getEVNT_bem(i);
      h_DeltaBe[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),deltbeta[i]);
      h_BeVSp[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),myDataPERP->getEVNT_bem(i));
    }

    //------------------ Delta Beta with Cuts ---------------//
    
    for (int i=0;i<myDataPERP->getNum_chargedtracks();i++){
      deltbetacut[i]=myDataPERP->getEVNT_track(i).Beta()-myDataPERP->getEVNT_bem(i);
      if(deltbetacut[0] > 0.05  || deltbetacut[0] < -0.045) continue;
      if(deltbetacut[1] > 0.025 || deltbetacut[1] < -0.025) continue;
      if(deltbetacut[2] > 0.05  || deltbetacut[2] < -0.05) continue;
      h_DeltaBecut[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),deltbeta[i]);
      h_BeVSpcut[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),myDataPERP->getEVNT_bem(i));
    }

    if(deltbetacut[2] > 0.05  || deltbetacut[2] < -0.05) continue;    //Cut from Delta Beta vs Missingmass,Missing momentum, Invariantmass


    
   
    //------------------Correlation Theta-Phi, -----------------------/
    
    h_ThePhi_proton->Fill(myDataPERP->getEVNT_track(0).Phi()*TMath::RadToDeg(), myDataPERP->getEVNT_track(0).Theta()*TMath::RadToDeg());  
    h_ThePhi_kaon->Fill(myDataPERP->getEVNT_track(1).Phi()*TMath::RadToDeg(), myDataPERP->getEVNT_track(1).Theta()*TMath::RadToDeg());
    h_ThePhi_pion->Fill(myDataPERP->getEVNT_track(2).Phi()*TMath::RadToDeg(), myDataPERP->getEVNT_track(2).Theta()*TMath::RadToDeg());

    
     //-----Cuts due to detector geometry (Fiduciary cuts)---------//
    //------------------------------------------------------------//

    Double_t phiproton_cut;
    phiproton_cut= myDataPERP->getEVNT_track(0).Phi()*TMath::RadToDeg();

    if( (phiproton_cut <= -25)   && (phiproton_cut >= -35)  ) continue;
    if( (phiproton_cut <= -85)   && (phiproton_cut >= -95)  ) continue;
    if( (phiproton_cut <= -145)  && (phiproton_cut >= -155) ) continue;
    if( (phiproton_cut >= 25)    && (phiproton_cut <= 35)   ) continue;
    if( (phiproton_cut >= 85)    && (phiproton_cut <= 95)   ) continue;
    if( (phiproton_cut >= 145)   && (phiproton_cut <= 155)  ) continue;
     
    h_ThePhi_protoncut->Fill(phiproton_cut, myDataPERP->getEVNT_track(0).Theta()*TMath::RadToDeg());
    
    //kaon cuts
    
    Double_t phikaon_cut;
    phikaon_cut= myDataPERP->getEVNT_track(1).Phi()*TMath::RadToDeg();
    
    if( (phikaon_cut <= -25)   && (phikaon_cut >= -35)  ) continue;
    if( (phikaon_cut <= -85)   && (phikaon_cut >= -95)  ) continue;
    if( (phikaon_cut <= -145)  && (phikaon_cut >= -155) ) continue;
    if( (phikaon_cut >= 25)    && (phikaon_cut <= 35)   ) continue;
    if( (phikaon_cut >= 85)    && (phikaon_cut <= 95)   ) continue;
    if( (phikaon_cut >= 145)   && (phikaon_cut <= 155)  ) continue;
    
    h_ThePhi_kaoncut->Fill(phikaon_cut, myDataPERP->getEVNT_track(1).Theta()*TMath::RadToDeg());
    
    //Pion cuts
    
    Double_t phiPion_cut;
    phiPion_cut= myDataPERP->getEVNT_track(2).Phi()*TMath::RadToDeg();
    
    if( (phiPion_cut <= -25)   && (phiPion_cut >= -35)  ) continue;
    if( (phiPion_cut <= -85)   && (phiPion_cut >= -95)  ) continue;
    if( (phiPion_cut <= -145)  && (phiPion_cut >= -155) ) continue;
    if( (phiPion_cut >= 25)    && (phiPion_cut <= 35)   ) continue;
    if( (phiPion_cut >= 85)    && (phiPion_cut <= 95)   ) continue;
    if( (phiPion_cut >= 145)   && (phiPion_cut <= 155)  ) continue;
    
    h_ThePhi_pioncut->Fill(phiPion_cut, myDataPERP->getEVNT_track(2).Theta()*TMath::RadToDeg());
    
    //------------------ Photons, Delta T  ------------------ // 
    
    for (int i=0;i<myDataPERP->getNum_photons();i++){
      //LOOP OVER ALL PHOTONS
      if (fabs(deltbeta[1])<0.025){
	h_DeltaTall[0]->Fill(myDataPERP->getDelt_t_k(i));
	h_DeltaTallvsp[0]->Fill(myDataPERP->getEVNT_track(1).Rho(),myDataPERP->getDelt_t_k(i));
      }
      if (fabs(deltbeta[2])<0.05){
	h_DeltaTall[1]->Fill(myDataPERP->getDelt_t_pi(i));
	h_DeltaTallvsp[1]->Fill(myDataPERP->getEVNT_track(2).Rho(),myDataPERP->getDelt_t_pi(i));
      }
    }
    
    //------------Delta T with Cuts ----------- //
    
    if (myDataPERP->getNumph_k()==1)
      h_DeltaT[0]->Fill(myDataPERP->getDelt_t_k(myDataPERP->getIndex_k(0)));
    if(myDataPERP->getNumph_pi()==1)
      h_DeltaT[1]->Fill(myDataPERP->getDelt_t_pi(myDataPERP->getIndex_pi(0)));
    
    if (myDataPERP->getNumph_k()!=1) continue;
    //if (myDataPERP->getNumph_pi()!=1) continue;

    //--------------- Energy loss ----------- //
    
    h_eloss[0]->Fill(1000.0*(myDataPERP->geteloss_track(0).Rho() - myDataPERP->getEVNT_track(0).Rho() ));
    h_eloss[1]->Fill(1000.0*(myDataPERP->geteloss_track(1).Rho() - myDataPERP->getEVNT_track(1).Rho() ));
    h_eloss[2]->Fill(1000.0*(myDataPERP->geteloss_track(2).Rho() - myDataPERP->getEVNT_track(2).Rho() ));


    //-------------- Reconstruction --------- //
    
    TLorentzVector photon, deuteron, kaon, kaonpion, proton, pion, Wneutron_kaon, Wneutron_pion, Sigma, Lambda, Neutron;
    photon.SetXYZM(0,0,myDataPERP->getTAGR_epho(myDataPERP->getIndex_k(0)),0);
    deuteron.SetXYZM(0,0,0,1.8756);
    double Px_kaonpion = myDataPERP->getEVNT_track(1).Rho()* sin(myDataPERP->getEVNT_track(1).Theta())* cos(myDataPERP->getEVNT_track(1).Phi());
    double Py_kaonpion = myDataPERP->getEVNT_track(1).Rho()* sin(myDataPERP->getEVNT_track(1).Theta())* sin(myDataPERP->getEVNT_track(1).Phi());
    double Pz_kaonpion = myDataPERP->getEVNT_track(1).Rho()*(myDataPERP->getEVNT_track(1).CosTheta());
    kaonpion.SetXYZM(Px_kaonpion, Py_kaonpion, Pz_kaonpion, 0.139);       // This mass is of Pion-, because we need remove the background of Pion-


    proton = myDataPERP->geteloss_track(0);
    kaon = myDataPERP->geteloss_track(1);
    pion = myDataPERP->geteloss_track(2);
    Wneutron_kaon = photon + deuteron - proton - kaon - pion;
    Neutron.SetXYZM(Wneutron_kaon.Px(), Wneutron_kaon.Py(), Wneutron_kaon.Pz(), 0.939);
    Sigma = pion + Neutron;
    Lambda = pion + proton;
    
    Wneutron_pion = photon + deuteron - proton - kaonpion - pion;         // This missing mass is with the Pion-
    

    double radx=0.04, rady=0.02, offsetx=0.9395601, offsety=1.052416, angle=45*TMath::DegToRad();
    
      Double_t El = TMath::Power((Wneutron_kaon.M()-offsetx)*cos(angle)+(Wneutron_pion.M()-offsety)*sin(angle),2)/TMath::Power(radx,2)+TMath::Power((Wneutron_kaon.M()-offsetx)*sin(angle)-(Wneutron_pion.M()-offsety)*cos(angle),2)/TMath::Power(rady,2);

    if(El > 1) continue;
    
    
    h_MissingMass->Fill(Wneutron_kaon.M());
    h_MissingMass_kaonpion->Fill(Wneutron_pion.M());
    
    //h_MissingP->Fill(W.P());
    //h_MissingPvsMass->Fill(W.P(), Wneutron_kaon.M());
    //h_MissingMassvsSigmaMass->Fill(Sigma.M(), Wneutron_kaon.M());
    //h_MissingPvsSigmaMass->Fill(Sigma.M(), Wneutron_kaon.P());
    h_MissingMass_vsMissingMasskaonpion->Fill(Wneutron_kaon.M(), Wneutron_pion.M());
    
    if( ((Lambda.M()<1.05) || (Lambda.M()>1.20)) && (Wneutron_kaon.M()>0.9) && (Wneutron_kaon.M()<1.0) && (Wneutron_kaon.P()>0.2))
    h_InvariantMass->Fill(Sigma.M());

    //if( ((Lambda.M()<1.05) || (Lambda.M()>1.20)) && (W.M()>0.9) && (W.M()<1.0) )
    //if( ((Lambda.M()<1.05) || (Lambda.M()>1.20)) && (W.M()>0.9) )
    
    
    //if(Sigma.M()<1.10 || Sigma.M()>1.40)
    h_LambdaMass->Fill(Lambda.M());

    // Pion - 
    h_DeltaBVSInvariantMass->Fill(Sigma.M(),deltbeta[2]);
    h_DeltaBVSMissingMass->Fill(Wneutron_kaon.M(),deltbeta[2]);
    h_DeltaBVSMissingMomentum->Fill(Wneutron_kaon.P(),deltbeta[2]);
    
  }
  cout<<endl;
  
  //TF1 *f1 = new TF1("f1","gaus",-0.1,0.1);
  //f1->SetParameters(0.,1.0);
  // h_MissingMass->Fit(f1);
   //f1->SetRange(0.1,1.7,-0.045,0.05);
  // f1->Draw("same");
  
   DoCanvas();
}

#endif
