#include "Histograms.h"
#include "../DataEvent.cpp"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TEllipse.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include <sstream>

#ifndef CODECUTS_H
#define CODECUTS_H

class Clase : public Histograms {

private:

public:
  Clase(){DoHistograms();}
  void CodeCuts();

};

void Clase::CodeCuts(){

  
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

    if(myDataPERP->getEVNT_vertex(1).Z()<-39.0  || myDataPERP->getEVNT_vertex(1).Z()>-1.0) continue;
    for (int i=0;i<myDataPERP->getNum_chargedtracks();i++){
      deltbeta[i]=myDataPERP->getEVNT_track(i).Beta()-myDataPERP->getEVNT_bem(i);
      if(deltbeta[0] > 0.05  || deltbeta[0] < -0.045) continue;
      if(deltbeta[2] > 0.05  || deltbeta[2] < -0.05) continue;
      if(deltbeta[1] > 0.025 || deltbeta[1] < -0.025) continue;
      h_DeltaBe[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),deltbeta[i]);
      h_BeVSp[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),myDataPERP->getEVNT_bem(i));
    }
    
    //------------------ Photons, Delta T  ------------------ // 
    
    for (int i=0;i<myDataPERP->getNum_photons();i++){
      //LOOP OVER ALL PHOTONS
      if (fabs(deltbeta[1])<0.025){
	h_DeltaTall[0]->Fill(myDataPERP->getDelt_t_k(i));
	h_DeltaTallvsp[0]->Fill(myDataPERP->getEVNT_track(1).Rho(), myDataPERP->getDelt_t_k(i));
      }
      if (fabs(deltbeta[2])<0.05){
	h_DeltaTall[1]->Fill(myDataPERP->getDelt_t_pi(i));
	h_DeltaTallvsp[1]->Fill(myDataPERP->getEVNT_track(2).Rho(), myDataPERP->getDelt_t_pi(i));
      }
    }
    
    //------------Delta T with Cuts ----------- //
    
    if (myDataPERP->getNumph_k()==1)
      h_DeltaT[0]->Fill(myDataPERP->getDelt_t_k(myDataPERP->getIndex_k(0)));
    
    if (myDataPERP->getNumph_k()!=1) continue;

    //--------------- Energy loss ----------- //
    
    h_eloss[0]->Fill(1000.0*(myDataPERP->geteloss_track(0).Rho() - myDataPERP->getEVNT_track(0).Rho() ));
    h_eloss[1]->Fill(1000.0*(myDataPERP->geteloss_track(1).Rho() - myDataPERP->getEVNT_track(1).Rho() ));
    h_eloss[2]->Fill(1000.0*(myDataPERP->geteloss_track(2).Rho() - myDataPERP->getEVNT_track(2).Rho() ));

    //-------------- Reconstruction --------- //
    
    TLorentzVector photon, deuteron, kaon, kaonpion, proton, pion, W, W2, Sigma, Lambda;
    photon.SetXYZM(0,0,myDataPERP->getTAGR_epho(myDataPERP->getIndex_k(0)),0);
    deuteron.SetXYZM(0,0,0,1.8756);
    double Px_kaonpion = myDataPERP->getEVNT_track(1).Rho()* sin(myDataPERP->getEVNT_track(1).Theta())* cos(myDataPERP->getEVNT_track(1).Phi());
    double Py_kaonpion = myDataPERP->getEVNT_track(1).Rho()* sin(myDataPERP->getEVNT_track(1).Theta())* sin(myDataPERP->getEVNT_track(1).Phi());
    double Pz_kaonpion = myDataPERP->getEVNT_track(1).Rho()*(myDataPERP->getEVNT_track(1).CosTheta());
    kaonpion.SetXYZM(Px_kaonpion, Py_kaonpion, Pz_kaonpion, 0.139);


    proton = myDataPERP->geteloss_track(0);
    kaon = myDataPERP->geteloss_track(1);
    pion = myDataPERP->geteloss_track(2);
    W = photon + deuteron - proton - kaon - pion;
    Sigma = pion + W;
    Lambda = pion + proton;
    
    W2 = photon + deuteron - proton - kaonpion - pion;
    
    if(W2.M() > 0.98)
      {
	h_MissingMass_kaonpion->Fill(W2.M());
	h_MissingMass_vsMissingMasskaonpion->Fill(W.M(), W2.M());
	h_MissingMass->Fill(W.M());
	h_MissingP->Fill(W.P());
	h_MissingPvsMass->Fill(W.P(), W.M());
	h_MissingMassvsSigmaMass->Fill(Sigma.M(), W.M());
	h_MissingPvsSigmaMass->Fill(Sigma.M(), W.P());
      }
    
    if( ((Lambda.M()<1.05) || (Lambda.M()>1.20)) && (W.M()>0.9) && (W.M()<1.0) && (W.P()>0.2))
      h_InvariantMass->Fill(Sigma.M());
    //if( ((Lambda.M()<1.05) || (Lambda.M()>1.20)) && (W.M()>0.9) && (W.M()<1.0) )
    //if( ((Lambda.M()<1.05) || (Lambda.M()>1.20)) && (W.M()>0.9) )
    
    
    if(Sigma.M()<1.10 || Sigma.M()>1.40)
      h_LambdaMass->Fill(Lambda.M());
    
    
  }
  cout<<endl;
  
  
  DoCanvas();
}

#endif
