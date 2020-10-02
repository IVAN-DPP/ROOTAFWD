#include "./DataEvent.cpp"
#include <iostream>
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


void Analysis(){
  string treeName="g13b";
  //For 1 root file use constructor below
  //string fileNamePERP="/home/emuneva/Analysis/PARA/skim2_55065.root";
  //DataEvent *myDataPERP=new DataEvent(fileNamePERP,treeName);
  string fileNamePERP="List1.txt";
  DataEvent *myDataPERP=new DataEvent(fileNamePERP,treeName, 35);
  
  ///Define histograms here
  TH2F *h_DeltaBe[3];
  h_DeltaBe[0]=new TH2F("h_DeltaBe_0","Proton ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  h_DeltaBe[1]=new TH2F("h_DeltaBe_1","Kaon ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.15, 0.15);
  h_DeltaBe[2]=new TH2F("h_DeltaBe_2","Pion ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);

  TH2F *h_DeltaBe2[3];
  h_DeltaBe2[0]=new TH2F("h_DeltaBe2_0","Proton ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  h_DeltaBe2[1]=new TH2F("h_DeltaBe2_1","Kaon ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.15, 0.15);
  h_DeltaBe2[2]=new TH2F("h_DeltaBe2_2","Pion ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);

  TH2F *h_DeltaBe3[3];
  h_DeltaBe3[0]=new TH2F("h_DeltaBe3_0","Proton ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  h_DeltaBe3[1]=new TH2F("h_DeltaBe3_1","Kaon ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.15, 0.15);
  h_DeltaBe3[2]=new TH2F("h_DeltaBe3_2","Pion ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);

  TH2F *h_DeltaBe4[3];
  h_DeltaBe4[0]=new TH2F("h_DeltaBe4_0","Proton ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  h_DeltaBe4[1]=new TH2F("h_DeltaBe4_1","Kaon ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.15, 0.15);
  h_DeltaBe4[2]=new TH2F("h_DeltaBe4_2","Pion ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);

  TH2F *h_DeltaBe5[3];
  h_DeltaBe5[0]=new TH2F("h_DeltaBe5_0","Proton ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);
  h_DeltaBe5[1]=new TH2F("h_DeltaBe5_1","Kaon ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.15, 0.15);
  h_DeltaBe5[2]=new TH2F("h_DeltaBe5_2","Pion ;p [GeV/c];#Delta #beta;",200, 0, 3, 200, -0.2, 0.2);

  
  TH2F *h_BeVSp[3];
  h_BeVSp[0]=new TH2F("h_BeVSp_0","Proton ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp[1]=new TH2F("h_BeVSp_1","Kaon ;p [GeV/c];#beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp[2]=new TH2F("h_BeVSp_2","Pion ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);

  TH2F *h_BeVSp2[3];
  h_BeVSp2[0]=new TH2F("h_BeVSp2_0","Proton ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp2[1]=new TH2F("h_BeVSp2_1","Kaon ;p [GeV/c];#beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp2[2]=new TH2F("h_BeVSp2_2","Pion ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  
  TH2F *h_BeVSp3[3];
  h_BeVSp3[0]=new TH2F("h_BeVSp3_0","Proton ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp3[1]=new TH2F("h_BeVSp3_1","Kaon ;p [GeV/c];#beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp3[2]=new TH2F("h_BeVSp3_2","Pion ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  
  TH2F *h_BeVSp4[3];
  h_BeVSp4[0]=new TH2F("h_BeVSp4_0","Proton ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp4[1]=new TH2F("h_BeVSp4_1","Kaon ;p [GeV/c];#beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp4[2]=new TH2F("h_BeVSp4_2","Pion ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  
  TH2F *h_BeVSp5[3];
  h_BeVSp5[0]=new TH2F("h_BeVSp5_0","Proton ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp5[1]=new TH2F("h_BeVSp5_1","Kaon ;p [GeV/c];#beta;",200, 0, 3, 200, 0, 1);
  h_BeVSp5[2]=new TH2F("h_BeVSp5_2","Pion ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);


  TH1F *h_eloss[3];
  h_eloss[0]=new TH1F("h_eloss_0","Proton; eloss",50, 0, 100);
  h_eloss[1]=new TH1F("h_eloss_1","Kaon; eloss",50, 0, 20);
  h_eloss[2]=new TH1F("h_eloss_2","Pion; eloss",50, 0, 20);

  
  TH2F *h_ThVSphi[3];
  h_ThVSphi[0]=new TH2F("h_ThVSphi_0","Proton ;p [GeV/c]; #beta;",200, -3.4, 3.4, 200, 0, 3);
  h_ThVSphi[1]=new TH2F("h_ThVSphi_1","Kaon ;p [GeV/c];#beta;",200, -3.4, 3.4, 200, 0, 3);
  h_ThVSphi[2]=new TH2F("h_ThVSphi_2","Pion ;p [GeV/c]; #beta;",200, -3.4, 3.4, 200, 0, 3);
  
  TH1F *h_CosTheta[3];
  h_CosTheta[0]=new TH1F("h_CosTheta_0","Proton ;p [GeV/c];#Delta #beta;",200, -1, 1);
  h_CosTheta[1]=new TH1F("h_CosTheta_1","Kaon ;p [GeV/c];#Delta #beta;",200, -1, 1);
  h_CosTheta[2]=new TH1F("h_CosTheta_2","Pion ;p [GeV/c];#Delta #beta;",200, -1, 1);
  
  TH1F *h_Theta[3];
  h_Theta[0]=new TH1F("h_Theta_0","Proton ;p [GeV/c];#Delta #beta;",200, 0, 3);
  h_Theta[1]=new TH1F("h_Theta_1","Kaon ;p [GeV/c];#Delta #beta;",200, 0, 3);
  h_Theta[2]=new TH1F("h_Theta_2","Pion ;p [GeV/c];#Delta #beta;",200, 0, 3);
  
  TH1F *h_phi[3];
  h_phi[0]=new TH1F("h_phi_0","Proton ;p [GeV/c]; #beta;",200, -3.4, 3.4);
  h_phi[1]=new TH1F("h_phi_1","Kaon ;p [GeV/c];#beta;",200, -3.4, 3.4);
  h_phi[2]=new TH1F("h_phi_2","Pion ;p [GeV/c]; #beta;",200, -3.4, 3.4);
  
  TH1F *h_DBe[3];
  h_DBe[0]=new TH1F("h_DBe_0","Proton ;p [GeV/c];#Delta #beta;",200, -0.2, 0.2);
  h_DBe[1]=new TH1F("h_DBe_1","Kaon ;p [GeV/c];#Delta #beta;",200, -0.2, 0.2);
  h_DBe[2]=new TH1F("h_DBe_2","Pion ;p [GeV/c];#Delta #beta;",200, -0.2, 0.2);
  
  TH1F *h_p[3];
  h_p[0]=new TH1F("h_p_0","Proton ;p [GeV/c]; #beta;",200, 0, 3);
  h_p[1]=new TH1F("h_p_1","Kaon ;p [GeV/c];#beta;",200, 0, 3);
  h_p[2]=new TH1F("h_p_2","Pion ;p [GeV/c]; #beta;",200, 0, 3);
  
  TH1F *h_DeltaTall[2];
  h_DeltaTall[0]=new TH1F("h_DeltaTall_0","Kaon ;#Delta t [ns];counts;", 200, -10, 10);
  h_DeltaTall[1]=new TH1F("h_DeltaTall_1","Pion ;#Delta t [ns];counts;", 200, -10, 10);
  
  TH2F *h_DeltaTallvsp[2];
  h_DeltaTallvsp[0]=new TH2F("h_DeltaTallvsp_0","Kaon ;#Delta t [ns];counts;", 200, 0, 3, 200, -10, 10);
  h_DeltaTallvsp[1]=new TH2F("h_DeltaTallvsp_1","Pion ;#Delta t [ns];counts;", 200, 0, 3, 200, -10, 10);


  TH1F *h_DeltaT[1];
  h_DeltaT[0]=new TH1F("h_DeltaT_0","Kaon ;#Delta t [ns];counts;", 200, -10, 10);

  TH1F *h_MissingMass = new TH1F("h_missingmass","Missing mass;", 100, 0.7, 1.2);
  TH1F *h_MissingMass_kaonpion = new TH1F("h_missingmass_kaonpion","Missing mass Kaon_pion;", 100, 0.7, 1.2);

  TH2F *h_MissingMass_vsMissingMasskaonpion = new TH2F("MissingMass_correlation", "",100, 0.7, 1.2, 100, 0.7, 1.2);


  TH1F *h_MissingP = new TH1F("h_missingp","Missing momentum;", 100, 0.0, 1.5);
  TH2F *h_MissingPvsMass = new TH2F("h_missingpvsm","Missing;", 100, 0.0, 1.5, 100, 0.7, 1.2);
  TH2F *h_MissingMassvsSigmaMass = new TH2F("h_missingmvsSigma","Missing;",100,1.0,1.5, 100, 0.7, 1.2);
  TH2F *h_MissingPvsSigmaMass = new TH2F("h_missingPvsSigma","Missing;",100,1.0,1.5, 100, 0.0, 1.5);
 
  TH1F *h_InvariantMass = new TH1F("h_InvariantMass","Invariant mass;", 100, 1.0, 1.5);
  TH1F *h_LambdaMass = new TH1F("h_LambdaMass","Invariant mass;", 100, 1.08, 1.16);

  
  //Analysis
  while (myDataPERP->getEntry()<myDataPERP->getEntries()){
    myDataPERP->getNextEntry();
    if (myDataPERP->getEntry() % 1000 == 0){
      fprintf (stderr, "Looped over PERP %.2f percent \r", myDataPERP->getEntry()*100.0/myDataPERP->getEntries());
      fflush (stderr);
    }
    
    double deltbeta[3];
    for (int i=0;i<myDataPERP->getNum_chargedtracks();i++){
      deltbeta[i]=myDataPERP->getEVNT_track(i).Beta()-myDataPERP->getEVNT_bem(i);
      h_DeltaBe[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),deltbeta[i]);
      h_BeVSp[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),myDataPERP->getEVNT_bem(i));
      //h_ThVSphi[i]->Fill(myDataPERP->getEVNT_track(i).Phi(),myDataPERP->getEVNT_track(i).Theta());
      //h_CosTheta[i]->Fill(myDataPERP->getEVNT_track(i).CosTheta());
      //h_Theta[i]->Fill(myDataPERP->getEVNT_track(i).Theta());
      //h_phi[i]->Fill(myDataPERP->getEVNT_track(i).Phi());
      //h_DBe[i]->Fill(deltbeta[i]);
      //h_p[i]->Fill(myDataPERP->getEVNT_track(i).Rho());
    }

    //Cuts
    //
    //z-vertex for K+
    //cout << "Vertz1: " << myDataPERP->getEVNT_vertex(1).Z() << endl;
    if(myDataPERP->getEVNT_vertex(1).Z()<-39.0  || myDataPERP->getEVNT_vertex(1).Z()>-1.0) continue;
    for (int i=0;i<myDataPERP->getNum_chargedtracks();i++){
      h_DeltaBe2[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),deltbeta[i]);
      h_BeVSp2[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),myDataPERP->getEVNT_bem(i));
    }
    //cout << "Vertz2: " << myDataPERP->getEVNT_vertex(1).Z() << endl;
    
    if( (deltbeta[0] > 0.05) || (deltbeta[0] < -0.045) ) continue;
    for (int i=0;i<myDataPERP->getNum_chargedtracks();i++){
      h_DeltaBe3[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),deltbeta[i]);
      h_BeVSp3[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),myDataPERP->getEVNT_bem(i));
    }
    
    if( (deltbeta[2] > 0.05) || (deltbeta[2] < -0.05) ) continue;
    for (int i=0;i<myDataPERP->getNum_chargedtracks();i++){
      h_DeltaBe4[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),deltbeta[i]);
      h_BeVSp4[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),myDataPERP->getEVNT_bem(i));
    }

    if( (deltbeta[1] > 0.025) || (deltbeta[1] < -0.025) ) continue;
    for (int i=0;i<myDataPERP->getNum_chargedtracks();i++){
      h_DeltaBe5[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),deltbeta[i]);
      h_BeVSp5[i]->Fill(myDataPERP->getEVNT_track(i).Rho(),myDataPERP->getEVNT_bem(i));
    }

    
    
    for (int ii=0;ii<myDataPERP->getNum_photons();ii++){
      //LOOP OVER ALL PHOTONS
      if (fabs(deltbeta[1])<0.025){
	h_DeltaTall[0]->Fill(myDataPERP->getDelt_t_k(ii));
	h_DeltaTallvsp[0]->Fill(myDataPERP->getEVNT_track(1).Rho(), myDataPERP->getDelt_t_k(ii));
      }
      if (fabs(deltbeta[2])<0.05){
	h_DeltaTall[1]->Fill(myDataPERP->getDelt_t_pi(ii));
	h_DeltaTallvsp[1]->Fill(myDataPERP->getEVNT_track(2).Rho(), myDataPERP->getDelt_t_pi(ii));
      }
    }

    if (myDataPERP->getNumph_k()==1)
      h_DeltaT[0]->Fill(myDataPERP->getDelt_t_k(myDataPERP->getIndex_k(0)));

    if (myDataPERP->getNumph_k()!=1) continue;
    

    //cout << "P1: " << myDataPERP->getEVNT_track(0).Rho() << endl;
    //cout << "P2: " << myDataPERP->geteloss_track(0).Rho() << endl;

 
    
    h_eloss[0]->Fill(1000.0*(myDataPERP->geteloss_track(0).Rho() - myDataPERP->getEVNT_track(0).Rho() ));
    h_eloss[1]->Fill(1000.0*(myDataPERP->geteloss_track(1).Rho() - myDataPERP->getEVNT_track(1).Rho() ));
    h_eloss[2]->Fill(1000.0*(myDataPERP->geteloss_track(2).Rho() - myDataPERP->getEVNT_track(2).Rho() ));

    
    //Missing mass
    TLorentzVector photon, deuteron, kaon, kaonpion, proton, pion, W, W2, Sigma, Lambda;
    photon.SetXYZM(0,0,myDataPERP->getTAGR_epho(myDataPERP->getIndex_k(0)),0);
    deuteron.SetXYZM(0,0,0,1.8756);
    double Px_kaonpion = myDataPERP->getEVNT_track(1).Rho()* sin(myDataPERP->getEVNT_track(1).Theta())* cos(myDataPERP->getEVNT_track(1).Phi());
    double Py_kaonpion = myDataPERP->getEVNT_track(1).Rho()* sin(myDataPERP->getEVNT_track(1).Theta())* sin(myDataPERP->getEVNT_track(1).Phi());
    double Pz_kaonpion = myDataPERP->getEVNT_track(1).Rho()*(myDataPERP->getEVNT_track(1).CosTheta());
    kaonpion.SetXYZM(Px_kaonpion, Py_kaonpion, Pz_kaonpion, 0.139);

    //double E_kaonpion = sqrt(
    //kaonpion.SetPxPyPzE(Px_kaonpion, Py_kaonpion, Pz_kaonpion, 0.139);


    //cout << "Px_kpi: " << Px_kaonpion << endl;
    //cout << "Py_kpi: " << Py_kaonpion << endl;
    //cout << "Pz_kpi: " << Pz_kaonpion << endl;
    //cout << endl;
    
    //proton = myDataPERP->getEVNT_track(0);
    //kaon = myDataPERP->getEVNT_track(1);
    //pion = myDataPERP->getEVNT_track(2);
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
  
  TCanvas *c0=new TCanvas("c0","My plots", 900, 500);
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


  TCanvas *c1=new TCanvas("c1","My plots", 900, 500);
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


  TCanvas *c3=new TCanvas("c3","My plots", 900, 500);
  c3->Divide(3,1);
  c3->cd(1);
  h_MissingPvsMass->Draw("colz");
  c3->cd(2);
  h_MissingMassvsSigmaMass->Draw("colz");
  c3->cd(3);
  h_MissingPvsSigmaMass->Draw("colz");

  TCanvas *c2=new TCanvas("c2","My plots", 900, 500);
  c2->Divide(2,3);
  c2->cd(1);
  h_MissingMass->Draw();
  c2->cd(2);
  h_MissingP->Draw();
  c2->cd(3);
  h_InvariantMass->Draw();
  c2->cd(4);
  h_LambdaMass->Draw();
  c2->cd(5);
  h_MissingMass_kaonpion->Draw();
  c2->cd(6);
  h_MissingMass_vsMissingMasskaonpion->Draw();

  
  TCanvas *c4=new TCanvas("c4","My plots", 900, 500);
  c4->Divide(3,1);
  c4->cd(1);
  h_eloss[0]->Draw();
  c4->cd(2);
  h_eloss[1]->Draw();
  c4->cd(3);
  h_eloss[2]->Draw();


  
  /*
  TCanvas *c1=new TCanvas("c1","My plots", 900, 500);
  c1->Divide(2,3);
  c1->cd(1);
  h_DeltaBe2[0]->Draw("colz");
  c1->cd(3);
  h_DeltaBe2[1]->Draw("colz");
  c1->cd(5);
  h_DeltaBe2[2]->Draw("colz");
  
  c1->cd(2);
  h_BeVSp2[0]->Draw("colz");
  c1->cd(4);
  h_BeVSp2[1]->Draw("colz");
  c1->cd(6);
  h_BeVSp2[2]->Draw("colz");

  TCanvas *c2=new TCanvas("c2","My plots", 900, 500);
  c2->Divide(2,3);
  c2->cd(1);
  h_DeltaBe3[0]->Draw("colz");
  c2->cd(3);
  h_DeltaBe3[1]->Draw("colz");
  c2->cd(5);
  h_DeltaBe3[2]->Draw("colz");
  
  c2->cd(2);
  h_BeVSp3[0]->Draw("colz");
  c2->cd(4);
  h_BeVSp3[1]->Draw("colz");
  c2->cd(6);
  h_BeVSp3[2]->Draw("colz");


  TCanvas *c3=new TCanvas("c3","My plots", 900, 500);
  c3->Divide(2,3);
  c3->cd(1);
  h_DeltaBe4[0]->Draw("colz");
  c3->cd(3);
  h_DeltaBe4[1]->Draw("colz");
  c3->cd(5);
  h_DeltaBe4[2]->Draw("colz");
  
  c3->cd(2);
  h_BeVSp4[0]->Draw("colz");
  c3->cd(4);
  h_BeVSp4[1]->Draw("colz");
  c3->cd(6);
  h_BeVSp4[2]->Draw("colz");
  */
  /*
  TCanvas *c4=new TCanvas("c4","My plots", 900, 500);
  c4->Divide(2,3);
  c4->cd(1);
  h_DeltaBe5[0]->Draw("colz");
  c4->cd(3);
  h_DeltaBe5[1]->Draw("colz");
  c4->cd(5);
  h_DeltaBe5[2]->Draw("colz");
  
  c4->cd(2);
  h_BeVSp5[0]->Draw("colz");
  c4->cd(4);
  h_BeVSp5[1]->Draw("colz");
  c4->cd(6);
  h_BeVSp5[2]->Draw("colz");
  */

  
  /*
  TCanvas *c1=new TCanvas("c1","My plots 2", 900, 500);
  c1->Divide(2,3);
  c1->cd(1);
  h_Theta[0]->Draw();
  //h_ThVSphi[0]->Draw("colz");
  c1->cd(3);
  h_Theta[1]->Draw();
  //h_ThVSphi[1]->Draw("colz");
  c1->cd(5);
  h_Theta[2]->Draw();
  //h_ThVSphi[2]->Draw("colz");
  
  c1->cd(2);
  h_phi[0]->Draw();
  //h_CosTheta[0]->Draw();
  c1->cd(4);
  h_phi[1]->Draw();
  //h_CosTheta[1]->Draw();
  c1->cd(6);
  h_phi[2]->Draw();
  //h_CosTheta[2]->Draw();
  
  TCanvas *c2 = new TCanvas("c2","My plots 3", 900, 500);
  c2->Divide(2,3);
  c2->cd(1);
  h_DBe[0]->Draw();
  //h_ThVSphi[0]->Draw("colz");
  c2->cd(3);
  h_DBe[1]->Draw();
  //h_ThVSphi[1]->Draw("colz");
  c2->cd(5);
  h_DBe[2]->Draw();
  //h_ThVSphi[2]->Draw("colz");
  
  c2->cd(2);
  h_p[0]->Draw();
  //h_CosTheta[0]->Draw();
  c2->cd(4);
  h_p[1]->Draw();
  //h_CosTheta[1]->Draw();
  c2->cd(6);
  h_p[2]->Draw();
  //h_CosTheta[2]->Draw();
  */


}




/* DESCRIPTION OF METHODS
 There are three particles in each event. The first particle is a proton, the second is a positive kaon and the third is a negative pion
 You can access the following info
 TLorentzVector getEVNT_track(int i){return loc_EVNT_track->at(i);} //TLorentzVector for proton (i=0) kaon (i=1) and negative pion (i=2). The TlorentzVector has the nominal masses of the tracks.
 int getEVNT_q(int i){return loc_EVNT_q->at(i);} //charge for track i
 int getEVNT_scsec(int i){return loc_EVNT_scsec->at(i);} //SC sector for track i
 int getEVNT_scpad(int i){return loc_EVNT_scpad->at(i);} //SC paddle for track i
 int getEVNT_schit(int i){return loc_EVNT_schit->at(i);} //SC hit for track i
 int getEVNT_stsec(int i){return loc_EVNT_stsec->at(i);} //ST sector for track i
 int getEVNT_sthit(int i){return loc_EVNT_sthit->at(i);} //ST hit for track i
 int getTAGR_eid(int i){return loc_tagr_eid->at(i);} //TAGR eid for photon i
 int getTAGR_tid(int i){return loc_tagr_tid->at(i);} //TAGR tid for photon i
 int getTAGR_stat(int i){return loc_tagr_stat->at(i);} //TAGR status for photon i
 int getSTPB_sthid(int i){return loc_STPB_sthid->at(i);} //STBP hitd for track i
 int getSCPB_ScPdHt(int i){return loc_SCPB_ScPdHt->at(i);} //SCBP ScPdHt for track i
 float getEVNT_bem(int i){return loc_EVNT_bem->at(i);} //beta measured for track i
 float getEVNT_sc_t(int i){return loc_EVNT_sc_t->at(i);} //sc time for track i
 float getEVNT_sc_d(int i){return loc_EVNT_sc_d->at(i);} //sc d for track i
 float getEVNT_st_t(int i){return loc_EVNT_st_t->at(i);} //st time for track i
 float getEVNT_st_d(int i){return loc_EVNT_st_d->at(i);} //st d for track i
 float getEVNT_sc_e(int i){return loc_EVNT_sc_e->at(i);} //sc energy for track i
 float getTAGR_epho(int i){return loc_tagr_epho->at(i);} //TAGR epho for photon i
 float getTAGR_tpho(int i){return loc_tagr_tpho->at(i);} //TAGR tpho for photon i
 int getNum_photons(){return loc_tagr_tpho->size();} //number of photons
 int getTrip_flag(){return loc_trip_flag;} //trip flag
 int getNumofpart(){return loc_numofpart;} //number of particles (including neutrals)
 int getNum_pos(){return loc_num_pos;} //number of positive
 int getNum_chargedtracks(){return loc_EVNT_track->size();} //number of charged
 int getNum_neg(){return loc_num_neg;} //number of negative
 int getNum_neu(){return loc_num_neu;} //number of neutrals
 int getHEAD_eventnum(){return loc_head_eventnum;} //HEAD event number
 int getHEAD_runnum(){return loc_head_runnum;} //HEAD run number
 int getNum_deuterons(){return loc_num_deuterons;} //number of deuterons
 int getNum_protons(){return loc_num_protons;} //number of protons
 int getNum_poskaons(){return loc_num_poskaons;} //number of postive kaons
 int getNum_pospions(){return loc_num_pospions;} //number of positive pions
 int getNum_negkaons(){return loc_num_negkaons;} //number of negative kaons
 int getNum_negpions(){return loc_num_negpions;} ////number of negative pions
 float getCoh_edge(){return loc_coh_edge;} //Coherent Edge from EPICS
 float getBeam_en(){return loc_beam_en;} //Beam energy from EPICS
 float getCoh_edge_nom(){return loc_coh_edge_nom;} //Nominal Coherent Edge
 int getCoh_plan_db(){return loc_coh_plan_db;} //Coherent plane from database
 int getCoh_radi(){return loc_coh_radi;} //Coherent radiator from EPICS
 int getCoh_plan(){return loc_coh_plan;} //Coherent radiator from EPICS
 float getDelt_t_k(int i){return loc_delt_t_k->at(i);} //Photon i coincidence time with kaon
 float getDelt_t_pi(int i){return loc_delt_t_pi->at(i);} //Photon i coincidence time with pion
 int getNumph_k(){return loc_numph_k;} //Number of photons within 1ns when looking at the photon-kaon coincidence time
 int getNumph_pi(){return loc_numph_pi;} //Number of photons within 1ns when looking at the photon-pion coincidence time
 int getIndex_k(int i){return loc_index_k->at(i);} //Photon index when sorted using the photon-kaon coincidence time
 int getIndex_pi(int ip){return loc_index_pi->at(i);} //Photon index when sorted using the photon-pion coincidence time. getIndex_pi(0) returns the photon position that produces the smallest coincidence time
 TLorentzVector getKpi_mm(int i){return loc_kpi_mm->at(i);} // 4-vector g n ->K+ pi- X
 TLorentzVector getK_mm(int i){return loc_k_mm->at(i);} // 4-vector g n ->K+ X
 TLorentzVector getPpi_mm(int i){return loc_ppi_mm->at(i);} // 4-vector g n ->p pi- X when kaon is given proton mass
 TLorentzVector getPipi_mm(int i){return loc_pipi_mm->at(i);} // 4-vector g n ->pi+ pi- X when kaon is given pion mass
 TLorentzVector getDKppi_mm(int i){return loc_d_kppi_mm->at(i);} // 4-vector g d ->p K+ pi- X
 TLorentzVector getDKp_mm(int i){return loc_d_kp_mm->at(i);} // 4-vector g d ->p K+ X
 TVector3 getEVNT_vertex(int i){return loc_EVNT_vertex->at(i);} //EVNT vertex of i track
 TVector3 getMVRT_vertex(){return *loc_MVRT_vertex;} //MVRT vertex
 int getNextEntry();
 int getEntry(){return eventno;}
 int getEntries();

*/
