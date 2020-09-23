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

    TH2F *h_BeVSp[3];
    h_BeVSp[0]=new TH2F("h_BeVSp_0","Proton ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);
    h_BeVSp[1]=new TH2F("h_BeVSp_1","Kaon ;p [GeV/c];#beta;",200, 0, 3, 200, 0, 1);
    h_BeVSp[2]=new TH2F("h_BeVSp_2","Pion ;p [GeV/c]; #beta;",200, 0, 3, 200, 0, 1);

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
    

    TH1F *h_phot[3];
    h_phot[0]=new TH1F("h_phot_0","Proton ;p [GeV/c]; #beta;",200, 0, 3);
    h_phot[1]=new TH1F("h_phot_1","Kaon ;p [GeV/c];#beta;",200, 0, 3);
    h_phot[2]=new TH1F("h_phot_2","Pion ;p [GeV/c]; #beta;",200, 0, 3);


    Int_t conteo=0;
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
	    h_ThVSphi[i]->Fill(myDataPERP->getEVNT_track(i).Phi(),myDataPERP->getEVNT_track(i).Theta());
	    h_CosTheta[i]->Fill(myDataPERP->getEVNT_track(i).CosTheta());
	    h_Theta[i]->Fill(myDataPERP->getEVNT_track(i).Theta());
	    h_phi[i]->Fill(myDataPERP->getEVNT_track(i).Phi());
	    h_DBe[i]->Fill(deltbeta[i]);
	    h_p[i]->Fill(myDataPERP->getEVNT_track(i).Rho());
	}

	
        // for (int ii=0;ii<myDataPERP->getNum_photons();ii++){
	
        // }

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
