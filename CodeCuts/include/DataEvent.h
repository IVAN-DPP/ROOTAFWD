#ifndef DATAEVENT_H
#define DATAEVENT_H

#include "Libraries.h"

#pragma link C++ class vector<TLorentzVector>+;
#pragma link C++ class vector<TVector3>+;

using namespace std;

class DataEvent{


private:
  vector<TLorentzVector>
    *loc_EVNT_track,
    *loc_eloss_track,
    *loc_kpi_mm,
    *loc_k_mm,
    *loc_ppi_mm,
    *loc_pipi_mm,
    *loc_d_kppi_mm,
    *loc_d_kp_mm;
  
  vector<int>
    *loc_EVNT_q,
    *loc_EVNT_scsec,
    *loc_EVNT_scpad,
    *loc_EVNT_schit,
    *loc_EVNT_stsec,
    *loc_EVNT_sthit,
    *loc_tagr_eid,
    *loc_tagr_tid,
    *loc_tagr_stat,
    //    *loc_STPB_stat,
    //    *loc_SCPB_stat,
    //    *loc_DCPB_stat,
    //    *loc_ECPB_stat,
    //    *loc_EVNT_stat,
    *loc_STPB_sthid,
    *loc_SCPB_ScPdHt,
    *loc_index_k,
    *loc_index_pi;
  
  vector<float>
    *loc_EVNT_bem,
    *loc_EVNT_sc_t,
    *loc_EVNT_sc_d,
    *loc_EVNT_st_t,
    *loc_EVNT_st_d,
    *loc_EVNT_sc_e,
    *loc_delt_t_k,
    *loc_delt_t_pi,
    *loc_tagr_epho,
    *loc_tagr_tpho;
  
  int
    loc_trip_flag,
    loc_numofpart,
    loc_num_pos,
    loc_num_neg,
    loc_num_neu,
    loc_head_eventnum,
    loc_head_runnum,
    loc_num_deuterons,
    loc_num_protons,
    loc_num_poskaons,
    loc_num_pospions,
    loc_num_negkaons,
    loc_num_negpions,
    loc_coh_plan,
    loc_coh_radi,
    loc_coh_plan_db,
    loc_numph_k,
    loc_numph_pi;
  
  float loc_coh_edge,
    loc_beam_en,
    loc_coh_edge_nom;
  
  
  vector<TVector3> *loc_EVNT_vertex;
  TVector3 *loc_MVRT_vertex;
  UInt_t trigger_pattern;
  Int_t RC26_bit[32];
  TFile *loc_file;
  TChain *loc_Tree;
  int eventno;
  int numofRootFiles;
 public:
  
  DataEvent(string fileName, string treeName);
  DataEvent(string filesfileName, string treeName, int numOfFiles);
  void SetUpTree(string fileName, string treeName);
  void SetUpChain(string filesfileName, string treeName);
  void SetUpBranches();
  
  TLorentzVector getEVNT_track(int ip){return loc_EVNT_track->at(ip);} //Track 4-vector
  TLorentzVector geteloss_track(int ip){return loc_eloss_track->at(ip);}//Eloss corrected Track 4-vector
  int getEVNT_q(int ip){return loc_EVNT_q->at(ip);} //charge
  int getEVNT_scsec(int ip){return loc_EVNT_scsec->at(ip);} //SC sector
  int getEVNT_scpad(int ip){return loc_EVNT_scpad->at(ip);} //SC paddle
  int getEVNT_schit(int ip){return loc_EVNT_schit->at(ip);} //SC hit
  int getEVNT_stsec(int ip){return loc_EVNT_stsec->at(ip);} //ST sector
  int getEVNT_sthit(int ip){return loc_EVNT_sthit->at(ip);} //ST hit
  int getTAGR_eid(int ip){return loc_tagr_eid->at(ip);} //TAGR eid
  int getTAGR_tid(int ip){return loc_tagr_tid->at(ip);} //TAGR tid
  int getTAGR_stat(int ip){return loc_tagr_stat->at(ip);} //TAGR status
  int getSTPB_sthid(int ip){return loc_STPB_sthid->at(ip);} //STBP hitd
  int getSCPB_ScPdHt(int ip){return loc_SCPB_ScPdHt->at(ip);} //SCBP ScPdHt
  float getEVNT_bem(int ip){return loc_EVNT_bem->at(ip);} //beta measured
  float getEVNT_sc_t(int ip){return loc_EVNT_sc_t->at(ip);} //sc time
  float getEVNT_sc_d(int ip){return loc_EVNT_sc_d->at(ip);} //sc d
  float getEVNT_st_t(int ip){return loc_EVNT_st_t->at(ip);} //st time
  float getEVNT_st_d(int ip){return loc_EVNT_st_d->at(ip);} //st d
  float getEVNT_sc_e(int ip){return loc_EVNT_sc_e->at(ip);} //sc energy
  float getTAGR_epho(int ip){return loc_tagr_epho->at(ip);} //TAGR epho
  float getTAGR_tpho(int ip){return loc_tagr_tpho->at(ip);} //TAGR tpho
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
  float getCoh_edge(){return loc_coh_edge;}
  float getBeam_en(){return loc_beam_en;}
  float getCoh_edge_nom(){return loc_coh_edge_nom;}
  int getCoh_plan_db(){return loc_coh_plan_db;}
  int getCoh_radi(){return loc_coh_radi;}
  int getCoh_plan(){return loc_coh_plan;}
  float getDelt_t_k(int ip){return loc_delt_t_k->at(ip);}
  float getDelt_t_pi(int ip){return loc_delt_t_pi->at(ip);}
  int getNumph_k(){return loc_numph_k;}
  int getNumph_pi(){return loc_numph_pi;}
  int getIndex_k(int ip){return loc_index_k->at(ip);}
  int getIndex_pi(int ip){return loc_index_pi->at(ip);}
  TLorentzVector getKpi_mm(int ip){return loc_kpi_mm->at(ip);} //Track 4-vector
  TLorentzVector getK_mm(int ip){return loc_k_mm->at(ip);} //Track 4-vector
  TLorentzVector getPpi_mm(int ip){return loc_ppi_mm->at(ip);} //Track 4-vector
  TLorentzVector getPipi_mm(int ip){return loc_pipi_mm->at(ip);} //Track 4-vector
  TLorentzVector getDKppi_mm(int ip){return loc_d_kppi_mm->at(ip);} //Track 4-vector
  TLorentzVector getDKp_mm(int ip){return loc_d_kp_mm->at(ip);} //Track 4-vector
  TVector3 getEVNT_vertex(int ip){return loc_EVNT_vertex->at(ip);} //EVNT vertex
  TVector3 getMVRT_vertex(){return *loc_MVRT_vertex;} //MVRT vertex
  int getNextEntry();
  int getEntry(){return eventno;}
  int getEntries();
};
#endif
