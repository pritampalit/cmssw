//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 22 21:43:10 2022 by ROOT version 6.22/09
// from TTree eventTree/Event tree
// found on file: bstophiphiFromttbarinput_200evt.root
//////////////////////////////////////////////////////////

#ifndef bstophiphi_phase2_h
#define bstophiphi_phase2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class bstophiphi_phase2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *trk_pt;
   vector<float>   *trk_eta;
   vector<float>   *trk_phi;
   vector<float>   *trk_phi_local;
   vector<float>   *trk_d0;
   vector<float>   *trk_z0;
   vector<float>   *trk_chi2;
   vector<float>   *trk_chi2dof;
   vector<float>   *trk_chi2rphi;
   vector<float>   *trk_chi2rz;
   vector<float>   *trk_bendchi2;
   vector<float>   *trk_MVA1;
   vector<int>     *trk_nstub;
   vector<int>     *trk_lhits;
   vector<int>     *trk_dhits;
   vector<int>     *trk_seed;
   vector<int>     *trk_hitpattern;
   vector<unsigned int> *trk_phiSector;
   vector<int>     *trk_genuine;
   vector<int>     *trk_loose;
   vector<int>     *trk_unknown;
   vector<int>     *trk_combinatoric;
   vector<int>     *trk_fake;
   vector<int>     *trk_matchtp_pdgid;
   vector<float>   *trk_matchtp_pt;
   vector<float>   *trk_matchtp_eta;
   vector<float>   *trk_matchtp_phi;
   vector<float>   *trk_matchtp_z0;
   vector<float>   *trk_matchtp_dxy;
   vector<float>   *trk_gtt_pt;
   vector<float>   *trk_gtt_eta;
   vector<float>   *trk_gtt_phi;
   vector<int>     *trk_gtt_selected_index;
   vector<int>     *trk_gtt_selected_emulation_index;
   vector<int>     *trk_poskaon_gtt_selected_index;
   vector<int>     *trk_poskaon_gtt_selected_emulation_index;
   vector<int>     *trk_negkaon_gtt_selected_index;
   vector<int>     *trk_negkaon_gtt_selected_emulation_index;
   vector<float>   *trkExt_pt;
   vector<float>   *trkExt_eta;
   vector<float>   *trkExt_phi;
   vector<float>   *trkExt_phi_local;
   vector<float>   *trkExt_d0;
   vector<float>   *trkExt_z0;
   vector<float>   *trkExt_chi2;
   vector<float>   *trkExt_chi2dof;
   vector<float>   *trkExt_chi2rphi;
   vector<float>   *trkExt_chi2rz;
   vector<float>   *trkExt_bendchi2;
   vector<float>   *trkExt_MVA;
   vector<int>     *trkExt_nstub;
   vector<int>     *trkExt_lhits;
   vector<int>     *trkExt_dhits;
   vector<int>     *trkExt_seed;
   vector<int>     *trkExt_hitpattern;
   vector<unsigned int> *trkExt_phiSector;
   vector<int>     *trkExt_genuine;
   vector<int>     *trkExt_loose;
   vector<int>     *trkExt_unknown;
   vector<int>     *trkExt_combinatoric;
   vector<int>     *trkExt_fake;
   vector<int>     *trkExt_matchtp_pdgid;
   vector<float>   *trkExt_matchtp_pt;
   vector<float>   *trkExt_matchtp_eta;
   vector<float>   *trkExt_matchtp_phi;
   vector<float>   *trkExt_matchtp_z0;
   vector<float>   *trkExt_matchtp_dxy;
   vector<float>   *trkExt_gtt_pt;
   vector<float>   *trkExt_gtt_eta;
   vector<float>   *trkExt_gtt_phi;
   vector<int>     *trkExt_gtt_selected_index;
   vector<int>     *trkExt_gtt_selected_emulation_index;
   vector<float>   *tp_pt;
   vector<float>   *tp_eta;
   vector<float>   *tp_phi;
   vector<float>   *tp_dxy;
   vector<float>   *tp_d0;
   vector<float>   *tp_z0;
   vector<float>   *tp_d0_prod;
   vector<float>   *tp_z0_prod;
   vector<int>     *tp_pdgid;
   vector<int>     *tp_nmatch;
   vector<int>     *tp_nstub;
   vector<int>     *tp_eventid;
   vector<int>     *tp_charge;
   vector<float>   *matchtrk_pt;
   vector<float>   *matchtrk_eta;
   vector<float>   *matchtrk_phi;
   vector<float>   *matchtrk_z0;
   vector<float>   *matchtrk_d0;
   vector<float>   *matchtrk_chi2;
   vector<float>   *matchtrk_chi2dof;
   vector<float>   *matchtrk_chi2rphi;
   vector<float>   *matchtrk_chi2rz;
   vector<float>   *matchtrk_bendchi2;
   vector<float>   *matchtrk_MVA1;
   vector<int>     *matchtrk_nstub;
   vector<int>     *matchtrk_lhits;
   vector<int>     *matchtrk_dhits;
   vector<int>     *matchtrk_seed;
   vector<int>     *matchtrk_hitpattern;
   vector<float>   *matchtrkExt_pt;
   vector<float>   *matchtrkExt_eta;
   vector<float>   *matchtrkExt_phi;
   vector<float>   *matchtrkExt_z0;
   vector<float>   *matchtrkExt_d0;
   vector<float>   *matchtrkExt_chi2;
   vector<float>   *matchtrkExt_chi2dof;
   vector<float>   *matchtrkExt_chi2rphi;
   vector<float>   *matchtrkExt_chi2rz;
   vector<float>   *matchtrkExt_bendchi2;
   vector<float>   *matchtrkExt_MVA;
   vector<int>     *matchtrkExt_nstub;
   vector<int>     *matchtrkExt_lhits;
   vector<int>     *matchtrkExt_dhits;
   vector<int>     *matchtrkExt_seed;
   vector<int>     *matchtrkExt_hitpattern;
   vector<float>   *pv_L1reco;
   vector<float>   *pv_L1reco_sum;
   vector<float>   *pv_L1reco_emu;
   vector<int>     *MC_lep;
   vector<float>   *pv_MC;
   vector<float>   *gen_pt;
   vector<float>   *gen_phi;
   vector<float>   *gen_pdgid;
   vector<float>   *gen_z0;
   vector<float>   *trkfastjet_eta;
   vector<float>   *trkfastjet_vz;
   vector<float>   *trkfastjet_p;
   vector<float>   *trkfastjet_pt;
   vector<float>   *trkfastjet_phi;
   vector<int>     *trkfastjet_ntracks;
   vector<float>   *trkfastjet_truetp_sumpt;
   vector<float>   *trkjet_eta;
   vector<float>   *trkjet_vz;
   vector<float>   *trkjet_p;
   vector<float>   *trkjet_pt;
   vector<float>   *trkjet_phi;
   vector<int>     *trkjet_ntracks;
   vector<int>     *trkjet_nDisplaced;
   vector<int>     *trkjet_nTight;
   vector<int>     *trkjet_nTightDisplaced;
   vector<float>   *trkjetem_eta;
   vector<float>   *trkjetem_pt;
   vector<float>   *trkjetem_phi;
   vector<float>   *trkjetem_z;
   vector<int>     *trkjetem_ntracks;
   vector<int>     *trkjetem_nxtracks;
   vector<float>   *trkfastjetExt_eta;
   vector<float>   *trkfastjetExt_vz;
   vector<float>   *trkfastjetExt_p;
   vector<float>   *trkfastjetExt_pt;
   vector<float>   *trkfastjetExt_phi;
   vector<int>     *trkfastjetExt_ntracks;
   vector<float>   *trkfastjetExt_truetp_sumpt;
   vector<float>   *trkjetExt_eta;
   vector<float>   *trkjetExt_vz;
   vector<float>   *trkjetExt_p;
   vector<float>   *trkjetExt_pt;
   vector<float>   *trkjetExt_phi;
   vector<int>     *trkjetExt_ntracks;
   vector<int>     *trkjetExt_nDisplaced;
   vector<int>     *trkjetExt_nTight;
   vector<int>     *trkjetExt_nTightDisplaced;
   vector<float>   *trkjetemExt_eta;
   vector<float>   *trkjetemExt_pt;
   vector<float>   *trkjetemExt_phi;
   vector<float>   *trkjetemExt_z;
   vector<int>     *trkjetemExt_ntracks;
   vector<int>     *trkjetemExt_nxtracks;
   Float_t         trueMET;
   Float_t         trueTkMET;
   Float_t         trkMET;
   Float_t         trkMETEmu;
   Float_t         trkMHT;
   Float_t         trkHT;
   Float_t         trkMHTEmu;
   Float_t         trkMHTEmuPhi;
   Float_t         trkHTEmu;
   Float_t         trkMETExt;
   Float_t         trkMHTExt;
   Float_t         trkHTExt;
   Float_t         trkMHTEmuExt;
   Float_t         trkMHTEmuPhiExt;
   Float_t         trkHTEmuExt;
   vector<float>   *trkphicands_eta;
   vector<float>   *trkphicands_mass;
   vector<float>   *trkphicands_pt;
   vector<float>   *trkphicands_phi;
   vector<float>   *trkphicandsemulation_eta;
   vector<float>   *trkphicandsemulation_mass;
   vector<float>   *trkphicandsemulation_pt;
   vector<float>   *trkphicandsemulation_phi;
   vector<float>   *trkphicandsExt_eta;
   vector<float>   *trkphicandsExt_mass;
   vector<float>   *trkphicandsExt_pt;
   vector<float>   *trkphicandsExt_phi;
   vector<float>   *trkphicandsemulationExt_eta;
   vector<float>   *trkphicandsemulationExt_mass;
   vector<float>   *trkphicandsemulationExt_pt;
   vector<float>   *trkphicandsemulationExt_phi;
   vector<float>   *trkbscands_eta;
   vector<float>   *trkbscands_mass;
   vector<float>   *trkbscands_pt;
   vector<float>   *trkbscands_phi;
   vector<float>   *trkbscandsemulation_eta;
   vector<float>   *trkbscandsemulation_mass;
   vector<float>   *trkbscandsemulation_pt;
   vector<float>   *trkbscandsemulation_phi;
   vector<float>   *trkbscandsExt_eta;
   vector<float>   *trkbscandsExt_mass;
   vector<float>   *trkbscandsExt_pt;
   vector<float>   *trkbscandsExt_phi;
   vector<float>   *trkbscandsemulationExt_eta;
   vector<float>   *trkbscandsemulationExt_mass;
   vector<float>   *trkbscandsemulationExt_pt;
   vector<float>   *trkbscandsemulationExt_phi;

   // List of branches
   TBranch        *b_trk_pt;   //!
   TBranch        *b_trk_eta;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trk_phi_local;   //!
   TBranch        *b_trk_d0;   //!
   TBranch        *b_trk_z0;   //!
   TBranch        *b_trk_chi2;   //!
   TBranch        *b_trk_chi2dof;   //!
   TBranch        *b_trk_chi2rphi;   //!
   TBranch        *b_trk_chi2rz;   //!
   TBranch        *b_trk_bendchi2;   //!
   TBranch        *b_trk_MVA1;   //!
   TBranch        *b_trk_nstub;   //!
   TBranch        *b_trk_lhits;   //!
   TBranch        *b_trk_dhits;   //!
   TBranch        *b_trk_seed;   //!
   TBranch        *b_trk_hitpattern;   //!
   TBranch        *b_trk_phiSector;   //!
   TBranch        *b_trk_genuine;   //!
   TBranch        *b_trk_loose;   //!
   TBranch        *b_trk_unknown;   //!
   TBranch        *b_trk_combinatoric;   //!
   TBranch        *b_trk_fake;   //!
   TBranch        *b_trk_matchtp_pdgid;   //!
   TBranch        *b_trk_matchtp_pt;   //!
   TBranch        *b_trk_matchtp_eta;   //!
   TBranch        *b_trk_matchtp_phi;   //!
   TBranch        *b_trk_matchtp_z0;   //!
   TBranch        *b_trk_matchtp_dxy;   //!
   TBranch        *b_trk_gtt_pt;   //!
   TBranch        *b_trk_gtt_eta;   //!
   TBranch        *b_trk_gtt_phi;   //!
   TBranch        *b_trk_gtt_selected_index;   //!
   TBranch        *b_trk_gtt_selected_emulation_index;   //!
   TBranch        *b_trk_poskaon_gtt_selected_index;   //!
   TBranch        *b_trk_poskaon_gtt_selected_emulation_index;   //!
   TBranch        *b_trk_negkaon_gtt_selected_index;   //!
   TBranch        *b_trk_negkaon_gtt_selected_emulation_index;   //!
   TBranch        *b_trkExt_pt;   //!
   TBranch        *b_trkExt_eta;   //!
   TBranch        *b_trkExt_phi;   //!
   TBranch        *b_trkExt_phi_local;   //!
   TBranch        *b_trkExt_d0;   //!
   TBranch        *b_trkExt_z0;   //!
   TBranch        *b_trkExt_chi2;   //!
   TBranch        *b_trkExt_chi2dof;   //!
   TBranch        *b_trkExt_chi2rphi;   //!
   TBranch        *b_trkExt_chi2rz;   //!
   TBranch        *b_trkExt_bendchi2;   //!
   TBranch        *b_trkExt_MVA;   //!
   TBranch        *b_trkExt_nstub;   //!
   TBranch        *b_trkExt_lhits;   //!
   TBranch        *b_trkExt_dhits;   //!
   TBranch        *b_trkExt_seed;   //!
   TBranch        *b_trkExt_hitpattern;   //!
   TBranch        *b_trkExt_phiSector;   //!
   TBranch        *b_trkExt_genuine;   //!
   TBranch        *b_trkExt_loose;   //!
   TBranch        *b_trkExt_unknown;   //!
   TBranch        *b_trkExt_combinatoric;   //!
   TBranch        *b_trkExt_fake;   //!
   TBranch        *b_trkExt_matchtp_pdgid;   //!
   TBranch        *b_trkExt_matchtp_pt;   //!
   TBranch        *b_trkExt_matchtp_eta;   //!
   TBranch        *b_trkExt_matchtp_phi;   //!
   TBranch        *b_trkExt_matchtp_z0;   //!
   TBranch        *b_trkExt_matchtp_dxy;   //!
   TBranch        *b_trkExt_gtt_pt;   //!
   TBranch        *b_trkExt_gtt_eta;   //!
   TBranch        *b_trkExt_gtt_phi;   //!
   TBranch        *b_trkExt_gtt_selected_index;   //!
   TBranch        *b_trkExt_gtt_selected_emulation_index;   //!
   TBranch        *b_tp_pt;   //!
   TBranch        *b_tp_eta;   //!
   TBranch        *b_tp_phi;   //!
   TBranch        *b_tp_dxy;   //!
   TBranch        *b_tp_d0;   //!
   TBranch        *b_tp_z0;   //!
   TBranch        *b_tp_d0_prod;   //!
   TBranch        *b_tp_z0_prod;   //!
   TBranch        *b_tp_pdgid;   //!
   TBranch        *b_tp_nmatch;   //!
   TBranch        *b_tp_nstub;   //!
   TBranch        *b_tp_eventid;   //!
   TBranch        *b_tp_charge;   //!
   TBranch        *b_matchtrk_pt;   //!
   TBranch        *b_matchtrk_eta;   //!
   TBranch        *b_matchtrk_phi;   //!
   TBranch        *b_matchtrk_z0;   //!
   TBranch        *b_matchtrk_d0;   //!
   TBranch        *b_matchtrk_chi2;   //!
   TBranch        *b_matchtrk_chi2dof;   //!
   TBranch        *b_matchtrk_chi2rphi;   //!
   TBranch        *b_matchtrk_chi2rz;   //!
   TBranch        *b_matchtrk_bendchi2;   //!
   TBranch        *b_matchtrk_MVA1;   //!
   TBranch        *b_matchtrk_nstub;   //!
   TBranch        *b_matchtrk_lhits;   //!
   TBranch        *b_matchtrk_dhits;   //!
   TBranch        *b_matchtrk_seed;   //!
   TBranch        *b_matchtrk_hitpattern;   //!
   TBranch        *b_matchtrkExt_pt;   //!
   TBranch        *b_matchtrkExt_eta;   //!
   TBranch        *b_matchtrkExt_phi;   //!
   TBranch        *b_matchtrkExt_z0;   //!
   TBranch        *b_matchtrkExt_d0;   //!
   TBranch        *b_matchtrkExt_chi2;   //!
   TBranch        *b_matchtrkExt_chi2dof;   //!
   TBranch        *b_matchtrkExt_chi2rphi;   //!
   TBranch        *b_matchtrkExt_chi2rz;   //!
   TBranch        *b_matchtrkExt_bendchi2;   //!
   TBranch        *b_matchtrkExt_MVA;   //!
   TBranch        *b_matchtrkExt_nstub;   //!
   TBranch        *b_matchtrkExt_lhits;   //!
   TBranch        *b_matchtrkExt_dhits;   //!
   TBranch        *b_matchtrkExt_seed;   //!
   TBranch        *b_matchtrkExt_hitpattern;   //!
   TBranch        *b_pv_L1reco;   //!
   TBranch        *b_pv_L1reco_sum;   //!
   TBranch        *b_pv_L1reco_emu;   //!
   TBranch        *b_MC_lep;   //!
   TBranch        *b_pv_MC;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_pdgid;   //!
   TBranch        *b_gen_z0;   //!
   TBranch        *b_trkfastjet_eta;   //!
   TBranch        *b_trkfastjet_vz;   //!
   TBranch        *b_trkfastjet_p;   //!
   TBranch        *b_trkfastjet_pt;   //!
   TBranch        *b_trkfastjet_phi;   //!
   TBranch        *b_trkfastjet_ntracks;   //!
   TBranch        *b_trkfastjet_truetp_sumpt;   //!
   TBranch        *b_trkjet_eta;   //!
   TBranch        *b_trkjet_vz;   //!
   TBranch        *b_trkjet_p;   //!
   TBranch        *b_trkjet_pt;   //!
   TBranch        *b_trkjet_phi;   //!
   TBranch        *b_trkjet_ntracks;   //!
   TBranch        *b_trkjet_nDisplaced;   //!
   TBranch        *b_trkjet_nTight;   //!
   TBranch        *b_trkjet_nTightDisplaced;   //!
   TBranch        *b_trkjetem_eta;   //!
   TBranch        *b_trkjetem_pt;   //!
   TBranch        *b_trkjetem_phi;   //!
   TBranch        *b_trkjetem_z;   //!
   TBranch        *b_trkjetem_ntracks;   //!
   TBranch        *b_trkjetem_nxtracks;   //!
   TBranch        *b_trkfastjetExt_eta;   //!
   TBranch        *b_trkfastjetExt_vz;   //!
   TBranch        *b_trkfastjetExt_p;   //!
   TBranch        *b_trkfastjetExt_pt;   //!
   TBranch        *b_trkfastjetExt_phi;   //!
   TBranch        *b_trkfastjetExt_ntracks;   //!
   TBranch        *b_trkfastjetExt_truetp_sumpt;   //!
   TBranch        *b_trkjetExt_eta;   //!
   TBranch        *b_trkjetExt_vz;   //!
   TBranch        *b_trkjetExt_p;   //!
   TBranch        *b_trkjetExt_pt;   //!
   TBranch        *b_trkjetExt_phi;   //!
   TBranch        *b_trkjetExt_ntracks;   //!
   TBranch        *b_trkjetExt_nDisplaced;   //!
   TBranch        *b_trkjetExt_nTight;   //!
   TBranch        *b_trkjetExt_nTightDisplaced;   //!
   TBranch        *b_trkjetemExt_eta;   //!
   TBranch        *b_trkjetemExt_pt;   //!
   TBranch        *b_trkjetemExt_phi;   //!
   TBranch        *b_trkjetemExt_z;   //!
   TBranch        *b_trkjetemExt_ntracks;   //!
   TBranch        *b_trkjetemExt_nxtracks;   //!
   TBranch        *b_trueMET;   //!
   TBranch        *b_trueTkMET;   //!
   TBranch        *b_trkMET;   //!
   TBranch        *b_trkMETEmu;   //!
   TBranch        *b_trkMHT;   //!
   TBranch        *b_trkHT;   //!
   TBranch        *b_trkMHTEmu;   //!
   TBranch        *b_trkMHTEmuPhi;   //!
   TBranch        *b_trkHTEmu;   //!
   TBranch        *b_trkMETExt;   //!
   TBranch        *b_trkMHTExt;   //!
   TBranch        *b_trkHTExt;   //!
   TBranch        *b_trkMHTEmuExt;   //!
   TBranch        *b_trkMHTEmuPhiExt;   //!
   TBranch        *b_trkHTEmuExt;   //!
   TBranch        *b_trkphicands_eta;   //!
   TBranch        *b_trkphicands_mass;   //!
   TBranch        *b_trkphicands_pt;   //!
   TBranch        *b_trkphicands_phi;   //!
   TBranch        *b_trkphicandsemulation_eta;   //!
   TBranch        *b_trkphicandsemulation_mass;   //!
   TBranch        *b_trkphicandsemulation_pt;   //!
   TBranch        *b_trkphicandsemulation_phi;   //!
   TBranch        *b_trkphicandsExt_eta;   //!
   TBranch        *b_trkphicandsExt_mass;   //!
   TBranch        *b_trkphicandsExt_pt;   //!
   TBranch        *b_trkphicandsExt_phi;   //!
   TBranch        *b_trkphicandsemulationExt_eta;   //!
   TBranch        *b_trkphicandsemulationExt_mass;   //!
   TBranch        *b_trkphicandsemulationExt_pt;   //!
   TBranch        *b_trkphicandsemulationExt_phi;   //!
   TBranch        *b_trkbscands_eta;   //!
   TBranch        *b_trkbscands_mass;   //!
   TBranch        *b_trkbscands_pt;   //!
   TBranch        *b_trkbscands_phi;   //!
   TBranch        *b_trkbscandsemulation_eta;   //!
   TBranch        *b_trkbscandsemulation_mass;   //!
   TBranch        *b_trkbscandsemulation_pt;   //!
   TBranch        *b_trkbscandsemulation_phi;   //!
   TBranch        *b_trkbscandsExt_eta;   //!
   TBranch        *b_trkbscandsExt_mass;   //!
   TBranch        *b_trkbscandsExt_pt;   //!
   TBranch        *b_trkbscandsExt_phi;   //!
   TBranch        *b_trkbscandsemulationExt_eta;   //!
   TBranch        *b_trkbscandsemulationExt_mass;   //!
   TBranch        *b_trkbscandsemulationExt_pt;   //!
   TBranch        *b_trkbscandsemulationExt_phi;   //!

   bstophiphi_phase2(TTree *tree=0);
   virtual ~bstophiphi_phase2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef bstophiphi_phase2_cxx
bstophiphi_phase2::bstophiphi_phase2(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("bstophiphiFromttbarinput_200evt.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("bstophiphiFromttbarinput_200evt.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("bstophiphiFromttbarinput_200evt.root:/L1TrackNtuple");
      dir->GetObject("eventTree",tree);

   }
   Init(tree);
}

bstophiphi_phase2::~bstophiphi_phase2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t bstophiphi_phase2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t bstophiphi_phase2::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void bstophiphi_phase2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trk_pt = 0;
   trk_eta = 0;
   trk_phi = 0;
   trk_phi_local = 0;
   trk_d0 = 0;
   trk_z0 = 0;
   trk_chi2 = 0;
   trk_chi2dof = 0;
   trk_chi2rphi = 0;
   trk_chi2rz = 0;
   trk_bendchi2 = 0;
   trk_MVA1 = 0;
   trk_nstub = 0;
   trk_lhits = 0;
   trk_dhits = 0;
   trk_seed = 0;
   trk_hitpattern = 0;
   trk_phiSector = 0;
   trk_genuine = 0;
   trk_loose = 0;
   trk_unknown = 0;
   trk_combinatoric = 0;
   trk_fake = 0;
   trk_matchtp_pdgid = 0;
   trk_matchtp_pt = 0;
   trk_matchtp_eta = 0;
   trk_matchtp_phi = 0;
   trk_matchtp_z0 = 0;
   trk_matchtp_dxy = 0;
   trk_gtt_pt = 0;
   trk_gtt_eta = 0;
   trk_gtt_phi = 0;
   trk_gtt_selected_index = 0;
   trk_gtt_selected_emulation_index = 0;
   trk_poskaon_gtt_selected_index = 0;
   trk_poskaon_gtt_selected_emulation_index = 0;
   trk_negkaon_gtt_selected_index = 0;
   trk_negkaon_gtt_selected_emulation_index = 0;
   trkExt_pt = 0;
   trkExt_eta = 0;
   trkExt_phi = 0;
   trkExt_phi_local = 0;
   trkExt_d0 = 0;
   trkExt_z0 = 0;
   trkExt_chi2 = 0;
   trkExt_chi2dof = 0;
   trkExt_chi2rphi = 0;
   trkExt_chi2rz = 0;
   trkExt_bendchi2 = 0;
   trkExt_MVA = 0;
   trkExt_nstub = 0;
   trkExt_lhits = 0;
   trkExt_dhits = 0;
   trkExt_seed = 0;
   trkExt_hitpattern = 0;
   trkExt_phiSector = 0;
   trkExt_genuine = 0;
   trkExt_loose = 0;
   trkExt_unknown = 0;
   trkExt_combinatoric = 0;
   trkExt_fake = 0;
   trkExt_matchtp_pdgid = 0;
   trkExt_matchtp_pt = 0;
   trkExt_matchtp_eta = 0;
   trkExt_matchtp_phi = 0;
   trkExt_matchtp_z0 = 0;
   trkExt_matchtp_dxy = 0;
   trkExt_gtt_pt = 0;
   trkExt_gtt_eta = 0;
   trkExt_gtt_phi = 0;
   trkExt_gtt_selected_index = 0;
   trkExt_gtt_selected_emulation_index = 0;
   tp_pt = 0;
   tp_eta = 0;
   tp_phi = 0;
   tp_dxy = 0;
   tp_d0 = 0;
   tp_z0 = 0;
   tp_d0_prod = 0;
   tp_z0_prod = 0;
   tp_pdgid = 0;
   tp_nmatch = 0;
   tp_nstub = 0;
   tp_eventid = 0;
   tp_charge = 0;
   matchtrk_pt = 0;
   matchtrk_eta = 0;
   matchtrk_phi = 0;
   matchtrk_z0 = 0;
   matchtrk_d0 = 0;
   matchtrk_chi2 = 0;
   matchtrk_chi2dof = 0;
   matchtrk_chi2rphi = 0;
   matchtrk_chi2rz = 0;
   matchtrk_bendchi2 = 0;
   matchtrk_MVA1 = 0;
   matchtrk_nstub = 0;
   matchtrk_lhits = 0;
   matchtrk_dhits = 0;
   matchtrk_seed = 0;
   matchtrk_hitpattern = 0;
   matchtrkExt_pt = 0;
   matchtrkExt_eta = 0;
   matchtrkExt_phi = 0;
   matchtrkExt_z0 = 0;
   matchtrkExt_d0 = 0;
   matchtrkExt_chi2 = 0;
   matchtrkExt_chi2dof = 0;
   matchtrkExt_chi2rphi = 0;
   matchtrkExt_chi2rz = 0;
   matchtrkExt_bendchi2 = 0;
   matchtrkExt_MVA = 0;
   matchtrkExt_nstub = 0;
   matchtrkExt_lhits = 0;
   matchtrkExt_dhits = 0;
   matchtrkExt_seed = 0;
   matchtrkExt_hitpattern = 0;
   pv_L1reco = 0;
   pv_L1reco_sum = 0;
   pv_L1reco_emu = 0;
   MC_lep = 0;
   pv_MC = 0;
   gen_pt = 0;
   gen_phi = 0;
   gen_pdgid = 0;
   gen_z0 = 0;
   trkfastjet_eta = 0;
   trkfastjet_vz = 0;
   trkfastjet_p = 0;
   trkfastjet_pt = 0;
   trkfastjet_phi = 0;
   trkfastjet_ntracks = 0;
   trkfastjet_truetp_sumpt = 0;
   trkjet_eta = 0;
   trkjet_vz = 0;
   trkjet_p = 0;
   trkjet_pt = 0;
   trkjet_phi = 0;
   trkjet_ntracks = 0;
   trkjet_nDisplaced = 0;
   trkjet_nTight = 0;
   trkjet_nTightDisplaced = 0;
   trkjetem_eta = 0;
   trkjetem_pt = 0;
   trkjetem_phi = 0;
   trkjetem_z = 0;
   trkjetem_ntracks = 0;
   trkjetem_nxtracks = 0;
   trkfastjetExt_eta = 0;
   trkfastjetExt_vz = 0;
   trkfastjetExt_p = 0;
   trkfastjetExt_pt = 0;
   trkfastjetExt_phi = 0;
   trkfastjetExt_ntracks = 0;
   trkfastjetExt_truetp_sumpt = 0;
   trkjetExt_eta = 0;
   trkjetExt_vz = 0;
   trkjetExt_p = 0;
   trkjetExt_pt = 0;
   trkjetExt_phi = 0;
   trkjetExt_ntracks = 0;
   trkjetExt_nDisplaced = 0;
   trkjetExt_nTight = 0;
   trkjetExt_nTightDisplaced = 0;
   trkjetemExt_eta = 0;
   trkjetemExt_pt = 0;
   trkjetemExt_phi = 0;
   trkjetemExt_z = 0;
   trkjetemExt_ntracks = 0;
   trkjetemExt_nxtracks = 0;
   trkphicands_eta = 0;
   trkphicands_mass = 0;
   trkphicands_pt = 0;
   trkphicands_phi = 0;
   trkphicandsemulation_eta = 0;
   trkphicandsemulation_mass = 0;
   trkphicandsemulation_pt = 0;
   trkphicandsemulation_phi = 0;
   trkphicandsExt_eta = 0;
   trkphicandsExt_mass = 0;
   trkphicandsExt_pt = 0;
   trkphicandsExt_phi = 0;
   trkphicandsemulationExt_eta = 0;
   trkphicandsemulationExt_mass = 0;
   trkphicandsemulationExt_pt = 0;
   trkphicandsemulationExt_phi = 0;
   trkbscands_eta = 0;
   trkbscands_mass = 0;
   trkbscands_pt = 0;
   trkbscands_phi = 0;
   trkbscandsemulation_eta = 0;
   trkbscandsemulation_mass = 0;
   trkbscandsemulation_pt = 0;
   trkbscandsemulation_phi = 0;
   trkbscandsExt_eta = 0;
   trkbscandsExt_mass = 0;
   trkbscandsExt_pt = 0;
   trkbscandsExt_phi = 0;
   trkbscandsemulationExt_eta = 0;
   trkbscandsemulationExt_mass = 0;
   trkbscandsemulationExt_pt = 0;
   trkbscandsemulationExt_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
   fChain->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trk_phi_local", &trk_phi_local, &b_trk_phi_local);
   fChain->SetBranchAddress("trk_d0", &trk_d0, &b_trk_d0);
   fChain->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
   fChain->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
   fChain->SetBranchAddress("trk_chi2dof", &trk_chi2dof, &b_trk_chi2dof);
   fChain->SetBranchAddress("trk_chi2rphi", &trk_chi2rphi, &b_trk_chi2rphi);
   fChain->SetBranchAddress("trk_chi2rz", &trk_chi2rz, &b_trk_chi2rz);
   fChain->SetBranchAddress("trk_bendchi2", &trk_bendchi2, &b_trk_bendchi2);
   fChain->SetBranchAddress("trk_MVA1", &trk_MVA1, &b_trk_MVA1);
   fChain->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
   fChain->SetBranchAddress("trk_lhits", &trk_lhits, &b_trk_lhits);
   fChain->SetBranchAddress("trk_dhits", &trk_dhits, &b_trk_dhits);
   fChain->SetBranchAddress("trk_seed", &trk_seed, &b_trk_seed);
   fChain->SetBranchAddress("trk_hitpattern", &trk_hitpattern, &b_trk_hitpattern);
   fChain->SetBranchAddress("trk_phiSector", &trk_phiSector, &b_trk_phiSector);
   fChain->SetBranchAddress("trk_genuine", &trk_genuine, &b_trk_genuine);
   fChain->SetBranchAddress("trk_loose", &trk_loose, &b_trk_loose);
   fChain->SetBranchAddress("trk_unknown", &trk_unknown, &b_trk_unknown);
   fChain->SetBranchAddress("trk_combinatoric", &trk_combinatoric, &b_trk_combinatoric);
   fChain->SetBranchAddress("trk_fake", &trk_fake, &b_trk_fake);
   fChain->SetBranchAddress("trk_matchtp_pdgid", &trk_matchtp_pdgid, &b_trk_matchtp_pdgid);
   fChain->SetBranchAddress("trk_matchtp_pt", &trk_matchtp_pt, &b_trk_matchtp_pt);
   fChain->SetBranchAddress("trk_matchtp_eta", &trk_matchtp_eta, &b_trk_matchtp_eta);
   fChain->SetBranchAddress("trk_matchtp_phi", &trk_matchtp_phi, &b_trk_matchtp_phi);
   fChain->SetBranchAddress("trk_matchtp_z0", &trk_matchtp_z0, &b_trk_matchtp_z0);
   fChain->SetBranchAddress("trk_matchtp_dxy", &trk_matchtp_dxy, &b_trk_matchtp_dxy);
   fChain->SetBranchAddress("trk_gtt_pt", &trk_gtt_pt, &b_trk_gtt_pt);
   fChain->SetBranchAddress("trk_gtt_eta", &trk_gtt_eta, &b_trk_gtt_eta);
   fChain->SetBranchAddress("trk_gtt_phi", &trk_gtt_phi, &b_trk_gtt_phi);
   fChain->SetBranchAddress("trk_gtt_selected_index", &trk_gtt_selected_index, &b_trk_gtt_selected_index);
   fChain->SetBranchAddress("trk_gtt_selected_emulation_index", &trk_gtt_selected_emulation_index, &b_trk_gtt_selected_emulation_index);
   fChain->SetBranchAddress("trk_poskaon_gtt_selected_index", &trk_poskaon_gtt_selected_index, &b_trk_poskaon_gtt_selected_index);
   fChain->SetBranchAddress("trk_poskaon_gtt_selected_emulation_index", &trk_poskaon_gtt_selected_emulation_index, &b_trk_poskaon_gtt_selected_emulation_index);
   fChain->SetBranchAddress("trk_negkaon_gtt_selected_index", &trk_negkaon_gtt_selected_index, &b_trk_negkaon_gtt_selected_index);
   fChain->SetBranchAddress("trk_negkaon_gtt_selected_emulation_index", &trk_negkaon_gtt_selected_emulation_index, &b_trk_negkaon_gtt_selected_emulation_index);
   fChain->SetBranchAddress("trkExt_pt", &trkExt_pt, &b_trkExt_pt);
   fChain->SetBranchAddress("trkExt_eta", &trkExt_eta, &b_trkExt_eta);
   fChain->SetBranchAddress("trkExt_phi", &trkExt_phi, &b_trkExt_phi);
   fChain->SetBranchAddress("trkExt_phi_local", &trkExt_phi_local, &b_trkExt_phi_local);
   fChain->SetBranchAddress("trkExt_d0", &trkExt_d0, &b_trkExt_d0);
   fChain->SetBranchAddress("trkExt_z0", &trkExt_z0, &b_trkExt_z0);
   fChain->SetBranchAddress("trkExt_chi2", &trkExt_chi2, &b_trkExt_chi2);
   fChain->SetBranchAddress("trkExt_chi2dof", &trkExt_chi2dof, &b_trkExt_chi2dof);
   fChain->SetBranchAddress("trkExt_chi2rphi", &trkExt_chi2rphi, &b_trkExt_chi2rphi);
   fChain->SetBranchAddress("trkExt_chi2rz", &trkExt_chi2rz, &b_trkExt_chi2rz);
   fChain->SetBranchAddress("trkExt_bendchi2", &trkExt_bendchi2, &b_trkExt_bendchi2);
   fChain->SetBranchAddress("trkExt_MVA", &trkExt_MVA, &b_trkExt_MVA);
   fChain->SetBranchAddress("trkExt_nstub", &trkExt_nstub, &b_trkExt_nstub);
   fChain->SetBranchAddress("trkExt_lhits", &trkExt_lhits, &b_trkExt_lhits);
   fChain->SetBranchAddress("trkExt_dhits", &trkExt_dhits, &b_trkExt_dhits);
   fChain->SetBranchAddress("trkExt_seed", &trkExt_seed, &b_trkExt_seed);
   fChain->SetBranchAddress("trkExt_hitpattern", &trkExt_hitpattern, &b_trkExt_hitpattern);
   fChain->SetBranchAddress("trkExt_phiSector", &trkExt_phiSector, &b_trkExt_phiSector);
   fChain->SetBranchAddress("trkExt_genuine", &trkExt_genuine, &b_trkExt_genuine);
   fChain->SetBranchAddress("trkExt_loose", &trkExt_loose, &b_trkExt_loose);
   fChain->SetBranchAddress("trkExt_unknown", &trkExt_unknown, &b_trkExt_unknown);
   fChain->SetBranchAddress("trkExt_combinatoric", &trkExt_combinatoric, &b_trkExt_combinatoric);
   fChain->SetBranchAddress("trkExt_fake", &trkExt_fake, &b_trkExt_fake);
   fChain->SetBranchAddress("trkExt_matchtp_pdgid", &trkExt_matchtp_pdgid, &b_trkExt_matchtp_pdgid);
   fChain->SetBranchAddress("trkExt_matchtp_pt", &trkExt_matchtp_pt, &b_trkExt_matchtp_pt);
   fChain->SetBranchAddress("trkExt_matchtp_eta", &trkExt_matchtp_eta, &b_trkExt_matchtp_eta);
   fChain->SetBranchAddress("trkExt_matchtp_phi", &trkExt_matchtp_phi, &b_trkExt_matchtp_phi);
   fChain->SetBranchAddress("trkExt_matchtp_z0", &trkExt_matchtp_z0, &b_trkExt_matchtp_z0);
   fChain->SetBranchAddress("trkExt_matchtp_dxy", &trkExt_matchtp_dxy, &b_trkExt_matchtp_dxy);
   fChain->SetBranchAddress("trkExt_gtt_pt", &trkExt_gtt_pt, &b_trkExt_gtt_pt);
   fChain->SetBranchAddress("trkExt_gtt_eta", &trkExt_gtt_eta, &b_trkExt_gtt_eta);
   fChain->SetBranchAddress("trkExt_gtt_phi", &trkExt_gtt_phi, &b_trkExt_gtt_phi);
   fChain->SetBranchAddress("trkExt_gtt_selected_index", &trkExt_gtt_selected_index, &b_trkExt_gtt_selected_index);
   fChain->SetBranchAddress("trkExt_gtt_selected_emulation_index", &trkExt_gtt_selected_emulation_index, &b_trkExt_gtt_selected_emulation_index);
   fChain->SetBranchAddress("tp_pt", &tp_pt, &b_tp_pt);
   fChain->SetBranchAddress("tp_eta", &tp_eta, &b_tp_eta);
   fChain->SetBranchAddress("tp_phi", &tp_phi, &b_tp_phi);
   fChain->SetBranchAddress("tp_dxy", &tp_dxy, &b_tp_dxy);
   fChain->SetBranchAddress("tp_d0", &tp_d0, &b_tp_d0);
   fChain->SetBranchAddress("tp_z0", &tp_z0, &b_tp_z0);
   fChain->SetBranchAddress("tp_d0_prod", &tp_d0_prod, &b_tp_d0_prod);
   fChain->SetBranchAddress("tp_z0_prod", &tp_z0_prod, &b_tp_z0_prod);
   fChain->SetBranchAddress("tp_pdgid", &tp_pdgid, &b_tp_pdgid);
   fChain->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
   fChain->SetBranchAddress("tp_nstub", &tp_nstub, &b_tp_nstub);
   fChain->SetBranchAddress("tp_eventid", &tp_eventid, &b_tp_eventid);
   fChain->SetBranchAddress("tp_charge", &tp_charge, &b_tp_charge);
   fChain->SetBranchAddress("matchtrk_pt", &matchtrk_pt, &b_matchtrk_pt);
   fChain->SetBranchAddress("matchtrk_eta", &matchtrk_eta, &b_matchtrk_eta);
   fChain->SetBranchAddress("matchtrk_phi", &matchtrk_phi, &b_matchtrk_phi);
   fChain->SetBranchAddress("matchtrk_z0", &matchtrk_z0, &b_matchtrk_z0);
   fChain->SetBranchAddress("matchtrk_d0", &matchtrk_d0, &b_matchtrk_d0);
   fChain->SetBranchAddress("matchtrk_chi2", &matchtrk_chi2, &b_matchtrk_chi2);
   fChain->SetBranchAddress("matchtrk_chi2dof", &matchtrk_chi2dof, &b_matchtrk_chi2dof);
   fChain->SetBranchAddress("matchtrk_chi2rphi", &matchtrk_chi2rphi, &b_matchtrk_chi2rphi);
   fChain->SetBranchAddress("matchtrk_chi2rz", &matchtrk_chi2rz, &b_matchtrk_chi2rz);
   fChain->SetBranchAddress("matchtrk_bendchi2", &matchtrk_bendchi2, &b_matchtrk_bendchi2);
   fChain->SetBranchAddress("matchtrk_MVA1", &matchtrk_MVA1, &b_matchtrk_MVA1);
   fChain->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);
   fChain->SetBranchAddress("matchtrk_lhits", &matchtrk_lhits, &b_matchtrk_lhits);
   fChain->SetBranchAddress("matchtrk_dhits", &matchtrk_dhits, &b_matchtrk_dhits);
   fChain->SetBranchAddress("matchtrk_seed", &matchtrk_seed, &b_matchtrk_seed);
   fChain->SetBranchAddress("matchtrk_hitpattern", &matchtrk_hitpattern, &b_matchtrk_hitpattern);
   fChain->SetBranchAddress("matchtrkExt_pt", &matchtrkExt_pt, &b_matchtrkExt_pt);
   fChain->SetBranchAddress("matchtrkExt_eta", &matchtrkExt_eta, &b_matchtrkExt_eta);
   fChain->SetBranchAddress("matchtrkExt_phi", &matchtrkExt_phi, &b_matchtrkExt_phi);
   fChain->SetBranchAddress("matchtrkExt_z0", &matchtrkExt_z0, &b_matchtrkExt_z0);
   fChain->SetBranchAddress("matchtrkExt_d0", &matchtrkExt_d0, &b_matchtrkExt_d0);
   fChain->SetBranchAddress("matchtrkExt_chi2", &matchtrkExt_chi2, &b_matchtrkExt_chi2);
   fChain->SetBranchAddress("matchtrkExt_chi2dof", &matchtrkExt_chi2dof, &b_matchtrkExt_chi2dof);
   fChain->SetBranchAddress("matchtrkExt_chi2rphi", &matchtrkExt_chi2rphi, &b_matchtrkExt_chi2rphi);
   fChain->SetBranchAddress("matchtrkExt_chi2rz", &matchtrkExt_chi2rz, &b_matchtrkExt_chi2rz);
   fChain->SetBranchAddress("matchtrkExt_bendchi2", &matchtrkExt_bendchi2, &b_matchtrkExt_bendchi2);
   fChain->SetBranchAddress("matchtrkExt_MVA", &matchtrkExt_MVA, &b_matchtrkExt_MVA);
   fChain->SetBranchAddress("matchtrkExt_nstub", &matchtrkExt_nstub, &b_matchtrkExt_nstub);
   fChain->SetBranchAddress("matchtrkExt_lhits", &matchtrkExt_lhits, &b_matchtrkExt_lhits);
   fChain->SetBranchAddress("matchtrkExt_dhits", &matchtrkExt_dhits, &b_matchtrkExt_dhits);
   fChain->SetBranchAddress("matchtrkExt_seed", &matchtrkExt_seed, &b_matchtrkExt_seed);
   fChain->SetBranchAddress("matchtrkExt_hitpattern", &matchtrkExt_hitpattern, &b_matchtrkExt_hitpattern);
   fChain->SetBranchAddress("pv_L1reco", &pv_L1reco, &b_pv_L1reco);
   fChain->SetBranchAddress("pv_L1reco_sum", &pv_L1reco_sum, &b_pv_L1reco_sum);
   fChain->SetBranchAddress("pv_L1reco_emu", &pv_L1reco_emu, &b_pv_L1reco_emu);
   fChain->SetBranchAddress("MC_lep", &MC_lep, &b_MC_lep);
   fChain->SetBranchAddress("pv_MC", &pv_MC, &b_pv_MC);
   fChain->SetBranchAddress("gen_pt", &gen_pt, &b_gen_pt);
   fChain->SetBranchAddress("gen_phi", &gen_phi, &b_gen_phi);
   fChain->SetBranchAddress("gen_pdgid", &gen_pdgid, &b_gen_pdgid);
   fChain->SetBranchAddress("gen_z0", &gen_z0, &b_gen_z0);
   fChain->SetBranchAddress("trkfastjet_eta", &trkfastjet_eta, &b_trkfastjet_eta);
   fChain->SetBranchAddress("trkfastjet_vz", &trkfastjet_vz, &b_trkfastjet_vz);
   fChain->SetBranchAddress("trkfastjet_p", &trkfastjet_p, &b_trkfastjet_p);
   fChain->SetBranchAddress("trkfastjet_pt", &trkfastjet_pt, &b_trkfastjet_pt);
   fChain->SetBranchAddress("trkfastjet_phi", &trkfastjet_phi, &b_trkfastjet_phi);
   fChain->SetBranchAddress("trkfastjet_ntracks", &trkfastjet_ntracks, &b_trkfastjet_ntracks);
   fChain->SetBranchAddress("trkfastjet_truetp_sumpt", &trkfastjet_truetp_sumpt, &b_trkfastjet_truetp_sumpt);
   fChain->SetBranchAddress("trkjet_eta", &trkjet_eta, &b_trkjet_eta);
   fChain->SetBranchAddress("trkjet_vz", &trkjet_vz, &b_trkjet_vz);
   fChain->SetBranchAddress("trkjet_p", &trkjet_p, &b_trkjet_p);
   fChain->SetBranchAddress("trkjet_pt", &trkjet_pt, &b_trkjet_pt);
   fChain->SetBranchAddress("trkjet_phi", &trkjet_phi, &b_trkjet_phi);
   fChain->SetBranchAddress("trkjet_ntracks", &trkjet_ntracks, &b_trkjet_ntracks);
   fChain->SetBranchAddress("trkjet_nDisplaced", &trkjet_nDisplaced, &b_trkjet_nDisplaced);
   fChain->SetBranchAddress("trkjet_nTight", &trkjet_nTight, &b_trkjet_nTight);
   fChain->SetBranchAddress("trkjet_nTightDisplaced", &trkjet_nTightDisplaced, &b_trkjet_nTightDisplaced);
   fChain->SetBranchAddress("trkjetem_eta", &trkjetem_eta, &b_trkjetem_eta);
   fChain->SetBranchAddress("trkjetem_pt", &trkjetem_pt, &b_trkjetem_pt);
   fChain->SetBranchAddress("trkjetem_phi", &trkjetem_phi, &b_trkjetem_phi);
   fChain->SetBranchAddress("trkjetem_z", &trkjetem_z, &b_trkjetem_z);
   fChain->SetBranchAddress("trkjetem_ntracks", &trkjetem_ntracks, &b_trkjetem_ntracks);
   fChain->SetBranchAddress("trkjetem_nxtracks", &trkjetem_nxtracks, &b_trkjetem_nxtracks);
   fChain->SetBranchAddress("trkfastjetExt_eta", &trkfastjetExt_eta, &b_trkfastjetExt_eta);
   fChain->SetBranchAddress("trkfastjetExt_vz", &trkfastjetExt_vz, &b_trkfastjetExt_vz);
   fChain->SetBranchAddress("trkfastjetExt_p", &trkfastjetExt_p, &b_trkfastjetExt_p);
   fChain->SetBranchAddress("trkfastjetExt_pt", &trkfastjetExt_pt, &b_trkfastjetExt_pt);
   fChain->SetBranchAddress("trkfastjetExt_phi", &trkfastjetExt_phi, &b_trkfastjetExt_phi);
   fChain->SetBranchAddress("trkfastjetExt_ntracks", &trkfastjetExt_ntracks, &b_trkfastjetExt_ntracks);
   fChain->SetBranchAddress("trkfastjetExt_truetp_sumpt", &trkfastjetExt_truetp_sumpt, &b_trkfastjetExt_truetp_sumpt);
   fChain->SetBranchAddress("trkjetExt_eta", &trkjetExt_eta, &b_trkjetExt_eta);
   fChain->SetBranchAddress("trkjetExt_vz", &trkjetExt_vz, &b_trkjetExt_vz);
   fChain->SetBranchAddress("trkjetExt_p", &trkjetExt_p, &b_trkjetExt_p);
   fChain->SetBranchAddress("trkjetExt_pt", &trkjetExt_pt, &b_trkjetExt_pt);
   fChain->SetBranchAddress("trkjetExt_phi", &trkjetExt_phi, &b_trkjetExt_phi);
   fChain->SetBranchAddress("trkjetExt_ntracks", &trkjetExt_ntracks, &b_trkjetExt_ntracks);
   fChain->SetBranchAddress("trkjetExt_nDisplaced", &trkjetExt_nDisplaced, &b_trkjetExt_nDisplaced);
   fChain->SetBranchAddress("trkjetExt_nTight", &trkjetExt_nTight, &b_trkjetExt_nTight);
   fChain->SetBranchAddress("trkjetExt_nTightDisplaced", &trkjetExt_nTightDisplaced, &b_trkjetExt_nTightDisplaced);
   fChain->SetBranchAddress("trkjetemExt_eta", &trkjetemExt_eta, &b_trkjetemExt_eta);
   fChain->SetBranchAddress("trkjetemExt_pt", &trkjetemExt_pt, &b_trkjetemExt_pt);
   fChain->SetBranchAddress("trkjetemExt_phi", &trkjetemExt_phi, &b_trkjetemExt_phi);
   fChain->SetBranchAddress("trkjetemExt_z", &trkjetemExt_z, &b_trkjetemExt_z);
   fChain->SetBranchAddress("trkjetemExt_ntracks", &trkjetemExt_ntracks, &b_trkjetemExt_ntracks);
   fChain->SetBranchAddress("trkjetemExt_nxtracks", &trkjetemExt_nxtracks, &b_trkjetemExt_nxtracks);
   fChain->SetBranchAddress("trueMET", &trueMET, &b_trueMET);
   fChain->SetBranchAddress("trueTkMET", &trueTkMET, &b_trueTkMET);
   fChain->SetBranchAddress("trkMET", &trkMET, &b_trkMET);
   fChain->SetBranchAddress("trkMETEmu", &trkMETEmu, &b_trkMETEmu);
   fChain->SetBranchAddress("trkMHT", &trkMHT, &b_trkMHT);
   fChain->SetBranchAddress("trkHT", &trkHT, &b_trkHT);
   fChain->SetBranchAddress("trkMHTEmu", &trkMHTEmu, &b_trkMHTEmu);
   fChain->SetBranchAddress("trkMHTEmuPhi", &trkMHTEmuPhi, &b_trkMHTEmuPhi);
   fChain->SetBranchAddress("trkHTEmu", &trkHTEmu, &b_trkHTEmu);
   fChain->SetBranchAddress("trkMETExt", &trkMETExt, &b_trkMETExt);
   fChain->SetBranchAddress("trkMHTExt", &trkMHTExt, &b_trkMHTExt);
   fChain->SetBranchAddress("trkHTExt", &trkHTExt, &b_trkHTExt);
   fChain->SetBranchAddress("trkMHTEmuExt", &trkMHTEmuExt, &b_trkMHTEmuExt);
   fChain->SetBranchAddress("trkMHTEmuPhiExt", &trkMHTEmuPhiExt, &b_trkMHTEmuPhiExt);
   fChain->SetBranchAddress("trkHTEmuExt", &trkHTEmuExt, &b_trkHTEmuExt);
   fChain->SetBranchAddress("trkphicands_eta", &trkphicands_eta, &b_trkphicands_eta);
   fChain->SetBranchAddress("trkphicands_mass", &trkphicands_mass, &b_trkphicands_mass);
   fChain->SetBranchAddress("trkphicands_pt", &trkphicands_pt, &b_trkphicands_pt);
   fChain->SetBranchAddress("trkphicands_phi", &trkphicands_phi, &b_trkphicands_phi);
   fChain->SetBranchAddress("trkphicandsemulation_eta", &trkphicandsemulation_eta, &b_trkphicandsemulation_eta);
   fChain->SetBranchAddress("trkphicandsemulation_mass", &trkphicandsemulation_mass, &b_trkphicandsemulation_mass);
   fChain->SetBranchAddress("trkphicandsemulation_pt", &trkphicandsemulation_pt, &b_trkphicandsemulation_pt);
   fChain->SetBranchAddress("trkphicandsemulation_phi", &trkphicandsemulation_phi, &b_trkphicandsemulation_phi);
   fChain->SetBranchAddress("trkphicandsExt_eta", &trkphicandsExt_eta, &b_trkphicandsExt_eta);
   fChain->SetBranchAddress("trkphicandsExt_mass", &trkphicandsExt_mass, &b_trkphicandsExt_mass);
   fChain->SetBranchAddress("trkphicandsExt_pt", &trkphicandsExt_pt, &b_trkphicandsExt_pt);
   fChain->SetBranchAddress("trkphicandsExt_phi", &trkphicandsExt_phi, &b_trkphicandsExt_phi);
   fChain->SetBranchAddress("trkphicandsemulationExt_eta", &trkphicandsemulationExt_eta, &b_trkphicandsemulationExt_eta);
   fChain->SetBranchAddress("trkphicandsemulationExt_mass", &trkphicandsemulationExt_mass, &b_trkphicandsemulationExt_mass);
   fChain->SetBranchAddress("trkphicandsemulationExt_pt", &trkphicandsemulationExt_pt, &b_trkphicandsemulationExt_pt);
   fChain->SetBranchAddress("trkphicandsemulationExt_phi", &trkphicandsemulationExt_phi, &b_trkphicandsemulationExt_phi);
   fChain->SetBranchAddress("trkbscands_eta", &trkbscands_eta, &b_trkbscands_eta);
   fChain->SetBranchAddress("trkbscands_mass", &trkbscands_mass, &b_trkbscands_mass);
   fChain->SetBranchAddress("trkbscands_pt", &trkbscands_pt, &b_trkbscands_pt);
   fChain->SetBranchAddress("trkbscands_phi", &trkbscands_phi, &b_trkbscands_phi);
   fChain->SetBranchAddress("trkbscandsemulation_eta", &trkbscandsemulation_eta, &b_trkbscandsemulation_eta);
   fChain->SetBranchAddress("trkbscandsemulation_mass", &trkbscandsemulation_mass, &b_trkbscandsemulation_mass);
   fChain->SetBranchAddress("trkbscandsemulation_pt", &trkbscandsemulation_pt, &b_trkbscandsemulation_pt);
   fChain->SetBranchAddress("trkbscandsemulation_phi", &trkbscandsemulation_phi, &b_trkbscandsemulation_phi);
   fChain->SetBranchAddress("trkbscandsExt_eta", &trkbscandsExt_eta, &b_trkbscandsExt_eta);
   fChain->SetBranchAddress("trkbscandsExt_mass", &trkbscandsExt_mass, &b_trkbscandsExt_mass);
   fChain->SetBranchAddress("trkbscandsExt_pt", &trkbscandsExt_pt, &b_trkbscandsExt_pt);
   fChain->SetBranchAddress("trkbscandsExt_phi", &trkbscandsExt_phi, &b_trkbscandsExt_phi);
   fChain->SetBranchAddress("trkbscandsemulationExt_eta", &trkbscandsemulationExt_eta, &b_trkbscandsemulationExt_eta);
   fChain->SetBranchAddress("trkbscandsemulationExt_mass", &trkbscandsemulationExt_mass, &b_trkbscandsemulationExt_mass);
   fChain->SetBranchAddress("trkbscandsemulationExt_pt", &trkbscandsemulationExt_pt, &b_trkbscandsemulationExt_pt);
   fChain->SetBranchAddress("trkbscandsemulationExt_phi", &trkbscandsemulationExt_phi, &b_trkbscandsemulationExt_phi);
   Notify();
}

Bool_t bstophiphi_phase2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void bstophiphi_phase2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t bstophiphi_phase2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef bstophiphi_phase2_cxx
