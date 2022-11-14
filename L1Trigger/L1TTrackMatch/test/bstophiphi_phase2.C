#define bstophiphi_phase2_cxx
#include "bstophiphi_phase2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cmath>
#include <math.h>
void bstophiphi_phase2::Loop()
{
//   In a ROOT session, you can do:
//      root> .L bstophiphi_phase2.C
//      root> bstophiphi_phase2 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   std::cout << "Total number of entries: " << nentries << std::endl;
   
   TH1F *hPosKaonPt  = new TH1F("hPosKaonPt","PosKaonPt", 25, 0., 10.);
   TH1F *hPosKaonEmuPt  = new TH1F("hPosKaonEmuPt","PosKaonEmuPt", 25, 0., 10.);

   TH1F *hPosKaonPhi  = new TH1F("hPosKaonPhi","PosKaonPhi", 40, -3.2, 3.2);
   TH1F *hPosKaonEmuPhi  = new TH1F("hPosKaonEmuPhi","PosKaonEmuPhi", 40, -3.2, 3.2);

   TH1F *hPosKaonEta  = new TH1F("hPosKaonEta","PosKaonEta", 40, -3.2, 3.2);
   TH1F *hPosKaonEmuEta  = new TH1F("hPosKaonEmuEta","PosKaonEmuEta", 40, -3.2, 3.2);

   TH1F *hNegKaonPt  = new TH1F("hNegKaonPt","NegKaonPt", 25, 0., 10.);
   TH1F *hNegKaonEmuPt  = new TH1F("hNegKaonEmuPt","NegKaonEmuPt", 25, 0., 10.);

    TH1F *hNegKaonPhi  = new TH1F("hNegKaonPhi","NegKaonPhi", 40, -3.2, 3.2);
   TH1F *hNegKaonEmuPhi  = new TH1F("hNegKaonEmuPhi","NegKaonEmuPhi", 40, -3.2, 3.2);

   TH1F *hNegKaonEta  = new TH1F("hNegKaonEta","NegKaonEta", 40, -3.2, 3.2);
   TH1F *hNegKaonEmuEta  = new TH1F("hNegKaonEmuEta","NegKaonEmuEta", 40, -3.2, 3.2);

   TH1F *hPhiMass  = new TH1F("hPhiMass","PhiMass", 20, 0.98, 1.13625); // calculated bin number by Alexx is 20
   TH1F *hPhiEmuMass  = new TH1F("hPhiEmuMass","PhiEmuMass", 20, 0.98, 1.13625);

   TH1F *hPhiPt  = new TH1F("hPhiPt","PhiPt", 25, 0., 10.);
   TH1F *hPhiEmuPt  = new TH1F("hPhiEmuPt","PhiEmuPt", 25, 0., 10.);

   TH1F *hPhiPhi  = new TH1F("hPhiPhi","PhiPhi", 40, -3.2, 3.2);
   TH1F *hPhiEmuPhi  = new TH1F("hPhiEmuPhi","PhiEmuPhi", 40, -3.2, 3.2);

   TH1F *hPhiEta  = new TH1F("hPhiEta","PhiEta", 40, -3.2, 3.2);
   TH1F *hPhiEmuEta  = new TH1F("hPhiEmuEta","PhiEmuEta", 40, -3.2, 3.2);

   TH1F *hPhiDeltaR  = new TH1F("hPhiDeltaR","PhiDeltaR", 30, 0., 3.0);
   TH1F *hPhiEmuDeltaR  = new TH1F("hPhiEmuDeltaR","PhiEmuDeltaR", 30, 0., 3.0);
   TH1F *hBsMass  = new TH1F("hBsMass","BsMass", 100, 5., 5.78125);
   TH1F *hBsEmuMass  = new TH1F("hBsEmuMass","BsEmuMass", 100, 5., 5.78125);

   TH1F *hTrkPairMass  = new TH1F("hTrkPairMass","TrkPairMass", 20, 0.98, 1.13625);
   TH1F *hTrkPairdR  = new TH1F("hTrkPairdR","TrkPairdR", 30, 0., 6.);
   TH1F *hTrkPairdEta  = new TH1F("hTrkPairdEta","TrkPairdEta", 20, 0., 4.0);
   TH1F *hTrkPairdPhi  = new TH1F("hTrkPairdPhi","TrkPairdPhi", 20, 0., 4.0);

   TH1F *hTrkPairEmuMass  = new TH1F("hTrkPairEmuMass","TrkPairEmuMass", 20, 0.98, 1.13625);
   TH1F *hTrkPairEmudR  = new TH1F("hTrkPairEmudR","TrkPairEmudR", 30, 0., 6.);
   TH1F *hTrkPairEmudEta  = new TH1F("hTrkPairEmudEta","TrkPairEmudEta", 20, 0., 4.0);
   TH1F *hTrkPairEmudPhi  = new TH1F("hTrkPairEmudPhi","TrkPairEmudPhi", 20, 0., 4.0);

   TH1F *hPhiPairMass  = new TH1F("hPhiPairMass","PhiPairMass", 100, 5., 5.78125);
   TH1F *hPhiPairdR  = new TH1F("hPhiPairdR","PhiPairdR", 20, 0., 5.);
   TH1F *hPhiPairdEta  = new TH1F("hPhiPairdEta","PhiPairdEta", 20, 0., 4.0);
   TH1F *hPhiPairdPhi  = new TH1F("hPhiPairdPhi","PhiPairdPhi", 20, 0., 4.0);

   TH1F *hPhiPairEmuMass  = new TH1F("hPhiPairEmuMass","PhiPairEmuMass", 100, 5., 5.78125);
   TH1F *hPhiPairEmudR  = new TH1F("hPhiPairEmudR","PhiPairEmudR", 20, 0., 5.);
   TH1F *hPhiPairEmudEta  = new TH1F("hPhiPairEmudEta","PhiPairEmudEta", 20, 0., 4.0);
   TH1F *hPhiPairEmudPhi  = new TH1F("hPhiPairEmudPhi","PhiPairEmudPhi", 20, 0., 4.0);

   TFile f("demo.root","recreate");


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //      std::cout << "Total number of bs's : " << trkbscands_eta->size() << std::endl;
      //std::cout << "Total number of bs  emulation's : " << trkbscandsemulation_eta->size() << std::endl;
      
      int Nposkaon = 0;
      TLorentzVector PosKaon;
      std::vector<TLorentzVector> PosKaonVec;
      for (int i = 0; i < trk_poskaon_gtt_selected_index->size(); i++){
	if (trk_gtt_selected_index->at(i) == -1) continue;
	else {
	  if (trk_poskaon_gtt_selected_index->at(i) == -1) continue;
	  else {
	    Nposkaon++;
 	    PosKaon.SetPtEtaPhiM(trk_pt->at(i), trk_eta->at(i), trk_phi->at(i), 0.494);
	    PosKaonVec.push_back(PosKaon);
	    hPosKaonPt->Fill(PosKaon.Pt());
	    hPosKaonPhi->Fill(PosKaon.Phi());
	    hPosKaonEta->Fill(PosKaon.Eta());
	  }
	}
      }

      int Nnegkaon = 0;
      TLorentzVector NegKaon;
      std::vector<TLorentzVector> NegKaonVec;
      for (int i = 0; i < trk_negkaon_gtt_selected_index->size(); i++){
	if (trk_gtt_selected_index->at(i) == -1) continue;
	else {
	  if (trk_negkaon_gtt_selected_index->at(i) == -1) continue;
	  else {
	    Nnegkaon++;
 	    NegKaon.SetPtEtaPhiM(trk_pt->at(i), trk_eta->at(i), trk_phi->at(i), 0.494);
	    NegKaonVec.push_back(NegKaon);
	    hNegKaonPt->Fill(NegKaon.Pt());
	    hNegKaonPhi->Fill(NegKaon.Phi());
	    hNegKaonEta->Fill(NegKaon.Eta());
	  }
	}
      }

      for (int i = 0; i < PosKaonVec.size() ; i++){
	for (int j = 0; j < NegKaonVec.size() ; j++){
	  TLorentzVector trkpairvec = PosKaonVec[i] + NegKaonVec[j];
	  hTrkPairMass->Fill(trkpairvec.M());
	  hTrkPairdR->Fill(PosKaonVec[i].DeltaR(NegKaonVec[j]));
	  double delEta = sqrt(pow(PosKaonVec[i].DeltaR(NegKaonVec[j]),2) - pow(PosKaonVec[i].DeltaPhi(NegKaonVec[j]),2));
	  hTrkPairdEta->Fill(std::fabs(delEta));
	  hTrkPairdPhi->Fill(std::fabs(PosKaonVec[i].DeltaPhi(NegKaonVec[j])));
	}
      }


      int NposkaonEmu = 0;
      TLorentzVector PosKaonEmu;
      std::vector<TLorentzVector> PosKaonEmuVec;
      for (int i = 0; i < trk_poskaon_gtt_selected_emulation_index->size(); i++){
	if (trk_gtt_selected_index->at(i) == -1) continue;
	else {
	  if (trk_poskaon_gtt_selected_emulation_index->at(i) == -1) continue;
	  else {
	    NposkaonEmu++;
 	    PosKaonEmu.SetPtEtaPhiM(trk_pt_emu->at(i), trk_eta_emu->at(i), trk_phi_emu->at(i), 0.494);
	    PosKaonEmuVec.push_back(PosKaonEmu);
	    hPosKaonEmuPt->Fill(PosKaonEmu.Pt());
	    hPosKaonEmuPhi->Fill(PosKaonEmu.Phi());
	    hPosKaonEmuEta->Fill(PosKaonEmu.Eta());
	  }
	}
      }

      int NnegkaonEmu = 0;
      TLorentzVector NegKaonEmu;
      std::vector<TLorentzVector> NegKaonEmuVec;
      for (int i = 0; i < trk_negkaon_gtt_selected_emulation_index->size(); i++){
	if (trk_gtt_selected_index->at(i) == -1) continue;
	else {
	  if (trk_negkaon_gtt_selected_emulation_index->at(i) == -1) continue;
	  else {
	    NnegkaonEmu++;
 	    NegKaonEmu.SetPtEtaPhiM(trk_pt_emu->at(i), trk_eta_emu->at(i), trk_phi_emu->at(i), 0.494);
	    NegKaonEmuVec.push_back(NegKaonEmu);
	    hNegKaonEmuPt->Fill(NegKaonEmu.Pt());
	    hNegKaonEmuPhi->Fill(NegKaonEmu.Phi());
	    hNegKaonEmuEta->Fill(NegKaonEmu.Eta());
	  }
	}
      }

      for (int i = 0; i < PosKaonEmuVec.size() ; i++){
	for (int j = 0; j < NegKaonEmuVec.size() ; j++){
	  TLorentzVector trkpairemuvec = PosKaonEmuVec[i] + NegKaonEmuVec[j];
	  hTrkPairEmuMass->Fill(trkpairemuvec.M());
	  hTrkPairEmudR->Fill(PosKaonEmuVec[i].DeltaR(NegKaonEmuVec[j]));
	  double delEta = sqrt(pow(PosKaonEmuVec[i].DeltaR(NegKaonEmuVec[j]),2) - pow(PosKaonEmuVec[i].DeltaPhi(NegKaonEmuVec[j]),2));
	  hTrkPairEmudEta->Fill(std::fabs(delEta));
	  hTrkPairEmudPhi->Fill(std::fabs(PosKaonEmuVec[i].DeltaPhi(NegKaonEmuVec[j])));
	}
      }


      //      std::cout << "Nposkaon : " << Nposkaon << "\tNposkaonEmu : " << NposkaonEmu << "\tNnegkaon : " << Nnegkaon << "\tNnegkaonEmu : " << Nnegkaon << std::endl; 

      //std::cout << "Pos track pt : " << PosKaon.Pt() << "\t phi : " << PosKaon.Phi() << "\t emu pt : " << PosKaonEmu.Pt() << "\t emu phi : " << PosKaonEmu.Phi() << std::endl;
      //std::cout << "Neg track pt : " << NegKaon.Pt() << "\t phi : " << NegKaon.Phi() << "\t emu pt : " << NegKaonEmu.Pt() << "\t emu phi : " << NegKaonEmu.Phi() << std::endl;


      //std::cout << "Total number of phi's : " << trkphicands_eta->size() << std::endl;
      //std::cout << "Total number of phi emulation's : " << trkphicandsemulation_eta->size() << std::endl;
      vector<TLorentzVector> PhiCandVec, PhiCandEmuVec;

      //      if (trkphicands_eta->size() > 0) {
      for (int i = 0; i < trkphicands_eta->size(); i++){
        TLorentzVector PhiCand;
	std::cout << "trk phi cands eta : " << trkphicands_eta->at(i) << "\ttrk phi cands phi : " << trkphicands_phi->at(i) << "\ttrk phi cands pt : " << trkphicands_pt->at(i) << "\ttrk phi cands mass : " << trkphicands_mass->at(i) << std::endl;
	/////std::cout << "trk phi cands mass from ntuple : " << trkphicands_mass->at(i) << std::endl;
        PhiCand.SetPtEtaPhiM(trkphicands_pt->at(i),trkphicands_eta->at(i),trkphicands_phi->at(i), trkphicands_mass->at(i));
        PhiCandVec.push_back(PhiCand);
	std::cout << "trk phi cands mass : " << PhiCand.M() << std::endl;
        hPhiMass->Fill(PhiCand.M());
	hPhiPt->Fill(PhiCand.Pt());
	hPhiEta->Fill(PhiCand.Eta());
	hPhiPhi->Fill(PhiCand.Phi());
	std::cout << "one phi loop ended" << std::endl;
      }

      for (int i = 0; i < PhiCandVec.size() ; i++){
	for (int j = i+1; j < PhiCandVec.size() ; j++){
	  TLorentzVector phipairvec = PhiCandVec[i] + PhiCandVec[j];
	  hPhiPairMass->Fill(phipairvec.M());
	  hPhiPairdR->Fill(PhiCandVec[i].DeltaR(PhiCandVec[j]));
	  double delEta = sqrt(pow(PhiCandVec[i].DeltaR(PhiCandVec[j]),2) - pow(PhiCandVec[i].DeltaPhi(PhiCandVec[j]),2));
	  hPhiPairdEta->Fill(std::fabs(delEta));
	  hPhiPairdPhi->Fill(std::fabs(PhiCandVec[i].DeltaPhi(PhiCandVec[j])));
	}
      }


      for (int i = 0; i < trkphicandsemulation_eta->size(); i++){
        TLorentzVector PhiCandEmu;
	//	std::cout << "trk phi cands emulation eta : " << trkphicandsemulation_eta->at(i) << std::endl;
	std::cout << "trk phi cands emu eta : " << trkphicandsemulation_eta->at(i) << "\ttrk phi cands phi : " << trkphicandsemulation_phi->at(i) << "\ttrk phi cands pt : " << trkphicandsemulation_pt->at(i) << "\ttrk phi cands mass : " << trkphicandsemulation_mass->at(i) << std::endl;
	PhiCandEmu.SetPtEtaPhiM(trkphicandsemulation_pt->at(i),trkphicandsemulation_eta->at(i),trkphicandsemulation_phi->at(i),trkphicandsemulation_mass->at(i));
        PhiCandEmuVec.push_back(PhiCandEmu);
	std::cout << "trk phi cands emulation mass : " << PhiCandEmu.M() << std::endl;
        hPhiEmuMass->Fill(PhiCandEmu.M());
	hPhiEmuPt->Fill(PhiCandEmu.Pt());
	hPhiEmuEta->Fill(PhiCandEmu.Eta());
	hPhiEmuPhi->Fill(PhiCandEmu.Phi());
	std::cout << "one emulation phi loop ended" << std::endl;
      }

      std::cout << "phicand vec size : " << PhiCandVec.size() << "\tphican emu vec size : " << PhiCandEmuVec.size() << std::endl;
      if (PhiCandVec.size() > 0 && PhiCandEmuVec.size() < 1) {
	for (int iphi = 0 ; iphi < PhiCandVec.size() ; iphi++){
	  std::cout << "extra selected Phi cand sim pt : " << PhiCandVec[iphi].Pt() << "\teta : " << PhiCandVec[iphi].Eta() << "\tphi : " << PhiCandVec[iphi].Phi() << std::endl;
	}
      }


      for (int i = 0; i < PhiCandEmuVec.size() ; i++){
	for (int j = i+1; j < PhiCandEmuVec.size() ; j++){
	  TLorentzVector phipairemuvec = PhiCandEmuVec[i] + PhiCandEmuVec[j];
	  hPhiPairEmuMass->Fill(phipairemuvec.M());
	  hPhiPairEmudR->Fill(PhiCandEmuVec[i].DeltaR(PhiCandEmuVec[j]));
	  double delEta = sqrt(pow(PhiCandEmuVec[i].DeltaR(PhiCandEmuVec[j]),2) - pow(PhiCandEmuVec[i].DeltaPhi(PhiCandEmuVec[j]),2));
	  hPhiPairEmudEta->Fill(std::fabs(delEta));
	  hPhiPairEmudPhi->Fill(std::fabs(PhiCandEmuVec[i].DeltaPhi(PhiCandEmuVec[j])));
	}
      }

      //      FillDeltaR(phicandVec , hPhiDeltaR);
      //std::cout << "end of simulation deltaR" << std::endl;
      //FillDeltaR(phicandEmuVec , hPhiEmuDeltaR);

   }
   
   hPosKaonPt->Write();
   hPosKaonEmuPt->Write();   
   hNegKaonPt->Write();
   hNegKaonEmuPt->Write();
   hPosKaonEta->Write();
   hPosKaonEmuEta->Write();   
   hNegKaonEta->Write();
   hNegKaonEmuEta->Write();
   hPosKaonPhi->Write();
   hPosKaonEmuPhi->Write();
   hNegKaonPhi->Write();
   hNegKaonEmuPhi->Write();
   hPhiMass->Write();
   hPhiEmuMass->Write();
   hPhiPt->Write();
   hPhiEmuPt->Write();   
   hPhiPhi->Write();
   hPhiEmuPhi->Write();
   hPhiEta->Write();
   hPhiEmuEta->Write();
   hPhiDeltaR->Write();
   hPhiEmuDeltaR->Write();
   hBsMass->Write();
   hBsEmuMass->Write();

   hTrkPairMass->Write();
   hTrkPairdR->Write();
   hTrkPairdEta->Write();
   hTrkPairdPhi->Write();
   hTrkPairEmuMass->Write();
   hTrkPairEmudR->Write();
   hTrkPairEmudEta->Write();
   hTrkPairEmudPhi->Write();
   
   hPhiPairMass->Write();
   hPhiPairdR->Write();
   hPhiPairdEta->Write();
   hPhiPairdPhi->Write();
   hPhiPairEmuMass->Write();
   hPhiPairEmudR->Write();
   hPhiPairEmudEta->Write();
   hPhiPairEmudPhi->Write();



}

void bstophiphi_phase2::FillDeltaR(std::vector<TLorentzVector> CandVec, TH1F* h){
  if(CandVec.size() > 1) {
    for (int i = 0; i < CandVec.size() ; i ++){
      for (int j = i+1 ; j < CandVec.size(); j++){
	//std::cout << "DeltaR value is : " << CandVec[i].DeltaR(CandVec[j]) << std::endl;
	h->Fill(CandVec[i].DeltaR(CandVec[j]));
      }
    }  
  }
}

double bstophiphi_phase2::PhiRangeConv(double glbphi){
  if (glbphi < 0.) {
    glbphi = glbphi + 2*M_PI; 
  }
  else if (glbphi > 2*M_PI) {
    glbphi = glbphi - 2*M_PI;
  }
  return glbphi;
}
