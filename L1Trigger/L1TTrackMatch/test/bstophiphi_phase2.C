#define bstophiphi_phase2_cxx
#include "bstophiphi_phase2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

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
   TH1F *hPhiMass  = new TH1F("hPhiMass","hPhiMass", 50, 0.95, 1.05);
   TH1F *hPhiEmuMass  = new TH1F("hPhiEmuMass","hPhiEmuMass", 50, 0.95, 1.05);
   TFile f("demo.root","recreate");


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //      std::cout << "Total number of bs's : " << trkbscands_eta->size() << std::endl;
      //std::cout << "Total number of bs  emulation's : " << trkbscandsemulation_eta->size() << std::endl;

      std::cout << "Total number of phi's : " << trkphicands_eta->size() << std::endl;
      std::cout << "Total number of phi emulation's : " << trkphicandsemulation_eta->size() << std::endl;
      vector<TLorentzVector> phicandVec, phicandEmuVec;

      //      if (trkphicands_eta->size() > 0) {
      for (int i = 0; i < trkphicands_eta->size(); i++){
        TLorentzVector PhiCand;
	std::cout << "trk phi cands eta : " << trkphicands_eta->at(i) << std::endl;
	std::cout << "trk phi cands mass from ntuple : " << trkphicands_mass->at(i) << std::endl;
        PhiCand.SetPtEtaPhiM(trkphicands_pt->at(i),trkphicands_eta->at(i),trkphicands_phi->at(i), trkphicands_mass->at(i));
        phicandVec.push_back(PhiCand);
	std::cout << "trk phi cands mass : " << PhiCand.M() << std::endl;
        hPhiMass->Fill(PhiCand.M());
	std::cout << "one phi loop ended" << std::endl;
      }

      for (int i = 0; i < trkphicandsemulation_eta->size(); i++){
        TLorentzVector PhiCandEmu;
	std::cout << "trk phi cands emulation eta : " << trkphicandsemulation_eta->at(i) << std::endl;
	std::cout << "one phi emulation loop ended" << std::endl;
        PhiCandEmu.SetPtEtaPhiM(trkphicandsemulation_pt->at(i),trkphicandsemulation_eta->at(i),trkphicandsemulation_phi->at(i),trkphicandsemulation_mass->at(i));
        phicandEmuVec.push_back(PhiCandEmu);
	std::cout << "trk phi cands emulation mass : " << PhiCandEmu.M() << std::endl;
        hPhiEmuMass->Fill(PhiCandEmu.M());
	std::cout << "one emulation phi loop ended" << std::endl;
      }
   }
   
   hPhiMass->Write();
   hPhiEmuMass->Write();

}
