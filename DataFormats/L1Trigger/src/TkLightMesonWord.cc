////////
//
// class to store the 96-bit track word produced by the L1 Track Trigger.  Intended to be inherited by L1 TTTrack.
// packing scheme given below.
//
// author:      Alexx Perloff
// created:     March 17, 2021
//
///////

#include "DataFormats/L1Trigger/interface/TkLightMesonWord.h"

namespace l1t {

  TkLightMesonWord::TkLightMesonWord(valid_t valid, pt_t pt, glbphi_t phi, glbeta_t eta, z0_t z0, mass_t mass, type_t type, ntracks_t ntracks, unassigned_t unassigned){
    // convert directly to AP types
    valid_t valid_ap = valid;
    pt_t pt_ap = pt;
    glbphi_t phi_ap = phi;
    glbeta_t eta_ap = eta;
    z0_t z0_ap = z0;
    mass_t mass_ap = mass;
    type_t type_ap = type;
    ntracks_t ntracks_ap = ntracks;
    unassigned_t unassigned_ap = unassigned;

    setTkLightMesonWord(valid_ap, pt_ap, phi_ap, eta_ap, z0_ap, mass_ap, type_ap, ntracks_ap, unassigned_ap);
  }

  void TkLightMesonWord::setTkLightMesonWord(valid_t valid, pt_t pt, glbphi_t phi,  glbeta_t eta, z0_t z0, mass_t mass, type_t type, ntracks_t ntracks, unassigned_t unassigned){
    // pack the TkLightMesonWord
    unsigned int offset = 0;
    for (unsigned int b = offset; b < (offset + TkLightMesonBitWidths::kValidSize); b++) {
      tkLightMesonWord_.set(b, valid[b - offset]);
    }
    offset += TkLightMesonBitWidths::kValidSize;
    
    for (unsigned int b = offset; b < (offset + TkLightMesonBitWidths::kPtSize); b++) {
      tkLightMesonWord_.set(b, z0[b - offset]);
    }
    offset += TkLightMesonBitWidths::kPtSize;

    for (unsigned int b = offset; b < (offset + TkLightMesonBitWidths::kGlbPhiSize); b++) {
      tkLightMesonWord_.set(b, z0[b - offset]);
    }
    offset += TkLightMesonBitWidths::kGlbPhiSize;
    
    for (unsigned int b = offset; b < (offset + TkLightMesonBitWidths::kGlbEtaSize); b++) {
      tkLightMesonWord_.set(b, z0[b - offset]);
    }
    offset += TkLightMesonBitWidths::kGlbEtaSize;

    for (unsigned int b = offset; b < (offset + TkLightMesonBitWidths::kZ0Size); b++) {
      tkLightMesonWord_.set(b, z0[b - offset]);
    }
    offset += TkLightMesonBitWidths::kZ0Size;

    for (unsigned int b = offset; b < (offset + TkLightMesonBitWidths::kMassSize); b++) {
      tkLightMesonWord_.set(b, z0[b - offset]);
    }
    offset += TkLightMesonBitWidths::kMassSize;

    for (unsigned int b = offset; b < (offset + TkLightMesonBitWidths::kTypeSize); b++) {
      tkLightMesonWord_.set(b, z0[b - offset]);
    }
    offset += TkLightMesonBitWidths::kTypeSize;

    for (unsigned int b = offset; b < (offset + TkLightMesonBitWidths::kNtracksSize); b++) {
      tkLightMesonWord_.set(b, z0[b - offset]);
    }
    offset += TkLightMesonBitWidths::kNtracksSize;

    for (unsigned int b = offset; b < (offset + TkLightMesonBitWidths::kUnassignedSize); b++) {
      tkLightMesonWord_.set(b, unassigned[b - offset]);
    }
  }

}  // namespace l1t
