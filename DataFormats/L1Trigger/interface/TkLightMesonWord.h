#ifndef FIRMWARE_TkLightMesonWord_h
#define FIRMWARE_TkLightMesonWord_h

#include <vector>
#include <ap_int.h>
#include <cassert>
#include <cmath>
#include <bitset>
#include <string>

namespace l1t {

  class TkLightMesonWord {
  public:
    // ----------constants, enums and typedefs ---------
    int INTPHI_PI = 720;
    int INTPHI_TWOPI = 2 * INTPHI_PI;
    float INTPT_LSB = 1 >> 5;
    double ETAPHI_LSB = M_PI / (1 << 12);
    double Z0_LSB = 0.05;

    enum TkLightMesonTypes{
      kPhiType = 1,
      kRhoType = 2,
      kBsType = 3,
    };

    enum TkLightMesonBitWidths {
      kValidSize = 1,
      kPtSize = 16,
      kPtMagSize = 11,
      kGlbPhiSize = 13,
      //kGlbPhiMagSize = 3, // signed  // DEl
      kGlbEtaSize = 14,
      //kGlbEtaMagSize = 4, // signed // del
      kZ0Size = 10,
      //kZ0MagSize = 6, //signed // del
      kMassSize = 10,
      kMassMagSize = 3,
      kTypeSize = 2,
      kNtracksSize = 3,
      kUnassignedSize = 27,
      kTkLightMesonWordSize = kValidSize + kPtSize + kGlbPhiSize + kGlbEtaSize + kZ0Size + kMassSize + kTypeSize + kNtracksSize + kUnassignedSize,
    };

    enum TkLightMesonBitLocations {
      kValidLSB = 0,
      kValidMSB = kValidLSB + TkLightMesonBitWidths::kValidSize - 1,
      kPtLSB = kValidMSB + 1,
      kPtMSB = kPtLSB + TkLightMesonBitWidths::kPtSize - 1,
      kGlbPhiLSB = kPtMSB + 1,
      kGlbPhiMSB = kGlbPhiLSB + TkLightMesonBitWidths::kGlbPhiSize - 1,
      kGlbEtaLSB = kGlbPhiMSB + 1,
      kGlbEtaMSB = kGlbEtaLSB + TkLightMesonBitWidths::kGlbEtaSize - 1,
      kZ0LSB = kGlbEtaMSB + 1,
      kZ0MSB = kZ0LSB + TkLightMesonBitWidths::kZ0Size - 1,
      kMassLSB = kZ0MSB + 1,
      kMassMSB = kMassLSB + TkLightMesonBitWidths::kMassSize - 1,
      kTypeLSB = kMassMSB + 1,
      kTypeMSB = kTypeLSB + TkLightMesonBitWidths::kTypeSize - 1,
      kNtracksLSB = kTypeMSB + 1,
      kNtracksMSB = kNtracksLSB + TkLightMesonBitWidths::kNtracksSize - 1,
      kUnassignedLSB = kNtracksMSB + 1,
      kUnassignedMSB = kUnassignedLSB + TkLightMesonBitWidths::kUnassignedSize - 1,
    };

    typedef ap_uint<TkLightMesonBitWidths::kValidSize> valid_t;
    typedef ap_ufixed<TkLightMesonBitWidths::kPtSize, TkLightMesonBitWidths::kPtMagSize, AP_RND_CONV, AP_SAT> pt_t;
    typedef ap_int<TkLightMesonBitWidths::kGlbPhiSize> glbphi_t;
    typedef ap_int<TkLightMesonBitWidths::kGlbEtaSize> glbeta_t;
    typedef ap_int<TkLightMesonBitWidths::kZ0Size> z0_t;     // 40cm / 0.1
    typedef ap_ufixed<TkLightMesonBitWidths::kMassSize, TkLightMesonBitWidths::kMassMagSize, AP_RND_CONV, AP_SAT> mass_t;
    typedef ap_uint<TkLightMesonBitWidths::kTypeSize> type_t;       //type of meson
    typedef ap_uint<TkLightMesonBitWidths::kNtracksSize> ntracks_t;                                       //number of tracks
    typedef ap_uint<TkLightMesonBitWidths::kUnassignedSize> unassigned_t;  // Unassigned bits
    typedef std::bitset<TkLightMesonBitWidths::kTkLightMesonWordSize> tklightmesonword_bs_t;
    typedef ap_uint<TkLightMesonBitWidths::kTkLightMesonWordSize> tklightmesonword_t;

  public:
    // ----------Constructors --------------------------
    TkLightMesonWord() {}
    TkLightMesonWord(valid_t valid, pt_t pt, glbphi_t phi, glbeta_t eta, z0_t z0, mass_t mass, type_t type, ntracks_t ntracks, unassigned_t unassigned); //{
      /*      std::string word = "";
      word.append(TkLightMesonBitWidths::kUnassignedSize - (unassigned.to_string().length() - 2), '0');
      word = word + (unassigned.to_string().substr(2, unassigned.to_string().length() - 2));
      word.append(TkLightMesonBitWidths::kNtracksSize - (ntracks.to_string().length() - 2), '0');
      word = word + (ntracks.to_string().substr(2, ntracks.to_string().length() - 2));
      word.append(TkLightMesonBitWidths::kTypeSize - (type.to_string().length() - 2), '0');
      word = word + (type.to_string().substr(2, type.to_string().length() - 2));
      word.append(TkLightMesonBitWidths::kMassSize - (mass.to_string().length() - 2), '0');
      word = word + (mass.to_string().substr(2, mass.to_string().length() - 2));
      word.append(TkLightMesonBitWidths::kZ0Size - (z0.to_string().length() - 2), '0');
      word = word + (z0.to_string().substr(2, z0.to_string().length() - 2));
      word.append(TkLightMesonBitWidths::kGlbEtaSize - (eta.to_string().length() - 2), '0');
      word = word + (eta.to_string().substr(2, eta.to_string().length() - 2));
      word.append(TkLightMesonBitWidths::kGlbPhiSize - (phi.to_string().length() - 2), '0');
      word = word + (phi.to_string().substr(2, phi.to_string().length() - 2));
      ap_ufixed<kPtSize + 5, kPtMagSize + 5, AP_TRN, AP_SAT> pt_2 = pt;
      ap_uint<kPtSize> pt_temp = pt_2 << 5;
      word.append(TkLightMesonBitWidths::kPtSize - (pt_temp.to_string().length() - 2), '0');
      word = word + (pt_temp.to_string().substr(2, pt_temp.to_string().length() - 2));
      word.append(TkLightMesonBitWidths::kValidSize - (valid.to_string().length() - 2), '0');
      word = word + (valid.to_string().substr(2, valid.to_string().length() - 2));*/

      //      tklightmesonword_bs_t tmp(word);
      //tkLightMesonWord_ = tmp;
      // }

    ~TkLightMesonWord() {}

    // ----------copy constructor ----------------------
    TkLightMesonWord(const TkLightMesonWord& word) { tkLightMesonWord_ = word.tkLightMesonWord_; }

    // ----------operators -----------------------------
    TkLightMesonWord& operator=(const TkLightMesonWord& word) {
      tkLightMesonWord_ = word.tkLightMesonWord_;
      return *this;
    }

    // ----------member functions (getters) ------------
    // These functions return arbitarary precision words (lists of bits) for each quantity
    valid_t validWord() const { return tkLightMesonWord()(TkLightMesonBitLocations::kValidMSB, TkLightMesonBitLocations::kValidLSB); }
    pt_t ptWord() const {
      pt_t ret;
      ret.V = tkLightMesonWord()(TkLightMesonBitLocations::kPtMSB, TkLightMesonBitLocations::kPtLSB);
      return ret;
    }
    glbphi_t glbPhiWord() const {
      glbphi_t ret;
      ret.V = tkLightMesonWord()(TkLightMesonBitLocations::kGlbPhiMSB, TkLightMesonBitLocations::kGlbPhiLSB);
      return ret;
    }
    glbeta_t glbEtaWord() const {
      glbeta_t ret;
      ret.V = tkLightMesonWord()(TkLightMesonBitLocations::kGlbEtaMSB, TkLightMesonBitLocations::kGlbEtaLSB);
      return ret;
    }
    z0_t z0Word() const {
      z0_t ret;
      ret.V = tkLightMesonWord()(TkLightMesonBitLocations::kZ0MSB, TkLightMesonBitLocations::kZ0LSB);
      return ret;
    }
    mass_t massWord() const {
      mass_t ret;
      ret.V = tkLightMesonWord()(TkLightMesonBitLocations::kMassMSB, TkLightMesonBitLocations::kMassLSB);
      return ret;
    }
    type_t typeWord() const {
      type_t ret;
      ret.V = tkLightMesonWord()(TkLightMesonBitLocations::kTypeMSB, TkLightMesonBitLocations::kTypeLSB);
      return ret;
    }
    ntracks_t ntracksWord() const {
      ntracks_t ret;
      ret.V = tkLightMesonWord()(TkLightMesonBitLocations::kNtracksMSB, TkLightMesonBitLocations::kNtracksLSB);
      return ret;
    }
    unassigned_t unassignedWord() const {
      return tkLightMesonWord()(TkLightMesonBitLocations::kUnassignedMSB, TkLightMesonBitLocations::kUnassignedLSB);
    }
    tklightmesonword_t tkLightMesonWord() const { return tklightmesonword_t(tkLightMesonWord_.to_string().c_str(), 2); }

    // These functions return the packed bits in integer format for each quantity
    // Signed quantities have the sign enconded in the left-most bit.
    unsigned int validBits() const { return validWord().to_uint(); }
    unsigned int ptBits() const { return ptWord().to_uint(); }
    unsigned int glbPhiBits() const { return glbPhiWord().to_uint(); }
    unsigned int glbEtaBits() const { return glbEtaWord().to_uint(); }
    unsigned int z0Bits() const { return z0Word().to_uint(); }
    unsigned int massBits() const { return massWord().to_uint(); }
    unsigned int typeBits() const { return typeWord().to_uint(); }
    unsigned int ntracksBits() const { return ntracksWord().to_uint(); }
    unsigned int unassignedBits() const { return unassignedWord().to_uint(); }

    // These functions return the unpacked and converted values
    // These functions return real numbers converted from the digitized quantities by unpacking the 64-bit vertex word
    bool valid() const { return validWord().to_bool(); }
    float pt() const { return ptWord().to_float(); }
    float glbphi() const { return glbPhiWord().to_float() * ETAPHI_LSB; }
    float glbeta() const { return glbEtaWord().to_float() * ETAPHI_LSB; }
    float z0() const { return z0Word().to_float() * Z0_LSB; }
    float mass() const { return massWord().to_float(); }
    unsigned int type() const { return typeWord().to_uint(); }
    unsigned int ntracks() const { return ntracksWord().to_uint(); }
    unsigned int unassigned() const { return unassignedWord().to_uint(); }

    // ----------member functions (setters) ------------
    void setTkLightMesonWord(valid_t valid, pt_t pt, glbphi_t phi,  glbeta_t eta, z0_t z0, mass_t mass, type_t type, ntracks_t ntracks, unassigned_t unassigned);

  private:
    // ----------private member functions --------------
    double unpackSignedValue(unsigned int bits, unsigned int nBits, double lsb) const {
      int isign = 1;
      unsigned int digitized_maximum = (1 << nBits) - 1;
      if (bits & (1 << (nBits - 1))) {  // check the sign
        isign = -1;
        bits = (1 << (nBits + 1)) - bits;  // if negative, flip everything for two's complement encoding
      }
      return (double(bits & digitized_maximum) + 0.5) * lsb * isign;
    }

    //    double dRTrkPair() const;
    //double dmass() const;
    //double dxyTrkPair() const;
    //double dzTrkPair() const;
    //double d0TrkPair() const;

    // ----------member data ---------------------------
    tklightmesonword_bs_t tkLightMesonWord_;
  };

  typedef std::vector<l1t::TkLightMesonWord> TkLightMesonWordCollection;

}  // namespace l1t

#endif
