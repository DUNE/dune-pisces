#pragma once

#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/IPrediction.h"

#include <Eigen/Dense>

namespace pisces {

  using namespace ana;

  // define const arrays - kinda inelegant but it works
  namespace oscchan { 
    const std::array<std::string, 6> kFlavNames
    {
      "nuetonue",
      "numutonumu",
      "nuetonumu",
      "nuetonutau",
      "numutonue",
      "numutonutau"
    };
    const std::array<Flavors::Flavors_t, 6> kFlavours
    {
      Flavors::kNuEToNuE,
      Flavors::kNuMuToNuMu,
      Flavors::kNuEToNuMu,
      Flavors::kNuEToNuTau,
      Flavors::kNuMuToNuE,
      Flavors::kNuMuToNuTau
    };
    const std::array<Loaders::SwappingConfig, 6> kSwapConfigs
    {
      Loaders::kNonSwap,
      Loaders::kNonSwap,
      Loaders::kNueSwap,
      Loaders::kNuTauSwap,
      Loaders::kNueSwap,
      Loaders::kNuTauSwap
    };
    const std::array<Cut, 6> kFlavCuts
    {
      kIsBeamNue,
      kIsNumuCC,
      kIsNumuApp,
      kIsTauFromE,
      kIsSig,
      kIsTauFromMu
    };
    const std::array<int, 6> kInitFlav { 12, 14, 12, 12, 14, 14 };
    const std::array<int, 6> kFinalFlav { 12, 14, 14, 16, 12, 16 };

    const std::array<std::string, 2> kSignNames { "nu", "nubar" };
    const std::array<Sign::Sign_t, 2> kSigns { Sign::kNu, Sign::kAntiNu };
    const std::array<Cut, 2> kSignCuts { !kIsAntiNu, kIsAntiNu };
  } // namespace oscchan

  class OscChannel {

  public:

    OscChannel(const std::string& name);
    ~OscChannel() {};

    /// Get name of oscillation channel
    std::string Name() const { return fName; };
    /// Get flavour of oscillation channel
    Flavors::Flavors_t Flav() const { return fFlavour; };
    /// Get current of oscillation channel
    Current::Current_t Curr() const { return fCurrent; };
    /// Get sign of oscillation channel
    Sign::Sign_t Sign() const { return fSign; };
    /// Get cut for current oscillation channel
    const Cut TruthCut() const { return fCut; };
    /// Get swap config for current oscillation channel
    Loaders::SwappingConfig Config() const { return fConfig; };
    /// Get initial neutrino flavour PDG
    int From() const { return fFrom; };
    /// Get final neutrino flavour PDG
    int To() const { return fTo; };

    bool operator <(const OscChannel& rhs) const { return Name() < rhs.Name(); };

    virtual void SaveTo(TDirectory* dir, const std::string& name) const;
    static std::unique_ptr<OscChannel> LoadFrom(TDirectory* dir,
                                                const std::string& name);

  protected:

    static std::string CCName(size_t f, size_t s) { return "cc_"+oscchan::kSignNames[s]+"_"+oscchan::kFlavNames[f]; };
    static std::string NCName() { return "nc"; };

    std::string                fName;
    Flavors::Flavors_t         fFlavour;
    Current::Current_t         fCurrent;
    Sign::Sign_t               fSign;
    Cut                        fCut;
    Loaders::SwappingConfig    fConfig;
    int                        fFrom;
    int                        fTo;

    friend class Sample;

  }; // class OscChannel

} // namespace pisces
