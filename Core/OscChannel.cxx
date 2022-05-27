#include "Core/OscChannel.h"

#include "TDirectory.h"
#include "TObjString.h"

namespace pisces {

  using namespace ana;

  //---------------------------------------------------------------------------
  OscChannel::OscChannel(const std::string& name)
    : fName(name), fFlavour(Flavors::kAll), fCurrent(Current::kBoth),
      fSign(Sign::kBoth), fCut(kNoCut)
  {
    // Add each combination of flavour and sign
    for (size_t i = 0; i < oscchan::kFlavNames.size(); ++i) {
      for (size_t j = 0; j < oscchan::kSignNames.size(); ++j) {
        if (name == OscChannel::CCName(i,j)) {
          fFlavour = oscchan::kFlavours[i];
          fCurrent = Current::kCC;
          fSign = oscchan::kSigns[j];
          fCut = oscchan::kFlavCuts[i] && oscchan::kSignCuts[j];
          fConfig = oscchan::kSwapConfigs[i];
          fFrom = pow(-1, j) * oscchan::kInitFlav[i];
          fTo = pow(-1, j) * oscchan::kFinalFlav[i];
          return;
        } // if name matches
      } // for flavour
    } // for sign
    if (name == OscChannel::NCName()) {
      fFlavour = Flavors::kAll;
      fCurrent = Current::kNC;
      fSign = Sign::kBoth;
      fCut = kIsNC;
      fConfig = Loaders::kNonSwap;
      fFrom = 12;
      fTo = 0;
      return;
    }
    assert(false && ("Sample "+name+" not recognised!").c_str());
  } // OscChannel constructor

  //---------------------------------------------------------------------------
  void OscChannel::SaveTo(TDirectory* dir, const std::string& name) const
  {
    TDirectory* tmp = gDirectory;
    dir = dir->mkdir(name.c_str());
    dir->cd();

    TObjString("OscChannel").Write("type");
    TObjString(fName.c_str()).Write("name");
    dir->Write();

    delete dir;
    tmp->cd();
  } // function OscChannel::SaveTo

  //---------------------------------------------------------------------------
  std::unique_ptr<OscChannel> OscChannel::LoadFrom(TDirectory* dir,
                                                   const std::string& name)
  {
    dir = dir->GetDirectory(name.c_str());
    assert(dir);

    [[maybe_unused]] TObjString* type = (TObjString*)dir->Get("type");
    assert(type);
    assert(type->GetString() == "OscChannel");

    TObjString* n = (TObjString*)dir->Get("name");
    assert(n);

    return std::make_unique<OscChannel>(n->GetString().Data());
  } // function OscChannal::LoadFrom

} // namespace pisces
