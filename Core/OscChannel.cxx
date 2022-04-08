#include "CAFAna/PISCES/OscChannel.h"

#include "TDirectory.h"
#include "TObjString.h"

namespace pisces {

  using std::make_unique;
  using std::string;
  using std::unique_ptr;
  using std::vector;
  using Eigen::ArrayXd;
  using namespace ana;
  using namespace oscchan;

  //-------------------------------------------------------------------------
  OscChannel::OscChannel(string const& name)
    : fName(name), fFlavour(Flavors::kAll), fCurrent(Current::kBoth),
      fSign(Sign::kBoth), fCut(kNoCut)
  {
    // Add each combination of flavour and sign
    for (size_t i = 0; i < flavNames.size(); ++i) {
      for (size_t j = 0; j < signNames.size(); ++j) {
        if (name == OscChannel::CCName(i,j)) {
          fFlavour = flavours[i];
          fCurrent = Current::kCC;
          fSign = signs[j];
          fCut = flavCuts[i] && signCuts[j];
          fConfig = swapConfigs[i];
          fFrom = pow(-1, j) * initFlav[i];
          fTo = pow(-1, j) * finalFlav[i];
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

  //-------------------------------------------------------------------------
  vector<OscChannel> OscChannel::OscChannels(Sample const& sample)
  {
    vector<OscChannel> ret;
    size_t nFlavs = sample.IsFD() ? 6 : 2;
    // Add each combination of flavour and sign
    for (size_t i = 0; i < nFlavs; ++i)
      for (size_t j = 0; j < 2; ++j)
        ret.push_back(OscChannel(CCName(i,j)));
    // Add the neutral currents
    ret.push_back(OscChannel(NCName()));
    return ret;
  } // function OscChannel::GetFromName
  
  //-------------------------------------------------------------------------
  ArrayXd OscChannel::Predict(Sample const& sample, osc::IOscCalc* osc,
                              SystShifts shifts) const
  {
    Spectrum spec = sample.PredictComponentSyst(osc, sample.Shifts(shifts), Flav(), Curr(), Sign());
    ArrayXd arr = spec.GetEigen(sample.POT());
    return arr.segment(1, sample.GetBinning().NBins());
  } // function OscChannel::Predict

  //-------------------------------------------------------------------------
  void OscChannel::SaveTo(TDirectory* dir, string const& name) const
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

  //-------------------------------------------------------------------------
  unique_ptr<OscChannel> OscChannel::LoadFrom(TDirectory* dir,
                                              string const& name)
  {
    dir = dir->GetDirectory(name.c_str());
    assert(dir);

    TObjString* type = (TObjString*)dir->Get("type");
    assert(type);
    assert(type->GetString() == "OscChannel");

    TObjString* n = (TObjString*)dir->Get("name");
    assert(n);

    return make_unique<OscChannel>(n->GetString().Data());
  } // function OscChannal::LoadFrom

} // namespace pisces
