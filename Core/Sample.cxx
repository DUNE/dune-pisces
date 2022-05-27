#include "Core/Sample.h"
#include "Core/OscChannel.h"

#include <iostream>

namespace pisces {

  using namespace ana;

  //-------------------------------------------------------------------------  
  void Error(const std::string& msg)
  {
    std::cerr << std::endl << "  error: " << msg << std::endl << std::endl;
    abort();
  } // function Error

  //-------------------------------------------------------------------------
  void Usage(char** argv, const std::string& msg)
  {
    std::cerr << std::endl << "  usage: " << argv[0];
    if (msg.empty()) std::cerr << " [no args]";
    else std::cerr << " " << msg;
    std::cerr << std::endl << std::endl;
    abort();
  } // function Usage

  //-------------------------------------------------------------------------
  Sample::Sample(Selection s, Polarity p, Detector d)
    : fSel(s), fPol(p), fDet(d)
  {} // Sample constructor

  //-------------------------------------------------------------------------
  Sample::Sample(unsigned int id)
  {
    size_t offset = nBitsDet+nBitsPol, val = pow(nBitsSel,2)-1;
    fSel = (Selection)((id & (val << offset)) >> offset);
    offset = nBitsDet, val = pow(nBitsPol,2)-1;
    fPol = (Polarity)((id & (val << offset)) >> offset);
    val = pow(nBitsDet,2)-1;
    fDet = (Detector)(id & val);
  } // Sample constructor

  //-------------------------------------------------------------------------
  const HistAxis Sample::GetAxis() const
  {
    if (fAxis.GetLabels().empty()) Error("axis not set in sample "+Name());
    return fAxis;
  } // function Sample::GetAxis

  //-------------------------------------------------------------------------
  const Cut Sample::GetCut() const
  {
    if (fCut.ID() == kNoCut.ID())  Error("cut not set in sample "+Name());
    return fCut;
  } // function Sample::GetCut

  //-------------------------------------------------------------------------
  double Sample::POT() const
  {
    if (HasData()) return fData.POT();
    if (fPOT == -1)  Error("POT not set in sample "+Name());
    return fPOT;
  } // function Sample::POT

  //-------------------------------------------------------------------------
  double Sample::Livetime() const
  {
    if (HasData()) return fData.Livetime();
    if (fLivetime == -1)  Error("livetime not set in sample "+Name());
    return fLivetime;
  } // function Sample::Livetime

  //---------------------------------------------------------------------------
  std::vector<OscChannel> Sample::AllChannels() const
  {
    std::vector<OscChannel> ret;
    size_t nFlavs = IsFD() ? 6 : 2;
    // Add each combination of flavour and sign
    for (size_t i = 0; i < nFlavs; ++i)
      for (size_t j = 0; j < 2; ++j)
        ret.push_back(OscChannel(OscChannel::CCName(i,j)));
    // Add the neutral currents
    ret.push_back(OscChannel(OscChannel::NCName()));
    return ret;
  } // function Sample::AllChannels

  //---------------------------------------------------------------------------
  bool Sample::IsSignal(const OscChannel& c) const
  {
    if (IsNC() && c.Curr() == Current::kNC)
      return true;
    if (IsNumu() && c.Curr() == Current::kCC && c.Flav() == Flavors::kNuMuToNuMu)
      return true;
    if (IsNue() && c.Curr() == Current::kCC && c.Flav() == Flavors::kNuMuToNuE)
      return true;
    return false;
  } // function Sample::IsSignal

  //---------------------------------------------------------------------------
  std::vector<OscChannel> Sample::SignalChannels() const
  {
    std::vector<OscChannel> ret;
    for (const OscChannel& c : AllChannels())
      if (IsSignal(c))
        ret.push_back(c);
    return ret;
  } // function Sample::SignalChannels

  //---------------------------------------------------------------------------
  std::vector<OscChannel> Sample::BackgroundChannels() const
  {
    std::vector<OscChannel> ret;
    for (const OscChannel& c : AllChannels())
      if (!IsSignal(c))
        ret.push_back(c);
    return ret;
  } // function Sample::BackgroundChannels

  //---------------------------------------------------------------------------
  Spectrum Sample::Predict(osc::IOscCalc* calc,
                           const SystShifts& shifts) const
  {
    return Prediction()->PredictSyst(calc, Shifts(shifts));
  } // function Sample::Predict

  //---------------------------------------------------------------------------
  Spectrum Sample::PredictComponent(Flavors::Flavors_t flav,
                                    Current::Current_t curr,
                                    Sign::Sign_t sign,
                                    osc::IOscCalc* calc,
                                    const SystShifts& shifts) const
  {
    return Prediction()->PredictComponentSyst(calc, Shifts(shifts), flav, curr, sign);
  } // function Sample::PredictComponent

  //---------------------------------------------------------------------------
  Spectrum Sample::PredictChannel(const OscChannel& channel,
                                  osc::IOscCalc* calc,
                                  const SystShifts& shifts) const
  {
    return PredictComponent(channel.Flav(), channel.Curr(), channel.Sign(), calc, Shifts(shifts));
  } // function Sample::PredictChannel

  //---------------------------------------------------------------------------
  Spectrum Sample::PredictSignal(osc::IOscCalc* calc,
                                 const SystShifts& shifts) const
  {
    Spectrum ret = Spectrum::Uninitialized();
    for (const OscChannel& c : SignalChannels()) {
      Spectrum spec = PredictChannel(c, calc, shifts);
      if (!ret.NDimensions()) ret = spec;
      else ret += spec;
    } // for signal channel
    return ret;
  } // function Sample::PredictSignal

  //---------------------------------------------------------------------------
  Spectrum Sample::PredictBackground(osc::IOscCalc* calc,
                                     const SystShifts& shifts) const
  {
    Spectrum ret = Spectrum::Uninitialized();
    for (const OscChannel& c : BackgroundChannels()) {
      Spectrum spec = PredictChannel(c, calc, shifts);
      if (!ret.NDimensions()) ret = spec;
      else ret += spec;
    } // for background channel
    return ret;
  } // function Sample::PredictBackground

  //-------------------------------------------------------------------------
  Spectrum Sample::Data() const
  {
    if (!HasData()) Error("data spectrum not set in sample "+Name());
    return fData;
  } // function Sample::Data

  //-------------------------------------------------------------------------
  Spectrum Sample::Cosmic() const
  {
    if (!HasCosmic()) Error("cosmic spectrum not set in sample "+Name());
    return fCosmic;
  } // function Sample::Cosmic

  //---------------------------------------------------------------------------
  Spectrum Sample::NewSpectrum(const Eigen::ArrayXd& arr) const
  {
    Eigen::ArrayXd tmp = Eigen::ArrayXd::Zero(GetBinning().NBins()+2);
    tmp.segment(1, GetBinning().NBins()) = arr;
    Spectrum ret(std::move(tmp), GetAxis(), POT(), Livetime());
    return ret;
  } // function Sample::NewSpectrum

  //-------------------------------------------------------------------------
  SystShifts Sample::Shifts(SystShifts shifts) const
  {
    SystShifts ret;
    for (const ISyst* syst : shifts.ActiveSysts()) {
      if (fSystMap.count(syst)) {
        if (fSystMap.at(syst)) { // If there's a nullptr here, skip it
          ret.SetShift(fSystMap.at(syst), shifts.GetShift(syst));
        }
      } else {
        ret.SetShift(syst, shifts.GetShift(syst));
      }
    } // for syst
    return ret;
  } // function Sample::GetSystShifts

  //-------------------------------------------------------------------------
  std::vector<const ISyst*> Sample::Systs(std::vector<const ISyst*> systs) const
  {
    std::vector<const ISyst*> ret;
    for (const ISyst* syst : systs) {
      if (!fSystMap.count(syst)) ret.push_back(syst);
      else if (fSystMap.at(syst)) ret.push_back(fSystMap.at(syst));
    }
    return ret;
  } // function Sample::Systs

  //-------------------------------------------------------------------------
  unsigned int Sample::GetID() const
  {
    unsigned int id = fSel;
    id <<= nBitsPol;
    id += fPol;
    id <<= nBitsDet;
    id += fDet;
    return id;
  } // function Sample::GetID

  //-------------------------------------------------------------------------
  std::string Sample::EnsembleID(const std::vector<Sample>& samples)
  {
    std::ostringstream oss;
    oss << "id";
    for (const Sample& s : samples)
      oss << "_" << s.GetID();
    return oss.str();
  } // function Sample::EnsembleID

  //-------------------------------------------------------------------------
  std::vector<Sample> Sample::FromEnsembleID(const std::string& id)
  {
    std::vector<Sample> ret;
    size_t start = 3;
    while (true) {
      size_t end = id.find("_", start);
      unsigned int val = std::stoi(id.substr(start, end));
      ret.push_back(Sample(val));
      if (end == std::string::npos) break;
      start = end + 1;
    }
    return ret;
  } // function Sample::FromEnsembleID

  //-------------------------------------------------------------------------
  std::vector<Sample> Sample::All()
  {
    std::vector<Sample> ret;
    for (const auto& [sel, tmp1] : kSelNames)
      for (const auto& [pol, tmp2] : kPolNames)
        for (const auto& [det, tmp3] : kDetNames)
          ret.push_back(Sample(sel, pol, det));
    return ret;
  } // function Sample::All

  //-------------------------------------------------------------------------
  bool Sample::IsNC() const
  {
    return fSel == kNCOld || fSel == kNCRes10 || fSel == kNCRes20 || fSel == kNCRes30;
  } // function Sample::IsNC

  //-------------------------------------------------------------------------
  bool Sample::IsNumu() const
  {
    return fSel == kCCNumu || fSel == kCCNumuQ1 || fSel == kCCNumuQ2 || fSel == kCCNumuQ3 || fSel == kCCNumuQ4;
  } // function Sample::IsNumu

  //-------------------------------------------------------------------------
  bool Sample::IsNue() const
  {
    return fSel == kCCNue;
  } // function Sample::IsNue

  //-------------------------------------------------------------------------
  bool Sample::IsFHC() const
  {
    return fPol == kFHC;
  } // function Sample::IsFHC

  //-------------------------------------------------------------------------
  bool Sample::IsRHC() const
  {
    return fPol == kRHC;
  } // function Sample::IsRHC

  //-------------------------------------------------------------------------
  bool Sample::IsND() const
  {
    return fDet == kNearDet;
  } // function Sample::IsND

  //-------------------------------------------------------------------------
  bool Sample::IsFD() const
  {
    return fDet == kFarDet;
  } // function Sample::IsFD

  //-------------------------------------------------------------------------
  std::shared_ptr<IPrediction> Sample::Prediction() const
  {
    if (!HasPrediction()) Error("prediction not set in sample "+Name());
    return fPred;
  } // function Sample::Prediction

} // namespace pisces
