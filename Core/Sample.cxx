#include "CAFAna/PISCES/Sample.h"

#include <cassert>

namespace pisces {

  using std::string;
  using std::vector;
  using namespace ana;

  //-------------------------------------------------------------------------
  Sample::Sample(Selection s, Polarity p, Detector d)
    : fSel(s), fPol(p), fDet(d), fAxis(kEmptyAxis), fCut(kNoCut), fPOT(-1),
      fLivetime(-1), fPred(nullptr), fData(Spectrum::Uninitialized()),
      fCosmic(Spectrum::Uninitialized()), fIsAux(false)
  {} // Sample constructor

  //-------------------------------------------------------------------------
  Sample::Sample(unsigned int id)
    : fAxis(kEmptyAxis), fCut(kNoCut), fPOT(-1), fLivetime(-1),
      fPred(nullptr), fData(Spectrum::Uninitialized()),
      fCosmic(Spectrum::Uninitialized()), fIsAux(false)
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
    if (fAxis.GetLabels().empty())
      assert(false && "Axis not set in sample!");
    return fAxis;
  } // function Sample::GetAxis

  //-------------------------------------------------------------------------
  const Cut Sample::GetCut() const
  {
    if (fCut.ID() == kNoCut.ID()) assert(false && "Cut not set in sample!");
    return fCut;
  } // function Sample::GetCut

  //-------------------------------------------------------------------------
  double Sample::POT() const
  {
    if (HasData()) return fData.POT();
    if (fPOT == -1) assert(false && "POT not set in sample!");
    return fPOT;
  } // function Sample::POT

  //-------------------------------------------------------------------------
  double Sample::Livetime() const
  {
    if (HasData()) return fData.Livetime();
    if (fLivetime == -1) assert(false && "Livetime not set in sample!");
    return fLivetime;
  } // function Sample::Livetime

  //-------------------------------------------------------------------------
  Spectrum Sample::Predict(osc::IOscCalc* calc) const
  {
    assert(HasPrediction() && "Prediction not set in sample!");
    return fPred->Predict(calc);
  } // function Sample::Predict

  //-------------------------------------------------------------------------
  Spectrum Sample::PredictComponent(osc::IOscCalc* calc,
                                    Flavors::Flavors_t flav,
                                    Current::Current_t curr,
                                    Sign::Sign_t sign) const
  {
    assert(HasPrediction() && "Prediction not set in sample!");
    return fPred->PredictComponent(calc, flav, curr, sign);
  } // function Sample::PredictComponent

  //-------------------------------------------------------------------------
  Spectrum Sample::PredictSyst(osc::IOscCalc* calc,
                               const SystShifts& syst) const
  {
    assert(HasPrediction() && "Prediction not set in sample!");
    return fPred->PredictSyst(calc, syst);
  } // function Sample::PredictSyst

  //-------------------------------------------------------------------------
  Spectrum Sample::PredictComponentSyst(osc::IOscCalc* calc,
                                        const SystShifts& syst,
                                        Flavors::Flavors_t flav,
                                        Current::Current_t curr,
                                        Sign::Sign_t sign) const
  {
    assert(HasPrediction() && "Prediction not set in sample!");
    return fPred->PredictComponentSyst(calc, syst, flav, curr, sign);
  } // function Sample::PredictComponentSyst

  //-------------------------------------------------------------------------
  Spectrum Sample::Data() const
  {
    if (!HasData()) assert(false && "Data spectrum not set in sample!");
    return fData;
  } // function Sample::Data

  //-------------------------------------------------------------------------
  Spectrum Sample::Cosmic() const
  {
    if (!HasCosmic()) assert(false && "Cosmic spectrum not set in sample!");
    return fCosmic;
  } // function Sample::Cosmic

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
  vector<const ISyst*> Sample::Systs(vector<const ISyst*> systs) const
  {
    vector<const ISyst*> ret;
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
  string Sample::EnsembleID(const vector<Sample>& samples)
  {
    std::ostringstream oss;
    oss << "id";
    for (const Sample& s : samples)
      oss << "_" << s.GetID();
    return oss.str();
  } // function Sample::EnsembleID

  //-------------------------------------------------------------------------
  vector<Sample> Sample::FromEnsembleID(string const& id)
  {
    vector<Sample> ret;
    size_t start = 3;
    while (true) {
      size_t end = id.find("_", start);
      unsigned int val = stoi(id.substr(start, end));
      ret.push_back(Sample(val));
      if (end == string::npos) break;
      start = end + 1;
    }
    return ret;
  } // function Sample::FromEnsembleID

  //-------------------------------------------------------------------------
  vector<Sample> Sample::All()
  {
    vector<Sample> ret;
    for (auto const& [sel, tmp1] : kSelNames)
      for (auto const& [pol, tmp2] : kPolNames)
        for (auto const& [det, tmp3] : kDetNames)
          ret.push_back(Sample(sel, pol, det));
    return ret;
  } // function Sample::All

  //-------------------------------------------------------------------------
  bool Sample::IsNC() const
  {
    return fSel == kNCOld || fSel == kNCRes10 || fSel == kNCRes20 ||
           fSel == kNCRes30;
  } // function Sample::IsNC

  //-------------------------------------------------------------------------
  bool Sample::IsNumu() const
  {
    return fSel == kCCNumu || fSel == kCCNumuQ1 || fSel == kCCNumuQ2 ||
           fSel == kCCNumuQ3 || fSel == kCCNumuQ4;
  } // return Sample::IsNumu

  //-------------------------------------------------------------------------
  bool Sample::IsNue() const
  {
    return fSel == kCCNue;
  } // return Sample::IsNue

  //-------------------------------------------------------------------------
  bool Sample::IsFHC() const
  {
    return fPol == kFHC;
  } // return Sample::IsFHC

  //-------------------------------------------------------------------------
  bool Sample::IsRHC() const
  {
    return fPol == kRHC;
  } // return Sample::IsRHC

  //-------------------------------------------------------------------------
  bool Sample::IsND() const
  {
    return fDet == kNearDet;
  } // return Sample::IsND

  //-------------------------------------------------------------------------
  bool Sample::IsFD() const
  {
    return fDet == kFarDet;
  } // return Sample::IsFD

} // namespace pisces
