#pragma once

#include "CAFAna/Core/Cut.h"
#include "CAFAna/Core/HistAxis.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SystShifts.h"
#include "CAFAna/Prediction/IPrediction.h"

#include <vector>
#include <string>

namespace pisces {

  class OscChannel;

  using namespace ana;

  void Error(const std::string& msg); // utility function for throwing an error
  void Usage(char** argv, const std::string& args=""); // utility function for printing arguments

  enum Selection
  {
    kCCNumu,
    kCCNumuQ1,
    kCCNumuQ2,
    kCCNumuQ3,
    kCCNumuQ4,
    kCCNue,
    kNCOld,
    kNCRes10,
    kNCRes20,
    kNCRes30
  };
  enum Polarity { kFHC, kRHC };
  enum Detector { kNearDet, kFarDet };

  namespace {
    const HistAxis kEmptyAxis({}, {}, {});

    // define size of each enum so we can bit shift
    const size_t nBitsSel = 6, nBitsPol = 2, nBitsDet = 2;

    const std::map<Selection, std::pair<std::string, std::string>> kSelNames
    {
      { kCCNumu,   { "numusel",    "CC $\\nu_{\\mu}$"    } },
      { kCCNumuQ1, { "numuq1sel",  "CC $\\nu_{\\mu}$ Q1" } },
      { kCCNumuQ2, { "numuq2sel",  "CC $\\nu_{\\mu}$ Q2" } },
      { kCCNumuQ3, { "numuq3sel",  "CC $\\nu_{\\mu}$ Q3" } },
      { kCCNumuQ4, { "numuq4sel",  "CC $\\nu_{\\mu}$ Q4" } },
      { kCCNue,    { "nuesel",     "CC $\\nu_{e}$"       } },
      { kNCOld,    { "ncoldsel",   "NC (old)"            } },
      { kNCRes10,  { "ncres10sel", "NC (10% res)"        } },
      { kNCRes20,  { "ncres20sel", "NC (20% res)"        } },
      { kNCRes30,  { "ncres30sel", "NC (30% res)"        } }
    };

    const std::map<Polarity, std::pair<std::string, std::string>> kPolNames
    {
      { kFHC, { "fhc", "FHC" } },
      { kRHC, { "rhc", "RHC" } }
    };

    const std::map<Detector, std::pair<std::string, std::string>> kDetNames
    {
      { kNearDet, { "neardet", "ND" } },
      { kFarDet, { "fardet", "FD" } }
    };

  } // anonymous namespace

  class Sample {

  public:

    // Constructors
    Sample(Selection s, Polarity p, Detector d);
    Sample(unsigned int id);

    std::string SelStr() const { return kSelNames.at(fSel).first; }
    std::string PolStr() const { return kPolNames.at(fPol).first; }
    std::string DetStr() const { return kDetNames.at(fDet).first; }

    std::string Name() const { return SelStr()+" "+PolStr()+" "+DetStr(); }
    std::string Tag() const { return SelStr()+"_"+PolStr()+"_"+DetStr(); }
    std::string LatexName() const
    {
      std::stringstream oss;
      oss << kSelNames.at(fSel).second << " "
          << kPolNames.at(fPol).second << " "
          << kDetNames.at(fDet).second;
      return oss.str();
    }

    void SetAxis(HistAxis a)   { fAxis     = a; }
    void SetCut(Cut c)         { fCut      = c; }
    void SetPOT(double d)      { fPOT      = d; }
    void SetLivetime(double l) { fLivetime = l; }

    void SetPrediction(std::unique_ptr<IPrediction>& p) { fPred = std::move(p); }
    void SetCosmic(Spectrum s) { fCosmic = s; }
    void SetData(Spectrum s)   { fData   = s; }

    void SetSystAlias(const ISyst* key, const ISyst* val) { fSystMap[key] = val; }
    void SetAuxiliary(bool val) { fIsAux = val; };

    Selection Sel() const { return fSel; }
    Polarity  Pol() const { return fPol; }
    Detector  Det() const { return fDet; }

    const HistAxis GetAxis()    const;
    const Binning  GetBinning() const { return GetAxis().GetBins1D(); }
    const Var      GetVar()     const { return GetAxis().GetVar1D();  }
    const Cut      GetCut()     const;
    double         POT()        const;
    double         Livetime()   const;

    std::vector<OscChannel> AllChannels() const;
    bool IsSignal(const OscChannel& c) const;
    std::vector<OscChannel> SignalChannels() const;
    std::vector<OscChannel> BackgroundChannels() const;

    Spectrum Predict(osc::IOscCalc* calc,
                     const SystShifts& shifts=kNoShift) const;
    Spectrum PredictComponent(Flavors::Flavors_t flav,
                              Current::Current_t curr,
                              Sign::Sign_t sign,
                              osc::IOscCalc* calc,
                              const SystShifts& shifts=kNoShift) const;
    Spectrum PredictChannel(const OscChannel& channel,
                            osc::IOscCalc* calc,
                            const SystShifts& shifts=kNoShift) const;
    Spectrum PredictSignal(osc::IOscCalc* calc,
                           const SystShifts& shifts=kNoShift) const;
    Spectrum PredictBackground(osc::IOscCalc* calc,
                               const SystShifts& shifts=kNoShift) const;

    Spectrum Data()   const;
    Spectrum Cosmic() const;

    Spectrum NewSpectrum(const Eigen::ArrayXd& arr) const;

    SystShifts Shifts(SystShifts shifts) const;
    std::vector<const ISyst*> Systs(std::vector<const ISyst*> systs) const;

    bool IsAuxiliary() const { return fIsAux; };

    bool HasPrediction() const { return bool(fPred);                  }
    bool HasData()       const { return fData.NDimensions();          }
    bool HasCosmic()     const { return fCosmic.NDimensions();        }
    void ResetPrediction()     { fPred.reset();                       }
    void ResetData()           { fData = Spectrum::Uninitialized();   }
    void ResetCosmic()         { fCosmic = Spectrum::Uninitialized(); }

    unsigned int GetID() const;
    static std::string EnsembleID(const std::vector<Sample>& samples);
    static std::vector<Sample> FromEnsembleID(const std::string& id);

    bool operator< (const Sample& lhs) const { return lhs.GetID() < this->GetID(); }
    bool operator==(const Sample& lhs) const { return lhs.GetID() == this->GetID(); }
    bool operator!=(const Sample& lhs) const { return lhs.GetID() != this->GetID(); }

    static std::vector<Sample> All();

    bool IsNC  () const;
    bool IsNumu() const;
    bool IsNue () const;
    bool IsFHC () const;
    bool IsRHC () const;
    bool IsND  () const;
    bool IsFD  () const;

  protected:

    std::shared_ptr<IPrediction> Prediction() const;

    Selection fSel;
    Polarity  fPol;
    Detector  fDet;

    HistAxis fAxis = kEmptyAxis;
    Cut      fCut = kNoCut;
    double   fPOT = -1;
    double   fLivetime = -1;

    std::shared_ptr<IPrediction> fPred;

    Spectrum fData = Spectrum::Uninitialized();
    Spectrum fCosmic = Spectrum::Uninitialized();

    std::map<const ISyst*, const ISyst*> fSystMap;

    bool fIsAux = false;

    friend class Ensemble;

  }; // class Sample

} // namespace pisces
