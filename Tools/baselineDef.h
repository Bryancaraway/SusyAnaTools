#ifndef ANABASELINE_BASELINEDEF_H
#define ANABASELINE_BASELINEDEF_H

#include "NTupleReader.h"
#include "customize.h"
#include "EventListFilter.h"

#include "Math/VectorUtil.h"

#include <memory>
#include <iostream>

class TFile;
class TF1;

enum AK8Flag : unsigned
{
  NoTag,
  TopTag,
  WTag,
  WinTopTag,
  WAloneTag,
  NoTagLooseb,
  NoTagMediumb,
  NoTagNob,
};

class BaselineVessel
{
private:
    std::string spec;
    bool isfastsim;
    std::string firstSpec;
    bool printOnce;

    NTupleReader *tr;

//~~~~~~~~~~~~~~~~~~ TopTagger ~~~~~~~~~~~~~~~~~~
    std::shared_ptr<TopTagger> ttPtr;

    //  container
    TLorentzVector metLVec; 
    std::vector<TLorentzVector> *jetsLVec_forTagger;
    std::vector<float> *recoJetsBtag_forTagger;
    std::vector<float> *qgLikelihood_forTagger;
    std::vector<TLorentzVector> *vBjs;
    std::vector<unsigned int> *vBidxs;
    std::vector<TLorentzVector> *vTops;
    std::vector<TLorentzVector> *vWs;
    std::vector<TLorentzVector> *vResolvedTops;
    std::map<int, std::vector<TLorentzVector> > *mTopJets;
    std::vector<unsigned> * vAK8Flag;

    std::vector<TLorentzVector> GetAK4NoSubjet(Constituent &ak8, 
    std::vector<TLorentzVector> &ak4jets) const;

public:

//~~~~~~~~~~~~~~~~~~ Member ~~~~~~~~~~~~~~~~~~
    int  bToFake;
    bool debug;
    bool printConfig;
    bool incZEROtop;
    bool UseLeptonCleanJet;
    bool UseDRLeptonCleanJet;
    bool UseDRPhotonCleanJet;
    bool UseDeepTagger;
    bool UseDeepCSV;

    std::string eraLabel;
    std::string jetVecLabel;
    std::string CSVVecLabel;
    std::string METLabel;
    std::string METPhiLabel;
    std::string jetVecLabelAK8;
    std::string muonsFlagIDLabel;
    std::string elesFlagIDLabel;
    std::string qgLikehoodLabel;
    std::string toptaggerCfgFile;
    
    bool doLeptonVeto; 
    bool doEleVeto;
    bool doMuonVeto;
    bool doIsoTrkVeto;
    bool doMET;
    bool dodPhis;
    bool passBaselineLowDM;
    bool passBaselineHighDM;


    BaselineVessel(NTupleReader &tr, const std::string specialization = "", const std::string filterString = "");
    BaselineVessel(const std::string specialization = "", const std::string filterString = "");
    ~BaselineVessel();

    inline std::string UseCleanedJetsVar(std::string varname) const;
    void PassBaseline();
    void PassTrigger();
    bool PrintoutConfig() const;
    bool CompCommonVar();
    bool passNoiseEventFilterFunc();
    bool passQCDHighMETFilterFunc();
    bool passFastsimEventFilterFunc();
    bool PredefineSpec();
    bool UseCleanedJets();
    bool OpenWMassCorrFile();

    bool FlagDeepAK8Jets();
    void operator()(NTupleReader& tr);
    

//~~~~~~~~~~~~~~~~~~ TopTagger ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bool SetupTopTagger(std::string CfgFile_ = "TopTagger.cfg");
    bool PassTopTagger();
    bool GetMHT() const;
    bool GetLeptons() const;
    bool GetRecoZ( const int zMassMin, const int zMassMax) const;
    bool GetRecoZ(const std::string leptype, const std::string lepchg, std::vector<TLorentzVector>* recoZVec,
        std::map<unsigned int, std::pair<unsigned int, unsigned int> > *TypeZLepIdx, 
        const int zMassMin, const int zMassMax) const;
    void prepareTopTagger();
    void prepareDeepTopTagger();
    std::shared_ptr<TopTagger> GetTopTaggerPtr() const {return ttPtr;};
    int GetnTops() const;
    bool GetTopCombs() const;
    float CalcMT2() const;
    float coreMT2calc(const TLorentzVector & fatJet1LVec, const TLorentzVector & fatJet2LVec) const;
    bool FlagAK8Jets();
    AK8Flag FlagAK8FromCSV(Constituent &ak8) const;
    AK8Flag FlagAK8FromTagger(Constituent &ak8 );
    AK8Flag FlagAK8DeepFromCSV(unsigned int AK8index) const;
    bool CombDeepCSV();
    bool GetISRJet() const;
    bool GetWAlone() const;
    bool GetSoftbJets();
};

inline void passBaselineFunc(NTupleReader &tr, std::string filterstring)
{
  BaselineVessel blv(tr, "", filterstring);
  blv.PassBaseline();
  blv.GetMHT();
}

#endif
