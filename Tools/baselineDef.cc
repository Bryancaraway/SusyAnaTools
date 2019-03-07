#include "baselineDef.h"

#include "TFile.h"
#include "TF1.h"

#include "lester_mt2_bisect.h"

//**************************************************************************//
//                              BaselineVessel                              //
//**************************************************************************//

BaselineVessel::BaselineVessel(NTupleReader &tr, const std::string specialization, const std::string filterString) : 
  tr(&tr), spec(specialization), ttPtr(NULL), WMassCorFile(NULL)
{
  debug                 = false;
  printConfig           = false;
  incZEROtop            = false;
  UseLeptonCleanJet     = false;
  UseDRLeptonCleanJet   = false;
  UseDRPhotonCleanJet   = false;
  UseDeepTagger         = true;
  UseDeepCSV            = true;
  eraLabel              = "2016MC";
  jetVecLabel           = "JetTLV";
  CSVVecLabel           = "Jet_btagDeepB";
  METLabel              = "MET_pt";
  METPhiLabel           = "MET_phi";
  jetVecLabelAK8        = "FatJetTLV";
  qgLikehoodLabel       = "qgLikelihood";
  muonsFlagIDLabel      = "Muon_Stop0l"; 
  elesFlagIDLabel       = "Electron_Stop0l";
  toptaggerCfgFile      = "TopTagger.cfg";
  doLeptonVeto          = true;
  doEleVeto             = true;
  doMuonVeto            = true;
  doIsoTrkVeto          = true;
  doMET                 = true;
  dodPhis               = true;
  passBaselineLowDM     = false;
  passBaselineHighDM    = false;
  metLVec.SetPtEtaPhiM(0, 0, 0, 0);
  if (UseDeepCSV)
    CSVVecLabel           = "Jet_btagDeepB";
  if (UseDeepTagger)
    toptaggerCfgFile      = "TopTagger.cfg";
    //toptaggerCfgFile      = "TopTagger_DeepCombined.cfg";

  if(filterString.compare("fastsim") ==0) isfastsim = true; else isfastsim = false; 

  //Check if simplified tagger is called for
  std::string taggerLabel = "";
  const std::string aggBinLabel = "AggregatedBins";
  size_t loc = spec.find(aggBinLabel);
  if(loc != std::string::npos)
  {
    toptaggerCfgFile = "TopTagger_Simplified.cfg";
    taggerLabel = "AggBins";
    //Remove aggBinLabel from spec
    spec.erase(loc, aggBinLabel.size());
    //Strip any white space ledt in spec
    spec.erase(spec.begin(), std::find_if(spec.begin(), spec.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    spec.erase(std::find_if(spec.rbegin(), spec.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), spec.end());
  }

  if( !spec.empty() ){
    TString stripT = spec;
    TObjArray * objArr = stripT.Tokenize(" ");
    TObjString* firstObj = dynamic_cast<TObjString*>(objArr->At(0));
    firstSpec = firstObj->GetString().Data();
  }
  firstSpec += taggerLabel;

  //TODO: not updated yet, plan to remove
  printOnce = false;
  PredefineSpec();

  //SetupTopTagger(toptaggerCfgFile);
}

// constructor without nullptr as argument
BaselineVessel::BaselineVessel(const std::string specialization, const std::string filterString) : BaselineVessel(*static_cast<NTupleReader*>(nullptr), specialization, filterString) {}

// ===  FUNCTION  ============================================================
//         Name:  BaselineVessel::UseCleanedJets
//  Description:  By default no Lep clean in Jets. Call this function to
//  switch input labels
// ===========================================================================
bool BaselineVessel::UseCleanedJets() 
{
  std::string prefix = "";
  std::string suffix = "";
  if      (UseLeptonCleanJet)   prefix = "prodJetsNoLep_";
  if      (UseDRPhotonCleanJet) suffix = "_drPhotonCleaned";
  else if (UseDRLeptonCleanJet) suffix = "_drLeptonCleaned";
  jetVecLabel     = prefix + "JetTLV"        + suffix;
  CSVVecLabel     = prefix + "Jet_btagDeepB" + suffix;
  qgLikehoodLabel = prefix + "qgLikelihood"  + suffix;
  jetVecLabelAK8  = prefix + "FatJetTLV"     + suffix;
  if (UseDeepCSV)
  {
    // Note that DeepCSVcomb is a derived variable... but it is derived with cleaned variables 
    CSVVecLabel   = prefix + "Jet_btagDeepB" + suffix;
  }
  return true;
}       // -----  end of function BaselineVessel::UseCleanedJets  -----

bool BaselineVessel::OpenWMassCorrFile()
{
  std::string puppiCorr = "puppiCorr.root";
  WMassCorFile = TFile::Open(puppiCorr.c_str(),"READ");
  if (!WMassCorFile)
    std::cout << "W mass correction file not found w mass!!!!!!! " << puppiCorr <<" Will not correct W mass" << std::endl;
  else{
    puppisd_corrGEN      = (TF1*)WMassCorFile->Get("puppiJECcorr_gen");
    puppisd_corrRECO_cen = (TF1*)WMassCorFile->Get("puppiJECcorr_reco_0eta1v3");
    puppisd_corrRECO_for = (TF1*)WMassCorFile->Get("puppiJECcorr_reco_1v3eta2v5");
  }
  return true;
}       // -----  end of function BaselineVessel::OpenWMassCorrFile  -----

std::string BaselineVessel::UseCleanedJetsVar(std::string varname) const
{
  std::string prefix = "";
  std::string suffix = "";
  if      (UseLeptonCleanJet)   prefix = "prodJetsNoLep_";
  if      (UseDRPhotonCleanJet) suffix = "_drPhotonCleaned";
  else if (UseDRLeptonCleanJet) suffix = "_drLeptonCleaned";
  return prefix + varname + suffix;
}       // -----  end of function BaselineVessel::UseCleanedJetsVar  -----

bool BaselineVessel::SetupTopTagger(std::string CfgFile_)
{
  toptaggerCfgFile = CfgFile_;

  ttPtr.reset(new TopTagger);
  ttPtr->setCfgFile(toptaggerCfgFile);
  OpenWMassCorrFile();
  
  return true;
}       // -----  end of function BaselineVessel::SetupTopTagger  -----

void BaselineVessel::prepareDeepTopTagger()
{
  *jetsLVec_forTagger     = tr->getVec<TLorentzVector>(jetVecLabel);
  *recoJetsBtag_forTagger = tr->getVec<float>(CSVVecLabel);
  *qgLikelihood_forTagger = tr->getVec<float>(qgLikehoodLabel);

  // ----- AK4 Jets -----
  
  // AK4 jet variables
  std::map<std::string, std::string> AK4Variables;
  //            name for top tagger                     name in ntuple
  AK4Variables["qgPtD"]                               = "qgPtD"; 
  AK4Variables["qgAxis1"]                             = "qgAxis1"; 
  AK4Variables["qgAxis2"]                             = "qgAxis2"; 
  AK4Variables["recoJetschargedHadronEnergyFraction"] = "recoJetschargedHadronEnergyFraction"; 
  AK4Variables["recoJetschargedEmEnergyFraction"]     = "recoJetschargedEmEnergyFraction"; 
  AK4Variables["recoJetsneutralEmEnergyFraction"]     = "recoJetsneutralEmEnergyFraction"; 
  AK4Variables["recoJetsmuonEnergyFraction"]          = "recoJetsmuonEnergyFraction"; 
  AK4Variables["recoJetsHFHadronEnergyFraction"]      = "recoJetsHFHadronEnergyFraction";
  AK4Variables["recoJetsHFEMEnergyFraction"]          = "recoJetsHFEMEnergyFraction"; 
  AK4Variables["recoJetsneutralEnergyFraction"]       = "recoJetsneutralEnergyFraction"; 
  AK4Variables["PhotonEnergyFraction"]                = "PhotonEnergyFraction"; 
  AK4Variables["ElectronEnergyFraction"]              = "ElectronEnergyFraction";
  AK4Variables["ChargedHadronMultiplicity"]           = "ChargedHadronMultiplicity";
  AK4Variables["NeutralHadronMultiplicity"]           = "NeutralHadronMultiplicity";
  AK4Variables["PhotonMultiplicity"]                  = "PhotonMultiplicity"; 
  AK4Variables["ElectronMultiplicity"]                = "ElectronMultiplicity"; 
  AK4Variables["MuonMultiplicity"]                    = "MuonMultiplicity"; 
  AK4Variables["DeepCSVb"]                            = "DeepCSVb"; 
  AK4Variables["DeepCSVc"]                            = "DeepCSVc"; 
  AK4Variables["DeepCSVl"]                            = "DeepCSVl"; 
  AK4Variables["DeepCSVbb"]                           = "DeepCSVbb";     
  AK4Variables["DeepCSVcc"]                           = "DeepCSVcc";      
  
  ttUtility::ConstAK4Inputs<float> AK4Inputs(*jetsLVec_forTagger, *recoJetsBtag_forTagger);
  
  // convert qgMult to float and then add to AK4Inputs
  const std::vector<int>   &qgMult_i = tr->getVec<int>(UseCleanedJetsVar("qgMult"));
  const std::vector<float> qgMult_f(qgMult_i.begin(), qgMult_i.end());
  AK4Inputs.addSupplamentalVector("qgMult", qgMult_f);
  
  // loop over variables and add to AK4Inputs
  for (const auto& variable : AK4Variables)
  {
    //                              first: name for top tagger                          second: name in ntuple
    AK4Inputs.addSupplamentalVector(variable.first, tr->getVec<float>(UseCleanedJetsVar(variable.second)));
  }

  // ----- AK8 Jets -----

  const std::vector<TLorentzVector>              &AK8JetLV           = tr->getVec<TLorentzVector>(jetVecLabelAK8);
  const std::vector<float>                       &AK8JetSoftdropMass = tr->getVec<float>(UseCleanedJetsVar("puppisoftDropMass"));
  const std::vector<float>                       &AK8JetDeepAK8Top   = tr->getVec<float>(UseCleanedJetsVar("deepAK8btop"));
  const std::vector<float>                       &AK8JetDeepAK8W     = tr->getVec<float>(UseCleanedJetsVar("deepAK8bW"));
  const std::vector<std::vector<TLorentzVector>> &AK8SubjetLV        = tr->getVec<std::vector<TLorentzVector>>(UseCleanedJetsVar("puppiAK8SubjetLVec"));

  //Create AK8 inputs object
  ttUtility::ConstAK8Inputs<float> AK8Inputs(
      AK8JetLV,
      AK8JetDeepAK8Top,
      AK8JetDeepAK8W,
      AK8JetSoftdropMass,
      AK8SubjetLV
      );

  //Create jets constituents list combining AK4 and AK8 jets, these are used to construct top candiates
  //The vector of input constituents can also be constructed "by hand"
  std::vector<Constituent> constituents = ttUtility::packageConstituents(AK4Inputs, AK8Inputs);

  //run the top tagger
  ttPtr->runTagger(constituents);

}


BaselineVessel::~BaselineVessel() 
{
}

bool BaselineVessel::PredefineSpec()
{

  if( spec.compare("noIsoTrksVeto") == 0)
  {
    doIsoTrkVeto = false;
  }
  else if( spec.compare("incZEROtop") == 0)
  {
    incZEROtop = true;
  }
  else if( spec.compare("hadtau") == 0)
  {
    doMuonVeto = false;
    doIsoTrkVeto = false;
    METLabel = "met_hadtau";
    METPhiLabel = "metphi_hadtau";
    jetVecLabel = "jetsLVec_hadtau";
    CSVVecLabel = "recoJetsCSVv2_hadtau";
  }
  else if( spec.compare("lostlept") == 0)
  {
    doLeptonVeto = false;
    doEleVeto    = false;
    doMuonVeto   = false;
    doIsoTrkVeto = false;
  }
  else if (spec.compare("NoVeto") == 0)
  {
    METLabel    = "cleanMetPt";
    METPhiLabel = "cleanMetPhi";
    
    UseDeepCSV          = false; // broken in CMSSW8028_2016 ntuples 
    UseLeptonCleanJet   = false;
    UseDRPhotonCleanJet = false;
    UseDRLeptonCleanJet = false;
    doLeptonVeto = false;
    doEleVeto    = false;
    doMuonVeto   = false;
    doIsoTrkVeto = false;
    dodPhis = false;
  }
  else if (spec.compare("PFLeptonCleaned") == 0)
  {
    METLabel    = "cleanMetPt";
    METPhiLabel = "cleanMetPhi";
    
    UseDeepCSV          = false;
    UseLeptonCleanJet   = true;
    UseDRPhotonCleanJet = false;
    UseDRLeptonCleanJet = false;
    doLeptonVeto = false;
    doEleVeto    = false;
    doMuonVeto   = false;
    doIsoTrkVeto = false;
    dodPhis = false;
  }
  // Z invisible Z to LL control region
  else if (spec.compare("DRLeptonCleaned") == 0)
  {
    METLabel    = "cleanMetPt";
    METPhiLabel = "cleanMetPhi";
    
    UseDeepCSV          = true;
    UseLeptonCleanJet   = false;
    UseDRPhotonCleanJet = false;
    UseDRLeptonCleanJet = true;
    doLeptonVeto = false;
    doEleVeto    = false;
    doMuonVeto   = false;
    doIsoTrkVeto = false;
    dodPhis = true;
  }
  // Z invisible photon control region
  else if (spec.compare("DRPhotonCleaned") == 0)
  {
    METLabel    = "metWithPhoton";
    METPhiLabel = "metphiWithPhoton";
    
    UseDeepCSV          = false;
    UseLeptonCleanJet   = false;
    UseDRPhotonCleanJet = true;
    UseDRLeptonCleanJet = false;
    doLeptonVeto = true;
    doEleVeto    = true;
    doMuonVeto   = true;
    doIsoTrkVeto = true;
    dodPhis = true;
  }
  else if(spec.compare("Zinv") == 0 || spec.compare("Zinv1b") == 0 || spec.compare("Zinv2b") == 0 || spec.compare("Zinv3b") == 0 || spec.compare("ZinvJEUUp") == 0 || spec.compare("ZinvJEUDn") == 0 || spec.compare("ZinvMEUUp") == 0 || spec.compare("ZinvMEUDn") == 0) 
  {
    UseDeepCSV          = true;
    UseLeptonCleanJet   = false;
    UseDRPhotonCleanJet = true;
    UseDRLeptonCleanJet = false;
    doLeptonVeto = false;
    doEleVeto    = false;
    doMuonVeto   = false;
    doIsoTrkVeto = false;
    dodPhis             = true;
    
    if(spec.compare("Zinv1b") == 0)
    {
      CSVVecLabel = "cleanJetpt30ArrBTag1fake";
      bToFake = 1;
    }
    else if(spec.compare("Zinv2b") == 0)
    {
      CSVVecLabel = "cleanJetpt30ArrBTag2fake";
      bToFake = 1; // This is not a typo
    }
    else if(spec.compare("Zinv3b") == 0)
    {
      CSVVecLabel = "cleanJetpt30ArrBTag3fake";
      bToFake = 1; // This is not a typo
    }
    else if(spec.compare("ZinvJEUUp") == 0)
    {
      jetVecLabel = "jetLVecUp";
    }
    else if(spec.compare("ZinvJEUDn") == 0)
    {
      jetVecLabel = "jetLVecDn";
    }
    else if(spec.compare("ZinvMEUUp") == 0)
    {
      METLabel    = "metMEUUp";
    }
    else if(spec.compare("ZinvMEUDn") == 0)
    {
      METLabel    = "metMEUDn";
    }
  }
  else if(spec.compare("QCD") == 0)
  {
    doMET = false;
    dodPhis = false;
  }else if( spec.find("jecUp") != std::string::npos || spec.find("jecDn") != std::string::npos || spec.find("metMagUp") != std::string::npos || spec.find("metMagDn") != std::string::npos || spec.find("metPhiUp") != std::string::npos || spec.find("metPhiDn") != std::string::npos ){
    if( spec.find("jecUp") != std::string::npos ){
      jetVecLabel = "jetLVec_jecUp";
      CSVVecLabel = "recoJetsBtag_jecUp";
    }else if(spec.find("jecDn") != std::string::npos ){
      jetVecLabel = "jetLVec_jecDn";
      CSVVecLabel = "recoJetsBtag_jecDn";
    }else if(spec.find("metMagUp") != std::string::npos ){
      METLabel = "met_metMagUp";
    }else if(spec.find("metMagDn") != std::string::npos ){
      METLabel = "met_metMagDn";
    }else if(spec.find("metPhiUp") != std::string::npos ){
      METPhiLabel = "metphi_metPhiUp";
    }else if(spec.find("metPhiDn") != std::string::npos ){
      METPhiLabel = "metphi_metPhiDn";
    }
    if( spec.find("usegenmet") != std::string::npos ){
      METLabel = "genmet";
      METPhiLabel = "genmetphi";
    } 
  }else if( spec.compare("usegenmet") == 0 ){
    METLabel = "genmet";
    METPhiLabel = "genmetphi";
  }

  if( !printOnce ){
    printOnce = true;
  }  
  

  return true;
}

bool BaselineVessel::PassTopTagger()
{
  int nTopCandSortedCnt = -1;
  int nWs = -1;
  int nResolvedTops = -1;
  bool passTagger = false;
  vTops = new std::vector<TLorentzVector>();
  vWs = new std::vector<TLorentzVector>();
  vResolvedTops = new std::vector<TLorentzVector>();
  mTopJets = new std::map<int, std::vector<TLorentzVector> >();

  nTopCandSortedCnt = GetnTops();
  nWs = vWs->size();
  nResolvedTops = vResolvedTops->size();
  passTagger = (incZEROtop || nTopCandSortedCnt >= AnaConsts::low_nTopCandSortedSel); 

  tr->registerDerivedVar("nTopCandSortedCnt" + firstSpec, nTopCandSortedCnt);
  tr->registerDerivedVar("nWs"+firstSpec, nWs);
  tr->registerDerivedVar("nResolvedTops"+firstSpec, nResolvedTops);
  tr->registerDerivedVec("vTops"+firstSpec, vTops);
  tr->registerDerivedVec("vWs"+firstSpec, vWs);
  tr->registerDerivedVec("mTopJets"+firstSpec, mTopJets);

  return passTagger;
}       // -----  end of function BaselineVessel::PassTopTagger  -----


bool BaselineVessel::PrintoutConfig() const
{
  if (!tr->isFirstEvent()) return false;
  
  std::cout << "=== Current Config ===" << std::endl;
  std::cout << "    Specialization : " << spec             << std::endl;
  std::cout << "    Era Label      : " << eraLabel         << std::endl;
  std::cout << "    AK4Jet Label   : " << jetVecLabel      << std::endl;
  std::cout << "    b-tag Label    : " << CSVVecLabel      << std::endl;
  std::cout << "    top-tag config : " << toptaggerCfgFile << std::endl;
  std::cout << "    MET Label      : " << METLabel         << std::endl;
  std::cout << "======================" << std::endl;
  return true;
}       // -----  end of function BaselineVessel::PrintoutConfig  -----

void BaselineVessel::PassBaseline()
{
  if (printConfig) PrintoutConfig();
  
  // Get jet collection
  const auto& jet_vec = tr->getVec_LVFromNano<float>(jetVecLabel);

  // Create TLorentzVector for MET
  metLVec.SetPtEtaPhiM(tr->getVar<float>(METLabel), 0, tr->getVar<float>(METPhiLabel), 0);
  float met = metLVec.Pt();
  float metphi = metLVec.Phi();

  // Pass_LeptonVeto
  int Pass_LeptonVeto = tr->getVar<bool>("Pass_LeptonVeto");
  
  // TODO: this is wrong... you need to count number of leptons passing veto selection (Electron_Stop0l, Muon_Stop0l, and IsoTrack_Stop0l bool flags)
  int nElectrons = 0;
  for(const auto& elecPass : tr->getVec<unsigned char>("Electron_Stop0l")) ++nElectrons;
  int nMuons     = 0;
  for(const auto& muonPass : tr->getVec<unsigned char>("Muon_Stop0l")) ++nMuons;
  int nIsoTrks   = 0;
  for(const auto& isoTrkPass : tr->getVec<unsigned char>("IsoTrack_Stop0l")) ++nIsoTrks;

  // Calculate number of jets and b-tagged jets
  int cntCSVS = 0;
  vBidxs = new std::vector<unsigned int>();
  vBjs = new std::vector<TLorentzVector>();
  // TODO: Move the cut value to map
  // 2016 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
  // 2017 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X
  
  // use Stop0l_nbtags
  cntCSVS = tr->getVar<int>("Stop0l_nbtags");  
    
  // Getting the b-jets. Sorted by pt by default
  for(auto idx : *vBidxs)
    vBjs->push_back(jet_vec.at(idx));

  // Pass lepton veto?
  bool passMuonVeto = (nMuons == AnaConsts::nMuonsSel), passEleVeto = (nElectrons == AnaConsts::nElectronsSel), passIsoTrkVeto = (nIsoTrks == AnaConsts::nIsoTrksSel);
  //bool passIsoLepTrkVeto = (nIsoLepTrks == AnaConsts::nIsoTrksSel), passIsoPionTrkVeto = (nIsoPionTrks == AnaConsts::nIsoTrksSel);
  //bool passLeptVeto = passMuonVeto && passEleVeto && passIsoTrkVeto;
  bool passLeptVeto = Pass_LeptonVeto;
  
  // Pass the baseline MET requirement?
  bool passMET = (met >= AnaConsts::defaultMETcut);

  // Pass the HT cut for trigger?
  float HT = tr->getVar<float>("Stop0l_HT");
  bool passHT = (HT >= AnaConsts::defaultHTcut);

  //const auto& nTops = tr->getVar<unsigned int>("nResolvedTopCandidate");
  const auto& nMergedTops   = tr->getVar<int>("Stop0l_nTop");
  const auto& nResolvedTops = tr->getVar<int>("Stop0l_nResolved");
  bool passTagger = (incZEROtop || nMergedTops >= AnaConsts::low_nTopCandSortedSel); 

  bool passNoiseEventFilter = passNoiseEventFilterFunc();
  bool passQCDHighMETFilter = true;
  bool passFastsimEventFilter = true;

  // low dm and high dm baselines from Hui Wang, branch hui_new_tagger
  // https://github.com/susy2015/SusyAnaTools/blob/hui_new_tagger/Tools/tupleRead.C#L629-L639
  // https://github.com/susy2015/SusyAnaTools/blob/5e4f54e1aa985daff90f1ad7a220b8d17e4b7290/Tools/tupleRead.C#L629-L639
  
  // variables for passBaselineLowDM and passBaselineHighDM
  // get ISR jet
  TLorentzVector ISRJet;
  const auto& FatJets = tr->getVec_LVFromNano<float>(jetVecLabelAK8);
  const auto& Stop0l_ISRJetIdx = tr->getVar<int>("Stop0l_ISRJetIdx");
  if (Stop0l_ISRJetIdx < FatJets.size()) ISRJet = FatJets[Stop0l_ISRJetIdx];

  int cntNJetsPt20Eta24 = AnaFunctions::countJets(jet_vec, AnaConsts::pt20Eta24Arr);

  const auto& ISRpt        = tr->getVar<float>("Stop0l_ISRJetPt");
  const auto& mtb          = tr->getVar<float>("Stop0l_Mtb");
  const auto& ptb          = tr->getVar<float>("Stop0l_Ptb");
  const auto& nWs          = tr->getVar<int>("Stop0l_nW");
  const auto& nBottoms     = cntCSVS;
  const auto& nSoftBottoms = tr->getVar<int>("Stop0l_nSoftb");;
  const auto& nJets        = cntNJetsPt20Eta24;
  float S_met              = met / sqrt(HT);

  //SUS-16-049, low dm, ISR cut
  bool pass_ISR = (
                   ISRpt > 200
                   && fabs(ISRJet.Eta()) < 2.4
                   && fabs(ROOT::Math::VectorUtil::DeltaPhi(ISRJet, metLVec)) > 2
                  );
  
  //SUS-16-049, low dm, mtb cut
  bool pass_mtb_lowdm = (nBottoms == 0 || (nBottoms > 0 && mtb < 175));  

  //SUS-16-049, low dm, dphi(met, j1) > 0.5, dphi(met, j23) > 0.15
  bool passdphi_lowdm = ( 
                          (nJets == 2 && fabs(ROOT::Math::VectorUtil::DeltaPhi(jet_vec[0], metLVec)) > 0.5 && fabs(ROOT::Math::VectorUtil::DeltaPhi(jet_vec[1], metLVec)) > 0.15) ||
                          (nJets  > 2 && fabs(ROOT::Math::VectorUtil::DeltaPhi(jet_vec[0], metLVec)) > 0.5 && fabs(ROOT::Math::VectorUtil::DeltaPhi(jet_vec[1], metLVec)) > 0.15 && fabs(ROOT::Math::VectorUtil::DeltaPhi(jet_vec[2], metLVec)) > 0.15)
                        );
  //SUS-16-049, high dm, dphi(met, jet1234) > 0.5
  bool passdphi_highdm = (
                           nJets >= 4 
                           && fabs(ROOT::Math::VectorUtil::DeltaPhi(jet_vec[0], metLVec)) > 0.5 
                           && fabs(ROOT::Math::VectorUtil::DeltaPhi(jet_vec[1], metLVec)) > 0.5 
                           && fabs(ROOT::Math::VectorUtil::DeltaPhi(jet_vec[2], metLVec)) > 0.5 
                           && fabs(ROOT::Math::VectorUtil::DeltaPhi(jet_vec[3], metLVec)) > 0.5
                         );
  
  //baseline for SUS-16-049 low dm
  passBaselineLowDM = (
                           nMergedTops == 0
                        && nWs == 0
                        && pass_ISR
                        && S_met > 10
                        && passdphi_lowdm
                        && pass_mtb_lowdm
                        && passMET
                        && passNoiseEventFilter
                        && nJets >= 2
                      );      
  
  //baseline for SUS-16-049 high dm plus HT cut
  passBaselineHighDM = (
                            passMET
                         && passNoiseEventFilter
                         && nJets >= 5
                         && passdphi_highdm
                         && nBottoms >= 1
                         && passHT
                       );      
  
  if (doLeptonVeto)
  {
      passBaselineLowDM  = passBaselineLowDM  && passLeptVeto;
      passBaselineHighDM = passBaselineHighDM && passLeptVeto;
  }


  // Register all the calculated variables
  tr->registerDerivedVec("vBjs" + firstSpec, vBjs);
  tr->registerDerivedVar("ISRJet" + firstSpec, ISRJet);
  tr->registerDerivedVar("nSoftBottoms" + firstSpec, nSoftBottoms);
  tr->registerDerivedVar("nMergedTops" + firstSpec, nMergedTops);
  tr->registerDerivedVar("nResolvedTops" + firstSpec, nResolvedTops);
  tr->registerDerivedVar("nBottoms" + firstSpec, nBottoms);
  tr->registerDerivedVar("nWs" + firstSpec, nWs);
  tr->registerDerivedVar("nJets" + firstSpec, nJets);
  tr->registerDerivedVar("ptb" + firstSpec, ptb);
  tr->registerDerivedVar("mtb" + firstSpec, mtb);
  tr->registerDerivedVar("passLeptVeto" + firstSpec, passLeptVeto);
  tr->registerDerivedVar("passMuonVeto" + firstSpec, passMuonVeto);
  tr->registerDerivedVar("passEleVeto" + firstSpec, passEleVeto);
  tr->registerDerivedVar("passIsoTrkVeto" + firstSpec, passIsoTrkVeto);
  tr->registerDerivedVar("passMET" + firstSpec, passMET);
  tr->registerDerivedVar("passHT" + firstSpec, passHT);
  tr->registerDerivedVar("passTagger" + firstSpec, passTagger);
  tr->registerDerivedVar("passNoiseEventFilter" + firstSpec, passNoiseEventFilter);
  tr->registerDerivedVar("passQCDHighMETFilter" + firstSpec, passQCDHighMETFilter);
  tr->registerDerivedVar("passFastsimEventFilter" + firstSpec, passFastsimEventFilter);
  tr->registerDerivedVar("HT" + firstSpec, HT);
  tr->registerDerivedVar("passBaselineLowDM"  + firstSpec, passBaselineLowDM);
  tr->registerDerivedVar("passBaselineHighDM" + firstSpec, passBaselineHighDM);
} 


int BaselineVessel::GetnTops() const
{
  //get output of tagger
  const TopTaggerResults& ttr = ttPtr->getResults();
  //Use result for top var
  std::vector<TopObject*> Ntop = ttr.getTops();  
  unsigned int topidx = 0;

  for(unsigned int it=0; it<Ntop.size(); it++)
  {
    TopObject::Type  type = Ntop.at(it)->getType() ;
    if ( type == TopObject::Type::MERGED_TOP
    || type   == TopObject::Type::SEMIMERGEDWB_TOP
    || type   == TopObject::Type::RESOLVED_TOP)
    {
      vTops->push_back(Ntop.at(it)->P());
      std::vector<TLorentzVector> temp;
      for(auto j : Ntop.at(it)->getConstituents())
      {
        temp.push_back(j->P());
      }
      mTopJets->insert(std::make_pair(topidx++, temp));
    }
    if ( type == TopObject::Type::MERGED_W ) vWs->push_back(Ntop.at(it)->P());
    if ( type == TopObject::Type::RESOLVED_TOP ) vResolvedTops->push_back(Ntop.at(it)->P());
  }

  int nTops = vTops->size();
  return nTops;
}       // -----  end of function VarPerEvent::GetnTops  -----

// ===  FUNCTION  ============================================================
//         Name:  BaselineVessel::FlagAK8Jets
//  Description:  Provide a flag for each AK8 jets
//  Frist digit: Notag, top tag, W in top, W alone
//  Second digit: Nobjet, Has Loose b-jet, Has medium bjet
// ===========================================================================
bool BaselineVessel::FlagAK8Jets()
{
  // AK8 + Ak4 for W + jet
  ttUtility::ConstAK8Inputs<float> myConstAK8Inputs(
      tr->getVec<TLorentzVector>(UseCleanedJetsVar("puppiJetsLVec")),
      tr->getVec<float>(UseCleanedJetsVar("puppitau1")),
      tr->getVec<float>(UseCleanedJetsVar("puppitau2")),
      tr->getVec<float>(UseCleanedJetsVar("puppitau3")),
      tr->getVec<float>(UseCleanedJetsVar("puppisoftDropMass")),
      tr->getVec<std::vector<TLorentzVector>>(UseCleanedJetsVar("puppiAK8SubjetLVec")));

  std::vector<Constituent> AK8constituents;
  myConstAK8Inputs.packageConstituents(AK8constituents);

  vAK8Flag = new std::vector<unsigned>();

  for(auto ak8_ : AK8constituents)
  {
    unsigned flag = FlagAK8FromTagger(ak8_);
    if (flag == NoTag)
    {
      flag = FlagAK8FromCSV(ak8_);
    }
    vAK8Flag->push_back(flag);
  }


  GetWAlone();
  GetISRJet();
  return true;
}

bool BaselineVessel::FlagDeepAK8Jets()
{
  const std::vector<TLorentzVector> &ak8s =  tr->getVec<TLorentzVector>(UseCleanedJetsVar("puppiJetsLVec"));
  const std::vector<float> &btops =  tr->getVec<float>(UseCleanedJetsVar("deepAK8btop"));
  const std::vector<float> &bWs =  tr->getVec<float>(UseCleanedJetsVar("deepAK8bW"));
  vAK8Flag = new std::vector<unsigned>(ak8s.size(), NoTag);
  
  std::vector<TLorentzVector> topjets;
  // Cross check with top tagger
  for(auto tops : *mTopJets)
  {
    if (tops.second.size() == 1)
      topjets.push_back(tops.second.front());
  }

  for (unsigned int i = 0; i < ak8s.size(); ++i)
  {
    bool matched = false;
    for(auto t : topjets)
    {
      if (t == ak8s.at(i))
      {
        vAK8Flag->at(i) = TopTag;
        matched = true;
        break;
      }
    }
    if (matched) continue;
    for(auto t : *vWs)
    {
      if (t == ak8s.at(i))
      {
        vAK8Flag->at(i) = WTag;
        matched = true;
        break;
      }
    }
  }

  for (unsigned int i = 0; i < ak8s.size(); ++i)
  {
    if (vAK8Flag->at(i) != NoTag) continue;
    vAK8Flag->at(i) = FlagAK8DeepFromCSV(i);
  }

  GetISRJet();
  return true;
}

AK8Flag BaselineVessel::FlagAK8FromCSV(Constituent &ak8) const
{
  unsigned loosebcnt =0 ;
  unsigned mediumbcnt = 0;
  const std::vector<TLorentzVector> &jets = tr->getVec<TLorentzVector>(jetVecLabel);
  const std::vector<float> &CSV = tr->getVec<float>(CSVVecLabel);

  for(auto sub : ak8.getSubjets())
  {
    for(unsigned int ij=0; ij<jets.size(); ij++)
    {
      if (sub.p().DeltaR(jets.at(ij)) < 0.4)
      {
        if (jets.at(ij).Pt() < 20 || fabs(jets.at(ij).Eta()) > 2.4) continue;
        if (CSV.at(ij) > AnaConsts::cutCSVS ) mediumbcnt ++;
        if (CSV.at(ij) > AnaConsts::cutCSVL ) loosebcnt ++;
      }
    }
    
  }

  if (mediumbcnt > 0 ) return NoTagMediumb;
  if (loosebcnt > 0 ) return NoTagLooseb;
  return NoTagNob;
}

AK8Flag BaselineVessel::FlagAK8DeepFromCSV(unsigned int AK8index) const
{
  unsigned loosebcnt =0 ;
  unsigned mediumbcnt = 0;

  const std::vector<std::vector<TLorentzVector> > &subjets = tr->getVec<std::vector<TLorentzVector> >(UseCleanedJetsVar("puppiAK8SubjetLVec"));
  const std::vector<TLorentzVector> &jets = tr->getVec<TLorentzVector>(jetVecLabel);
  const std::vector<float> &CSV = tr->getVec<float>(CSVVecLabel);

  for(auto sub : subjets.at(AK8index))
  {
    for(unsigned int ij=0; ij<jets.size(); ij++)
    {
      if (sub.DeltaR(jets.at(ij)) < 0.4)
      {
        if (jets.at(ij).Pt() < 20 || fabs(jets.at(ij).Eta()) > 2.4) continue;
        if (UseDeepCSV)
        {
          if (CSV.at(ij) > AnaConsts::DeepCSV.at(eraLabel).at("cutM")) mediumbcnt ++;
          if (CSV.at(ij) > AnaConsts::DeepCSV.at(eraLabel).at("cutL")) loosebcnt ++;
        }
        else{
          if (CSV.at(ij) > AnaConsts::CSVv2.at(eraLabel).at("cutM")) mediumbcnt ++;
          if (CSV.at(ij) > AnaConsts::CSVv2.at(eraLabel).at("cutM")) loosebcnt ++;
        }
      }
    }
    
  }

  if (mediumbcnt > 0 ) return NoTagMediumb;
  if (loosebcnt > 0 ) return NoTagLooseb;
  return NoTagNob;
}

AK8Flag BaselineVessel::FlagAK8FromTagger(Constituent &ak8 )
{

  for(auto t : *mTopJets)
  {
    // Find the mono top jet from tagger
    if (t.second.size() == 1)
    {
      if (t.second.at(0).Pt() == ak8.P().Pt()
          && t.second.at(0).Eta() == ak8.P().Eta()
          && t.second.at(0).Phi() == ak8.P().Phi())
      {
        return TopTag;
      }
    }
    // Find the W jet from tagger
    if (t.second.size() == 2)
    {
      for(auto w : t.second)
      {
        if (w.Pt() == ak8.P().Pt()
            && w.Eta() == ak8.P().Eta()
            && w.Phi() == ak8.P().Phi())
        {
          return WinTopTag;
        }
      }
    }
    // Find the W jet from 3jet tagger
    if (t.second.size() == 3)
    {
      for(auto tri : t.second)
      {
        for(auto sub : ak8.getSubjets())
        {
          if (tri.DeltaR(sub.p())<0.4)
          {
            return WinTopTag;
          }
        }
      }
    }
  }

  // Looking for stand alone W tagger
  float corrSDMass = ak8.getSoftDropMass() * ak8.getWMassCorr();
  float tau21 = ak8.getTau2()/ak8.getTau1();
  if ( corrSDMass > 65 && corrSDMass < 100 &&
      tau21 < 0.6 && ak8.p().Pt() > 200)
  {
    return WAloneTag;
  }

  return NoTag;
}

bool BaselineVessel::GetWAlone() const
{
  const std::vector<TLorentzVector> &AK8 = tr->getVec<TLorentzVector>(UseCleanedJetsVar("puppiJetsLVec"));
  std::vector<TLorentzVector> *WAlone= new std::vector<TLorentzVector>();
  for(unsigned int i=0; i < vAK8Flag->size(); ++i)
  {
    if (vAK8Flag->at(i) == WAloneTag)
    {
      WAlone->push_back(AK8.at(i));
    }
  }
  tr->registerDerivedVec("vWAlone"+spec, WAlone);
  return true;
}

// ===  FUNCTION  ============================================================
//         Name:  BaselineVessel::GetISRJet
//  Description:  Following the ISR defition from SUS-16-049
//  The hardest of large-R jets with pT > 200GeV
//  Failed loose b-tagging, nor top/W tagging
// ===========================================================================
bool BaselineVessel::GetISRJet() const
{
  const std::vector<TLorentzVector> &AK8 = tr->getVec<TLorentzVector>(UseCleanedJetsVar("puppiJetsLVec"));
  std::vector<TLorentzVector> *ISRJet = new std::vector<TLorentzVector>();
  std::map<float, TLorentzVector> ISRJetAll;

  for(unsigned int i=0; i < vAK8Flag->size(); ++i)
  {
    if (vAK8Flag->at(i) == NoTagNob)
    {
      float pt = AK8.at(i).Pt();
      if (pt > 200 && fabs(AK8.at(i).Eta()) < 2.4
          && fabs(TVector2::Phi_mpi_pi(AK8.at(i).Phi()- tr->getVar<float>(METPhiLabel))) > 2)
        ISRJetAll[pt] = AK8.at(i);
    }
  }
  if (!ISRJetAll.empty())
  {
    ISRJet->push_back(ISRJetAll.rbegin()->second);
  }

  tr->registerDerivedVec("vISRJet"+spec, ISRJet);

  return true;
}

bool BaselineVessel::GetTopCombs() const
{
  std::vector<TLorentzVector> *vCombs = new std::vector<TLorentzVector>();
  std::map<int, std::vector<TLorentzVector> > *vCombJets = new std::map<int, std::vector<TLorentzVector> >();

  std::vector<float> *vDeepResCand = new std::vector<float>();
  std::vector<float> *vDeepTopCand = new std::vector<float>();
  std::vector<float> *vDeepWCand = new std::vector<float>();
  //get output of tagger
  //Only MVA combs so far
  const TopTaggerResults& ttr = ttPtr->getResults();
  int i = 0;
  std::vector<TLorentzVector> temp;
  for(auto topr : ttr.getTopCandidates() )
  {
    vCombs->push_back(topr.P());
    temp.clear();
    const std::vector<const Constituent*>& constituents =  topr.getConstituents();
    for(auto cons : constituents)
    {
      //std::cout << " top score? "  << cons->getTopDisc() << std::endl;
      temp.push_back(cons->P());
    }
    if (constituents.size() == 3)
    {
      vDeepResCand->push_back( topr.getDiscriminator() );
    }
    if (constituents.size() == 1)
    {
      vDeepTopCand->push_back(constituents.front()->getTopDisc());
    }
    vCombJets->insert(std::make_pair(i, temp));
    i++;
  }


  tr->registerDerivedVec("vCombs"+spec, vCombs);
  tr->registerDerivedVec("mCombJets"+spec, vCombJets);
  tr->registerDerivedVec("vDeepResCand"+spec, vDeepResCand);
  tr->registerDerivedVec("vDeepTopCand"+spec, vDeepTopCand);

  return true;
}

std::vector<TLorentzVector>  BaselineVessel::GetAK4NoSubjet(Constituent &ak8, std::vector<TLorentzVector> &ak4jets) const
{
  std::vector<TLorentzVector>  temp;
  for(auto ak4 : ak4jets)
  {
    bool ismatched = false;
    for(auto sub : ak8.getSubjets())
    {
      if (ak4.DeltaR(sub.p())<0.4)
      {
        ismatched = true;
        break;
      }
    }

    if (!ismatched)
    {
      temp.push_back(ak4);
    }
  }

  return temp;
}       // -----  end of function BaselineVessel::GetAK4NoSubjet  -----

bool BaselineVessel::passNoiseEventFilterFunc()
{
  // According to https://twiki.cern.ch/twiki/bin/view/CMS/SUSRecommendationsICHEP16#Filters_to_be_applied,
  // "Do not apply filters to signal monte carlo (fastsim)"
  if( isfastsim ) return true;

  try
  {
    bool cached_rethrow = tr->getReThrow();
    tr->setReThrow(false);

    bool passDataSpec = true;
    if( tr->getVar<unsigned int>("run") >= 100000 ){ // hack to know if it's data or MC...
      int goodVerticesFilter = true;
      if (tr->checkBranch("goodVerticesFilter"))
      {
          goodVerticesFilter = tr->getVar<int>("goodVerticesFilter");
      }
      // new filters
      //const int & globalTightHalo2016Filter = tr->getVar<int>("globalTightHalo2016Filter");
      //bool passGlobalTightHalo2016Filter = (&globalTightHalo2016Filter) != nullptr? tr->getVar<int>("globalTightHalo2016Filter") !=0 : true;
      bool globalTightHalo2016Filter = true;
      if (tr->checkBranch("globalTightHalo2016Filter"))
      {
          globalTightHalo2016Filter = tr->getVar<bool>("globalTightHalo2016Filter");
      }
      bool passGlobalTightHalo2016Filter = globalTightHalo2016Filter;

      //int eeBadScFilter = tr->getVar<int>("eeBadScFilter");
      int eeBadScFilter = true;

      passDataSpec = goodVerticesFilter && eeBadScFilter && passGlobalTightHalo2016Filter;
    }

    bool hbheNoiseFilter = isfastsim? 1:tr->getVar<bool>("Flag_HBHENoiseFilter");
    bool hbheIsoNoiseFilter = isfastsim? 1:tr->getVar<bool>("Flag_HBHENoiseIsoFilter");
    bool ecalTPFilter = tr->getVar<bool>("Flag_EcalDeadCellTriggerPrimitiveFilter");

    //bool jetIDFilter = isfastsim? 1:tr->getVar<bool>("AK4looseJetID");
    //bool jetIDFilter = true;
    // new filters

    bool BadPFMuonFilter = true;
    std::string Flag_BadPFMuonFilter_type;
    tr->getType("Flag_BadPFMuonFilter", Flag_BadPFMuonFilter_type);
    if (Flag_BadPFMuonFilter_type.compare("bool") == 0 || Flag_BadPFMuonFilter_type.compare("Bool_t") == 0)
    {
        BadPFMuonFilter = tr->getVar<bool>("Flag_BadPFMuonFilter");
    }
    else if (Flag_BadPFMuonFilter_type.compare("char") == 0 || Flag_BadPFMuonFilter_type.compare("UChar_t") == 0)
    {
        BadPFMuonFilter = bool(tr->getVar<char>("Flag_BadPFMuonFilter"));
    }

    //const auto& BadPFMuonFilter = bool(tr->getVar<char>("Flag_BadPFMuonFilter"));
    //const auto& BadPFMuonFilter = bool(tr->getVar<bool>("Flag_BadPFMuonFilter"));
    bool passBadPFMuonFilter = &BadPFMuonFilter != nullptr ? BadPFMuonFilter : true;

    bool BadChargedCandidateFilter = true;
    std::string Flag_BadChargedCandidateFilter_type;
    tr->getType("Flag_BadChargedCandidateFilter", Flag_BadChargedCandidateFilter_type);
    if (Flag_BadChargedCandidateFilter_type.compare("bool") == 0 || Flag_BadChargedCandidateFilter_type.compare("Bool_t") == 0)
    {
        BadChargedCandidateFilter = tr->getVar<bool>("Flag_BadChargedCandidateFilter");
    }
    else if (Flag_BadChargedCandidateFilter_type.compare("char") == 0 || Flag_BadChargedCandidateFilter_type.compare("UChar_t") == 0)
    {
        BadChargedCandidateFilter = bool(tr->getVar<char>("Flag_BadChargedCandidateFilter"));
    }
    //const auto& BadChargedCandidateFilter = bool(tr->getVar<char>("Flag_BadChargedCandidateFilter"));
    bool passBadChargedCandidateFilter = &BadChargedCandidateFilter != nullptr ? BadChargedCandidateFilter : true;

    //bool passMETratioFilter = tr->getVar<float>("calomet")!=0 ? tr->getVar<float>("met")/tr->getVar<float>("calomet") < 5 : true;

    tr->setReThrow(cached_rethrow);
    //return passDataSpec && hbheNoiseFilter && hbheIsoNoiseFilter && ecalTPFilter && jetIDFilter && passBadPFMuonFilter && passBadChargedCandidateFilter && passMETratioFilter;
    return passDataSpec && hbheNoiseFilter && hbheIsoNoiseFilter && ecalTPFilter && passBadPFMuonFilter && passBadChargedCandidateFilter;
  }
  catch (std::string var)
  {
    if(tr->isFirstEvent()) 
    {
      printf("NTupleReader::getTupleObj(const std::string var):  Variable not found: \"%s\"!!!\n", var.c_str());
      printf("Running with PHYS14 Config\n");
    }
  }
  return true;
}

bool BaselineVessel::passQCDHighMETFilterFunc()
{
  std::vector<TLorentzVector> jetsLVec = tr->getVec<TLorentzVector>(jetVecLabel);
  std::vector<float> recoJetsmuonEnergyFraction = tr->getVec<float>("recoJetsmuonEnergyFraction");
  float metphi = tr->getVar<float>("metphi");

  int nJetsLoop = recoJetsmuonEnergyFraction.size();
  std::vector<float> dPhisVec = AnaFunctions::calcDPhi( jetsLVec, metphi, nJetsLoop, AnaConsts::dphiArr);

  for(int i=0; i<nJetsLoop ; i++)
  {
    float thisrecoJetsmuonenergy = recoJetsmuonEnergyFraction.at(i)*(jetsLVec.at(i)).Pt();
    if( (recoJetsmuonEnergyFraction.at(i)>0.5) && (thisrecoJetsmuonenergy>200) && (std::abs(dPhisVec.at(i)-3.1416)<0.4) ) return false;
  }

  return true;
}

bool BaselineVessel::passFastsimEventFilterFunc()
{
  bool passFilter = true;

  if( isfastsim ){
    bool cached_rethrow = tr->getReThrow();
    tr->setReThrow(false);
    const std::vector<TLorentzVector> & genjetsLVec = tr->getVec<TLorentzVector>("genjetsLVec");
    const std::vector<TLorentzVector> & recoJetsLVec = tr->getVec<TLorentzVector>("jetsLVec");
    const std::vector<float> & recoJetschargedHadronEnergyFraction = tr->getVec<float>("recoJetschargedHadronEnergyFraction");

    if( recoJetschargedHadronEnergyFraction.size() != recoJetsLVec.size() ) std::cout<<"\nWARNING ... Non-equal recoJetschargedHadronEnergyFraction.size : "<<recoJetschargedHadronEnergyFraction.size()<<"  recoJetsLVec.size : "<<recoJetsLVec.size()<<std::endl<<std::endl;

    if( !recoJetsLVec.empty() && (&genjetsLVec) != nullptr ){
      for(unsigned int ij=0; ij<recoJetsLVec.size(); ij++){
        //                if( !AnaFunctions::jetPassCuts(recoJetsLVec[ij], AnaConsts::pt20Eta25Arr) ) continue;
        if( !AnaFunctions::jetPassCuts(recoJetsLVec[ij], AnaConsts::pt30Eta24Arr) ) continue;
        float mindeltaR = 999.0;
        int matchedgenJetsIdx = -1;
        for(unsigned int ig=0; ig<genjetsLVec.size(); ig++){
          float dR = recoJetsLVec[ij].DeltaR(genjetsLVec[ig]);
          if( mindeltaR > dR ){ mindeltaR = dR; matchedgenJetsIdx = (int)ig; }
        }
        if( matchedgenJetsIdx != -1 && mindeltaR > 0.3 && recoJetschargedHadronEnergyFraction[ij] < 0.1 ) passFilter = false;
      }
    }
    tr->setReThrow(cached_rethrow);
  }
  return passFilter;
}

float BaselineVessel::CalcMT2() const
{

  //~~~~ Initial the input ~~~~~
  TLorentzVector fatJet1LVec(0, 0, 0,0);
  TLorentzVector fatJet2LVec(0, 0, 0,0);
  //get output of tagger
  const TopTaggerResults& ttr = ttPtr->getResults();
  //Use result for top var
  const std::vector<TopObject*> &Ntop = ttr.getTops();  

  if (Ntop.size() == 0)
  {
    return 0.0;
  }

  if (Ntop.size() == 1)
  {
    fatJet1LVec = Ntop.at(0)->P();
    fatJet2LVec = ttr.getRsys().P();
     
    return coreMT2calc(fatJet1LVec, fatJet2LVec);
  }

  if (Ntop.size() >= 2)
  {
    std::vector<float> cachedMT2vec;
    for(unsigned int it=0; it<Ntop.size(); it++)
    {
       for(unsigned int jt=it+1; jt<Ntop.size(); jt++)
       {
          cachedMT2vec.push_back(coreMT2calc(Ntop.at(it)->P(), Ntop.at(jt)->P()));
       } 
    }
    std::sort(cachedMT2vec.begin(), cachedMT2vec.end());

    return cachedMT2vec.front();
//    return cachedMT2vec.back();
  }

  return 0.0;
}
void BaselineVessel::operator()(NTupleReader& tr_)
{
  tr = &tr_;
  UseCleanedJets();

  PassBaseline();
  PassTrigger();

  GetSoftbJets();
  GetMHT();
}

void BaselineVessel::PassTrigger()
{
    bool passMuTrigger     = false;
    bool passPhotonTrigger = false;
    
    if( tr->getVar<bool>("HLT_Photon175") ||
        tr->getVar<bool>("HLT_Photon75") ||
        tr->getVar<bool>("HLT_Photon90_CaloIdL_PFHT500") ||
        tr->getVar<bool>("HLT_Photon90")
      )
    {
        passPhotonTrigger = true;
    }

    if( tr->getVar<bool>("HLT_IsoMu24") ||
        tr->getVar<bool>("HLT_IsoTkMu24") ||
        tr->getVar<bool>("HLT_Mu50") ||
        tr->getVar<bool>("HLT_Mu55")
      )
    {
        passMuTrigger = true;
    }
    
    tr->registerDerivedVar("passMuTrigger",     passMuTrigger);
    tr->registerDerivedVar("passPhotonTrigger", passPhotonTrigger);
}


bool BaselineVessel::CombDeepCSV()
{
  std::vector<float> *DeepCSVcomb = new std::vector<float>();
  const std::vector<float> &DeepCSVb = tr->getVec<float>(UseCleanedJetsVar("DeepCSVb"));
  const std::vector<float> &DeepCSVbb = tr->getVec<float>(UseCleanedJetsVar("DeepCSVbb"));
  for (int i = 0; i < DeepCSVb.size(); ++i)
  {
    DeepCSVcomb->push_back(DeepCSVb.at(i) + DeepCSVbb.at(i));
  }

  tr->registerDerivedVec(UseCleanedJetsVar("DeepCSVcomb"), DeepCSVcomb);
  return true;
}

bool BaselineVessel::GetSoftbJets() 
{
  if (!tr->checkBranch("svLVec"))
  {
    return false;
  }

  std::vector<TLorentzVector> *SoftbLVec = new std::vector<TLorentzVector>();


  const std::vector<TLorentzVector> &svLVec   = tr->getVec<TLorentzVector>("svLVec");
  const std::vector<TLorentzVector> &jetsLVec   = tr->getVec<TLorentzVector>("jetsLVec");
  const std::vector<float> &svPT   = tr->getVec<float>("svPT");
  const std::vector<float> &svDXY   = tr->getVec<float>("svDXY");
  const std::vector<float> &svD3D   = tr->getVec<float>("svD3D");
  const std::vector<float> &svD3Derr   = tr->getVec<float>("svD3Derr");
  const std::vector<float> &svNTracks   = tr->getVec<float>("svNTracks");
  const std::vector<float> &svCosThetaSVPS   = tr->getVec<float>("svCosThetaSVPS");

  for(unsigned int i =0; i<svLVec.size(); i++)
  {
    if(svPT.at(i) >= 20.0 ) continue;

    if(svDXY.at(i) >= 3.0) continue;
    if(svD3D.at(i)/svD3Derr.at(i) <= 4.0) continue;
    if(svCosThetaSVPS.at(i) <= 0.98) continue;
    if(svNTracks.at(i)<3) continue;

    bool failDR = false;
    for(unsigned int j=0; j<jetsLVec.size(); j++)
    {
      if (jetsLVec.at(j).Pt() < 20 || fabs(jetsLVec.at(j).Eta()) > 2.4) continue;

      if( ROOT::Math::VectorUtil::DeltaR( svLVec.at(i), jetsLVec.at(j) ) <= 0.4 ) {
        failDR = true;
        break;
      }
    }
    if(failDR) continue;

    SoftbLVec->push_back(svLVec.at(i));
  }

  tr->registerDerivedVec("softbLVec"+firstSpec, SoftbLVec);
  return true;
}       // -----  end of function BaselineVessel::GetSoftbJets  -----

bool BaselineVessel::GetMHT() const
{
  // Calculate MHT
  TLorentzVector MHT(0, 0, 0, 0);
  float SumHT = 0.0; //Using jet > 30 , |eta| < 5
  for(auto &jet : tr->getVec<TLorentzVector>(jetVecLabel))
  {
    if (jet.Pt() >= 30)
    {
      MHT -= jet;
      SumHT += jet.Pt();
    }
  }
  tr->registerDerivedVar("MHT"+firstSpec, static_cast<float>(MHT.Pt()));
  tr->registerDerivedVar("MHTPhi"+firstSpec, static_cast<float>(MHT.Phi()));
  tr->registerDerivedVar("MHTSig"+firstSpec, static_cast<float>(MHT.Pt()/ sqrt(SumHT)));

  tr->registerDerivedVar("METSig"+firstSpec, static_cast<float>(tr->getVar<float>(METLabel)/ sqrt(SumHT)));
  return true;
}

bool BaselineVessel::GetLeptons() const
{
  std::vector<TLorentzVector> *vMuons = new std::vector<TLorentzVector> ();
  std::vector<TLorentzVector> *vEles = new std::vector<TLorentzVector> ();
  std::vector<int> *vMuonChg = new std::vector<int> ();
  std::vector<int> *vEleChg = new std::vector<int> ();

  const std::vector<TLorentzVector> &MuonTLV   = tr->getVec<TLorentzVector>("MuonTLV");
  const std::vector<float>          &muonsRelIso = tr->getVec<float>("muonsMiniIso");
  const std::vector<float>          &muonsMtw    = tr->getVec<float>("muonsMtw");
  const std::vector<int>            &muonsFlagID = tr->getVec<int>(muonsFlagIDLabel.c_str());
  const std::vector<float>          &muonsCharge = tr->getVec<float>("muonsCharge");
  for(unsigned int im=0; im<MuonTLV.size(); im++){
    if(AnaFunctions::passMuon(MuonTLV[im], muonsRelIso[im], muonsMtw[im], muonsFlagID[im], AnaConsts::muonsMiniIsoArr))
    {
      if (!vMuons->empty()) // Making sure the vMuons are sorted in Pt
        assert(MuonTLV.at(im).Pt() <= vMuons->back().Pt());
      vMuons->push_back(MuonTLV.at(im));
      vMuonChg->push_back(muonsCharge.at(im));
    }
  }

  const std::vector<TLorentzVector> &electronsLVec   = tr->getVec<TLorentzVector>("ElectronTLV");
  const std::vector<float> &electronsRelIso         = tr->getVec<float>("elesMiniIso");
  const std::vector<float> &electronsMtw            = tr->getVec<float>("elesMtw");
  const std::vector<int> &electronsFlagID            = tr->getVec<int>(elesFlagIDLabel.c_str());
  const std::vector<float>         &electronsCharge = tr->getVec<float>("elesCharge");
  for(unsigned int ie=0; ie<electronsLVec.size(); ie++){
    if(AnaFunctions::passElectron(electronsLVec[ie], electronsRelIso[ie], electronsMtw[ie], electronsFlagID[ie], AnaConsts::elesMiniIsoArr)) 
    {
      if (!vEles->empty()) // Making sure the vEles are sorted in Pt
        assert(electronsLVec.at(ie).Pt() <= vEles->back().Pt());
      vEles->push_back(electronsLVec.at(ie));
      vEleChg->push_back(electronsCharge.at(ie));

    }
  }

  tr->registerDerivedVar("cutMuID"+firstSpec, muonsFlagIDLabel);
  tr->registerDerivedVar("cutEleID"+firstSpec, elesFlagIDLabel);
  tr->registerDerivedVec("cutMuVec"+firstSpec, vMuons);
  tr->registerDerivedVec("cutEleVec"+firstSpec, vEles);
  tr->registerDerivedVec("cutMuCharge"+firstSpec, vMuonChg);
  tr->registerDerivedVec("cutEleCharge"+firstSpec, vEleChg);

  return true;
}

bool BaselineVessel::GetRecoZ( const int zMassMin, const int zMassMax) const
{
  std::vector<TLorentzVector>* recoZVec = new std::vector<TLorentzVector>();
  std::map<unsigned int, std::pair<unsigned int, unsigned int> > *ZLepIdx = 
    new std::map<unsigned int, std::pair<unsigned int, unsigned int> >();

  GetRecoZ("cutMuVec"+firstSpec,"cutMuCharge"+firstSpec, recoZVec, ZLepIdx, zMassMin, zMassMax );
  GetRecoZ("cutEleVec"+firstSpec,"cutEleCharge"+firstSpec, recoZVec, ZLepIdx, zMassMin, zMassMax );

  tr->registerDerivedVec("recoZVec"+spec, recoZVec);
  tr->registerDerivedVec("ZLepIdx"+spec, ZLepIdx);
  return true;
}       // -----  end of function BaselineVessel::GetRecoZ  -----


bool BaselineVessel::GetRecoZ(const std::string leptype, const std::string lepchg, 
    std::vector<TLorentzVector>* recoZVec ,
    std::map<unsigned int, std::pair<unsigned int, unsigned int> > *ZLepIdx,
    const int zMassMin, const int zMassMax) const
{
  
  const std::vector<TLorentzVector> & LepVec = tr->getVec<TLorentzVector>(leptype);
  const std::vector<int> & LepChg = tr->getVec<int>(lepchg);

  for(unsigned int i = 0; i < LepVec.size(); ++i)
  {
    for(unsigned int j = i; j < LepVec.size(); ++j)
    {
      float zm = (LepVec.at(i) + LepVec.at(j)).M();
      int sumchg = LepChg.at(i) + LepChg.at(j); 
      if (sumchg == 0 && zm > zMassMin && zm < zMassMax)
      {
        recoZVec->push_back((LepVec.at(i) + LepVec.at(j)));
        if (leptype.find("Mu") != std::string::npos)
        {
          ZLepIdx->insert(std::make_pair( recoZVec->size(), 
                std::make_pair(i, j)));
        }
        if (leptype.find("Ele") != std::string::npos) // offset by 100
        {
          ZLepIdx->insert(std::make_pair( recoZVec->size(), 
                std::make_pair(i+100, j+100)));
        }
      }
    }
  }
  return true;
}

bool BaselineVessel::CompCommonVar()
{
  const std::vector<float>  &bdisc = tr->getVec<float>(CSVVecLabel);
  const std::vector<TLorentzVector> &jets = tr->getVec<TLorentzVector>(jetVecLabel);
  std::map<float, unsigned> discmap;
  for(auto idx : *vBidxs)
  {
    discmap[bdisc[idx]] = idx;
  }
  float mtb = 99999;
  float ptb = 0;
  unsigned cnt = 0;

  for (auto iter = discmap.rbegin(); iter != discmap.rend(); ++iter)
  {
    float temp = sqrt(2*jets.at(iter->second).Pt()*tr->getVar<float>(METLabel)
        *(1-cos(jets.at(iter->second).Phi() - tr->getVar<float>(METPhiLabel))));
    mtb = std::min(mtb, temp);
    cnt ++;
    if (cnt == 2) break;
  }
  if (mtb == 99999) mtb=0;

  for (unsigned i = 0; i < vBjs->size() && i < 2; ++i)
  {
    ptb += vBjs->at(i).Pt();
  }

  tr->registerDerivedVar("mtb"+firstSpec, mtb);
  tr->registerDerivedVar("ptb"+firstSpec, ptb);

  return true;
}

