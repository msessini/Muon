#include <TFile.h>
#include "../interface/MuonMVAEstimator.h"
#include <cmath>
#include <vector>

using namespace std;

#ifndef STANDALONE
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
using namespace reco;
#endif

//--------------------------------------------------------------------------------------------------
MuonMVAEstimator::MuonMVAEstimator() :
  fMethodname("BDTG method"),
  fisInitialized(kFALSE),
  fPrintMVADebug(kFALSE),
  fMVAType(kIDIsoRingsCombined),
  fUseBinnedVersion(kTRUE),
  fNMVABins(0)
{
  // Constructor.  
  fTMVAReader = std::vector<TMVA::Reader*>(0);
}

//--------------------------------------------------------------------------------------------------
MuonMVAEstimator::~MuonMVAEstimator(){
  for (uint i=0;i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void MuonMVAEstimator::initialize( std::string methodName,
				   std::string weightsfile,
				   MuonMVAEstimator::MVAType type){
  
  std::vector<std::string> tempWeightFileVector;
  tempWeightFileVector.push_back(weightsfile);
  initialize(methodName,type,kFALSE,tempWeightFileVector);
}


//--------------------------------------------------------------------------------------------------
void MuonMVAEstimator::initialize( std::string methodName,
				   MuonMVAEstimator::MVAType type,
				   Bool_t useBinnedVersion,
				   std::vector<std::string> weightsfiles
				   ) {
  
  //clean up first
  for (uint i=0;i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
  fTMVAReader.clear();

  //initialize
  fisInitialized = kTRUE;
  fMVAType = type;
  fMethodname = methodName;
  fUseBinnedVersion = useBinnedVersion;

  //Define expected number of bins
  UInt_t ExpectedNBins = 0;
  if (!fUseBinnedVersion) {
    ExpectedNBins = 1;
  } else if (type == kIDIsoRingsCombined) {
    ExpectedNBins = 5;
  } else if (type == kIsoRings)  {
    ExpectedNBins = 4;
  } else if (type == kIsoDeltaR) {
    ExpectedNBins = 4;
  }
  fNMVABins = ExpectedNBins;
  
  //Check number of weight files given
  if (fNMVABins != weightsfiles.size() ) {
    std::cout << "Error: Expected Number of bins = " << fNMVABins << " does not equal to weightsfiles.size() = " 
              << weightsfiles.size() << std::endl;
    assert(fNMVABins == weightsfiles.size());
  }


  //Loop over all bins
  for (uint i=0;i<fNMVABins; ++i) {
  
    TMVA::Reader *tmpTMVAReader = new TMVA::Reader( "!Color:!Silent:Error" );  
    tmpTMVAReader->SetVerbose(kTRUE);

    if (type == kIDIsoRingsCombined) {
      // Pure tracking variables
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p0To0p1",         &fMVAVar_ChargedIso_DR0p0To0p1        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p1To0p2",         &fMVAVar_ChargedIso_DR0p1To0p2        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p2To0p3",         &fMVAVar_ChargedIso_DR0p2To0p3        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p3To0p4",         &fMVAVar_ChargedIso_DR0p3To0p4        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p4To0p5",         &fMVAVar_ChargedIso_DR0p4To0p5        );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p0To0p1",           &fMVAVar_GammaIso_DR0p0To0p1          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p1To0p2",           &fMVAVar_GammaIso_DR0p1To0p2          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p2To0p3",           &fMVAVar_GammaIso_DR0p2To0p3          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p3To0p4",           &fMVAVar_GammaIso_DR0p3To0p4          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p4To0p5",           &fMVAVar_GammaIso_DR0p4To0p5          );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p0To0p1",   &fMVAVar_NeutralHadronIso_DR0p0To0p1  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p1To0p2",   &fMVAVar_NeutralHadronIso_DR0p1To0p2  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p2To0p3",   &fMVAVar_NeutralHadronIso_DR0p2To0p3  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p3To0p4",   &fMVAVar_NeutralHadronIso_DR0p3To0p4  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p4To0p5",   &fMVAVar_NeutralHadronIso_DR0p4To0p5  );
      tmpTMVAReader->AddVariable( "TkNchi2",                       &fMVAVar_MuTkNchi2               );
      if (i != 4) {
        tmpTMVAReader->AddVariable( "GlobalNchi2",                   &fMVAVar_MuGlobalNchi2           );
        tmpTMVAReader->AddVariable( "NValidHits",                    &fMVAVar_MuNValidHits            );
      }
      tmpTMVAReader->AddVariable( "NTrackerHits",                  &fMVAVar_MuNTrackerHits          );
      tmpTMVAReader->AddVariable( "NPixelHits",                    &fMVAVar_MuNPixelHits            );
      tmpTMVAReader->AddVariable( "NMatches",                      &fMVAVar_MuNMatches              );
      tmpTMVAReader->AddVariable( "TrkKink",                       &fMVAVar_MuTrkKink               );      
      tmpTMVAReader->AddVariable( "SegmentCompatibility",          &fMVAVar_MuSegmentCompatibility  );      
      tmpTMVAReader->AddVariable( "CaloCompatibility",             &fMVAVar_MuCaloCompatibility     );      
      tmpTMVAReader->AddVariable( "HadEnergy",                     &fMVAVar_MuHadEnergy       );      
      tmpTMVAReader->AddVariable( "EmEnergy",                      &fMVAVar_MuEmEnergy        );      
      tmpTMVAReader->AddVariable( "HadS9Energy",                   &fMVAVar_MuHadS9Energy     );      
      tmpTMVAReader->AddVariable( "EmS9Energy",                    &fMVAVar_MuEmS9Energy      );      
      tmpTMVAReader->AddSpectator("eta",                           &fMVAVar_MuEta);
      tmpTMVAReader->AddSpectator("pt",                            &fMVAVar_MuPt);
      tmpTMVAReader->AddSpectator("typeBits",                      &fMVAVar_MuTypeBits);
    }
  
    if (type == kIsoRings) {
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p0To0p1",         &fMVAVar_ChargedIso_DR0p0To0p1        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p1To0p2",         &fMVAVar_ChargedIso_DR0p1To0p2        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p2To0p3",         &fMVAVar_ChargedIso_DR0p2To0p3        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p3To0p4",         &fMVAVar_ChargedIso_DR0p3To0p4        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p4To0p5",         &fMVAVar_ChargedIso_DR0p4To0p5        );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p0To0p1",           &fMVAVar_GammaIso_DR0p0To0p1          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p1To0p2",           &fMVAVar_GammaIso_DR0p1To0p2          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p2To0p3",           &fMVAVar_GammaIso_DR0p2To0p3          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p3To0p4",           &fMVAVar_GammaIso_DR0p3To0p4          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p4To0p5",           &fMVAVar_GammaIso_DR0p4To0p5          );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p0To0p1",   &fMVAVar_NeutralHadronIso_DR0p0To0p1  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p1To0p2",   &fMVAVar_NeutralHadronIso_DR0p1To0p2  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p2To0p3",   &fMVAVar_NeutralHadronIso_DR0p2To0p3  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p3To0p4",   &fMVAVar_NeutralHadronIso_DR0p3To0p4  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p4To0p5",   &fMVAVar_NeutralHadronIso_DR0p4To0p5  );
      tmpTMVAReader->AddSpectator("eta",                           &fMVAVar_MuEta);
      tmpTMVAReader->AddSpectator("pt",                            &fMVAVar_MuPt);
    }
    
    if (type == kIsoDeltaR) {
      tmpTMVAReader->AddVariable("DZ",                            &fMVAVar_MuDZ              );
      tmpTMVAReader->AddVariable("IP2d",                          &fMVAVar_MuIP2d            );
      tmpTMVAReader->AddVariable("PFCharged",                     &fMVAVar_MuRelIsoPFCharged );
      tmpTMVAReader->AddVariable("PFNeutral",                     &fMVAVar_MuRelIsoPFNeutral );
      tmpTMVAReader->AddVariable("PFPhotons",                     &fMVAVar_MuRelIsoPFPhotons );
      tmpTMVAReader->AddVariable("DeltaR",                        &fMVAVar_MuDeltaRSum       );
      tmpTMVAReader->AddVariable("DeltaRMean",                    &fMVAVar_MuDeltaRMean      );
      tmpTMVAReader->AddVariable("Density",                       &fMVAVar_MuDensity         );
    }
    tmpTMVAReader->BookMVA(fMethodname , weightsfiles[i]);
    std::cout << "MVABin " << i << " : MethodName = " << fMethodname 
              << " , type == " << type << " , "
              << "Load weights file : " << weightsfiles[i] 
              << std::endl;
    fTMVAReader.push_back(tmpTMVAReader);
  }
  std::cout << "Muon ID MVA Completed\n";
}


//--------------------------------------------------------------------------------------------------
UInt_t MuonMVAEstimator::GetMVABin( double eta, double pt, Bool_t isGlobal, Bool_t isTrackerMuon) const {
  
    //Default is to return the first bin
    uint bin = 0;

    if (fMVAType == MuonMVAEstimator::kIsoRings) {
      if (pt < 10 && fabs(eta) < 1.479) bin = 0;
      if (pt < 10 && fabs(eta) >= 1.479) bin = 1;
      if (pt >= 10 && fabs(eta) < 1.479) bin = 2;
      if (pt >= 10 && fabs(eta) >= 1.479) bin = 3;
    }

    if (fMVAType == MuonMVAEstimator::kIDIsoRingsCombined ) {
      bin = 0;
      if (isGlobal && isTrackerMuon) {
        if (pt < 10 && fabs(eta) < 1.479) bin = 0;
        if (pt < 10 && fabs(eta) >= 1.479) bin = 1;
        if (pt >= 10 && fabs(eta) < 1.479) bin = 2;
        if (pt >= 10 && fabs(eta) >= 1.479) bin = 3;        
      } else if (!isGlobal && isTrackerMuon) {
        bin = 4;
      } else {
        cout << "Warning: Muon is not a tracker muon. Such muons are not supported. \n";
        bin = 0;
      }
    }
    if (fMVAType == MuonMVAEstimator::kIsoDeltaR){
      if (pt <  20 && fabs(eta) <  1.479) bin = 0;
      if (pt <  20 && fabs(eta) >= 1.479) bin = 1;
      if (pt >= 20 && fabs(eta) <  1.479) bin = 2;
      if (pt >= 20 && fabs(eta) >= 1.479) bin = 3;
    }
    return bin;
}

// //--------------------------------------------------------------------------------------------------
// Double_t MuonMVAEstimator::mvaValue(Double_t fbrem, 
// 					Double_t kfchi2,
// 					Int_t    kfhits,
// 					Double_t gsfchi2,
// 					Double_t deta,
// 					Double_t dphi,
// 					Double_t detacalo,
// 					//Double_t dphicalo,
// 					Double_t see,
// 					Double_t spp,
// 					Double_t etawidth,
// 					Double_t phiwidth,
// 					Double_t e1x5e5x5,
// 					Double_t R9,
// 					//Int_t    nbrems,
// 					Double_t HoE,
// 					Double_t EoP,
// 					Double_t IoEmIoP,
// 					Double_t eleEoPout,
// 					Double_t PreShowerOverRaw,
// 					//Double_t EoPout,
// 					Double_t eta,
// 					Double_t pt,
// 					Bool_t printDebug) {
  
//   if (!fisInitialized) { 
//     std::cout << "Error: MuonMVAEstimator not properly initialized.\n"; 
//     return -9999;
//   }

//   fMVAVar_fbrem           = fbrem; 
//   fMVAVar_kfchi2          = kfchi2;
//   fMVAVar_kfhits          = float(kfhits);   // BTD does not support int variables
//   fMVAVar_gsfchi2         = gsfchi2;

//   fMVAVar_deta            = deta;
//   fMVAVar_dphi            = dphi;
//   fMVAVar_detacalo        = detacalo;
//   // fMVAVar_dphicalo        = dphicalo;


//   fMVAVar_see             = see;
//   fMVAVar_spp             = spp;
//   fMVAVar_etawidth        = etawidth;
//   fMVAVar_phiwidth        = phiwidth;
//   fMVAVar_e1x5e5x5        = e1x5e5x5;
//   fMVAVar_R9              = R9;
//   //fMVAVar_nbrems          = float(nbrems);   // BTD does not support int variables


//   fMVAVar_HoE             = HoE;
//   fMVAVar_EoP             = EoP;
//   fMVAVar_IoEmIoP         = IoEmIoP;
//   fMVAVar_eleEoPout       = eleEoPout;
//   fMVAVar_PreShowerOverRaw= PreShowerOverRaw;
//   //fMVAVar_EoPout          = EoPout; 

//   fMVAVar_eta             = eta;
//   fMVAVar_pt              = pt;


//   bindVariables();
//   Double_t mva = -9999;  
//   if (fUseBinnedVersion) {
//     mva = fTMVAReader[GetMVABin(fMVAVar_eta,fMVAVar_pt)]->EvaluateMVA(fMethodname);
//   } else {
//     mva = fTMVAReader[0]->EvaluateMVA(fMethodname);
//   }



//   if(printDebug) {
//     cout << " *** Inside the class fMethodname " << fMethodname << endl;
//     cout << " fbrem " <<  fMVAVar_fbrem  
//       	 << " kfchi2 " << fMVAVar_kfchi2  
// 	 << " mykfhits " << fMVAVar_kfhits  
// 	 << " gsfchi2 " << fMVAVar_gsfchi2  
// 	 << " deta " <<  fMVAVar_deta  
// 	 << " dphi " << fMVAVar_dphi  
//       	 << " detacalo " << fMVAVar_detacalo  
//       // << " dphicalo " << fMVAVar_dphicalo  
// 	 << " see " << fMVAVar_see  
// 	 << " spp " << fMVAVar_spp  
// 	 << " etawidth " << fMVAVar_etawidth  
// 	 << " phiwidth " << fMVAVar_phiwidth  
// 	 << " e1x5e5x5 " << fMVAVar_e1x5e5x5  
// 	 << " R9 " << fMVAVar_R9  
//       // << " mynbrems " << fMVAVar_nbrems  
// 	 << " HoE " << fMVAVar_HoE  
// 	 << " EoP " << fMVAVar_EoP  
// 	 << " IoEmIoP " << fMVAVar_IoEmIoP  
// 	 << " eleEoPout " << fMVAVar_eleEoPout  
//       //<< " EoPout " << fMVAVar_EoPout  
// 	 << " PreShowerOverRaw " << fMVAVar_PreShowerOverRaw  
// 	 << " eta " << fMVAVar_eta  
// 	 << " pt " << fMVAVar_pt << endl;
//     cout << " ### MVA " << mva << endl;
//   }


//   return mva;
// }


//--------------------------------------------------------------------------------------------------
#ifndef STANDALONE
Double_t MuonMVAEstimator::mvaValue(const reco::Muon& mu, 
				    const reco::Vertex& vertex, 
				    const reco::PFCandidateCollection &PFCandidates,
				    double Rho,
				    MuonEffectiveArea::MuonEffectiveAreaTarget EATarget,
				    const reco::GsfElectronCollection &IdentifiedElectrons,
				    const reco::MuonCollection &IdentifiedMuons) {
  
  if (!fisInitialized) { 
    std::cout << "Error: MuonMVAEstimator not properly initialized.\n"; 
    return -9999;
  }

  if (fMVAType==2) return -9999;
  
  TrackRef muTrk = mu.track();
  if (muTrk.isNull()) {
    muTrk = mu.standAloneMuon();
  }
  if (muTrk.isNull()) {
    //if muon is not standalone either, then return -9999
    return -9999;
  }
  
  double muNchi2 = 0.0; 
  if (mu.combinedMuon().isNonnull()) { 
    muNchi2 = mu.combinedMuon()->chi2() / (Double_t)mu.combinedMuon()->ndof(); 
  } else if (mu.standAloneMuon().isNonnull()) { 
    muNchi2 = mu.standAloneMuon()->chi2() / (Double_t)mu.standAloneMuon()->ndof(); 
  } else if (mu.track().isNonnull()) { 
    muNchi2 = mu.track()->chi2() / (Double_t)mu.track()->ndof(); 
  }

  double rho = 0;
  if (!(isnan(float(Rho)) || isinf(float(Rho)))) rho = Rho;

  // Spectators
  fMVAVar_MuEta             =  muTrk->eta();         
  fMVAVar_MuPt              =  muTrk->pt();                          
  
  //set all input variables
  fMVAVar_MuTkNchi2              = muTrk->chi2() / (Double_t)muTrk->ndof();
  fMVAVar_MuGlobalNchi2          = muNchi2;
  fMVAVar_MuNValidHits           = mu.globalTrack().isNonnull() ? mu.globalTrack()->hitPattern().numberOfValidMuonHits() : 0;
  fMVAVar_MuNTrackerHits         = muTrk->numberOfValidHits();
  fMVAVar_MuNPixelHits           = muTrk->hitPattern().numberOfValidPixelHits();
  fMVAVar_MuNMatches             = mu.numberOfMatches();
  fMVAVar_MuTrkKink              = mu.combinedQuality().trkKink;
  fMVAVar_MuSegmentCompatibility = muon::segmentCompatibility(mu, reco::Muon::SegmentAndTrackArbitration);
  fMVAVar_MuCaloCompatibility    = mu.caloCompatibility();
  fMVAVar_MuHadEnergy             = mu.calEnergy().had;
  fMVAVar_MuEmEnergy       = mu.calEnergy().em;
  fMVAVar_MuHadS9Energy    = mu.calEnergy().hadS9;
  fMVAVar_MuEmS9Energy     = mu.calEnergy().emS9;
  
  //**********************************************************
  //Isolation variables
  //**********************************************************
  Double_t tmpChargedIso_DR0p0To0p1  = 0;
  Double_t tmpChargedIso_DR0p1To0p2  = 0;
  Double_t tmpChargedIso_DR0p2To0p3  = 0;
  Double_t tmpChargedIso_DR0p3To0p4  = 0;
  Double_t tmpChargedIso_DR0p4To0p5  = 0;
  Double_t tmpGammaIso_DR0p0To0p1  = 0;
  Double_t tmpGammaIso_DR0p1To0p2  = 0;
  Double_t tmpGammaIso_DR0p2To0p3  = 0;
  Double_t tmpGammaIso_DR0p3To0p4  = 0;
  Double_t tmpGammaIso_DR0p4To0p5  = 0;
  Double_t tmpNeutralHadronIso_DR0p0To0p1  = 0;
  Double_t tmpNeutralHadronIso_DR0p1To0p2  = 0;
  Double_t tmpNeutralHadronIso_DR0p2To0p3  = 0;
  Double_t tmpNeutralHadronIso_DR0p3To0p4  = 0;
  Double_t tmpNeutralHadronIso_DR0p4To0p5  = 0;
  
  double muonTrackZ = 0;
  if (mu.track().isNonnull()) {
    muonTrackZ = mu.track()->dz(vertex.position());
  }

  for (reco::PFCandidateCollection::const_iterator iP = PFCandidates.begin(); 
       iP != PFCandidates.end(); ++iP) {
      
    //exclude the muon itself
    if(iP->trackRef().isNonnull() && mu.innerTrack().isNonnull() &&
       refToPtr(iP->trackRef()) == refToPtr(mu.innerTrack())) continue;

    //************************************************************
    // New Isolation Calculations
    //************************************************************
    double dr = sqrt(pow(iP->eta() - mu.eta(),2) + pow(acos(cos(iP->phi() - mu.phi())),2));

    if (dr < 1.0) {
      Bool_t IsLeptonFootprint = kFALSE;
      //************************************************************
      // Lepton Footprint Removal
      //************************************************************   
      for (reco::GsfElectronCollection::const_iterator iE = IdentifiedElectrons.begin(); 
           iE != IdentifiedElectrons.end(); ++iE) {
	//if pf candidate matches an electron passing ID cuts, then veto it
	if(iP->gsfTrackRef().isNonnull() && iE->gsfTrack().isNonnull() &&
	   refToPtr(iP->gsfTrackRef()) == refToPtr(iE->gsfTrack())) IsLeptonFootprint = kTRUE;
        if(iP->trackRef().isNonnull() && iE->closestCtfTrackRef().isNonnull() &&
           refToPtr(iP->trackRef()) == refToPtr(iE->closestCtfTrackRef())) IsLeptonFootprint = kTRUE;

	//if pf candidate lies in veto regions of electron passing ID cuts, then veto it
        double tmpDR = sqrt(pow(iP->eta() - iE->eta(),2) + pow(acos(cos(iP->phi() - iE->phi())),2));
	if(iP->trackRef().isNonnull() && fabs(iE->superCluster()->eta()) >= 1.479 
           && tmpDR < 0.015) IsLeptonFootprint = kTRUE;
	if(iP->particleId() == reco::PFCandidate::gamma && fabs(iE->superCluster()->eta()) >= 1.479 
           && tmpDR < 0.08) IsLeptonFootprint = kTRUE;
      }
      for (reco::MuonCollection::const_iterator iM = IdentifiedMuons.begin(); 
           iM != IdentifiedMuons.end(); ++iM) {
	//if pf candidate matches an muon passing ID cuts, then veto it
	if(iP->trackRef().isNonnull() && iM->innerTrack().isNonnull() &&
	   refToPtr(iP->trackRef()) == refToPtr(iM->innerTrack())) IsLeptonFootprint = kTRUE;

	//if pf candidate lies in veto regions of muon passing ID cuts, then veto it
        double tmpDR = sqrt(pow(iP->eta() - iM->eta(),2) + pow(acos(cos(iP->phi() - iM->phi())),2));
	if(iP->trackRef().isNonnull() && tmpDR < 0.01) IsLeptonFootprint = kTRUE;
      }

     if (!IsLeptonFootprint) {
	Bool_t passVeto = kTRUE;
	//Charged
	 if(iP->trackRef().isNonnull()) {	  	   
	   if (!(fabs(iP->trackRef()->dz(vertex.position()) - muonTrackZ) < 0.2)) passVeto = kFALSE;
	   //************************************************************
	   // Veto any PFmuon, or PFEle
	   if (iP->particleId() == reco::PFCandidate::e || iP->particleId() == reco::PFCandidate::mu) passVeto = kFALSE;
	   //************************************************************
	   //************************************************************
	   // Footprint Veto
	   if (fabs(fMVAVar_MuEta) > 1.479 && dr < 0.01) passVeto = kFALSE;
	   //************************************************************
	   if (passVeto) {
	     if (dr < 0.1) tmpChargedIso_DR0p0To0p1 += iP->pt();
	     if (dr >= 0.1 && dr < 0.2) tmpChargedIso_DR0p1To0p2 += iP->pt();
	     if (dr >= 0.2 && dr < 0.3) tmpChargedIso_DR0p2To0p3 += iP->pt();
	     if (dr >= 0.3 && dr < 0.4) tmpChargedIso_DR0p3To0p4 += iP->pt();
	     if (dr >= 0.4 && dr < 0.5) tmpChargedIso_DR0p4To0p5 += iP->pt();
	   } //pass veto	   
	 }
	 //Gamma
	 else if (iP->particleId() == reco::PFCandidate::gamma) {
	     if (dr < 0.1) tmpGammaIso_DR0p0To0p1 += iP->pt();
	     if (dr >= 0.1 && dr < 0.2) tmpGammaIso_DR0p1To0p2 += iP->pt();
	     if (dr >= 0.2 && dr < 0.3) tmpGammaIso_DR0p2To0p3 += iP->pt();
	     if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += iP->pt();
	     if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += iP->pt();
	 }
	 //NeutralHadron
	 else {
           if (dr < 0.1) tmpNeutralHadronIso_DR0p0To0p1 += iP->pt();
           if (dr >= 0.1 && dr < 0.2) tmpNeutralHadronIso_DR0p1To0p2 += iP->pt();
           if (dr >= 0.2 && dr < 0.3) tmpNeutralHadronIso_DR0p2To0p3 += iP->pt();
           if (dr >= 0.3 && dr < 0.4) tmpNeutralHadronIso_DR0p3To0p4 += iP->pt();
           if (dr >= 0.4 && dr < 0.5) tmpNeutralHadronIso_DR0p4To0p5 += iP->pt();
	 }
      } //not lepton footprint
    } //in 1.0 dr cone
  } //loop over PF candidates

  fMVAVar_ChargedIso_DR0p0To0p1   = TMath::Min((tmpChargedIso_DR0p0To0p1)/mu.pt(), 2.5);
  fMVAVar_ChargedIso_DR0p1To0p2   = TMath::Min((tmpChargedIso_DR0p1To0p2)/mu.pt(), 2.5);
  fMVAVar_ChargedIso_DR0p2To0p3 = TMath::Min((tmpChargedIso_DR0p2To0p3)/mu.pt(), 2.5);
  fMVAVar_ChargedIso_DR0p3To0p4 = TMath::Min((tmpChargedIso_DR0p3To0p4)/mu.pt(), 2.5);
  fMVAVar_ChargedIso_DR0p4To0p5 = TMath::Min((tmpChargedIso_DR0p4To0p5)/mu.pt(), 2.5); 
  fMVAVar_GammaIso_DR0p0To0p1 = TMath::Max(TMath::Min((tmpGammaIso_DR0p0To0p1 - rho*MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaIsoDR0p0To0p1, fMVAVar_MuEta, EATarget))/mu.pt(), 2.5), 0.0);
  fMVAVar_GammaIso_DR0p1To0p2 = TMath::Max(TMath::Min((tmpGammaIso_DR0p1To0p2 - rho*MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaIsoDR0p1To0p2, fMVAVar_MuEta, EATarget))/mu.pt(), 2.5), 0.0);
  fMVAVar_GammaIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpGammaIso_DR0p2To0p3 - rho*MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaIsoDR0p2To0p3, fMVAVar_MuEta, EATarget))/mu.pt(), 2.5), 0.0);
  fMVAVar_GammaIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpGammaIso_DR0p3To0p4 - rho*MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaIsoDR0p3To0p4, fMVAVar_MuEta, EATarget))/mu.pt(), 2.5), 0.0);
  fMVAVar_GammaIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpGammaIso_DR0p4To0p5 - rho*MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaIsoDR0p4To0p5, fMVAVar_MuEta, EATarget))/mu.pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p0To0p1 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p0To0p1 - rho*MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuNeutralHadronIsoDR0p0To0p1, fMVAVar_MuEta, EATarget))/mu.pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p1To0p2 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p1To0p2 - rho*MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuNeutralHadronIsoDR0p1To0p2, fMVAVar_MuEta, EATarget))/mu.pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p2To0p3 - rho*MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuNeutralHadronIsoDR0p2To0p3, fMVAVar_MuEta, EATarget))/mu.pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p3To0p4 - rho*MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuNeutralHadronIsoDR0p3To0p4, fMVAVar_MuEta, EATarget))/mu.pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p4To0p5 - rho*MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuNeutralHadronIsoDR0p4To0p5, fMVAVar_MuEta, EATarget))/mu.pt(), 2.5), 0.0);
    
  if(fPrintMVADebug) {
    cout << fUseBinnedVersion << " -> BIN: " << fMVAVar_MuEta << " " << fMVAVar_MuPt << " "
         << "isGlobalMuon=" << mu.isGlobalMuon() << " " 
         << "isTrackerMuon=" << mu.isTrackerMuon() << " "
         << " : " << GetMVABin(fMVAVar_MuEta,fMVAVar_MuPt, mu.isGlobalMuon() , mu.isTrackerMuon()) << endl;
  }

  // evaluate
  bindVariables();
  Double_t mva = -9999; 
   
  if (fUseBinnedVersion) {    
    mva = fTMVAReader[GetMVABin(fMVAVar_MuEta,fMVAVar_MuPt, mu.isGlobalMuon() , mu.isTrackerMuon() )]->EvaluateMVA(fMethodname);
  } else {
    mva = fTMVAReader[0]->EvaluateMVA(fMethodname);
  }

  if(fPrintMVADebug) {
    cout << " *** Inside the class fMethodname " << fMethodname << " fMVAType " << fMVAType << endl;
    cout << " MuTkNchi2 " << fMVAVar_MuTkNchi2              
         << " MuGlobalNchi2 " << fMVAVar_MuGlobalNchi2          
         << " MuNValidHits " << fMVAVar_MuNValidHits           
         << " MuNTrackerHits " << fMVAVar_MuNTrackerHits         
         << " MuNPixelHits " << fMVAVar_MuNPixelHits           
         << " MuNMatches " << fMVAVar_MuNMatches             
         << " MuTrkKink " << fMVAVar_MuTrkKink              
         << " MuSegmentCompatibility " << fMVAVar_MuSegmentCompatibility 
         << " MuCaloCompatibility " << fMVAVar_MuCaloCompatibility    
         << " MuHadEnergy " << fMVAVar_MuHadEnergy      
         << " MuEmEnergy " << fMVAVar_MuEmEnergy       
         << " MuHadS9Energy " << fMVAVar_MuHadS9Energy    
         << " MuEmS9Energy " << fMVAVar_MuEmS9Energy     
	 << " eta " << fMVAVar_MuEta  
	 << " pt " << fMVAVar_MuPt << endl;
    cout  << "ChargedIso ( 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 ): " 
          << fMVAVar_ChargedIso_DR0p0To0p1   << " "
          << fMVAVar_ChargedIso_DR0p1To0p2   << " "
          << fMVAVar_ChargedIso_DR0p2To0p3 << " "
          << fMVAVar_ChargedIso_DR0p3To0p4 << " "
          << fMVAVar_ChargedIso_DR0p4To0p5 << endl;
    cout  << "PF Gamma Iso ( 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 ): " 
          << fMVAVar_GammaIso_DR0p0To0p1 << " "
          << fMVAVar_GammaIso_DR0p1To0p2 << " "
          << fMVAVar_GammaIso_DR0p2To0p3 << " "
          << fMVAVar_GammaIso_DR0p3To0p4 << " "
          << fMVAVar_GammaIso_DR0p4To0p5 << endl;
    cout  << "PF Neutral Hadron Iso ( 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 ): " 
          << fMVAVar_NeutralHadronIso_DR0p0To0p1 << " "
          << fMVAVar_NeutralHadronIso_DR0p1To0p2 << " "
          << fMVAVar_NeutralHadronIso_DR0p2To0p3 << " "
          << fMVAVar_NeutralHadronIso_DR0p3To0p4 << " "
          << fMVAVar_NeutralHadronIso_DR0p4To0p5 << " "
          << endl;
    cout << " ### MVA " << mva << endl;
  }
  

  return mva;
}
// THIS METHOD TAKES THE PfNoPileUp collection
Double_t MuonMVAEstimator::mvaValue(const reco::Muon& mu, 
				    const reco::Vertex& vertex, 
				    const reco::PFCandidateCollection &PFCandidates,
				    double Rho,
				    MuonEffectiveArea::MuonEffectiveAreaTarget EATarget){
  
  if (!fisInitialized) { 
    std::cout << "Error: MuonMVAEstimator not properly initialized.\n"; 
    return -9999;
  }

  if (fMVAType!=2) return -9999;

  TrackRef muTrk = mu.track();
  if (muTrk.isNull()) {
    muTrk = mu.standAloneMuon();
  }
  if (muTrk.isNull()) {
    //if muon is not standalone either, then return -9999
    return -9999;
  }
  
  double rho = 0;
  if (!(isnan(float(Rho)) || isinf(float(Rho)))) rho = Rho;

  Double_t tmpMuRelIsoPFCharged = 0;
  Double_t tmpMuRelIsoPFNeutral = 0;
  Double_t tmpMuRelIsoPFPhotons = 0;
  Double_t tmpMuDeltaRMean = 0;
  Double_t tmpMuDeltaRSum = 0;
  Double_t tmpMuDensity = 0;
  Double_t tmpMuDZ = 0;
  Double_t tmpMuIP2d = 0;
  Double_t tmpMuNPFCand = 0;
  
  if (mu.track().isNonnull()) {
    tmpMuDZ = mu.track()->dz(vertex.position());
    tmpMuIP2d = mu.track()->dxy(vertex.position());
  }
  
  for (reco::PFCandidateCollection::const_iterator iP = PFCandidates.begin(); iP != PFCandidates.end(); ++iP) {
    //exclude the muon itself
    if(iP->trackRef().isNonnull() && mu.innerTrack().isNonnull() &&
       refToPtr(iP->trackRef()) == refToPtr(mu.innerTrack())) continue;

    //************************************************************
    // New Isolation Calculations
    //************************************************************
    double dr = sqrt(pow(iP->eta() - mu.eta(),2) + pow(acos(cos(iP->phi() - mu.phi())),2));

    if (dr > 0.5)  continue;
    if (dr < 0.01) continue; 
    
    //Charged
    if(iP->trackRef().isNonnull()) {	  	   
      //************************************************************
      // Veto any PFmuon, or PFEle
      if (iP->particleId() == reco::PFCandidate::e || iP->particleId() == reco::PFCandidate::mu) continue;
      //************************************************************
      tmpMuRelIsoPFCharged += iP->pt(); 
    }
    //Gamma
    else if (iP->particleId() == reco::PFCandidate::gamma && iP->pt() > 1.0) {
      tmpMuRelIsoPFPhotons += iP->pt();
    }
    //NeutralHadron
    else if (iP->pt() > 1.0) {
      tmpMuRelIsoPFNeutral += iP->pt();
    } 
    
    tmpMuNPFCand++;
    tmpMuDeltaRMean += dr;
    tmpMuDeltaRSum  += dr;
    tmpMuDensity    += iP->pt() / dr;
  } //loop over PF candidates
  
  fMVAVar_MuRelIsoPFCharged = (tmpMuRelIsoPFCharged) / mu.pt();  
  fMVAVar_MuRelIsoPFNeutral = TMath::Max((tmpMuRelIsoPFNeutral - rho*MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuNeutralIso05, fMVAVar_MuEta, EATarget))/mu.pt(),0.0);
  fMVAVar_MuRelIsoPFPhotons = TMath::Max((tmpMuRelIsoPFPhotons - rho*MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaIso05,   fMVAVar_MuEta, EATarget))/mu.pt(),0.0);
  fMVAVar_MuDeltaRMean      = tmpMuDeltaRMean/TMath::Max(1.0,tmpMuNPFCand);
  fMVAVar_MuDeltaRSum       = tmpMuDeltaRSum;
  fMVAVar_MuDensity         = tmpMuDensity;
  fMVAVar_MuDZ              = tmpMuDZ;
  fMVAVar_MuIP2d            = tmpMuIP2d;

  if(fPrintMVADebug) {
    cout << fUseBinnedVersion << " -> BIN: " << fMVAVar_MuEta << " " << fMVAVar_MuPt << " "
         << "isGlobalMuon=" << mu.isGlobalMuon() << " " 
         << "isTrackerMuon=" << mu.isTrackerMuon() << " "
         << " : " << GetMVABin(fMVAVar_MuEta,fMVAVar_MuPt, mu.isGlobalMuon() , mu.isTrackerMuon()) << endl;
  }
  
  // evaluate
  bindVariables();
  Double_t mva = -9999; 
  
  if (fUseBinnedVersion) {    
    mva = fTMVAReader[GetMVABin(fMVAVar_MuEta,fMVAVar_MuPt, mu.isGlobalMuon() , mu.isTrackerMuon() )]->EvaluateMVA(fMethodname);
  } else {
    mva = fTMVAReader[0]->EvaluateMVA(fMethodname);
  }
  
  return mva;
}
#endif
void MuonMVAEstimator::bindVariables() {
  
  // this binding is needed for variables that sometime diverge. 
  
  
  return;
}








