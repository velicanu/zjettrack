#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TArrow.h"

#include <iostream>
#include <vector>

const double pi = 3.14159265358979323846;

using namespace std;

float getDphi(float phi1, float phi2) {
  float dphi = phi1 - phi2;
  if ( dphi > pi )
    dphi = dphi - 2.*pi;
  if ( dphi <= -pi ) 
    dphi = dphi + 2.*pi;
  return dphi;
}

//variables branches in ggNtuple
   UInt_t          _run;
   ULong64_t       _event;
   UInt_t          _lumis;
   Bool_t          _isData;
   Int_t           _nEle;
   vector<int>     *_eleCharge;
   vector<int>     *_eleChargeConsistent;
   vector<float>   *_eleEn;
   vector<float>   *_eleD0;
   vector<float>   *_eleDz;
   vector<float>   *_eleD0Err;
   vector<float>   *_eleDzErr;
   vector<float>   *_eleTrkPt;
   vector<float>   *_eleTrkEta;
   vector<float>   *_eleTrkPhi;
   vector<int>     *_eleTrkCharge;
   vector<float>   *_eleTrkChi2;
   vector<float>   *_eleTrkNdof;
   vector<float>   *_eleTrkNormalizedChi2;
   vector<int>     *_eleTrkValidHits;
   vector<int>     *_eleTrkLayers;
   vector<float>   *_elePt;
   vector<float>   *_eleEta;
   vector<float>   *_elePhi;
   vector<float>   *_eleSCEn;
   vector<float>   *_eleESEn;
   vector<float>   *_eleSCEta;
   vector<float>   *_eleSCPhi;
   vector<float>   *_eleSCRawEn;
   vector<float>   *_eleSCEtaWidth;
   vector<float>   *_eleSCPhiWidth;
   vector<float>   *_eleHoverE;
   vector<float>   *_eleEoverP;
   vector<float>   *_eleEoverPInv;
   vector<float>   *_eleBrem;
   vector<float>   *_eledEtaAtVtx;
   vector<float>   *_eledPhiAtVtx;
   vector<float>   *_eleSigmaIEtaIEta;
   vector<float>   *_eleSigmaIEtaIEta_2012;
   vector<float>   *_eleSigmaIPhiIPhi;
   vector<int>     *_eleMissHits;
   vector<float>   *_eleESEffSigmaRR;
   vector<float>   *_elePFChIso;
   vector<float>   *_elePFPhoIso;
   vector<float>   *_elePFNeuIso;
   vector<float>   *_elePFPUIso;
   vector<float>   *_elePFChIso03;
   vector<float>   *_elePFPhoIso03;
   vector<float>   *_elePFNeuIso03;
   vector<float>   *_elePFChIso04;
   vector<float>   *_elePFPhoIso04;
   vector<float>   *_elePFNeuIso04;
   vector<float>   *_eleBC1E;
   vector<float>   *_eleBC1Eta;
   vector<float>   *_eleBC2E;
   vector<float>   *_eleBC2Eta;
   Int_t           _nPho;
   vector<float>   *_phoE;
   vector<float>   *_phoEt;
   vector<float>   *_phoEta;
   vector<float>   *_phoPhi;
   vector<float>   *_phoSCE;
   vector<float>   *_phoSCRawE;
   vector<float>   *_phoESEn;
   vector<float>   *_phoSCEta;
   vector<float>   *_phoSCPhi;
   vector<float>   *_phoSCEtaWidth;
   vector<float>   *_phoSCPhiWidth;
   vector<float>   *_phoSCBrem;
   vector<int>     *_phohasPixelSeed;
   vector<float>   *_phoR9;
   vector<float>   *_phoHoverE;
   vector<float>   *_phoSigmaIEtaIEta;
   vector<float>   *_phoE1x3;
   vector<float>   *_phoE2x2;
   vector<float>   *_phoE3x3;
   vector<float>   *_phoE2x5Max;
   vector<float>   *_phoE1x5;
   vector<float>   *_phoE2x5;
   vector<float>   *_phoE5x5;
   vector<float>   *_phoMaxEnergyXtal;
   vector<float>   *_phoSigmaEtaEta;
   vector<float>   *_phoR1x5;
   vector<float>   *_phoR2x5;
   vector<float>   *_phoESEffSigmaRR;
   vector<float>   *_phoSigmaIEtaIEta_2012;
   vector<float>   *_phoSigmaIEtaIPhi_2012;
   vector<float>   *_phoSigmaIPhiIPhi_2012;
   vector<float>   *_phoE1x3_2012;
   vector<float>   *_phoE2x2_2012;
   vector<float>   *_phoE3x3_2012;
   vector<float>   *_phoE2x5Max_2012;
   vector<float>   *_phoE5x5_2012;
   vector<float>   *_phoBC1E;
   vector<float>   *_phoBC1Eta;
   vector<float>   *_phoBC2E;
   vector<float>   *_phoBC2Eta;
   vector<float>   *_pho_ecalClusterIsoR2;
   vector<float>   *_pho_ecalClusterIsoR3;
   vector<float>   *_pho_ecalClusterIsoR4;
   vector<float>   *_pho_ecalClusterIsoR5;
   vector<float>   *_pho_hcalRechitIsoR1;
   vector<float>   *_pho_hcalRechitIsoR2;
   vector<float>   *_pho_hcalRechitIsoR3;
   vector<float>   *_pho_hcalRechitIsoR4;
   vector<float>   *_pho_hcalRechitIsoR5;
   vector<float>   *_pho_trackIsoR1PtCut20;
   vector<float>   *_pho_trackIsoR2PtCut20;
   vector<float>   *_pho_trackIsoR3PtCut20;
   vector<float>   *_pho_trackIsoR4PtCut20;
   vector<float>   *_pho_trackIsoR5PtCut20;
   vector<float>   *_pho_swissCrx;
   vector<float>   *_pho_seedTime;
   vector<float>   *_pfcIso1;
   vector<float>   *_pfcIso2;
   vector<float>   *_pfcIso3;
   vector<float>   *_pfcIso4;
   vector<float>   *_pfcIso5;
   vector<float>   *_pfpIso1;
   vector<float>   *_pfpIso2;
   vector<float>   *_pfpIso3;
   vector<float>   *_pfpIso4;
   vector<float>   *_pfpIso5;
   vector<float>   *_pfnIso1;
   vector<float>   *_pfnIso2;
   vector<float>   *_pfnIso3;
   vector<float>   *_pfnIso4;
   vector<float>   *_pfnIso5;
   Int_t				   _nMC;
	 vector<int>     *_mcPID;
   vector<int>     *_mcStatus;
	 vector<float>   *_mcPt;
   vector<float>   *_mcEta;
   vector<float>   *_mcPhi;
	 vector<float>   *_mcMomPt;
   vector<float>   *_mcMomEta;
   vector<float>   *_mcMomPhi;
   vector<int>     *_mcMomPID;

   Int_t           _nMu;
   vector<float>   *_muPt;
   vector<float>   *_muEta;
   vector<float>   *_muPhi;
   vector<int>     *_muCharge;
   vector<int>     *_muType;
   vector<int>     *_muIsGood;
   vector<float>   *_muD0;
   vector<float>   *_muDz;
   vector<float>   *_muChi2NDF;
   vector<float>   *_muInnerD0;
   vector<float>   *_muInnerDz;
   vector<int>     *_muTrkLayers;
   vector<int>     *_muPixelLayers;
   vector<int>     *_muPixelHits;
   vector<int>     *_muMuonHits;
   vector<int>     *_muTrkQuality;
   vector<int>     *_muStations;
   vector<float>   *_muIsoTrk;
   vector<float>   *_muPFChIso;
   vector<float>   *_muPFPhoIso;
   vector<float>   *_muPFNeuIso;
   vector<float>   *_muPFPUIso;

void initggTree(TTree *tree) {

   _eleCharge = 0;
   _eleChargeConsistent = 0;
   _eleEn = 0;
   _eleD0 = 0;
   _eleDz = 0;
   _eleD0Err = 0;
   _eleDzErr = 0;
   _eleTrkPt = 0;
   _eleTrkEta = 0;
   _eleTrkPhi = 0;
   _eleTrkCharge = 0;
   _eleTrkChi2 = 0;
   _eleTrkNdof = 0;
   _eleTrkNormalizedChi2 = 0;
   _eleTrkValidHits = 0;
   _eleTrkLayers = 0;
   _elePt = 0;
   _eleEta = 0;
   _elePhi = 0;
   _eleSCEn = 0;
   _eleESEn = 0;
   _eleSCEta = 0;
   _eleSCPhi = 0;
   _eleSCRawEn = 0;
   _eleSCEtaWidth = 0;
   _eleSCPhiWidth = 0;
   _eleHoverE = 0;
   _eleEoverP = 0;
   _eleEoverPInv = 0;
   _eleBrem = 0;
   _eledEtaAtVtx = 0;
   _eledPhiAtVtx = 0;
   _eleSigmaIEtaIEta = 0;
   _eleSigmaIEtaIEta_2012 = 0;
   _eleSigmaIPhiIPhi = 0;
   _eleMissHits = 0;
   _eleESEffSigmaRR = 0;
   _elePFChIso = 0;
   _elePFPhoIso = 0;
   _elePFNeuIso = 0;
   _elePFPUIso = 0;
   _elePFChIso03 = 0;
   _elePFPhoIso03 = 0;
   _elePFNeuIso03 = 0;
   _elePFChIso04 = 0;
   _elePFPhoIso04 = 0;
   _elePFNeuIso04 = 0;
   _eleBC1E = 0;
   _eleBC1Eta = 0;
   _eleBC2E = 0;
   _eleBC2Eta = 0;
   _phoE = 0;
   _phoEt = 0;
   _phoEta = 0;
   _phoPhi = 0;
   _phoSCE = 0;
   _phoSCRawE = 0;
   _phoESEn = 0;
   _phoSCEta = 0;
   _phoSCPhi = 0;
   _phoSCEtaWidth = 0;
   _phoSCPhiWidth = 0;
   _phoSCBrem = 0;
   _phohasPixelSeed = 0;
   _phoR9 = 0;
   _phoHoverE = 0;
   _phoSigmaIEtaIEta = 0;
   _phoE1x3 = 0;
   _phoE2x2 = 0;
   _phoE3x3 = 0;
   _phoE2x5Max = 0;
   _phoE1x5 = 0;
   _phoE2x5 = 0;
   _phoE5x5 = 0;
   _phoMaxEnergyXtal = 0;
   _phoSigmaEtaEta = 0;
   _phoR1x5 = 0;
   _phoR2x5 = 0;
   _phoESEffSigmaRR = 0;
   _phoSigmaIEtaIEta_2012 = 0;
   _phoSigmaIEtaIPhi_2012 = 0;
   _phoSigmaIPhiIPhi_2012 = 0;
   _phoE1x3_2012 = 0;
   _phoE2x2_2012 = 0;
   _phoE3x3_2012 = 0;
   _phoE2x5Max_2012 = 0;
   _phoE5x5_2012 = 0;
   _phoBC1E = 0;
   _phoBC1Eta = 0;
   _phoBC2E = 0;
   _phoBC2Eta = 0;
   _pho_ecalClusterIsoR2 = 0;
   _pho_ecalClusterIsoR3 = 0;
   _pho_ecalClusterIsoR4 = 0;
   _pho_ecalClusterIsoR5 = 0;
   _pho_hcalRechitIsoR1 = 0;
   _pho_hcalRechitIsoR2 = 0;
   _pho_hcalRechitIsoR3 = 0;
   _pho_hcalRechitIsoR4 = 0;
   _pho_hcalRechitIsoR5 = 0;
   _pho_trackIsoR1PtCut20 = 0;
   _pho_trackIsoR2PtCut20 = 0;
   _pho_trackIsoR3PtCut20 = 0;
   _pho_trackIsoR4PtCut20 = 0;
   _pho_trackIsoR5PtCut20 = 0;
   _pho_swissCrx = 0;
   _pho_seedTime = 0;
   _pfcIso1 = 0;
   _pfcIso2 = 0;
   _pfcIso3 = 0;
   _pfcIso4 = 0;
   _pfcIso5 = 0;
   _pfpIso1 = 0;
   _pfpIso2 = 0;
   _pfpIso3 = 0;
   _pfpIso4 = 0;
   _pfpIso5 = 0;
   _pfnIso1 = 0;
   _pfnIso2 = 0;
   _pfnIso3 = 0;
   _pfnIso4 = 0;
   _pfnIso5 = 0;
	 _nMC = 0;
	 _mcPID = 0;
   _mcStatus = 0;
	 _mcPt = 0;
   _mcEta = 0;
   _mcPhi = 0;
	 _mcMomPt = 0;
   _mcMomEta = 0;
   _mcMomPhi = 0;
   _mcMomPID = 0;
   _muPt = 0;
   _muEta = 0;
   _muPhi = 0;
   _muCharge = 0;
   _muType = 0;
   _muIsGood = 0;
   _muD0 = 0;
   _muDz = 0;
   _muChi2NDF = 0;
   _muInnerD0 = 0;
   _muInnerDz = 0;
   _muTrkLayers = 0;
   _muPixelLayers = 0;
   _muPixelHits = 0;
   _muMuonHits = 0;
   _muTrkQuality = 0;
   _muStations = 0;
   _muIsoTrk = 0;
   _muPFChIso = 0;
   _muPFPhoIso = 0;
   _muPFNeuIso = 0;
   _muPFPUIso = 0;

   tree->SetBranchAddress("run", &_run);
   tree->SetBranchAddress("event", &_event);
   tree->SetBranchAddress("lumis", &_lumis);
   tree->SetBranchAddress("isData", &_isData);
   tree->SetBranchAddress("nEle", &_nEle);
   tree->SetBranchAddress("eleCharge", &_eleCharge);
   tree->SetBranchAddress("eleChargeConsistent", &_eleChargeConsistent);
   tree->SetBranchAddress("eleEn", &_eleEn);
   tree->SetBranchAddress("eleD0", &_eleD0);
   tree->SetBranchAddress("eleDz", &_eleDz);
   tree->SetBranchAddress("eleD0Err", &_eleD0Err);
   tree->SetBranchAddress("eleDzErr", &_eleDzErr);
   tree->SetBranchAddress("eleTrkPt", &_eleTrkPt);
   tree->SetBranchAddress("eleTrkEta", &_eleTrkEta);
   tree->SetBranchAddress("eleTrkPhi", &_eleTrkPhi);
   tree->SetBranchAddress("eleTrkCharge", &_eleTrkCharge);
   tree->SetBranchAddress("eleTrkChi2", &_eleTrkChi2);
   tree->SetBranchAddress("eleTrkNdof", &_eleTrkNdof);
   tree->SetBranchAddress("eleTrkNormalizedChi2", &_eleTrkNormalizedChi2);
   tree->SetBranchAddress("eleTrkValidHits", &_eleTrkValidHits);
   tree->SetBranchAddress("eleTrkLayers", &_eleTrkLayers);
   tree->SetBranchAddress("elePt", &_elePt);
   tree->SetBranchAddress("eleEta", &_eleEta);
   tree->SetBranchAddress("elePhi", &_elePhi);
   tree->SetBranchAddress("eleSCEn", &_eleSCEn);
   tree->SetBranchAddress("eleESEn", &_eleESEn);
   tree->SetBranchAddress("eleSCEta", &_eleSCEta);
   tree->SetBranchAddress("eleSCPhi", &_eleSCPhi);
   tree->SetBranchAddress("eleSCRawEn", &_eleSCRawEn);
   tree->SetBranchAddress("eleSCEtaWidth", &_eleSCEtaWidth);
   tree->SetBranchAddress("eleSCPhiWidth", &_eleSCPhiWidth);
   tree->SetBranchAddress("eleHoverE", &_eleHoverE);
   tree->SetBranchAddress("eleEoverP", &_eleEoverP);
   tree->SetBranchAddress("eleEoverPInv", &_eleEoverPInv);
   tree->SetBranchAddress("eleBrem", &_eleBrem);
   tree->SetBranchAddress("eledEtaAtVtx", &_eledEtaAtVtx);
   tree->SetBranchAddress("eledPhiAtVtx", &_eledPhiAtVtx);
   tree->SetBranchAddress("eleSigmaIEtaIEta", &_eleSigmaIEtaIEta);
   tree->SetBranchAddress("eleSigmaIEtaIEta_2012", &_eleSigmaIEtaIEta_2012);
   tree->SetBranchAddress("eleSigmaIPhiIPhi", &_eleSigmaIPhiIPhi);
   tree->SetBranchAddress("eleMissHits", &_eleMissHits);
   tree->SetBranchAddress("eleESEffSigmaRR", &_eleESEffSigmaRR);
   tree->SetBranchAddress("elePFChIso", &_elePFChIso);
   tree->SetBranchAddress("elePFPhoIso", &_elePFPhoIso);
   tree->SetBranchAddress("elePFNeuIso", &_elePFNeuIso);
   tree->SetBranchAddress("elePFPUIso", &_elePFPUIso);
   tree->SetBranchAddress("elePFChIso03", &_elePFChIso03);
   tree->SetBranchAddress("elePFPhoIso03", &_elePFPhoIso03);
   tree->SetBranchAddress("elePFNeuIso03", &_elePFNeuIso03);
   tree->SetBranchAddress("elePFChIso04", &_elePFChIso04);
   tree->SetBranchAddress("elePFPhoIso04", &_elePFPhoIso04);
   tree->SetBranchAddress("elePFNeuIso04", &_elePFNeuIso04);
   tree->SetBranchAddress("eleBC1E", &_eleBC1E);
   tree->SetBranchAddress("eleBC1Eta", &_eleBC1Eta);
   tree->SetBranchAddress("eleBC2E", &_eleBC2E);
   tree->SetBranchAddress("eleBC2Eta", &_eleBC2Eta);
   tree->SetBranchAddress("nPho", &_nPho);
   tree->SetBranchAddress("phoE", &_phoE);
   tree->SetBranchAddress("phoEt", &_phoEt);
   tree->SetBranchAddress("phoEta", &_phoEta);
   tree->SetBranchAddress("phoPhi", &_phoPhi);
   tree->SetBranchAddress("phoSCE", &_phoSCE);
   tree->SetBranchAddress("phoSCRawE", &_phoSCRawE);
   tree->SetBranchAddress("phoESEn", &_phoESEn);
   tree->SetBranchAddress("phoSCEta", &_phoSCEta);
   tree->SetBranchAddress("phoSCPhi", &_phoSCPhi);
   tree->SetBranchAddress("phoSCEtaWidth", &_phoSCEtaWidth);
   tree->SetBranchAddress("phoSCPhiWidth", &_phoSCPhiWidth);
   tree->SetBranchAddress("phoSCBrem", &_phoSCBrem);
   tree->SetBranchAddress("phohasPixelSeed", &_phohasPixelSeed);
   tree->SetBranchAddress("phoR9", &_phoR9);
   tree->SetBranchAddress("phoHoverE", &_phoHoverE);
   tree->SetBranchAddress("phoSigmaIEtaIEta", &_phoSigmaIEtaIEta);
   tree->SetBranchAddress("phoE1x3", &_phoE1x3);
   tree->SetBranchAddress("phoE2x2", &_phoE2x2);
   tree->SetBranchAddress("phoE3x3", &_phoE3x3);
   tree->SetBranchAddress("phoE2x5Max", &_phoE2x5Max);
   tree->SetBranchAddress("phoE1x5", &_phoE1x5);
   tree->SetBranchAddress("phoE2x5", &_phoE2x5);
   tree->SetBranchAddress("phoE5x5", &_phoE5x5);
   tree->SetBranchAddress("phoMaxEnergyXtal", &_phoMaxEnergyXtal);
   tree->SetBranchAddress("phoSigmaEtaEta", &_phoSigmaEtaEta);
   tree->SetBranchAddress("phoR1x5", &_phoR1x5);
   tree->SetBranchAddress("phoR2x5", &_phoR2x5);
   tree->SetBranchAddress("phoESEffSigmaRR", &_phoESEffSigmaRR);
   tree->SetBranchAddress("phoSigmaIEtaIEta_2012", &_phoSigmaIEtaIEta_2012);
   tree->SetBranchAddress("phoSigmaIEtaIPhi_2012", &_phoSigmaIEtaIPhi_2012);
   tree->SetBranchAddress("phoSigmaIPhiIPhi_2012", &_phoSigmaIPhiIPhi_2012);
   tree->SetBranchAddress("phoE1x3_2012", &_phoE1x3_2012);
   tree->SetBranchAddress("phoE2x2_2012", &_phoE2x2_2012);
   tree->SetBranchAddress("phoE3x3_2012", &_phoE3x3_2012);
   tree->SetBranchAddress("phoE2x5Max_2012", &_phoE2x5Max_2012);
   tree->SetBranchAddress("phoE5x5_2012", &_phoE5x5_2012);
   tree->SetBranchAddress("phoBC1E", &_phoBC1E);
   tree->SetBranchAddress("phoBC1Eta", &_phoBC1Eta);
   tree->SetBranchAddress("phoBC2E", &_phoBC2E);
   tree->SetBranchAddress("phoBC2Eta", &_phoBC2Eta);
   tree->SetBranchAddress("pho_ecalClusterIsoR2", &_pho_ecalClusterIsoR2);
   tree->SetBranchAddress("pho_ecalClusterIsoR3", &_pho_ecalClusterIsoR3);
   tree->SetBranchAddress("pho_ecalClusterIsoR4", &_pho_ecalClusterIsoR4);
   tree->SetBranchAddress("pho_ecalClusterIsoR5", &_pho_ecalClusterIsoR5);
   tree->SetBranchAddress("pho_hcalRechitIsoR1", &_pho_hcalRechitIsoR1);
   tree->SetBranchAddress("pho_hcalRechitIsoR2", &_pho_hcalRechitIsoR2);
   tree->SetBranchAddress("pho_hcalRechitIsoR3", &_pho_hcalRechitIsoR3);
   tree->SetBranchAddress("pho_hcalRechitIsoR4", &_pho_hcalRechitIsoR4);
   tree->SetBranchAddress("pho_hcalRechitIsoR5", &_pho_hcalRechitIsoR5);
   tree->SetBranchAddress("pho_trackIsoR1PtCut20", &_pho_trackIsoR1PtCut20);
   tree->SetBranchAddress("pho_trackIsoR2PtCut20", &_pho_trackIsoR2PtCut20);
   tree->SetBranchAddress("pho_trackIsoR3PtCut20", &_pho_trackIsoR3PtCut20);
   tree->SetBranchAddress("pho_trackIsoR4PtCut20", &_pho_trackIsoR4PtCut20);
   tree->SetBranchAddress("pho_trackIsoR5PtCut20", &_pho_trackIsoR5PtCut20);
   tree->SetBranchAddress("pho_swissCrx", &_pho_swissCrx);
   tree->SetBranchAddress("pho_seedTime", &_pho_seedTime);
   tree->SetBranchAddress("pfcIso1", &_pfcIso1);
   tree->SetBranchAddress("pfcIso2", &_pfcIso2);
   tree->SetBranchAddress("pfcIso3", &_pfcIso3);
   tree->SetBranchAddress("pfcIso4", &_pfcIso4);
   tree->SetBranchAddress("pfcIso5", &_pfcIso5);
   tree->SetBranchAddress("pfpIso1", &_pfpIso1);
   tree->SetBranchAddress("pfpIso2", &_pfpIso2);
   tree->SetBranchAddress("pfpIso3", &_pfpIso3);
   tree->SetBranchAddress("pfpIso4", &_pfpIso4);
   tree->SetBranchAddress("pfpIso5", &_pfpIso5);
   tree->SetBranchAddress("pfnIso1", &_pfnIso1);
   tree->SetBranchAddress("pfnIso2", &_pfnIso2);
   tree->SetBranchAddress("pfnIso3", &_pfnIso3);
   tree->SetBranchAddress("pfnIso4", &_pfnIso4);
   tree->SetBranchAddress("pfnIso5", &_pfnIso5);
	 tree->SetBranchAddress("nMC", &_nMC);
	 tree->SetBranchAddress("mcPID", &_mcPID);
   tree->SetBranchAddress("mcStatus", &_mcStatus);
	 tree->SetBranchAddress("mcPt", &_mcPt);
   tree->SetBranchAddress("mcEta", &_mcEta);
   tree->SetBranchAddress("mcPhi", &_mcPhi);
	 tree->SetBranchAddress("mcMomPt", &_mcMomPt);
   tree->SetBranchAddress("mcMomEta", &_mcMomEta);
   tree->SetBranchAddress("mcMomPhi", &_mcMomPhi);
   tree->SetBranchAddress("mcMomPID", &_mcMomPID);
   tree->SetBranchAddress("nMu", &_nMu);
   tree->SetBranchAddress("muPt", &_muPt);
   tree->SetBranchAddress("muEta", &_muEta);
   tree->SetBranchAddress("muPhi", &_muPhi);
   tree->SetBranchAddress("muCharge", &_muCharge);
   tree->SetBranchAddress("muType", &_muType);
   tree->SetBranchAddress("muIsGood", &_muIsGood);
   tree->SetBranchAddress("muD0", &_muD0);
   tree->SetBranchAddress("muDz", &_muDz);
   tree->SetBranchAddress("muChi2NDF", &_muChi2NDF);
   tree->SetBranchAddress("muInnerD0", &_muInnerD0);
   tree->SetBranchAddress("muInnerDz", &_muInnerDz);
   tree->SetBranchAddress("muTrkLayers", &_muTrkLayers);
   tree->SetBranchAddress("muPixelLayers", &_muPixelLayers);
   tree->SetBranchAddress("muPixelHits", &_muPixelHits);
   tree->SetBranchAddress("muMuonHits", &_muMuonHits);
   tree->SetBranchAddress("muTrkQuality", &_muTrkQuality);
   tree->SetBranchAddress("muStations", &_muStations);
   tree->SetBranchAddress("muIsoTrk", &_muIsoTrk);
   tree->SetBranchAddress("muPFChIso", &_muPFChIso);
   tree->SetBranchAddress("muPFPhoIso", &_muPFPhoIso);
   tree->SetBranchAddress("muPFNeuIso", &_muPFNeuIso);
   tree->SetBranchAddress("muPFPUIso", &_muPFPUIso);

}