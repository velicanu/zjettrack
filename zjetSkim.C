#include "ggTree.h"
#include "jetTree.h"
#include "trackTree.h"
#include "SkimTree.h"
#include "L2L3ResidualWFits.h"
#include "getTrkCorr.h"

//! cuts defined here , {veto, loose, medium, tight}
float cuts_barrel_pbpb_eleSigmaIEtaIEta_2012[4] = {0.0111,0.0108,0.0106,0.0097};
float cuts_barrel_pbpb_eleEoverPInv[4] = {0.2805,0.1240,0.0148,0.0148};
float cuts_barrel_pbpb_eledEtaAtVtx[4] = {0.0158,0.0157,0.0149,0.0147};
float cuts_barrel_pbpb_eledPhiAtVtx[4] = {0.1572,0.0544,0.0340,0.0257};
float cuts_barrel_pbpb_eleHoverE[4]    = {0.0885,0.0359,0.0333,0.0325};
float cuts_barrel_pbpb_eleD0[4]        = {0.0522,0.0260,0.0160,0.0159};
float cuts_barrel_pbpb_eleDz[4]        = {0.1300,0.1267,0.0694,0.0653};
float cuts_barrel_pbpb_eleMissHits[4]  = {1,1,1,1};

float cuts_endcap_pbpb_eleSigmaIEtaIEta_2012[4] = {0.03488 , 0.0312 , 0.0299 , 0.0294};
float cuts_endcap_pbpb_eleEoverPInv[4] = {0.1867 , 0.1628 , 0.1262 , 0.0596};
float cuts_endcap_pbpb_eledEtaAtVtx[4] = {0.0171 , 0.0158 , 0.0086 , 0.0085};
float cuts_endcap_pbpb_eledPhiAtVtx[4] = {0.3554 , 0.1907 , 0.1224 , 0.0350};
float cuts_endcap_pbpb_eleHoverE[4]    = {0.1228 , 0.1014 , 0.0927 , 0.0377};
float cuts_endcap_pbpb_eleD0[4]        = {0.1909 , 0.1060 , 0.0872 , 0.0800};
float cuts_endcap_pbpb_eleDz[4]        = {0.2641 , 0.2238 , 0.2032 , 0.1852};
float cuts_endcap_pbpb_eleMissHits[4]  = {1,1,1,1};

float cuts_barrel_pp_eleSigmaIEtaIEta_2012[4] = {0.0114,0.0103,0.0101,0.0101};
float cuts_barrel_pp_eleEoverPInv[4] = {0.207,0.102,0.0174,0.012};
float cuts_barrel_pp_eledEtaAtVtx[4] = {0.0152,0.0105,0.0103,0.00926};
float cuts_barrel_pp_eledPhiAtVtx[4] = {0.216,0.115,0.0336,0.0336};
float cuts_barrel_pp_eleHoverE[4]    = {0.181,0.104,0.0876,0.0597};
float cuts_barrel_pp_eleD0[4]        = {0.0564,0.0261,0.0118,0.0111};
float cuts_barrel_pp_eleDz[4]        = {0.472,0.41,0.373,0.0466};
float cuts_barrel_pp_eleMissHits[4]  = {2,2,2,2};

float cuts_endcap_pp_eleSigmaIEtaIEta_2012[4] = {0.0352,0.0301,0.0283,0.0279};
float cuts_endcap_pp_eleEoverPInv[4] = {0.174,0.126,0.0898,0.00999};
float cuts_endcap_pp_eledEtaAtVtx[4] = {0.0113,0.00814,0.00733,0.00724};
float cuts_endcap_pp_eledPhiAtVtx[4] = {0.237,0.182,0.114,0.0918};
float cuts_endcap_pp_eleHoverE[4]    = {0.116,0.0897,0.0678,0.0615};
float cuts_endcap_pp_eleD0[4]        = {0.222,0.118,0.0739,0.0351};
float cuts_endcap_pp_eleDz[4]        = {0.921,0.822,0.602,0.417};
float cuts_endcap_pp_eleMissHits[4]  = {3,1,1,1};
int cut_type_pbpb = 1;
int cut_type_pp = 2;


bool goodJet(int i) {
  if(	_neutralSum[i]/rawpt[i] < 0.9
      && _chargedSum[i]/rawpt[i] > 0.01
      && chargedN[i]+photonN[i]+neutralN[i]+eN[i]+muN[i] > 0
      && chargedMax[i]/rawpt[i] < 0.99
      && photonSum[i]/rawpt[i] < 0.99
      && _eSum[i]/rawpt[i] < 0.99
      ) return true;
  else return false;
}


float getTrkWeight(TrkCorr * trkCorr, int itrk, int hiBin)
{
  float rmin = 999;
  for(int k = 0; k<nref; k++)
  {
    if(jtpt[k]<50) break;
    if(!goodJet(k)) continue;
    if(TMath::Abs(jteta[k]>2)) continue;//jet quality cut
    float R = TMath::Power(jteta[k]-trkEta_[itrk],2)+TMath::Power(TMath::ACos(TMath::Cos(jtphi[k]-trkPhi_[itrk])),2);
    if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
  }
  return trkCorr->getTrkCorr(trkPt_[itrk],trkEta_[itrk],trkPhi_[itrk],hiBin,rmin);
}


bool goodMuon(int imu) {
  if(_muChi2NDF->at(imu)<10
      && _muMuonHits->at(imu)>0
      && _muStations->at(imu)>1
      && _muTrkLayers->at(imu)>5
      && _muPixelHits->at(imu)>0
      && fabs(_muInnerD0->at(imu))<0.2
      && fabs(_muInnerDz->at(imu))<0.5
      ) return true;
  else return false;
}

//! Electron can be pp or PbPb which have different cuts
bool goodElectron(int i, bool is_pp) {
  if(is_pp)
  {
    if(fabs(_eleSCEta->at(i))<1.479) {
      if(  fabs(_eledEtaAtVtx->at(i))<cuts_barrel_pp_eledEtaAtVtx[cut_type_pp]
          && fabs(_eledPhiAtVtx->at(i))<cuts_barrel_pp_eledPhiAtVtx[cut_type_pp]
          && _eleSigmaIEtaIEta->at(i)<cuts_barrel_pp_eleSigmaIEtaIEta_2012[cut_type_pp]
          && _eleHoverE->at(i)<cuts_barrel_pp_eleHoverE[cut_type_pp]
          && fabs(_eleD0->at(i))<cuts_barrel_pp_eleD0[cut_type_pp]
          && fabs(_eleDz->at(i))<cuts_barrel_pp_eleDz[cut_type_pp]
          && fabs(_eleEoverPInv->at(i))<cuts_barrel_pp_eleEoverPInv[cut_type_pp]
          && _eleMissHits->at(i) <= cuts_barrel_pp_eleMissHits[cut_type_pp]
          ) return true;
      else return false;
    }
    if(fabs(_eleSCEta->at(i))>1.479) {
      if(  fabs(_eledEtaAtVtx->at(i))<cuts_endcap_pp_eledEtaAtVtx[cut_type_pp]
          && fabs(_eledPhiAtVtx->at(i))<cuts_endcap_pp_eledPhiAtVtx[cut_type_pp]
          && _eleSigmaIEtaIEta->at(i)<cuts_endcap_pp_eleSigmaIEtaIEta_2012[cut_type_pp]
          && _eleHoverE->at(i)<cuts_endcap_pp_eleHoverE[cut_type_pp]
          && fabs(_eleD0->at(i))<cuts_endcap_pp_eleD0[cut_type_pp]
          && fabs(_eleDz->at(i))<cuts_endcap_pp_eleDz[cut_type_pp]
          && fabs(_eleEoverPInv->at(i))<cuts_endcap_pp_eleEoverPInv[cut_type_pp]
          && _eleMissHits->at(i) <= cuts_endcap_pp_eleMissHits[cut_type_pp]
          ) return true;
      else return false;
    }
  } else {
    if(fabs(_eleSCEta->at(i))<1.479) {
      if(  fabs(_eledEtaAtVtx->at(i))<cuts_barrel_pbpb_eledEtaAtVtx[cut_type_pp]
          && fabs(_eledPhiAtVtx->at(i))<cuts_barrel_pbpb_eledPhiAtVtx[cut_type_pp]
          && _eleSigmaIEtaIEta->at(i)<cuts_barrel_pbpb_eleSigmaIEtaIEta_2012[cut_type_pp]
          && _eleHoverE->at(i)<cuts_barrel_pbpb_eleHoverE[cut_type_pp]
          && fabs(_eleD0->at(i))<cuts_barrel_pbpb_eleD0[cut_type_pp]
          && fabs(_eleDz->at(i))<cuts_barrel_pbpb_eleDz[cut_type_pp]
          && fabs(_eleEoverPInv->at(i))<cuts_barrel_pbpb_eleEoverPInv[cut_type_pp]
          && _eleMissHits->at(i) <= cuts_barrel_pbpb_eleMissHits[cut_type_pp]
          ) return true;
      else return false;
    }
    if(fabs(_eleSCEta->at(i))>1.479) {
      if(  fabs(_eledEtaAtVtx->at(i))<cuts_endcap_pbpb_eledEtaAtVtx[cut_type_pp]
          && fabs(_eledPhiAtVtx->at(i))<cuts_endcap_pbpb_eledPhiAtVtx[cut_type_pp]
          && _eleSigmaIEtaIEta->at(i)<cuts_endcap_pbpb_eleSigmaIEtaIEta_2012[cut_type_pp]
          && _eleHoverE->at(i)<cuts_endcap_pbpb_eleHoverE[cut_type_pp]
          && fabs(_eleD0->at(i))<cuts_endcap_pbpb_eleD0[cut_type_pp]
          && fabs(_eleDz->at(i))<cuts_endcap_pbpb_eleDz[cut_type_pp]
          && fabs(_eleEoverPInv->at(i))<cuts_endcap_pbpb_eleEoverPInv[cut_type_pp]
          && _eleMissHits->at(i) <= cuts_endcap_pbpb_eleMissHits[cut_type_pp]
          ) return true;
      else return false;
    }
  }

  return false;
}


void gammajetSkim(TString infilename="HiForest.root", TString outfilename="Zevents.root", string jetname="ak4PFJetAnalyzer", int i_is_pp = 0 ) {

  bool is_pp = (i_is_pp == 1) ;
  TrkCorr* trkCorr;
  if(is_pp) trkCorr = new TrkCorr("TrkCorr_Mar15_Iterative_PbPb/");
  else trkCorr = new TrkCorr("TrkCorr_Mar15_Iterative_pp/");
  L2L3Residual * jetcorr = new L2L3Residual(3);
  TFile *fin = TFile::Open(infilename);

  TFile *fout = new TFile(outfilename,"recreate");

  TH1::SetDefaultSumw2();


  int Ztype; //type 1 muon, type 2 electron
  float Zmass, Zpt, Zeta, Zrapidity, Zphi;
  float Zlepton1Pt, Zlepton2Pt, Zlepton1Eta, Zlepton2Eta, Zlepton1Phi, Zlepton2Phi;
  float weight = 0 , vz = -99;
  int hiBin = -99;
  int Zcharge;
  float leptonptcut = 20;


  int njet;
  float jetpt[200], jeteta[200], jetphi[200];
  float chargedSum[200], neutralSum[200], eSum[200];
  int jetID[200], subid[200];

  const int nPtBins = 13;
  const double PtBins[nPtBins+1]={0,2.5,5.0,7.5,10.0,12.5,15.0,20,30,40,50,70,100,150};
  const int nYBins = 13;

  TH1F *mmMass = new TH1F("mmMass",";m^{#mu#mu} (GeV/c^{2});Events",30,60,120);
  TH1F *mmMassSS = new TH1F("mmMassSS",";m^{#mu#mu} (GeV/c^{2});Events",30,60,120);

  TH1F *mmPt = new TH1F("mmPt",";p_{T}^{#mu#mu} (GeV/c);dN/dp_{T}",nPtBins,PtBins);
  TH1F *mmY = new TH1F("mmY",";y^{#mu#mu};Events",nYBins,-2.6,2.6);
  TH1F *mmPhi = new TH1F("mmPhi",";#phi^{#mu#mu};Events",20,-pi,pi);

  TH1F *eeMass = new TH1F("eeMass",";m^{ee} (GeV/c^{2});Events",30,60,120);
  TH1F *eeMassSS = new TH1F("eeMassSS",";m^{ee} (GeV/c^{2});Events",30,60,120);

  TH1F *eePt = new TH1F("eePt",";p_{T}^{ee} (GeV/c);dN/dp_{T}",nPtBins,PtBins);
  TH1F *eeY = new TH1F("eeY",";y^{ee};Events",nYBins,-2.6,2.6);
  TH1F *eePhi = new TH1F("eePhi",";#phi^{ee};Events",20,-pi,pi);

  int nMuPair = 0;
  int nElPair = 0;


  Int_t           nTrk;
  Int_t           run;
  Int_t           event;
  Int_t           lumis;
  Float_t         trkPt[100000];   //[nTrk]
  Float_t         trkPtError[100000];   //[nTrk]
  UChar_t         trkNHit[100000];   //[nTrk]
  UChar_t         trkNlayer[100000];   //[nTrk]
  Float_t         trkEta[100000];   //[nTrk]
  Float_t         trkPhi[100000];   //[nTrk]
  Int_t           trkCharge[100000];   //[nTrk]
  UChar_t         trkNVtx[100000];   //[nTrk]
  Bool_t          highPurity[100000];   //[nTrk]
  Bool_t          tight[100000];   //[nTrk]
  Bool_t          loose[100000];   //[nTrk]
  Float_t         trkChi2[100000];   //[nTrk]
  UChar_t         trkNdof[100000];   //[nTrk]
  Float_t         trkDxy1[100000];   //[nTrk]
  Float_t         trkDxyError1[100000];   //[nTrk]
  Float_t         trkDz1[100000];   //[nTrk]
  Float_t         trkDzError1[100000];   //[nTrk]
  Bool_t          trkFake[100000];   //[nTrk]
  UChar_t         trkAlgo[100000];   //[nTrk]
  UChar_t         trkOriginalAlgo[100000];   //[nTrk]
  Float_t         trkMVA[100000];   //[nTrk]
  Int_t           pfType[100000];   //[nTrk]
  Float_t         pfCandPt[100000];   //[nTrk]
  Float_t         pfEcal[100000];   //[nTrk]
  Float_t         pfHcal[100000];   //[nTrk]
  Float_t         trkWeight[100000];   //[nTrk]
  Int_t    nPho;
  Float_t  phoE[100];   //_nPho
  Float_t  phoEt[100];   //_nPho
  Float_t  phoEta[100];   //_nPho
  Float_t  phoPhi[100];   //_nPho
  Float_t  phoSCE[100];   //_nPho
  Float_t  phoSCRawE[100];   //_nPho
  Float_t  phoESEn[100];   //_nPho
  Float_t  phoSCEta[100];   //_nPho
  Float_t  phoSCPhi[100];   //_nPho
  Float_t  phoSCEtaWidth[100];   //_nPho
  Float_t  phoSCPhiWidth[100];   //_nPho
  Float_t  phoSCBrem[100];   //_nPho
  Int_t    phohasPixelSeed[100];   //_nPho
  Float_t  phoR9[100];   //_nPho
  Float_t  phoHoverE[100];   //_nPho
  Float_t  phoSigmaIEtaIEta[100];   //_nPho
  Float_t  pho_ecalClusterIsoR2[100];   //_nPho
  Float_t  pho_ecalClusterIsoR3[100];   //_nPho
  Float_t  pho_ecalClusterIsoR4[100];   //_nPho
  Float_t  pho_ecalClusterIsoR5[100];   //_nPho
  Float_t  pho_hcalRechitIsoR1[100];   //_nPho
  Float_t  pho_hcalRechitIsoR2[100];   //_nPho
  Float_t  pho_hcalRechitIsoR3[100];   //_nPho
  Float_t  pho_hcalRechitIsoR4[100];   //_nPho
  Float_t  pho_hcalRechitIsoR5[100];   //_nPho
  Float_t  pho_trackIsoR1PtCut20[100];   //_nPho
  Float_t  pho_trackIsoR2PtCut20[100];   //_nPho
  Float_t  pho_trackIsoR3PtCut20[100];   //_nPho
  Float_t  pho_trackIsoR4PtCut20[100];   //_nPho
  Float_t  pho_trackIsoR5PtCut20[100];   //_nPho
  Float_t  pho_swissCrx[100];   //_nPho
  Float_t  pho_seedTime[100];   //_nPho
  Float_t  pfcIso1[100];   //_nPho
  Float_t  pfcIso2[100];   //_nPho
  Float_t  pfcIso3[100];   //_nPho
  Float_t  pfcIso4[100];   //_nPho
  Float_t  pfcIso5[100];   //_nPho
  Float_t  pfpIso1[100];   //_nPho
  Float_t  pfpIso2[100];   //_nPho
  Float_t  pfpIso3[100];   //_nPho
  Float_t  pfpIso4[100];   //_nPho
  Float_t  pfpIso5[100];   //_nPho
  Float_t  pfnIso1[100];   //_nPho
  Float_t  pfnIso2[100];   //_nPho
  Float_t  pfnIso3[100];   //_nPho
  Float_t  pfnIso4[100];   //_nPho
  Float_t  pfnIso5[100];   //_nPho
	Int_t    nMC;
	Int_t    mcPID[100];
	Int_t    mcStatus[100];
	Float_t  mcPt[100];
	Float_t  mcEta[100];
	Float_t  mcPhi[100];
	Float_t  mcMomPt[100];
	Float_t  mcMomEta[100];
	Float_t  mcMomPhi[100];
	Int_t    mcMomPID[100];

  TTree *ztree = new TTree("ztree","Jet track tree");
  ztree->Branch("run",	&run,	"run/I");
  ztree->Branch("event",	&event,	"event/I");
  ztree->Branch("lumis",	&lumis,	"lumis/I");
  ztree->Branch("hiBin", &hiBin, "hiBin/I");
  ztree->Branch("vz", &vz, "vz/F");
  ztree->Branch("Ztype",	&Ztype,	"Ztype/I");
  ztree->Branch("Zmass",	&Zmass,	"Zmass/F");
  ztree->Branch("Zpt",	&Zpt,	"Zpt/F");
  ztree->Branch("Zeta",	&Zeta,	"Zeta/F");
  ztree->Branch("Zphi",	&Zphi,	"Zphi/F");
  ztree->Branch("Zrapidity",	&Zrapidity,	"Zrapidity/F");
  ztree->Branch("Zlepton1Pt",	&Zlepton1Pt,	"Zlepton1Pt/F");
  ztree->Branch("Zlepton2Pt",	&Zlepton2Pt,	"Zlepton2Pt/F");
  ztree->Branch("Zlepton1Eta",	&Zlepton1Eta,	"Zlepton1Eta/F");
  ztree->Branch("Zlepton2Eta",	&Zlepton2Eta,	"Zlepton2Eta/F");
  ztree->Branch("Zlepton1Phi",	&Zlepton1Phi,	"Zlepton1Phi/F");
  ztree->Branch("Zlepton2Phi",	&Zlepton2Phi,	"Zlepton2Phi/F");
  ztree->Branch("Zcharge",	&Zcharge,	"Zcharge/I");
  ztree->Branch("njet",	&njet,	"njet/I");
  ztree->Branch("jetpt",	&jetpt,	"jetpt[njet]/F");
  ztree->Branch("jeteta",	&jeteta,	"jeteta[njet]/F");
  ztree->Branch("jetphi",	&jetphi,	"jetphi[njet]/F");
  ztree->Branch("jetID",	&jetID,	"jetID[njet]/I");
  ztree->Branch("subid",	&subid,	"subid[njet]/I");
  ztree->Branch("chargedSum",	&chargedSum,	"chargedSum[njet]/F");
  ztree->Branch("neutralSum",	&neutralSum,	"neutralSum[njet]/F");
  ztree->Branch("eSum",	&eSum,	"eSum[njet]/F");
  ztree->Branch("nTrk",	&nTrk,	"nTrk/I");
  ztree->Branch("trkPt",	&trkPt,	"trkPt[nTrk]/F");
  ztree->Branch("trkPtError", &trkPtError,"trkPtError[nTrk]/F");
  ztree->Branch("trkNHit", &trkNHit,"trkNHit[nTrk]/b");
  ztree->Branch("trkNlayer", &trkNlayer,"trkNlayer[nTrk]/b");
  ztree->Branch("trkEta", &trkEta,"trkEta[nTrk]/F");
  ztree->Branch("trkPhi", &trkPhi,"trkPhi[nTrk]/F");
  ztree->Branch("trkCharge", &trkCharge,"trkCharge[nTrk]/I");
  ztree->Branch("trkNVtx", &trkNVtx,"trkNVtx[nTrk]/b");
  ztree->Branch("highPurity", &highPurity,"highPurity[nTrk]/O");
  ztree->Branch("tight", &tight,"tight[nTrk]/O");
  ztree->Branch("loose", &loose,"loose[nTrk]/O");
  ztree->Branch("trkChi2", &trkChi2,"trkChi2[nTrk]/F");
  ztree->Branch("trkNdof", &trkNdof,"trkNdof[nTrk]/b");
  ztree->Branch("trkDxy1", &trkDxy1,"trkDxy1[nTrk]/F");
  ztree->Branch("trkDxyError1", &trkDxyError1,"trkDxyError1[nTrk]/F");
  ztree->Branch("trkDz1", &trkDz1,"trkDz1[nTrk]/F");
  ztree->Branch("trkDzError1", &trkDzError1,"trkDzError1[nTrk]/F");
  ztree->Branch("trkFake", &trkFake,"trkFake[nTrk]/O");
  ztree->Branch("trkAlgo", &trkAlgo,"trkAlgo[nTrk]/b");
  ztree->Branch("trkOriginalAlgo", &trkOriginalAlgo,"trkOriginalAlgo[nTrk]/b");
  ztree->Branch("trkMVA", &trkMVA,"trkMVA[nTrk]/F");
  ztree->Branch("pfType", &pfType,"pfType[nTrk]/I");
  ztree->Branch("pfCandPt", &pfCandPt,"pfCandPt[nTrk]/F");
  ztree->Branch("pfEcal", &pfEcal,"pfEcal[nTrk]/F");
  ztree->Branch("pfHcal", &pfHcal,"pfHcal[nTrk]/F");
  ztree->Branch("trkWeight", &trkWeight,"trkWeight[nTrk]/F");
  ztree->Branch("weight", &weight,"weight/F");

  ztree->Branch("nPho",&nPho,"nPho/I");
  ztree->Branch("phoE",&phoE,"phoE[nPho]/F");
  ztree->Branch("phoEt",&phoEt,"phoEt[nPho]/F");
  ztree->Branch("phoEta",&phoEta,"phoEta[nPho]/F");
  ztree->Branch("phoPhi",&phoPhi,"phoPhi[nPho]/F");
  ztree->Branch("phoSCE",&phoSCE,"phoSCE[nPho]/F");
  ztree->Branch("phoSCRawE",&phoSCRawE,"phoSCRawE[nPho]/F");
  ztree->Branch("phoESEn",&phoESEn,"phoESEn[nPho]/F");
  ztree->Branch("phoSCEta",&phoSCEta,"phoSCEta[nPho]/F");
  ztree->Branch("phoSCPhi",&phoSCPhi,"phoSCPhi[nPho]/F");
  ztree->Branch("phoSCEtaWidth",&phoSCEtaWidth,"phoSCEtaWidth[nPho]/F");
  ztree->Branch("phoSCPhiWidth",&phoSCPhiWidth,"phoSCPhiWidth[nPho]/F");
  ztree->Branch("phoSCBrem",&phoSCBrem,"phoSCBrem[nPho]/F");
  ztree->Branch("phohasPixelSeed",&phohasPixelSeed,"phohasPixelSeed[nPho]/I");
  ztree->Branch("phoR9",&phoR9,"phoR9[nPho]/F");
  ztree->Branch("phoHoverE",&phoHoverE,"phoHoverE[nPho]/F");
  ztree->Branch("phoSigmaIEtaIEta",&phoSigmaIEtaIEta,"phoSigmaIEtaIEta[nPho]/F");
  ztree->Branch("pho_ecalClusterIsoR2",&pho_ecalClusterIsoR2,"pho_ecalClusterIsoR2[nPho]/F");
  ztree->Branch("pho_ecalClusterIsoR3",&pho_ecalClusterIsoR3,"pho_ecalClusterIsoR3[nPho]/F");
  ztree->Branch("pho_ecalClusterIsoR4",&pho_ecalClusterIsoR4,"pho_ecalClusterIsoR4[nPho]/F");
  ztree->Branch("pho_ecalClusterIsoR5",&pho_ecalClusterIsoR5,"pho_ecalClusterIsoR5[nPho]/F");
  ztree->Branch("pho_hcalRechitIsoR1",&pho_hcalRechitIsoR1,"pho_hcalRechitIsoR1[nPho]/F");
  ztree->Branch("pho_hcalRechitIsoR2",&pho_hcalRechitIsoR2,"pho_hcalRechitIsoR2[nPho]/F");
  ztree->Branch("pho_hcalRechitIsoR3",&pho_hcalRechitIsoR3,"pho_hcalRechitIsoR3[nPho]/F");
  ztree->Branch("pho_hcalRechitIsoR4",&pho_hcalRechitIsoR4,"pho_hcalRechitIsoR4[nPho]/F");
  ztree->Branch("pho_hcalRechitIsoR5",&pho_hcalRechitIsoR5,"pho_hcalRechitIsoR5[nPho]/F");
  ztree->Branch("pho_trackIsoR1PtCut20",&pho_trackIsoR1PtCut20,"pho_trackIsoR1PtCut20[nPho]/F");
  ztree->Branch("pho_trackIsoR2PtCut20",&pho_trackIsoR2PtCut20,"pho_trackIsoR2PtCut20[nPho]/F");
  ztree->Branch("pho_trackIsoR3PtCut20",&pho_trackIsoR3PtCut20,"pho_trackIsoR3PtCut20[nPho]/F");
  ztree->Branch("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20,"pho_trackIsoR4PtCut20[nPho]/F");
  ztree->Branch("pho_trackIsoR5PtCut20",&pho_trackIsoR5PtCut20,"pho_trackIsoR5PtCut20[nPho]/F");
  ztree->Branch("pho_swissCrx",&pho_swissCrx,"pho_swissCrx[nPho]/F");
  ztree->Branch("pho_seedTime",&pho_seedTime,"pho_seedTime[nPho]/F");
  ztree->Branch("pfcIso1",&pfcIso1,"pfcIso1[nPho]/F");
  ztree->Branch("pfcIso2",&pfcIso2,"pfcIso2[nPho]/F");
  ztree->Branch("pfcIso3",&pfcIso3,"pfcIso3[nPho]/F");
  ztree->Branch("pfcIso4",&pfcIso4,"pfcIso4[nPho]/F");
  ztree->Branch("pfcIso5",&pfcIso5,"pfcIso5[nPho]/F");
  ztree->Branch("pfpIso1",&pfpIso1,"pfpIso1[nPho]/F");
  ztree->Branch("pfpIso2",&pfpIso2,"pfpIso2[nPho]/F");
  ztree->Branch("pfpIso3",&pfpIso3,"pfpIso3[nPho]/F");
  ztree->Branch("pfpIso4",&pfpIso4,"pfpIso4[nPho]/F");
  ztree->Branch("pfpIso5",&pfpIso5,"pfpIso5[nPho]/F");
  ztree->Branch("pfnIso1",&pfnIso1,"pfnIso1[nPho]/F");
  ztree->Branch("pfnIso2",&pfnIso2,"pfnIso2[nPho]/F");
  ztree->Branch("pfnIso3",&pfnIso3,"pfnIso3[nPho]/F");
  ztree->Branch("pfnIso4",&pfnIso4,"pfnIso4[nPho]/F");
  ztree->Branch("pfnIso5",&pfnIso5,"pfnIso5[nPho]/F");
	ztree->Branch("nMC", &nMC,"nMC/I");
	ztree->Branch("mcPID", &mcPID,"mcPID[nMC]/I");
	ztree->Branch("mcStatus", &mcStatus,"mcStatus[nMC]/I");
	ztree->Branch("mcPt", &mcPt,"mcPt[nMC]/F");
	ztree->Branch("mcEta", &mcEta,"mcEta[nMC]/F");
	ztree->Branch("mcPhi", &mcPhi,"mcPhi[nMC]/F");
	ztree->Branch("mcMomPt", &mcMomPt,"mcMomPt[nMC]/F");
	ztree->Branch("mcMomEta", &mcMomEta,"mcMomEta[nMC]/F");
	ztree->Branch("mcMomPhi", &mcMomPhi,"mcMomPhi[nMC]/F");
	ztree->Branch("mcMomPID", &mcMomPID,"mcMomPID[nMC]/I");
  // ztree->Branch("nMu",&nMu,"nMu/I");
  // ztree->Branch("muPt",&muPt,"muPt[nMu]/F");
  // ztree->Branch("muEta",&muEta,"[nMu]/F");
  // ztree->Branch("muPhi",&muPhi,"muPhi[nMu]/F");
  // ztree->Branch("muCharge",&muCharge,"muCharge[nMu]/I");
  // ztree->Branch("muType",&muType,"muType[nMu]/I");
  // ztree->Branch("muIsGood",&muIsGood,"muIsGood[nMu]/I");
  // ztree->Branch("muD0",&muD0,"muD0[nMu]/F");
  // ztree->Branch("muDz",&muDz,"muDz[nMu]/F");
  // ztree->Branch("muChi2NDF",&muChi2NDF,"muChi2NDF[nMu]/F");
  // ztree->Branch("muInnerD0",&muInnerD0,"muInnerD0[nMu]/F");
  // ztree->Branch("muInnerDz",&muInnerDz,"muInnerDz[nMu]/F");
  // ztree->Branch("muTrkLayers",&muTrkLayers,"muTrkLayers[nMu]/I");
  // ztree->Branch("muPixelLayers",&muPixelLayers,"muPixelLayers[nMu]/I");
  // ztree->Branch("muPixelHits",&muPixelHits,"muPixelHits[nMu]/I");
  // ztree->Branch("muMuonHits",&muMuonHits,"muMuonHits[nMu]/I");
  // ztree->Branch("muTrkQuality",&muTrkQuality,"muTrkQuality[nMu]/I");
  // ztree->Branch("muStations",&muStations,"muStations[nMu]/I");
  // ztree->Branch("muIsoTrk",&muIsoTrk,"muIsoTrk[nMu]/F");
  // ztree->Branch("muPFChIso",&muPFChIso,"muPFChIso[nMu]/F");
  // ztree->Branch("muPFPhoIso",&muPFPhoIso,"muPFPhoIso[nMu]/F");
  // ztree->Branch("muPFNeuIso",&muPFNeuIso,"muPFNeuIso[nMu]/F");
  // ztree->Branch("muPFPUIso",&muPFPUIso,"muPFPUIso[nMu]/F");


  // const int nPtBins = 13;
  // const double PtBins[nPtBins+1]={0,2.5,5.0,7.5,10.0,12.5,15.0,20,30,40,50,70,100,150};
  // const int nYBins = 13;

  // TH1F *mmMass = new TH1F("mmMass",";m^{#mu#mu} (GeV/c^{2});Events",30,60,120);
  // TH1F *mmMassSS = new TH1F("mmMassSS",";m^{#mu#mu} (GeV/c^{2});Events",30,60,120);

  // TH1F *mmPt = new TH1F("mmPt",";p_{T}^{#mu#mu} (GeV/c);dN/dp_{T}",nPtBins,PtBins);
  // TH1F *mmY = new TH1F("mmY",";y^{#mu#mu};Events",nYBins,-2.6,2.6);
  // TH1F *mmPhi = new TH1F("mmPhi",";#phi^{#mu#mu};Events",20,-pi,pi);

  // TH1F *eeMass = new TH1F("eeMass",";m^{ee} (GeV/c^{2});Events",30,60,120);
  // TH1F *eeMassSS = new TH1F("eeMassSS",";m^{ee} (GeV/c^{2});Events",30,60,120);

  // TH1F *eePt = new TH1F("eePt",";p_{T}^{ee} (GeV/c);dN/dp_{T}",nPtBins,PtBins);
  // TH1F *eeY = new TH1F("eeY",";y^{ee};Events",nYBins,-2.6,2.6);
  // TH1F *eePhi = new TH1F("eePhi",";#phi^{ee};Events",20,-pi,pi);

  // int nMuPair = 0;
  // int nElPair = 0;


  TTree *inggTree = (TTree*)fin->Get("ggHiNtuplizer/EventTree");
  if(!inggTree){
    cout<<"Could not access gg tree!"<<endl;
    return;
  }
  initggTree(inggTree);
  // akPu4CaloJetAnalyzer
  TTree *injetTree = (TTree*)fin->Get(Form("%s/t",jetname.data()));
  // TTree *injetTree = (TTree*)fin->Get("ak3PFJetAnalyzer/t");
  // if(!injetTree) injetTree = (TTree*)fin->Get("akPu3PFJetAnalyzer/t");
  if(!injetTree){
    cout<<"Could not access jet tree!"<<endl;
    return;
  }
  initjetTree(injetTree);

  TTree *evttree = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
  if(!evttree){
    cout<<"Could not access event tree!"<<endl;
    return;
  }
  evttree->SetBranchAddress("weight", &weight);
  evttree->SetBranchAddress("hiBin", &hiBin);
  evttree->SetBranchAddress("vz", &vz);

  TTree * tracktree_                     = (TTree*) fin->Get("anaTrack/trackTree");
  if( tracktree_ == 0 ) tracktree_        = (TTree*) fin->Get("ppTrack/trackTree");
  if(!tracktree_){
    cout<<"Could not access track tree!"<<endl;
    return;
  }
  initTrackTree(tracktree_);

  TTree * skimTree                     = (TTree*) fin->Get("skimanalysis/HltTree");
  if( skimTree == 0 )
  {
    cout<<"Could not access skim tree!"<<endl;
    return;
  }
  initSkimTree(skimTree);



  // int nEv = inggTree->GetEntries();
  int nEv = evttree->GetEntries();

  for (int j=0; j<nEv; j++) {
    Zlepton1Pt=-99;
    Zlepton2Pt=-99;
    Zlepton1Eta=-99;
    Zlepton2Eta=-99;
    Zlepton1Phi=-99;
    Zlepton2Phi=-99;

    Zlepton1Pt=-99;
    Zlepton2Pt=-99;
    Zlepton1Eta=-99;
    Zlepton2Eta=-99;
    Zlepton1Phi=-99;
    Zlepton2Phi=-99;

    skimTree->GetEntry(j);
    // if(!(HBHENoiseFilterResultRun2Loose && pPAprimaryVertexFilter && pBeamScrapingFilter)) continue;
    evttree->GetEntry(j);
    if(fabs(vz)>15) continue;
    injetTree->GetEntry(j);
    if(j%10000 == 0) { cout << "Processing event: " << j << endl; }

    njet = 0;
    for(int ij=0; ij<nref; ij++) {
      // if(jtpt[ij]>1)
      // if(goodJet(ij))
      if(jtpt[ij]>1 && goodJet(ij))
      {
        jetpt[njet] = jetcorr->get_corrected_pt(jtpt[ij],jteta[ij]);
        jetpt[njet] = jtpt[ij];
        jeteta[njet] = jteta[ij];
        jetphi[njet] = jtphi[ij];
        jetID[njet] = goodJet(ij);
        subid[njet] = _subid[ij];
        chargedSum[njet] = _chargedSum[ij];
        neutralSum[njet] = _neutralSum[ij];
        eSum[njet] = _eSum[ij];
        njet++;
      }
    } //end of jet loop
    // if(njet==0) continue;
    // cout<<njet<<endl;

    inggTree->GetEntry(j);
		int nmcphoton = 0;
    for(int imc = 0 ; imc < _nMC ; ++imc)
		{
			// if(abs((*_mcMomPID)[imc])!=22) continue; //only signal photons if mc
			nmcphoton++;
			mcPID[imc] = (*_mcPID)[imc];
			mcStatus[imc] = (*_mcStatus)[imc];
			mcPt[imc] = (*_mcPt)[imc];
			mcEta[imc] = (*_mcEta)[imc];
			mcPhi[imc] = (*_mcPhi)[imc];
			mcMomPt[imc] = (*_mcMomPt)[imc];
			mcMomEta[imc] = (*_mcMomEta)[imc];
			mcMomPhi[imc] = (*_mcMomPhi)[imc];
			mcMomPID[imc] = (*_mcMomPID)[imc];
		}
		nMC = nmcphoton;
		// if(!(_nMC==0 || nmcphoton!=0)) continue;

    int nphoton = 0;
		for(int ipho = 0 ; ipho < _nPho ; ++ipho)
    {

      // if((*_phoHoverE)[ipho]<0.1 && (*_pho_swissCrx)[ipho]<0.9 && abs((*_pho_seedTime)[ipho])<3.0 && ((*_pho_ecalClusterIsoR4)[ipho] + (*_pho_hcalRechitIsoR4)[ipho] + (*_pho_trackIsoR4PtCut20)[ipho]) < 1.0 && (*_phoSigmaIEtaIEta_2012)[ipho]<0.01 && (*_phoR9)[ipho]>0.3 && _phoEt->at(ipho)>40 ) //photon selection
      if(_phoEt->at(ipho)>40 ) //photon selection
      {
        float phopt = _phoEt->at(ipho);
        float phoptphopt = phopt*phopt;
        if(fabs(_pho_seedTime->at(ipho)) > 3.0) continue;
        if(fabs(_pho_swissCrx->at(ipho)) > 0.9) continue;
        if(!(_pfcIso4->at(ipho) < 0.76  && _pfnIso4->at(ipho) < (0.97 + 0.014*phopt + 0.000019*phoptphopt) && _pfpIso4->at(ipho) < (0.08 + 0.0053*phopt))) continue;
        if(_phoHoverE->at(ipho)>0.05) continue;
        if(_phoSigmaIEtaIEta->at(ipho) > 0.0100 ) continue;
        phoE[nphoton] = (*_phoE)[ipho];
        phoEt[nphoton] = (*_phoEt)[ipho];
        phoEta[nphoton] = (*_phoEta)[ipho];
        phoPhi[nphoton] = (*_phoPhi)[ipho];
        phoSCE[nphoton] = (*_phoSCE)[ipho];
        phoSCRawE[nphoton] = (*_phoSCRawE)[ipho];
        phoESEn[nphoton] = (*_phoESEn)[ipho];
        phoSCEta[nphoton] = (*_phoSCEta)[ipho];
        phoSCPhi[nphoton] = (*_phoSCPhi)[ipho];
        phoSCEtaWidth[nphoton] = (*_phoSCEtaWidth)[ipho];
        phoSCPhiWidth[nphoton] = (*_phoSCPhiWidth)[ipho];
        phoSCBrem[nphoton] = (*_phoSCBrem)[ipho];
        phohasPixelSeed[nphoton] = (*_phohasPixelSeed)[ipho];
        phoR9[nphoton] = (*_phoR9)[ipho];
        phoHoverE[nphoton] = (*_phoHoverE)[ipho];
        phoSigmaIEtaIEta[nphoton] = (*_phoSigmaIEtaIEta)[ipho];
        pho_ecalClusterIsoR2[nphoton] = (*_pho_ecalClusterIsoR2)[ipho];
        pho_ecalClusterIsoR3[nphoton] = (*_pho_ecalClusterIsoR3)[ipho];
        pho_ecalClusterIsoR4[nphoton] = (*_pho_ecalClusterIsoR4)[ipho];
        pho_ecalClusterIsoR5[nphoton] = (*_pho_ecalClusterIsoR5)[ipho];
        pho_hcalRechitIsoR1[nphoton] = (*_pho_hcalRechitIsoR1)[ipho];
        pho_hcalRechitIsoR2[nphoton] = (*_pho_hcalRechitIsoR2)[ipho];
        pho_hcalRechitIsoR3[nphoton] = (*_pho_hcalRechitIsoR3)[ipho];
        pho_hcalRechitIsoR4[nphoton] = (*_pho_hcalRechitIsoR4)[ipho];
        pho_hcalRechitIsoR5[nphoton] = (*_pho_hcalRechitIsoR5)[ipho];
        pho_trackIsoR1PtCut20[nphoton] = (*_pho_trackIsoR1PtCut20)[ipho];
        pho_trackIsoR2PtCut20[nphoton] = (*_pho_trackIsoR2PtCut20)[ipho];
        pho_trackIsoR3PtCut20[nphoton] = (*_pho_trackIsoR3PtCut20)[ipho];
        pho_trackIsoR4PtCut20[nphoton] = (*_pho_trackIsoR4PtCut20)[ipho];
        pho_trackIsoR5PtCut20[nphoton] = (*_pho_trackIsoR5PtCut20)[ipho];
        pho_swissCrx[nphoton] = (*_pho_swissCrx)[ipho];
        pho_seedTime[nphoton] = (*_pho_seedTime)[ipho];
        pfcIso1[nphoton] = (*_pfcIso1)[ipho];
        pfcIso2[nphoton] = (*_pfcIso2)[ipho];
        pfcIso3[nphoton] = (*_pfcIso3)[ipho];
        pfcIso4[nphoton] = (*_pfcIso4)[ipho];
        pfcIso5[nphoton] = (*_pfcIso5)[ipho];
        pfpIso1[nphoton] = (*_pfpIso1)[ipho];
        pfpIso2[nphoton] = (*_pfpIso2)[ipho];
        pfpIso3[nphoton] = (*_pfpIso3)[ipho];
        pfpIso4[nphoton] = (*_pfpIso4)[ipho];
        pfpIso5[nphoton] = (*_pfpIso5)[ipho];
        pfnIso1[nphoton] = (*_pfnIso1)[ipho];
        pfnIso2[nphoton] = (*_pfnIso2)[ipho];
        pfnIso3[nphoton] = (*_pfnIso3)[ipho];
        pfnIso4[nphoton] = (*_pfnIso4)[ipho];
        pfnIso5[nphoton] = (*_pfnIso5)[ipho];
        nphoton++;
      }
    }
    nPho = nphoton;
    // if(nphoton==0) continue;

    bool flagMu = 0; bool flagEle = 0;
    if(flagMu || flagEle) cout<<"nothing"<<endl;

    TLorentzVector muon1, muon2;
    TLorentzVector ele1, ele2;

    for(int i1 = 0; i1 < _nMu; i1++) {

      if(_muPt->at(i1)>(leptonptcut-10) && fabs(_muEta->at(i1))<2.4 && goodMuon(i1)) {

        for(int i2 = i1+1; i2 < _nMu; i2++) {

          if(_muPt->at(i2)>(leptonptcut-10) && fabs(_muEta->at(i2))<2.4 && goodMuon(i2)) {

            muon1.SetPtEtaPhiM(_muPt->at(i1), _muEta->at(i1), _muPhi->at(i1), 0.105658);
            muon2.SetPtEtaPhiM(_muPt->at(i2), _muEta->at(i2), _muPhi->at(i2), 0.105658);
            TLorentzVector pair = muon1 + muon2;

            if(pair.M()>60 && pair.M()<120) {

              Ztype = 1;
              Zmass = pair.M();
              Zpt = pair.Pt();
              Zrapidity = pair.Rapidity();
              Zeta = pair.Eta();
              Zphi = pair.Phi();
              Zcharge = _muCharge->at(i1) + _muCharge->at(i2);
              flagMu = 1;
              Zlepton1Pt = muon1.Pt();
              Zlepton2Pt = muon2.Pt();
              Zlepton1Eta = muon1.Eta();
              Zlepton2Eta = muon2.Eta();
              Zlepton1Phi = muon1.Phi();
              Zlepton2Phi = muon2.Phi();
              if(_muCharge->at(i1) != _muCharge->at(i2)) {
                mmMass->Fill(pair.M());
                mmPt->Fill(pair.Pt());
                mmY->Fill(pair.Rapidity());
                mmPhi->Fill(pair.Phi());
                nMuPair++;
              }
              else {
                mmMassSS->Fill(pair.M());
              }

            }
          }
        }
      }
    } //end of muon loop

    for(int i1 = 0; i1 < _nEle; i1++) {

      if(_elePt->at(i1)>leptonptcut && fabs(_eleSCEta->at(i1))<2.5 && goodElectron(i1,is_pp) && (fabs(_eleSCEta->at(i1))<1.4442 || fabs(_eleSCEta->at(i1))>1.566)) {

        for(int i2 = i1+1; i2 < _nEle; i2++) {

          if(_elePt->at(i2)>leptonptcut && fabs(_eleSCEta->at(i2))<2.5 && goodElectron(i2,is_pp) && (fabs(_eleSCEta->at(i2))<1.4442 || fabs(_eleSCEta->at(i2))>1.566)) {

            ele1.SetPtEtaPhiM(_elePt->at(i1), _eleEta->at(i1), _elePhi->at(i1), 0.000511);
            ele2.SetPtEtaPhiM(_elePt->at(i2), _eleEta->at(i2), _elePhi->at(i2), 0.000511);
            TLorentzVector pair = ele1 + ele2;

            if(pair.M()>60 && pair.M()<120) {

              Ztype = 2;
              Zmass = pair.M();
              Zpt = pair.Pt();
              Zrapidity = pair.Rapidity();
              Zeta = pair.Eta();
              Zphi = pair.Phi();
              Zcharge = _eleCharge->at(i1) + _eleCharge->at(i2);
              flagEle = 1;
              Zlepton1Pt = ele1.Pt();
              Zlepton2Pt = ele2.Pt();
              Zlepton1Eta = ele1.Eta();
              Zlepton2Eta = ele2.Eta();
              Zlepton1Phi = ele1.Phi();
              Zlepton2Phi = ele2.Phi();

              if(_eleCharge->at(i1) != _eleCharge->at(i2)) {
                eeMass->Fill(pair.M());
                eePt->Fill(pair.Pt());
                eeY->Fill(pair.Rapidity());
                eePhi->Fill(pair.Phi());
                nElPair++;
              }
              else {
                eeMassSS->Fill(pair.M());
              }

            }
          }
        }
      }
    } //end of electron loop

    if( flagEle==0 && flagMu==0 ) continue;


    tracktree_->GetEntry(j);

    int ntracks = 0;
    // std::cout<<nTrk_<<std::endl;
    for(int itrk = 0 ; itrk < nTrk_ ; ++itrk)
    {
      // if((trkMVA_[itrk]<0.5 && trkMVA_[itrk]!=-99) || (int)trkNHit_[itrk]<8 || trkPtError_[itrk]/trkPt_[itrk]>0.3 || fabs(trkDz1_[itrk])/trkDzError1_[itrk]>3 || fabs(trkDxy1_[itrk])/trkDxyError1_[itrk]>3) continue;
      // if((Zlepton1Pt!=-99&&sqrt(pow(Zlepton1Phi- trkPhi_[itrk],2) + pow(Zlepton1Eta- trkEta_[itrk],2))<0.006)) continue; // reject z leptons
      // if((Zlepton2Pt!=-99&&sqrt(pow(Zlepton2Phi- trkPhi_[itrk],2) + pow(Zlepton2Eta- trkEta_[itrk],2))<0.006)) continue; // reject z leptons
      // if(!highPurity_[itrk]) continue;
      // std::cout<<"here"<<std::endl;
      if(trkPt_[itrk]<1 || trkPt_[itrk]>300 || fabs(trkEta_[itrk])>2.4 ) continue;
      if(highPurity_[itrk]!=1) continue;
      if(trkPtError_[itrk]/trkPt_[itrk]>0.1 || TMath::Abs(trkDz1_[itrk]/trkDzError1_[itrk])>3 || TMath::Abs(trkDxy1_[itrk]/trkDxyError1_[itrk])>3) continue;
      if(trkChi2_[itrk]/(float)trkNdof_[itrk]/(float)trkNlayer_[itrk]>0.15) continue;
      if(trkNHit_[itrk]<11 && trkPt_[itrk]>0.7) continue;
      if((maxJetPt>50 && trkPt[itrk]>maxJetPt) || (maxJetPt<50 && trkPt[itrk]>50)) continue;

      float Et = (pfHcal[itrk]+pfEcal[itrk])/TMath::CosH(trkEta[itrk]);
      if(!(trkPt[itrk]<20 || (Et>0.5*trkPt[itrk]))) continue;

      float trkweight = 0;
      if(is_pp) trkweight = getTrkWeight(trkCorr,i,0);
      else trkweight = getTrkWeight(trkCorr,i,hiBin);
      trkPt[ntracks] = trkPt_[itrk];   //[nTrk]
      trkPtError[ntracks] = trkPtError_[itrk];   //[nTrk]
      trkNHit[ntracks] = trkNHit_[itrk];   //[nTrk]
      trkNlayer[ntracks] = trkNlayer_[itrk];   //[nTrk]
      trkEta[ntracks] = trkEta_[itrk];   //[nTrk]
      trkPhi[ntracks] = trkPhi_[itrk];   //[nTrk]
      trkCharge[ntracks] = trkCharge_[itrk];   //[nTrk]
      trkNVtx[ntracks] = trkNVtx_[itrk];   //[nTrk]
      highPurity[ntracks] = highPurity_[itrk];   //[nTrk]
      tight[ntracks] = tight_[itrk];   //[nTrk]
      loose[ntracks] = loose_[itrk];   //[nTrk]
      trkChi2[ntracks] = trkChi2_[itrk];   //[nTrk]
      trkNdof[ntracks] = trkNdof_[itrk];   //[nTrk]
      trkDxy1[ntracks] = trkDxy1_[itrk];   //[nTrk]
      trkDxyError1[ntracks] = trkDxyError1_[itrk];   //[nTrk]
      trkDz1[ntracks] = trkDz1_[itrk];   //[nTrk]
      trkDzError1[ntracks] = trkDzError1_[itrk];   //[nTrk]
      trkFake[ntracks] = trkFake_[itrk];   //[nTrk]
      trkAlgo[ntracks] = trkAlgo_[itrk];   //[nTrk]
      trkOriginalAlgo[ntracks] = trkOriginalAlgo_[itrk];   //[nTrk]
      trkMVA[ntracks] = trkMVA_[itrk];   //[nTrk]
      pfType[ntracks] = pfType_[itrk];   //[nTrk]
      pfCandPt[ntracks] = pfCandPt_[itrk];   //[nTrk]
      pfEcal[ntracks] = pfEcal_[itrk];   //[nTrk]
      pfHcal[ntracks] = pfHcal_[itrk];   //[nTrk]
      trkWeight[ntracks] = trkweight;
      ntracks++;
      //if((trkPt[itrk]-2*trkPtError[itrk])*TMath::CosH(trkEta[itrk])>15 && (trkPt[itrk]-2*trkPtError[itrk])*TMath::CosH(trkEta[itrk])>pfHcal[itrk]+pfEcal[itrk]) continue;} //Calo Matching
    }
    nTrk=ntracks;

    ztree->Fill();

  } //end of loop over events


  fout->cd();
  ztree->Write();
  // injetTree->Write();
  // trackTree->Write();
  fout->Write();
  fout->Close();

}

int main(int argc, char *argv[])
{
  if((argc < 3))
  {
    std::cout << "Usage: ./gammajetSkim.exe <inputfile> <outputfile> [jetname] [is_pp]" << std::endl;
    return 1;
  }
  if(argc==3)  gammajetSkim(argv[1], argv[2]);
  if(argc==4)  gammajetSkim(argv[1], argv[2], argv[3]);
  if(argc==5)  gammajetSkim(argv[1], argv[2], argv[3], std::atoi(argv[4]));
  return 0;
}
