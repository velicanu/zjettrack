   const static int NTRK = 100000 ; // 87441
   const static int NVTX = 10 ; // 53122
   // Declaration of leaf types

      //if((trkPt[i]-2*trkPtError[i])*TMath::CosH(trkEta[i])>15 && (trkPt[i]-2*trkPtError[i])*TMath::CosH(trkEta[i])>pfHcal[i]+pfEcal[i]) continue;} //Calo Matching
   
   Int_t           nEv_;
   Int_t           nLumi_;
   Int_t           nBX_;
   Int_t           nRun_;
   Int_t           N_;
   Int_t           nVtx_;
   Int_t           nTrk_;
   Int_t           maxPtVtx_;
   Int_t           maxMultVtx_;
   Int_t           * nTrkVtx_ = (Int_t*) malloc(sizeof(Int_t)*NVTX);   //[nVtx]
   Float_t         * normChi2Vtx_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtx]
   Float_t         * sumPtVtx_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtx]
   Float_t         * xVtx_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtx]
   Float_t         * yVtx_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtx]
   Float_t         * zVtx_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtx]
   Float_t         * xVtxErr_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtx]
   Float_t         * yVtxErr_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtx]
   Float_t         * zVtxErr_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtx]
   Float_t         * vtxDist2DErr_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtx]
   Float_t         * vtxDist2DSig_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtx]
   Float_t         * vtxDist2D_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtx]
   Float_t         * vtxDist3DErr_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtx]
   Float_t         * vtxDist3DSig_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtx]
   Int_t           nVtxSim_;
   Float_t         * xVtxSim_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtxSim]
   Float_t         * yVtxSim_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtxSim]
   Float_t         * zVtxSim_ = (Float_t*) malloc(sizeof(Float_t)*NVTX);   //[nVtxSim]
   Float_t         * trkPt_ = (Float_t*) malloc(sizeof(Float_t)*NTRK);   //[nTrk]
   Float_t         * trkPtError_ = (Float_t*) malloc(sizeof(Float_t)*NTRK);   //[nTrk]
   UChar_t         * trkNHit_ = (UChar_t*) malloc(sizeof(UChar_t)*NTRK);   //[nTrk]
   UChar_t         * trkNlayer_ = (UChar_t*) malloc(sizeof(UChar_t)*NTRK);   //[nTrk]
   Float_t         * trkEta_ = (Float_t*) malloc(sizeof(Float_t)*NTRK);   //[nTrk]
   Float_t         * trkPhi_ = (Float_t*) malloc(sizeof(Float_t)*NTRK);   //[nTrk]
   Int_t           * trkCharge_ = (Int_t*) malloc(sizeof(Int_t)*NTRK);   //[nTrk]
   UChar_t         * trkNVtx_ = (UChar_t*) malloc(sizeof(UChar_t)*NTRK);   //[nTrk]
   Int_t           nTrkTimesnVtx_;
   Bool_t          * trkAssocVtx_ = (Bool_t*) malloc(sizeof(Bool_t)*NTRK);   //[nTrkTimesnVtx]
   Bool_t          * highPurity_ = (Bool_t*) malloc(sizeof(Bool_t)*NTRK);   //[nTrk]
   Bool_t          * tight_ = (Bool_t*) malloc(sizeof(Bool_t)*NTRK);   //[nTrk]
   Bool_t          * loose_ = (Bool_t*) malloc(sizeof(Bool_t)*NTRK);   //[nTrk]
   Float_t         * trkChi2_ = (Float_t*) malloc(sizeof(Float_t)*NTRK);   //[nTrk]
   UChar_t         * trkNdof_ = (UChar_t*) malloc(sizeof(UChar_t)*NTRK);   //[nTrk]
   Float_t         * trkDxy1_ = (Float_t*) malloc(sizeof(Float_t)*NTRK);   //[nTrk]
   Float_t         * trkDxyError1_ = (Float_t*) malloc(sizeof(Float_t)*NTRK);   //[nTrk]
   Float_t         * trkDz1_ = (Float_t*) malloc(sizeof(Float_t)*NTRK);   //[nTrk]
   Float_t         * trkDzError1_ = (Float_t*) malloc(sizeof(Float_t)*NTRK);   //[nTrk]
   Bool_t          * trkFake_ = (Bool_t*) malloc(sizeof(Bool_t)*NTRK);   //[nTrk]
   UChar_t         * trkAlgo_ = (UChar_t*) malloc(sizeof(UChar_t)*NTRK);   //[nTrk]
   UChar_t         * trkOriginalAlgo_ = (UChar_t*) malloc(sizeof(UChar_t)*NTRK);   //[nTrk]
   Float_t         * trkMVA_ = (Float_t*) malloc(sizeof(Float_t)*NTRK);   //[nTrk]
   Int_t           * pfType_ = (Int_t*) malloc(sizeof(Int_t)*NTRK);   //[nTrk]
   Float_t         * pfCandPt_ = (Float_t*) malloc(sizeof(Float_t)*NTRK);   //[nTrk]
   Float_t         * pfEcal_ = (Float_t*) malloc(sizeof(Float_t)*NTRK);   //[nTrk]
   Float_t         * pfHcal_ = (Float_t*) malloc(sizeof(Float_t)*NTRK);   //[nTrk]


   
void initTrackTree(TTree *tree)
{
   tree->SetBranchAddress("nEv", &nEv_);
   tree->SetBranchAddress("nLumi", &nLumi_);
   tree->SetBranchAddress("nBX", &nBX_);
   tree->SetBranchAddress("nRun", &nRun_);
   tree->SetBranchAddress("N", &N_);
   tree->SetBranchAddress("nVtx", &nVtx_);
   tree->SetBranchAddress("nTrk", &nTrk_);
   tree->SetBranchAddress("maxPtVtx", &maxPtVtx_);
   tree->SetBranchAddress("maxMultVtx", &maxMultVtx_);
   tree->SetBranchAddress("nTrkVtx", nTrkVtx_);
   tree->SetBranchAddress("normChi2Vtx", normChi2Vtx_);
   tree->SetBranchAddress("sumPtVtx", sumPtVtx_);
   tree->SetBranchAddress("xVtx", xVtx_);
   tree->SetBranchAddress("yVtx", yVtx_);
   tree->SetBranchAddress("zVtx", zVtx_);
   tree->SetBranchAddress("xVtxErr", xVtxErr_);
   tree->SetBranchAddress("yVtxErr", yVtxErr_);
   tree->SetBranchAddress("zVtxErr", zVtxErr_);
   tree->SetBranchAddress("vtxDist2D", vtxDist2D_);
   tree->SetBranchAddress("vtxDist2DErr", vtxDist2DErr_);
   tree->SetBranchAddress("vtxDist2DSig", vtxDist2DSig_);
   tree->SetBranchAddress("vtxDist3DErr", vtxDist3DErr_);
   tree->SetBranchAddress("vtxDist3DSig", vtxDist3DSig_);
   tree->SetBranchAddress("nVtxSim", &nVtxSim_);
   tree->SetBranchAddress("xVtxSim", &xVtxSim_);
   tree->SetBranchAddress("yVtxSim", &yVtxSim_);
   tree->SetBranchAddress("zVtxSim", &zVtxSim_);
   tree->SetBranchAddress("trkPt", trkPt_);
   tree->SetBranchAddress("trkPtError", trkPtError_);
   tree->SetBranchAddress("trkNHit", trkNHit_);
   tree->SetBranchAddress("trkNlayer", trkNlayer_);
   tree->SetBranchAddress("trkEta", trkEta_);
   tree->SetBranchAddress("trkPhi", trkPhi_);
   tree->SetBranchAddress("trkCharge", trkCharge_);
   tree->SetBranchAddress("trkNVtx", trkNVtx_);
   tree->SetBranchAddress("nTrkTimesnVtx", &nTrkTimesnVtx_);
   tree->SetBranchAddress("trkAssocVtx", trkAssocVtx_);
   tree->SetBranchAddress("highPurity", highPurity_);
   tree->SetBranchAddress("tight", tight_);
   tree->SetBranchAddress("loose", loose_);
   tree->SetBranchAddress("trkChi2", trkChi2_);
   tree->SetBranchAddress("trkNdof", trkNdof_);
   tree->SetBranchAddress("trkDxy1", trkDxy1_);
   tree->SetBranchAddress("trkDxyError1", trkDxyError1_);
   tree->SetBranchAddress("trkDz1", trkDz1_);
   tree->SetBranchAddress("trkDzError1", trkDzError1_);
   tree->SetBranchAddress("trkFake", trkFake_);
   tree->SetBranchAddress("trkAlgo", trkAlgo_);
   tree->SetBranchAddress("trkOriginalAlgo", trkOriginalAlgo_);
   tree->SetBranchAddress("trkMVA", trkMVA_);
   tree->SetBranchAddress("pfType", pfType_);
   tree->SetBranchAddress("pfCandPt", pfCandPt_);
   tree->SetBranchAddress("pfEcal", pfEcal_);
   tree->SetBranchAddress("pfHcal", pfHcal_);
}

