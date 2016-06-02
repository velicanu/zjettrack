#ifndef GETTRKCORR
#define GETTRKCORR

#include "TMath.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TrkSettings.h"
#include <iostream>
#include <vector>

using namespace std;

class TrkCorr{
  public:
    double getTrkCorr(float pt, float eta, float phi, float hiBin, float rmin=99, int nev=0, float jtpt=0, int correction=0);
    TrkCorr(std::string inputDirectory = "trkCorrections/");
    ~TrkCorr();    

  // private:
    int nFiles;
    int nSteps;

    std::vector<std::vector<TH1D*> > eff;
    std::vector<std::vector<TH1D*> > fake;
    std::vector<std::vector<TH2D*> > eff2;
    std::vector<std::vector<TH2D*> > fake2;
    std::vector<TH2D*> secondary;
    std::vector<TH1D*> multiple;

    TrkSettings * s;
};

TrkCorr::TrkCorr(std::string inputDirectory)
{
  std::cout << "Initializing tracking correction files..." << std::endl;
  nFiles = 0;
  nSteps = 0;

  s = new TrkSettings(Form("%sTrkCorrInputFile.txt",inputDirectory.c_str()));
  for(int i = 0; i<s->nPtBinCoarse; i++)
  {
    nFiles += s->nCentPUBinCoarse[i];
  }
  nSteps = s->nStep;
  std::cout << "Correction tables reading " << nFiles << " Files, with " << nSteps << " iterations per file." << std::endl;

  TFile * f;
  for(int i = 0; i<nFiles; i++)
  {
    std::string isPP = "";
    if(s->nPb==0) isPP = "pp_";
    f = TFile::Open(Form("%s%scorrHists_job%d.root",inputDirectory.c_str(),isPP.c_str(),i),"read");
    std::vector<TH1D*> tempTH1EffVec;
    std::vector<TH2D*> tempTH2EffVec;
    std::vector<TH1D*> tempTH1FakeVec;
    std::vector<TH2D*> tempTH2FakeVec;
    for(int j = 0; j<nSteps; j++)
    {
      if(s->stepOrder.at(j)!=1 && s->stepOrder.at(j)!=7)
      {
	tempTH1EffVec.push_back((TH1D*)f->Get(Form("finalEff_type%d",j)));
        tempTH1EffVec.back()->SetDirectory(0);
        tempTH1FakeVec.push_back((TH1D*)f->Get(Form("finalFake_type%d",j)));
        tempTH1FakeVec.back()->SetDirectory(0);
      }
      else 
      {
	tempTH2EffVec.push_back((TH2D*)f->Get(Form("finalEff_type%d",j)));
        tempTH2EffVec.back()->SetDirectory(0);
        tempTH2FakeVec.push_back((TH2D*)f->Get(Form("finalFake_type%d",j)));
        tempTH2FakeVec.back()->SetDirectory(0);
      }
    }
    eff.push_back(tempTH1EffVec);
    fake.push_back(tempTH1FakeVec);
    eff2.push_back(tempTH2EffVec);
    fake2.push_back(tempTH2FakeVec);

    secondary.push_back((TH2D*)f->Get("SecondaryRate"));
    secondary.back()->SetDirectory(0);
    multiple.push_back((TH1D*)f->Get("MultipleRecoRate"));
    multiple.back()->SetDirectory(0);


    f->Close();
  }
  std::cout << "Initialization complete." << std::endl; 
}

//correction=0 is total, 1 is eff, 2 is fake, 3 is second, 4 is mult
double TrkCorr::getTrkCorr(float pt, float eta, float phi, float hiBin, float rmin, int nev, float jtpt, int correction)
{
  //  std::cout << "pre crash " <<pt << " " << eta << " " << phi  << " " <<rmin<<" "<<hiBin<< std::endl;
  if(pt<0.5 || pt>=400){  std::cout << "\nPt of " << pt << " less than 500 MeV or > 300 GeV, please place a cut to prevent this. Returning a correction of 1" << std::endl; return 1;}
  if(eta<-2.4 || eta>2.4){  std::cout << "\nEta outside of |2.4|, please place a cut to prevent this. Returning a correction of 1" << std::endl; return 1;}
  if(hiBin<0 || hiBin>199){  std::cout << "\nhiBin not within 0 to 200, please place a cut to prevent this.  Returning a correction of 1" << std::endl; return 1;}
  
  //calculating what file to take corrections out of 
  int coarseBin = 0;
  float cent = hiBin;
  if(s->nPb==2) cent = cent/2.0;
  // std::cout<<"pre first for loop"<<endl;
  for(int i = 0; i<s->nPtBinCoarse; i++)
  {
    if(pt >= s->ptBinCoarse[i+1]) coarseBin+=s->nCentPUBinCoarse[i];
    else
    {
      for(int j = 0; j<s->nCentPUBinCoarse[i]; j++)
      {
        if(cent >= s->centPUBinCoarse[i][j+1]) coarseBin++;
      }
      break;
    }
  }
  //end bin calculation
 
  float netEff = 1;
  float netFake = 1;
  float netSec = 0;
  float netMult = 0; 

  int th1indx = 0;
  int th2indx = 0;
  //std::cout<<"pre second for loop"<< coarseBin<<" "<<th1indx<<endl;
  for(int j = 0; j<nSteps; j++)
  {
    //    std::cout<<"in for "<<j<<std::endl;
    if(th1indx==0 && coarseBin==10){
        
      // std::cout<<"would break "<<0<<" "<<coarseBin<<" "<<th1indx<<" "<<pt<<" "<<nev<<std::endl;
      continue;
    }
    if(s->stepOrder.at(j)==0)
    {
      // std::cout<<"in if "<<0<<" "<<coarseBin<<" "<<th1indx<<" "<<pt<<std::endl;
      // std::cout<<"in if "<<fake[coarseBin][th1indx]->GetBinContent(fake[coarseBin][th1indx]->FindBin(pt))<<std::endl;
      if(eff[coarseBin][th1indx]->GetBinContent(eff[coarseBin][th1indx]->FindBin(pt))>0)
        netEff *= eff[coarseBin][th1indx]->GetBinContent(eff[coarseBin][th1indx]->FindBin(pt));
      if(fake[coarseBin][th1indx]->GetBinContent(fake[coarseBin][th1indx]->FindBin(pt))>0)
        netFake *= fake[coarseBin][th1indx]->GetBinContent(fake[coarseBin][th1indx]->FindBin(pt));
    }
    if(s->stepOrder.at(j)==1)
    {
      // std::cout<<"in if "<<1<<std::endl;
      if(eff2[coarseBin][th2indx]->GetBinContent(eff2[coarseBin][th2indx]->GetXaxis()->FindBin(eta),eff2[coarseBin][th2indx]->GetYaxis()->FindBin(phi))>0)
        netEff *= eff2[coarseBin][th2indx]->GetBinContent(eff2[coarseBin][th2indx]->GetXaxis()->FindBin(eta),eff2[coarseBin][th2indx]->GetYaxis()->FindBin(phi));
      if(fake2[coarseBin][th2indx]->GetBinContent(fake2[coarseBin][th2indx]->GetXaxis()->FindBin(eta),fake2[coarseBin][th2indx]->GetYaxis()->FindBin(phi)))
        netFake *= fake2[coarseBin][th2indx]->GetBinContent(fake2[coarseBin][th2indx]->GetXaxis()->FindBin(eta),fake2[coarseBin][th2indx]->GetYaxis()->FindBin(phi));
    }
    if(s->stepOrder.at(j)==2)
    {
      // std::cout<<"in if "<<2<<std::endl;
      if(eff[coarseBin][th1indx]->GetBinContent(eff[coarseBin][th1indx]->FindBin(pt))>0)
      netEff *= eff[coarseBin][th1indx]->GetBinContent(eff[coarseBin][th1indx]->FindBin(cent));
      if(fake[coarseBin][th1indx]->GetBinContent(fake[coarseBin][th1indx]->FindBin(pt))>0)
      netFake *= fake[coarseBin][th1indx]->GetBinContent(fake[coarseBin][th1indx]->FindBin(cent));
    }
    if(s->stepOrder.at(j)==3)
    {
      // std::cout<<"in if "<<3<<std::endl;
      if(eff[coarseBin][th1indx]->GetBinContent(eff[coarseBin][th1indx]->FindBin(pt))>0)
      netEff *= eff[coarseBin][th1indx]->GetBinContent(eff[coarseBin][th1indx]->FindBin(jtpt));
      if(fake[coarseBin][th1indx]->GetBinContent(fake[coarseBin][th1indx]->FindBin(pt))>0)
      netFake *= fake[coarseBin][th1indx]->GetBinContent(fake[coarseBin][th1indx]->FindBin(jtpt));
    }
    if(s->stepOrder.at(j)==4)
    {
      // std::cout<<"in if "<<4<<std::endl;
      if(eff[coarseBin][th1indx]->GetBinContent(eff[coarseBin][th1indx]->FindBin(pt))>0)
      netEff *= eff[coarseBin][th1indx]->GetBinContent(eff[coarseBin][th1indx]->FindBin(eta));
      if(fake[coarseBin][th1indx]->GetBinContent(fake[coarseBin][th1indx]->FindBin(pt))>0)
      netFake *= fake[coarseBin][th1indx]->GetBinContent(fake[coarseBin][th1indx]->FindBin(eta));
    }
    if(s->stepOrder.at(j)==5)
    {
      // std::cout<<"in if "<<5<<std::endl;
      if(eff[coarseBin][th1indx]->GetBinContent(eff[coarseBin][th1indx]->FindBin(pt))>0)
      netEff *= eff[coarseBin][th1indx]->GetBinContent(eff[coarseBin][th1indx]->FindBin(rmin));
      if(fake[coarseBin][th1indx]->GetBinContent(fake[coarseBin][th1indx]->FindBin(pt))>0)
      netFake *= fake[coarseBin][th1indx]->GetBinContent(fake[coarseBin][th1indx]->FindBin(rmin));
    }
    if(s->stepOrder.at(j)==7)
    {
      // std::cout<<"in if "<<7<<std::endl;
      if(eff[coarseBin][th1indx]->GetBinContent(eff[coarseBin][th1indx]->FindBin(pt))>0)
      netEff *= eff2[coarseBin][th2indx]->GetBinContent(eff2[coarseBin][th2indx]->GetXaxis()->FindBin(eta),eff2[coarseBin][th2indx]->GetYaxis()->FindBin(pt));
      if(fake[coarseBin][th1indx]->GetBinContent(fake[coarseBin][th1indx]->FindBin(pt))>0)
      netFake *= fake2[coarseBin][th2indx]->GetBinContent(fake2[coarseBin][th2indx]->GetXaxis()->FindBin(eta),fake2[coarseBin][th2indx]->GetYaxis()->FindBin(pt));
    }

    if(s->stepOrder.at(j)==1 || s->stepOrder.at(j)==7) th2indx++;
    else                      th1indx++;
  }
  //std::cout<<"post second for loop"<<endl;
  if(th1indx==0 && coarseBin==10)
  {
    // std::cout<<"would break post"<<" "<<coarseBin<<" "<<th1indx<<" "<<pt<<" "<<nev<<std::endl;
    // continue;
  }
  else
  {
    netMult = multiple[coarseBin]->GetBinContent(multiple[coarseBin]->FindBin(pt));
    netSec  = secondary[coarseBin]->GetBinContent(secondary[coarseBin]->GetXaxis()->FindBin(pt),secondary[coarseBin]->GetYaxis()->FindBin(eta));
  }
/*
  netEff *= eff[coarseBin][0]->GetBinContent(eff[coarseBin][0]->FindBin(pt));
  netEff *= eff[coarseBin][1]->GetBinContent(eff[coarseBin][1]->FindBin(cent));
  netEff *= eff2[coarseBin][2]->GetBinContent(eff2[coarseBin][2]->GetXaxis()->FindBin(eta),eff2[coarseBin][2]->GetYaxis()->FindBin(phi));
  netEff *= eff[coarseBin][3]->GetBinContent(eff[coarseBin][3]->FindBin(rmin));
  
  netFake *= fake[coarseBin][0]->GetBinContent(fake[coarseBin][0]->FindBin(pt));
  netFake *= fake[coarseBin][1]->GetBinContent(fake[coarseBin][1]->FindBin(cent));
  netFake *= fake2[coarseBin][2]->GetBinContent(fake2[coarseBin][2]->GetXaxis()->FindBin(eta),fake2[coarseBin][2]->GetYaxis()->FindBin(phi));
  netFake *= fake[coarseBin][3]->GetBinContent(fake[coarseBin][3]->FindBin(rmin));*/

  if(netFake<1) netFake = 1;
  if(netEff>1)  netEff = 1;

//  std::cout << coarseBin << std::endl;
//  std::cout << "pt: " << pt << " cent: " << cent << " eta: " << eta << " phi: " << phi   << std::endl;
//  std::cout << "Efficiency: " << netEff << " (pt: " << eff[coarseBin][0]->GetBinContent(eff[coarseBin][0]->FindBin(pt)) << " cent: " <<  eff[coarseBin][1]->GetBinContent(eff[coarseBin][1]->FindBin(cent)) << " eta/phi: " << eff2[coarseBin][2]->GetBinContent(eff2[coarseBin][2]->GetXaxis()->FindBin(eta),eff2[coarseBin][2]->GetYaxis()->FindBin(phi)) << " rmin: " << eff[coarseBin][3]->GetBinContent(eff[coarseBin][3]->FindBin(rmin))<< std::endl;
//   std::cout << "Fake Rate: " << (1-1./netFake) << " (pt: " << 1-1./fake[coarseBin][0]->GetBinContent(fake[coarseBin][0]->FindBin(pt)) << " cent: " << 1-1./fake[coarseBin][1]->GetBinContent(fake[coarseBin][1]->FindBin(cent)) << " eta/phi: " << 1-1./fake2[coarseBin][2]->GetBinContent(fake2[coarseBin][2]->GetXaxis()->FindBin(eta),fake2[coarseBin][2]->GetYaxis()->FindBin(phi)) << "rmin: " << fake[coarseBin][3]->GetBinContent(fake[coarseBin][3]->FindBin(rmin)) << std::endl;
//  std::cout << "Secondary Rate: " <<  netSec << std::endl;
//  std::cout << "Multiple Reco Rate: " << netMult << "\nTotal Correction: " << (1.0-netSec)/(netEff*netFake*(1+netMult)) << std::endl;


  if(1/netEff>100000) std::cout << "problem here! " << netEff <<  " " <<pt << " " << eta << " " << phi << " " << " " << coarseBin << " " <<rmin<<" "<<hiBin<< std::endl;

  if(correction==1) return 1/(netEff);
  else if(correction==2) return 1/(netFake);
  else if(correction==3) return 1-netSec;
  else if(correction==4) return 1/(1+netMult);
  else 
  {
    if(s->nPb==0) return (1-netSec)/(netEff*netFake);
    if(s->nPb==2) return 1.0/(netEff*netFake);
  }

  //return -1 (if there is a problem and the code gets to here)
  return -1;
//  else return 1.0/(netEff*netFake*(1+netMult)); */
}

TrkCorr::~TrkCorr()
{
  delete s;
}
#endif
