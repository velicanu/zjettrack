#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

using namespace std;
class L2L3Residual
{
 private:
  static const int neta=16;
  double lower_pt_cut;
  double higher_pt_cut;
  int radius;
  TFile *correction_file;
  TF1 * fits[neta];
  TString algo_corr;
  public:
  
   double eta_min[neta];
   double eta_max[neta];
   void reset()
   { 
    for(int ieta=0;ieta<neta;ieta++){
     fits[ieta]=NULL;
    }

    correction_file=NULL;
   
    eta_min[0] = -3;
    eta_max[0] = eta_min[1] = -2.500;
    eta_max[1] = eta_min[2] = -2.172;
    eta_max[2] = eta_min[3] = -1.740;
    eta_max[3] = eta_min[4] = -1.392;
    eta_max[4] = eta_min[5] = -1.044;
    eta_max[5] = eta_min[6] = -0.696;
    eta_max[6] = eta_min[7] = -0.348;
    eta_max[7] = eta_min[8] = 0.000;
    eta_max[8] = eta_min[9] = 0.348;
    eta_max[9] = eta_min[10] = 0.696;
    eta_max[10] = eta_min[11] = 1.044;
    eta_max[11] = eta_min[12] = 1.392;
    eta_max[12] = eta_min[13] = 1.740;
    eta_max[13] = eta_min[14] = 2.172;
    eta_max[14] = eta_min[15] = 2.500;
    eta_max[15] = 3;
   }
  
  L2L3Residual(int radius=3)
  {
   reset();
   this->radius=radius;
   algo_corr=Form("ak%dPF",radius);   
   correction_file = new TFile(Form("L2L3VsPtEtaBinned_%s.root",algo_corr.Data()));
   for(int ieta=0;ieta<neta;ieta++){
    fits[ieta] = (TF1*)correction_file->Get(Form("fit%d",ieta));
   } 
   lower_pt_cut = 20;
   higher_pt_cut = 400;
  }
   
  
  double get_corrected_pt(double jetpt, double jeteta)
  {
   double correction = 1;
   if( abs(jeteta)> 3) return correction*jetpt;
   if(jetpt < lower_pt_cut || jetpt > higher_pt_cut) return correction*jetpt;
   
   int etaindex = 0;
   for(int ieta = 0; ieta < neta; ieta++){
    if(eta_min[ieta] > jeteta ) continue;
	else etaindex = ieta;
   }
   return fits[etaindex]->Eval(jetpt)*jetpt;
  }
};