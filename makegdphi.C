{
  string thistype = "mm";
  TFile * ppfile = TFile::Open("/export/d00/scratch/dav2105/ztrees/g.pp-photonHLTFilter-v0-HiForest.root");
  TFile * pbpbfile = TFile::Open("/export/d00/scratch/dav2105/ztrees/g.newtrkcuts.HIPhoton40AndZ-PbPb-photonHLTFilter-v3.root");

  TTree * ppztree = (TTree*) ppfile->Get("ztree");
  TTree * pbpbztree = (TTree*) pbpbfile->Get("ztree");

  const int ntrkcuts = 6;

  int trkptcuts [ntrkcuts] = {4,8};
  float trkptmin [ntrkcuts] = {0.5,1,2,3,4,8};
  float trkptmax [ntrkcuts] = {1,2,3,4,8,1000};
  int zptcuts   [2] = {200, 220};
  float ymins [2] = {0,0};
  float subymins [2] = {-1,-1};
  float ymaxs [2] = {6,3};
  float subymaxs [2] = {5,2};
  int histcolors [ntrkcuts] = {TColor::GetColor("#9999FF"),TColor::GetColor("#FFFF99"),TColor::GetColor("#FF9933"),TColor::GetColor("#9933CC"),TColor::GetColor("#006600"),TColor::GetColor("#CC0000")};
  TH1D * hppzdphi [2][ntrkcuts];
  TH1D * hsubppzdphi [2][ntrkcuts];
  TH1D * hpbpbzdphi [2][ntrkcuts];
  TH1D * hsubpbpbzdphi [2][ntrkcuts];
  TF1 * flat[2][ntrkcuts][2]; // [izcut][itrkcut][pp or PbPb]
  string collision[2] = {"pp","PbPb"};
  int linestyle[2] = {2,1};
  TCanvas * canvases1D [2][ntrkcuts];
  TCanvas * canvases [2][ntrkcuts];
  TCanvas * subcanvases [2][ntrkcuts];
  TCanvas * stackedcanvasepp [2];
  TCanvas * stackedcanvasepbpb [2];
  THStack * hspp[2];
  THStack * hspbpb[2];
  int itrkstart = 1;
  for(int izcut = 0 ; izcut < 2 ; ++izcut)
  {
    for(int itrkcut = itrkstart ; itrkcut < ntrkcuts ; ++itrkcut)
    {

      hppzdphi[izcut][itrkcut] = new TH1D(Form("hppzdphi_%d_%d_%d",zptcuts[izcut],int(trkptmin[itrkcut]),int(trkptmax[itrkcut])),";#gamma-trk #Delta#phi",12,0,3.1415);
      hpbpbzdphi[izcut][itrkcut] = new TH1D(Form("hpbpbzdphi_%d_%d_%d",zptcuts[izcut],int(trkptmin[itrkcut]),int(trkptmax[itrkcut])),";#gamma-trk #Delta#phi",12,0,3.1415);
      int nppzpassingcuts = ppztree->Draw("phoEt[0]",Form("phoEt[0]>%d",zptcuts[izcut]),"goff");
      int npbpbzpassingcuts = pbpbztree->Draw("phoEt[0]",Form("phoEt[0]>%d",zptcuts[izcut]),"goff");
      float ppbinwidth   = hppzdphi[izcut][itrkcut]->GetBinWidth(1);
      float pbpbbinwidth = hpbpbzdphi[izcut][itrkcut]->GetBinWidth(1);
      // ppztree->Draw(Form("acos(cos(phoPhi[0] - trkPhi))>>hppzdphi_%d_%d_%d",zptcuts[izcut],int(trkptmin[itrkcut]),int(trkptmax[itrkcut])),Form("(phoEt[0]>%d && trkPt>%2.1f && trkPt<%2.1f)*trkPt*highPurity",zptcuts[izcut],trkptmin[itrkcut],trkptmax[itrkcut]),"goff");
      // pbpbztree->Draw(Form("acos(cos(phoPhi[0] - trkPhi))>>hpbpbzdphi_%d_%d_%d",zptcuts[izcut],int(trkptmin[itrkcut]),int(trkptmax[itrkcut])),Form("(phoEt[0]>%d && trkPt>%2.1f && trkPt<%2.1f)*trkPt*highPurity",zptcuts[izcut],trkptmin[itrkcut],trkptmax[itrkcut]),"goff");
      ppztree->Draw(Form("acos(cos(phoPhi[0] - trkPhi))>>hppzdphi_%d_%d_%d",zptcuts[izcut],int(trkptmin[itrkcut]),int(trkptmax[itrkcut])),Form("(phoEt[0]>%d && trkPt>%2.1f && trkPt<%2.1f)*trkPt*highPurity*trkWeight",zptcuts[izcut],trkptmin[itrkcut],trkptmax[itrkcut]),"goff");
      pbpbztree->Draw(Form("acos(cos(phoPhi[0] - trkPhi))>>hpbpbzdphi_%d_%d_%d",zptcuts[izcut],int(trkptmin[itrkcut]),int(trkptmax[itrkcut])),Form("(phoEt[0]>%d && trkPt>%2.1f && trkPt<%2.1f)*trkPt*highPurity*trkWeight",zptcuts[izcut],trkptmin[itrkcut],trkptmax[itrkcut]),"goff");
      hppzdphi[izcut][itrkcut]->Sumw2();
      hppzdphi[izcut][itrkcut]->SetMarkerStyle(24);
      hpbpbzdphi[izcut][itrkcut]->Sumw2();

      hppzdphi[izcut][itrkcut]->Scale(1.0/(float)nppzpassingcuts/ppbinwidth);
      hpbpbzdphi[izcut][itrkcut]->Scale(1.0/(float)npbpbzpassingcuts/pbpbbinwidth);

      hsubppzdphi[izcut][itrkcut] = (TH1D*) hppzdphi[izcut][itrkcut]->Clone(Form("hsubppzdphi_%d_%d_%d",zptcuts[izcut],int(trkptmin[itrkcut]),int(trkptmax[itrkcut])));
      hsubpbpbzdphi[izcut][itrkcut] = (TH1D*) hpbpbzdphi[izcut][itrkcut]->Clone(Form("hsubpbpbzdphi_%d_%d_%d",zptcuts[izcut],int(trkptmin[itrkcut]),int(trkptmax[itrkcut])));


      canvases[izcut][itrkcut] = new TCanvas();
      float maxvalue = hppzdphi[izcut][itrkcut]->GetMaximum();
      if( hpbpbzdphi[izcut][itrkcut]->GetMaximum() > maxvalue ) maxvalue = hpbpbzdphi[izcut][itrkcut]->GetMaximum();
      float minvalue = hppzdphi[izcut][itrkcut]->GetMinimum();
      if( hpbpbzdphi[izcut][itrkcut]->GetMinimum() < minvalue ) minvalue = hpbpbzdphi[izcut][itrkcut]->GetMinimum();

      TH2D * dummy = new TH2D("dummy",";#gamma-trk #Delta#phi;#frac{1}{N_{#gamma}} #frac{dN}{d#Delta#phi}",1,0,3.1415,1,minvalue-1,maxvalue+1);
      dummy->GetXaxis()->CenterTitle();
      dummy->GetYaxis()->CenterTitle();
      dummy->GetYaxis()->SetTitleOffset(1.25);
      TH2D * subdummy = new TH2D("subdummy",";#gamma-trk #Delta#phi;subtracted #frac{1}{N_{#gamma}} #frac{dN}{d#Delta#phi}",1,0,3.1415,1,-4,13);
      subdummy->GetXaxis()->CenterTitle();
      subdummy->GetYaxis()->CenterTitle();
      subdummy->GetYaxis()->SetTitleOffset(1.25);

      //canvases1D[izcut][itrkcut] = new TCanvas(Form("canvas1D_%d_%d",izcut,itrkcut));
      dummy->Draw();
      hppzdphi[izcut][itrkcut]->Draw("same");
      hpbpbzdphi[izcut][itrkcut]->Draw("same");
      // break;

      for(int icoll = 0 ; icoll < 2 ; ++icoll)
      {
        flat[izcut][itrkcut][icoll] = new TF1(Form("flat_%s",collision[icoll].data()),"[0]", 0, TMath::Pi());
        flat[izcut][itrkcut][icoll]->SetParameter(0,0);
        if(icoll==0)
        {
          hppzdphi[izcut][itrkcut]->Fit(Form("flat_%s",collision[icoll].data()),"0","",0, TMath::Pi()/2.0);
        }
        else
        {
          hpbpbzdphi[izcut][itrkcut]->Fit(Form("flat_%s",collision[icoll].data()),"0","",0, TMath::Pi()/2.0);
        }
        flat[izcut][itrkcut][icoll]->SetLineStyle(linestyle[icoll]);
        // flat[izcut][itrkcut][icoll]->Draw("same");
      }

      TLegend * leg = new TLegend(0.45,0.58,0.75,0.9);
      leg->AddEntry(hpbpbzdphi[izcut][itrkcut],"#sqrt{s_{NN}}=5TeV PbPb","lp");
      leg->AddEntry(hppzdphi[izcut][itrkcut],"#sqrt{s}   =5TeV pp","lp");
      leg->AddEntry(hppzdphi[izcut][itrkcut],Form("#gamma p_{T}>%d GeV",zptcuts[izcut]),"");
      leg->AddEntry(hppzdphi[izcut][itrkcut],Form("%2.1f<trk p_{T}<%2.1f GeV",trkptmin[itrkcut],trkptmax[itrkcut]),"");
      // leg->AddEntry(hppzdphi[izcut][itrkcut],Form("80<Zmass<110 GeV"),"");
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->Draw();


      TLatex * lint = new TLatex(.22,0.95,"CMS Preliminary");
      lint->SetNDC();
      lint->Draw();
      TLatex * lztype = new TLatex(.52,0.95,Form("%d pp, %d PbPb #gamma",nppzpassingcuts,npbpbzpassingcuts));
      // TLatex * lztype = new TLatex(.52,0.95,Form("%d pp, %d PbPb Z-#mu#mu",nppzpassingcuts,npbpbzpassingcuts));
      lztype->SetNDC();
      lztype->Draw();

      canvases[izcut][itrkcut]->SaveAs(Form("gtrkdphi_%d_%d_%d_%s.png",zptcuts[izcut],int(trkptmin[itrkcut]),int(trkptmax[itrkcut]),thistype.data()));

      subcanvases[izcut][itrkcut] = new TCanvas();
      for(int ibin = 0 ; ibin < hsubppzdphi[izcut][itrkcut]->GetNbinsX()+1 ; ++ibin)
      {
        hsubppzdphi[izcut][itrkcut]->SetBinContent(ibin,hppzdphi[izcut][itrkcut]->GetBinContent(ibin) - flat[izcut][itrkcut][0]->Eval(1.0));
        hsubpbpbzdphi[izcut][itrkcut]->SetBinContent(ibin,hpbpbzdphi[izcut][itrkcut]->GetBinContent(ibin) - flat[izcut][itrkcut][1]->Eval(1.0));
      }

      subdummy->Draw();
      hsubppzdphi[izcut][itrkcut]->Draw("same");
      hsubpbpbzdphi[izcut][itrkcut]->Draw("same");
      leg->Draw();
      lint->Draw();
      lztype->Draw();
      TLine * lzero = new TLine(0,0,TMath::Pi(),0);
      lzero->SetLineStyle(9);
      lzero->Draw();
      subcanvases[izcut][itrkcut]->SaveAs(Form("sub_gtrkdphi_%d_%d_%d_%s.png",zptcuts[izcut],int(trkptmin[itrkcut]),int(trkptmax[itrkcut]),thistype.data()));
    }
    // break;
    stackedcanvasepp[izcut] = new TCanvas();
    hspp[izcut] = new THStack(Form("hspp_%d",zptcuts[izcut]),"");

    TLegend * leg = new TLegend(0.25,0.45,0.55,0.9);
    leg->SetFillColor(0);
    leg->SetTextSize(0.05);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);


    TLatex * lint = new TLatex(.22,0.95,"CMS Preliminary");
    lint->SetNDC();
    TLatex * lztype = new TLatex(.52,0.95,Form("pp #gamma"));
    // TLatex * lztype = new TLatex(.52,0.95,Form("%d pp, %d PbPb Z-#mu#mu",nppzpassingcuts,npbpbzpassingcuts));
    lztype->SetNDC();

    for(int itrkcut = itrkstart ; itrkcut < ntrkcuts ; ++itrkcut)
    {
      hsubppzdphi[izcut][itrkcut]->SetFillColor(histcolors[itrkcut]);
      hspp[izcut]->Add(hsubppzdphi[izcut][itrkcut]);
      leg->AddEntry(hsubppzdphi[izcut][itrkcut],Form("%2.1f<trk p_{T}<%2.1f GeV",trkptmin[itrkcut],trkptmax[itrkcut]),"fl");
    }
    leg->AddEntry(hppzdphi[izcut][itrkstart],Form("#gamma p_{T}>%d GeV",zptcuts[izcut]),"");
    // leg->AddEntry(hppzdphi[izcut][itrkstart],Form("80<Zmass<110 GeV"),"");

    // hspp[izcut]->GetXaxis()->SetTitle("Z-trk #Delta#phi");
    hspp[izcut]->Draw("hist");
    hspp[izcut]->Draw("pe same");
    lint->Draw();
    leg->Draw();
    lztype->Draw();
    stackedcanvasepp[izcut]->SaveAs(Form("sub_gtrkdphi_stacked_pp_%d.png",zptcuts[izcut]));

    stackedcanvasepbpb[izcut] = new TCanvas();

    hspbpb[izcut] = new THStack(Form("hspbpb_%d",zptcuts[izcut]),"");
    TLegend * legpbpb = new TLegend(0.25,0.45,0.55,0.9);
    legpbpb->SetFillColor(0);
    legpbpb->SetTextSize(0.05);
    legpbpb->SetFillStyle(0);
    legpbpb->SetTextFont(42);
    TLatex * lztypepbpb = new TLatex(.52,0.95,Form("PbPb #gamma"));
    // TLatex * lztype = new TLatex(.52,0.95,Form("%d pbpb, %d PbPb Z-#mu#mu",npbpbzpassingcuts,npbpbzpassingcuts));
    lztypepbpb->SetNDC();

    for(int itrkcut = itrkstart ; itrkcut < ntrkcuts ; ++itrkcut)
    {
      hsubpbpbzdphi[izcut][itrkcut]->SetFillColor(histcolors[itrkcut]);
      hspbpb[izcut]->Add(hsubpbpbzdphi[izcut][itrkcut]);
      legpbpb->AddEntry(hsubpbpbzdphi[izcut][itrkcut],Form("%2.1f<trk p_{T}<%2.1f GeV",trkptmin[itrkcut],trkptmax[itrkcut]),"fl");
    }
    legpbpb->AddEntry(hpbpbzdphi[izcut][itrkstart],Form("#gamma p_{T}>%d GeV",zptcuts[izcut]),"");
    //legpbpb->AddEntry(hpbpbzdphi[izcut][itrkstart],Form("80<Zmass<110 GeV"),"");

    hspbpb[izcut]->Draw("hist");
    hspbpb[izcut]->Draw("pe same");
    lint->Draw();
    legpbpb->Draw();
    lztypepbpb->Draw();
    stackedcanvasepbpb[izcut]->SaveAs(Form("sub_gtrkdphi_stacked_pbpb_%d.png",zptcuts[izcut]));
    // break;
  }
}
