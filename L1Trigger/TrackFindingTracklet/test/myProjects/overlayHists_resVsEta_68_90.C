//--------------------------------------------------------------------------------------------------------
// Opens root file in another directory ("dir") which contains histograms and overlays them.
// In this case, it overlays some property ("prop") of three different root files like "eff_eta" or something.
// Saves plots to a chosen directory "saveDir"
//
// (SetPlotStyle function comes from Louise Skinnari's codes. It's an ATLAS plot style macro)
//
// By David Abrahamyan July 2023
//--------------------------------------------------------------------------------------------------------

//#include "plotstyle.h"

void SetPlotStyle(); // Sets plot style to some CMS thing I think idk Emily recommended I add this

double GetMaxHists(vector<TH1*> hists); // Gets maximum value of three histograms

void overlayHists_resVsEta_68_90 (){
    
    SetPlotStyle();

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";
    
    // Property you want to plot
    TString prop = "resVsEta_"; 

    // muon
    TFile *muonFile= new TFile(dir + "output_SingleMuon_PU0_D88.root");
    TH1F *muonEtaHist68= (TH1F*)muonFile->Get(prop + "eta_68");
    TH1F *muonPtHist68= (TH1F*)muonFile->Get(prop + "ptRel_68");
    TH1F *muonPhiHist68= (TH1F*)muonFile->Get(prop + "phi_68");
    TH1F *muonZ0Hist68= (TH1F*)muonFile->Get(prop + "z0_68");

    TH1F *muonEtaHist90= (TH1F*)muonFile->Get(prop + "eta_90");
    TH1F *muonPtHist90= (TH1F*)muonFile->Get(prop + "ptRel_90");
    TH1F *muonPhiHist90= (TH1F*)muonFile->Get(prop + "phi_90");
    TH1F *muonZ0Hist90= (TH1F*)muonFile->Get(prop + "z0_90");
    
    // electron
    TFile *electronFile= new TFile(dir + "output_SingleElectronPU0D88.root");
    TH1F *electronEtaHist68= (TH1F*)electronFile->Get(prop + "eta_68");
    TH1F *electronPtHist68= (TH1F*)electronFile->Get(prop + "ptRel_68");
    TH1F *electronPhiHist68= (TH1F*)electronFile->Get(prop + "phi_68");
    TH1F *electronZ0Hist68= (TH1F*)electronFile->Get(prop + "z0_68");

    TH1F *electronEtaHist90= (TH1F*)electronFile->Get(prop + "eta_90");
    TH1F *electronPtHist90= (TH1F*)electronFile->Get(prop + "ptRel_90");
    TH1F *electronPhiHist90= (TH1F*)electronFile->Get(prop + "phi_90");
    TH1F *electronZ0Hist90= (TH1F*)electronFile->Get(prop + "z0_90");

    // TTbarPU200
    TFile *TTbarPU200File= new TFile(dir + "output_TTbar_PU200_D88.root");
    TH1F *TTbarPU200EtaHist68= (TH1F*)TTbarPU200File->Get(prop + "eta_68");
    TH1F *TTbarPU200PtHist68= (TH1F*)TTbarPU200File->Get(prop + "ptRel_68");
    TH1F *TTbarPU200PhiHist68= (TH1F*)TTbarPU200File->Get(prop + "phi_68");
    TH1F *TTbarPU200Z0Hist68= (TH1F*)TTbarPU200File->Get(prop + "z0_68");

    TH1F *TTbarPU200EtaHist90= (TH1F*)TTbarPU200File->Get(prop + "eta_90");
    TH1F *TTbarPU200PtHist90= (TH1F*)TTbarPU200File->Get(prop + "ptRel_90");
    TH1F *TTbarPU200PhiHist90= (TH1F*)TTbarPU200File->Get(prop + "phi_90");
    TH1F *TTbarPU200Z0Hist90= (TH1F*)TTbarPU200File->Get(prop + "z0_90");

    // TTbarPU0
    TFile *TTbarPU0File= new TFile(dir + "output_TTbar_PU0_D88.root");
    TH1F *TTbarPU0EtaHist68= (TH1F*)TTbarPU0File->Get(prop + "eta_68");
    TH1F *TTbarPU0PtHist68= (TH1F*)TTbarPU0File->Get(prop + "ptRel_68");
    TH1F *TTbarPU0PhiHist68= (TH1F*)TTbarPU0File->Get(prop + "phi_68");
    TH1F *TTbarPU0Z0Hist68= (TH1F*)TTbarPU0File->Get(prop + "z0_68");

    TH1F *TTbarPU0EtaHist90= (TH1F*)TTbarPU0File->Get(prop + "eta_90");
    TH1F *TTbarPU0PtHist90= (TH1F*)TTbarPU0File->Get(prop + "ptRel_90");
    TH1F *TTbarPU0PhiHist90= (TH1F*)TTbarPU0File->Get(prop + "phi_90");
    TH1F *TTbarPU0Z0Hist90= (TH1F*)TTbarPU0File->Get(prop + "z0_90");

    

    // Make vectors for ease of use (maybe? lets see)
    vector<TH1*> allHists = {muonEtaHist68, muonPtHist68, muonPhiHist68, muonZ0Hist68,
    electronEtaHist68, electronPtHist68, electronPhiHist68, electronZ0Hist68,
    TTbarPU200EtaHist68, TTbarPU200PtHist68, TTbarPU200PhiHist68, TTbarPU200Z0Hist68,
    TTbarPU0EtaHist68, TTbarPU0PtHist68, TTbarPU0PhiHist68, TTbarPU0Z0Hist68,
    muonEtaHist90, muonPtHist90, muonPhiHist90, muonZ0Hist90,
    electronEtaHist90, electronPtHist90, electronPhiHist90, electronZ0Hist90,
    TTbarPU200EtaHist90, TTbarPU200PtHist90, TTbarPU200PhiHist90, TTbarPU200Z0Hist90,
    TTbarPU0EtaHist90, TTbarPU0PtHist90, TTbarPU0PhiHist90, TTbarPU0Z0Hist90};

    vector<TH1*> muonHists {muonEtaHist68, muonPtHist68, muonPhiHist68, muonZ0Hist68,
    muonEtaHist90, muonPtHist90, muonPhiHist90, muonZ0Hist90};

    vector<TH1*> muonHists68 {muonEtaHist68, muonPtHist68, muonPhiHist68, muonZ0Hist68};
    vector<TH1*> muonHists90 {muonEtaHist90, muonPtHist90, muonPhiHist90, muonZ0Hist90};

    vector<TH1*> electronHists {electronEtaHist68, electronPtHist68, electronPhiHist68, electronZ0Hist68,
    electronEtaHist90, electronPtHist90, electronPhiHist90, electronZ0Hist90};

    vector<TH1*> electronHists68 {electronEtaHist68, electronPtHist68, electronPhiHist68, electronZ0Hist68};
    vector<TH1*> electronHists90 {electronEtaHist90, electronPtHist90, electronPhiHist90, electronZ0Hist90};

    vector<TH1*> TTbarPU200Hists {TTbarPU200EtaHist68, TTbarPU200PtHist68, TTbarPU200PhiHist68, TTbarPU200Z0Hist68,
    TTbarPU200EtaHist90, TTbarPU200PtHist90, TTbarPU200PhiHist90, TTbarPU200Z0Hist90};

    vector<TH1*> TTbarPU200Hists68 {TTbarPU200EtaHist68, TTbarPU200PtHist68, TTbarPU200PhiHist68, TTbarPU200Z0Hist68};
    vector<TH1*> TTbarPU200Hists90 {TTbarPU200EtaHist90, TTbarPU200PtHist90, TTbarPU200PhiHist90, TTbarPU200Z0Hist90};

    vector<TH1*> TTbarPU0Hists {TTbarPU0EtaHist68, TTbarPU0PtHist68, TTbarPU0PhiHist68, TTbarPU0Z0Hist68,
    TTbarPU0EtaHist90, TTbarPU0PtHist90, TTbarPU0PhiHist90, TTbarPU0Z0Hist90};

    vector<TH1*> TTbarPU0Hists68 {TTbarPU0EtaHist68, TTbarPU0PtHist68, TTbarPU0PhiHist68, TTbarPU0Z0Hist68};
    vector<TH1*> TTbarPU0Hists90 {TTbarPU0EtaHist90, TTbarPU0PtHist90, TTbarPU0PhiHist90, TTbarPU0Z0Hist90};

    vector<vector<TH1*>> hists68 {muonHists68, electronHists68, TTbarPU200Hists68, TTbarPU0Hists68};
    vector<vector<TH1*>> hists90 {muonHists90, electronHists90, TTbarPU200Hists90, TTbarPU0Hists90};

    vector<TH1*> etaHists {muonEtaHist68, electronEtaHist68, TTbarPU200EtaHist68, TTbarPU0EtaHist68,
    muonEtaHist90, electronEtaHist90, TTbarPU200EtaHist90, TTbarPU0EtaHist90};

    vector<TH1*> ptHists {muonPtHist68, electronPtHist68, TTbarPU200PtHist68, TTbarPU0PtHist68,
    muonPtHist90, electronPtHist90, TTbarPU200PtHist90, TTbarPU0PtHist90};

    vector<TH1*> phiHists {muonPhiHist68, electronPhiHist68, TTbarPU200PhiHist68, TTbarPU0PhiHist68,
    muonPhiHist90, electronPhiHist90, TTbarPU200PhiHist90, TTbarPU0PhiHist90};

    vector<TH1*> z0Hists {muonZ0Hist68, electronZ0Hist68, TTbarPU200Z0Hist68, TTbarPU0Z0Hist68,
    muonZ0Hist90, electronZ0Hist90, TTbarPU200Z0Hist90, TTbarPU0Z0Hist90};
    

    for (int i = 0; i < allHists.size(); i++) {
      allHists[i]->Sumw2();
    }


    // Set colors of points and error lines
    int i_muon = muonHists.size();
    for (int i = 0; i < i_muon; i++) {
        muonHists[i]->SetMarkerColor(8);
        muonHists[i]->SetLineColor(8);
        electronHists[i]->SetMarkerColor(4);
        electronHists[i]->SetLineColor(4);
        TTbarPU200Hists[i]->SetMarkerColor(1);
        TTbarPU200Hists[i]->SetLineColor(1);
        TTbarPU0Hists[i]->SetMarkerColor(2);
        TTbarPU0Hists[i]->SetLineColor(2);
    }

    // Set shapes of plots
    for (int i = 0; i < muonHists68.size(); i++) {
        muonHists68[i]->SetMarkerStyle(kFullSquare);
        electronHists68[i]->SetMarkerStyle(kFullCircle);
        TTbarPU200Hists68[i]->SetMarkerStyle(kFullTriangleUp);
        TTbarPU0Hists68[i]->SetMarkerStyle(kFullDiamond);
        muonHists90[i]->SetMarkerStyle(kOpenSquare);
        electronHists90[i]->SetMarkerStyle(kOpenCircle);
        TTbarPU200Hists90[i]->SetMarkerStyle(kOpenTriangleUp);
        TTbarPU0Hists90[i]->SetMarkerStyle(kOpenDiamond);
    }

    // Draw and Print Histograms to pdf
    TString saveDir = "plotPDFs/";

    TCanvas c;
    double max;

    // eta
    TLegend* etaLeg = new TLegend(0.25, 0.70, 0.54, 0.90);
    etaLeg->AddEntry(muonEtaHist68, "Single Muon PU=0, Res=68%");
    etaLeg->AddEntry(TTbarPU200EtaHist68, "TTbar PU=200, Res=68%");
    etaLeg->AddEntry(TTbarPU0EtaHist68, "TTbar PU=0, Res=68%");
    etaLeg->AddEntry(electronEtaHist68, "Single Electron PU=0, Res=68%");
    etaLeg->AddEntry(muonEtaHist90, "Single Muon PU=0, Res=90%");
    etaLeg->AddEntry(TTbarPU200EtaHist90, "TTbar PU=200, Res=90%");
    etaLeg->AddEntry(TTbarPU0EtaHist90, "TTbar PU=0, Res=90%");
    etaLeg->AddEntry(electronEtaHist90, "Single Electron PU=0, Res=90%");

    etaHists[0]->Draw("pe");
    for (int i = 1; i < etaHists.size(); i++) {
        etaHists[i]->Draw("pe same");
    }

    for (int i = 0; i < etaHists.size(); i++) {
        int nbins = etaHists[i]->GetNbinsX();
        for (int ii = 0; ii < nbins; ii++) {
          etaHists[i]->SetBinError(ii, 0);
        }
    }

    etaLeg->Draw();
    max = GetMaxHists(etaHists);
    etaHists[0]->GetYaxis()->SetRangeUser(0, 0.01); // MANUAL FIX CAUSE GetMaximum() NOT WORKING
    c.SaveAs(saveDir + prop + "eta_60_90_MuonElectronTTbar.pdf");

    // pT
    TLegend* ptLeg = new TLegend(0.2, 0.68, 0.5, 0.88);
    ptLeg->AddEntry(muonPtHist68, "Single Muon PU=0, Res=68%");
    ptLeg->AddEntry(TTbarPU200PtHist68, "TTbar PU=200, Res=68%");
    ptLeg->AddEntry(TTbarPU0PtHist68, "TTbar PU=0, Res=68%");
    ptLeg->AddEntry(electronPtHist68, "Single Electron PU=0, Res=68%");
    ptLeg->AddEntry(muonPtHist90, "Single Muon PU=0, Res=90%");
    ptLeg->AddEntry(TTbarPU200PtHist90, "TTbar PU=200, Res=90%");
    ptLeg->AddEntry(TTbarPU0PtHist90, "TTbar PU=0, Res=90%");
    ptLeg->AddEntry(electronPtHist90, "Single Electron PU=0, Res=90%");

    ptHists[0]->Draw("pe");
    for (int i = 1; i < ptHists.size(); i++) {
        ptHists[i]->Draw("pe same");
    }

    for (int i = 0; i < ptHists.size(); i++) {
        int nbins = ptHists[i]->GetNbinsX();
        for (int ii = 0; ii < nbins; ii++) {
          ptHists[i]->SetBinError(ii, 0);
        }
    }

    ptLeg->Draw();
    max = GetMaxHists(ptHists);
    ptHists[0]->GetYaxis()->SetRangeUser(0, 1); // MANUAL FIX CAUSE GetMaximum() NOT WORKING
    c.SaveAs(saveDir + prop + "pt_60_90_MuonElectronTTbar.pdf");


    // phi
    TLegend* phiLeg = new TLegend(0.2, 0.68, 0.5, 0.88);
    phiLeg->AddEntry(muonPhiHist68, "Single Muon PU=0, Res=68%");
    phiLeg->AddEntry(TTbarPU200PhiHist68, "TTbar PU=200, Res=68%");
    phiLeg->AddEntry(TTbarPU0PhiHist68, "TTbar PU=0, Res=68%");
    phiLeg->AddEntry(electronPhiHist68, "Single Electron PU=0, Res=68%");
    phiLeg->AddEntry(muonPhiHist90, "Single Muon PU=0, Res=90%");
    phiLeg->AddEntry(TTbarPU200PhiHist90, "TTbar PU=200, Res=90%");
    phiLeg->AddEntry(TTbarPU0PhiHist90, "TTbar PU=0, Res=90%");
    phiLeg->AddEntry(electronPhiHist90, "Single Electron PU=0, Res=90%");

    phiHists[0]->Draw("pe");
    for (int i = 1; i < phiHists.size(); i++) {
        phiHists[i]->Draw("pe same");
    }

    for (int i = 0; i < phiHists.size(); i++) {
        int nbins = phiHists[i]->GetNbinsX();
        for (int ii = 0; ii < nbins; ii++) {
          phiHists[i]->SetBinError(ii, 0);
        }
    }

    phiLeg->Draw();
    max = GetMaxHists(phiHists);
    phiHists[0]->GetYaxis()->SetRangeUser(0, 0.016); // MANUAL FIX CAUSE GetMaximum() NOT WORKING
    c.SaveAs(saveDir + prop + "phi_60_90_MuonElectronTTbar.pdf");

    // // z0
    TLegend* z0Leg = new TLegend(0.2, 0.68, 0.5, 0.88);
    z0Leg->AddEntry(muonZ0Hist68, "Single Muon PU=0, Res=68%");
    z0Leg->AddEntry(TTbarPU200Z0Hist68, "TTbar PU=200, Res=68%");
    z0Leg->AddEntry(TTbarPU0Z0Hist68, "TTbar PU=0, Res=68%");
    z0Leg->AddEntry(electronZ0Hist68, "Single Electron PU=0, Res=68%");
    z0Leg->AddEntry(muonZ0Hist90, "Single Muon PU=0, Res=90%");
    z0Leg->AddEntry(TTbarPU200Z0Hist90, "TTbar PU=200, Res=90%");
    z0Leg->AddEntry(TTbarPU0Z0Hist90, "TTbar PU=0, Res=90%");
    z0Leg->AddEntry(electronZ0Hist90, "Single Electron PU=0, Res=90%");

    z0Hists[0]->Draw("pe");
    for (int i = 1; i < z0Hists.size(); i++) {
        z0Hists[i]->Draw("pe same");
    }

    for (int i = 0; i < z0Hists.size(); i++) {
        int nbins = z0Hists[i]->GetNbinsX();
        for (int ii = 0; ii < nbins; ii++) {
          z0Hists[i]->SetBinError(ii, 0);
        }
    }

    z0Leg->Draw();
    max = GetMaxHists(z0Hists);
    z0Hists[0]->GetYaxis()->SetRangeUser(0, 1.6); // MANUAL FIX CAUSE GetMaximum() NOT WORKING
    c.SaveAs(saveDir + prop + "z0_60_90_MuonElectronTTbar.pdf");

    // delete all pointers
    for (int i = 0; i < allHists.size(); i++) {
      delete allHists[i];
    }
    delete muonFile;
    delete electronFile;
    delete TTbarPU200File;
    delete TTbarPU0File;
    delete etaLeg, ptLeg, phiLeg, z0Leg;
}


double GetMaxHists(vector<TH1*> hists) {
  double maxValue = 0;
  double currMax = 0;

  for(int i=0; i<hists.size(); i++) {
    currMax = hists[i]->GetMaximum();

    if(currMax > maxValue) {
      maxValue = currMax;
    }
  }

  return maxValue;
}


void SetPlotStyle() {
  // from ATLAS plot style macro

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);

  gStyle->SetPalette(1);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20, 26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42, "x");
  gStyle->SetTitleFont(42, "x");
  gStyle->SetLabelFont(42, "y");
  gStyle->SetTitleFont(42, "y");
  gStyle->SetLabelFont(42, "z");
  gStyle->SetTitleFont(42, "z");
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");
  gStyle->SetTitleSize(0.05, "y");
  gStyle->SetLabelSize(0.05, "z");
  gStyle->SetTitleSize(0.05, "z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2, "[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}