//--------------------------------------------------------------------------------------------------------
// Opens root file in another directory ("dir") which contains histograms and overlays them.
// In this case, it overlays some property ("prop") of three different root files like "eff_eta" or something.
// Saves plots to a chosen directory "saveDir"
//
// ALSO it normalizes the plots that it overlays
//
// (SetPlotStyle function comes from Louise Skinnari's codes. It's an ATLAS plot style macro)
//
// By David Abrahamyan August 2023
//--------------------------------------------------------------------------------------------------------

#include "plotstyle.h"

//void SetPlotStyle(); // Sets plot style to some CMS thing I think idk Emily recommended I add this

double GetMaxHists(vector<TH1*> hists); // Gets maximum value of three histograms

void overlayHists_bendchi2 (){
    
    SetPlotStyle();

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";
    
    // Property you want to plot
    TString prop = "trk_bendchi2";

    // muon
    TFile *muonFile= new TFile(dir + "output_SingleMuon_PU0_D88.root");
    TH1F *muonHist= (TH1F*)muonFile->Get(prop);
    
    // electron
    TFile *electronFile= new TFile(dir + "output_SingleElectronPU0D88.root");
    TH1F *electronHist= (TH1F*)electronFile->Get(prop);

    // TTbarPU200
    TFile *TTbarPU200File= new TFile(dir + "output_TTbar_PU200_D88.root");
    TH1F *TTbarPU200Hist= (TH1F*)TTbarPU200File->Get(prop);

    // // TTbarPU0
    TFile *TTbarPU0File= new TFile(dir + "output_TTbar_PU0_D88.root");
    TH1F *TTbarPU0Hist= (TH1F*)TTbarPU0File->Get(prop);

    // Set colors of points and error lines
    muonHist->SetMarkerColor(8);
    muonHist->SetLineColor(8);
    electronHist->SetMarkerColor(4);
    electronHist->SetLineColor(4);
    TTbarPU200Hist->SetMarkerColor(1);
    TTbarPU200Hist->SetLineColor(1);
    TTbarPU0Hist->SetMarkerColor(2);
    TTbarPU0Hist->SetLineColor(2);
    muonHist->SetMarkerStyle(kFullSquare);
    TTbarPU200Hist->SetMarkerStyle(kFullTriangleUp);
    TTbarPU0Hist->SetMarkerStyle(kFullDiamond);

    // Draw and Print Histograms to pdf
    TString saveDir = "plotPDFs/";

    TCanvas c;
    double max;

    // Normalizing histos
    vector<TH1*> allHists{muonHist, electronHist, TTbarPU200Hist, TTbarPU0Hist};

    int histNum = allHists.size();
    for (int i = 0; i < histNum; i++) {
      allHists[i]->Scale(1./allHists[i]->Integral(), "width");
    }
    
    // Plot
    TLegend* leg = new TLegend(0.43, 0.4, 0.68, 0.6);
    leg->AddEntry(muonHist, "Single Muon PU=0");
    leg->AddEntry(TTbarPU200Hist, "TTbar PU=200");
    leg->AddEntry(TTbarPU0Hist, "TTbar PU=0");
    leg->AddEntry(electronHist, "Single Electron PU=0");

    muonHist->Draw();
    electronHist->Draw("same");
    TTbarPU200Hist->Draw("same");
    TTbarPU0Hist->Draw("same");
    leg->Draw();
    max = GetMaxHists(allHists);
    muonHist->GetYaxis()->SetRangeUser(0, max*1.1);
    muonHist->GetXaxis()->SetRangeUser(0, 30);
    c.SaveAs(saveDir + prop + "_MuonElectronTTbar.pdf");

    // ---------------Debug code-----------------
    
    // ------------------------------------------
}


double GetMaxHists(vector<TH1*> hists) {
  double maxValue = 0;

  for(int i=0; i<hists.size(); i++) {

    double currMax = hists[i]->GetMaximum();

    if(currMax > maxValue) {
      maxValue = currMax;
    }
  }

  return maxValue;
}


// void SetPlotStyle() {
//   // from ATLAS plot style macro

//   // use plain black on white colors
//   gStyle->SetFrameBorderMode(0);
//   gStyle->SetFrameFillColor(0);
//   gStyle->SetCanvasBorderMode(0);
//   gStyle->SetCanvasColor(0);
//   gStyle->SetPadBorderMode(0);
//   gStyle->SetPadColor(0);
//   gStyle->SetStatColor(0);
//   gStyle->SetHistLineColor(1);

//   gStyle->SetPalette(1);

//   // set the paper & margin sizes
//   gStyle->SetPaperSize(20, 26);
//   gStyle->SetPadTopMargin(0.05);
//   gStyle->SetPadRightMargin(0.05);
//   gStyle->SetPadBottomMargin(0.16);
//   gStyle->SetPadLeftMargin(0.16);

//   // set title offsets (for axis label)
//   gStyle->SetTitleXOffset(1.4);
//   gStyle->SetTitleYOffset(1.4);

//   // use large fonts
//   gStyle->SetTextFont(42);
//   gStyle->SetTextSize(0.05);
//   gStyle->SetLabelFont(42, "x");
//   gStyle->SetTitleFont(42, "x");
//   gStyle->SetLabelFont(42, "y");
//   gStyle->SetTitleFont(42, "y");
//   gStyle->SetLabelFont(42, "z");
//   gStyle->SetTitleFont(42, "z");
//   gStyle->SetLabelSize(0.05, "x");
//   gStyle->SetTitleSize(0.05, "x");
//   gStyle->SetLabelSize(0.05, "y");
//   gStyle->SetTitleSize(0.05, "y");
//   gStyle->SetLabelSize(0.05, "z");
//   gStyle->SetTitleSize(0.05, "z");

//   // use bold lines and markers
//   gStyle->SetMarkerStyle(20);
//   gStyle->SetMarkerSize(1.2);
//   gStyle->SetHistLineWidth(2.);
//   gStyle->SetLineStyleString(2, "[12 12]");

//   // get rid of error bar caps
//   gStyle->SetEndErrorSize(0.);

//   // do not display any of the standard histogram decorations
//   gStyle->SetOptTitle(0);
//   gStyle->SetOptStat(0);
//   gStyle->SetOptFit(0);

//   // put tick marks on top and RHS of plots
//   gStyle->SetPadTickX(1);
//   gStyle->SetPadTickY(1);
// }