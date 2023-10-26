//--------------------------------------------------------------------------------------------------------
// Opens root file in another directory ("dir") which contains histograms and overlays them.
// In this case, it overlays some property ("prop") of three different root files like "eff_eta" or something.
// Saves plots to a chosen directory "saveDir"
//
// (SetPlotStyle function comes from Louise Skinnari's codes. It's an ATLAS plot style macro)
//
// By David Abrahamyan August 2023
//--------------------------------------------------------------------------------------------------------

#include "plotstyle.h"

void mySmallText(Double_t x, Double_t y, Color_t color, char* text);

double GetMaxHists(vector<TH1*> hists);

void overlayTwoHists_pT_superLow (){
    
    SetPlotStyle();

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";

    // directory to save plots to
    TString saveDir = "hybrid_vs_newkf_plots/";

    // parameter to compare
    TString param = "eff_pt_superLow";

    // files to compare
    TString file1 = "output_hybridvsnewkf_TTbar_PU200_D88_HYBRID.root";
    TString file2 = "output_hybridvsnewkf_TTbar_PU200_D88_NEWKF.root";

    // LABELS //
    TString dataSetLabel = "TTbar PU=200";

    TCanvas c;
    char ctxt[500];

    TFile *hybridFile= new TFile(dir + file1);
    TFile *newkfFile= new TFile(dir + file2);

    TH1F *hybridHist= (TH1F*)hybridFile->Get(param);
    TH1F *newkfHist= (TH1F*)newkfFile->Get(param);

    // Set colors, markers, etc.
    hybridHist->SetMarkerColor(kRed-3);
    hybridHist->SetLineColor(kRed-3);
    hybridHist->SetMarkerStyle(kFullCircle);
    newkfHist->SetMarkerStyle(kOpenCircle);

    TLegend* leg = new TLegend(0.23, 0.25, 0.35, 0.4);
    leg->AddEntry(hybridHist, "Hybrid");
    leg->AddEntry(newkfHist, "NewKF");
            
    // Draw and save histos
    hybridHist->Draw("p");
    newkfHist->Draw("p same");

    sprintf(ctxt, dataSetLabel); // Add label saying 
    mySmallText(0.47, 0.25, 1, ctxt); // which data set it is
    leg->Draw();
    c.SaveAs(saveDir + "eff_pt_superLow.pdf");

    /////////////////////// DEBUG CODE ////////////////////////////////////
    // cout << "Entries in Old KF " + labels[iDataSet] + ": " << hybridHist->GetEntries() << endl;
    // cout << "Entries in New KF " + labels[iDataSet] + ": " << newkfHist->GetEntries() << endl;
    ////////////////////////////////////////////////////////////////////////
}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.050;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
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

