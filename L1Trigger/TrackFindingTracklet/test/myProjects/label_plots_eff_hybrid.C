//--------------------------------------------------------------------------------------------------------
// Opens root file in another directory ("dir") which contains histograms and overlays them.
// In this case, it overlays some property ("prop").
// It does this for sets of algorithms, parameters, truncations, etc that you can customize.
// Saves plots to a chosen directory "saveDir"
//
// (SetPlotStyle function comes from Louise Skinnari's codes. It's an ATLAS plot style macro)
//
// By David Abrahamyan September 2023
//--------------------------------------------------------------------------------------------------------

#include "plotstyle.h"

void mySmallText(Double_t x, Double_t y, Color_t color, char* text);

void label_plots_eff_hybrid (){
    
    SetPlotStyle();

    // filename
    TString fileName = "output_retryTest_TTbar_PU200_D88_HYBRID_Comb_LatestDev_2023_10_20";
    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";
    // Property you want to plot
    TString prop = "eff_";
    // directory to save plots to
    TString saveDir = "latestDev_2023_10_20/";
    

    /////////////////////// CYCLING THROUGH FILES ////////////////////
    // separate plots for params within algos
    vector <TString> params = {"eta", "pt", "phi", "z0"};
    
    TCanvas c;
    char ctxt[500];
    char ctxt2[500];

    TFile *file = new TFile(dir + fileName + ".root");

    for (int iParam = 0; iParam < params.size(); iParam++) {
        TH1F *hist= (TH1F*)file->Get(prop + params[iParam]);

        gPad->SetGridx();
        gPad->SetGridy();

        hist->Draw();

        sprintf(ctxt, "TTbar PU=200"); // Add label saying 
        mySmallText(0.5, 0.31, 1, ctxt); // which data set it is

        sprintf(ctxt2, "Hybrid Prompt Tracking"); // Add label saying 
        mySmallText(0.5, 0.36, 1, ctxt2); // which data set it is

        c.SaveAs(saveDir + "retryTest_TTbar_PU200_D88_HYBRID_Comb_eff_" + params[iParam] + ".pdf");

        delete hist;
    }
}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.043;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}