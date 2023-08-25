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

void overlayHists_eff_allJets (){
    
    SetPlotStyle();

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";
    // Property you want to plot
    TString prop = "eff_";
    // directory to save plots to
    TString saveDir = "truncation_examination/";

    /////////////// Which trunc do you want to test? ////////////////////
    TString trunc = "noIRTrunc";

    /////////////////////// CYCLING THROUGH FILES ////////////////////
    // separate plots for algos
    vector <TString> algos = {"HYBRID", "NEWKF"};
    // separate plots for params within algos
    vector <TString> params = {"eta", "pt", "phi", "z0"};
    // plot all these on one plot
    vector <TString> jets = {"","_injet","_injet_highpt","_injet_vhighpt"};
   
    /////////////////////// STYLIN ////////////////////////
    // Labels for plots
    vector <TString> algoLabels = {"Hybrid", "Hybrid_NewKF"};
    vector <TString> jetLabels = {"All Tracks", "In Jet (>30 GeV)", "In Jet (>100 GeV)", "In Jet (>200 GeV)"};
    // Markers n lines
    vector <int> colors = {1, 2, 3, 9};
    vector <int> markers = {21, 20, 22, 33};

    TCanvas c;
    char ctxt[500];

    // For two alogrithms (HYBRID and NEWKF)
    for (int iAlgo = 0; iAlgo < algos.size(); iAlgo++) {
      // open all necessary files
      vector <TFile*> files;
      for (int iJet = 0; iJet < jets.size(); iJet++) {
        files.push_back(new TFile(dir + "output_TTbar_PU200_D88_" + algos[iAlgo] + "_" + trunc + jets[iJet] + ".root"));
      }

      // for params (eta, pt, phi, z0)
      for (int iParam = 0; iParam < params.size(); iParam++) {
        vector<TH1*> hists;
        TLegend* leg = new TLegend(0.6, 0.3, 0.75, 0.4);

        // for different truncation settings
        for (int iJet = 0; iJet < jets.size(); iJet++) {
          hists.push_back((TH1F*)files[iJet]->Get(prop + params[iParam]));
          hists[iJet]->SetMarkerColor(colors[iJet]);
          hists[iJet]->SetLineColor(colors[iJet]);
          hists[iJet]->SetMarkerStyle(markers[iJet]);

          leg->AddEntry(hists[iJet], jetLabels[iJet]);

          if (iJet == 0) {
            hists[iJet]->Draw();
          }
          else {
            hists[iJet]->Draw("same");
          }

        }
        sprintf(ctxt, algoLabels[iAlgo]); // Add label saying 
        mySmallText(0.25, 0.3, 1, ctxt); // which data set it is

        leg->Draw();
        c.SaveAs(saveDir + "test_jets_" + algos[iAlgo] + "_" + prop + params[iParam] + ".pdf");

        delete leg;
        hists.clear();
      }

      // delete dynamically allocated memory
      for (int iJet = 0; iJet < jets.size(); iJet++) {
        delete files[iJet];
      }
      files.clear();
    }

}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.050;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}