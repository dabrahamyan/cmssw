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

void overlayHists_eff_allTruncs (){
    
    SetPlotStyle();

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";
    // Property you want to plot
    TString prop = "eff_";
    // directory to save plots to
    TString saveDir = "truncation_examination_assertsOn/";

    /////////////////////// CYCLING THROUGH FILES ////////////////////
    // separate plots for algos
    vector <TString> algos = {"HYBRID", "NEWKF"};
    // separate plots for params within algos
    vector <TString> params = {"eta", "pt", "phi", "z0"};
    // plot all these on one plot
    vector <TString> combined = {"fullTrunc", "noTrunc", "noMCTrunc"}; 


    /////////////////////// STYLIN ////////////////////////
    // Labels for plots
    vector <TString> algoLabels = {"Hybrid", "Hybrid_NewKF"};
    // Markers n lines
    vector <int> colors = {1, 2, 30, 9, 6, 15, 28};
    vector <int> markers = {21, 20, 33, 22, 25, 24, 26};

    TCanvas c;
    char ctxt[500];

    // For two alogrithms (HYBRID and NEWKF)
    for (int iAlgo = 0; iAlgo < algos.size(); iAlgo++) {
      // open all necessary files
      vector <TFile*> files;
      for (int iTrunc = 0; iTrunc < truncs.size(); iTrunc++) {
        files.push_back(new TFile(dir + "output_TTbar_PU200_D88_" + algos[iAlgo] + "_" + truncs[iTrunc] + "_assertsOn_injet_vhighpt.root"));
      }

      // for params (eta, pt, phi, z0)
      for (int iParam = 0; iParam < params.size(); iParam++) {
        vector<TH1*> hists;
        TLegend* leg = new TLegend(0.65, 0.2, 0.9, 0.45);

        // for different truncation settings
        for (int iTrunc = 0; iTrunc < truncs.size(); iTrunc++) {
          hists.push_back((TH1F*)files[iTrunc]->Get(prop + params[iParam]));
          hists[iTrunc]->SetMarkerColor(colors[iTrunc]);
          hists[iTrunc]->SetLineColor(colors[iTrunc]);
          hists[iTrunc]->SetMarkerStyle(markers[iTrunc]);

          leg->AddEntry(hists[iTrunc], truncs[iTrunc]);

          if (iTrunc == 0) {
            hists[iTrunc]->Draw();
          }
          else {
            hists[iTrunc]->Draw("same");
          }

        }


        sprintf(ctxt, algoLabels[iAlgo] + " (p_{T} Jet>200 GeV)"); // Add label saying 
        mySmallText(0.2, 0.3, 1, ctxt); // which data set it is

        leg->Draw();
        c.SaveAs(saveDir + "combinedVsNotCombined_" + algos[iAlgo] + "_" + prop + params[iParam] + "_assertsOn.pdf");

        delete leg;
        hists.clear();
      }

      // delete dynamically allocated memory
      for (int iTrunc = 0; iTrunc < truncs.size(); iTrunc++) {
        delete files[iTrunc];
      }
      files.clear();
    }

}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.045;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}