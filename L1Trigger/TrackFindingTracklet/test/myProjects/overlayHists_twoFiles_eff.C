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

void overlayHists_twoFiles_eff (){
    
    SetPlotStyle();

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/DisplacedCombinedBugFix/";
    // Property you want to plot
    TString prop = "eff_";
    // directory to save plots to
    TString saveDir = "displaced_combined_vs_not/";

    /////////////////////// CYCLING THROUGH FILES ////////////////////
    // separate plots for params within algos
    vector <TString> params = {"eta", "pt", "phi", "z0", "absd0", "absd0_eta2"};
    // combined or not
    vector <TString> fileNames = {"output_DisplacedMuon_PU0_D88_DISPLACED_Uncomb_BugFix_2023_11_8", "output_DisplacedMuon_PU0_D88_DISPLACED_Comb_BugFix_2023_11_8"};

    /////////////////////// STYLIN ////////////////////////
    // Labels for plots
    TString algoLabel = "Hybrid_Displaced"; // "Hybrid_NewKF"};
    TString dataLabel = "Displaced Muon PU=0";
    vector <TString> legLabels = {"#splitline{Uncombined}{Modules}", "#splitline{Combined}{Modules}"};
    // Markers n lines
    vector <short> colors = {kRed-3, kBlack};  //{1, 2, 38, 9, 6, 15, 28};
    vector <short> markers = {kFullCircle, kOpenCircle}; //, 22, 25, 24, 26};

    TCanvas c;
    char ctxt[500];
    char ctxt2[500];
    char ctxt3[500];

        // open all necessary files
        vector <TFile*> files;
        for (int iFile = 0; iFile < fileNames.size(); iFile++) {
          files.push_back(new TFile(dir + fileNames[iFile] + ".root"));
        }

        // for params (eta, pt)
        for (int iParam = 0; iParam < params.size(); iParam++) {
          vector<TH1*> hists;
          TLegend* leg = new TLegend(0.72, 0.19, 0.89, 0.37);

          // for different files
          for (int iFile = 0; iFile < fileNames.size(); iFile++) {
            hists.push_back((TH1F*)files[iFile]->Get(prop + params[iParam]));
            hists[iFile]->SetMarkerColor(colors[iFile]);
            hists[iFile]->SetLineColor(colors[iFile]);
            hists[iFile]->SetMarkerStyle(markers[iFile]);

            leg->AddEntry(hists[iFile], legLabels[iFile]);

            if (iFile == 0) {
              hists[iFile]->Draw();
            }
            else {
              hists[iFile]->Draw("same");
            }
          }

          sprintf(ctxt, algoLabel); // Add label
          mySmallText(0.2, 0.3, 1, ctxt);
          sprintf(ctxt2, dataLabel); // Add label
          mySmallText(0.2, 0.35, 1, ctxt2);
          if (params[iParam] == "absd0_eta2") {
            sprintf(ctxt3, "#eta < 2.0"); // Add label
            mySmallText(0.2, 0.25, 1, ctxt3);
          }

          leg->SetTextSize(0.032);
          leg->Draw();
          c.SaveAs(saveDir + "DisplacedMuon_combVsNot_bugFix_2023_11_8_" + prop + params[iParam] + ".pdf");

          delete leg;
          hists.clear();
        }

        // delete dynamically allocated memory
        for (int iFile = 0; iFile < fileNames.size(); iFile++) {
          delete files[iFile];
        }
        files.clear();
      
    
}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.042;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}