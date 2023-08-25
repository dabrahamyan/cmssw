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
    TString saveDir = "truncation_examination/";

    /////////////////////// CYCLING THROUGH FILES ////////////////////
    // separate plots for algos
    vector <TString> algos = {"HYBRID", "NEWKF"};
    // separate plots for params within algos
    vector <TString> params = {"eta", "pt", "phi", "z0"};
    // plot all these on one plot
    vector <TString> truncs = {"fullTrunc", "noTrunc", "noIRTrunc"}; //, "noVMRTrunc", "noTETrunc",
    //"noTCTrunc", "noPRTrunc", "noMETrunc", "noMCTrunc", "noTBTrunc", "noDRTrunc"};
   
    /////////////////////// STYLIN ////////////////////////
    // Labels for plots
    vector <TString> algoLabels = {"Hybrid", "Hybrid_NewKF"};
    // Markers n lines
    vector <int> colors = {1, 2, 3}; //, "4", "6", "8", "9", "15", "28", "41", "29"};
    vector <int> markers = {21, 20, 22};

    TCanvas c;
    char ctxt[500];

    // For two alogrithms (HYBRID and NEWKF)
    for (int iAlgo = 0; iAlgo < algos.size(); iAlgo++) {
      // open all necessary files
      vector <TFile*> files;
      for (int iTrunc = 0; iTrunc < truncs.size(); iTrunc++) {
        files.push_back(new TFile(dir + "output_TTbar_PU200_D88_" + algos[iAlgo] + "_" + truncs[iTrunc] + "_injet_vhighpt.root"));
      }

      // for params (eta, pt, phi, z0)
      for (int iParam = 0; iParam < params.size(); iParam++) {
        vector<TH1*> hists;
        TLegend* leg = new TLegend(0.6, 0.3, 0.75, 0.4);

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
        sprintf(ctxt, algoLabels[iAlgo]); // Add label saying 
        mySmallText(0.25, 0.3, 1, ctxt); // which data set it is

        leg->Draw();
        c.SaveAs(saveDir + "test_truncs_" + algos[iAlgo] + "_" + prop + params[iParam] + ".pdf");

        delete leg;
        hists.clear();
      }

      // delete dynamically allocated memory
      for (int iTrunc = 0; iTrunc < truncs.size(); iTrunc++) {
        delete files[iTrunc];
      }
      files.clear();
    }


    // // For data sets like SingleMuon, SingleElectron
    // for (int iDataSet = 0; iDataSet < dataSets.size(); iDataSet++) {
    //     // Load hybrid and newkf root files
    //     TFile *hybridFile= new TFile(dir + "output_TTbar_PU200_D88" + dataSets[iDataSet] + "_HYBRID.root");
    //     TFile *newkfFile= new TFile(dir + "output_" + dataSets[iDataSet] + "_NEWKF.root");

    //     // For parameters like eta, phi, ...
    //     for (int iParam = 0; iParam < params.size(); iParam++) {
    //         // copy hybrid and newkf hists for param
    //         TH1F *hybridHist= (TH1F*)hybridFile->Get(prop + params[iParam]);
    //         TH1F *newkfHist= (TH1F*)newkfFile->Get(prop + params[iParam]);

    //         // Set colors, markers, etc.
    //         hybridHist->SetMarkerColor(1);
    //         hybridHist->SetLineColor(1);
    //         newkfHist->SetMarkerColor(kGreen);
    //         newkfHist->SetLineColor(kGreen);
    //         newkfHist->SetMarkerStyle(kFullTriangleUp);

    //         // Set max height as a bit higher to fit in label
    //         hybridHist->GetYaxis()->SetRangeUser(0,1.3);

    //         // make legend
    //         TLegend* leg = new TLegend(0.5, 0.3, 0.65, 0.4);
    //         leg->AddEntry(hybridHist, "Hybrid");
    //         leg->AddEntry(newkfHist, "NewKF");
            

    //         // Draw and save histos
    //         hybridHist->Draw();
    //         newkfHist->Draw("same");
    //         sprintf(ctxt, labels[iDataSet]); // Add label saying 
    //         mySmallText(0.47, 0.83, 1, ctxt); // which data set it is
    //         leg->Draw();
    //         c.SaveAs(saveDir + "HYBRIDvsNEWKF" + "_" + prop + params[iParam] + "_" + dataSets[iDataSet] + ".pdf");

    //         // delete pointers you want to remake
    //         delete hybridHist, newkfHist, leg;
    //     }

    //     // delete pointers you're gonna remake
    //     delete hybridFile, newkfFile;
    // }

}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.050;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}