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

void overlayHists_combinedvsNot (){
    
    SetPlotStyle();

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";
    // Property you want to plot
    TString prop = "eff_";
    // directory to save plots to
    TString saveDir = "combined_vs_not/";
    // truncation
    TString trunc = "fullTrunc";

    /////////////////////// CYCLING THROUGH FILES ////////////////////
    // separate plots for algos
    vector <TString> algos = {"HYBRID"}; // "NEWKF"};
    // separate plots for params within algos
    vector <TString> params = {"eta", "pt"};
    // What jet pT??
    vector <TString> pTJet = {"", "_injet_vhighpt"};
    // combined or not
    vector <TString> combinedOrNo = {"_uncombined", "_combined"};

    // largest effect: PR, ME, MC

    // same as noTrunc no asserts: IR, VMR, DR, TB
    // same as noTrunc asserts on: IR, VMR, DR, TB

    /////////////////////// STYLIN ////////////////////////
    // Labels for plots
    vector <TString> algoLabels = {"Hybrid"}; // "Hybrid_NewKF"};
    vector <TString> combinedLabels = {"#splitline{Uncombined}{Modules}", "#splitline{Combined}{Modules}"};
    vector <TString> jetLabels = {"All Particles", "In Jets with p_{T} > 200 GeV"};
    // Markers n lines
    vector <short> colors = {kRed-3, kBlack};  //{1, 2, 38, 9, 6, 15, 28};
    vector <short> markers = {kFullCircle, kOpenCircle}; //, 22, 25, 24, 26};

    TCanvas c;
    char ctxt[500];

    for (int iJetPt = 0; iJetPt < pTJet.size(); iJetPt++) {
      // For two alogrithms (HYBRID and NEWKF)
      for (int iAlgo = 0; iAlgo < algos.size(); iAlgo++) {
        // open all necessary files
        vector <TFile*> files;
        for (int iComb = 0; iComb < combinedOrNo.size(); iComb++) {
          files.push_back(new TFile(dir + "output_SingleMuon_PU0_D88_" + algos[iAlgo] + combinedOrNo[iComb] + pTJet[iJetPt] + ".root"));
        }
        // output_SingleMuon_PU0_D88_HYBRID_combined.root
        // output_TTbar_PU200_D88_HYBRID_fullTrunc_assertsOn_oneTrunc.root
        // output_TTbar_PU200_D88_HYBRID_fullTrunc_assertsOn_oneTrunc_injet_highpt.root
        // output_TTbar_PU200_D88_NEWKF_combined_fullTrunc_assertsOn_oneTrunc.root

        // for params (eta, pt)
        for (int iParam = 0; iParam < params.size(); iParam++) {
          vector<TH1*> hists;
          TLegend* leg = new TLegend(0.72, 0.19, 0.89, 0.37);

          // for different truncation settings
          for (int iComb = 0; iComb < combinedOrNo.size(); iComb++) {
            hists.push_back((TH1F*)files[iComb]->Get(prop + params[iParam]));
            hists[iComb]->SetMarkerColor(colors[iComb]);
            hists[iComb]->SetLineColor(colors[iComb]);
            hists[iComb]->SetMarkerStyle(markers[iComb]);

            leg->AddEntry(hists[iComb], combinedLabels[iComb]);

            if (iComb == 0) {
              hists[iComb]->Draw();
            }
            else {
              hists[iComb]->Draw("same");
            }
          }

          sprintf(ctxt, algoLabels[iAlgo] + " (" + jetLabels[iJetPt] + ")"); // Add label saying 
          mySmallText(0.2, 0.3, 1, ctxt); // which data set it is
          
          leg->SetTextSize(0.032);
          //leg->SetHeader("#splitline{Modules}{Truncated}", "C");
          leg->Draw();
          c.SaveAs(saveDir + "Muon_combinedVsNot_" + algos[iAlgo] + "_" + prop + params[iParam] + pTJet[iJetPt] + ".pdf");

          delete leg;
          hists.clear();
        }

        // delete dynamically allocated memory
        for (int iComb = 0; iComb < combinedOrNo.size(); iComb++) {
          delete files[iComb];
        }
        files.clear();
      }
    }
}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.042;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}