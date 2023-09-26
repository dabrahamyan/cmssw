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

void overlayHists_eff_allTruncs_separate (){
    
    SetPlotStyle();

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";
    // Property you want to plot
    TString prop = "eff_";
    // directory to save plots to
    TString saveDir = "truncation_examination_assertsOn_oneTrunc/";
    

    /////////////////////// CYCLING THROUGH FILES ////////////////////
    // separate plots for algos
    vector <TString> algos = {"HYBRID", "NEWKF"};
    // separate plots for params within algos
    vector <TString> params = {"eta", "pt"};
    // plot all these on one plot
    vector <TString> truncs = {"no", "full", "placeholder"}; //, "VMRTrunc", "TETrunc", "TCTrunc", "PRTrunc"}; //"noTETrunc", "noTCTrunc", "noPRTrunc", "noMETrunc", "noMCTrunc"}; // "noMETrunc", "noMCTrunc"}; //"noMETrunc", "noTETrunc",
    //"noTCTrunc", "noPRTrunc", "noMETrunc", "noMCTrunc"};
    vector <TString> truncsToTest = {"TP", "MP"}; // "VMR", "TE", "TC", "PR", "ME", "MC"
    // What jet pT??
    vector <TString> pTJet = {"", "_injet", "_injet_highpt", "_injet_vhighpt"};

    // largest effect: PR, ME, MC

    // same as noTrunc no asserts: IR, VMR, DR, TB
    // same as noTrunc asserts on: IR, VMR, DR, TB

    /////////////////////// STYLIN ////////////////////////
    // Labels for plots
    vector <TString> algoLabels = {"Hybrid", "Hybrid_NewKF"};
    vector <TString> truncLabels = {"None", "All", "placeholder"};
    vector <TString> jetLabels = {"All Particles", "In Jets with p_{T} > 30 GeV", "In Jets with p_{T} > 100 GeV", "In Jets with p_{T} > 200 GeV"};
    // Markers n lines
    vector <short> colors = {kBlack, kAzure+7, kOrange-3};  //{1, 2, 38, 9, 6, 15, 28};
    vector <short> markers = {kFullSquare, kFullCircle, kFullDiamond}; //, 22, 25, 24, 26};

    TCanvas c;
    char ctxt[500];

    for (int iJetPt = 0; iJetPt < pTJet.size(); iJetPt++) {
        for (int iTruncsToTest = 0; iTruncsToTest < truncsToTest.size(); iTruncsToTest++) {
            truncs[2] = truncsToTest[iTruncsToTest];
            truncLabels[2] = truncsToTest[iTruncsToTest];
            // For two alogrithms (HYBRID and NEWKF)
            for (int iAlgo = 0; iAlgo < algos.size(); iAlgo++) {
                // open all necessary files
                vector <TFile*> files;
                for (int iTrunc = 0; iTrunc < truncs.size(); iTrunc++) {
                    files.push_back(new TFile(dir + "output_TTbar_PU200_D88_" + algos[iAlgo] + "_combined_" + truncs[iTrunc] + "Trunc_assertsOn_oneTrunc" + pTJet[iJetPt] + ".root"));
                }

                // for params (eta, pt, phi, z0)
                for (int iParam = 0; iParam < params.size(); iParam++) {
                    vector<TH1*> hists;
                    TLegend* leg = new TLegend(0.73, 0.2, 0.88, 0.44);

                    
                    // for different truncation settings
                    for (int iTrunc = 0; iTrunc < truncs.size(); iTrunc++) {
                        hists.push_back((TH1F*)files[iTrunc]->Get(prop + params[iParam]));
                        hists[iTrunc]->SetMarkerColor(colors[iTrunc]);
                        hists[iTrunc]->SetLineColor(colors[iTrunc]);
                        hists[iTrunc]->SetMarkerStyle(markers[iTrunc]);

                        leg->AddEntry(hists[iTrunc], truncLabels[iTrunc]);

                        if (iTrunc == 0) {
                        hists[iTrunc]->Draw();
                        }
                        else {
                        hists[iTrunc]->Draw("same");
                        }

                    }

                    
                    sprintf(ctxt, algoLabels[iAlgo] + " (" + jetLabels[iJetPt] + ")"); // Add label saying 
                    mySmallText(0.2, 0.3, 1, ctxt); // which data set it is

                    leg->SetTextSize(0.032);
                    leg->SetHeader("#splitline{Modules}{Truncated}", "C");
                    leg->Draw();
                    c.SaveAs(saveDir + "truncs_combined_" + truncsToTest[iTruncsToTest] + "_" + algos[iAlgo] + "_" + prop + params[iParam] + pTJet[iJetPt] + ".pdf");

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