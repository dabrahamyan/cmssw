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

void overlayHists_oneParam_combinedvsNot (){
    
    SetPlotStyle();

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";
    
    // Property you want to plot
    TString prop = "trk_";

    // directory to save plots to
    TString saveDir = "truncation_examination_assertsOn_oneTrunc/";

    // data set
    TString dataUsed = "TTbar_PU200_D88_";

    // data set to compare
    vector <TString> dataSets = {"combined_", ""};
    vector <TString> params = {"eta"};
    vector <TString> labels = {"Combined Modules", "Not Combined Modules"};

    // output_TTbar_PU200_D88_HYBRID_fullTrunc_assertsOn_oneTrunc.root
    // output_TTbar_PU200_D88_HYBRID_combined_fullTrunc_assertsOn_oneTrunc.root

    // doing this manually bc max does not work w resVsEta plots for some reason
    // vector <float> maxVals = {0.035, 0.22, 0.002, 2.4, 0.035, 1.1, 0.025, 2.4, 0.035, 0.12, 0.006, 2.4};

    TCanvas c;
    char ctxt[500];
    double max;

    int maxCounter = 0;


    // For data sets like SingleMuon, SingleElectron
    for (int iDataSet = 0; iDataSet < dataSets.size(); iDataSet++) {
        // Load hybrid and newkf root files
        TFile *hybridFile= new TFile(dir + "output_" + dataUsed + "HYBRID_" + dataSets[iDataSet] + "fullTrunc_assertsOn_oneTrunc.root");
        TFile *newkfFile= new TFile(dir + "output_" + dataUsed + "HYBRID_" + dataSets[iDataSet] + "fullTrunc_assertsOn_oneTrunc.root");

        // For parameters like eta, phi, ...
        for (int iParam = 0; iParam < params.size(); iParam++) {
            // copy hybrid and newkf hists for param
            TH1F *hybridHist= (TH1F*)hybridFile->Get(prop + params[iParam]);
            TH1F *newkfHist= (TH1F*)newkfFile->Get(prop + params[iParam]);


            // Set colors, markers, etc.
            newkfHist->SetMarkerColor(kRed-3);
            newkfHist->SetLineColor(kRed-3);
            newkfHist->SetMarkerStyle(kFullDiamond);


            // Set max height done manually bc max of resVsEta is being weird
            //vector<TH1*> histList {hybridHist68, newkfHist68, hybridHist90, newkfHist90};
            //max = GetMaxHists(histList);
            //hybridHist68->GetYaxis()->SetRangeUser(0, maxVals[maxCounter]);
            maxCounter++;

            // make legend
            TLegend* leg = new TLegend(0.23, 0.25, 0.35, 0.4);
            leg->AddEntry(hybridHist, "Hybrid");
            leg->AddEntry(newkfHist, "NewKF");
            
            // Draw and save histos
            hybridHist->Draw("p");
            newkfHist->Draw("p same");

            sprintf(ctxt, labels[iDataSet]); // Add label saying 
            mySmallText(0.47, 0.25, 1, ctxt); // which data set it is
            leg->Draw();
            c.SaveAs(saveDir + "combinedvsnot" + "_" + prop + params[iParam] + "_" + dataUsed + "HYBRID.pdf");

            /////////////////////// DEBUG CODE ////////////////////////////////////
            cout << "Entries in Old KF " + labels[iDataSet] + ": " << hybridHist->GetEntries() << endl;
            cout << "Entries in New KF " + labels[iDataSet] + ": " << newkfHist->GetEntries() << endl;
            ////////////////////////////////////////////////////////////////////////

            // delete pointers you want to remake
            delete hybridHist, newkfHist, leg;
        }

        // delete pointers you're gonna remake
        delete hybridFile, newkfFile;
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

