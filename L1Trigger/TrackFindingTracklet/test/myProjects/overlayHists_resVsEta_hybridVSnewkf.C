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

//void SetPlotStyle(); // Sets plot style to some CMS thing I think idk Emily recommended I add this

void overlayHists_resVsEta_hybridVSnewkf (){
    
    SetPlotStyle();

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";
    
    // Property you want to plot
    TString prop = "resVsEta_";

    // directory to save plots to
    TString saveDir = "hybrid_vs_newkf_plots/";

    // data set to compare
    vector <TString> dataSets = {"SingleMuon_PU0_D88", "SingleElectron_PU0_D88", "TTbar_PU200_D88"};
    vector <TString> params = {"eta", "ptRel", "phi", "z0"};
    vector <TString> labels = {"Single Muon PU=0", "Single Electron PU=0", "TTbar PU=200"};
    
    // doing this manually bc max does not work w resVsEta plots for some reason
    vector <float> maxVals = {0.035, 0.22, 0.002, 2.4, 0.035, 1.1, 0.025, 2.4, 0.035, 0.12, 0.006, 2.4};

    TCanvas c;
    char ctxt[500];
    double max;

    int maxCounter = 0;


    // For data sets like SingleMuon, SingleElectron
    for (int iDataSet = 0; iDataSet < dataSets.size(); iDataSet++) {
        // Load hybrid and newkf root files
        TFile *hybridFile= new TFile(dir + "output_" + dataSets[iDataSet] + "_HYBRID.root");
        TFile *newkfFile= new TFile(dir + "output_" + dataSets[iDataSet] + "_NEWKF.root");

        // For parameters like eta, phi, ...
        for (int iParam = 0; iParam < params.size(); iParam++) {
            // copy hybrid and newkf hists for param
            TH1F *hybridHist68= (TH1F*)hybridFile->Get(prop + params[iParam] + "_68");
            TH1F *newkfHist68= (TH1F*)newkfFile->Get(prop + params[iParam] + "_68");
            TH1F *hybridHist90= (TH1F*)hybridFile->Get(prop + params[iParam] + "_90");
            TH1F *newkfHist90= (TH1F*)newkfFile->Get(prop + params[iParam] + "_90");

            // Set colors, markers, etc.
            newkfHist68->SetMarkerColor(kRed-3);
            newkfHist68->SetLineColor(kRed-3);
            newkfHist68->SetMarkerStyle(kFullDiamond);

            newkfHist90->SetMarkerColor(kRed-3);
            newkfHist90->SetLineColor(kRed-3);
            newkfHist90->SetMarkerStyle(kOpenDiamond);
            hybridHist90->SetMarkerStyle(kOpenCircle);

            // Set max height done manually bc max of resVsEta is being weird
            //vector<TH1*> histList {hybridHist68, newkfHist68, hybridHist90, newkfHist90};
            //max = GetMaxHists(histList);
            hybridHist68->GetYaxis()->SetRangeUser(0, maxVals[maxCounter]);
            maxCounter++;

            // make legend
            TLegend* leg = new TLegend(0.2, 0.65, 0.45, 0.85);
            leg->AddEntry(hybridHist68, "Hybrid, Res=68%"); 
            leg->AddEntry(newkfHist68, "NewKF, Res=68%");
            leg->AddEntry(hybridHist90, "Hybrid, Res=90%");  
            leg->AddEntry(newkfHist90, "NewKF, Res=90%");
            
            // Draw and save histos
            hybridHist68->Draw("p");
            newkfHist68->Draw("p same");
            hybridHist90->Draw("p same");
            newkfHist90->Draw("p same");
            sprintf(ctxt, labels[iDataSet]); // Add label saying 
            mySmallText(0.47, 0.85, 1, ctxt); // which data set it is
            leg->Draw();
            gStyle->SetOptStat(1);
            c.SaveAs(saveDir + "HYBRIDvsNEWKF" + "_" + prop + params[iParam] + "_" + dataSets[iDataSet] + ".pdf");
            gStyle->SetOptStat(0);

            /////////////////////// DEBUG CODE ////////////////////////////////////
            cout << "Entries in Old KF " + labels[iDataSet] + " 68%: " << hybridHist68->GetEntries() << endl;
            cout << "Entries in New KF " + labels[iDataSet] + " 68%: " << newkfHist68->GetEntries() << endl;
            cout << "Entries in Old KF " + labels[iDataSet] + " 90%: " << hybridHist90->GetEntries() << endl;
            cout << "Entries in New KF " + labels[iDataSet] + " 90%: " << newkfHist90->GetEntries() << endl;
            ////////////////////////////////////////////////////////////////////////

            // delete pointers you want to remake
            delete hybridHist68, newkfHist68, hybridHist90, newkfHist90, leg;
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