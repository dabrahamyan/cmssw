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

//void SetPlotStyle(); // Sets plot style to some CMS thing I think idk Emily recommended I add this

void overlayHists_eff_hybridVSnewkf (){
    
    SetPlotStyle();

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";
    
    // Property you want to plot
    TString prop = "eff_";

    // directory to save plots to
    TString saveDir = "hybrid_vs_newkf_plots/";

    // data set to compare
    vector <TString> dataSets = {"SingleMuon_PU0_D88", "SingleElectron_PU0_D88", "TTbar_PU200_D88"};
    vector <TString> params = {"eta", "pt", "phi", "z0"};
    //vector <TString> leg?  
    vector <TString> labels = {"Single Muon PU=0", "Single Electron PU=0", "TTbar PU=200"};

    TCanvas c;
    char ctxt[500];

    // For data sets like SingleMuon, SingleElectron
    for (int iDataSet = 0; iDataSet < dataSets.size(); iDataSet++) {
        // Load hybrid and newkf root files
        TFile *hybridFile= new TFile(dir + "output_" + dataSets[iDataSet] + "_HYBRID.root");
        TFile *newkfFile= new TFile(dir + "output_hybridvsnewkf_" + dataSets[iDataSet] + "_NEWKF.root");

        // For parameters like eta, phi, ...
        for (int iParam = 0; iParam < params.size(); iParam++) {
            // copy hybrid and newkf hists for param
            TH1F *hybridHist= (TH1F*)hybridFile->Get(prop + params[iParam]);
            TH1F *newkfHist= (TH1F*)newkfFile->Get(prop + params[iParam]);

            // Set colors, markers, etc.
            hybridHist->SetMarkerColor(kBlack);
            hybridHist->SetLineColor(kBlack);
    
            newkfHist->SetMarkerColor(kRed-3);
            newkfHist->SetLineColor(kRed-3);
            newkfHist->SetMarkerStyle(kFullDiamond);

            // Set max height as a bit higher to fit in label
            hybridHist->GetYaxis()->SetRangeUser(0,1.3);

            // make legend
            TLegend* leg = new TLegend(0.5, 0.3, 0.65, 0.4);
            leg->AddEntry(hybridHist, "Hybrid");
            leg->AddEntry(newkfHist, "NewKF");
            

            // Draw and save histos
            hybridHist->Draw();
            newkfHist->Draw("same");
            sprintf(ctxt, labels[iDataSet]); // Add label saying 
            mySmallText(0.47, 0.83, 1, ctxt); // which data set it is
            leg->Draw();
            c.SaveAs(saveDir + "rerun_HYBRIDvsNEWKF" + "_" + prop + params[iParam] + "_" + dataSets[iDataSet] + ".pdf");

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