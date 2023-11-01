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

void overlayHists_resVsEta_newkfVSnewkf (){
    
    SetPlotStyle();

    // DR version
    TString DRsetting = "z0";

    // Load in directory and files where the root files containing histograms are stored
    TString dir = "../LocalChecks/";
    TString OGNewkfFileName = "output_newkfDebug_SingleMuon_DR_off";
    TString newkfFileName = "output_newkfDebug_SingleMuon_DR_off_ownDR_" + DRsetting;

    // save settings
    TString saveDir = "hybrid_vs_newkf_plots/";
    TString saveFileStart = "newkfDubug_ownDR_" + DRsetting;
    
    // Property you want to plot
    TString prop = "resVsEta_";

    // data set to compare
    vector <TString> dataSets = {"SingleMuon_PU0_D88"};
    vector <TString> params = {"eta", "ptRel", "phi", "z0"};
    vector <TString> labels = {"Single Muon PU=0"};
    
    // doing this manually bc max does not work w resVsEta plots for some reason
    vector <float> maxVals = {0.035, 0.22, 0.0012, 2.4};

    TCanvas c;
    char ctxt[500];
    char ctxt2[500];
    char ctxt3[500];
    double max;

    int maxCounter = 0;


    // For data sets like SingleMuon, SingleElectron
    for (int iDataSet = 0; iDataSet < dataSets.size(); iDataSet++) {
        // Load OGNewkf and newkf root files
        TFile *OGNewkfFile= new TFile(dir + OGNewkfFileName + ".root");
        TFile *newkfFile= new TFile(dir + newkfFileName + ".root");

        // For parameters like eta, phi, ...
        for (int iParam = 0; iParam < params.size(); iParam++) {
            // copy OGNewkf and newkf hists for param
            TH1F *OGNewkfHist68= (TH1F*)OGNewkfFile->Get(prop + params[iParam] + "_68");
            TH1F *OGNewkfHist90= (TH1F*)OGNewkfFile->Get(prop + params[iParam] + "_90");
            TH1F *newkfHist68= (TH1F*)newkfFile->Get(prop + params[iParam] + "_68");
            TH1F *newkfHist90= (TH1F*)newkfFile->Get(prop + params[iParam] + "_90");

            // Set colors, markers, etc.
            newkfHist68->SetMarkerStyle(kOpenSquare);

            newkfHist90->SetMarkerStyle(kOpenCircle);

            OGNewkfHist90->SetMarkerStyle(kFullCircle);
            OGNewkfHist90->SetMarkerColor(kRed-3);

            OGNewkfHist68->SetMarkerColor(kRed-3);
            OGNewkfHist68->SetMarkerStyle(kFullSquare);

            // Set max height done manually bc max of resVsEta is being weird
            //vector<TH1*> histList {OGNewkfHist68, newkfHist68, OGNewkfHist90, newkfHist90};
            //max = GetMaxHists(histList);
            OGNewkfHist68->GetYaxis()->SetRangeUser(0, maxVals[maxCounter]);
            maxCounter++;

            // make legend
            TLegend* leg = new TLegend(0.2, 0.65, 0.55, 0.93);
            leg->AddEntry(OGNewkfHist68, "No special selection, Res=68%", "p"); 
            leg->AddEntry(newkfHist68, "Lowest " + DRsetting + " Residual, Res=68%", "p");
            leg->AddEntry(OGNewkfHist90, "No special selection, Res=90%", "p");  
            leg->AddEntry(newkfHist90, "Lowest " + DRsetting + " Residual, Res=90%", "p");
            
            // Draw and save histos
            OGNewkfHist68->Draw("p");
            OGNewkfHist90->Draw("p same");
            newkfHist68->Draw("p same");
            newkfHist90->Draw("p same");
            sprintf(ctxt, labels[iDataSet]); // Add label saying 
            mySmallText(0.2, 0.57, 1, ctxt); // which data set it is
            sprintf(ctxt2, "No DR");
            mySmallText(0.2, 0.52, 1, ctxt2);
            // sprintf(ctxt3, "Own DR eta");
            // mySmallText(0.2, 0.47, 1, ctxt3);
            leg->Draw();
            gStyle->SetOptStat(1);
            c.SaveAs(saveDir + saveFileStart + "_" + dataSets[iDataSet] + "_" + prop + params[iParam] + ".pdf");
            gStyle->SetOptStat(0);

            /////////////////////// DEBUG CODE ////////////////////////////////////
            // cout << "Entries in Old KF " + labels[iDataSet] + " 68%: " << OGNewkfHist68->GetEntries() << endl;
            // cout << "Entries in New KF " + labels[iDataSet] + " 68%: " << newkfHist68->GetEntries() << endl;
            // cout << "Entries in Old KF " + labels[iDataSet] + " 90%: " << OGNewkfHist90->GetEntries() << endl;
            // cout << "Entries in New KF " + labels[iDataSet] + " 90%: " << newkfHist90->GetEntries() << endl;
            ////////////////////////////////////////////////////////////////////////

            // delete pointers you want to remake
            delete OGNewkfHist68, newkfHist68, OGNewkfHist90, newkfHist90, leg;
        }

        // delete pointers you're gonna remake
        delete OGNewkfFile, newkfFile;
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