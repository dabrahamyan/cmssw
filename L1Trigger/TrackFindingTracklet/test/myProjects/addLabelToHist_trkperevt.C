//--------------------------------------------------------------------------------------------------------
// Opens root file in another directory ("dir") which contains histogram and adds a label to it
// Saves plots to a chosen directory "saveDir"
//
// (SetPlotStyle function comes from Louise Skinnari's codes. It's an ATLAS plot style macro)
//
// By David Abrahamyan July 2023
//--------------------------------------------------------------------------------------------------------


void SetPlotStyle(); // Sets plot style to some CMS thing I think idk Emily recommended I add this

void mySmallText(Double_t x, Double_t y, Color_t color, char* text);

void addLabelToHist_trkperevt (){
    
    SetPlotStyle();

  

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";

    // focus on this file and then get histograms from it
    // muon
    TFile *muonFile= new TFile(dir + "output_SingleMuon_PU0_D88.root");
    TH1F *muonHist= (TH1F*)muonFile->Get("ntrk_tot");

    // electron
    TFile *electronFile= new TFile(dir + "output_SingleElectronPU0D88.root");
    TH1F *electronHist= (TH1F*)electronFile->Get("ntrk_tot");
    
    // TTbar
    TFile *TTbarFile= new TFile(dir + "output_TTbar_PU200_D88.root");
    TH1F *TTbarHist= (TH1F*)TTbarFile->Get("ntrk_tot");

    // Draw and Print Histograms to pdf
    TString saveDir = "plotPDFs/";

    TCanvas c;
    char ctxt[500];

    muonHist->Draw(); 
    sprintf(ctxt, "Single Muon PU0");
    mySmallText(0.45, 0.76, 1, ctxt);    
    c.SaveAs(saveDir + "muon_trkperevt.pdf");

    electronHist->Draw(); 
    sprintf(ctxt, "Single Electron PU0");
    mySmallText(0.45, 0.76, 1, ctxt);    
    c.SaveAs(saveDir + "electron_trkperevt.pdf");

    TTbarHist->Draw(); 
    sprintf(ctxt, "TTbar PU200");
    mySmallText(0.25, 0.76, 1, ctxt);    
    c.SaveAs(saveDir + "TTbar_trkperevt.pdf");

  
    // ---------------Debug code-----------------
    // muonEtaHist->Print("all");
    // electronEtaHist->Print("all");
    // TTbarEtaHist->Print();
    // ------------------------------------------
}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.050;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}

void SetPlotStyle() {
  // from ATLAS plot style macro

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);

  gStyle->SetPalette(1);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20, 26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42, "x");
  gStyle->SetTitleFont(42, "x");
  gStyle->SetLabelFont(42, "y");
  gStyle->SetTitleFont(42, "y");
  gStyle->SetLabelFont(42, "z");
  gStyle->SetTitleFont(42, "z");
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");
  gStyle->SetTitleSize(0.05, "y");
  gStyle->SetLabelSize(0.05, "z");
  gStyle->SetTitleSize(0.05, "z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2, "[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}