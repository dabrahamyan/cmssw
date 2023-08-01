//--------------------------------------------------------------------------------------------------------
// Opens root file in another directory ("dir") which contains histograms and overlays them.
// In this case, it overlays some property ("prop") of three different root files like "eff_eta" or something.
// Saves plots to a chosen directory "saveDir"
//
// (SetPlotStyle function comes from Louise Skinnari's codes. It's an ATLAS plot style macro)
//
// By David Abrahamyan July 2023
//--------------------------------------------------------------------------------------------------------


void SetPlotStyle(); // Sets plot style to some CMS thing I think idk Emily recommended I add this

double GetMaxHists(vector<TH1*> hists); // Gets maximum value of three histograms

void overlayHists_resVsEta (){
    
    SetPlotStyle();

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";
    
    // Property you want to plot
    TString prop = "resVsEta_";

    // focus on this file and then get histograms from it
    TFile *muonFile= new TFile(dir + "output_SingleMuon_PU0_D88.root");
    TH1F *muonEtaHist= (TH1F*)muonFile->Get(prop + "eta");
    TH1F *muonPtHist= (TH1F*)muonFile->Get(prop + "ptRel");
    TH1F *muonPhiHist= (TH1F*)muonFile->Get(prop + "phi");
    TH1F *muonZ0Hist= (TH1F*)muonFile->Get(prop + "z0");

    //----------Test Code-------------
    //muonFile->ls();
    //eff_eta->Draw();
    //muonEtaHist->Draw();
    //muonEtaHist->Print();
    //--------------------------------

    // rinse and repeat for two more root files
    
    // electron
    TFile *electronFile= new TFile(dir + "output_SingleElectronPU0D88.root");
    TH1F *electronEtaHist= (TH1F*)electronFile->Get(prop + "eta");
    TH1F *electronPtHist= (TH1F*)electronFile->Get(prop + "ptRel");
    TH1F *electronPhiHist= (TH1F*)electronFile->Get(prop + "phi");
    TH1F *electronZ0Hist= (TH1F*)electronFile->Get(prop + "z0");

    // TTbar
    TFile *TTbarFile= new TFile(dir + "output_TTbar_PU200_D88.root");
    TH1F *TTbarEtaHist= (TH1F*)TTbarFile->Get(prop + "eta");
    TH1F *TTbarPtHist= (TH1F*)TTbarFile->Get(prop + "ptRel");
    TH1F *TTbarPhiHist= (TH1F*)TTbarFile->Get(prop + "phi");
    TH1F *TTbarZ0Hist= (TH1F*)TTbarFile->Get(prop + "z0");

    // Set colors of points and error lines
    muonEtaHist->SetMarkerColor(8);
    muonEtaHist->SetLineColor(8);
    electronEtaHist->SetMarkerColor(4);
    electronEtaHist->SetLineColor(4);
    TTbarEtaHist->SetMarkerColor(1);
    TTbarEtaHist->SetLineColor(1);
    muonEtaHist->SetMarkerStyle(kFullSquare);
    TTbarEtaHist->SetMarkerStyle(kFullTriangleUp);

    muonPtHist->SetMarkerColor(8);
    muonPtHist->SetLineColor(8);
    electronPtHist->SetMarkerColor(4);
    electronPtHist->SetLineColor(4);
    TTbarPtHist->SetMarkerColor(1);
    TTbarPtHist->SetLineColor(1);
    muonPtHist->SetMarkerStyle(kFullSquare);
    TTbarPtHist->SetMarkerStyle(kFullTriangleUp);

    muonPhiHist->SetMarkerColor(8);
    muonPhiHist->SetLineColor(8);
    electronPhiHist->SetMarkerColor(4);
    electronPhiHist->SetLineColor(4);
    TTbarPhiHist->SetMarkerColor(1);
    TTbarPhiHist->SetLineColor(1);
    muonPhiHist->SetMarkerStyle(kFullSquare);
    TTbarPhiHist->SetMarkerStyle(kFullTriangleUp);

    muonZ0Hist->SetMarkerColor(8);
    muonZ0Hist->SetLineColor(8);
    electronZ0Hist->SetMarkerColor(4);
    electronZ0Hist->SetLineColor(4);
    TTbarZ0Hist->SetMarkerColor(1);
    TTbarZ0Hist->SetLineColor(1);
    muonZ0Hist->SetMarkerStyle(kFullSquare);
    TTbarZ0Hist->SetMarkerStyle(kFullTriangleUp);

    // Draw and Print Histograms to pdf
    TString saveDir = "plotPDFs/";

    TCanvas c;
    double max;

    // eta
    TLegend* etaLeg = new TLegend(0.6, 0.25, 0.89, 0.45);
    etaLeg->AddEntry(muonEtaHist, "Single Muon PU0");
    etaLeg->AddEntry(TTbarEtaHist, "TTbar PU200");
    etaLeg->AddEntry(electronEtaHist, "Single Electron PU0");

    muonEtaHist->Draw();
    electronEtaHist->Draw("same");
    TTbarEtaHist->Draw("same");
    etaLeg->Draw();
    vector<TH1*> etaHists{muonEtaHist, electronEtaHist, TTbarEtaHist};
    max = GetMaxHists(etaHists);
    muonEtaHist->GetYaxis()->SetRangeUser(0, max*1.1);
    c.SaveAs(saveDir + prop + "eta_MuonElectronTTbar.pdf");

    // pT
    TLegend* ptLeg = new TLegend(0.4, 0.45, 0.69, 0.65);
    ptLeg->AddEntry(muonEtaHist, "Single Muon PU0");
    ptLeg->AddEntry(TTbarEtaHist, "TTbar PU200");
    ptLeg->AddEntry(electronEtaHist, "Single Electron PU0");

    muonPtHist->Draw();
    electronPtHist->Draw("same");
    TTbarPtHist->Draw("same");
    ptLeg->Draw();
    vector<TH1*> ptHists{muonPtHist, electronPtHist, TTbarPtHist};
    max = GetMaxHists(ptHists);
    muonPtHist->GetYaxis()->SetRangeUser(0, max*1.1);
    c.SaveAs(saveDir + prop + "ptRel_MuonElectronTTbar.pdf");

    // phi
    TLegend* phiLeg = new TLegend(0.25, 0.3, 0.55, 0.5);
    phiLeg->AddEntry(muonEtaHist, "Single Muon PU0");
    phiLeg->AddEntry(TTbarEtaHist, "TTbar PU200");
    phiLeg->AddEntry(electronEtaHist, "Single Electron PU0");

    muonPhiHist->Draw();
    electronPhiHist->Draw("same");
    TTbarPhiHist->Draw("same");
    phiLeg->Draw();
    vector<TH1*> phiHists{muonPhiHist, electronPhiHist, TTbarPhiHist};
    max = GetMaxHists(phiHists);
    muonPhiHist->GetYaxis()->SetRangeUser(0, max*1.1);
    c.SaveAs(saveDir + prop + "phi_MuonElectronTTbar.pdf");

    // z0
    TLegend* z0Leg = new TLegend(0.25, 0.65, 0.55, 0.85);
    z0Leg->AddEntry(muonEtaHist, "Single Muon PU0");
    z0Leg->AddEntry(TTbarEtaHist, "TTbar PU200");
    z0Leg->AddEntry(electronEtaHist, "Single Electron PU0");

    muonZ0Hist->Draw();
    electronZ0Hist->Draw("same");
    TTbarZ0Hist->Draw("same");
    z0Leg->Draw();
    vector<TH1*> z0Hists{muonZ0Hist, electronZ0Hist, TTbarZ0Hist};
    max = GetMaxHists(z0Hists);
    muonZ0Hist->GetYaxis()->SetRangeUser(0, max*1.1);
    c.SaveAs(saveDir + prop + "z0_MuonElectronTTbar.pdf");

    // ---------------Debug code-----------------
    // muonEtaHist->Print("all");
    // electronEtaHist->Print("all");
    // TTbarEtaHist->Print();
    // ------------------------------------------
}


double GetMaxHists(vector<TH1*> hists) {
  double maxValue = 0;

  for(int i=0; i<hists.size(); i++) {

    double currMax = hists[i]->GetMaximum();

    if(currMax > maxValue) {
      maxValue = currMax;
    }
  }

  return maxValue;
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