//--------------------------------------------------------------------------------------------------------
// Opens root file in another directory ("dir") which contains histograms and overlays them.
// In this case, it overlays some property ("prop") of three different root files like "eff_eta" or something.
// Saves plots to a chosen directory "saveDir"
//
// ALSO it normalizes the plots that it overlays
//
// (SetPlotStyle function comes from Louise Skinnari's codes. It's an ATLAS plot style macro)
//
// By David Abrahamyan August 2023
//--------------------------------------------------------------------------------------------------------


void SetPlotStyle(); // Sets plot style to some CMS thing I think idk Emily recommended I add this

double GetMaxHists(vector<TH1*> hists); // Gets maximum value of three histograms

void overlayHists_trk_params (){
    
    SetPlotStyle();

    // Load in directory where the root files containing histograms are stored
    TString dir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";
    
    // Property you want to plot
    TString prop = "trk_";

    // muon
    TFile *muonFile= new TFile(dir + "output_SingleMuon_PU0_D88.root");
    TH1F *muonEtaHist= (TH1F*)muonFile->Get(prop + "eta");
    TH1F *muonPtHist= (TH1F*)muonFile->Get(prop + "pt");
    TH1F *muonPhiHist= (TH1F*)muonFile->Get(prop + "phi");
    TH1F *muonZ0Hist= (TH1F*)muonFile->Get(prop + "z0");
    
    // electron
    TFile *electronFile= new TFile(dir + "output_SingleElectronPU0D88.root");
    TH1F *electronEtaHist= (TH1F*)electronFile->Get(prop + "eta");
    TH1F *electronPtHist= (TH1F*)electronFile->Get(prop + "pt");
    TH1F *electronPhiHist= (TH1F*)electronFile->Get(prop + "phi");
    TH1F *electronZ0Hist= (TH1F*)electronFile->Get(prop + "z0");

    // TTbarPU200
    TFile *TTbarPU200File= new TFile(dir + "output_TTbar_PU200_D88.root");
    TH1F *TTbarPU200EtaHist= (TH1F*)TTbarPU200File->Get(prop + "eta");
    TH1F *TTbarPU200PtHist= (TH1F*)TTbarPU200File->Get(prop + "pt");
    TH1F *TTbarPU200PhiHist= (TH1F*)TTbarPU200File->Get(prop + "phi");
    TH1F *TTbarPU200Z0Hist= (TH1F*)TTbarPU200File->Get(prop + "z0");

    // // TTbarPU0
    TFile *TTbarPU0File= new TFile(dir + "output_TTbar_PU0_D88.root");
    TH1F *TTbarPU0EtaHist= (TH1F*)TTbarPU0File->Get(prop + "eta");
    TH1F *TTbarPU0PtHist= (TH1F*)TTbarPU0File->Get(prop + "pt");
    TH1F *TTbarPU0PhiHist= (TH1F*)TTbarPU0File->Get(prop + "phi");
    TH1F *TTbarPU0Z0Hist= (TH1F*)TTbarPU0File->Get(prop + "z0");

    // Set colors of points and error lines
    muonEtaHist->SetMarkerColor(8);
    muonEtaHist->SetLineColor(8);
    electronEtaHist->SetMarkerColor(4);
    electronEtaHist->SetLineColor(4);
    TTbarPU200EtaHist->SetMarkerColor(1);
    TTbarPU200EtaHist->SetLineColor(1);
    TTbarPU0EtaHist->SetMarkerColor(2);
    TTbarPU0EtaHist->SetLineColor(2);
    muonEtaHist->SetMarkerStyle(kFullSquare);
    TTbarPU200EtaHist->SetMarkerStyle(kFullTriangleUp);
    TTbarPU0EtaHist->SetMarkerStyle(kFullDiamond);

    muonPtHist->SetMarkerColor(8);
    muonPtHist->SetLineColor(8);
    electronPtHist->SetMarkerColor(4);
    electronPtHist->SetLineColor(4);
    TTbarPU200PtHist->SetMarkerColor(1);
    TTbarPU200PtHist->SetLineColor(1);
    TTbarPU0PtHist->SetMarkerColor(2);
    TTbarPU0PtHist->SetLineColor(2);
    muonPtHist->SetMarkerStyle(kFullSquare);
    TTbarPU200PtHist->SetMarkerStyle(kFullTriangleUp);
    TTbarPU0PtHist->SetMarkerStyle(kFullDiamond);

    muonPhiHist->SetMarkerColor(8);
    muonPhiHist->SetLineColor(8);
    electronPhiHist->SetMarkerColor(4);
    electronPhiHist->SetLineColor(4);
    TTbarPU200PhiHist->SetMarkerColor(1);
    TTbarPU200PhiHist->SetLineColor(1);
    TTbarPU0PhiHist->SetMarkerColor(2);
    TTbarPU0PhiHist->SetLineColor(2);
    muonPhiHist->SetMarkerStyle(kFullSquare);
    TTbarPU200PhiHist->SetMarkerStyle(kFullTriangleUp);
    TTbarPU0PhiHist->SetMarkerStyle(kFullDiamond);

    muonZ0Hist->SetMarkerColor(8);
    muonZ0Hist->SetLineColor(8);
    electronZ0Hist->SetMarkerColor(4);
    electronZ0Hist->SetLineColor(4);
    TTbarPU200Z0Hist->SetMarkerColor(1);
    TTbarPU200Z0Hist->SetLineColor(1);
    TTbarPU0Z0Hist->SetMarkerColor(2);
    TTbarPU0Z0Hist->SetLineColor(2);
    muonZ0Hist->SetMarkerStyle(kFullSquare);
    TTbarPU200Z0Hist->SetMarkerStyle(kFullTriangleUp);
    TTbarPU0Z0Hist->SetMarkerStyle(kFullDiamond);

    // Draw and Print Histograms to pdf
    TString saveDir = "plotPDFs/";

    TCanvas c;
    double max;

    // Normalizing histos
    vector<TH1*> allHists{muonEtaHist, muonPtHist, muonPhiHist, muonZ0Hist, 
    electronEtaHist, electronPtHist, electronPhiHist, electronZ0Hist, 
    TTbarPU200EtaHist, TTbarPU200PtHist, TTbarPU200PhiHist, TTbarPU200Z0Hist, 
    TTbarPU0EtaHist, TTbarPU0PtHist, TTbarPU0PhiHist, TTbarPU0Z0Hist};

    int histNum = allHists.size();
    for (int i = 0; i < histNum; i++) {
      allHists[i]->Scale(1./allHists[i]->Integral(), "width");
    }

    
    // eta
    TLegend* etaLeg = new TLegend(0.43, 0.18, 0.63, 0.37);
    etaLeg->AddEntry(muonEtaHist, "Single Muon PU=0");
    etaLeg->AddEntry(TTbarPU200EtaHist, "TTbar PU=200");
    etaLeg->AddEntry(TTbarPU0EtaHist, "TTbar PU=0");
    etaLeg->AddEntry(electronEtaHist, "Single Electron PU=0");

    muonEtaHist->Draw();
    electronEtaHist->Draw("same");
    TTbarPU200EtaHist->Draw("same");
    TTbarPU0EtaHist->Draw("same");
    etaLeg->Draw();
    vector<TH1*> etaHists{muonEtaHist, electronEtaHist, TTbarPU200EtaHist, TTbarPU0EtaHist}; // 
    max = GetMaxHists(etaHists);
    muonEtaHist->GetYaxis()->SetRangeUser(0, max*1.1);
    muonEtaHist->GetYaxis()->SetTitle("Tracks Probability Distribution / 0.026");
    c.SaveAs(saveDir + prop + "eta_MuonElectronTTbar.pdf");

    // pT
    TLegend* ptLeg = new TLegend(0.4, 0.45, 0.69, 0.65);
    ptLeg->AddEntry(muonEtaHist, "Single Muon PU=0");
    ptLeg->AddEntry(TTbarPU200EtaHist, "TTbar PU=200");
    ptLeg->AddEntry(TTbarPU0EtaHist, "TTbar PU=0");
    ptLeg->AddEntry(electronEtaHist, "Single Electron PU=0");

    muonPtHist->Draw();
    electronPtHist->Draw("same");
    TTbarPU200PtHist->Draw("same");
    TTbarPU0PtHist->Draw("same");
    ptLeg->Draw();
    vector<TH1*> ptHists{muonPtHist, electronPtHist, TTbarPU200PtHist, TTbarPU0PtHist}; // 
    max = GetMaxHists(ptHists);
    muonPtHist->GetYaxis()->SetRangeUser(0, max*1.1);
    muonPtHist->GetYaxis()->SetTitle("Tracks Probability Distribution / 0.5 GeV");
    c.SaveAs(saveDir + prop + "pt_MuonElectronTTbar.pdf");

    // zoom in on pT 
    muonPtHist->GetYaxis()->SetRangeUser(0, 0.1);
    c.SaveAs(saveDir + prop + "pt_MuonElectronTTbar_zoom.pdf");

    // phi
    TLegend* phiLeg = new TLegend(0.25, 0.3, 0.55, 0.5);
    phiLeg->AddEntry(muonEtaHist, "Single Muon PU=0");
    phiLeg->AddEntry(TTbarPU200EtaHist, "TTbar PU=200");
    phiLeg->AddEntry(TTbarPU0EtaHist, "TTbar PU=0");
    phiLeg->AddEntry(electronEtaHist, "Single Electron PU=0");

    muonPhiHist->Draw();
    electronPhiHist->Draw("same");
    TTbarPU200PhiHist->Draw("same");
    TTbarPU0PhiHist->Draw("same");
    phiLeg->Draw();
    vector<TH1*> phiHists{muonPhiHist, electronPhiHist, TTbarPU200PhiHist, TTbarPU0PhiHist}; // 
    max = GetMaxHists(phiHists);
    muonPhiHist->GetYaxis()->SetRangeUser(0, max*1.1);
    muonPhiHist->GetYaxis()->SetTitle("Tracks Probability Distribution / 0.0314");
    c.SaveAs(saveDir + prop + "phi_MuonElectronTTbar.pdf");

    // z0
    TLegend* z0Leg = new TLegend(0.20, 0.65, 0.50, 0.85);
    z0Leg->AddEntry(muonEtaHist, "Single Muon PU=0");
    z0Leg->AddEntry(TTbarPU200EtaHist, "TTbar PU=200");
    z0Leg->AddEntry(TTbarPU0EtaHist, "TTbar PU=0");
    z0Leg->AddEntry(electronEtaHist, "Single Electron PU=0");

    muonZ0Hist->Draw();
    electronZ0Hist->Draw("same");
    TTbarPU200Z0Hist->Draw("same");
    TTbarPU0Z0Hist->Draw("same");
    z0Leg->Draw();
    vector<TH1*> z0Hists{muonZ0Hist, electronZ0Hist, TTbarPU200Z0Hist, TTbarPU0Z0Hist}; 
    max = GetMaxHists(z0Hists);
    muonZ0Hist->GetYaxis()->SetRangeUser(0, max*1.1);
    muonZ0Hist->GetYaxis()->SetTitle("Tracks Probability Distribution / 0.35");
    c.SaveAs(saveDir + prop + "z0_MuonElectronTTbar.pdf");

    // ---------------Debug code-----------------
    
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