// ----------------------------------------------------------------------------------------------------
// Makes a double bar graph with custom x axis labels
// Manually input data 
// referenced here: https://root.cern/doc/master/classTHistPainter.html#HP100 under "the bar chart option"
//
// By David Abrahamyan September 11, 2023
// -----------------------------------------------------------------------------------------------------

#include <fstream>

void mySmallText(Double_t x, Double_t y, Color_t color, char* text);

void efficiencyComparisonBarChart () {

   const Int_t nx = 9;
   string os_X[nx]   = {"IR","VMR","TE","TC","PR","ME","MC","TB", "DR"};
   // All tp
   float hybrid[nx] = {96.1977, 96.1178, 96.1791, 96.1605, 95.7747, 95.7747, 95.6393, 96.1948, 96.1977};
   float newkf[nx] = {95.8282, 95.8593, 95.4922, 95.4816, 95.4922, 95.4922, 95.541, 95.877, 95.877};
 
   // // In Jet pT > 30 GeV
   // float hybrid[nx] = {96.2799, 96.1962, 96.2554, 96.2288, 95.7401, 95.3698, 95.5297, 96.2765, 96.2799};
   // float newkf[nx] = {95.8983, 95.8442, 95.8727, 95.3976, 95.3907, 95.3976, 95.4536, 95.8983, 95.8983};

   // In Jet pT > 100 GeV
   // float hybrid[nx] = {96.1617, 96.0559, 96.1003, 96.0305, 94.9048, 93.8278, 94.1092, 96.1511, 96.1617};
   // float newkf[nx] = {95.4892, 95.4028, 95.4306, 94.2334, 94.3476, 94.2334, 94.2396, 95.4892, 95.4892};

   // // In Jet pT > 200 GeV
   // float hybrid[nx] = {95.5925, 95.3802, 95.4676, 95.1679, 91.8716, 88.875, 88.8875, 95.58, 95.5925};
   // float newkf[nx] = {94.9319, 94.8245, 94.8603, 90.8309, 92.0487, 90.8309, 90.7951, 94.9319, 94.9319};

   auto cb = new TCanvas("cb","cb",600,400);
   cb->SetGrid();
 
   gStyle->SetHistMinimumZero();
 
   auto h1b = new TH1F("h1b","",nx,0,nx);
   h1b->SetFillColor(4);
   h1b->SetBarWidth(0.4);
   h1b->SetBarOffset(0.1);
   h1b->SetStats(0);
   h1b->SetMinimum(83);
   h1b->SetMaximum(100);

   //h1b->SetTitle(";Truncated Module;Efficiency");
   h1b->GetXaxis()->SetTitle("Truncated Module");
   h1b->GetXaxis()->SetTitleSize(.045);
   h1b->GetXaxis()->CenterTitle(true);

   h1b->GetYaxis()->SetTitle("Efficiency (%)");
   h1b->GetYaxis()->SetTitleSize(.045);
   h1b->GetYaxis()->CenterTitle(true);

   h1b->GetXaxis()->SetLabelSize(0.06);
   h1b->GetYaxis()->SetLabelSize(0.04);
   //h1b->GetXaxis()->SetTitleSize(0.2);
   //h1b->GetXaxis()->SetTitleSize(0.2);
 
   int i;
   for (i=1; i<=nx; i++) {
      h1b->SetBinContent(i, hybrid[i-1]);
      h1b->GetXaxis()->SetBinLabel(i,os_X[i-1].c_str());
   }
 
   //cb->Update();

   

   
   h1b->Draw("b");
 
   auto h2b = new TH1F("h2b","h2b",nx,0,nx);
   h2b->SetFillColor(38);
   h2b->SetBarWidth(0.4);
   h2b->SetBarOffset(0.5);
   h2b->SetStats(0);
   for (i=1;i<=nx;i++) h2b->SetBinContent(i, newkf[i-1]);
   
   TLegend* leg = new TLegend(0.65, 0.75, 0.87, 0.87);
   leg->AddEntry(h1b, "Hybrid");
   leg->AddEntry(h2b, "Hybrid New KF");

   h2b->Draw("b same");
   leg->Draw();

   cb->Update();

   // line values for full and no trunc
   // All Particles
   double hybridFull = 94.8894;
   double hybridNo = 96.1977;
   double newkfFull = 94.9644;
   double newkfNo = 95.877;

   // In Jet > 30 GeV
   // double hybridFull = 94.572;
   // double hybridNo = 96.2799;
   // double newkfFull = 94.712;
   // double newkfNo = 95.8983;

   // In Jet > 100 GeV
   // double hybridFull = 91.9827;
   // double hybridNo = 96.1617;
   // double newkfFull = 92.5735;
   // double newkfNo = 95.4892;

   // In Jet > 200 GeV
   // double hybridFull = 84.1304;
   // double hybridNo = 95.5925;
   // double newkfFull = 86.712;
   // double newkfNo = 94.9319;

   TLine *hybrid_full = new TLine(h1b->GetXaxis()->GetXmin(), hybridFull, h1b->GetXaxis()->GetXmax(), hybridFull);
   TLine *hybrid_no = new TLine(h1b->GetXaxis()->GetXmin(), hybridNo, h1b->GetXaxis()->GetXmax(), hybridNo);
   TLine *newkf_full = new TLine(h1b->GetXaxis()->GetXmin(), newkfFull, h1b->GetXaxis()->GetXmax(), newkfFull);
   TLine *newkf_no = new TLine(h1b->GetXaxis()->GetXmin(), newkfNo, h1b->GetXaxis()->GetXmax(), newkfNo);

   hybrid_full->SetLineColor(kOrange-3);
   hybrid_full->SetLineWidth(3);
   hybrid_full->Draw("same");

   hybrid_no->SetLineColor(kOrange-3);
   hybrid_no->SetLineWidth(3);
   hybrid_no->Draw("same");

   newkf_full->SetLineColor(kGreen+1);
   newkf_full->SetLineWidth(3);
   newkf_full->Draw("same");

   newkf_no->SetLineColor(kGreen+1);
   newkf_no->SetLineWidth(3);
   newkf_no->Draw("same");

   char ctxt1[500];
   sprintf(ctxt1, "Hybrid Full Trunc"); // Add label saying 
   mySmallText(0.13, 0.58, kOrange-3, ctxt1); // 

   char ctxt2[500];
   sprintf(ctxt2, "Hybrid No Trunc"); // Add label saying 
   mySmallText(0.13, 0.73, kOrange-3, ctxt2); // 

   char ctxt3[500];
   sprintf(ctxt3, "NewKF Full Trunc"); // Add label saying 
   mySmallText(0.4, 0.58, kGreen+1, ctxt3); // which data set it is

   char ctxt4[500];
   sprintf(ctxt4, "NewKF No Trunc"); // Add label saying 
   mySmallText(0.4, 0.73, kGreen+1, ctxt4); // 



   cb->Update();
}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.043;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}