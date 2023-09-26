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
   // dIRECTORY TO SAVE PLOTS TO
   TString saveDir = "truncation_examination_assertsOn_oneTrunc/";

   // set this to true to do the combined plots
   bool combinedOn = "true";
   //Set this to 0 (all tp), 1 (injet > 30 GeV), 2 (injet > 100 GeV), 3 (injet > 200 GeV)
   int inJetSetting = 0;

   int nx;
   vector <string>binLabel;
   vector <double>hybrid;
   vector <double>newkf;

   double hybridFull;
   double hybridNo;
   double newkfFull;
   double newkfNo;

   int minPlotY;
   int maxPlotY;

   if (!combinedOn) {
      nx = 9;
      for (string i : {"IR","VMR","TE","TC","PR","ME","MC","TB","DR"}) {
         binLabel.push_back(i);
      }
      
      switch(inJetSetting) {
         case 0:
            for (double i : {96.1977, 96.1178, 96.1791, 96.1605, 95.7747, 95.7747, 95.6393, 96.1948, 96.1977}) {
               hybrid.push_back(i);
            }
            for (double i : {95.877, 95.8282, 95.8593, 95.4922, 95.4816, 95.4922, 95.541, 95.877, 95.877}) {
               newkf.push_back(i);
            }
            hybridFull = 94.8894;
            hybridNo = 96.1977;
            newkfFull = 94.9644;
            newkfNo = 95.877;
            cout << "true 0" << endl;
            break;

         case 1:
            for (double i : {96.2799, 96.1962, 96.2554, 96.2288, 95.7401, 95.3698, 95.5297, 96.2765, 96.2799}) {
               hybrid.push_back(i);
            }
            for (double i : {95.8983, 95.8442, 95.8727, 95.3976, 95.3907, 95.3976, 95.4536, 95.8983, 95.8983}) {
               newkf.push_back(i);
            }
            hybridFull = 94.572;
            hybridNo = 96.2799;
            newkfFull = 94.712;
            newkfNo = 95.8983;
            cout << "true 1" << endl;
            break;

         case 2:
            for (double i : {96.1617, 96.0559, 96.1003, 96.0305, 94.9048, 93.8278, 94.1092, 96.1511, 96.1617}) {
               hybrid.push_back(i);
            }
            for (double i : {95.4892, 95.4028, 95.4306, 94.2334, 94.3476, 94.2334, 94.2396, 95.4892, 95.4892}) {
               newkf.push_back(i);
            }
            hybridFull = 91.9827;
            hybridNo = 96.1617;
            newkfFull = 92.5735;
            newkfNo = 95.4892;
            break;

         case 3:
            for (double i : {95.5925, 95.3802, 95.4676, 95.1679, 91.8716, 88.875, 88.8875, 95.58, 95.5925}) {
               hybrid.push_back(i);
            }
            for (double i : {94.9319, 94.8245, 94.8603, 90.8309, 92.0487, 90.8309, 90.7951, 94.9319, 94.9319}) {
               newkf.push_back(i);
            }
            hybridFull = 84.1304;
            hybridNo = 95.5925;
            newkfFull = 86.712;
            newkfNo = 94.9319;
            break;
                  
         default:
            cout << "whaaaaa your number doesnt even make sense man" << endl << "you gotta pick a better number man" << endl << "idk smth like 0, 1, 2, or 3 maybe" << endl << "idk man get your stuff together" << endl;
            throw;
      }
      minPlotY = 83;
      maxPlotY = 100;
   }

   // combined off
   else {
      nx = 2;
      for (string i : {"TP", "MP"}) {
         binLabel.push_back(i);
      }

      switch (inJetSetting) {
         case 0:
            // All tp
            for (double i : {91.1722, 89.6291}) {
               hybrid.push_back(i);
            }
            for (double i : {93.1033, 92.0389}) {
               newkf.push_back(i);
            }
            hybridFull = 89.4804;
            hybridNo = 91.3847;
            newkfFull = 91.8994;
            newkfNo = 93.266;
            break;

         case 1:
            // In Jet pT > 30 GeV
            for (double i : {91.5496, 89.5644}) {
               hybrid.push_back(i);
            }
            for (double i : {93.4352, 92.0701}) {
               newkf.push_back(i);
            }
            hybridFull = 89.402;
            hybridNo = 91.8126;
            newkfFull = 91.9232;
            newkfNo = 93.6394;
            break;
            
         case 2:
            // In Jet pT > 100 GeV
            for (double i : {91.7952, 87.1632}) {
               hybrid.push_back(i);
            }
            for (double i : {93.2006, 90.0678}) {
               newkf.push_back(i);
            }
            hybridFull = 86.9909;
            hybridNo = 92.35;
            newkfFull = 89.8553;
            newkfNo = 93.601;
            break;

         case 3:
            // In Jet pT > 200 GeV
            for (double i : {90.4326, 78.2077}) {
               hybrid.push_back(i);
            }
            for (double i : {92.0918, 83.0038}) {
               newkf.push_back(i);
            }
            hybridFull = 78.2324;
            hybridNo = 91.7058;
            newkfFull = 83.1792;
            newkfNo = 93.2079;
            break;

         default:
            cout << "whaaaaa your number doesnt even make sense man" << endl << "you gotta pick a better number man" << endl << "idk smth like 0, 1, 2, or 3 maybe" << endl << "idk man get your stuff together" << endl;
            throw;
      }
      minPlotY = 75;
      maxPlotY = 100;
   }
   auto cb = new TCanvas("cb","cb",600,400);
   cb->SetGrid();
 
   gStyle->SetHistMinimumZero();
 
   auto h1b = new TH1F("h1b","",nx,0,nx);
   h1b->SetMinimum(minPlotY);
   h1b->SetMaximum(maxPlotY); 
   h1b->SetFillColor(4);
   h1b->SetBarWidth(0.4);
   h1b->SetBarOffset(0.1);
   h1b->SetStats(0);

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
      h1b->GetXaxis()->SetBinLabel(i,binLabel[i-1].c_str());
   }
   
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
   mySmallText(0.32, 0.515, kOrange-3, ctxt1); // 
   // y vals: 
   // combined: 0- 1-0.515 2-0.44 3-0.22

   char ctxt2[500];
   sprintf(ctxt2, "Hybrid No Trunc"); // Add label saying 
   mySmallText(0.32, 0.637, kOrange-3, ctxt2); // 
   // y vals: 
   // combined: 0- 1-0.595 2-0.615 3-0.57

   char ctxt3[500];
   sprintf(ctxt3, "NewKF Full Trunc"); // Add label saying 
   mySmallText(0.55, 0.60, kGreen+1, ctxt3); // which data set it is  
   // y vals: 
   // combined: 0- 1-0.595 2-0.53 3-0.37

   char ctxt4[500];
   sprintf(ctxt4, "NewKF No Trunc"); // Add label saying 
   mySmallText(0.55, 0.695, kGreen+1, ctxt4); // 
   // y vals: 
   // combined: 0- 1-0.705 2-0.705 3-0.69

   cb->Update();

   TString inJetString = to_string(inJetSetting);
   cb->SaveAs(saveDir + "effComparison_barChart_combined_" + inJetString + ".pdf");
}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.043;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}