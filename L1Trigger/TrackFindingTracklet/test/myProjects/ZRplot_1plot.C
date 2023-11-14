////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Info: 
// tp        - all tracking particles
// matchtrk  - *L1 track* properties, for tracking particles matched to an L1 track
// trk       - all L1 tracks
// tp_nmatch - number of matched tracks for a tracking particle
// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void mySmallText(Double_t x, Double_t y, Color_t color, char* text);

void ZRplot_1plot() {
        // File and directory to get Ntuple from -------------------------------------------------------------------------
        TString file = "newkfDebug_SingleMuon_DR_off_stubs"; // without ".root"
        TString fileDir = "../LocalChecks/";

        // save settings
        TString saveDir = "hybrid_vs_newkf_plots/";
        TString saveName = "ZRplot.pdf";

        // read ntuples ---------------------------------------------------------------------------------------------------
        TChain* tree = new TChain("L1TrackNtuple/eventTree");
        tree->Add(fileDir + file + ".root");

        if (tree->GetEntries() == 0) {
        cout << "File doesn't exist or is empty, returning..."
        << endl;
        return; 
        }

        // -----------------------------------------------------------------------------------------------
        // define leaves and branches 

        vector<float>* tp_z0;
        vector<float>* matchtrk_z0;
        vector<float>* trk_z0;
        vector<float>* tp_eta;
        vector<int>* tp_nmatch;
        vector<int>* matchtrk_nstub;
        vector<int>* trk_nstub;
        vector<float>* tp_phi;
        vector<float>* matchtrk_phi;
        vector<float>* trk_phi;
        vector<float>* trk_eta;
        vector<float>* matchtrk_eta;
        vector<float>* allstub_z;
        vector<float>* allstub_x;
        vector<float>* allstub_y;
        
        TBranch* b_tp_z0;
        TBranch* b_matchtrk_z0;
        TBranch* b_trk_z0;
        TBranch* b_tp_eta;
        TBranch* b_tp_nmatch;
        TBranch* b_matchtrk_nstub;
        TBranch* b_trk_nstub;
        TBranch* b_tp_phi;
        TBranch* b_matchtrk_phi;
        TBranch* b_trk_phi;
        TBranch* b_trk_eta;
        TBranch* b_matchtrk_eta;
        TBranch* b_allstub_z;
        TBranch* b_allstub_x;
        TBranch* b_allstub_y;

        tp_z0 = 0;
        matchtrk_z0 = 0;
        trk_z0 = 0;
        tp_phi = 0;
        matchtrk_phi = 0;
        trk_phi = 0;
        tp_eta = 0;
        tp_nmatch = 0;
        matchtrk_nstub = 0;
        trk_nstub = 0;
        trk_eta = 0;
        matchtrk_eta = 0;
        allstub_z = 0;
        allstub_x = 0;
        allstub_y = 0;

        tree->SetBranchAddress("tp_z0", &tp_z0, &b_tp_z0);
        tree->SetBranchAddress("matchtrk_z0", &matchtrk_z0, &b_matchtrk_z0);
        tree->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
        tree->SetBranchAddress("tp_phi", &tp_phi, &b_tp_phi);
        tree->SetBranchAddress("matchtrk_phi", &matchtrk_phi, &b_matchtrk_phi);
        tree->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
        tree->SetBranchAddress("tp_eta", &tp_eta, &b_tp_eta);
        tree->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
        tree->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);
        tree->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
        tree->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
        tree->SetBranchAddress("matchtrk_eta", &matchtrk_eta, &b_matchtrk_eta);
        tree->SetBranchAddress("allstub_z", &allstub_z, &b_allstub_z);
        tree->SetBranchAddress("allstub_x", &allstub_x, &b_allstub_x);
        tree->SetBranchAddress("allstub_y", &allstub_y, &b_allstub_y);

        TCanvas c;
        char ctxt[500];
        char ctxt2[500];
        char ctxt3[500];

        
        // choose event and tp that has the bad eta residual
        int eventNum = 38;
        int tpNum = 0;

        tree->GetEntry(eventNum,0);

        // -----------------------------------------------------------------------------------------------
        // My own Duplicate Removal

        // Choosing best trk based on eta and assigning it to matchtrk_eta
        float etaResidual = 100; // set high because we are finding minimum
        int minResidualInd = 100; // put high so that error shows up if not reassigned
        float tempResidual;
        float residual;
        for (int itrk = 0; itrk < (int)trk_eta->size(); itrk++) {
                tempResidual = std::abs(tp_eta->at(tpNum) - trk_eta->at(itrk));

                if (tempResidual < etaResidual) {
                etaResidual = tempResidual;
                minResidualInd = itrk;
                }
        }

        // effectively, if minResidualInd was assigned to anything at all
        if (minResidualInd < 99) {
                matchtrk_eta->at(tpNum) = trk_eta->at(minResidualInd);
        };

        // --------------------------------------------------------------------------------------------
        // getting stub and line data

        // is tp_eta positive or negative (we can then know which stubs to use)
        bool zpos = false;
        bool zneg = false;
        bool zunclear = false;

        if (tp_eta->at(tpNum) > 0.1)
                zpos = true;
        else if (tp_eta->at(tpNum) < -0.1)
                zneg = true;
        else
                zunclear = true; // if its less than 0.1, some of the stubs might 
                                 // be back and forth between positive and negative z for one tp

        // get stub indicies
        vector<int> stubInd;
        if (zpos) {
                for (int istub = 0; istub < (int)allstub_z->size(); istub++) {
                        if (allstub_z->at(istub) > 0) {
                                stubInd.push_back(istub);
                        }
                }
        }
        else if (zneg) {
                for (int istub = 0; istub < (int)allstub_z->size(); istub++) {
                        if (allstub_z->at(istub) < 0) {
                                stubInd.push_back(istub);
                        }
                }
        }
        else if (zunclear) {
                throw std::invalid_argument("too close to eta=0 so it's unclear which stubs to use");
        }

        // get stub coordinates
        int nstub = stubInd.size();
        float z_stub[nstub], r_stub[nstub];
        for (int i = 0; i < nstub; i++) {
                z_stub[i] = allstub_z->at(stubInd[i]);
                r_stub[i] = sqrt(pow(allstub_x->at(stubInd[i]),2) + pow(allstub_y->at(stubInd[i]),2));
                // Debug stuff
                // cout << "z: " << z_stub[i] << endl;
                // cout << "r: " << r_stub[i] << endl;
        }

        // Get tp and matchmatchtrk lines
        float tp_x[2], tp_y[2], matchtrk_x[2], matchtrk_y[2];
        float matchtrk_theta = 2 * atan(exp(-1 * matchtrk_eta->at(tpNum))); // theta is in radians
        float tp_theta = 2 * atan(exp(-1 * tp_eta->at(tpNum)));
        float matchtrk_slope = cos(matchtrk_theta)/sin(matchtrk_theta);
        float tp_slope = cos(tp_theta)/sin(tp_theta);
        // first point
        tp_x[0] = 0;
        matchtrk_x[0] = 0;
        tp_y[0] = tp_z0->at(tpNum);
        matchtrk_y[0] = matchtrk_z0->at(tpNum);
        // second point
        tp_x[1] = 500;
        matchtrk_x[1] = 500;
        tp_y[1] = tp_y[0] + tp_x[1] * tp_slope;
        matchtrk_y[1] = matchtrk_y[0] + matchtrk_x[1] * matchtrk_slope;

        cout << "tp_slope: " << tp_slope << endl;
        cout << "matchtrk_slope: " << matchtrk_slope << endl;
        cout << "tp_z0: " << tp_z0->at(tpNum) << endl;
        cout << "matchtrk_z0: " << matchtrk_z0->at(tpNum) << endl;

        // ---------------------------------------------------------------------------------------------
        // Draw things 

        // Draw Stubs
        TGraph *stubGraph = new TGraph(nstub, r_stub, z_stub);
        stubGraph->SetTitle(";R;Z");
        stubGraph->Draw("A*");
        //stubGraph->SetMinimum(0); FAIL: that sets minimum of Y-axis. I want to set minimum of x-axis

        // Draw tp and matchtrk
        TGraph *tpGraph = new TGraph(2, tp_x, tp_y);
        TGraph *matchtrkGraph = new TGraph(2, matchtrk_x, matchtrk_y);
        tpGraph->SetLineColor(kGreen);
        matchtrkGraph->SetLineColor(kBlue);
        tpGraph->Draw("same");
        matchtrkGraph->Draw("same");

        // ----------------------------------------------------------------------------------------
        // dZ error bars


        // ----------------------------------------------------------------------------------------
        // legend and labels

        // make legend
        TLegend* leg = new TLegend(0.7, 0.4, 0.85, 0.55);
        leg->AddEntry(stubGraph, "stubs"); 
        leg->AddEntry(tpGraph, "tp");
        leg->AddEntry(matchtrkGraph, "matchtrk"); 
        leg->Draw();

        // make labels 
        TString eventLabel = "Event: " + to_string(eventNum);
        sprintf(ctxt, eventLabel);
        mySmallText(0.15, 0.52, 1, ctxt);
        TString tpLabel = "tp: " + to_string(tpNum);
        sprintf(ctxt2, tpLabel);
        mySmallText(0.15, 0.47, 1, ctxt2);
        TString residualLabel = "#eta Residual: " + to_string(etaResidual); 
        sprintf(ctxt3, residualLabel);
        mySmallText(0.15, 0.42, 1, ctxt3);

        c.SaveAs(saveDir + saveName);
}

void mySmallText(Double_t x, Double_t y, Color_t color, char* text) {
  Double_t tsize = 0.04;
  TLatex l;
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}