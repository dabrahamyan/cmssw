////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Info: 
// tp        - all tracking particles
// matchtrk  - *L1 track* properties, for tracking particles matched to an L1 track
// trk       - all L1 tracks
// tp_nmatch - number of matched tracks for a tracking particle
// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


void ZRplot_badEta_finder() {
        // File and directory to get Ntuple from -------------------------------------------------------------------------
        TString file = "newkfDebug_SingleMuon_DR_off"; // without ".root"
        TString fileDir = "../LocalChecks/";

        // read ntuples ---------------------------------------------------------------------------------------------------
        TChain* tree = new TChain("L1TrackNtuple/eventTree");
        tree->Add(fileDir + file + ".root");

        if (tree->GetEntries() == 0) {
        cout << "File doesn't exist or is empty, returning..."
        << endl;
        return; 
        }

        // define leaves and branches ----------------------------------------------------------------------------------------------------------------
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


        float etaResidual;
        int minResidualInd;
        float tempResidual;
        float residual;
        float absEta;

        int nevt = tree->GetEntries();
        // event loop -------------------------------------------------------------------------------------------
        int entries[] = {4006, 4045, 4056, 2047, 2012};
        for (int i = 0; i < nevt; i++) { 
                tree->GetEntry(i, 0);
                
                // cout << "----------------------------------------------------------------" << endl;
                // cout << "Event " << i << ":" << endl;

                // tracking particle loop -----------------------------------------------------------------------
                for (int it = 0; it < (int)tp_eta->size(); it++) {
                        // if 1.3 < eta 1.4 (bad bin)
                        if (std::abs(tp_eta->at(it)) < 1.9 || std::abs(tp_eta->at(it)) > 2.0) { //1.3-1.4  1.9-2.0
                                continue;
                        }
                        // Choosing best trk based on eta and assigning it to matchtrk_eta
                        etaResidual = 100; // set high because we are finding minimum
                        minResidualInd = 100; // put high so that error shows up if not reassigned
                        tempResidual;

                        for (int itrk = 0; itrk < (int)trk_eta->size(); itrk++) {
                                tempResidual = std::abs(tp_eta->at(it) - trk_eta->at(itrk));

                                if (tempResidual < etaResidual) {
                                etaResidual = tempResidual;
                                minResidualInd = itrk;
                                }
                        }

                        // effectively, if minResidualInd was assigned to anything at all
                        if (minResidualInd < 99) {
                        matchtrk_eta->at(it) = trk_eta->at(minResidualInd);
                        //matchtrk_pt->at(it) = trk_pt->at(minResidualInd);
                        matchtrk_phi->at(it) = trk_phi->at(minResidualInd);
                        matchtrk_z0->at(it) = trk_z0->at(minResidualInd);
                        };

                        // if matchtrk - tp eta residual is close to 0.1, print 
                        residual = std::abs(tp_eta->at(it) - matchtrk_eta->at(it));
                        if (residual > 0.07 && residual < 0.1) { // residual > 0.035 && residual < 0.1
                                cout << "Event: " << i << "\ttp: " << it << "\teta residual: " << residual << endl;
                        }
                }
        }
}