////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Info: 
// tp        - all tracking particles
// matchtrk  - *L1 track* properties, for tracking particles matched to an L1 track
// trk       - all L1 tracks
// tp_nmatch - number of matched tracks for a tracking particle
// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


void showTpMatchtrkTrk() {
        // File and directory to get Ntuple from -------------------------------------------------------------------------
        TString file = "hybridvsnewkf_SingleMuonPU0D88_NEWKF_1_numEvent10000"; // without ".root"
        TString fileDir = "/eos/user/d/dabraham/L1NtupleTrackExamples/";

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

        int nevt = tree->GetEntries();
        // event loop -------------------------------------------------------------------------------------------
        int entries[] = {4006, 4045, 4056, 2047, 2012};
        for (int i : entries) { //int i = 4006; i < 4007; i++
                tree->GetEntry(i, 0);
                
                cout << "----------------------------------------------------------------" << endl;
                cout << "Event " << i << ":" << endl;
                //cout << "tp\t\tmatchtrk\ttrk" << endl;
                // tracking particle loop -----------------------------------------------------------------------
                cout << "tp_nmatch:" << endl;
                for (int it = 0; it < (int)tp_phi->size(); it++) {
                        cout << tp_nmatch->at(it) << endl;
                }
                cout << endl;

                cout << "tp_phi\tmatchtrk_phi\ttrk_phi" << endl;
                for (int it = 0; it < (int)tp_phi->size(); it++) {
                        cout << tp_phi->at(it) << "\t" << matchtrk_phi->at(it) << "\t" << trk_phi->at(it) << endl;
                }
                cout << endl;

                cout << "trk_z0\ttrk_nstub" << endl;
                for (int it = 0; it < (int)trk_phi->size(); it++) {
                        //if (trk_nstub->at(it) < 4)
                        //continue;
                        cout << trk_phi->at(it) << "\t" << trk_nstub->at(it) << endl;
                }
                cout << endl;

                cout << "tp_eta" << endl;
                for (int it = 0; it < (int)tp_eta->size(); it++) {
                        cout << tp_eta->at(it) << endl;
                }
        }
}