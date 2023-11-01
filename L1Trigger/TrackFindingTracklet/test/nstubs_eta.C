///////////////////////////////////////////////////////////////////////////////
// This code will make a graph of nstubs vs eta
// It reads in a root file with an Ntuple and then makes a plot :D
//
// By David Abrahamyan Oct 11 2023

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include <TError.h>
#include "TSystem.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;


void nstubs_eta () {
    // Read in root file and ntuple
    TChain* tree = new TChain("L1TrackNtuple/eventTree");
    tree->Add("combinedDebug_SingleMuon_set1_HYBRID_jan27commit.root"); // root file here

    if (tree->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;  //cout's kept in this file as it is an example standalone plotting script, not running in central CMSSW
    return;
    }

        vector<int>* trk_nstub;
        vector<double>* trk_eta;

        TBranch* b_trk_nstub;
        TBranch* b_trk_eta;

        trk_nstub = 0;
        trk_eta = 0;

        tree->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
        tree->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);

    int nevt = tree->GetEntries();
    cout << "number of events = " << nevt << endl;

    // ----------------------------------------------------------------------------------------------------------------
    // event loop
    for (int i = 0; i < nevt; i++) {
        tree->GetEntry(i, 0);
        
        for (int i=0; i<(int)trk_eta->size();i++){
            cout << trk_eta->at(i) << endl;     
        }

        double etaArr[100], nstubArr[100]; // create empty arrays
        double binStartMid = -2.35; // middle of bin all th eway to left of graph
        double binWidth = 0.1;

        // create eta array
        int n; // # of bins (will go into TGraph)
        n = 0;
        for (double b = binStartMid; b < 2.4; b = b + binWidth) {
            etaArr[n] = b;
            n++;
        }
    }

    //for ()

    // create nstubs array


    
    // create graph
    // auto g = new TGraph(n, etaArr, nstubArr);
    // g->SetTitle("; #eta; nstubs");
    // g->Draw("AC*");
    
    //delete tree, trk_nstub, trk_eta, b_trk_eta, b_trk_nstub;



}
