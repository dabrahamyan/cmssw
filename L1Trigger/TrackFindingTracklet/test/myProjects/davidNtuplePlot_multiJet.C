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
#include "plotstyle.h"
#include "TROOT.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

// FAILED                                  **************************************************************
// **********************************************************************************

void davidNtuplePlot_multiJet () {
    //gROOT->ProcessLine(".L davidNtuplePlot.C");
    gROOT->LoadMacro("davidNtuplePlot.C");
    davidNtuplePlot("TTbar_PU200_D88_HYBRID_noMETrunc", "/eos/user/d/dabraham/L1NtupleTrackExamples/", "", 0);
}