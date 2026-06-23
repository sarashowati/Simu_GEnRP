#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TF1.h>
#include <TH2F.h>
#include <TFile.h>
#include <iomanip>
#include <time.h>
#include <TVector3.h>
#include <TMath.h>
#include "TCutG.h"
#include <TLine.h>
#include <TLorentzVector.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <TGraphErrors.h>
#include "TLegend.h"
#include "TVirtualPad.h"
#include "TObjString.h"
#include "TVirtualHistPainter.h"
#include "THLimitsFinder.h"
#include "TFitResult.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>


void Multiplicity_AA() {
    TFile *f = TFile::Open("/w/halla-scshelf2102/sbs/saru/Simu_GEnRP/genrp_agc_sim4.root");
    TTree *tree = (TTree*)f->Get("T");

    TTreeReader reader(tree);
    TTreeReaderValue<int> ActAna_hit_nhits(reader, "Harm.ActAnScint.hit.nhits");

  TH1D *h_AA_mult = new TH1D("h_AA_mult","Active Analyzer Multiplicity;AA hit multiplicity;Counts",   32, 0.5, 32.5);

    while(reader.Next()) {
        h_AA_mult->Fill(*ActAna_hit_nhits);
    }
    TCanvas *c1 = new TCanvas("c1","AA Multiplicity",800,600);

    h_AA_mult->Draw();
//    gPad->SetLogy(); // optional
}
