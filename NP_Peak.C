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

void NP_Peak(){
//void hcal_analysis(char const *input_rootfile){ //alternate to tchain
    gROOT->Reset();
    using namespace std;

//####### tchain adding #####
    TChain *tchnT = new TChain("T");
      tchnT->Add(Form("/w/halla-scshelf2102/sbs/saru/Simu_GEnRP/RootFile_Simu/genrp_Saru_New.root"));


    cout<<"processing input files "<<endl;
    if(!tchnT || tchnT->IsZombie()){cout <<"no file found"<<endl; exit(-1);}

    cout<<"##### Following files are added to TChain for this analysis #####"<<endl;
    TObjArray *files = tchnT->GetListOfFiles();
    TIter next(files);
    TChainElement *chEl=0;
    while(( chEl=(TChainElement*)next()  )){
        cout<< chEl->GetTitle()  <<endl;
    }

    TTreeReader tree(tchnT);
    TFile *fout = new TFile("./SimuSaru.root","RECREATE");

TTreeReaderArray<int>OTrack_PID(tree,"OTrack.PID");
TTreeReaderArray<int>OTrack_MID(tree,"OTrack.MID");
TTreeReaderArray<int>OTrack_TID(tree,"OTrack.TID");
TTreeReaderArray<double>OTrack_posx(tree,"OTrack.posx");
TTreeReaderArray<double>OTrack_posy(tree,"OTrack.posy");
TTreeReaderArray<double>OTrack_posz(tree,"OTrack.posz");
TTreeReaderArray<double>OTrack_momx(tree,"OTrack.momx");
TTreeReaderArray<double>OTrack_momy(tree,"OTrack.momy");
TTreeReaderArray<double>OTrack_momz(tree,"OTrack.momz");
TTreeReaderArray<double>OTrack_Etot(tree,"OTrack.Etot");

TTreeReaderValue<int> ActAn_nhits(tree, "Harm.ActAnScint.hit.nhits");
TTreeReaderArray<double> ActAn_t(tree, "Harm.ActAnScint.hit.tmin");
TTreeReaderArray<int> ActAn_cell(tree, "Harm.ActAnScint.hit.cell");
TTreeReaderArray<int> ActAn_sdtridx(tree, "Harm.ActAnScint.hit.sdtridx");
TTreeReaderArray<int> ActAn_ptridx(tree, "Harm.ActAnScint.hit.ptridx");
TTreeReaderArray<int> ActAn_otridx(tree, "Harm.ActAnScint.hit.otridx");
TTreeReaderArray<double> ActAn_esum(tree, "Harm.ActAnScint.hit.sumedep");
TTreeReaderArray<double> ActAn_xhitg(tree,"Harm.ActAnScint.hit.xhitg");
TTreeReaderArray<double> ActAn_yhitg(tree,"Harm.ActAnScint.hit.yhitg");
TTreeReaderArray<double> ActAn_zhitg(tree,"Harm.ActAnScint.hit.zhitg");

TTreeReaderValue<int> HCal_nhits(tree, "Harm.HCalScint.hit.nhits");
TTreeReaderArray<double> HCal_t(tree, "Harm.HCalScint.hit.tmin");
TTreeReaderArray<int> HCal_cell(tree, "Harm.HCalScint.hit.cell");
TTreeReaderArray<int> HCal_sdtridx(tree, "Harm.HCalScint.hit.sdtridx");
TTreeReaderArray<int> HCal_ptridx(tree, "Harm.HCalScint.hit.ptridx");
TTreeReaderArray<int> HCal_otridx(tree, "Harm.HCalScint.hit.otridx");

TTreeReaderArray<int> ActAn_pid(tree, "Harm.ActAnScint.pid");
TTreeReaderArray<int> ActAn_mid(tree, "Harm.ActAnScint.mid");
TTreeReaderArray<int> HCal_pid(tree, "Harm.HCalScint.pid");
TTreeReaderArray<int> HCal_mid(tree, "Harm.HCalScint.mid");

TTreeReaderValue<int> ev_nucl(tree,"ev.nucl");

// ===== Kinematics =====
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(0001);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleSize(0.045,"TXYZ");
    gStyle->SetTitleSize(0.045,"T");
    gStyle->SetTitleXOffset(1.1);
    gStyle->SetTitleYOffset(1.3);
    gStyle->SetTitleFont(62,"TXYZ");
    gStyle->SetLabelSize(0.045,"XYZ");
    gStyle->SetStatY(0.90);
    gStyle->SetStatX(0.90);
    gStyle->SetStatW(0.22);
    gStyle->SetStatH(0.23);
    gStyle->SetPadLeftMargin(0.05); //best setting previously for (3,2)
    gStyle->SetPadRightMargin(0.05); ////best setting previously for (3,2)
    gStyle->SetPadBottomMargin(0.10);
    gStyle->SetStatStyle(0); //set statbox transparent
    gStyle->SetTitleStyle(0);
    TGaxis::SetMaxDigits(4);
    gROOT->ForceStyle(true);


TH1D *h_AA_Etot_neutron = new TH1D("h_AA_Etot_neutron","AA E_{tot};E_{tot} (GeV);Events",200, 2, 4);
TH1D *h_AA_Etot_proton = new TH1D("h_AA_Etot_proton", "AA E_{tot};E_{tot} (GeV);Events", 200, 2, 4);
TH1D *h_HCal_Etot_neutron = new TH1D("h_HCal_Etot_neutron","HCal E_{tot};E_{tot} (GeV);Events",200, -1, 4);
TH1D *h_HCal_Etot_proton = new TH1D("h_HCal_Etot_proton", "HCal E_{tot};E_{tot} (GeV);Events", 200, -1, 4);




int neutron_events_AA =0;
int proton_events_AA  = 0;
int neutron_events_HCal =0;
int proton_events_HCal  = 0;


while (tree.Next()) {
    bool neutron_event_AA = false;
    bool proton_event_AA  = false;
    for (int jj = 0; jj < *ActAn_nhits; jj++) {
        int idxx = ActAn_otridx[jj];
        double pid_AA  = OTrack_PID[idxx];
        double tid_AA  = OTrack_TID[idxx];
        double mid_AA  = OTrack_MID[idxx];
        double AA_Etot = OTrack_Etot[idxx];

       if (tid_AA == 2 && mid_AA == 0) {
            if (pid_AA == 2112) {          // neutron
                h_AA_Etot_neutron->Fill(AA_Etot);
                neutron_event_AA = true;
            }

            else if (pid_AA == 2212) {     // proton
                h_AA_Etot_proton->Fill(AA_Etot);
                proton_event_AA = true;
       }  
    }}
    if (neutron_event_AA)
        neutron_events_AA++;
    if (proton_event_AA)
        proton_events_AA++;

    ////////////////////////////////////////////////////////////////
    bool neutron_event_HCal = false;
    bool proton_event_HCal = false;
    for (int j = 0; j < *HCal_nhits; j++) {
        int idx = HCal_otridx[j];
        double pid_HCal  = OTrack_PID[idx];
        double tid_HCal  = OTrack_TID[idx];
        double mid_HCal  = OTrack_MID[idx];
        double HCal_Etot = OTrack_Etot[idx];

       if (tid_HCal != 2 && mid_HCal != 0) {
            if (pid_HCal == 2112) {          // neutron
                h_HCal_Etot_neutron->Fill(HCal_Etot);
                neutron_event_HCal = true;
            }

            else if (pid_HCal == 2212) {     // proton
                h_HCal_Etot_proton->Fill(HCal_Etot);
                proton_event_HCal = true;
            }

    }}
    if (neutron_event_HCal)
        neutron_events_HCal++;
    if (proton_event_HCal)
        proton_events_HCal++;



}

// Create one canvas with two pads
TCanvas *c = new TCanvas("c", "AA and HCal Energy", 1400, 800);
c->Divide(2, 1);
// Top: AA
c->cd(1);
// Histogram appearance
h_AA_Etot_neutron->SetLineColor(kBlue+1);
h_AA_Etot_neutron->SetLineWidth(2);
h_AA_Etot_proton->SetLineColor(kRed+1);
h_AA_Etot_proton->SetLineWidth(2);
// Set y-axis maximum
h_AA_Etot_neutron->SetMaximum(1000);
// Draw all three on the same plot
h_AA_Etot_neutron->Draw("HIST");
h_AA_Etot_proton->Draw("HIST SAME");
// Legend
TLegend *leg = new TLegend(0.60,0.65,0.88,0.88);
leg->AddEntry(h_AA_Etot_neutron, "Neutron (PID=2112)","l");
leg->AddEntry(h_AA_Etot_proton,"Proton (PID=2212)","l");
leg->Draw();
c->cd(2);
h_HCal_Etot_neutron->SetLineColor(kBlue+1);
h_HCal_Etot_neutron->SetLineWidth(2);
h_HCal_Etot_proton->SetLineColor(kRed+1);
h_HCal_Etot_proton->SetLineWidth(2);
// Set y-axis maximum
//h_HCal_Etot_neutron->SetMaximum(12000);
// Draw all three on the same plot
h_HCal_Etot_neutron->Draw("HIST");
h_HCal_Etot_proton->Draw("HIST SAME");
TLegend *leg1 = new TLegend(0.60,0.65,0.88,0.88);
leg1->AddEntry(h_HCal_Etot_neutron, "Neutron (PID=2112)","l");
leg1->AddEntry(h_HCal_Etot_proton,"Proton (PID=2212)","l");
leg1->Draw();

}
