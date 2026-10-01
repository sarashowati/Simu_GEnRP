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

void Neutron_AAHCAL_SDTrack(){
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


//HCal  and AA Hit_Timing Variables
TTreeReaderValue<int> HCal_hit_nhits(tree, "Harm.HCalScint.hit.nhits");
TTreeReaderArray<int> HCal_cell(tree, "Harm.HCalScint.hit.cell");
TTreeReaderArray<int> HCal_sdtridx(tree, "Harm.HCalScint.hit.sdtridx");
TTreeReaderArray<int> HCal_ptridx(tree, "Harm.HCalScint.hit.ptridx");
TTreeReaderArray<int> HCal_otridx(tree, "Harm.HCalScint.hit.otridx");
TTreeReaderArray<double> HCal_t(tree, "Harm.HCalScint.hit.tmin");
TTreeReaderArray<double> HCal_esum(tree, "Harm.HCalScint.hit.sumedep");
TTreeReaderArray<double> HCal_xhitg(tree,"Harm.HCalScint.hit.xhitg");
TTreeReaderArray<double> HCal_yhitg(tree,"Harm.HCalScint.hit.yhitg");
TTreeReaderArray<double> HCal_zhitg(tree,"Harm.HCalScint.hit.zhitg");

TTreeReaderArray<int>PTrack_PID(tree,"PTrack.PID");
TTreeReaderArray<int>PTrack_TID(tree,"PTrack.TID");
TTreeReaderArray<double>PTrack_posx(tree,"PTrack.posx");
TTreeReaderArray<double>PTrack_posy(tree,"PTrack.posy");
TTreeReaderArray<double>PTrack_posz(tree,"PTrack.posz");
TTreeReaderArray<double>PTrack_momx(tree,"PTrack.momx");
TTreeReaderArray<double>PTrack_momy(tree,"PTrack.momy");
TTreeReaderArray<double>PTrack_momz(tree,"PTrack.momz");
TTreeReaderArray<double>PTrack_Etot(tree,"PTrack.Etot");

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

TTreeReaderValue<int>SDTrack_ntracks(tree,"SDTrack.ntracks");
TTreeReaderArray<int>SDTrack_PID(tree,"SDTrack.PID");
TTreeReaderArray<int>SDTrack_MID(tree,"SDTrack.MID");
TTreeReaderArray<int>SDTrack_TID(tree,"SDTrack.TID");
TTreeReaderArray<double>SDTrack_posx(tree,"SDTrack.posx");
TTreeReaderArray<double>SDTrack_posy(tree,"SDTrack.posy");
TTreeReaderArray<double>SDTrack_posz(tree,"SDTrack.posz");
TTreeReaderArray<double>SDTrack_momx(tree,"SDTrack.momx");
TTreeReaderArray<double>SDTrack_momy(tree,"SDTrack.momy");
TTreeReaderArray<double>SDTrack_momz(tree,"SDTrack.momz");
TTreeReaderArray<double>SDTrack_Etot(tree,"SDTrack.Etot");
TTreeReaderArray<double>SDTrack_vx(tree,"SDTrack.vx");
TTreeReaderArray<double>SDTrack_vy(tree,"SDTrack.vy");
TTreeReaderArray<double>SDTrack_vz(tree,"SDTrack.vz");


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

TTreeReaderArray<int> ActAn_pid(tree, "Harm.ActAnScint.pid");
TTreeReaderArray<int> ActAn_mid(tree, "Harm.ActAnScint.mid");
TTreeReaderArray<int> HCal_pid(tree, "Harm.HCalScint.pid");
TTreeReaderArray<int> HCal_mid(tree, "Harm.HCalScint.mid");

TTreeReaderValue<int> ev_nucl(tree,"ev.nucl");

TTreeReaderValue<int> PRGEM_nhits(tree,"Harm.PRPolGEMFarSide.hit.nhits");
TTreeReaderArray<int> PRGEM_sdtridx(tree, "Harm.PRPolGEMFarSide.hit.sdtridx");


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

TH2D *h_2d_AA = new TH2D("h_2d_AA","AA xhitg vs AA zhitg;AA zhitg;HCal xhitg",100,4.6,5,100,-2.35,-2.05);
TH2D *h_2d_HCal = new TH2D("h_2d_HCal","HCal xhitg vs HCal zhitg;HCal zhitg;HCal xhitg",100,7.8,9.4,100,-5.1,-2.9);
TH1D *h_ev =new TH1D("h_ev","number of event generated",2,-0.5,1.5);
TH1D *h_ev_AA =new TH1D("h_ev_AA","number of event generated",2,-0.5,1.5);
TH1D *h_ev_HCalAA =new TH1D("h_ev_HCalAA","number of event generated",2,-0.5,1.5);

int Nentries = tree.GetEntries();
int neutron_AA_events =0;
int neutron_HCal_events=0;
int neutron_AA_HCal_events =0;

while (tree.Next()) {	
    bool condition_satisfied = false;
    for (int j = 0; j < *HCal_hit_nhits; j++) {
        int idx = HCal_sdtridx[j];
        double vx  = SDTrack_vx[idx];
        double vz  = SDTrack_vz[idx];
        double pid = SDTrack_PID[idx];
        double tid = SDTrack_TID[idx];
	double mid = SDTrack_MID[idx];
        double HCal_x = HCal_xhitg[idx];
        double HCal_z = HCal_zhitg[idx];
        if (pid == 2112 && tid == 2&&mid==0) {
	   h_2d_HCal->Fill(HCal_z, HCal_x);
            // Mark that this event has at least one qualifying hit
            condition_satisfied = true;
        }
    }
    // Fill ONLY ONCE per event
    if (condition_satisfied) {
	    neutron_AA_events++;
           h_ev->Fill(*ev_nucl);
    }

 bool condition_satisfied_AA = false;
for (int jj = 0; jj < *ActAn_nhits; jj++) {
	int idxx = ActAn_sdtridx[jj];
    double vx_AA = SDTrack_vx[idxx];
    double vz_AA = SDTrack_vz[idxx];
    double pid_AA =     SDTrack_PID[idxx];
     double tid_AA =     SDTrack_TID[idxx];
      double mid_AA =     SDTrack_MID[idxx];
     double AA_x = ActAn_xhitg[idxx];
        double AA_z = ActAn_zhitg[idxx];
  if (pid_AA ==2112&&tid_AA==2&&mid_AA==0){
 h_2d_AA->Fill(AA_z,AA_x);

     // Mark that this event has at least one qualifying hit
            condition_satisfied_AA = true;

}}
// Fill ONLY ONCE per event
    if (condition_satisfied_AA) {
	    neutron_AA_events++;
   h_ev_AA->Fill(*ev_nucl);
    }


bool condition_satisfied_AAHCal = false;
for (int jj = 0; jj < *ActAn_nhits; jj++) {
    int idxx_AA = ActAn_sdtridx[jj];
    double pid_AA = SDTrack_PID[idxx_AA];
    double tid_AA = SDTrack_TID[idxx_AA];
    double mid_AA = SDTrack_MID[idxx_AA];
    for (int kk = 0; kk < *HCal_hit_nhits; kk++) {
        int idx = HCal_sdtridx[kk];
        double pid_HCal = SDTrack_PID[idx];
        double tid_HCal = SDTrack_TID[idx];
        double mid_HCal = SDTrack_MID[idx];
        if (pid_AA == 2112 &&
            tid_AA == 2 &&
            mid_AA == 0 &&
            pid_HCal == 2112 &&
            tid_HCal == 2 &&
            mid_HCal == 0) {
            condition_satisfied_AAHCal = true;
        }
    }
}

if (condition_satisfied_AAHCal) {
    neutron_AA_HCal_events++;
    h_ev_HCalAA->Fill(*ev_nucl);

 //   std::cout << "AA-HCal neutron event: ev_nucl = "
   //           << *ev_nucl << std::endl;
}

}





TCanvas *c1 = new TCanvas("c1","AA and HCAl Neutron events",1400,1000);
c1->Divide(2,2);
c1->cd(1);
h_2d_HCal->Draw();
c1->cd(2);
h_ev->Draw();
c1->cd(3);
h_2d_AA->Draw();
c1->cd(4);
h_ev_AA->Draw();

TCanvas *c2 = new TCanvas("c2","AA-HCAl Neutron events",1000,600);
c2->Divide(1,1);
c2->cd(1);
h_ev_HCalAA->Draw();



}
