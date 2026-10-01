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

void TOF_Corr_with_BBPR(){
//void hcal_analysis(char const *input_rootfile){ //alternate to tchain
    gROOT->Reset();
    using namespace std;

//####### tchain adding #####
    TChain *tchnT = new TChain("T");
//    tchnT->Add(Form("/lustre24/expphy/volatile/halla/sbs/saru/genrp_replayed_1125_500k_event*root"));
      tchnT->Add(Form("/w/halla-scshelf2102/sbs/saru/Simu_GEnRP/genrp_agc_sim4.root"));


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
    TFile *fout = new TFile("./Simu_genrp_PR.root","RECREATE");


//Quasi elastic CUt variables
TTreeReaderValue<double>ev_W2(tree,"ev.W2");
TTreeReaderValue<double>ev_ep(tree,"ev.ep");
TTreeReaderValue<int>Earm_BBGEM_Track_ntracks(tree,"Earm.BBGEM.Track.ntracks");
TTreeReaderArray<double>Earm_BBPS_esum(tree,"Earm.BBPSTF1.det.esum");
TTreeReaderArray<double>Earm_BBSH_esum(tree,"Earm.BBSHTF1.det.esum");
TTreeReaderArray<int>Earm_BBPS_cell(tree,"Earm.BBPSTF1.hit.cell");
TTreeReaderArray<int>Earm_BBSH_cell(tree,"Earm.BBSHTF1.hit.cell");



 // ===== Harm HodoPR Variables =====
TTreeReaderValue<int> hodoPR_nhits(tree, "Harm.PRPolScintFarSide.hit.nhits");
TTreeReaderArray<double> hodoPR_hit_tmin(tree, "Harm.PRPolScintFarSide.hit.tmin");
TTreeReaderArray<double> hodoPR_hit_tmax(tree, "Harm.PRPolScintFarSide.hit.tmax");
TTreeReaderArray<double> hodoPR_t(tree, "Harm.PRPolScintFarSide.hit.tavg");
TTreeReaderArray<int> hodoPR_cell(tree, "Harm.PRPolScintFarSide.hit.cell");
TTreeReaderArray<double> hodoPR_esum(tree, "Harm.PRPolScintFarSide.hit.sumedep");



 // ===== Harm PRGEMs Variables =====
TTreeReaderArray<int> GEMPR_hit_plane(tree,"Harm.PRPolGEMFarSide.hit.plane");
TTreeReaderValue<int> GEMPR_hit_nhits(tree, "Harm.PRPolGEMFarSide.hit.nhits");
TTreeReaderArray<double> GEMPR_hit_t(tree, "Harm.PRPolGEMFarSide.hit.t");
TTreeReaderArray<double> GEMPR_hit_tmin(tree, "Harm.PRPolGEMFarSide.hit.tmin");
TTreeReaderArray<double> GEMPR_hit_tmax(tree, "Harm.PRPolGEMFarSide.hit.tmax");
TTreeReaderArray<int> GEMPR_hit_pid(tree,"Harm.PRPolGEMFarSide.hit.pid");

//HCal  and AA Hit_Timing Variables
TTreeReaderValue<int> HCal_hit_nhits(tree, "Harm.HCalScint.hit.nhits");
TTreeReaderArray<int> HCal_cell(tree, "Harm.HCalScint.hit.cell");
TTreeReaderArray<double> HCal_t(tree, "Harm.HCalScint.hit.tmin");
TTreeReaderArray<double> HCal_esum(tree, "Harm.HCalScint.hit.sumedep");


TTreeReaderArray<double> ActAna_esum(tree,"Harm.ActAnScint.hit.sumedep");
TTreeReaderValue<int> ActAna_nhits(tree, "Harm.ActAnScint.hit.nhits");
TTreeReaderArray<double> ActAna_t(tree, "Harm.ActAnScint.hit.tmin");
TTreeReaderArray<int> ActAna_cell(tree, "Harm.ActAnScint.hit.cell");

//BB Side hit timinig variables
TTreeReaderArray<int>Earm_BBGEM_pid(tree,"Earm.BBGEM.hit.pid");
TTreeReaderValue<int>Earm_BBGEM_nhits(tree,"Earm.BBGEM.hit.nhits");
TTreeReaderValue<int>Earm_BBPS_nhits(tree,"Earm.BBPSTF1.hit.nhits");
TTreeReaderValue<int>Earm_BBSH_nhits(tree,"Earm.BBSHTF1.hit.nhits");
TTreeReaderArray<double>Earm_BBPS_t(tree,"Earm.BBPSTF1.hit.tavg");
TTreeReaderArray<double>Earm_BBSH_t(tree,"Earm.BBSHTF1.hit.tavg");
TTreeReaderValue<int>Earm_BBHodo_nhits(tree,"Earm.BBHodoScint.hit.nhits");
TTreeReaderArray<double>Earm_BBHodo_t(tree,"Earm.BBHodoScint.hit.tavg");
TTreeReaderArray<double>Earm_BBHodo_esum(tree,"Earm.BBHodoScint.hit.sumedep");


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



//TOF    1D plots
TH1D *h_TOF = new TH1D("h_TOF","BBHodo time-hodoPR time;BBHodo time-hodoPR time;counts",200,-100,100);
TH1D *h_TOF_BBHCal = new TH1D("h_TOF_BBHCal","HCal time -BBhodo_time;HCal Time- BBhodo Time;count",200,-100,100);
TH1D *h_TOF_HodoHCal = new TH1D("h_TOF_HodoHCal","HCal Time- HodoPR Time;HCal Time-HodoPR Time;counts",200,-100,100);
TH1D *h_TOF_AA = new TH1D("h_TOF_AA","BBHodo time-AA time;BBHodo time-AA time;counts",200,-100,100);
TH1D *h_TOF_AAHCal = new TH1D("h_TOF_AAHCal","AA time - HCal time;AA time-HCal time;counts",200,-100,100);
TH1D *h_TOF_AAHodoPR = new TH1D("h_TOF_AAHodoPR","HodoPR time - AA time;HodoPR time - AA time;counts",200,-100,100);


// TOF 2D plots
TH2D *h_TOF_Corr = new TH2D("h_TOF_Corr","BBHodo time - hodoPR time vs HodoPR_Cell;HodoPR_Cell;BBHodo time - hodoPR time",28,-2,26,200,-100.0,100.0);
TProfile *ph_TOF_Corr = new TProfile( "ph_TOF_Corr","",28, -2, 26);
TH2D *h_TOF_Corr_HCal = new TH2D("h_TOF_Corr_HCal","BBHodo time - hcaltime vs Hcal_Cell;HCal_Cell; Hcal Time - BBHodo Time",300,-5,295,200,-100.0,100.0);
TProfile *ph_TOF_Corr_HCal = new TProfile( "ph_TOF_Corr_HCal","",300, 0, 295);
TH2D *h_TOF_Hodo_HCal = new TH2D("h_TOF_Hodo_HCal","HodoPR time - hcaltime vs Hcal_Cell;HCal_Cell; Hcal Time - HodoPR Time",300,-5,295,200,-100.0,100.0);
TH2D *h_TOF_HodoPR = new TH2D("h_TOF_HodoPR","Hcaltime - HodoPR time  vs HodoPR_Cell;HodoPR_Cell; Hcal time - HodoPR time",28,-2,26,200,-100.0,100.0);
TProfile *ph_TOF_HodoPR = new TProfile( "ph_TOF_HodoPR","",28,-2, 26);
TH2D *h_TOF_Corr_AA = new TH2D("h_TOF_Corr_AA","BBHodo time -AA time vs AA_Cell;BBSH_Cell;BBHodo time -AA time",35,-2,33,200,-100,100);
TH2D *h_TOF_Corr_AAHCal = new TH2D("h_TOF_Corr_AAHCal","AA time - HCal time vs HCal_Cell;HCal_Cell;AA time - HCal time",300,-5,295,200,-100,100);
TH2D *h_TOF_Corr_AAHodoPR = new TH2D("h_TOF_Corr_AAHodoPR","HodoPR time - AA time vs AA_Cell;AA_Cell;HodoPR time - AA time",35,-2,33,200,-100,100);



int Nentries = tree.GetEntries();
while (tree.Next()) {
/*
for (int i = 0; i < Earm_BBHodo_t.GetSize(); i++) {
    double bb_t = Earm_BBHodo_t[i];
    double bestDiff = 1e9;
    int best_j=-1;
    // --- find best PR match for THIS BB hit ---
    for (int j = 0; j < hodoPR_t.GetSize(); j++) {
        double diff = fabs(bb_t - hodoPR_t[j]);
//	h_TOF->Fill(diff);
        if (diff < bestDiff) {
            bestDiff = diff;
            best_j = j;
        }
    }
    // --- apply matching cut ---
    if (best_j < 0) continue;
    if (bestDiff > tcut) continue;
    double pr_t = hodoPR_t[best_j];
    double cell = hodoPR_cell[best_j];
    double dT = bb_t - pr_t;
    h_TOF->Fill(dT);
    h_TOF_Corr->Fill(cell, dT);
    ph_TOF_Corr->Fill(cell, dT);
}
*/
for (int i = 0; i < *Earm_BBHodo_nhits; i++) {
    double bb_t = Earm_BBHodo_t[i];
    double BBHodo_esum =Earm_BBHodo_esum[i];
     for (int j = 0; j < *hodoPR_nhits; j++) {
	     double HodoPR = hodoPR_t[j];
	      double cell = hodoPR_cell[j];
	      double HodoPR_esum= hodoPR_esum[j];
        double diff = (bb_t - HodoPR);
if(BBHodo_esum>0.002 && HodoPR_esum>0.003 && HodoPR<1000 && bb_t<1000){
		  h_TOF->Fill(diff);
    h_TOF_Corr->Fill(cell, diff);
    ph_TOF_Corr->Fill(cell, diff);
	}  

     }}


//For HCAL
for (int i = 0; i < *Earm_BBHodo_nhits; i++) {
    double bb_t = Earm_BBHodo_t[i];
    double BBHodo_esum =Earm_BBHodo_esum[i];
    for (int j = 0; j < *HCal_hit_nhits; j++) {
   	    double hcal_t = HCal_t[j];
	    double hcal_esum =HCal_esum[j];
    double diff = (bb_t - hcal_t);
// h_TOF_BBHCal->Fill(diff);
  double hcal_cell = HCal_cell[j];
 if(hcal_esum>0.008 && hcal_t<1000&& BBHodo_esum>0.002 && bb_t<1000 ){
	  h_TOF_BBHCal->Fill(diff);
 h_TOF_Corr_HCal->Fill(hcal_cell, diff);
    ph_TOF_Corr_HCal->Fill(hcal_cell, diff);

 }    }}
 

//For HodoPR_HCaaaaaaal
for (int i = 0; i < *hodoPR_nhits; i++) {
    double HodoPR = hodoPR_t[i];
     double HodoPR_esum= hodoPR_esum[i];
    for (int j = 0; j < *HCal_hit_nhits; j++) {
	    double hcal_t = HCal_t[j];
	     double hcal_esum =HCal_esum[j];
        double diff = (HodoPR - hcal_t);
//    h_TOF_HodoHCal->Fill(diff);
 double hcal_cell = HCal_cell[j];
 if (HodoPR_esum>0.003 && HodoPR<1000 && hcal_esum>0.008&&hcal_t<1000){
	  h_TOF_HodoHCal->Fill(diff);
 h_TOF_Hodo_HCal->Fill(hcal_cell, diff);}}
}

//HodoPR_cell vs HodoPR_t- HCal  time
for (int i = 0; i < *hodoPR_nhits; i++) {
    double HodoPR = hodoPR_t[i];
    double HodoPR_cell = hodoPR_cell[i];
    double HodoPR_esum= hodoPR_esum[i];
    for (int j = 0; j < *HCal_hit_nhits; j++) {
            double hcal_t = HCal_t[j];  
  double hcal_esum =HCal_esum[j];
        double diff = (HodoPR - hcal_t);
//    h_TOF_HodoHCal->Fill(diff);
 if(HodoPR_esum>0.003 &&  HodoPR<1000 && hcal_esum>0.008&&hcal_t<1000){
	  h_TOF_HodoHCal->Fill(diff);
 h_TOF_HodoPR->Fill(HodoPR_cell, diff);
 ph_TOF_HodoPR->Fill(HodoPR_cell,diff);
 }}}


//BBhodo to AA
for (int i = 0; i < *Earm_BBHodo_nhits; i++) {
    double bb_t = Earm_BBHodo_t[i];
    double BBHodo_esum =Earm_BBHodo_esum[i];
     for (int j = 0; j < *ActAna_nhits; j++) {
             double ActAna = ActAna_t[j];
              double cell = ActAna_cell[j];
              double Actana_esum= ActAna_esum[j];
        double diff = (bb_t - ActAna);
if(BBHodo_esum>0.002 && Actana_esum>0.004 && bb_t<1000&&ActAna<1000){
                  h_TOF_AA->Fill(diff);
    h_TOF_Corr_AA->Fill(cell, diff);}}}

    

//AA with HCal
     for (int i = 0; i< *ActAna_nhits; i++) {
             double ActAna = ActAna_t[i];
              double Actana_esum= ActAna_esum[i];
	      for (int j = 0; j < *HCal_hit_nhits; j++) {
            double hcal_t = HCal_t[j];
  double hcal_esum =HCal_esum[j];
  double diff = ( ActAna-hcal_t);
  double cell = HCal_cell[j];
  if( Actana_esum>0.004 && ActAna<1000 &&hcal_esum>0.008 && hcal_t<1000){
	  h_TOF_AAHCal->Fill(diff);
	  h_TOF_Corr_AAHCal->Fill(cell,diff);}}}
	  

//AA with HoDoPR
     for (int i = 0; i< *ActAna_nhits; i++) {
             double ActAna = ActAna_t[i];
              double Actana_esum= ActAna_esum[i];
	        double cell = ActAna_cell[i];
	      for (int j = 0; j < *hodoPR_nhits; j++) {
    double HodoPR = hodoPR_t[j];
    double HodoPR_esum= hodoPR_esum[j];
      double diff = ( HodoPR - ActAna);
    if (Actana_esum>0.004 && HodoPR_esum>0.003&& HodoPR<1000&&ActAna<1000){
	    h_TOF_AAHodoPR->Fill(diff);
	    h_TOF_Corr_AAHodoPR->Fill(cell,diff);}}}

}


TCanvas *c1 = new TCanvas("c1","1D of Time Difference plot",1200,400);
c1->Divide(3,1);
c1->cd(1);
h_TOF->Draw();
c1->cd(2);
h_TOF_BBHCal->Draw();
c1->cd(3);
h_TOF_HodoHCal->Draw();

TCanvas *c = new TCanvas("c","BBHodo with PRHodo",1200,600);
c->Divide(2,1);
c->cd(1);
gPad->SetRightMargin(0.15);
h_TOF_Corr->GetXaxis()->SetTitle("HodoPR_Cell");
h_TOF_Corr->GetYaxis()->SetTitle("BBHodo Time - HodoPR Time");
h_TOF_Corr->GetXaxis()->CenterTitle();
h_TOF_Corr->GetYaxis()->CenterTitle();
h_TOF_Corr->Draw("COLZ");
c->cd(2);
h_TOF_Corr_HCal->GetXaxis()->SetTitle("HCal_Cell");
h_TOF_Corr_HCal->GetYaxis()->SetTitle("BBHodo TIme - HCal Time");
h_TOF_Corr_HCal->GetXaxis()->CenterTitle();
h_TOF_Corr_HCal->GetYaxis()->CenterTitle();
h_TOF_Corr_HCal->Draw("COLZ");


TCanvas *c0 = new TCanvas("c0","HodoPR with HCal",1200,600);
c0->Divide(2,1);
c0->cd(1);
h_TOF_Hodo_HCal->GetXaxis()->SetTitle("HCal_Cell");
h_TOF_Hodo_HCal->GetYaxis()->SetTitle("HodoPR Time - HCal Time");
h_TOF_Hodo_HCal->GetXaxis()->CenterTitle();
h_TOF_Hodo_HCal->GetYaxis()->CenterTitle();
h_TOF_Hodo_HCal->Draw("COLZ");
c0->cd(2);
h_TOF_HodoPR->GetXaxis()->SetTitle("HodoPR_Cell");
h_TOF_HodoPR->GetYaxis()->SetTitle("HCal Time - HodoPR Time");
h_TOF_HodoPR->GetXaxis()->CenterTitle();
h_TOF_HodoPR->GetYaxis()->CenterTitle();
h_TOF_HodoPR->Draw("COLZ");

TCanvas *c11 = new TCanvas("c11","AA with BBhodo TOF plot",1200,600);
c11->Divide(2,1);
c11->cd(1);
h_TOF_AA->Draw();
c11->cd(2);
h_TOF_Corr_AA->GetXaxis()->SetTitle("AA_Cell");
h_TOF_Corr_AA->GetYaxis()->SetTitle("BBHodo Time - AA Time");
h_TOF_Corr_AA->GetXaxis()->CenterTitle();
h_TOF_Corr_AA->GetYaxis()->CenterTitle();
h_TOF_Corr_AA->Draw("COLZ");

TCanvas *c10 = new TCanvas("c10","AA with HCal TOF plot",1200,600);
c10->Divide(2,1);
c10->cd(1);
h_TOF_AAHCal->Draw();
c10->cd(2);
h_TOF_Corr_AAHCal->GetXaxis()->SetTitle("HCal_Cell");
h_TOF_Corr_AAHCal->GetYaxis()->SetTitle("AA Time-HCal Time");
h_TOF_Corr_AAHCal->GetXaxis()->CenterTitle();
h_TOF_Corr_AAHCal->GetYaxis()->CenterTitle();
h_TOF_Corr_AAHCal->Draw("COLZ");


TCanvas *c2 = new TCanvas("c2","AA with HodoPR TOF plot",1200,600);
c2->Divide(2,1);
c2->cd(1);
h_TOF_AAHodoPR->Draw();
c2->cd(2);
h_TOF_Corr_AAHodoPR->GetXaxis()->SetTitle("AA_Cell");
h_TOF_Corr_AAHodoPR->GetYaxis()->SetTitle("HodoPR Time - AA Time");
h_TOF_Corr_AAHodoPR->GetXaxis()->CenterTitle();
h_TOF_Corr_AAHodoPR->GetYaxis()->CenterTitle();
h_TOF_Corr_AAHodoPR->Draw("COLZ");


}

