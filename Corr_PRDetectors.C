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


using namespace std;
#define gemL1L2Tdiff 10.0
#define gemL1L2Pdiff 0.1


void Corr_PRDetectors(){
//void hcal_analysis(char const *input_rootfile){ //alternate to tchain
    gROOT->Reset();
    using namespace std;

//####### tchain adding #####
    TChain *tchnT = new TChain("T");
   tchnT->Add(Form("/w/halla-scshelf2102/sbs/saru/g4sbs/install/run_g4sbs_here/genrp_test_Saru.root"));
//      tchnT->Add(Form("/w/halla-scshelf2102/sbs/saru/Simu_GEnRP/genrp_agc_sim4.root"));


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

 // ===== Harm HodoPR Variables =====
TTreeReaderValue<int> hodoPR_hit_nhits(tree, "Harm.PRPolScintFarSide.hit.nhits");
TTreeReaderArray<double> hodoPR_hit_xhitg(tree, "Harm.PRPolScintFarSide.hit.xhitg");
TTreeReaderArray<double> hodoPR_hit_yhitg(tree, "Harm.PRPolScintFarSide.hit.yhitg");
TTreeReaderArray<double> hodoPR_hit_zhitg(tree, "Harm.PRPolScintFarSide.hit.zhitg");
TTreeReaderArray<double> hodoPR_hit_tmin(tree, "Harm.PRPolScintFarSide.hit.tmin");
TTreeReaderArray<double> hodoPR_hit_tmax(tree, "Harm.PRPolScintFarSide.hit.tmax");
TTreeReaderArray<double> hodoPR_hit_tavg(tree, "Harm.PRPolScintFarSide.hit.tavg");
/*
 // ===== Harm ActiveAna Variables =====
TTreeReaderValue<int> ActAna_hit_row(tree, "Harm.ActAna.hit.row");
TTreeReaderValue<int> ActAna_hit_col(tree, "Harm.ActAna.hit.col");
TTreeReaderValue<double> ActAna_hit_xhit(tree, "Harm.ActAna.hit.xhit");
TTreeReaderValue<double> ActAna_hit_yhit(tree, "Harm.ActAna.hit.yhit");
TTreeReaderValue<double> ActAna_hit_zhit(tree, "Harm.ActAna.hit.zhit");
TTreeReaderValue<double> ActAna_hit_t(tree, "Harm.ActAna.hit.t");
TTreeReaderValue<double> ActAna_hit_tmin(tree, "Harm.ActAna.hit.tmin");
TTreeReaderValue<double> ActAna_hit_tmax(tree, "Harm.ActAna.hit.tmax");

*/
 // ===== Harm PRGEMs Variables =====
TTreeReaderArray<int> GEMPR_hit_plane(tree,"Harm.PRPolGEMFarSide.hit.plane");
TTreeReaderValue<int> GEMPR_hit_nhits(tree, "Harm.PRPolGEMFarSide.hit.nhits");
TTreeReaderArray<double> GEMPR_hit_xg(tree, "Harm.PRPolGEMFarSide.hit.xg");
TTreeReaderArray<double> GEMPR_hit_yg(tree, "Harm.PRPolGEMFarSide.hit.yg");
TTreeReaderArray<double> GEMPR_hit_zg(tree, "Harm.PRPolGEMFarSide.hit.zg");
TTreeReaderArray<double> GEMPR_hit_t(tree, "Harm.PRPolGEMFarSide.hit.t");
TTreeReaderArray<double> GEMPR_hit_tmin(tree, "Harm.PRPolGEMFarSide.hit.tmin");
TTreeReaderArray<double> GEMPR_hit_tmax(tree, "Harm.PRPolGEMFarSide.hit.tmax");

//HCal Hit Variables
TTreeReaderValue<int> HCal_hit_nhits(tree, "Harm.HCalScint.hit.nhits");
TTreeReaderValue<int> ActAna_hit_nhits(tree, "Harm.ActAnScint.hit.nhits");

//Quasi elastic CUt variables
TTreeReaderValue<double>ev_W2(tree,"ev.W2");
TTreeReaderValue<double>ev_ep(tree,"ev.ep");
TTreeReaderValue<int>Earm_BBGEM_Track_ntracks(tree,"Earm.BBGEM.Track.ntracks");
TTreeReaderArray<double>Earm_BBPSTF1_det_esum(tree,"Earm.BBPSTF1.det.esum");
TTreeReaderArray<double>Earm_BBSHTF1_det_esum(tree,"Earm.BBSHTF1.det.esum");


TTreeReaderValue<int>Earm_BBPS_nhits(tree,"Earm.BBPSTF1.hit.nhits");
TTreeReaderValue<int>Earm_BBSH_nhits(tree,"Earm.BBSHTF1.hit.nhits");
TTreeReaderArray<double>Earm_BBPS(tree,"Earm.BBPSTF1.hit.tavg");
TTreeReaderArray<double>Earm_BBSH(tree,"Earm.BBSHTF1.hit.tavg");
TTreeReaderArray<double> ActAna_hit_t(tree, "Harm.ActAnScint.hit.tavg");
TTreeReaderValue<int>Earm_BBHodo_nhits(tree,"Earm.BBHodoScint.hit.nhits");
TTreeReaderArray<double>Earm_BBHodo_tavg(tree,"Earm.BBHodoScint.hit.tavg");

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
    gStyle->SetStatFont(62);
    gStyle->SetPadLeftMargin(0.05); //best setting previously for (3,2)
    gStyle->SetPadRightMargin(0.05); ////best setting previously for (3,2)
    gStyle->SetPadBottomMargin(0.10);
    gStyle->SetStatStyle(0); //set statbox transparent
    gStyle->SetTitleStyle(0);
    TGaxis::SetMaxDigits(4);
    gROOT->ForceStyle(true);


// ===== 1D Histograms =====
TH2D *h_xgxg_p1p2 = new TH2D("h_xgxg_p1p2","PR GEM Hit Xg vs Hit Xg for Plane12;Hitxg_p1;HitXg_p2",100,-3,-2.5,100,-3,-2.5);
TH2D *h_ygyg_p1p2 = new TH2D("h_ygyg_p1p2","PR GEM Hit Yg vs Hit Yg for Plane12;Hityg_p1;Hityg_p2",100,-1.25,1.25,100,-1.25,1.25);
TH2D *h_zgzg_p1p2 = new TH2D("h_zgzg_p1p2","PR GEM Hit Zg vs Hit Zg for Plane12;Hitzg_p1;Hitzg_p2",100,4.3,5.2,100,4.3,5.2);
//1d plot
TH1D *h_xgdelta = new TH1D("h_xgdelta", "X coordinate hit difference between GEM plane 1 GEM Plane2;GEM xg hit plane1-GEM xg hit plane2;counts",100,-1,1);
TH1D *h_ygdelta = new TH1D("h_ygdelta", "Y coordinate hit difference between GEM plane 1 GEM Plane2;GEM yg hit plane1-GEM yg hit plane2;counts",100,-1,1);
TH1D *h_zgdelta = new TH1D("h_zgdelta", "Z coordinate hit difference between GEM plane 1 GEM Plane2;GEM zg hit plane1-GEM zg hit plane2;counts",100,-1,1);

TH1D *h_T1T2 = new TH1D("h_T1T2","TIme difference between GEM plane 1 GEM Plane2;GEM Plane1 TIme-GEM  plane2 TIme;counts",100,-10,10);
TH2D *h_T1T2_p1p2 = new TH2D("h_T1T2_p1p2","PR GEMs Plane1 TIme vs PR GEM Plane2 TIme;GEM_Time_Plane1;GEM_Time_Plane2",100,0,1000,100,0,1000);


TH1D *h_NewCor_delta = new TH1D("h_NewCor_delta", "#sqrt{x_{g}^{2}+z_{g}^{2}};#sqrt{x_{g}^{2}+z_{g}^{2}} (m);Counts",100, -1, 1);
TH2D *h_NewCor_p1p2 = new TH2D("h_NewCor_p1p2", "Plane 1 vs Plane 2 #sqrt{x_{g}^{2}+z_{g}^{2}};Plane 1 #sqrt{x_{g}^{2}+z_{g}^{2}} (m);Plane 2 #sqrt{x_{g}^{2}+z_{g}^{2}} (m)",100, 5.1,5.9,100, 5.1, 5.9);
/*
TH2D *h_xgxg_p1p2 = new TH2D("h_xgxg_p1p2","PR GEM Hit Xg vs Hit Xg for Plane12;Hitx_p1;HitX_p2",100,-1.1,1.1,100,-1.1,-1.1);
TH2D *h_ygyg_p1p2 = new TH2D("h_ygyg_p1p2","PR GEM Hit Yg vs Hit Yg for Plane12;Hity_p1;Hity_p2",100,-0.60,-0.6015,100,-0.50,-0.5015);
TH2D *h_zgzg_p1p2 = new TH2D("h_zgzg_p1p2","PR GEM Hit Xg vs Hit Xg for Plane12;Hitz_p1;Hitz_p2",100,-0.4,0.4,100,-0.4,0.4);
*/

TH2D *h_yg_GemHodo = new TH2D("h_yg_GemHodo","PR GEMs Hit vs HodoPR Hit;GEMHityg;HodoPR_yg",100,-1.1,1.1, 100, -1.1,1.1);
TH2D *h_xg_GemHodo = new TH2D("h_xg_GemHodo","PR GEMs Hit vs HodoPR Hit;GEMHitxg;HodoPR_xg",100,-3.1,-2.6, 100, -3.1,-2.7);
TH2D *h_zg_GemHodo = new TH2D("h_zg_GemHodo","PR GEMs Hit vs HodoPR Hit;GEMHitzg;HodoPR_zg",100,4.3,5.2, 100, 4.3,5.2);

//TH2D *h_yg_GemHodo = new TH2D("h_yg_GemHodo","PR GEMs Hit vs HodoPR Hit;GEMHityg;HodoPR_yg",100,4.3,5.1, 100, 4.3,5.1);
TH2D *h_time_GemHodo = new TH2D("h_time_GemHodo", "GEMs Hit_t vs HodoPR _hit_tavg; GEM_hit_t; HodoPR_hit_tavg",100,0,1000,100,0,1000); 
TH1D *h_deltaT = new TH1D("h_deltaT","GEM time-hodo time;GEM time-hodo time;counts",100,-20,20);
TH1D *h_deltahit_yg = new TH1D("h_deltahit_yg","GEM yg hit position-hodo hit position;GEM yg hit position-hodo yg hit position;counts",100,-1.1,1.1);
TH1D *h_deltahit_xg = new TH1D("h_deltahit_xg","GEM xg hit position-hodo xg hit position;GEM xg hit position-hodo hit xg position;counts",100,-1.1,1.1);
TH1D *h_deltahit_zg = new TH1D("h_deltahit_zg","GEM zg hit position-hodo zg hit position;GEM zg hit position-hodo hit zg position;counts",100,-1.1,1.1);



// ===== 2D Histograms for each plane =====
const int NPLANES = 2;
// ===== Event Loop =====
int Nentries = tree.GetEntries();
int entry_num = 0;
while(tree.Next()) {
    double x1=-999, y1=-999, z1=-999;
    double x2=-999, y2=-999, z2=-999;

    double New_Cor1 = 0.0;
    double New_Cor2 = 0.0;

    double T1 = 0.0;
    double T2 = 0.0;

    bool hasP1 = false;
    bool hasP2 = false;
    // Loop over hits
    for (int i = 0; i < *GEMPR_hit_nhits; i++) {
        int plane = GEMPR_hit_plane[i];
        double x = GEMPR_hit_xg[i];
        double y = GEMPR_hit_yg[i];
        double z = GEMPR_hit_zg[i];
	double T = GEMPR_hit_t[i];
	double New_Cor = sqrt(x*x + z*z);
        // Plane 1
        if(plane == 1) {
            x1 = x;
            y1 = y;
            z1 = z;
	    T1=T;
	     New_Cor1 = New_Cor;
            hasP1 = true;
        }
        // Plane 2
        if(plane == 2) {
            x2 = x;
            y2 = y;
            z2 = z;
	    T2=T;
	     New_Cor2 = New_Cor;
            hasP2 = true; }
    }
    // Fill the hit correlations between two GEM planes only if both planes exist
    if(hasP1 && hasP2) {
	h_xgdelta->Fill(x1-x2);
        h_xgxg_p1p2->Fill(x1, x2);
	h_ygdelta->Fill(y1-y2);
        h_ygyg_p1p2->Fill(y1, y2);
        h_zgdelta->Fill(z1-z2);
        h_zgzg_p1p2->Fill(z1, z2);

	h_T1T2->Fill(T1-T2);
	h_T1T2_p1p2->Fill(T1,T2);

    h_NewCor_delta->Fill(New_Cor1 - New_Cor2);
    h_NewCor_p1p2->Fill(New_Cor1, New_Cor2);

    }



//FOr COrrelation Plot between  PRGEMS and HodoPR in terms of Hit Position
for (int i = 0; i<*GEMPR_hit_nhits; i++) {
	 int plane = GEMPR_hit_plane[i];
	double gemhit_yg = GEMPR_hit_yg[i];
for (int j= 0; j<*hodoPR_hit_nhits;j++){
	double hodoPR_yg = hodoPR_hit_yhitg[j];
	  // Plane 1
//        if(plane == 1) {

if(*ActAna_hit_nhits>0 && *HCal_hit_nhits>0 ){      
 if (fabs(*ev_W2 - 0.86) < 0.7 && Earm_BBPSTF1_det_esum[0] > 0.1){
    if (*Earm_BBGEM_Track_ntracks >= 0 && ((Earm_BBPSTF1_det_esum[0] + Earm_BBSHTF1_det_esum[0]) / (*ev_ep)) > 0.1){ 
    
	h_deltahit_yg->Fill(GEMPR_hit_yg[i] - hodoPR_hit_yhitg[j]);
	h_yg_GemHodo->Fill(gemhit_yg,hodoPR_yg);}}}}}

for (int i = 0; i<*GEMPR_hit_nhits; i++) {
	int plane = GEMPR_hit_plane[i];
        double gemhit_xg = GEMPR_hit_xg[i];
for (int j= 0; j<*hodoPR_hit_nhits;j++){
        double hodoPR_xg = hodoPR_hit_xhitg[j];
	 // Plane 1
    //    if(plane == 1){
    if(*ActAna_hit_nhits>0 && *HCal_hit_nhits>0 ){
 if (fabs(*ev_W2 - 0.86) < 0.7 && Earm_BBPSTF1_det_esum[0] > 0.1){
    if (*Earm_BBGEM_Track_ntracks >= 0 && ((Earm_BBPSTF1_det_esum[0] + Earm_BBSHTF1_det_esum[0]) / (*ev_ep)) > 0.1){
        h_deltahit_xg->Fill(GEMPR_hit_xg[i] - hodoPR_hit_xhitg[j]);
        h_xg_GemHodo->Fill(gemhit_xg,hodoPR_xg);
}}}}}

for (int i = 0; i<*GEMPR_hit_nhits; i++) {
	int plane = GEMPR_hit_plane[i];
        double gemhit_zg = GEMPR_hit_zg[i];
for (int j= 0; j<*hodoPR_hit_nhits;j++){
        double hodoPR_zg = hodoPR_hit_zhitg[j];
	 // Plane 1
     //   if(plane == 1) {
    if(*ActAna_hit_nhits>0 && *HCal_hit_nhits>0 ){
 if (fabs(*ev_W2 - 0.86) < 0.7 && Earm_BBPSTF1_det_esum[0] > 0.1){
    if (*Earm_BBGEM_Track_ntracks >= 0 && ((Earm_BBPSTF1_det_esum[0] + Earm_BBSHTF1_det_esum[0]) / (*ev_ep)) > 0.1){
        h_deltahit_zg->Fill(GEMPR_hit_zg[i] - hodoPR_hit_zhitg[j]);
        h_zg_GemHodo->Fill(gemhit_zg,hodoPR_zg);
}}}}}




//Timing correlation
for (int i = 0; i<*GEMPR_hit_nhits; i++) {
	 int plane = GEMPR_hit_plane[i];
	double gemhit_t = GEMPR_hit_t[i];
for (int j= 0; j<*hodoPR_hit_nhits;j++){
        double hodoPR_tavg = hodoPR_hit_tavg[j];
//	if(plane == 1) {
if(*ActAna_hit_nhits>0 && *HCal_hit_nhits>0 ){
 if (fabs(*ev_W2 - 0.86) < 0.7 && Earm_BBPSTF1_det_esum[0] < 0.1){
    if (*Earm_BBGEM_Track_ntracks <= 0 && ((Earm_BBPSTF1_det_esum[0] + Earm_BBSHTF1_det_esum[0]) / (*ev_ep)) < 0.1){
	h_deltaT->Fill(GEMPR_hit_t[i] - hodoPR_hit_tavg[j]);
        h_time_GemHodo->Fill(gemhit_t,hodoPR_tavg);
    }}}}}

}
/*
// Canvas for 1D histograms
TCanvas *c1 = new TCanvas("c1","PR GEMs hit vs HodoPR hit",1200,600);
c1->Divide(3,2);
c1->cd(1);
h_xgdelta->Draw();
c1->cd(2);
h_ygdelta->Draw();
c1->cd(3);
h_zgdelta->Draw();
c1->cd(4);
h_xgxg_p1p2->GetXaxis()->SetTitle("GEM x_{g} (m) for Plane1");
h_xgxg_p1p2->GetYaxis()->SetTitle("GEM x_{g} (m) for plane2");
// Optional title centering
h_xgxg_p1p2->GetXaxis()->CenterTitle();
h_xgxg_p1p2->GetYaxis()->CenterTitle();
h_xgxg_p1p2->Draw("COLZ");
c1->cd(5);
h_ygyg_p1p2->GetXaxis()->SetTitle("GEM y_{g} (m) for Plane1");
h_ygyg_p1p2->GetYaxis()->SetTitle("GEM y_{g} (m) for plane2");
// Optional title centering
h_ygyg_p1p2->GetXaxis()->CenterTitle();
h_ygyg_p1p2->GetYaxis()->CenterTitle();
h_ygyg_p1p2->Draw("COLZ");
c1->cd(6);
h_zgzg_p1p2->GetXaxis()->SetTitle("GEM z_{g} (m) for Plane1");
h_zgzg_p1p2->GetYaxis()->SetTitle("GEM z_{g} (m) for plane2");
// Optional title centering
h_zgzg_p1p2->GetXaxis()->CenterTitle();
h_zgzg_p1p2->GetYaxis()->CenterTitle();
h_zgzg_p1p2->Draw("COLZ");

// Canvas for 1D histograms
TCanvas *c0 = new TCanvas("c0","PR GEMs hits for Plane1 vs Plane2",1000,400);
c0->Divide(2,1);
c0->cd(1);
h_T1T2->Draw();
c0->cd(2);
h_T1T2_p1p2->GetXaxis()->SetTitle("GEM Time for Plane 1");
h_T1T2_p1p2->GetYaxis()->SetTitle("GEM Time for Plane 2");
// Optional title centering
h_T1T2_p1p2->GetXaxis()->CenterTitle();
h_T1T2_p1p2->GetYaxis()->CenterTitle();
h_T1T2_p1p2->Draw("COLZ");



// Canvas for 1D histograms
TCanvas *c11 = new TCanvas("c11","PR GEMs hits for Plane1 vs Plane2",1000,400);
c11->Divide(2,1);
c11->cd(1);
h_NewCor_delta->Draw();
c11->cd(2);
h_NewCor_p1p2->GetXaxis()->SetTitle("#sqrt{x_{g}^{2}+z_{g}^{2}} for Plane 1");
h_NewCor_p1p2->GetYaxis()->SetTitle("#sqrt{x_{g}^{2}+z_{g}^{2}} for Plane 2");
// Optional title centering
h_NewCor_p1p2->GetXaxis()->CenterTitle();
h_NewCor_p1p2->GetYaxis()->CenterTitle();
h_NewCor_p1p2->Draw("COLZ");



// Canvas for 1D histograms
TCanvas *c = new TCanvas("c","PR GEMs hit vs HodoPR hit", 1400,800);
c->Divide(3,2);
c->cd(1);
h_deltahit_xg->Draw();
c->cd(2);
h_deltahit_yg->Draw();
c->cd(3);
h_deltahit_zg->Draw();
//gPad->SetLogy();
c->cd(4);
// Axis titles
h_xg_GemHodo->GetXaxis()->SetTitle("GEM x_{g} (m)");
h_xg_GemHodo->GetYaxis()->SetTitle("HodoPR x_{hitg} (m)");
// Optional title centering
h_xg_GemHodo->GetXaxis()->CenterTitle();
h_xg_GemHodo->GetYaxis()->CenterTitle();
h_xg_GemHodo->Draw("COLZ");
c->cd(5);
h_yg_GemHodo->GetXaxis()->SetTitle("GEM y_{g} (m)");
h_yg_GemHodo->GetYaxis()->SetTitle("HodoPR y_{hitg} (m)");
// Optional title centering
h_yg_GemHodo->GetXaxis()->CenterTitle();
h_yg_GemHodo->GetYaxis()->CenterTitle();
h_yg_GemHodo->Draw("COLZ");
c->cd(6);
h_zg_GemHodo->GetXaxis()->SetTitle("GEM z_{g} (m)");
h_zg_GemHodo->GetYaxis()->SetTitle("HodoPR z_{hitg} (m)");
// Optional title centering
h_zg_GemHodo->GetXaxis()->CenterTitle();
h_zg_GemHodo->GetYaxis()->CenterTitle();
h_zg_GemHodo->Draw("COLZ");
*/


TCanvas *c3 = new TCanvas("c3","PR GEMs vs HodoPR Time COrrelation",800,400);
c3->Divide(2,1);
c3->cd(1);
h_deltaT->Draw();
c3->cd(2);
h_time_GemHodo->GetXaxis()->SetTitle("GEM hit_t");
h_time_GemHodo->GetYaxis()->SetTitle("HodoPR hit_tavg");
// Optional title centering
h_time_GemHodo->GetXaxis()->CenterTitle();
h_time_GemHodo->GetYaxis()->CenterTitle();
h_time_GemHodo->Draw("COLZ");


}
