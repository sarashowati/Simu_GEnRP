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

void xg_vs_zg(){
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

 // ===== Harm HodoPR Variables =====
TTreeReaderValue<int> hodoPR_hit_nhits(tree, "Harm.PRPolScintFarSide.hit.nhits");
TTreeReaderArray<int> hodoPR_hit_plane(tree,"Harm.PRPolScintFarSide.hit.plane");
TTreeReaderArray<double> hodoPR_hit_xhitg(tree, "Harm.PRPolScintFarSide.hit.xhitg");
TTreeReaderArray<double> hodoPR_hit_yhitg(tree, "Harm.PRPolScintFarSide.hit.yhitg");
TTreeReaderArray<double> hodoPR_hit_zhitg(tree, "Harm.PRPolScintFarSide.hit.zhitg");

 // ===== Harm PRGEMs Variables =====
TTreeReaderArray<int> GEMPR_hit_plane(tree,"Harm.PRPolGEMFarSide.hit.plane");
TTreeReaderValue<int> GEMPR_hit_nhits(tree, "Harm.PRPolGEMFarSide.hit.nhits");
TTreeReaderArray<double> GEMPR_hit_xg(tree, "Harm.PRPolGEMFarSide.hit.xg");
TTreeReaderArray<double> GEMPR_hit_yg(tree, "Harm.PRPolGEMFarSide.hit.yg");
TTreeReaderArray<double> GEMPR_hit_zg(tree, "Harm.PRPolGEMFarSide.hit.zg");

 // ===== Harm PRGEMs Variables =====
TTreeReaderArray<int> CEGEMs_hit_plane(tree,"Harm.CEPolRear.hit.plane");
TTreeReaderValue<int> CEGEMs_hit_nhits(tree, "Harm.CEPolRear.hit.nhits");
TTreeReaderArray<double> CEGEMs_hit_xg(tree, "Harm.CEPolRear.hit.xg");
TTreeReaderArray<double> CEGEMs_hit_yg(tree, "Harm.CEPolRear.hit.yg");
TTreeReaderArray<double> CEGEMs_hit_zg(tree, "Harm.CEPolRear.hit.zg");

TTreeReaderValue<int> ActAnScint_hit_nhits(tree, "Harm.ActAnScint.hit.nhits");
TTreeReaderArray<double> ActAnScint_hit_xhitg(tree, "Harm.ActAnScint.hit.xhitg");
TTreeReaderArray<double> ActAnScint_hit_yhitg(tree, "Harm.ActAnScint.hit.yhitg");
TTreeReaderArray<double> ActAnScint_hit_zhitg(tree, "Harm.ActAnScint.hit.zhitg");

//HCAl
TTreeReaderValue<int> HCal_hit_nhits(tree, "Harm.HCalScint.hit.nhits");
TTreeReaderArray<double> HCal_hit_xhitg(tree, "Harm.HCalScint.hit.xhitg");
TTreeReaderArray<double> HCal_hit_yhitg(tree, "Harm.HCalScint.hit.yhitg");
TTreeReaderArray<double> HCal_hit_zhitg(tree, "Harm.HCalScint.hit.zhitg");




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


int Nentries = tree.GetEntries();
int entry_num = 0;
// ============================
// Create TGraphs (NOT TH2D)
// ============================
TGraph *gr_p12   = new TGraph();
TGraph *gr_hodo = new TGraph();
TGraph *gr_ce = new TGraph();
TGraph *gr_aa   = new TGraph();
TGraph *gr_hcal   = new TGraph();


int n1=0, nh=0, nce=0,naa=0,nhcal=0;
// ============================
// Event loop
// ============================
while(tree.Next()) {
    // ---- GEM PR ----
    for (int i = 0; i < *GEMPR_hit_nhits; i++) {
        int plane = GEMPR_hit_plane[i];
        double x  = GEMPR_hit_xg[i];
        double z  = GEMPR_hit_zg[i];
//	if(plane==1){
            gr_p12->SetPoint(n1++, z, x);
    }
    // ---- HODO ----
    for (int i = 0; i < *hodoPR_hit_nhits; i++) {
        double x = hodoPR_hit_xhitg[i];
	int plane = hodoPR_hit_plane[i];
        double z = hodoPR_hit_zhitg[i];
        gr_hodo->SetPoint(nh++, z, x);
    }
    // ---- CE GEM ----
    for (int i = 0; i < *CEGEMs_hit_nhits; i++){
	     int plane = CEGEMs_hit_plane[i];
        double x  = CEGEMs_hit_xg[i];
        double z  = CEGEMs_hit_zg[i];
//if(plane==4){
            gr_ce->SetPoint(nce++, z, x);
    }


   for (int i = 0; i < *ActAnScint_hit_nhits; i++){
//             int plane = CEGEMs_hit_plane[i];
        double x  = ActAnScint_hit_xhitg[i];
        double z  = ActAnScint_hit_zhitg[i];
//if(plane==4){
            gr_aa->SetPoint(naa++, z, x);
    }


 for (int i = 0; i < *HCal_hit_nhits; i++){
        double x  = HCal_hit_xhitg[i];
        double z  = HCal_hit_zhitg[i];
            gr_hcal->SetPoint(nhcal++, z, x);
    }


}


// ============================
// Style
// ============================
gr_p12->SetMarkerColor(kRed);
gr_p12->SetMarkerStyle(20);
gr_hodo->SetMarkerColor(kGreen+2);
gr_hodo->SetMarkerStyle(20);
gr_ce->SetMarkerColor(kMagenta);
gr_ce->SetMarkerStyle(20);
gr_aa->SetMarkerColor(kGreen);
gr_aa->SetMarkerStyle(20);
gr_hcal->SetMarkerColor(kBlue+2);
gr_hcal->SetMarkerStyle(22);



// ============================
// Draw
// ============================
TCanvas *c1 = new TCanvas("c1","Xg vs Zg overlay",800,800);
TH2F *frame = new TH2F("frame","Zg vs Xg;Zg;Xg",100,3.5,9.5,100,-5.5,0.5); //if i add HCAL
//TH2F *frame = new TH2F("frame","Zg vs Xg;Zg;Xg",100,3.5,5.5,100,-3.2,-1.2);


frame->Draw();
// Overlay scatter
gr_p12->Draw("P SAME");
gr_hodo->Draw("P SAME");
gr_ce->Draw("P SAME");
gr_aa->Draw("P SAME");
gr_hcal->Draw("P SAME");


// Legend
TLegend *leg = new TLegend(0.70,0.70,0.90,0.90);
leg->AddEntry(gr_p12,"GEM PR Plane 12","p");
leg->AddEntry(gr_hodo,"Hodo","p");
leg->AddEntry(gr_ce,"CE GEMRear Plane ","p");
leg->AddEntry(gr_aa,"AA","p");
leg->AddEntry(gr_hcal,"HCal","p");
leg->Draw();
c1->Update();


}

