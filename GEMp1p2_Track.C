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



void GEMp1p2_Track()
{
     gROOT->Reset();
    using namespace std;
//####### tchain adding #####
    TChain *tchnT = new TChain("T");
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
TTreeReaderArray<double> hodoPR_hit_xhitg(tree, "Harm.PRPolScintFarSide.hit.xhitg");
TTreeReaderArray<double> hodoPR_hit_yhitg(tree, "Harm.PRPolScintFarSide.hit.yhitg");
TTreeReaderArray<double> hodoPR_hit_zhitg(tree, "Harm.PRPolScintFarSide.hit.zhitg");
TTreeReaderArray<double> hodoPR_hit_tmin(tree, "Harm.PRPolScintFarSide.hit.tmin");
TTreeReaderArray<double> hodoPR_hit_tmax(tree, "Harm.PRPolScintFarSide.hit.tmax");
TTreeReaderArray<double> hodoPR_hit_tavg(tree, "Harm.PRPolScintFarSide.hit.tavg");

 // ===== Harm PRGEMs Variables =====
TTreeReaderArray<int> GEMPR_hit_plane(tree,"Harm.PRPolGEMFarSide.hit.plane");
TTreeReaderValue<int> GEMPR_hit_nhits(tree, "Harm.PRPolGEMFarSide.hit.nhits");
TTreeReaderArray<double> GEMPR_hit_xg(tree, "Harm.PRPolGEMFarSide.hit.xg");
TTreeReaderArray<double> GEMPR_hit_yg(tree, "Harm.PRPolGEMFarSide.hit.yg");
TTreeReaderArray<double> GEMPR_hit_zg(tree, "Harm.PRPolGEMFarSide.hit.zg");
TTreeReaderArray<double> GEMPR_hit_t(tree, "Harm.PRPolGEMFarSide.hit.t");
TTreeReaderArray<double> GEMPR_hit_tmin(tree, "Harm.PRPolGEMFarSide.hit.tmin");
TTreeReaderArray<double> GEMPR_hit_tmax(tree, "Harm.PRPolGEMFarSide.hit.tmax");

//HCal and AA Hit Variables
TTreeReaderValue<int> HCal_hit_nhits(tree, "Harm.HCalScint.hit.nhits");
TTreeReaderValue<int> ActAna_hit_nhits(tree, "Harm.ActAnScint.hit.nhits");
TTreeReaderArray<double> ActAna_hit_xhit(tree, "Harm.ActAnScint.hit.xhitg");
TTreeReaderArray<double> ActAna_hit_yhit(tree, "Harm.ActAnScint.hit.yhitg");
TTreeReaderArray<double> ActAna_hit_zhit(tree, "Harm.ActAnScint.hit.zhitg");

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

 
TH1D *hResX = new TH1D("hResX","X residual;x_{fit}-x;counts",200,-2,2);
TH1D *hResY = new TH1D("hResY","Y residual;y_{fit}-y;counts",200,-2,2);

bool firstEvent = true;

// ===== Event Loop =====
int Nentries = tree.GetEntries();
int entry_num = 0;
while(tree.Next()){
    // =========================
    // GEM hit containers
    // =========================
    std::vector<double> x1, y1, z1;  // GEM1
    std::vector<double> x2, y2, z2;  // GEM2
    int N = *GEMPR_hit_nhits;
    if(N < 2) continue;
    // =========================
    // Split by plane
    // =========================
    for(int i=0;i<N;i++){
        if(GEMPR_hit_plane[i] == 1){
            x1.push_back(GEMPR_hit_xg[i]);
            y1.push_back(GEMPR_hit_yg[i]);
            z1.push_back(GEMPR_hit_zg[i]);
        }
        else if(GEMPR_hit_plane[i] == 2){
            x2.push_back(GEMPR_hit_xg[i]);
            y2.push_back(GEMPR_hit_yg[i]);
            z2.push_back(GEMPR_hit_zg[i]);
    }}

    // =========================
    // Require both planes
    // =========================
    if(x1.size() == 0 || x2.size() == 0) continue;
    // =========================
    // Merge GEM1 + GEM2 for LSQ fit
    // =========================
    std::vector<double> x, y, z;
    x.insert(x.end(), x1.begin(), x1.end());
    y.insert(y.end(), y1.begin(), y1.end());
    z.insert(z.end(), z1.begin(), z1.end());

    x.insert(x.end(), x2.begin(), x2.end());
    y.insert(y.end(), y2.begin(), y2.end());
    z.insert(z.end(), z2.begin(), z2.end());
    int n = x.size();
    if(n < 2) continue;
    // =========================
    // FIT X(Z)
    // =========================
    double Sz=0, Sx=0, Szz=0, Szx=0;
    for(int i=0;i<n;i++) {
        Sz  += z[i];
        Sx  += x[i];
        Szz += z[i]*z[i];
        Szx += z[i]*x[i];}

    double a =(n*Szx - Sz*Sx)/(n*Szz - Sz*Sz);
    double x0 =(Sx - a*Sz)/n;

    // =========================
    // FIT Y(Z)
    // =========================
    double Sy=0, Szy=0;
    for(int i=0;i<n;i++) {
        Sy  += y[i];
        Szy += z[i]*y[i];
    }
    double b = (n*Szy - Sz*Sy)/(n*Szz - Sz*Sz);
    double y0 =(Sy - b*Sz)/n;

    // =========================
    // RESIDUALS
    // =========================
    for(int i=0;i<n;i++){
        double xfit = x0 + a*z[i];
        double yfit = y0 + b*z[i];

        hResX->Fill(xfit - x[i]);
        hResY->Fill(yfit - y[i]);
    }

    // =========================
    // PLOT FIRST EVENT ONLY
    // =========================
    if(firstEvent){
        firstEvent = false;
        std::cout << "\nTRACK RESULT:\n";
        std::cout << "x(z) = " << x0 << " + " << a << " z\n";
        std::cout << "y(z) = " << y0 << " + " << b << " z\n";
        double zmin = *std::min_element(z.begin(), z.end());
        double zmax = *std::max_element(z.begin(), z.end());

        TCanvas *c1 = new TCanvas("c1","GEM Track",1000,500);
        c1->Divide(2,1);
        // =========================
        // X vs Z
        // =========================
        c1->cd(1);
        TGraph *gxz = new TGraph(n);

        for(int i=0;i<n;i++)
            gxz->SetPoint(i, z[i], x[i]);

        gxz->SetMarkerStyle(20);
        gxz->SetTitle("X vs Z;Z;X");
        gxz->Draw("AP");

        TF1 *fx = new TF1("fx","[0]+[1]*x", zmin, zmax);
        fx->SetParameters(x0, a);
        fx->SetLineWidth(2);
        fx->Draw("same");

        // =========================
        // Y vs Z
        // =========================
        c1->cd(2);
        TGraph *gyz = new TGraph(n);

        for(int i=0;i<n;i++)
            gyz->SetPoint(i, z[i], y[i]);

        gyz->SetMarkerStyle(20);
        gyz->SetTitle("Y vs Z;Z;Y");
        gyz->Draw("AP");

        TF1 *fy = new TF1("fy","[0]+[1]*x", zmin, zmax);
        fy->SetParameters(y0, b);
        fy->SetLineWidth(2);
        fy->Draw("same");
    }
}

// =========================
// Residual plots (AFTER LOOP)
// =========================
TCanvas *c2 = new TCanvas("c2","Residuals",1000,500);
c2->Divide(2,1);
c2->cd(1);
hResX->Draw();
c2->cd(2);
hResY->Draw();

} 
