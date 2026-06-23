#include <vector>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLine.h"

// ================= LSQ FIT =================
struct FitResult {
    double a;
    double b;
};

FitResult fitLSQ(const std::vector<double>& z,
                 const std::vector<double>& x)
{
    double S=0, Sz=0, Szz=0, Sx=0, Sxz=0;

    int n = z.size();
    for(int i=0;i<n;i++){
        S += 1;
        Sz += z[i];
        Szz += z[i]*z[i];
        Sx += x[i];
        Sxz += x[i]*z[i];
    }

    double denom = S*Szz - Sz*Sz;
    FitResult r;
    if(denom==0){ r.a=0; r.b=0; return r; }

    r.a = (S*Sxz - Sz*Sx)/denom;
    r.b = (Sx*Szz - Sz*Sxz)/denom;
    return r;
}

// ================= MAIN =================
void track_2D(const char* filename="/w/halla-scshelf2102/sbs/saru/Simu_GEnRP/genrp_agc_sim4.root")
{
    TFile *f = TFile::Open(filename);
    if(!f || f->IsZombie()){
        std::cout<<"File error"<<std::endl;
        return;
    }

    TTree *tree = (TTree*)f->Get("T");
    if(!tree){
        std::cout<<"Tree error"<<std::endl;
        return;
    }

    TTreeReader reader(tree);

    // ===== AA =====
    TTreeReaderValue<int> aa_n(reader,"Harm.ActAnScint.hit.nhits");
    TTreeReaderArray<double> aa_x(reader,"Harm.ActAnScint.hit.xhitg");
    TTreeReaderArray<double> aa_y(reader,"Harm.ActAnScint.hit.yhitg");
    TTreeReaderArray<double> aa_z(reader,"Harm.ActAnScint.hit.zhitg");

    // ===== HOD =====
    TTreeReaderValue<int> hod_n(reader,"Harm.PRPolScintFarSide.hit.nhits");
    TTreeReaderArray<double> hod_x(reader,"Harm.PRPolScintFarSide.hit.xhitg");
    TTreeReaderArray<double> hod_y(reader,"Harm.PRPolScintFarSide.hit.yhitg");
    TTreeReaderArray<double> hod_z(reader,"Harm.PRPolScintFarSide.hit.zhitg");

    // ===== GEM =====
    TTreeReaderValue<int> gem_n(reader,"Harm.PRPolGEMFarSide.hit.nhits");
    TTreeReaderArray<int> gem_plane(reader,"Harm.PRPolGEMFarSide.hit.plane");
    TTreeReaderArray<double> gem_x(reader,"Harm.PRPolGEMFarSide.hit.xg");
    TTreeReaderArray<double> gem_y(reader,"Harm.PRPolGEMFarSide.hit.yg");
    TTreeReaderArray<double> gem_z(reader,"Harm.PRPolGEMFarSide.hit.zg");

    // ================= HIST =================
TH1F *hAAx = new TH1F("hAAx","AA X residual;X_{AA} - X_{track} (m);Counts",200,-10,10);
TH1F *hAAy = new TH1F("hAAy","AA Y residual;Y_{AA} - Y_{track} (m);Counts",200,-10,10);
TH1F *hG0x = new TH1F("hG0x","GEM1 X residual;X_{GEM1} - X_{track} (m);Counts",200,-10,10);
TH1F *hG0y = new TH1F("hG0y","GEM1 Y residual;Y_{GEM1} - Y_{track} (m);Counts",200,-10,10);
TH1F *hG1x = new TH1F("hG1x","GEM2 X residual;X_{GEM2} - X_{track} (m);Counts",200,-10,10);
TH1F *hG1y = new TH1F("hG1y","GEM2 Y residual;Y_{GEM2} - Y_{track} (m);Counts",200,-10,10);
TH1F *hHOx = new TH1F("hHOx","HOD X residual;X_{HOD} - X_{track} (m);Counts",200,-10,10);
TH1F *hHOy = new TH1F("hHOy","HOD Y residual;Y_{HOD} - Y_{track} (m);Counts",200,-10,10);

TH1F *hChi2All = new TH1F("hChi2All","Global Chi2/NDF (All Detector Together)",100,-5,20);

    // FIRST EVENT STORAGE
    bool first=true;
    std::vector<double> fxz, fxx, fyz, fxy;
    // ================= LOOP =================
    while(reader.Next()) {
        // ---- AA ----
        if(*aa_n<=0) continue;
        double xAA=aa_x[0], yAA=aa_y[0], zAA=aa_z[0];

        // ---- HOD ----
        if(*hod_n<=0) continue;
        double xH=hod_x[0], yH=hod_y[0], zH=hod_z[0];

        // ---- GEM ----
        bool ok0=false, ok1=false;
        double xG0,yG0,zG0,xG1,yG1,zG1;
        for(int i=0;i<*gem_n;i++){
            if(gem_plane[i]==1 && !ok0){
                xG0=gem_x[i]; yG0=gem_y[i]; zG0=gem_z[i]; ok0=true;
            }
            if(gem_plane[i]==2 && !ok1){
                xG1=gem_x[i]; yG1=gem_y[i]; zG1=gem_z[i]; ok1=true;
            }
        }

        if(!ok0 || !ok1) continue;
        // ================= FIRST EVENT STORE =================
        if(first){
            first=false;
            fxz={zAA,zG0,zG1,zH};
            fxx={xAA,xG0,xG1,xH};
            fyz={zAA,zG0,zG1,zH};
            fxy={yAA,yG0,yG1,yH};
  }
{
    // ================= LOO FITS (ALL IN SAME SCOPE) =================

    std::vector<double> zAAfit = {zG0,zG1,zH};
    std::vector<double> xAAfit = {xG0,xG1,xH};
    auto fxAA = fitLSQ(zAAfit,xAAfit);

    std::vector<double> yAAfit = {yG0,yG1,yH};
    auto fyAA = fitLSQ(zAAfit,yAAfit);

    hAAx->Fill(xAA - (fxAA.a*zAA + fxAA.b));
    hAAy->Fill(yAA - (fyAA.a*zAA + fyAA.b));

    std::vector<double> zG0fit = {zAA,zG1,zH};
    std::vector<double> xG0fit = {xAA,xG1,xH};
    auto fxG0 = fitLSQ(zG0fit,xG0fit);

    std::vector<double> yG0fit = {yAA,yG1,yH};
    auto fyG0 = fitLSQ(zG0fit,yG0fit);

    hG0x->Fill(xG0 - (fxG0.a*zG0 + fxG0.b));
    hG0y->Fill(yG0 - (fyG0.a*zG0 + fyG0.b));

    std::vector<double> zG1fit = {zAA,zG0,zH};
    std::vector<double> xG1fit = {xAA,xG0,xH};
    auto fxG1 = fitLSQ(zG1fit,xG1fit);

    std::vector<double> yG1fit = {yAA,yG0,yH};
    auto fyG1 = fitLSQ(zG1fit,yG1fit);

    hG1x->Fill(xG1 - (fxG1.a*zG1 + fxG1.b));
    hG1y->Fill(yG1 - (fyG1.a*zG1 + fyG1.b));

    std::vector<double> zHfit = {zAA,zG0,zG1};
    std::vector<double> xHfit = {xAA,xG0,xG1};
    auto fxH = fitLSQ(zHfit,xHfit);

    std::vector<double> yHfit = {yAA,yG0,yG1};
    auto fyH = fitLSQ(zHfit,yHfit);

    hHOx->Fill(xH - (fxH.a*zH + fxH.b));
    hHOy->Fill(yH - (fyH.a*zH + fyH.b));

    // ================= GLOBAL CHI2 =================
    std::vector<double> zAll = {zAA, zG0, zG1, zH};
    std::vector<double> xAll = {xAA, xG0, xG1, xH};
    std::vector<double> yAll = {yAA, yG0, yG1, yH};

    double chi2 = 0;
    for(int i=0;i<4;i++){
        double pred_x = fxAA.a*zAll[i] + fxAA.b; // using AA-excluded track
        double pred_y = fyAA.a*zAll[i] + fyAA.b;

        chi2 += ((xAll[i] - pred_x)*(xAll[i] - pred_x)) + ((yAll[i] - pred_y)*(yAll[i] - pred_y)) ;
//        chi2 += (yAll[i] - pred_y)*(yAll[i] - pred_y);
    }
    double chi2_ndf = chi2 / 4.0;
    hChi2All->Fill(chi2_ndf);
}

}


    // ================= FIRST EVENT PLOTS =================
    TCanvas *c1=new TCanvas("c1","First XZ",800,600);
    TGraph *g1=new TGraph(fxz.size(),fxz.data(),fxx.data());
    g1->SetTitle("XZ Track; Z position; X position");  // <-- axis labels
    g1->SetMarkerStyle(29);
    g1->SetMarkerColor(kRed);   // <-- red points
    g1->SetMarkerSize(2);
    g1->Draw("AP");

    TCanvas *c2=new TCanvas("c2","First YZ",800,600);
    TGraph *g2=new TGraph(fyz.size(),fyz.data(),fxy.data());
    g2->SetTitle("YZ Track; Z position; Y position");  // <-- axis labels
    g2->SetMarkerStyle(29);
    g2->SetMarkerColor(kRed);   // <-- red points
    g2->SetMarkerSize(2);
    g2->Draw("AP");

    // ================= RESIDUAL PLOTS =================
    TCanvas *c3=new TCanvas("c3","X res",1200,800);
    c3->Divide(2,2);
    c3->cd(1); hAAx->Draw();
    c3->cd(2); hG0x->Draw();
    c3->cd(3); hG1x->Draw();
    c3->cd(4); hHOx->Draw();

    TCanvas *c4=new TCanvas("c4","Y res",1200,800);
    c4->Divide(2,2);
    c4->cd(1); hAAy->Draw();
    c4->cd(2); hG0y->Draw();
    c4->cd(3); hG1y->Draw();
    c4->cd(4); hHOy->Draw();

   
TCanvas *cChi2 = new TCanvas("cChi2","Global Chi2/NDF",800,600);
//cChi2->SetLogy();
hChi2All->SetLineWidth(2);
hChi2All->GetXaxis()->SetTitle("#chi^{2}/NDF");
hChi2All->GetYaxis()->SetTitle("Counts");

hChi2All->Draw();

}

