void hit_multi()
{
    TFile *f = TFile::Open("/w/halla-scshelf2102/sbs/saru/Simu_GEnRP/genrp_agc_sim4.root");
    TTree *T;
    f->GetObject("T", T);

    Long64_t nEvents = T->GetEntries();

    // hit counter
    Int_t nhFront, nhRear,nhGEM;

    std::vector<int> *front_pid = nullptr;
    std::vector<int> *rear_pid  = nullptr;
    std::vector<int> *gem_pid   = nullptr;


    T->SetBranchAddress("Harm.CEPolFront.hit.nhits", &nhFront);
    T->SetBranchAddress("Harm.CEPolRear.hit.nhits",  &nhRear);
    T->SetBranchAddress("Harm.PRPolGEMFarSide.hit.nhits", &nhGEM);

    T->SetBranchAddress("Harm.CEPolFront.hit.pid", &front_pid);
    T->SetBranchAddress("Harm.CEPolRear.hit.pid",  &rear_pid);
    T->SetBranchAddress("Harm.PRPolGEMFarSide.hit.pid", &gem_pid);


    // histograms
    TH1F *hFront = new TH1F("hFront","Front Hit Multiplicity;Hits;Events",80,0,80);
    TH1F *hRear  = new TH1F("hRear","Rear Hit Multiplicity;Hits;Events",80,0,80);
    TH1F *hGEM   = new TH1F("hGEM","PR Hit Multiplicity;Hits;Events",80,0,80);
    
for (Long64_t i = 0; i < nEvents; i++) {
    T->GetEntry(i);
    /*
    for (int j = 0; j < nhFront; j++) {
        if ((*front_pid)[j] == 2212)
            hFront->Fill(nhFront);
    }

    for (int j = 0; j < nhRear; j++) {
        if ((*rear_pid)[j] == 2212)
            hRear->Fill(nhRear);
    }

    for (int j = 0; j < nhGEM; j++) {
        if ((*gem_pid)[j] ==2212)
            hGEM->Fill(nhGEM);
    }
}
*/
for (int j = 0; j < nhFront; j++) {
    if ((*front_pid)[j] != 2212 && (*front_pid)[j] != 2112)
        hFront->Fill(nhFront);
}

for (int j = 0; j < nhRear; j++) {
    if ((*rear_pid)[j] != 2212 && (*rear_pid)[j] != 2112)
        hRear->Fill(nhRear);
}

for (int j = 0; j < nhGEM; j++) {
    if ((*gem_pid)[j] != 2212 && (*gem_pid)[j] != 2112)
        hGEM->Fill(nhGEM);
}
}

// Set line colors
hFront->SetLineColor(kRed);
hRear->SetLineColor(kBlue);
hGEM->SetLineColor(kGreen);

// Make lines thicker (optional)
hFront->SetLineWidth(2);
hRear->SetLineWidth(2);
hGEM->SetLineWidth(2);

// Turn off statistics box (optional)
gStyle->SetOptStat(0);

// Axis titles
hFront->SetTitle("Other particle Hit Multiplicity");
hFront->GetXaxis()->SetTitle("Hit Multiplicity");
hFront->GetYaxis()->SetTitle("Counts");

// Create one canvas
TCanvas *c1 = new TCanvas("c1","other particle Hit Multiplicity",800,600);

// Set y-axis maximum so all histograms are visible
double ymax = std::max({hFront->GetMaximum(),
                        hRear->GetMaximum(),
                        hGEM->GetMaximum()});

hFront->SetMaximum(1.1*ymax);

// Draw
hFront->Draw("HIST");
hRear->Draw("HIST SAME");
hGEM->Draw("HIST SAME");

// Legend
TLegend *leg = new TLegend(0.65,0.70,0.88,0.88);
leg->AddEntry(hFront,"CEPol Front","l");
leg->AddEntry(hRear,"CEPol Rear","l");
leg->AddEntry(hGEM,"PRPol GEM","l");
leg->Draw();

c1->Update();





}
