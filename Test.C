void Test()
{
    TFile *f = TFile::Open("/w/halla-scshelf2102/sbs/saru/Simu_GEnRP/genrp_agc_sim4.root");
    TTree *T = nullptr;
    f->GetObject("T", T);

    Int_t nhFront, nhRear, nhGEM;

    std::vector<int> *front_pid = nullptr;
    std::vector<int> *rear_pid  = nullptr;
    std::vector<int> *gem_pid   = nullptr;

    std::vector<int> *front_trid = nullptr;
    std::vector<int> *rear_trid  = nullptr;
    std::vector<int> *gem_trid   = nullptr;

    T->SetBranchAddress("Harm.CEPolFront.hit.nhits", &nhFront);
    T->SetBranchAddress("Harm.CEPolRear.hit.nhits",  &nhRear);
    T->SetBranchAddress("Harm.PRPolGEMFarSide.hit.nhits", &nhGEM);

    T->SetBranchAddress("Harm.CEPolFront.hit.pid", &front_pid);
    T->SetBranchAddress("Harm.CEPolRear.hit.pid",  &rear_pid);
    T->SetBranchAddress("Harm.PRPolGEMFarSide.hit.pid", &gem_pid);

    T->SetBranchAddress("Harm.CEPolFront.hit.trid", &front_trid);
    T->SetBranchAddress("Harm.CEPolRear.hit.trid",  &rear_trid);
    T->SetBranchAddress("Harm.PRPolGEMFarSide.hit.trid", &gem_trid);

    Long64_t nentries = T->GetEntries();
	
int num_event =0;
    for (Long64_t i = 0; i < nentries; i++) {
        T->GetEntry(i);

        if (!front_pid || !rear_pid || !gem_pid ||
            !front_trid || !rear_trid || !gem_trid) continue;

        // =========================
        // (B) PID CONDITIONS
        // =========================
        bool frontHasProton = false;
        bool rearHasNeutron = false;
        bool gemHasProton   = false;

        for (int pid : *front_pid)
            if (pid == 2212 || pid ==2112) frontHasProton = true;

        for (int pid : *rear_pid)
            if (pid == 2212 || pid ==2112) rearHasNeutron = true;

        for (int pid : *gem_pid)
            if (pid == 2212) gemHasProton = true;

        if (!(frontHasProton && rearHasNeutron && gemHasProton))
            continue;

        // =========================
        // (A) TRID MATCH CONDITION
        // =========================
        bool foundCommonTrack = false;

        for (int trid_r : *rear_trid) {

            bool inFront = false;
            bool inGEM   = false;

            for (int trid_f : *front_trid) {
                if (trid_f == trid_r) {
                    inFront = true;
                    break;
                }
            }

            for (int trid_g : *gem_trid) {
                if (trid_g == trid_r) {
                    inGEM = true;
                    break;
                }
            }

            if (inFront && inGEM) {
                foundCommonTrack = true;
                break;
            }
        }

        // =========================
        // FINAL SELECTION
        // =========================
        if (foundCommonTrack) {

		num_event++;
            std::cout << "Entry " << i
                      << " passes PID + common-trid condition"
                      << std::endl;
        }
    }

    std::cout<< "Number of event:" <<num_event <<std::endl;

    f->Close();
}
