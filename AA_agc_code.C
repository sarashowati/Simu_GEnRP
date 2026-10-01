// truth_scenarios_example.C
// -----------------------------------------------------------------------------
//   PTrack  = one row per primary generator particle (usually the scattered
//             electron AND the recoil nucleon), index is inconsistent.
//             Find the recoil by PID (2112=n, 2212=p).
//   OTrack  = the particle's species as of its most recent "rebirth". 
//             Only gets a new OTrack row when a new track is created in the target or analyzer volume.
//             So OTrack.PID vs PTrack.PID lets you compare the particle's identity at different stages, 
//             and OTrack.posz for position info.
//   SDTrack = snapshot at the boundary crossing for a given hit.
//
// Run:
//   root -l -b -q 'truth_scenarios_example.C+("your_sim_file.root")'
// -----------------------------------------------------------------------------

#include <TChain.h>
#include <TLeaf.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <vector>
#include <map>
#include <set>
#include <cstdio>
#include <cmath>

// ---- geometry: to define which region is "target" or "analyzer"
// (OTrack.posz is hall-frame, so rotate into the transport frame first), and hcal distances for dx
const double thSBS_deg = 24.7;
const double Eb = 4.3, Mtgt = 0.93827;          // beam energy, nucleon mass [GeV]
const double hcalDist = 5.6534, hcal_voffset = 0.75;   // HCal distance + real vertical offset [m]
const double zTgtMax = 0.5;                     // |transport z| < this -> born at the target
const double zAnaLo = 4.30, zAnaHi = 4.65;       // transport z window -> born in the analyzer
const double hcalHitThreshold = 0.007;           // GeV, per-block threshold for clustering

// ---- 7 scenarios ----
enum Scenario { SC_NN, SC_PP, SC_NP, SC_PN, SC_RESCAT_N, SC_RESCAT_P, SC_OTHER, SC_NA };
const char* scenName[8] = {"nn","pp","np","pn","rescat_n","rescat_p","other","n/a"};
const int   scenColor[8] = {kBlack, kRed+1, kBlue+1, kAzure+7, kGreen+2, kOrange+7, kGray+2, kMagenta};

// -- classifier: given OTrack/PTrack indices for one hit or track, decide what happened.
int classifyOrigin(int otIdx, int ptIdx,
                    const std::vector<int>& otPID, const std::vector<double>& otPosZ_transport,
                    const std::vector<int>& ptPID)
{
  if(otIdx<0 || ptIdx<0 || otIdx>=(int)otPID.size() || ptIdx>=(int)ptPID.size()) return SC_NA; // invalid indices

  const int otSpecies = otPID[otIdx];
  const double otZ    = otPosZ_transport[otIdx];
  const int ptSpecies = ptPID[ptIdx];

  const bool bornAtTarget   = std::fabs(otZ) < zTgtMax;
  const bool bornInAnalyzer = (otZ > zAnaLo && otZ < zAnaHi);

  if(bornAtTarget){
    // in target and still the same species it was generated as -- no charge exchange
    if(otSpecies==2112 && ptSpecies==2112) return SC_NN;
    if(otSpecies==2212 && ptSpecies==2212) return SC_PP;
    return SC_OTHER;
  }
  if(bornInAnalyzer){
    if(otSpecies==2112 && ptSpecies==2212) return SC_NP;   // p -> n in the analyzer
    if(otSpecies==2212 && ptSpecies==2112) return SC_PN;   // n -> p in the analyzer
    if(otSpecies==ptSpecies && (otSpecies==2112||otSpecies==2212))
      return (otSpecies==2112) ? SC_RESCAT_N : SC_RESCAT_P;   // same species, but re-scattered
    return SC_OTHER;   // analyzer-born, but not a clean n/p 
  }
  return SC_OTHER;   // born somewhere else
}

// ---- minimal HCal seeded-island clustering: group neighbouring fired blocks,
// keep the highest-energy island, and remember which raw hit was its seed (that
// seed hit's own otridx/ptridx is what we classify with).
struct HCalCluster { bool found=false; double x=0, y=0, edep=0; int seedRawHit=-1; };

HCalCluster clusterHCal(int nHits, const std::vector<double>& edep,
                         const std::vector<double>& xcell, const std::vector<double>& ycell,
                         const std::vector<int>& row, const std::vector<int>& col, const std::vector<int>& cell)
{
  std::map<int,int> cellToIndex;
  std::vector<double> cellEdep, cellX, cellY; std::vector<int> cellRow, cellCol, cellRawHit;
  for(int j=0;j<nHits;++j){
    if(edep[j]<=0) continue;
    auto it = cellToIndex.find(cell[j]);
    if(it==cellToIndex.end()){
      cellToIndex[cell[j]] = (int)cellEdep.size();
      cellEdep.push_back(edep[j]); cellX.push_back(xcell[j]); cellY.push_back(ycell[j]);
      cellRow.push_back(row[j]); cellCol.push_back(col[j]); cellRawHit.push_back(j);
    } else {
      const int k = it->second;
      if(edep[j] > edep[cellRawHit[k]]) cellRawHit[k] = j;
      cellEdep[k] += edep[j];
    }
  }
  std::set<int> unused;
  for(size_t k=0;k<cellEdep.size();++k) if(cellEdep[k]>=hcalHitThreshold) unused.insert((int)k);

  HCalCluster best;
  while(!unused.empty()){
    int seed=-1; double seedE=-1;
    for(int k : unused) if(cellEdep[k]>seedE){ seedE=cellEdep[k]; seed=k; }
    std::vector<int> island = {seed}; unused.erase(seed);
    double xsum=cellEdep[seed]*cellX[seed], ysum=cellEdep[seed]*cellY[seed], esum=cellEdep[seed];
    int seedRawHit = cellRawHit[seed];
    for(size_t idx=0; idx<island.size(); ++idx){
      const int r=cellRow[island[idx]], c=cellCol[island[idx]];
      const int dr[4]={0,0,-1,1}, dc[4]={-1,1,0,0};
      for(int d=0; d<4; ++d){
        const int rr=r+dr[d], cc=c+dc[d];
        if(rr<0||rr>=8||cc<0||cc>=4) continue;
        auto it = cellToIndex.find(cc + 4*rr); if(it==cellToIndex.end()) continue;
        const int k = it->second; if(unused.find(k)==unused.end()) continue;
        unused.erase(k); island.push_back(k);
        xsum += cellEdep[k]*cellX[k]; ysum += cellEdep[k]*cellY[k]; esum += cellEdep[k];
      }
    }
    if(esum > best.edep){ best.found=true; best.x=xsum/esum; best.y=ysum/esum; best.edep=esum; best.seedRawHit=seedRawHit; }
  }
  return best;
}

void AA_agc_code(const char* file)
{
  gStyle->SetOptStat(0);
  TChain chain("T");
  chain.Add(file);

  const double thSBS_rad = thSBS_deg*TMath::DegToRad();
  const TVector3 zHatTr(-std::sin(thSBS_rad),0,std::cos(thSBS_rad)), xHatTr(0,-1,0);
  const TVector3 yHatTr = zHatTr.Cross(xHatTr).Unit();
  const TLorentzVector beamP4(0,0,Eb,Eb), targetP4(0,0,0,Mtgt);
  auto toTransport = [&](const TVector3& p){ return TVector3(p.Dot(xHatTr),p.Dot(yHatTr),p.Dot(zHatTr)); }; // transform to transport frame
  
  // predicted positions to a plane in SBS from straight line to plane at z=planeZ, starting from vertexTr and recoilP4.
  auto projectToPlane = [&](TLorentzVector r, const TVector3& vTr, double planeZ, double& oX, double& oY){
    r.RotateY(thSBS_rad);
    oX = -(planeZ - vTr.Z())*r.Py()/r.Pz() + vTr.X();
    oY =  (planeZ - vTr.Z())*r.Px()/r.Pz() + vTr.Y();
  };

  // HCal hits, for clustering + the seed hit's own otridx/ptridx
  int hcNH=0;
  std::vector<double> *hcE=0,*hcX=0,*hcY=0; std::vector<int> *hcRow=0,*hcCol=0,*hcCell=0;
  std::vector<int> *hcOTIdx=0,*hcPTIdx=0;
  chain.SetBranchAddress("Harm.ActAnScint.hit.nhits",&hcNH);
  chain.SetBranchAddress("Harm.ActAnScint.hit.sumedep",&hcE);
  chain.SetBranchAddress("Harm.ActAnScint.hit.xcell",&hcX);
  chain.SetBranchAddress("Harm.ActAnScint.hit.ycell",&hcY);
  chain.SetBranchAddress("Harm.ActAnScint.hit.row",&hcRow);
  chain.SetBranchAddress("Harm.ActAnScint.hit.col",&hcCol);
  chain.SetBranchAddress("Harm.ActAnScint.hit.cell",&hcCell);
  chain.SetBranchAddress("Harm.ActAnScint.hit.otridx",&hcOTIdx);
  chain.SetBranchAddress("Harm.ActAnScint.hit.ptridx",&hcPTIdx);

  // CEPolRear = the analyzer's rear tracker (charged particles only). Need P too to pick the highest-momentum track.
  int nAnaTracks=0;
  std::vector<double> *anaP=0;
  std::vector<int> *anaOTIdx=0, *anaPTIdx=0;
  chain.SetBranchAddress("Harm.CEPolRear.Track.ntracks", &nAnaTracks);
  chain.SetBranchAddress("Harm.CEPolRear.Track.P", &anaP);
  chain.SetBranchAddress("Harm.CEPolRear.Track.otridx", &anaOTIdx);
  chain.SetBranchAddress("Harm.CEPolRear.Track.ptridx", &anaPTIdx);

  std::vector<int> *otPID=0;
  std::vector<double> *otPosX=0, *otPosY=0, *otPosZ=0;
  chain.SetBranchAddress("OTrack.PID", &otPID);
  chain.SetBranchAddress("OTrack.posx", &otPosX);
  chain.SetBranchAddress("OTrack.posy", &otPosY);
  chain.SetBranchAddress("OTrack.posz", &otPosZ);

  std::vector<int> *ptPID=0;
  chain.SetBranchAddress("PTrack.PID", &ptPID);

  // ev.rate is the per-event physics weight (sigma*Lumi*GenVol/Nthrown, in Hz) --
  // MC is generated flat over a wide kinematic range, not proportional to cross section, 
  // need to normalize by the event weight for a true representation.
  chain.LoadTree(0);
  TLeaf *leafTh=chain.GetLeaf("ev","th"), *leafPh=chain.GetLeaf("ev","ph");
  TLeaf *leafVx=chain.GetLeaf("ev","vx"), *leafVy=chain.GetLeaf("ev","vy"), *leafVz=chain.GetLeaf("ev","vz");
  TLeaf *leafRate=chain.GetLeaf("ev","rate");
  int curTree = chain.GetTreeNumber(); // track when the TChain switches to a new TTree, so we can re-load the leaves if needed

  // ---- histograms and counters
  double countHCal[8] = {0}, countAna[8] = {0};
  TH1D* hDx[7];
  TH1D* hDx_all = new TH1D("hDx_all","",80,-3.0,3.0); hDx_all->SetLineColor(kBlack); hDx_all->SetLineWidth(2);
  for(int s=0;s<7;++s){ hDx[s]=new TH1D(Form("hDx_%s",scenName[s]),"",80,-3.0,3.0); hDx[s]->SetLineColor(scenColor[s]); hDx[s]->SetLineWidth(2); }

  // ---- loop over events
  Long64_t n = chain.GetEntries();
  for(Long64_t i=0; i<n; i++){

    if (i% 10000 == 0) printf("\rProcessing event %lld / %lld", i, n), fflush(stdout);

    chain.GetEntry(i);
    if(chain.GetTreeNumber()!=curTree){ curTree=chain.GetTreeNumber(); // rebind the leaves if the TChain switched to a new TTree (will happen if the input file is a TChain of multiple files)
      leafTh=chain.GetLeaf("ev","th"); leafPh=chain.GetLeaf("ev","ph");
      leafVx=chain.GetLeaf("ev","vx"); leafVy=chain.GetLeaf("ev","vy"); leafVz=chain.GetLeaf("ev","vz");
      leafRate=chain.GetLeaf("ev","rate"); }
    const double w = leafRate->GetValue();

    // OTrack.posz is hall-frame -- rotate every row into transport-z once per event
    std::vector<double> otPosZTr(otPosX->size());
    for(size_t j=0; j<otPosX->size(); ++j)
      otPosZTr[j] = toTransport(TVector3(otPosX->at(j),otPosY->at(j),otPosZ->at(j))).Z();

    // ---- classify via the analyzer track (charged only); highest-momentum track
    if(nAnaTracks>0 && anaP && !anaP->empty()){
      int iBest=0; double pBest=-1;
      for(int k=0;k<(int)anaP->size();++k) if(anaP->at(k)>pBest){ iBest=k; pBest=anaP->at(k); }
      int scenario = classifyOrigin(anaOTIdx->at(iBest), anaPTIdx->at(iBest), *otPID, otPosZTr, *ptPID);
      countAna[scenario] += w;
    }

    // ---- classify via the HCal cluster; the seed hit's own otridx/ptridx is what we classify with (sees ALL species, not just charged)
    HCalCluster clus = clusterHCal(hcNH, *hcE,*hcX,*hcY,*hcRow,*hcCol,*hcCell);
    if(!clus.found) continue;
    int j = clus.seedRawHit;
    int scenario = classifyOrigin(hcOTIdx->at(j), hcPTIdx->at(j), *otPID, otPosZTr, *ptPID);
    countHCal[scenario] += w;

    // dx_HCal = actual - predicted 
    const double th=leafTh->GetValue(), ph=leafPh->GetValue();
    TVector3 eDir(std::sin(th)*std::cos(ph), std::sin(th)*std::sin(ph), std::cos(th));
    const double elecE = Eb/(1.0 + (Eb/Mtgt)*(1.0-std::cos(th)));
    TLorentzVector eP4; eP4.SetVectM(elecE*eDir,0.0);
    TLorentzVector rP4 = targetP4 + (beamP4 - eP4);
    TVector3 vTr = toTransport(TVector3(leafVx->GetValue(),leafVy->GetValue(),leafVz->GetValue()));
    double predX,predY; projectToPlane(rP4,vTr,hcalDist,predX,predY);
    const double measX = -(clus.y);
    const double dxHCal = measX - predX;
    hDx_all->Fill(dxHCal, w);
    if(scenario!=SC_NA) hDx[scenario]->Fill(dxHCal, w); // fill the corresponding scenario histogram for the dx_HCal distribution
  }

  printf("\n");
  printf("scenario rates [Hz], classified via the AA SDTrack:\n");
  for(int s=0; s<8; s++) printf("  %-10s %.4g\n", scenName[s], countHCal[s]);
  printf("\nscenario rates [Hz], classified via the analyzer GEM SDTrack --\n");
  for(int s=0; s<8; s++) printf("  %-10s %.4g\n", scenName[s], countAna[s]);

  // dx_HCal plot, one curve per scenario, to see the classification in action
  TCanvas* c = new TCanvas("c","dx_Active_Ana by scenario",800,600);
  double mx=0; for(int s=0;s<7;++s) mx=std::max(mx,hDx[s]->GetMaximum());
  auto* lg = new TLegend(0.7,0.55,0.9,0.88);
  lg->AddEntry(hDx_all,"All scenarios","l");
  hDx_all->SetLineColor(kGray+2); hDx_all->SetLineWidth(2); hDx_all->SetLineStyle(2);
  hDx_all->SetMaximum(1.3*mx); hDx_all->SetTitle("#Deltax_{AA} by scenario;#Deltax_{AA} [m];events");
  hDx_all->Draw(); 
  for(int s=0;s<7;++s){
    hDx[s]->SetMaximum(1.3*mx);
    hDx[s]->SetTitle("#Deltax_{AA} by scenario;#Deltax_{AA} [m];events");
    hDx[s]->Draw(s==0?"hist":"hist same");
    lg->AddEntry(hDx[s],scenName[s],"l");
  }
  lg->Draw();
//  c->SaveAs("truth_scenarios_example_dx.png");
}
