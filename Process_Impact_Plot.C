// Stephen JD Kay - University of York
// Stephen.Kay@york.ac.uk
// 18/06/26
// Script to process EpIC generator output in combination with simulated output (efficiencies per x/Q2/t bin)

using namespace ROOT::Math;

#include "TString.h"
#include "Process_Sim_Output.h"
#include "ePICStyle.C"
#include <vector>

void Process_Impact_Plot.C(TString InGenFile="", TString InSimOutputFile=""){

  if(CheckFile(InGenFile) == kFALSE){ // Check files exist, can be opened and contain tree with fn
    exit(1);
  }
  if(CheckFile(InSimOutputFile) == kFALSE){ // Check files exist, can be opened and contain tree with fn
    exit(2);
  }
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);

  // Binning for Q2/x/t - From Oliver's script
  double q2edges[9]  = {1., 1.78, 3.16, 5.62, 10., 18., 32., 56., 100};
  double xBedges[12] = {1e-4, 2.5e-4, 6.3e-4, 1e-3, 2.5e-3, 6.3e-3, 1e-2, 2.5e-2, 6.3e-2, 0.1, 0.25, 0.7};
  double tedges[15]  = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.4, 1.6};
  // Count no. of bins (nEdges-1)
  const int nQ2bins = sizeof(q2edges)/sizeof(q2edges[0]) - 1;
  const int nxBbins = sizeof(xBedges)/sizeof(xBedges[0]) - 1;
  const int ntbins  = sizeof(tedges)/sizeof(tedges[0]) - 1;
  TH1D* tmpHist1D;
  TH2D* tmpHist2D;
  // Integrated lumi (in fb-1) of the generated file
  double IntLumiGen = 2*0.12; // 2* because the filelist is BOTH helicities

  // Open generator level file and get xB/Q2/t for each event, fill relevant histogram
  TH1D* h1_tGen_Q2xB[nQ2bins][nxBbins]; // Full t dists for each x/Q2 bin

  for(int binq2{0}; binq2<nQ2bins; binq2++){
    for(int binxB{0}; binxB<nxBbins; binxB++){
      	h1_tGen_Q2xB[binq2][binxB] = new TH1D(Form("h1_tGen_Q2xB[%i][%i]",binq2,binxB),
						Form("%.1f<Q^{2}<%.1f GeV^{2}, %.2e<x_{B}<%.2e;|t| [GeV^{2}];",
						     q2edges[binq2],q2edges[binq2+1],
						     xBedges[binxB],xBedges[binxB+1]),
						20, 0., 2.);
    }// End of xB bin loop
  } // End of Q2 bin loop

  TChain *AnalysisChain = new TChain("hepmc3_tree");
  AnalysisChain->Add(InGenFile);
  TTreeReader tree_reader(AnalysisChain);

  TTreeReaderArray<int>  PID(tree_reader, "hepmc3_event.particles.pid");
  TTreeReaderArray<int>  Status(tree_reader, "hepmc3_event.particles.status");
  TTreeReaderArray<double>  Px(tree_reader, "hepmc3_event.particles.momentum.m_v1");
  TTreeReaderArray<double>  Py(tree_reader, "hepmc3_event.particles.momentum.m_v2");
  TTreeReaderArray<double>  Pz(tree_reader, "hepmc3_event.particles.momentum.m_v3");
  TTreeReaderArray<double>  E(tree_reader, "hepmc3_event.particles.momentum.m_v4");
  
  while(tree_reader.Next()){ // Loop over all events
    // Reset booleans
    eBeam = kFALSE; pBeam = kFALSE; eSc = kFALSE; pSc = kFALSE; Gamma = kFALSE; GoodEvent = kFALSE;
    for(unsigned int i = 0; i < Status.GetSize(); i++){ // Loop over all particles in event, assign to 4-vectors
      if(Status[i] == 4 && PID[i] == 11) {Vec_e_beam.SetPxPyPzE(Px[i], Py[i], Pz[i], E[i]); eBeam = kTRUE;}
      if(Status[i] == 4 && PID[i] == 2212) {Vec_p_beam.SetPxPyPzE(Px[i], Py[i], Pz[i], E[i]); pBeam = kTRUE;}
      if(Status[i] == 1 && PID[i] == 11) {Vec_eSc.SetPxPyPzE(Px[i], Py[i], Pz[i], E[i]); eSc = kTRUE;}
      if(Status[i] == 1 && PID[i] == 2212) {Vec_pSc.SetPxPyPzE(Px[i], Py[i], Pz[i], E[i]); pSc = kTRUE;}
      if(Status[i] == 1 && PID[i] == 22) {Vec_gamma.SetPxPyPzE(Px[i], Py[i], Pz[i], E[i]); Gamma = kTRUE;}
    }
    if (eBeam == kTRUE && pBeam == kTRUE && eSc == kTRUE && pSc == kTRUE && Gamma == kTRUE){GoodEvent = kTRUE;}
    if (GoodEvent != kTRUE){cout << "Beam and output particles not found in event" << endl; continue;} // Check all relevant particles found and assigned

    // Need to undo AB for each particle
    // Calculate x/Q2/t from truth info
    Vec_Q2 = (Vec_e_beam - Vec_eSc); // Virtual photon beam vector
    Q2 = -1*(Vec_Q2.mag2());
    Vec_t = (Vec_Q2 - Vec_gamma );
    t = -1*(Vec_t.mag2());
    y =(Vec_p_beam.Dot(Vec_Q2))/(Vec_p_beam.Dot(Vec_e_beam));
    xB = Q2/(4*Vec_e_beam.E()*Vec_p_beam.E()*y);
    // Fill hists
    h1_Q2->Fill(Q2);
    h1_t->Fill(t);
    h1_xB->Fill(xB);
    h1_y->Fill(y);
    // Need to loop over binning scheme and fill relevant histograms - probably a more efficiency way of doing this, but it will work
    for(int binq2{0}; binq2<nQ2bins; binq2++){
      if(q2edges[binq2] < Q2 && Q2 < q2edges[binq2+1]){
	for(int binxB{0}; binxB<nxBbins; binxB++){
	  if(xBedges[binxB] < xB && xB < xBedges[binxB+1]){
	    h1_tDiff_v2[binq2][binxB]->Fill(t); // Fill t dist for this x/Q2 bin
	  }
	} // End xB binning loop
      }
    } // End Q2 binning loop
  } // End event loop
  
}
