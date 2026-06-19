// Stephen JD Kay - University of York
// Stephen.Kay@york.ac.uk
// 18/06/26
// Script to process EpIC generator output in combination with simulated output (efficiencies per x/Q2/t bin)

using namespace ROOT::Math;
using ROOT::Math::VectorUtil::boost;

#include "TString.h"
#include "Process_Impact_Plot.h"
#include "ePICStyle.C"
#include <vector>

TH1::SetDefaultSumw2(kTRUE);

void Process_Impact_Plot(TString InGenFile="", TString InSimOutputFile=""){

  if(CheckFile_Gen(InGenFile) == kFALSE){ // Check files exist, can be opened and contain tree with fn
    exit(1);
  }
  if(CheckFile_Sim(InSimOutputFile) == kFALSE){ // Check files exist, can be opened and contain tree with fn
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
  TEfficiency *tmpEff;
  TLegend* Leg_Comp = new TLegend (0.1, 0.2, 0.9, 0.6);
  // Integrated lumi (in fb-1) of the generated file
  double IntLumiGen = 0.12; //  What matters is the int lumi in the processed HepMC3 file
  TDatime d;
  TFile *ofile = TFile::Open(Form("9x130_ImpactPlots_%d_0%d_%d.root", d.GetDay(), d.GetMonth(), d.GetYear()) ,"RECREATE");

  // Open generator level file and get xB/Q2/t for each event, fill relevant histogram
  TH1D* h1_tGen_Q2xB[nQ2bins][nxBbins]; // Full t dists for each x/Q2 bin
  TH1D* h1_tResult_Q2xB[nQ2bins][nxBbins]; // Full t dists for each x/Q2 bin - TEfficiency method
  TH1D* h1_tResult_Q2xB_v2[nQ2bins][nxBbins]; // Full t dists for each x/Q2 bin - TH1 method

  for(int binq2{0}; binq2<nQ2bins; binq2++){
    for(int binxB{0}; binxB<nxBbins; binxB++){
      	h1_tGen_Q2xB[binq2][binxB] = new TH1D(Form("h1_tGen_Q2xB[%i][%i]",binq2,binxB),
						Form("%.1f<Q^{2}<%.1f GeV^{2}, %.2e<x_{B}<%.2e;|t| [GeV^{2}];",
						     q2edges[binq2],q2edges[binq2+1],
						     xBedges[binxB],xBedges[binxB+1]),
						16, 0., 1.6);
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
    Vec_CoM_Boost = Vec_e_beam + Vec_p_beam;
    Boost_CoM.SetXYZ(-Vec_CoM_Boost.X()/Vec_CoM_Boost.E(), -Vec_CoM_Boost.Y()/Vec_CoM_Boost.E(), -Vec_CoM_Boost.Z()/Vec_CoM_Boost.E());
    Vec_e_beam = boost(Vec_e_beam, Boost_CoM);
    Vec_p_beam = boost(Vec_p_beam, Boost_CoM);
    RotX = RotationX(1.0*TMath::ATan2(Vec_p_beam.Y(), Vec_p_beam.Z()));
    RotY = RotationY(-1.0*TMath::ATan2(Vec_p_beam.X(), Vec_p_beam.Z()));
    Vec_p_beam = RotX*Vec_p_beam;
    Vec_p_beam = RotY*Vec_p_beam;
    Vec_e_beam = RotX*Vec_e_beam;
    Vec_e_beam = RotY*Vec_e_beam;
    Vec_HoF_Boost.SetPxPyPzE(0.,0., Vec_CoM_Boost.Z(), Vec_CoM_Boost.E());
    Boost_HoF.SetXYZ(Vec_HoF_Boost.X()/Vec_HoF_Boost.E(), Vec_HoF_Boost.Y()/Vec_HoF_Boost.E(), Vec_HoF_Boost.Z()/Vec_HoF_Boost.E());
    Vec_p_beam = boost(Vec_p_beam, Boost_HoF);
    Vec_e_beam = boost(Vec_e_beam, Boost_HoF);
    // Now for FS products, boost to CoM, rotate X and Y, boost to HoF
    Vec_eSc = boost(Vec_eSc, Boost_CoM);
    Vec_eSc = RotX*Vec_eSc;
    Vec_eSc = RotY*Vec_eSc;
    Vec_eSc = boost(Vec_eSc, Boost_HoF);
    Vec_pSc = boost(Vec_pSc, Boost_CoM);
    Vec_pSc = RotX*Vec_pSc;
    Vec_pSc = RotY*Vec_pSc;
    Vec_pSc = boost(Vec_pSc, Boost_HoF);
    Vec_gamma = boost(Vec_gamma, Boost_CoM);
    Vec_gamma = RotX*Vec_gamma;
    Vec_gamma = RotY*Vec_gamma;
    Vec_gamma = boost(Vec_gamma, Boost_HoF);
    
    // Calculate x/Q2/t from truth info
    Vec_Q2 = (Vec_e_beam - Vec_eSc); // Virtual photon beam vector
    Q2 = -1*(Vec_Q2.mag2());
    Vec_t = (Vec_Q2 - Vec_gamma );
    t = -1*(Vec_t.mag2());
    y =(Vec_p_beam.Dot(Vec_Q2))/(Vec_p_beam.Dot(Vec_e_beam));
    xB = Q2/(4*Vec_e_beam.E()*Vec_p_beam.E()*y);
    // Fill hists
    // h1_Q2->Fill(Q2);
    // h1_t->Fill(t);
    // h1_xB->Fill(xB);
    // h1_y->Fill(y);
    // Need to loop over binning scheme and fill relevant histograms - probably a more efficiency way of doing this, but it will work
    for(int binq2{0}; binq2<nQ2bins; binq2++){
      if(q2edges[binq2] < Q2 && Q2 < q2edges[binq2+1]){
	for(int binxB{0}; binxB<nxBbins; binxB++){
	  if(xBedges[binxB] < xB && xB < xBedges[binxB+1]){
	    h1_tGen_Q2xB[binq2][binxB]->Fill(t); // Fill t dist for this x/Q2 bin
	  }
	} // End xB binning loop
      }
    } // End Q2 binning loop
  } // End event loop

  // Scale by integrated lumi
  for(int binq2{0}; binq2<nQ2bins; binq2++){
    for(int binxB{0}; binxB<nxBbins; binxB++){
      h1_tResult_Q2xB[binq2][binxB] = (TH1D*)h1_tGen_Q2xB[binq2][binxB]->Clone(Form("h1_tResult_Q2xB[%i][%i]",binq2,binxB));
      h1_tResult_Q2xB[binq2][binxB]->Scale(1/IntLumiGen);
      h1_tResult_Q2xB_v2[binq2][binxB] = (TH1D*)h1_tGen_Q2xB[binq2][binxB]->Clone(Form("h1_tResult_Q2xB_v2[%i][%i]",binq2,binxB));
      h1_tResult_Q2xB_v2[binq2][binxB]->Scale(1/IntLumiGen);
      h1_tGen_Q2xB[binq2][binxB]->Scale(1/IntLumiGen); // For comparison later
    }
  }
  
  // Get bin by bin efficiencies and scale bins of histograms - TEfficiency method
  TFile *SimFile =  new TFile(InSimOutputFile);
  double tmp_content, tmp_error, tmp_EffVal, tmp_EffErr;
  for(int binq2{0}; binq2<nQ2bins; binq2++){
    for(int binxB{0}; binxB<nxBbins; binxB++){
      tmpEff = (TEfficiency*)(((TEfficiency*)SimFile->Get(Form("Q2xB_Binned_Dists/h1DEff_t_Q2xB_Eff[%i][%i]", binq2, binxB))));
      for(int i = 1; i <= h1_tResult_Q2xB[binq2][binxB]->GetNbinsX(); ++i){
	// Get initial values of bin and error
	tmp_content = h1_tResult_Q2xB[binq2][binxB]->GetBinContent(i);
	tmp_error = h1_tResult_Q2xB[binq2][binxB]->GetBinError(i);
	//Get corresponding efficiency value to scale by
	tmp_EffVal = tmpEff->GetEfficiency(i);
	// Take the average of the upper/lower efficiency error bars - symmetric error propagation
	tmp_EffErr =(tmpEff->GetEfficiencyErrorUp(i) + tmpEff->GetEfficiencyErrorLow(i))/2;
	// Set content and error
	h1_tResult_Q2xB[binq2][binxB]->SetBinContent(i, tmp_content * tmp_EffVal);
	h1_tResult_Q2xB[binq2][binxB]->SetBinError(i, TMath::Sqrt(TMath::Power(tmp_error*tmp_EffVal, 2) + TMath::Power(tmp_EffErr*tmp_content, 2)));
      }
    }
  }

  // Get bin by bin efficiencies and scale bins of histograms - TH1 method
  for(int binq2{0}; binq2<nQ2bins; binq2++){
    for(int binxB{0}; binxB<nxBbins; binxB++){
      tmpHist1D = (TH1D*)(((TH1D*)SimFile->Get(Form("Q2xB_Binned_Dists/h1_t_Q2xB_Eff[%i][%i]",binq2,binxB))));
      // Multiply hists together
      h1_tResult_Q2xB_v2[binq2][binxB]->Multiply(tmpHist1D); 
      // Manually multiply bin by bin
      // for(int i = 1; i <= h1_tResult_Q2xB_v2[binq2][binxB]->GetNbinsX(); ++i){
      // 	tmp_content = h1_tResult_Q2xB_v2[binq2][binxB]->GetBinContent(i);
      // 	tmp_error = h1_tResult_Q2xB_v2[binq2][binxB]->GetBinError(i);
      // 	tmp_EffVal = tmpHist1D->GetBinContent(i);
      // 	tmp_EffErr = tmpHist1D->GetBinError(i);
      // 	h1_tResult_Q2xB_v2[binq2][binxB]->SetBinContent(i, tmp_content * tmp_EffVal);
      // 	h1_tResult_Q2xB_v2[binq2][binxB]->SetBinError(i, TMath::Sqrt(TMath::Power(tmp_error*tmp_EffVal, 2) + TMath::Power(tmp_EffErr*tmp_content, 2))); 
      // }
    }
  }
  
  // TH1D* h1_tGen_Q2xB[nQ2bins][nxBbins]; // Full t dists for each x/Q2 bin
  // TH1D* h1_tResult_Q2xB[nQ2bins][nxBbins]; // Full t dists for each x/Q2 bin - TEfficiency method
  // TH1D* h1_tResult_Q2xB_v2[nQ2bins][nxBbins]; // Full t dists for each x/Q2 bin - TH1 method
  TCanvas* c_Q2xB_Results[9];
  TLatex *Q2_Range_Text[9];
  TString OutPdf = (Form("9x130_ImpactPlots_%d_0%d_%d_TEff_Ver.pdf", d.GetDay(), d.GetMonth(), d.GetYear()));
  for(int binq2{0}; binq2<nQ2bins; binq2++){
    c_Q2xB_Results[binq2] = new TCanvas(Form("c_Q2xB_Results_%i", binq2+1), Form("t Dists across xB bins, Q2 %i", binq2+1), 100, 0, 2560, 1920);
    c_Q2xB_Results[binq2]->Divide(4,3); 
    for(int binxB{0}; binxB<nxBbins; binxB++){
      c_Q2xB_Results[binq2]->cd(binxB+1);
      h1_tGen_Q2xB[binq2][binxB]->SetTitle(Form("%.2e<x_{B}<%.2e", xBedges[binxB],xBedges[binxB+1]));
      h1_tGen_Q2xB[binq2][binxB]->SetLineColor(kP6Blue);
      h1_tGen_Q2xB[binq2][binxB]->Draw("HISTERR");
      h1_tResult_Q2xB[binq2][binxB]->SetLineColor(kP6Yellow);
      h1_tResult_Q2xB[binq2][binxB]->Draw("SAMEHISTERR");
      if(binq2 == 0 && binxB == 0){
	Leg_Comp->AddEntry(h1_tGen_Q2xB[0][0], "Generated MC");
	Leg_Comp->AddEntry(h1_tResult_Q2xB[0][0], "ePIC Reco Scaled MC");
      }
      gPad->SetLogy(1);
    }
    c_Q2xB_Results[binq2]->cd(12);
    Leg_Comp->Draw();
    Q2_Range_Text[binq2] = new TLatex(0.2, 0.8, Form("%.1f<Q^{2}<%.1f GeV^{2}", q2edges[binq2],q2edges[binq2+1]));
    Q2_Range_Text[binq2]->Draw();
    if(binq2 == 0){
      c_Q2xB_Results[binq2]->Print(OutPdf + "(");
    }
    else if(binq2 == 7){
      c_Q2xB_Results[binq2]->Print(OutPdf + ")");
    }
    else{
      c_Q2xB_Results[binq2]->Print(OutPdf);
    }
  }

  TCanvas* c_Q2xB_Results_2[9];
  TString OutPdf2 = (Form("9x130_ImpactPlots_%d_0%d_%d_THist_Ver.pdf", d.GetDay(), d.GetMonth(), d.GetYear()));
  for(int binq2{0}; binq2<nQ2bins; binq2++){
    c_Q2xB_Results_2[binq2] = new TCanvas(Form("c_Q2xB_Results_2_%i", binq2+1), Form("t Dists across xB bins, Q2 %i", binq2+1), 100, 0, 2560, 1920);
    c_Q2xB_Results_2[binq2]->Divide(4,3); 
    for(int binxB{0}; binxB<nxBbins; binxB++){
      c_Q2xB_Results_2[binq2]->cd(binxB+1);
      h1_tGen_Q2xB[binq2][binxB]->SetTitle(Form("%.2e<x_{B}<%.2e", xBedges[binxB],xBedges[binxB+1]));
      h1_tGen_Q2xB[binq2][binxB]->SetLineColor(kP6Blue);
      h1_tGen_Q2xB[binq2][binxB]->Draw("HISTERR");
      h1_tResult_Q2xB[binq2][binxB]->SetLineColor(kP6Yellow);
      h1_tResult_Q2xB[binq2][binxB]->Draw("SAMEHISTERR");
      gPad->SetLogy(1);
    }
    c_Q2xB_Results_2[binq2]->cd(12);
    Leg_Comp->Draw();
    Q2_Range_Text[binq2]->Draw();
    if(binq2 == 0){
      c_Q2xB_Results_2[binq2]->Print(OutPdf2 + "(");
    }
    else if(binq2 == 7){
      c_Q2xB_Results_2[binq2]->Print(OutPdf2 + ")");
    }
    else{
      c_Q2xB_Results_2[binq2]->Print(OutPdf2);
    }
  }

  ofile->Write(); // Write histograms to file
  ofile->Close(); // Close output file
  
}
