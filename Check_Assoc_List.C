// Stephen JD Kay - University of York
// Stephen.Kay@york.ac.uk
// 12/06/26
// Script to process reco output and plot theta for Reconstructed Truth Seeded Charged Particles

using namespace ROOT::Math;

#include "TString.h"
#include "Check_Assoc_List.h"
#include "ePICStyle.C"
#include <vector>

void Check_Assoc_List(TString InFile = ""){

  if(CheckFile(InFile) == kFALSE){ // Check files exist, can be opened and contain tree with fn
    exit(1);
  }

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);

  TChain *AnalysisChain = new TChain("events");
  ifstream fstream (InFile);
  string line;
  int FileCount = 0;
  while(getline(fstream, line)){
    if ( FileCount != 0 && FileCount % 100 == 0 ) {
      cout << "Added " << FileCount << " files to chain" << endl;
    }
    //if (FileCount > 5) continue; // Stop loop after 5 files, comment out to read full file
    TString tmpFile{line};
    auto RootFile = TFile::Open(tmpFile);
    if(CheckRootFile(tmpFile) == kFALSE) cout << "!!!!! File either missing or broken !!!!!" << endl << tmpFile << endl << "!!!!! File either missing or broken !!!!!" << endl;
    else AnalysisChain->Add(tmpFile);
    FileCount++;
  }
  cout << "All files in file list - " << InFile << " - added to chain" << endl;

  TTreeReader tree_reader(AnalysisChain);
  TFile *ofile = TFile::Open(Form("%s_OutputHists.root", ((InFile.Remove(InFile.Length() - 5)).Remove(0, 24)).Data()),"RECREATE"); // 22 cuts the .eicrecon.edm4eic.root from the original input file name

  XYZVector Vec_pSc_MC;
  XYZVector Vec_pSc_Rec;
  XYZVector Vec_pSc_Rec_Rot;
  XYZVector Vec_pSc_RecAssoc;
  XYZVector Vec_pSc_RecAssoc_Rot;
  //XYZVector Vec_tmpRP;
  
  TLegend* Leg_Comp = new TLegend (0.1, 0.2, 0.9, 0.6);
  TLegend* Leg_Comp2 = new TLegend (0.1, 0.2, 0.9, 0.6);
  TLegend* Leg_Comp3 = new TLegend (0.1, 0.2, 0.9, 0.6);
  
  TH1D* h1_Theta_MC = new TH1D("h1_Theta_MC", "#theta_{MC}; #theta (mrad); Counts/0.2 mrad bin", 125, 0, 25);
  TH1D* h1_Theta_MC_RP = new TH1D("h1_Theta_MC_RP", "#theta_{MC} in RP; #theta (mrad); Counts/0.2 mrad bin", 125, 0, 25);
  TH1D* h1_Theta_MC_B0 = new TH1D("h1_Theta_MC_B0", "#theta_{MC} in B0; #theta (mrad); Counts/0.2 mrad bin", 125, 0, 25);
  TH1D* h1_Theta_RecAssoc = new TH1D("h1_Theta_RecAssoc", "#theta_{RecAssoc}; #theta (mrad); Counts/0.2 mrad bin", 125, 0, 25);
  TH1D* h1_Theta_Rec = new TH1D("h1_Theta_Rec", "#theta_{Rec}; #theta (mrad); Counts/0.2 mrad bin", 125, 0, 25);
  TH1D* h1_Theta_Rec_RP = new TH1D("h1_Theta_Rec_RP", "#theta_{Rec_RP}; #theta (mrad); Counts/0.2 mrad bin", 125, 0, 25);

  TH1D* h1_eta_MC = new TH1D("h1_eta_MC", "#eta_{MC}; #eta ; Counts/0.1 bin", 60, 4, 10);
  TH1D* h1_eta_MC_RP = new TH1D("h1_eta_MC_RP", "#eta_{MC} in RP; #eta ; Counts/0.1 bin", 60, 4, 10);
  TH1D* h1_eta_MC_B0 = new TH1D("h1_eta_MC_B0", "#eta_{MC} in B0; #eta ; Counts/0.1 bin", 60, 4, 10);
  TH1D* h1_eta_RecAssoc = new TH1D("h1_eta_RecAssoc", "#eta_{RecAssoc}; #eta ; Counts/0.1 bin", 60, 4, 10);
  TH1D* h1_eta_Rec = new TH1D("h1_eta_Rec", "#eta_{Rec}; #eta ; Counts/0.1 bin", 60, 4, 10);
  TH1D* h1_eta_Rec_RP = new TH1D("h1_eta_Rec_RP", "#eta_{Rec_RP}; #eta ; Counts/0.1 bin", 60, 4, 10);

  // Series of hists for efficiency plots
  TH1D* h1_Theta_pScRecMC = new TH1D("h1_Theta_MC_pScRec", "#theta_{MC}, p' Rec; #theta (mrad); Counts/0.2 mrad bin", 125, 0, 25);
  TH1D* h1_Theta_MC_pScRecRP = new TH1D("h1_Theta_MC_pScRecRP", "#theta_{MC}, p' in RP; #theta (mrad); Counts/0.2 mrad bin", 125, 0, 25);
  TH1D* h1_Theta_MC_pScRecB0 = new TH1D("h1_Theta_MC_pScRecB0", "#theta_{MC}, p' in B0; #theta (mrad); Counts/0.2 mrad bin", 125, 0, 25);
  TH1D* h1_eta_MC_pScRec = new TH1D("h1_eta_MC_pScRec", "#eta_{MC}, p' Rec; #eta ; Counts/0.1 bin", 60, 4, 10);
  TH1D* h1_eta_MC_pScRecRP = new TH1D("h1_eta_MC_pScRecRP", "#eta_{MC}, p' in RP; #eta ; Counts/0.1 bin", 60, 4, 10);
  TH1D* h1_eta_MC_pScRecB0 = new TH1D("h1_eta_MC_pScRecB0", "#eta_{MC}, p' in B0; #eta ; Counts/0.1 bin", 60, 4, 10);
  
  TTreeReaderArray<int>  MC_PID(tree_reader, "MCParticlesHeadOnFrameNoBeamFX.PDG");
  TTreeReaderArray<int>  MC_Status(tree_reader, "MCParticlesHeadOnFrameNoBeamFX.generatorStatus");
  TTreeReaderArray<double>  MC_Px(tree_reader, "MCParticlesHeadOnFrameNoBeamFX.momentum.x");
  TTreeReaderArray<double>  MC_Py(tree_reader, "MCParticlesHeadOnFrameNoBeamFX.momentum.y");
  TTreeReaderArray<double>  MC_Pz(tree_reader, "MCParticlesHeadOnFrameNoBeamFX.momentum.z");
  TTreeReaderArray<float>  TSCP_Px(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x");
  TTreeReaderArray<float>  TSCP_Py(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y");
  TTreeReaderArray<float>  TSCP_Pz(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z");
  TTreeReaderArray<int> TSCPRec_Assoc(tree_reader, "_ReconstructedTruthSeededChargedParticleAssociations_rec.index");
  TTreeReaderArray<int> TSCPSim_Assoc(tree_reader, "_ReconstructedTruthSeededChargedParticleAssociations_sim.index");
  TTreeReaderArray<float>  RP_Px(tree_reader, "ForwardRomanPotRecParticles.momentum.x");
  TTreeReaderArray<float>  RP_Py(tree_reader, "ForwardRomanPotRecParticles.momentum.y");
  TTreeReaderArray<float>  RP_Pz(tree_reader, "ForwardRomanPotRecParticles.momentum.z");

  cout << "Processing - " << nEntries << " events" << endl;
  
  while(tree_reader.Next()){ // Loop over all events
    EventCounter++;
    if ( EventCounter % ( nEntries / 20 ) == 0 ) {
      cout << "Processed " << setw(4) << ceil(((1.0*EventCounter)/(1.0*nEntries))*100.0) << " % of events - " << EventCounter << endl;
    }
    // Reset index and reset vector for each event
    int pSc_Index = 0;
    Bool_t B0_Part=kFALSE;
    // Find MC info
    for (unsigned int i = 0; i < MC_Status.GetSize(); i++){
      if(MC_Status[i] == 1 && MC_PID[i] == 2212){Vec_pSc_MC.SetXYZ(MC_Px[i],MC_Px[i],MC_Pz[i]); pSc_Index = i;}
    } // End MC loop'
    if(Vec_pSc_MC.Theta()*1000 < 0.5){continue;} // Skip events where p' MC theta is below 0.5mrad 
    h1_Theta_MC->Fill(Vec_pSc_MC.Theta()*1000);
    h1_eta_MC->Fill(Vec_pSc_MC.eta());
    if(0.005 > Vec_pSc_MC.Theta()*1000){
      h1_Theta_MC_RP->Fill(Vec_pSc_MC.Theta()*1000);
      h1_eta_MC_RP->Fill(Vec_pSc_MC.eta());      
    }
    if(0.0055 < Vec_pSc_MC.Theta()*1000 && 0.02 > Vec_pSc_MC.Theta()*1000){
      h1_Theta_MC_B0->Fill(Vec_pSc_MC.Theta()*1000);
      h1_eta_MC_B0->Fill(Vec_pSc_MC.eta());
    }
    for(unsigned int i = 0; i < TSCPSim_Assoc.GetSize(); i++){
      if (TSCPSim_Assoc[i] == pSc_Index){ // If matching track for p' found, assign reconstructed p' 3 vector
	Vec_pSc_RecAssoc.SetXYZ(TSCP_Px[TSCPRec_Assoc[i]], TSCP_Py[TSCPRec_Assoc[i]], TSCP_Pz[TSCPRec_Assoc[i]]);
	Vec_pSc_RecAssoc_Rot = rY*Vec_pSc_RecAssoc;
	h1_Theta_RecAssoc->Fill(Vec_pSc_RecAssoc_Rot.Theta()*1000); 
	h1_eta_RecAssoc->Fill(Vec_pSc_RecAssoc_Rot.eta());
	h1_Theta_Rec->Fill(Vec_pSc_Rec_Rot.Theta()*1000);
	h1_eta_Rec->Fill(Vec_pSc_Rec_Rot.eta());
      	// Fill true MC info for particles that have reconstructed in B0
	h1_Theta_pScRecMC->Fill(Vec_pSc_MC.Theta()*1000);
	h1_Theta_MC_pScRecB0->Fill(Vec_pSc_MC.Theta()*1000);
	h1_eta_MC_pScRec->Fill(Vec_pSc_MC.eta());
	h1_eta_MC_pScRecB0->Fill(Vec_pSc_MC.eta());
	B0_Part=kTRUE;
      }
    } // End association loop
    if(B0_Part == kFALSE){ // If not in B0, check RP - No associations for RP (?)
      for(unsigned int i = 0; i < RP_Px.GetSize(); i++){
	//Vec_tmpRP.SetXYZ(RP_Px[i], RP_Py[i], RP_Pz[i]);
	if(RP_Px.GetSize() == 1){ // Only fill if just 1 RP reconstructed track (for now), do something more sensible in future
	  Vec_pSc_Rec.SetXYZ(RP_Px[i], RP_Py[i], RP_Pz[i]);
	  if(Vec_pSc_Rec.Theta()*1000 < 25){// Only fill if in FF region
	    h1_Theta_Rec->Fill(Vec_pSc_Rec.Theta()*1000);
	    h1_Theta_Rec_RP->Fill(Vec_pSc_Rec.Theta()*1000);
	    h1_eta_Rec->Fill(Vec_pSc_Rec.eta());
	    h1_eta_Rec_RP->Fill(Vec_pSc_Rec.eta());
	    // Fill true MC info for particles that have reconstructed in B0
	    h1_Theta_pScRecMC->Fill(Vec_pSc_MC.Theta()*1000);
	    h1_Theta_MC_pScRecRP->Fill(Vec_pSc_MC.Theta()*1000);
	    h1_eta_MC_pScRec->Fill(Vec_pSc_MC.eta());
	    h1_eta_MC_pScRecRP->Fill(Vec_pSc_MC.eta());
	  }
	}
      }
    }
  } // End event while loop

  TH1D* h1_Eff_Theta = new TH1D("h1_Eff_Theta","Eff(#theta_{MC}); #theta_{MC} (mrad); Efficiency (%)", 125, 0, 25);
  h1_Eff_Theta->Divide(h1_Theta_pScRecMC, h1_Theta_MC, 1, 1, "b");
  TH1D* h1_Eff_Theta_RP = new TH1D("h1_Eff_Theta_RP", "Eff(#theta_{MC}), p' in RP; #theta_{MC} (mrad); Efficiency (%)", 125, 0, 25);
  h1_Eff_Theta_RP->Divide(h1_Theta_MC_pScRecRP, h1_Theta_MC, 1, 1, "b");
  TH1D* h1_Eff_Theta_B0 = new TH1D("h1_Eff_Theta_B0", "Eff(#theta_{MC}), p' in B0; #theta_{MC} (mrad); Efficiency (%)",  125, 0, 25);
  h1_Eff_Theta_B0->Divide(h1_Theta_MC_pScRecB0, h1_Theta_MC, 1, 1, "b");

  TH1D* h1_Eff_eta = new TH1D("h1_Eff_eta","Eff(#eta_{MC}); #eta_{MC}; Efficiency (%)", 60, 4, 10);
  h1_Eff_eta->Divide(h1_eta_MC_pScRec, h1_eta_MC, 1, 1, "b");
  TH1D* h1_Eff_eta_RP = new TH1D("h1_Eff_eta_RP", "Eff(#eta_{MC}), p' in RP; #eta_{MC}; Efficiency (%)", 60, 4, 10);
  h1_Eff_eta_RP->Divide(h1_eta_MC_pScRecRP, h1_eta_MC, 1, 1, "b");
  TH1D* h1_Eff_eta_B0 = new TH1D("h1_Eff_eta_B0", "Eff(#eta_{MC}), p' in B0; #theta_{MC}; Efficiency (%)",  60, 4, 10);
  h1_Eff_eta_B0->Divide(h1_eta_MC_pScRecB0, h1_eta_MC, 1, 1, "b");
  
  TCanvas* c_Results;
  TString OutPdf = Form("%s_pSc_Rec_Results.pdf", InFile.Data());
  c_Results = new TCanvas("c_Results","pSc_Rec_Results", 100, 0, 2560, 1920);
  c_Results->Divide(2,2); 
  c_Results->cd(1);
  h1_Theta_MC->SetLineColor(kBlack);
  h1_Theta_MC->SetTitle("#theta Comparison");
  h1_Theta_MC->Draw("HISTERR");
  h1_Theta_Rec->SetLineColor(kP6Blue);
  h1_Theta_Rec->Draw("SAMEHISTERR");
  h1_Theta_Rec_RP->SetLineColor(kP6Yellow);
  h1_Theta_Rec_RP->Draw("SAMEHISTERR");
  h1_Theta_RecAssoc->SetLineColor(kP6Red);
  h1_Theta_RecAssoc->Draw("SAMEHISTERR");
  Leg_Comp->AddEntry(h1_Theta_MC, "MC");
  Leg_Comp->AddEntry(h1_Theta_Rec, "Rec");
  Leg_Comp->AddEntry(h1_Theta_Rec_RP, "RecRP");
  Leg_Comp->AddEntry(h1_Theta_RecAssoc, "RecB0");
  gPad->SetLogy(1);
  c_Results->cd(2);
  h1_eta_MC->SetLineColor(kBlack);
  h1_eta_MC->SetTitle("#eta Comparison");
  h1_eta_MC->Draw("HISTERR");
  h1_eta_Rec->SetLineColor(kP6Blue);
  h1_eta_Rec->Draw("SAMEHISTERR");
  h1_eta_Rec_RP->SetLineColor(kP6Yellow);
  h1_eta_Rec_RP->Draw("SAMEHISTERR");
  h1_eta_RecAssoc->SetLineColor(kP6Red);
  h1_eta_RecAssoc->Draw("SAMEHISTERR");
  gPad->SetLogy(1);
  c_Results->cd(4);
  Leg_Comp->Draw("SAME");
  c_Results->Print(OutPdf + "(");

  TCanvas* c_Results2;
  c_Results2 = new TCanvas("c_Results2","pSc_Rec_Results2", 100, 0, 2560, 1920);
  c_Results2->Divide(3,2); 
  c_Results2->cd(1);
  h1_Theta_MC->SetLineColor(kBlack);
  h1_Theta_MC->SetTitle("#theta_{MC} Comparison");
  h1_Theta_MC->Draw("HISTERR");
  h1_Theta_pScRecMC->SetLineColor(kP6Blue);
  h1_Theta_pScRecMC->Draw("SAMEHISTERR");
  h1_Theta_MC_pScRecRP->SetLineColor(kP6Yellow);  
  h1_Theta_MC_pScRecRP->Draw("SAMEHISTERR"); 
  h1_Theta_MC_pScRecB0->SetLineColor(kP6Red);  
  h1_Theta_MC_pScRecB0->Draw("SAMEHISTERR"); 
  gPad->SetLogy(1);
  Leg_Comp2->AddEntry(h1_Theta_MC, "MC Truth");
  Leg_Comp2->AddEntry(h1_Theta_pScRecMC, "MC Truth - Reconstructed p'");
  Leg_Comp2->AddEntry(h1_Theta_MC_pScRecRP, "MC Truth - Reconstructed p' in RP");
  Leg_Comp2->AddEntry(h1_Theta_MC_pScRecB0, "MC Truth - Reconstructed p' in B0");
  c_Results2->cd(2);
  h1_eta_MC->SetLineColor(kBlack);
  h1_eta_MC->SetTitle("#eta_{MC} Comparison");
  h1_eta_MC->Draw("HISTERR");
  h1_eta_MC_pScRec->SetLineColor(kP6Blue);
  h1_eta_MC_pScRec->Draw("SAMEHISTERR");
  h1_eta_MC_pScRecRP->SetLineColor(kP6Yellow);  
  h1_eta_MC_pScRecRP->Draw("SAMEHISTERR"); 
  h1_eta_MC_pScRecB0->SetLineColor(kP6Red);  
  h1_eta_MC_pScRecB0->Draw("SAMEHISTERR"); 
  gPad->SetLogy(1);
  c_Results2->cd(3);
  Leg_Comp2->Draw("SAME");
  c_Results2->cd(4);
  h1_Eff_Theta->SetLineColor(kP6Blue);
  h1_Eff_Theta->Sumw2();
  h1_Eff_Theta->Scale(100);
  h1_Eff_Theta->Draw("HISTERR");
  h1_Eff_Theta_RP->SetLineColor(kP6Yellow);  
  h1_Eff_Theta_RP->Sumw2();
  h1_Eff_Theta_RP->Scale(100);
  h1_Eff_Theta_RP->Draw("SAMEHISTERR"); 
  h1_Eff_Theta_B0->SetLineColor(kP6Red);
  h1_Eff_Theta_B0->Sumw2();
  h1_Eff_Theta_B0->Scale(100);
  h1_Eff_Theta_B0->Draw("SAMEHISTERR"); 
  Leg_Comp3->AddEntry(h1_Eff_Theta, "Efficiency");
  Leg_Comp3->AddEntry(h1_Eff_Theta_RP, "Efficiency - Reconstructed p' in RP");
  Leg_Comp3->AddEntry(h1_Eff_Theta_B0, "Efficiency - Reconstructed p' in B0");
  c_Results2->cd(5);
  h1_Eff_eta->SetLineColor(kP6Blue);
  h1_Eff_eta->Sumw2();
  h1_Eff_eta->Scale(100);
  h1_Eff_eta->Draw("HISTERR");
  h1_Eff_eta_RP->SetLineColor(kP6Yellow);  
  h1_Eff_eta_RP->Sumw2();
  h1_Eff_eta_RP->Scale(100);
  h1_Eff_eta_RP->Draw("SAMEHISTERR"); 
  h1_Eff_eta_B0->SetLineColor(kP6Red);
  h1_Eff_eta_B0->Sumw2();
  h1_Eff_eta_B0->Scale(100);
  h1_Eff_eta_B0->Draw("SAMEHISTERR"); 
  c_Results2->cd(6);
  Leg_Comp3->Draw("SAME");
  c_Results2->Print(OutPdf + ")");
  
  ofile->Write();
  ofile->Close();
  
}
