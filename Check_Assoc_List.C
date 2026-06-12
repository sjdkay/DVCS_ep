// Stephen JD Kay - University of York
// Stephen.Kay@york.ac.uk
// 12/06/26
// Script to process reco output and plot theta for Reconstructed Truth Seeded Charged Particles

using namespace ROOT::Math;

#include "TString.h"
#include "Check_Assoc.h"
#include "ePICStyle.C"
#include <vector>

void Check_Assoc(TString InFile = ""){

  if(CheckFile(InFile) == kFALSE){ // Check files exist, can be opened and contain tree with fn
    exit(1);
  }

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
  cout << "All files in file list - " << FileList << " - added to chain" << endl;

  TTreeReader tree_reader(AnalysisChain);
  TFile *ofile = TFile::Open(Form("%s_OutputHists.root", ((InFile.Remove(InFile.Length() - 5)).Remove(0, 24)).Data()),"RECREATE"); // 22 cuts the .eicrecon.edm4eic.root from the original input file name

  XYZVector Vec_pSc_MC;
  XYZVector Vec_pSc_Rec;
  XYZVector Vec_pSc_Rec_Rot;
  XYZVector Vec_pSc_RecAssoc;
  XYZVector Vec_pSc_RecAssoc_Rot;
  
  TH1D* h1_Theta_MC = new TH1D("h1_Theta_MC", "#theta_{MC}; #theta (mRad); Counts/0.1 mRad bin", 250, 0, 0.025);
  TH1D* h1_Theta_MC_RP = new TH1D("h1_Theta_MC_RP", "#theta_{MC} in RP; #theta (mRad); Counts/0.1 mRad bin", 250, 0, 0.025);
  TH1D* h1_Theta_MC_B0 = new TH1D("h1_Theta_MC_B0", "#theta_{MC} in B0; #theta (mRad); Counts/0.1 mRad bin", 250, 0, 0.025);
  TH1D* h1_Theta_MCAssoc = new TH1D("h1_Theta_MCAssoc", "#theta_{MCAssoc}; #theta (mRad); Counts/0.1 mRad bin", 250, 0, 0.025);
  TH1D* h1_Theta_RecAssoc = new TH1D("h1_Theta_RecAssoc", "#theta_{RecAssoc}; #theta (mRad); Counts/0.1 mRad bin", 250, 0, 0.025);
  TH1D* h1_Theta_Rec = new TH1D("h1_Theta_Rec", "#theta_{Rec}; #theta (mRad); Counts/0.1 mRad bin", 250, 0, 0.025);
  TH1D* h1_Theta_Rec_RP = new TH1D("h1_Theta_Rec_RP", "#theta_{Rec_RP}; #theta (mRad); Counts/0.1 mRad bin", 250, 0, 0.025);

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

  while(tree_reader.Next()){ // Loop over all events
    // Reset index and reset vector for each event
    int pSc_Index = 0;
    Bool_t B0_Part=kFALSE;
    // Find MC info
    for (unsigned int i = 0; i < MC_Status.GetSize(); i++){
      if(MC_Status[i] == 1 && MC_PID[i] == 2212){Vec_pSc_MC.SetXYZ(MC_Px[i],MC_Px[i],MC_Pz[i]); pSc_Index = i;}
    } // End MC loop'
    h1_Theta_MC->Fill(Vec_pSc_MC.Theta());
    if(0.005 > Vec_pSc_MC.Theta()){
      h1_Theta_MC_RP->Fill(Vec_pSc_MC.Theta());
    }
    if(0.0055 < Vec_pSc_MC.Theta() && 0.02 > Vec_pSc_MC.Theta()){
      h1_Theta_MC_B0->Fill(Vec_pSc_MC.Theta());
    }
    for(unsigned int i = 0; i < TSCPSim_Assoc.GetSize(); i++){
      if (TSCPSim_Assoc[i] == pSc_Index){ // If matching track for p' found, assign reconstructed p' 3 vector
	Vec_pSc_RecAssoc.SetXYZ(TSCP_Px[TSCPRec_Assoc[i]], TSCP_Py[TSCPRec_Assoc[i]], TSCP_Pz[TSCPRec_Assoc[i]]);
	Vec_pSc_RecAssoc_Rot = rY*Vec_pSc_RecAssoc;
	h1_Theta_RecAssoc->Fill(Vec_pSc_RecAssoc_Rot.Theta()); 
	h1_Theta_MCAssoc->Fill(Vec_pSc_MC.Theta());
	B0_Part=kTRUE;
      }
    } // End association loop
    if(B0_Part == kFALSE){ // If not in B0, check RP - No associations for RP (?)
      for(unsigned int i = 0; i < RP_Px.GetSize(); i++){
	Vec_pSc_Rec.SetXYZ(RP_Px[i], RP_Py[i], RP_Pz[i]);
	if(Vec_pSc_Rec.Theta() < 0.025){// Only fill if in FF region
	  h1_Theta_Rec->Fill(Vec_pSc_Rec.Theta());
	  h1_Theta_Rec_RP->Fill(Vec_pSc_Rec.Theta());
	}
      }
    }
    // Fill all reconstructed particles, associated or not
    for(unsigned int i = 0; i < TSCP_Px.GetSize(); i++){
      Vec_pSc_Rec.SetXYZ(TSCP_Px[i], TSCP_Py[i], TSCP_Pz[i]);
      Vec_pSc_Rec_Rot = rY*Vec_pSc_Rec;
      if(Vec_pSc_Rec_Rot.Theta() < 0.025){// Only fill if in FF region
	h1_Theta_Rec->Fill(Vec_pSc_Rec_Rot.Theta());
      }
    }
  } // End event while loop
  
  ofile->Write();
  ofile->Close();
  
}
