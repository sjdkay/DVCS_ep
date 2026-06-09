// ePIC DVCS analysis class definition
#include "preLoadLib.hh"

// Data model headers
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/utils/vector_utils.h"
#include "edm4hep/utils/kinematics.h"
#include "edm4eic/ClusterCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "podio/Frame.h"
#include "podio/ROOTReader.h"

// ROOT Includes
#include <TSystem.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <Math/LorentzVector.h>
#include <TTree.h>
#include <TString.h>

#include <fstream>

// Class header include
#include "ePIC_DVCS_TASK.h"

//----------------------------------------------------
//----------------------------------------------------
//                    CONSTRUCTORS
//----------------------------------------------------
//----------------------------------------------------

// Default constructor
ePIC_DVCS_TASK::ePIC_DVCS_TASK(){
}

// Specific constructor
ePIC_DVCS_TASK::ePIC_DVCS_TASK(TString camp, TString energy, TString sett){
  // Set stored campaign attributes
  setDate(camp);
  setEnergy(energy);
  setSetting(sett);

  setBeamMomenta();
}


//----------------------------------------------------
//----------------------------------------------------
//                    SETTERS
//----------------------------------------------------
//----------------------------------------------------

// Set input file list
void ePIC_DVCS_TASK::setInFileList(TString name){
  sInList = name;
  if(gSystem->AccessPathName(sInList) == kTRUE){ // Check if provided file list is openable, if it isn't, exit
    cout << "File list - " << sInList << " not found." << endl;
    cout << "Check pathing and re-run" << endl;
    std::exit(EXIT_FAILURE);
  }
  else if(gSystem->AccessPathName(sInList) == kFALSE){
    std::cout<<"Input file list used: "<<sInList<<std::endl;
  }
}

// Set output file name and create new
void ePIC_DVCS_TASK::setOutFile(TString name){
  std::cout<<"Output ROOT file: "<<name<<std::endl;
  fOutFile = new TFile(name,"RECREATE");
}

void ePIC_DVCS_TASK::setBeamMomenta(){
  if(sEnergy == "5x41"){
    fPBeam_p=41.0;
    fPBeam_e=5.0;
  }
  else if(sEnergy == "10x100"){
    fPBeam_p=100.0;
    fPBeam_e=10.0;
  }
  else if(sEnergy == "10x130"){
    fPBeam_p=130.0;
    fPBeam_e=10.0;
  }
  else if(sEnergy == "10x250"){
    fPBeam_p=250.0;
    fPBeam_e=10.0;
  }
  else if(sEnergy == "18x275"){
    fPBeam_p=275.0;
    fPBeam_e=18.0;
  }
  else{
    fPBeam_p=100.0;
    fPBeam_e=10.0;
  }
}

void ePIC_DVCS_TASK::setMomCutFactors(Float_t factore = 1.,Float_t factorp = 1. ){
  fPMaxFactor_e = factore;
  fPMaxFactor_p = factorp;
}

//----------------------------------------------------
//----------------------------------------------------
//                    APPLY CUTS
//----------------------------------------------------
//----------------------------------------------------

// Single particle cuts - electron
Bool_t ePIC_DVCS_TASK::applyCuts_Electron(P3EVector beame, std::vector<P3EVector> scate){
   Bool_t passCuts{kTRUE};
   
   // EVENT CUTS
   // Require single particle in final state
   if(scate.size() != 1) passCuts = kFALSE;
   //if(scate.size() == 0) passCuts = kFALSE;
   // Return out of function if array is not filled
   if(!passCuts) return passCuts;

   // KINEMATIC CUTS
   // 1. Momentum
   if(scate[0].P() > (fPBeam_e*fPMaxFactor_e)) passCuts = kFALSE;
   // 2. Q2
   fQ2 = calcQ2_Elec(beame, scate[0]);
   if(fQ2 < fMinQ2) passCuts = kFALSE;

   return passCuts;
}

// Single particle cuts - photon
Bool_t ePIC_DVCS_TASK::applyCuts_Photon(std::vector<P3EVector> scatg){
   Bool_t passCuts{kTRUE};
   
   // EVENT CUTS
   // Require single particle in final state
   if(scatg.size() != 1) passCuts = kFALSE;
   //if(scatg.size() == 0) passCuts = kFALSE;
   // Return out of function if array is not filled
   if(!passCuts) return passCuts;

   //----------------------------------
   // INSERT ANY OTHER PHOTON CUTS HERE
   //----------------------------------

   //if(scatg[0].Eta() < 1.4 || scatg[0].Eta() > 3.0) passCuts = kFALSE;

   return passCuts;
}

// Single particle cuts - proton
Bool_t ePIC_DVCS_TASK::applyCuts_Proton(std::vector<P3EVector> scatp, TString sProtonDet="all"){
  Bool_t passCuts{kTRUE};
  
  // EVENT CUTS
  // Require single particle in final state
  if(scatp.size() != 1) passCuts = kFALSE;
  // Return out of function if array is not filled
  if(!passCuts) return passCuts;
  
  // KINEMATIC CUTS
  // 1. Momentum
  if(scatp[0].P() > (fPBeam_p*fPMaxFactor_p)) passCuts = kFALSE;
  
  // 2. Scattered proton theta (ensure within B0, Roman Pots or 'all')
  // If invalid detector name used, consider all
  if(sProtonDet != "B0" && sProtonDet != "RP" && sProtonDet != "all") sProtonDet="all";
  Float_t fMinPTheta{0.};
  Float_t fMaxPTheta{0.};
  // Need to know beam proton energy for minimum theta in RP
  Float_t beamP{0};
  if(sEnergy == "5x41"){
    beamP=41.0;
  }
  else if(sEnergy == "10x100"){
    beamP=100.0;
  }
  else if(sEnergy == "10x130"){
    beamP=130.0;
  }
  else if(sEnergy == "18x275"){
    beamP=275.0;
  }
  else{
    beamP=100.0;
  }
  // RP momentum acceptance < 200 MeV
  
  // B0 angular acceptance: 5.5 mrad - 20 mrad
  if(sProtonDet == "B0"){
    fMinPTheta = 0.0055;
    fMaxPTheta = 0.02;
  }
  // RP angluar acceptance: < 5.0 mrad
  else if(sProtonDet == "RP"){
    fMinPTheta = 0.;
    fMaxPTheta = 0.005;
  }
  // Full FF proton acceptance: < 20 mrad
  else if(sProtonDet == "all"){
    fMinPTheta = 0.;
    fMaxPTheta = 0.02;
  } 
  if(scatp[0].Theta()<fMinPTheta || scatp[0].Theta()>fMaxPTheta) passCuts = kFALSE;

  // 3. Longitudinal momentum
  //if(scatp[0].Z() < 0.98*beamP) passCuts = kFALSE;

  // 4. Transverse momentum (200 MeV minimum)
  //if(scatp[0].Pt() < 0.2) passCuts = kFALSE;

  return passCuts;
}

// Event-level cuts (DVCS kinematics)
Bool_t ePIC_DVCS_TASK::applyCuts_DVCS(TString sProtonDet="all"){
  Bool_t passCuts{kTRUE};

  // 1. MAXIMUM T CUT FOR ROMAN POTS
  if(sProtonDet != "B0" && sProtonDet != "RP" && sProtonDet != "all") sProtonDet="all";
  //if(sProtonDet == "RP" && ft > fMaxt_RP) passCuts = kFALSE;

  // 2. BJORKEN X CUT (removing tail from reconstructed histogram)
  //if(TMath::Log10(fxB) < fxB_Tail) passCuts = kFALSE;
  
  // 3. MAXIMUM MISSING MASS^2
  if(TMath::Abs(fM2miss) > fMax_M2miss) passCuts = kFALSE;

  return passCuts;
}

// Combination of all cuts
Bool_t ePIC_DVCS_TASK::applyCuts_All(P3EVector beame, P3EVector beamp, vector<P3EVector> scate, vector<P3EVector> scatp, vector<P3EVector> scatg, TString sProtonDet="all"){
  Bool_t passCuts{kTRUE};
  
  // 1. Electron cuts
  // Need to calculate Q2 first - set to zero if no detected electron
  if(scate.size() == 0) fQ2 = 0;
  else fQ2 = calcQ2_Elec(beame, scate[0]);
  passCuts = applyCuts_Electron(beame,scate);
  // Exit from function if failure
  if(!passCuts) return passCuts;

  // 2. Photon cuts
  passCuts = applyCuts_Photon(scatg);
  // Exit from function if failure
  if(!passCuts) return passCuts;

  // 3. Proton cuts
  passCuts = applyCuts_Proton(scatp, sProtonDet);
  // Exit from function if failure
  if(!passCuts) return passCuts;

  // 4. Event cuts
  // Need to calculate t, xB and MM2 first (e'p'y final state guaranteed by this point)
  fxB = calcX_Elec(beame, scate[0], beamp);
  fM2miss = calcM2Miss_3Body(beame, beamp, scate[0], scatp[0], scatg[0]);
  ft = calcT_BABE(beamp, scatp[0]);
  passCuts = applyCuts_DVCS(sProtonDet);

  return passCuts;
}

//----------------------------------------------------
//----------------------------------------------------
//            UNDO AFTERBURNER PROCEDURE
//----------------------------------------------------
//----------------------------------------------------

// Undo AB and calculate boost vectors - DO THIS FIRST FOR EACH EVENT
// USE BEAM VECTORS
void ePIC_DVCS_TASK::undoAfterburnAndCalc(P3EVector& p, P3EVector& k){
  // Holding vectors for beam - undoing crossing angle ONLY
  //P3EVector p_beam(fXAngle*p.E(), 0., p.E(), p.E());
  //P3EVector e_beam(0., 0., -k.E(), k.E());
  
  // Correction Feb. '26 - Use nomimal beam momenta for ab calculation
  P3EVector p_beam(fXAngle*fPBeam_p, 0., fPBeam_p, fPBeam_p);
  P3EVector e_beam(0., 0., -fPBeam_e, fPBeam_e);
    
  // Define boost vector to CoM frame
  P3EVector CoM_boost = p_beam+e_beam;
  vBoostToCoM.SetXYZ(-CoM_boost.X()/CoM_boost.E(), -CoM_boost.Y()/CoM_boost.E(), -CoM_boost.Z()/CoM_boost.E());
  
  // Apply boost to beam vectors
  p_beam = boost(p_beam, vBoostToCoM);
  e_beam = boost(e_beam, vBoostToCoM);
  
  // Calculate rotation angles and create rotation objects
  fRotY = -1.0*TMath::ATan2(p_beam.X(), p_beam.Z());
  fRotX = 1.0*TMath::ATan2(p_beam.Y(), p_beam.Z());

  rotAboutY = RotationY(fRotY);
  rotAboutX = RotationX(fRotX);

  // Apply rotation to beam vectors
  p_beam = rotAboutY(p_beam);
  p_beam = rotAboutX(p_beam);
  e_beam = rotAboutY(e_beam);
  e_beam = rotAboutX(e_beam);

  // Define boost vector back to head-on frame
  P3EVector HoF_boost(0., 0., CoM_boost.Z(), CoM_boost.E());
  vBoostToHoF.SetXYZ(HoF_boost.X()/HoF_boost.E(), HoF_boost.Y()/HoF_boost.E(), HoF_boost.Z()/HoF_boost.E());

  // Apply boost back to head on frame to beam vectors
  p_beam = boost(p_beam, vBoostToHoF);
  e_beam = boost(e_beam, vBoostToHoF);

  // Make changes to input vectors
  p.SetPxPyPzE(p_beam.X(), p_beam.Y(), p_beam.Z(), calcE(p_beam.Vect(),fMass_proton));
  k.SetPxPyPzE(e_beam.X(), e_beam.Y(), e_beam.Z(), calcE(e_beam.Vect(),fMass_electron));
}

// Undo afterburn procedure only
void ePIC_DVCS_TASK::undoAfterburn(P3EVector& a){
  Float_t mass = a.M();
  
  // Undo AB procedure for single vector, a^{mu}
  a = boost(a, vBoostToCoM); // BOOST TO COM FRAME
  a = rotAboutY(a);          // ROTATE TO Z-AXIS
  a = rotAboutX(a);          // ROTATE TO Z-AXIS
  a = boost(a, vBoostToHoF); // BOOST BACK TO HEAD ON FRAME

  a.SetPxPyPzE(a.X(), a.Y(), a.Z(), calcE(a.Vect(),mass));
}

//----------------------------------------------------
//----------------------------------------------------
//              KINEMATIC CALCULATIONS
//----------------------------------------------------
//----------------------------------------------------

// Calculate angle between hadronic and leptonic planes (Trento phi)
// Using planes defined by [k, q] and [q, p']
// Source: Bachetta, A. et al; Phys. Rev. D (2004); eq. 16
// Double_t ePIC_DVCS_TASK::calcTrentoPhi_qg(P3EVector k, P3EVector p, P3EVector kprime, P3EVector gprime){  
//   // Before calculating angle, boost into target rest frame
//   //MomVector vTgtRest = p.BoostToCM();

//   // Before calculating angle, boost into gamma*-p rest frame
//   // Calculate q in lab frame
//   P3EVector q = (k-kprime);
//   // Boost vector
//   MomVector vTgtRest = (p+q).BoostToCM();

//   k = boost(k,vTgtRest);
//   kprime = boost(kprime,vTgtRest);
//   gprime = boost(gprime,vTgtRest);

//   MomVector k3 = k.Vect();
//   MomVector kp3 = kprime.Vect();
//   MomVector gp3 = gprime.Vect();
//   MomVector qhat3 = (k3-kp3).Unit();

//   // Define leptonic plane using virtual photon and scattered electron
//   MomVector lNorm = qhat3.Cross(kp3);
//   lNorm /= lNorm.R();
//   // Define hadronic plane using q vector and scattered photon
//   MomVector hNorm = qhat3.Cross(gp3);
//   hNorm /= hNorm.R();

//   // Angle() function just returns magnitude of angle
//   // If photon vector has a component parallel to the leptonic normal, should be positive. If opposite, negative.
//   float phi = TMath::Sign(1.,gp3.Dot(lNorm))*Angle(lNorm,hNorm);

//   if (phi < 0) return phi+2*TMath::Pi();
//   else return phi;
// }

// Calculate angle between planes of qp and qg
Double_t ePIC_DVCS_TASK::calcPhiQPQG(P3EVector k, P3EVector p, P3EVector kprime, P3EVector gprime){
  MomVector p3 = p.Vect();
  MomVector k3 = k.Vect();
  MomVector kp3 = kprime.Vect();
  MomVector gp3 = gprime.Vect();
  MomVector q3 = k3-kp3;

  MomVector qp = q3.Cross(p3);
  MomVector qg = q3.Cross(gp3);

  return Angle(qp,qg);
}

// Calculation of cone angle (angle between measured photon and expected photon)
Double_t ePIC_DVCS_TASK::calcConeAngle(P3EVector k, P3EVector p, P3EVector kprime, P3EVector pprime, P3EVector gprime){
  // Initial state vectors
  MomVector p3 = p.Vect(); // Proton beam
  MomVector k3 = k.Vect(); // Electron beam
  MomVector pi = p3+k3;    // Total initial momenta
  // Final state vectors
  MomVector pp3 = pprime.Vect(); // Scattered proton
  MomVector kp3 = kprime.Vect(); // Scattered electron
  MomVector gp3 = gprime.Vect(); // Real photon
  MomVector pfkp = pp3+kp3;      // Sum of scattered proton and scattered electron momenta

  // Return angle between real photon and expected photon (calc. from missing momentum of ep system)
  MomVector gExpected = pi-pfkp;

  return Angle(gp3,gExpected);
}

//----------------------------------------------------
//----------------------------------------------------
//                     DO ANALYSIS
//----------------------------------------------------
//----------------------------------------------------

void ePIC_DVCS_TASK::doAnalysis(){

  //---------------------------------------------------------
  // Setup: Load input file list
  //---------------------------------------------------------
  // IF TESTING, LOAD TEST FILE LIST
  //sInList="./filelists/inputFileList_BGastestSingle.list";
  //sInList="./filelists/inputFileList_BGastest.list";

  ifstream fileListStream;
  fileListStream.open(sInList);
  string fileName;
  TFile* inputRootFile;

  int NumFiles = 0;  
  // Determine number of files in file list
  while(getline(fileListStream,fileName)){
    NumFiles++;
  }
  cout << "Proccessing - " << NumFiles << " total files in input file list" << endl;
  // Reset file
  fileListStream.clear();
  fileListStream.seekg(0, fileListStream.beg);

  //---------------------------------------------------------
  // Setup: Declare histograms
  //---------------------------------------------------------
  // Starting with TDR histograms
  // 1a) Eta - MC particles
  TH1D* h_eta_MCp   = new TH1D("eta_mcp",";#eta_{p'}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_MCe   = new TH1D("eta_mce",";#eta_{e'}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_MCg   = new TH1D("eta_mcg",";#eta_{#gamma}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_BGe   = new TH1D("eta_bge",";#eta_{e'}(Bkg.)", 275, -11.0, 11.0);
  TH1D* h_eta_BGg   = new TH1D("eta_bgg",";#eta_{#gamma}(Bkg.)", 275, -11.0, 11.0);
  // 1a) Eta - MC accepted
  TH1D* h_eta_MCAp   = new TH1D("eta_mcap",";#eta_{p'}(MCA)", 275, -11.0, 11.0);
  TH1D* h_eta_MCAe   = new TH1D("eta_mcae",";#eta_{e'}(MCA)", 275, -11.0, 11.0);
  TH1D* h_eta_MCAg   = new TH1D("eta_mcag",";#eta_{#gamma}(MCA)", 275, -11.0, 11.0);
  // 1b) Eta - reco. particles
  TH1D* h_eta_RPp   = new TH1D("eta_RPp",";#eta_{p'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_RPPp   = new TH1D("eta_RPPp",";#eta_{p'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_RPe   = new TH1D("eta_RPe",";#eta_{e'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_RPg   = new TH1D("eta_RPg",";#eta_{#gamma}(Reco)", 275, -11.0, 11.0);
  // 1c) E/eta coverage
  TH2D* h_2D_EvEta_g = new TH2D("2d_eveta_g",";#eta_{#gamma};E_{#gamma} [GeV]",200,-4.,4.,100,0.,50.);
  TH2D* h_2D_EvEta_e = new TH2D("2d_eveta_e",";#eta_{e'};E_{e'} [GeV]",200,-4.,4.,(Int_t)4*fPBeam_e, 0., 2.*fPBeam_e);
  TH2D* h_2D_EvEta_p = new TH2D("2d_eveta_p",";#eta_{p'};E_{p'} [GeV]",150,4.,10.,(Int_t)4.*fPBeam_p, 0., 2.*fPBeam_p);
  
  // 2) Energy of e', gamma
  TH1D* h_E_MCe = new TH1D("e_mce", ";E_{e'}(MC) [GeV]", (Int_t)4*fPBeam_e, 0.0, 2*fPBeam_e);
  TH1D* h_E_MCg = new TH1D("e_mcg", ";E_{#gamma}(MC) [GeV]", 50, 0.0, 50.);
  TH1D* h_E_BGe = new TH1D("e_bge", ";E_{e'}(Bkg.) [GeV]", (Int_t)4*fPBeam_e, 0.0, 2*fPBeam_e);
  TH1D* h_E_BGg = new TH1D("e_bgg", ";E_{#gamma}(Bkg.) [GeV]", 50, 0.0, 50.);
  TH1D* h_E_RPe = new TH1D("e_rpe", ";E_{e'}(Reco.) [GeV]", (Int_t)4*fPBeam_e, 0.0, 2*fPBeam_e);
  TH1D* h_E_RPg = new TH1D("e_rpg", ";E_{#gamma}(Reco.) [GeV]", 50, 0.0, 50.);
  // 2b) Theta of e'/gamma
  TH1D* h_Theta_MCe = new TH1D("theta_mce", ";#Theta_{e'}(MC) [deg]", 300, 0., 180.);
  TH1D* h_Theta_MCg = new TH1D("theta_mcg", ";#Theta_{#gamma}(MC) [deg]", 300, 0., 180.);
  TH1D* h_Theta_BGe = new TH1D("theta_bge", ";#Theta_{e'}(Bkg.) [deg]", 300, 0., 180.);
  TH1D* h_Theta_BGg = new TH1D("theta_bgg", ";#Theta_{#gamma}(Bkg.) [deg]", 300, 0., 180.);
  TH1D* h_Theta_RPe = new TH1D("theta_rpe", ";#Theta_{e'}(Reco.) [deg]", 300, 0., 180.);
  TH1D* h_Theta_RPg = new TH1D("theta_rpg", ";#Theta_{#gamma}(Reco.) [deg]", 300, 0., 180.);

  // 3) Mandelstam t
  TH1D* h_t_MC = new TH1D("t_mc", ";|t|(MC) [(GeV/c^{2})^{2}]", 40, 0., 2.);
  TH1D* h_t_RP = new TH1D("t_rp", ";|t|(Reco.) [(GeV/c^{2})^{2}]", 40, 0., 2.);

  // 4) Background species genStatus
  TH1D* h_genStat_MCe = new TH1D("genstat_mce",";generatorStatus",7500,0,7500);
  TH1D* h_genStat_MCg = new TH1D("genstat_mcg",";generatorStatus",7500,0,7500);

  // 5) Effect of cuts on reco. particles
  const int nCuts{13};
  TString cutname[13] = {"inc-sing_ele","inc-sing_pho","inc-sing_pro",
			 "e'p'#gamma-mult","e'p'#gamma-Q^{2}","e'p'#gamma-#theta_{p}","e'p'#gamma-NegTrackVeto","e'p'#gamma-FFVeto","e'p'#gamma-MM^{2}",
			 "e'#gamma-mult","e'#gamma-BarrelPosVeto","e'#gamma-NegTrackVeto","e'#gamma-FFVeto"};
  TH1D* hPassCuts_Rec = new TH1D("passcuts_rec",";;",nCuts,0,nCuts);
  for(int bin{1}; bin<=nCuts; bin++) hPassCuts_Rec->GetXaxis()->SetBinLabel(bin,cutname[bin-1]);

  // For testing - MC electrons reco. w/out charge
  TH1D* h_E_FakePhot = new TH1D("e_fakephot", ";E_{e'->\"#gamma\"}(MC) [GeV]", (Int_t)4*fPBeam_e, 0.0, 2*fPBeam_e);
  // For testing - MC matched BG elec/phot w/in reco. arrays
  TH1D* h_eta_BGeReco = new TH1D("eta_bgereco",";#eta_{e'}(Bkg.|Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_BGgReco = new TH1D("eta_bggreco",";#eta_{#gamma}(Bkg.|Reco)", 275, -11.0, 11.0);

  //---------------------------------------------------------
  // Loop over files in list
  //---------------------------------------------------------
  int fileCounter{0};

  // 4-vectors for beam particles - need these defined outside of file loop
  P3EVector beame4(0,0,0,-1);     // Beam electron (generated)
  P3EVector beamp4(0,0,0,-1);     // Beam proton (generated)
  
  // Start file loop
  while(getline(fileListStream,fileName)){
    //std::cout<<"Input file "<<fileCounter<<" : "<<fileName<<std::endl;
    
    // Open podio reader
    if ( fileCounter % ( NumFiles / 10 ) == 0 ) {
      cout << "Processed " << setw(4) << ceil(((1.0*fileCounter)/(1.0*NumFiles))*100.0) << " % of Files - " << fileCounter << endl;
    }
    // New reader for each file
    auto reader = podio::ROOTReader();
    reader.openFile(fileName);
    Int_t nEntries = reader.getEntries("events");
    //std::cout<<"File has "<<nEntries<<" events..."<<std::endl;
    
    // Booleans - check for vetoes
    bool kBarrelPos_Rec{false};    // Positive track outside FF region?
    bool kOtherFFTrack_Rec{false}; // Track in OMD/ZDC?
    bool kNonElecNeg_Rec{false};   // Non-electron -ve tracks in barrel

    // Case of taking average beams from file
    if(!kUseEventBeams){
      // MUST DO THIS FIRST
      // Full run over tree in first file before anything else
      // Calculate beams from average of individual event beam particles
      if(fileCounter==0){
	// Accumulator variables
	P3EVector beame4_acc(0,0,0,-1);
	P3EVector beamp4_acc(0,0,0,-1);
	
	for(size_t ev = 0; ev < reader.getEntries("events"); ev++){
	  const auto event = podio::Frame(reader.readNextEntry("events"));
	  
	  // Define holding particles
	  edm4hep::MCParticle beame_evt;
	  edm4hep::MCParticle beamp_evt;

	  // LOOP AND FIND PARTICLE OBJECTS
	  // MCParticles
	  auto& mcparts = event.get<edm4hep::MCParticleCollection>("MCParticlesHeadOnFrameNoBeamFX");
	  for(const auto& mcp : mcparts){
	    if(mcp.getPDG() == 11 && mcp.getGeneratorStatus() == 4) beame_evt = mcp;
	    if(mcp.getPDG() == 2212 && mcp.getGeneratorStatus() == 4) beamp_evt = mcp;
	  } // BEAM PARTICLES FOR EVENT FOUND
	  
	  // Add to accumulator
	  XYZVector e3vec_temp(beame_evt.getMomentum().x,beame_evt.getMomentum().y,beame_evt.getMomentum().z);
	  P3EVector e4vec_temp(e3vec_temp.X(), e3vec_temp.Y(), e3vec_temp.Z(), calcE(e3vec_temp,fMass_electron));
	  XYZVector p3vec_temp(beamp_evt.getMomentum().x,beamp_evt.getMomentum().y,beamp_evt.getMomentum().z);
	  P3EVector p4vec_temp(p3vec_temp.X(), p3vec_temp.Y(), p3vec_temp.Z(), calcE(p3vec_temp,fMass_proton));

	  beame4_acc += e4vec_temp;
	  beamp4_acc += p4vec_temp;	  
	} // End of event loop - FOR AVERAGED BEAMS
	
	// Divide by number of events in file
	beame4.SetCoordinates(beame4_acc.X()/nEntries, beame4_acc.Y()/nEntries, beame4_acc.Z()/nEntries, beame4_acc.E()/nEntries);
	beamp4.SetCoordinates(beamp4_acc.X()/nEntries, beamp4_acc.Y()/nEntries, beamp4_acc.Z()/nEntries, beamp4_acc.E()/nEntries);

	undoAfterburnAndCalc(beamp4,beame4);

	std::cout<<"First file - beams\n\te:"<<beame4<<"\n\tp:"<<beamp4<<std::endl;
      } // fi (fileCounter == 0)
      //else std::cout<<"Using beams from first file."<<std::endl;
    } // fi (!kUseEventBeams)

    // (Re)Run reader for main events
    for(size_t ev = 0; ev < reader.getEntries("events"); ev++){
      //for(size_t ev = 0; ev < 1; ev++){
      // Load next event
      //const auto event = podio::Frame(reader.readNextEntry("events"));
      const auto event = podio::Frame(reader.readEntry("events",ev));
      
      // 4-vectors for MC raw particles
      vector<P3EVector> scate4_gen;   // Scattered electron (generated)
      vector<P3EVector> scatp4_gen;   // Scattered proton (generated)
      vector<P3EVector> scatg4_gen;   // Scattered photon (generated)
      vector<P3EVector> bgasg4_gen;   // Beam gas photons
      vector<P3EVector> bgase4_gen;   // Background electrons
      // 4-vectors for associated MC particles (ONLY SCATTERED)
      vector<P3EVector> scate4_aso;   // Scattered electron (associated MC)
      vector<P3EVector> scatp4_aso;   // Scattered proton (associated MC)
      vector<P3EVector> scatg4_aso;   // Scattered photon (associated MC)
      // 4-vectors for reconstructed particles (SEPARATE PROTONS FOR B0 AND ROMAN POTS)
      vector<P3EVector> scate4_rec;   // Scattered electron (reconstructed)
      vector<P3EVector> scatp4_rec;   // Scattered proton (B0 reconstructed)
      vector<P3EVector> scatp4_rom;   // Scattered proton (Roman Pots reconstructed)
      vector<P3EVector> scatg4_rec;   // Scattered photon (reconstructed)     
      
      // MC truth
      // NEED TO USE BASE MCPARTICLES BRANCH - HoFNoBFX DOESN'T RETAIN GENERATOR STATUS CODES PROPERLY
      auto& mcparts = event.get<edm4hep::MCParticleCollection>("MCParticles");
      for(const auto& mcp : mcparts){
	// If using beams per event, look for generatorStatus 4
	if(kUseEventBeams){
	  if(mcp.getPDG() == 11 && mcp.getGeneratorStatus() == 4){
	    beame4.SetCoordinates(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, 
				  calcE(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, fMass_electron));
	  }
	  if(mcp.getPDG() == 2212 && mcp.getGeneratorStatus() == 4){
	    beamp4.SetCoordinates(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, 
				  calcE(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, fMass_proton));
	  }
	  undoAfterburnAndCalc(beamp4,beame4);
	} // fi (kUseEventBeams)
	
	// Then look for rest of MC particles
	if(mcp.getGeneratorStatus() != 4){
	  P3EVector temp(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, 
			 calcE(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, mcp.getMass()));
	  // Need to undo afterburn (not using HoFNoBFX)
	  undoAfterburn(temp);

	  if(mcp.getGeneratorStatus() == 1){
	    if(mcp.getPDG() == 11)   scate4_gen.push_back(temp);
	    if(mcp.getPDG() == 22)   scatg4_gen.push_back(temp);
	    if(mcp.getPDG() == 2212) scatp4_gen.push_back(temp);
	  } // fi (mcp primary physics particles)

	  if(mcp.getGeneratorStatus() > 100){
	    if(mcp.getPDG() == 11)   bgase4_gen.push_back(temp);
	    if(mcp.getPDG() == 22){
	      bgasg4_gen.push_back(temp);
	      //cout<<"Beam gas photon: "<<temp<<endl;
	    }
	  }// fi (mcp beam bkg)
	  
	} // fi (mcp.getGeneratorStatus() != 4)
	
	
	// Plot generator status codes for e/gamma (not expecting any "background" protons
	if(mcp.getPDG() == 11) h_genStat_MCe->Fill(mcp.getGeneratorStatus());
	if(mcp.getPDG() == 22) h_genStat_MCg->Fill(mcp.getGeneratorStatus());
	
      } // END OF MCPARTICLES LOOP

      // Reconstructed and associated particles (electrons/photons)
      const auto& assocReco = event.get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");
      for(const auto& mcreco : assocReco){
	// Declare holding vector variables
	P3EVector temp_mca(0,0,0,0);
	P3EVector temp_rec(0,0,0,0);
	
	// CASE 1: Using explicit MC matching
	// Only fill arrays if associated MC exists
	if(kUseExplicitMatch){
	  // Electrons - from sim.getPDG() flag
	  if(mcreco.getSim().getGeneratorStatus() == 1 && mcreco.getSim().getPDG() == 11){
	    // Simulated mass comes from MC truth - can trust sim.getMass()
	    temp_mca.SetCoordinates(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z,
				    calcE(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z, mcreco.getSim().getMass()));
	    // Reconstructed mass comes from PID hypothesis - cannot trust rec.getMass()
	    // Set mass by hand
	    temp_rec.SetCoordinates(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z,
				    calcE(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z, fMass_electron));
	  
	    // Undo afterburner
	    undoAfterburn(temp_mca);
	    undoAfterburn(temp_rec);
	    
	    // Add to particle arrays
	    scate4_aso.push_back(temp_mca); 
	    scate4_rec.push_back(temp_rec); 
	  }
	  // Photons - from sim.getPDG() flag
	  else if(mcreco.getSim().getGeneratorStatus() == 1 && mcreco.getSim().getPDG() == 22){
	    // Set LorentzVector coordinates by hand - as for electron case (but for massless reco.)
	    temp_mca.SetCoordinates(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z,
				    calcE(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z, mcreco.getSim().getMass()));
	    temp_rec.SetCoordinates(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z,
				    calcE(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z, 0.));
	  
	    // Undo afterburner
	    undoAfterburn(temp_mca);
	    undoAfterburn(temp_rec);
	 
	    // Add to particle arrays
	    scatg4_aso.push_back(temp_mca); 
	    scatg4_rec.push_back(temp_rec); 
	  }
	} // fi (kUseExplicitMatch)
	
	// Not using explicit matching?
	// Use reconstructed particle properties for PID
	else{
	  // CASE 2: Using ePIC PID
	  if(kUsePID){
	    // PID CODE HERE
	    // LOOKING FOR ELECTRONS AND PHOTONS
	  } // fi (kUsePID)
	  
	  // CASE 3: Using other properties for PID (charge, E/p, etc.)
	  else{
	    // Look for electrons - start from Q = -1
	    if(mcreco.getRec().getCharge() == -1){
	      // Skip if particle is missing ECAL clusters
	      if(!mcreco.getRec().getClusters()) continue;
	      
	      // Choose electrons from E/p 0.8 - 1.2
	      float clus_e{0};
	      for(auto& clust : mcreco.getRec().getClusters()) clus_e += clust.getEnergy();
	      float eoverp = clus_e/edm4hep::utils::magnitude(mcreco.getRec().getMomentum());

	      if(eoverp < 0.8 || eoverp > 1.2){
		// VETO - any negative tracks that aren't electrons
		kNonElecNeg_Rec = true;
		continue;
	      }
	      
	      // IF WANTING TO ADD CLUSTER ISOLATION - DO HERE
	      // OTHER ELECTRON FINDER CUTS
	      
	      // Create Lorentz vectors
	      // Simulated mass comes from MC truth - can trust sim.getMass()
	      temp_mca.SetCoordinates(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z,
				      calcE(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z, mcreco.getSim().getMass()));
	      // Reconstructed mass comes from PID hypothesis - cannot trust rec.getMass()
	      // Set mass by hand
	      temp_rec.SetCoordinates(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z,
				      calcE(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z, fMass_electron));
	      
	      // Undo afterburner
	      undoAfterburn(temp_mca);
	      undoAfterburn(temp_rec);
	      
	      // Add to particle arrays
	      scate4_aso.push_back(temp_mca); 
	      scate4_rec.push_back(temp_rec);
	      
	      // Looking at MC matched bkg electrons
	      if(mcreco.getSim().getGeneratorStatus()>100) h_eta_BGeReco->Fill(temp_rec.Eta());

	    } // End of electron finding
	    
	    // Look for photons - Q = 0 
	    // Will also catch cases where MC electron is missing track, but not clusters (clusters w/o track reconstructs as neutral)
	    if(mcreco.getRec().getCharge() == 0){
	      // Set LorentzVector coordinates by hand - as for electron case (but for massless reco.)
	      temp_mca.SetCoordinates(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z,
				      calcE(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z, mcreco.getSim().getMass()));
	      temp_rec.SetCoordinates(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z,
				      calcE(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z, 0.));
	      
	      // Undo afterburner
	      undoAfterburn(temp_mca);
	      undoAfterburn(temp_rec);
	      
	      // Add to particle arrays
	      scatg4_aso.push_back(temp_mca); 
	      scatg4_rec.push_back(temp_rec);

	      // TESTING - fake photons
	      // MC electrons which miss track
	      if(mcreco.getSim().getCharge()==-1) h_E_FakePhot->Fill(temp_rec.E());
	      // Looking at MC matched bkg photons
	      if(mcreco.getSim().getGeneratorStatus()>100) h_eta_BGgReco->Fill(temp_rec.Eta());

	    } // End of neutral finding
	  } // fi (!kUsePID)
	  
	} // fi (!kUseExplicitMatch)

	// VETOES
	// Positive tracks in barrel (B0 tracks ONLY in `TruthSeededCharged`)
	if(mcreco.getRec().getCharge() == 1) kBarrelPos_Rec = true;
      } // End of ReconstructedParticleAssociations

      // Now look for protons
      // Using ReconstructedTruthSeededChargedParticles
      const auto& assocTSReco = event.get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedTruthSeededChargedParticleAssociations");
      for(const auto& mcreco : assocTSReco){
	// Declare holding vector variables
	P3EVector temp_mca(0,0,0,0);
	P3EVector temp_rec(0,0,0,0);
      
	// CASE 1: Using explicit MC matching
	// Only fill arrays if associated MC exists
	if(kUseExplicitMatch){
	  // Select protons from sim.getPDG() flag
	  if(mcreco.getSim().getGeneratorStatus() == 1 && mcreco.getSim().getPDG() == 2212){
	    // Simulated mass comes from MC truth - can trust sim.getMass()
	    temp_mca.SetCoordinates(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z,
				    calcE(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z, mcreco.getSim().getMass()));
	    // Reconstructed mass comes from PID hypothesis - cannot trust rec.getMass()
	    // Set mass by hand
	    temp_rec.SetCoordinates(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z,
				    calcE(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z, fMass_proton));
	  
	    // Undo afterburner
	    undoAfterburn(temp_mca);
	    undoAfterburn(temp_rec);
	    
	    // Add to particle arrays
	    scatp4_aso.push_back(temp_mca); 
	    scatp4_rec.push_back(temp_rec); 
	  } // fi (sim.PDG() == 2212)
	} // fi (kUseExplicitMatch)

	// Not using explicit matching?
	// Use reconstructed particle properties for PID
	else{
	  // CASE 2: Using ePIC PID
	  if(kUsePID){
	    // PID CODE HERE
	    // LOOKING FOR ELECTRONS AND PHOTONS
	  } // fi (kUsePID)

	  // CASE 3: Using charge for PID
	  // Will ensure particle has associated track
	  else{
	    if(mcreco.getRec().getCharge() == 1){
	      // Simulated mass comes from MC truth - can trust sim.getMass()
	      temp_mca.SetCoordinates(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z,
				      calcE(mcreco.getSim().getMomentum().x, mcreco.getSim().getMomentum().y, mcreco.getSim().getMomentum().z, mcreco.getSim().getMass()));
	      // Reconstructed mass comes from PID hypothesis - cannot trust rec.getMass()
	      // Set mass by hand
	      temp_rec.SetCoordinates(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z,
				      calcE(mcreco.getRec().getMomentum().x, mcreco.getRec().getMomentum().y, mcreco.getRec().getMomentum().z, fMass_proton));
	  
	      // Undo afterburner
	      undoAfterburn(temp_mca);
	      undoAfterburn(temp_rec);

	      // Add to particle arrays
	      scatp4_aso.push_back(temp_mca); 
	      scatp4_rec.push_back(temp_rec); 
	    } // End of Q = 1
	  } // fi (!kUsePID)

	} // fi (!kUseExplicitMatch)
      } // End of ReconstructedTruthSeededChargedParticlesAssociations

      // Now get Roman Pot tracks
      // Branches from edm4eic::ReconstructedParticle without MC associations
      const auto& RPparts = event.get<edm4eic::ReconstructedParticleCollection>("ForwardRomanPotRecParticles");
      for(const auto& rpreco : RPparts){
	// Declare holding vector
	P3EVector temp_rec(0,0,0,0);
	
	// Assume all Roman Pot tracks are from reconstructed protons
	temp_rec.SetCoordinates(rpreco.getMomentum().x, rpreco.getMomentum().y, rpreco.getMomentum().z,
				calcE(rpreco.getMomentum().x, rpreco.getMomentum().y, rpreco.getMomentum().z, fMass_proton));

	// No need to undo afterburner - this is built in to the track reconstruction
	// Add to particle arrays
	scatp4_rom.push_back(temp_rec);
      } // End of ForwardRomanPotRecParticles
      
      // VETOES
      // OMD
      const auto& OMDparts = event.get<edm4eic::ReconstructedParticleCollection>("ForwardOffMRecParticles");
      for(const auto& omdreco : OMDparts){
	// Check if track exists
	if(omdreco.getTracks()) kOtherFFTrack_Rec = true;
      }
      // ZDC
      const auto& ZDCparts = event.get<edm4eic::ReconstructedParticleCollection>("ReconstructedFarForwardZDCNeutrals");
      for(const auto& zdcreco : ZDCparts){
	// Check for associated clusters
	if(zdcreco.getClusters()) kOtherFFTrack_Rec = true;
      }


      //---------------------------------------------------------
      // Fill histograms
      //---------------------------------------------------------
      // Single species - MC generated
      for(int ele_ind{0}; ele_ind<scate4_gen.size(); ele_ind++){
	h_eta_MCe->Fill(scate4_gen[ele_ind].Eta());
	h_2D_EvEta_e->Fill(scate4_gen[ele_ind].Eta(), scate4_gen[ele_ind].E());	
	h_E_MCe->Fill(scate4_gen[ele_ind].E());
	h_Theta_MCe->Fill(scate4_gen[ele_ind].Theta()*TMath::RadToDeg());
      }
      for(int pho_ind{0}; pho_ind<scatg4_gen.size(); pho_ind++){
	h_eta_MCg->Fill(scatg4_gen[pho_ind].Eta());
	h_2D_EvEta_g->Fill(scatg4_gen[pho_ind].Eta(), scatg4_gen[pho_ind].E());
	h_E_MCg->Fill(scatg4_gen[pho_ind].E());
	h_Theta_MCg->Fill(scatg4_gen[pho_ind].Theta()*TMath::RadToDeg());
      }
      for(int pro_ind{0}; pro_ind<scatp4_gen.size(); pro_ind++){
	h_eta_MCp->Fill(scatp4_gen[pro_ind].Eta());
	h_2D_EvEta_p->Fill(scatp4_gen[pro_ind].Eta(), scatp4_gen[pro_ind].E());
      }
      for(int ele_bkg{0}; ele_bkg<bgase4_gen.size(); ele_bkg++){
	h_eta_BGe->Fill(bgase4_gen[ele_bkg].Eta());
	h_2D_EvEta_e->Fill(bgase4_gen[ele_bkg].Eta(), bgase4_gen[ele_bkg].E());
	h_E_BGe->Fill(bgase4_gen[ele_bkg].E());
	h_Theta_BGe->Fill(bgase4_gen[ele_bkg].Theta()*TMath::RadToDeg());
      }
      for(int pho_bkg{0}; pho_bkg<bgasg4_gen.size(); pho_bkg++){
	h_eta_BGg->Fill(bgasg4_gen[pho_bkg].Eta());
	h_2D_EvEta_g->Fill(bgasg4_gen[pho_bkg].Eta(), bgasg4_gen[pho_bkg].E());
	h_E_BGg->Fill(bgasg4_gen[pho_bkg].E());
	h_Theta_BGg->Fill(bgasg4_gen[pho_bkg].Theta()*TMath::RadToDeg());
      }
      // Single species - MC accepted
      for(int ele_ind{0}; ele_ind<scate4_aso.size(); ele_ind++)	h_eta_MCAe->Fill(scate4_aso[ele_ind].Eta());
      for(int pho_ind{0}; pho_ind<scatg4_aso.size(); pho_ind++)	h_eta_MCAg->Fill(scatg4_aso[pho_ind].Eta());
      for(int pro_ind{0}; pro_ind<scatp4_aso.size(); pro_ind++)	h_eta_MCAp->Fill(scatp4_aso[pro_ind].Eta());
      // Single species - raw reco.
      for(int ele_ind{0}; ele_ind<scate4_rec.size(); ele_ind++){
	h_eta_RPe->Fill(scate4_rec[ele_ind].Eta());
	h_E_RPe->Fill(scate4_rec[ele_ind].E());
	h_Theta_RPe->Fill(scate4_rec[ele_ind].Theta()*TMath::RadToDeg());
      }
      for(int pho_ind{0}; pho_ind<scatg4_rec.size(); pho_ind++){
	h_eta_RPg->Fill(scatg4_rec[pho_ind].Eta());
	h_E_RPg->Fill(scatg4_rec[pho_ind].E());
	h_Theta_RPg->Fill(scatg4_rec[pho_ind].Theta()*TMath::RadToDeg());
      }
      for(int pro_ind{0}; pro_ind<scatp4_rec.size(); pro_ind++)	h_eta_RPp->Fill(scatp4_rec[pro_ind].Eta());
      for(int pro_ind{0}; pro_ind<scatp4_rom.size(); pro_ind++)	h_eta_RPPp->Fill(scatp4_rom[pro_ind].Eta());
      
      // Mandelstam t
      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all")){
	h_t_MC->Fill(calcT_BABE(beamp4,scatp4_gen[0]));
      }
      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "B0") && scatp4_rom.size() == 0){
	h_t_RP->Fill(calcT_BABE(beamp4,scatp4_rec[0]));
      }
      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rom, scatg4_rec, "RP") && scatp4_rec.size() == 0){
	h_t_RP->Fill(calcT_BABE(beamp4,scatp4_rom[0]));
      }
      
      // Count no. of events which pass cuts - reconstructed only
      // Single particle multiplicities
      if(scate4_rec.size() == 1) hPassCuts_Rec->Fill(0);
      if(scatg4_rec.size() == 1) hPassCuts_Rec->Fill(1);
      if(scatp4_rec.size() == 1 && scatp4_rom.size() == 0) hPassCuts_Rec->Fill(2);
      if(scatp4_rom.size() == 1 && scatp4_rec.size() == 0) hPassCuts_Rec->Fill(2);
      // Full DVCS multiplicity (using B0 protons)...
      if(scate4_rec.size() == 1 && scatg4_rec.size() == 1 && scatp4_rec.size() == 1 && scatp4_rom.size() == 0){
	hPassCuts_Rec->Fill(3);
	
	// ...and Q2 > 1...
	fQ2 = calcQ2_Elec(beame4, scate4_rec[0]);
	if(fQ2 >= 1.){
	  hPassCuts_Rec->Fill(4);
	  
	  // ...and proton track theta cut (B0 region)...
	  if(scatp4_rec[0].Theta() >= 0.0055 && scatp4_rec[0].Theta() <= 0.02){
	    hPassCuts_Rec->Fill(5);
	    
	    // ...and no extra -ve tracks
	    if(!kNonElecNeg_Rec){
	      hPassCuts_Rec->Fill(6);
	      
	      // ...and no ZDC/OMD particles
	      if(!kOtherFFTrack_Rec){
		hPassCuts_Rec->Fill(7);
		
		// ..and missing mass cut
		float mm2 = calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0]);
		if(TMath::Abs(mm2) <= 1) hPassCuts_Rec->Fill(8);
      	      } // fi (ZDC/OMD veto)
	    }   // fi (non-electron -ve particles)
	  }     // fi (p' theta cut)
	}       // fi (Q2 > 1)
      }         // fi (e'p'gamma)
      // Full DVCS multiplicity (using RP tracks)...
      if(scate4_rec.size() == 1 && scatg4_rec.size() == 1 && scatp4_rom.size() == 1 && scatp4_rec.size() == 0){
	hPassCuts_Rec->Fill(3);
	
	// ...and Q2 > 1...
	fQ2 = calcQ2_Elec(beame4, scate4_rec[0]);
	if(fQ2 >= 1.){
	  hPassCuts_Rec->Fill(4);
	  
	  // ...and proton track theta cut (RP region)...
	  if(scatp4_rom[0].Theta() > 0 && scatp4_rom[0].Theta() <= 0.005){
	    hPassCuts_Rec->Fill(5);
	    
	    // ...and no extra -ve tracks
	    if(!kNonElecNeg_Rec){
	      hPassCuts_Rec->Fill(6);
	      
	      // ...and no ZDC/OMD particles
	      if(!kOtherFFTrack_Rec){
		hPassCuts_Rec->Fill(7);
		
		// ..and missing mass cut
		float mm2 = calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0]);
		if(TMath::Abs(mm2) <= 1) hPassCuts_Rec->Fill(8);
      	      } // fi (ZDC/OMD veto)
	    }   // fi (non-electron -ve particles)
	  }     // fi (p' theta cut)
	}       // fi (Q2 > 1)
      }         // fi (e'p'gamma)
      // e'gamma final state (fakes for eXBE reco.)
      if(scate4_rec.size() == 1 && scatg4_rec.size() == 1){
	hPassCuts_Rec->Fill(9);
	
	// ...and barrel +ve track veto
	if(!kBarrelPos_Rec){
	  hPassCuts_Rec->Fill(10);
	  
	  // ...and non-elec -ve veto
	  if(!kNonElecNeg_Rec){
	    hPassCuts_Rec->Fill(11);
	    
	    // ...and ZDC/OMD veto
	    if(!kOtherFFTrack_Rec) hPassCuts_Rec->Fill(12);  // fi (ZDC/OMD veto)
	  } // fi (non-elec -ve veto)
	}   // fi (barrel +ve tracks)
      }     // fi (e'gamma)
      
      
    } // END OF EVENT LOOP - MAIN LOOP
    
    fileCounter++;
  } // END OF FILE LIST
  
  //------------------------------------------------------------
  // Write to output file
  //------------------------------------------------------------
  fOutFile->cd();
  
  h_eta_MCp->Write();
  h_eta_MCe->Write();
  h_eta_MCg->Write();
  h_eta_BGe->Write();
  h_eta_BGg->Write();
  h_eta_MCAp->Write();
  h_eta_MCAe->Write();
  h_eta_MCAg->Write();
  h_eta_RPp->Write();
  h_eta_RPPp->Write();
  h_eta_RPe->Write();
  h_eta_RPg->Write();
  h_2D_EvEta_g->Write();
  h_2D_EvEta_e->Write();
  h_2D_EvEta_p->Write();
  h_E_MCe->Write();
  h_E_MCg->Write();
  h_E_BGe->Write();
  h_E_BGg->Write();
  h_E_RPe->Write();
  h_E_RPg->Write();
  h_Theta_MCe->Write();
  h_Theta_MCg->Write();
  h_Theta_BGe->Write();
  h_Theta_BGg->Write();
  h_Theta_RPe->Write();
  h_Theta_RPg->Write();
  h_t_MC->Write();
  h_t_RP->Write();
  h_genStat_MCe->Write();
  h_genStat_MCg->Write();
  hPassCuts_Rec->Write();
  h_E_FakePhot->Write();
  h_eta_BGeReco->Write();
  h_eta_BGgReco->Write();

  return;
}
