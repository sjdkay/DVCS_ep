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
  std::cout<<"Input file list used: "<<name<<std::endl;
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

  // 3. Proton cuts (if not semi-inclusive)
  if(sProtonDet != "semi"){
    passCuts = applyCuts_Proton(scatp, sProtonDet);
    // Exit from function if failure
    if(!passCuts) return passCuts;
  }

  // 4. Event cuts
  if(sProtonDet == "semi"){
    // E-pz cut (general to fully exclusive and semi-inclusive reco.)
    fEmPz = (scate[0]+scatg[0]).E() - (scate[0]+scatg[0]).Pz();
    if(fEmPz < fMin_EmPz || fEmPz > fMax_EmPz) passCuts = false;
    // Exit from function if failure
    if(!passCuts) return passCuts;
  }
  else{
    // E-pz cut (general to fully exclusive and semi-inclusive reco.)
    fEmPz = (scate[0]+scatg[0]+scatp[0]).E() - (scate[0]+scatg[0]+scatp[0]).Pz();
    if(fEmPz < fMin_EmPz || fEmPz > fMax_EmPz) passCuts = false;
    // Exit from function if failure
    if(!passCuts) return passCuts;

    // Last cuts - fully exclusive DVCS event kinematics
    // Need to calculate t, xB and MM2 first (e'p'y final state guaranteed by this point)
    fxB = calcX_Elec(beame, scate[0], beamp);
    fM2miss = calcM2Miss_3Body(beame, beamp, scate[0], scatp[0], scatg[0]);
    ft = calcT_BABE(beamp, scatp[0]);
    passCuts = applyCuts_DVCS(sProtonDet);
  }

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
// Using planes defined by [k, q] and [q, g']
// Source: Bachetta, A. et al; Phys. Rev. D (2004); eq. 16
Double_t ePIC_DVCS_TASK::calcTrentoPhi_qg(P3EVector k, P3EVector p, P3EVector kprime, P3EVector gprime){  
  // Before calculating angle, boost into target rest frame
  //MomVector vTgtRest = p.BoostToCM();

  // Before calculating angle, boost into gamma*-p rest frame
  // Calculate q in lab frame
  P3EVector q = (k-kprime);
  // Boost vector
  MomVector vTgtRest = (p+q).BoostToCM();

  k = boost(k,vTgtRest);
  kprime = boost(kprime,vTgtRest);
  gprime = boost(gprime,vTgtRest);

  MomVector k3 = k.Vect();
  MomVector kp3 = kprime.Vect();
  MomVector gp3 = gprime.Vect();
  MomVector qhat3 = (k3-kp3).Unit();

  // Define leptonic plane using virtual photon and scattered electron
  MomVector lNorm = qhat3.Cross(kp3);
  lNorm /= lNorm.R();
  // Define hadronic plane using q vector and scattered photon
  MomVector hNorm = qhat3.Cross(gp3);
  hNorm /= hNorm.R();

  // Angle() function just returns magnitude of angle
  // If photon vector has a component parallel to the leptonic normal, should be positive. If opposite, negative.
  float phi = TMath::Sign(1.,gp3.Dot(lNorm))*Angle(lNorm,hNorm);

  if (phi < 0) return phi+2*TMath::Pi();
  else return phi;
}

// DIFF CALC. G.HILL PHD THESIS 2008
// EQ. 5.2 (p102)

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
  //if(sSett != "hiAcc" && sSett != "hiDiv") sInList="./filelists/inputFileList_test.list";
  //if(sSett != "hiAcc" && sSett != "hiDiv") sInList="./filelists/inputFileList_testDec.list";
  
  ifstream fileListStream;
  fileListStream.open(sInList);
  string fileName;
  TFile* inputRootFile;
  
  //---------------------------------------------------------
  // Setup: Declare histograms
  //---------------------------------------------------------
  // Starting with TDR histograms
  // 1a) Eta - MC particles
  TH1D* h_eta_MCp   = new TH1D("eta_MCp",";#eta_{p'}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_MCe   = new TH1D("eta_MCe",";#eta_{e'}(MC)", 275, -11.0, 11.0);
  TH1D* h_eta_MCg   = new TH1D("eta_MCg",";#eta_{#gamma}(MC)", 275, -11.0, 11.0);
  // 1a) Eta - MC particles
  TH1D* h_eta_MCAp   = new TH1D("eta_MCAp",";#eta_{p'}(MCA)", 275, -11.0, 11.0);
  TH1D* h_eta_MCAe   = new TH1D("eta_MCAe",";#eta_{e'}(MCA)", 275, -11.0, 11.0);
  TH1D* h_eta_MCAg   = new TH1D("eta_MCAg",";#eta_{#gamma}(MCA)", 275, -11.0, 11.0);
  // 1b) Eta - reco. particles
  TH1D* h_eta_RPp   = new TH1D("eta_RPp",";#eta_{p'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_RPPp   = new TH1D("eta_RPPp",";#eta_{p'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_RPe   = new TH1D("eta_RPe",";#eta_{e'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_eta_RPg   = new TH1D("eta_RPg",";#eta_{#gamma}(Reco)", 275, -11.0, 11.0);
  // 1c) E/p for electron candidates (before cuts applied)
  TH1D* h_EoverP_elec = new TH1D("eoverp_elec",";E/p",100,0.5,1.5);
  // 1d) E/eta coverage
  TH2D* h_2D_EvEta_g = new TH2D("2d_eveta_g",";#eta_{#gamma};E_{#gamma} [GeV]",200,-4.,4.,100,0.,50.);
  TH2D* h_2D_EvEta_e = new TH2D("2d_eveta_e",";#eta_{e'};E_{e'} [GeV]",200,-4.,4.,(Int_t)4*fPBeam_e, 0., 2.*fPBeam_e);
  TH2D* h_2D_EvEta_p = new TH2D("2d_eveta_p",";#eta_{p'};E_{p'} [GeV]",150,4.,10.,(Int_t)4.*fPBeam_p, 0., 2.*fPBeam_p);
  // DEBUGGING
  TH1D* h_eta_FakePhot = new TH1D("eta_fakephot",";#eta_{e^{-} #rightarrow #gamma}(MC)", 275, -11.0, 11.0);
  TH2D* h_eta_ElecGamma_MC = new TH2D("eta_elecgamma_mc",";#eta_{e'}(MC);#eta_{#gamma}(MC)", 275, -11.0, 11.0, 275, -11.0, 11.0);
  TH2D* h_eta_ElecGamma_RP = new TH2D("eta_elecgamma_rp",";#eta_{e'}(Reco);#eta_{#gamma}(Reco)", 275, -11.0, 11.0, 275, -11.0, 11.0);

  // 2) Photon theta resolution
  TH1D* h_PhotRes_theta = new TH1D("photres_theta",";#theta_{#gamma}(Reco)-#theta_{#gamma}(MC) [deg]",360,-90,90);
  TH2D* h_PhotRes2D_theta = new TH2D("photres2d_theta",";#theta_{#gamma, MC} [deg]; #Delta#theta_{#gamma} [deg]",370,0,185,360,-90,90);

  // 3) t distribution
  TH1D* h_t_Truth  = new TH1D("t_truth" ,";|t|(MC) [(GeV/#it{c}^{2})^{2}]"           , 20, 0., 2.);
  TH1D* h_t_B0Acc  = new TH1D("t_b0acc" ,";|t|(MC|Reco. - B0) [(GeV/#it{c}^{2})^{2}]", 20, 0., 2.);
  TH1D* h_t_RPAcc  = new TH1D("t_rpacc" ,";|t|(MC|Reco. - RP) [(GeV/#it{c}^{2})^{2}]", 20, 0., 2.);
  TH1D* h_t_LCAcc  = new TH1D("t_lcacc" ,";|t_{e'#gamma}|(MC|Reco.) [(GeV/#it{c}^{2})^{2}]", 20, 0., 2.);
  TH1D* h_t_B0Reco = new TH1D("t_b0reco",";|t|(Reco. - B0) [(GeV/#it{c}^{2})^{2}]"   , 20, 0., 2.);
  TH1D* h_t_RPReco = new TH1D("t_rpreco",";|t|(Reco. - RP) [(GeV/#it{c}^{2})^{2}]"   , 20, 0., 2.);
  TH1D* h_t_LCReco = new TH1D("t_lcreco",";|t_{e''gamma}|(Reco.) [(GeV/#it{c}^{2})^{2}]"   , 20, 0., 2.);

  // 4) t resolution - as absolute or as percentage (plot as preferred)
  TH2D* h_tResB0_2d = new TH2D("tresb0_2d",";|t|_{MC} [(GeV/#it{c})^{2}];#Deltat [(GeV/#it{c})^{2}]", 20, 0., 2., 500, -5., 5.);
  TH2D* h_tResRP_2d = new TH2D("tresrp_2d",";|t|_{MC} [(GeV/#it{c})^{2}];#Deltat [(GeV/#it{c})^{2}]", 20, 0., 2., 500, -5., 5.);
  TH2D* h_tResB0Pct_2d = new TH2D("tresb0pct_2d",";|t|_{MC} [(GeV/#it{c})^{2}];#Deltat/t_{MC}", 20, 0., 2., 200, -1., 1.);
  TH2D* h_tResRPPct_2d = new TH2D("tresrppct_2d",";|t|_{MC} [(GeV/#it{c})^{2}];#Deltat/t_{MC}", 20, 0., 2., 200, -1., 1.);
  TH2D* h_tResLC_2d    = new TH2D("treslc_2d",";|t_{e'#gamma}|_{MC} [(GeV/#it{c})^{2}];#Deltat [(GeV/#it{c})^{2}]", 20, 0., 2., 500, -5., 5.);
  TH2D* h_tResLCPct_2d = new TH2D("treslcpct_2d",";|t_{e'#gamma}|_{MC} [(GeV/#it{c})^{2}];#Deltat/t_{MC}", 20, 0., 2., 200, -1., 1.);

  // 5) Inclusive event kinematics - distributions, 2D response and resolution
  // 5a) 1D distributions
  TH1D* h_Q2_MC   = new TH1D("q2_mc"  , ";Q^{2}(MC) [GeV^{2}]"     , 5500, 0., 110.);
  TH1D* h_Q2_Acc  = new TH1D("q2_acc" , ";Q^{2}(MC|Reco) [GeV^{2}]", 5500, 0., 110.);
  TH1D* h_Q2_Reco = new TH1D("q2_reco", ";Q^{2}(Reco) [GeV^{2}]"   , 5500, 0., 110.);
  TH1D* h_xB_MC   = new TH1D("xb_mc"  , ";x_{B}(MC)"     , 10000, 0., 1.);
  TH1D* h_xB_Acc  = new TH1D("xb_acc" , ";x_{B}(MC|Reco)", 10000, 0., 1.);
  TH1D* h_xB_Reco = new TH1D("xb_reco", ";x_{B}(Reco)"   , 10000, 0., 1.);
  TH1D* h_y_MC   = new TH1D("y_mc"  , ";y(MC)"     , 100, 0., 1.);
  TH1D* h_y_Acc  = new TH1D("y_acc" , ";y(MC|Reco)", 100, 0., 1.);
  TH1D* h_y_Reco = new TH1D("y_reco", ";y(Reco)"   , 100, 0., 1.);
  // 5b) 2D distributions
  TH2D* h_Q2_2d = new TH2D("q2_2d",";Q^{2}(MC) [GeV^{2}];Q^{2}(Reco.) [GeV^{2}]", 5500, 0., 100., 5500, 0., 100.);
  TH2D* h_xB_2d = new TH2D("xb_2d",";x_{B}(MC);x_{B}(Reco.)", 10000, 0., 1., 10000, 0., 1.);
  TH2D* h_y_2d  = new TH2D("y_2d" ,";y(MC);y(Reco.)", 100, 0.,  1.,100, 0., 1.);
  // 5c) Resolutions
  TH2D* h_PctResQ2 = new TH2D("q2_pctres",";Q^{2} [GeV^{2}];#DeltaQ^{2}/Q^{2}", 5500, 0., 100., 200, -1., 1.);
  TH2D* h_PctResxB = new TH2D("xb_pctres",";x_{B};#Deltax_{B}/x_{B}", 1e4, 0., 1., 1000, -5., 5.);
  TH2D* h_PctResy  = new TH2D("y_pctres" ,";y;#Deltay/y", 100, 0.,  1., 200, -1., 1.);
  // 5d) Cross-variable coverages
  TH2D* h_2D_xVQ2_MC = new TH2D("2d_xvq2_mc",";x_{B,MC};Q^{2}_{MC} [GeV^{2}]",1e4,0.,1.,200,0.,100.);
  TH2D* h_2D_xVQ2_RP = new TH2D("2d_xvq2_rp",";x_{B,Reco};Q^{2}_{Reco} [GeV^{2}]",1e4,0.,1.,200,0.,100.);

  // 6) Cut variables
  // 6a) Proton track theta
  TH1D* h_theta_MCp = new TH1D("theta_mcp", ";#theta_{p'}(MC) [mrad]", 200, 0., 50.);
  TH1D* h_theta_B0p = new TH1D("theta_b0p", ";#theta_{p'}(Reco. - B0) [mrad]", 200, 0., 50.);
  TH1D* h_theta_RPp = new TH1D("theta_RPp", ";#theta_{p'}(Reco. - RP) [mrad]", 200, 0., 50.);
  // 6b) Full event missing mass (squared)
  TH1D* h_M2miss3_MC  = new TH1D("M2miss3_MC" ,";M^{2}_{miss}(MC) [(GeV/#it{c}^{2})^{2}]",600,-150,150);
  TH1D* h_M2miss3_MCA = new TH1D("M2miss3_MCA",";M^{2}_{miss}(MC|Reco) [(GeV/#it{c}^{2})^{2}]",600,-150,150);
  TH1D* h_M2miss3_RP  = new TH1D("M2miss3_RP" ,";M^{2}_{miss}(Reco) [(GeV/#it{c}^{2})^{2}]",600,-150,150);
  // 6c) E-pz (to reject radiative events)
  TH1D* h_EmPz3_MC = new TH1D("empz3_mc",";E-p_{z}(e'p'#gamma - MC) [GeV/#it{c}]",600,0.,150.);
  TH1D* h_EmPz3_RP = new TH1D("empz3_rp",";E-p_{z}(e'p'#gamma - Reco.) [GeV/#it{c}]",600,0.,150.);
  TH1D* h_EmPz2_MC = new TH1D("empz2_mc",";E-p_{z}(e'#gamma - MC) [GeV/#it{c}]",600,0.,150.);
  TH1D* h_EmPz2_RP = new TH1D("empz2_rp",";E-p_{z}(e'#gamma - MC) [GeV/#it{c}]",600,0.,150.);

  // 7) Trento phi
  TH1D* h_TPhi_MC    = new TH1D("tphi_mc"    ,";#phi_{h}(MC) [rad]"     , 160, 0, 6.4);
  TH1D* h_TPhi_B0Acc = new TH1D("tphi_b0acc" ,";#phi_{h}(MC|Reco) [rad]", 160, 0, 6.4);
  TH1D* h_TPhi_RPAcc = new TH1D("tphi_rpacc" ,";#phi_{h}(MC|Reco) [rad]", 160, 0, 6.4);
  TH1D* h_TPhi_B0Reco = new TH1D("tphi_b0reco",";#phi_{h}(Reco) [rad]"   , 160, 0, 6.4);
  TH1D* h_TPhi_RPReco = new TH1D("tphi_rpreco",";#phi_{h}(Reco) [rad]"   , 160, 0, 6.4);

  TH1D* h_TPhiRes = new TH1D("tphires",";#Delta#phi_{h}", 160, 0, 6.4);

  // 8) 3D differential distributions: Trento phi
  // [Q2][xB][|t|]
  double q2edges[9]  = {1., 1.78, 3.16, 5.62, 10., 18., 32., 56., 100};
  double xBedges[12] = {1e-4, 2.5e-4, 6.3e-4, 1e-3, 2.5e-3, 6.3e-3, 1e-2, 2.5e-2, 6.3e-2, 0.1, 0.25, 0.7};
  double tedges[15]  = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.4, 1.6};
  // Count no. of bins (nEdges-1)
  const int nQ2bins = sizeof(q2edges)/sizeof(q2edges[0]) - 1;
  const int nxBbins = sizeof(xBedges)/sizeof(xBedges[0]) - 1;
  const int ntbins  = sizeof(tedges)/sizeof(tedges[0]) - 1;
  
  // Q2/xB/t - for average values
  TH1D* h_Q2Diff[nQ2bins][nxBbins][ntbins];
  TH1D* h_xBDiff[nQ2bins][nxBbins][ntbins];
  TH1D* h_tDiff[nQ2bins][nxBbins][ntbins];

  // Phi - actual distributions
  TH1D* h_TPhiDiff_MC[nQ2bins][nxBbins][ntbins];
  TH1D* h_TPhiDiff_B0Acc[nQ2bins][nxBbins][ntbins];
  TH1D* h_TPhiDiff_RPAcc[nQ2bins][nxBbins][ntbins];
  TH1D* h_TPhiDiff_B0Reco[nQ2bins][nxBbins][ntbins];
  TH1D* h_TPhiDiff_RPReco[nQ2bins][nxBbins][ntbins];

  /*TH1D* h_TPhiResDiff_B0[nQ2bins][nxBbins][ntbins];
    TH1D* h_TPhiResDiff_RP[nQ2bins][nxBbins][ntbins];*/
  
  for(int binq2{0}; binq2<nQ2bins; binq2++){
    for(int binxB{0}; binxB<nxBbins; binxB++){
      for(int bint{0}; bint<ntbins; bint++){
	h_Q2Diff[binq2][binxB][bint] = new TH1D(Form("q2diff[%i][%i][%i]",binq2,binxB,bint),
						Form("%.2e<Q^{2}<%.2e GeV^{2}, %.2e<x_{B}<%.2e, %.2f<|t|<%.2f;Q^{2}(MCA) [GeV^{2}];",
						     q2edges[binq2],q2edges[binq2+1],
						     xBedges[binxB],xBedges[binxB+1],
						     tedges[bint],tedges[bint+1]),
						550, 0., 110.);
	h_xBDiff[binq2][binxB][bint] = new TH1D(Form("xbdiff[%i][%i][%i]",binq2,binxB,bint),
						Form("%.2e<Q^{2}<%.2e GeV^{2}, %.2e<x_{B}<%.2e, %.2f<|t|<%.2f;x_{B};",
						     q2edges[binq2],q2edges[binq2+1],
						     xBedges[binxB],xBedges[binxB+1],
						     tedges[bint],tedges[bint+1]),
						10000, 0., 1.);
	h_tDiff[binq2][binxB][bint] = new TH1D(Form("tdiff[%i][%i][%i]",binq2,binxB,bint),
					       Form("%.2e<Q^{2}<%.2e GeV^{2}, %.2e<x_{B}<%.2e, %.2f<|t|<%.2f;|t|(MCA) [GeV^{2}];",
						    q2edges[binq2],q2edges[binq2+1],
						    xBedges[binxB],xBedges[binxB+1],
						    tedges[bint],tedges[bint+1]),
					       20, 0., 2.);
	
	h_TPhiDiff_MC[binq2][binxB][bint] = new TH1D(Form("tphi_mc[%i][%i][%i]",binq2,binxB,bint), 
						     Form("%.2e<Q^{2}<%.2e GeV^{2}, %.2e<x_{B}<%.2e, %.2f<|t|<%.2f;#phi_{h}(MC) [rad];",
							  q2edges[binq2],q2edges[binq2+1],
							  xBedges[binxB],xBedges[binxB+1],
							  tedges[bint],tedges[bint+1]), 
						     10, 0, 6.3);
	h_TPhiDiff_B0Acc[binq2][binxB][bint] = new TH1D(Form("tphi_b0acc[%i][%i][%i]",binq2,binxB,bint) , 
							Form("%.2e<Q^{2}<%.2e GeV^{2}, %.2e<x_{B}<%.2e, %.2f<|t|<%.2f;#phi_{h}(MC|Reco) [rad];",
							     q2edges[binq2],q2edges[binq2+1],
							     xBedges[binxB],xBedges[binxB+1],
							     tedges[bint],tedges[bint+1]), 
							10, 0, 6.3);
	h_TPhiDiff_RPAcc[binq2][binxB][bint] = new TH1D(Form("tphi_rpacc[%i][%i][%i]",binq2,binxB,bint) , 
							Form("%.2e<Q^{2}<%.2e GeV^{2}, %.2e<x_{B}<%.2e, %.2f<|t|<%.2f;#phi_{h}(MC|Reco) [rad];",
							     q2edges[binq2],q2edges[binq2+1],
							     xBedges[binxB],xBedges[binxB+1],
							     tedges[bint],tedges[bint+1]), 
							10, 0, 6.3);
	h_TPhiDiff_B0Reco[binq2][binxB][bint] = new TH1D(Form("tphi_b0reco[%i][%i][%i]",binq2,binxB,bint), 
							 Form("%.2e<Q^{2}<%.2e GeV^{2}, %.2e<x_{B}<%.2e, %.2f<|t|<%.2f;#phi_{h}(Reco) [rad];",
							      q2edges[binq2],q2edges[binq2+1],
							      xBedges[binxB],xBedges[binxB+1],
							      tedges[bint],tedges[bint+1]), 
							 10, 0, 6.3);
	h_TPhiDiff_RPReco[binq2][binxB][bint] = new TH1D(Form("tphi_rpreco[%i][%i][%i]",binq2,binxB,bint), 
							 Form("%.2e<Q^{2}<%.2e GeV^{2}, %.2e<x_{B}<%.2e, %.2f<|t|<%.2f;#phi_{h}(Reco) [rad];",
							      q2edges[binq2],q2edges[binq2+1],
							      xBedges[binxB],xBedges[binxB+1],
							      tedges[bint],tedges[bint+1]), 
							 10, 0, 6.3);

	/*h_TPhiResDiff_B0[binq2][binxB][bint] = new TH1D(Form("tphi_rpreco[%i][%i][%i]",binq2,binxB,bint), 
							 Form("%.2e<Q^{2}<%.2e GeV^{2}, %.2e<x_{B}<%.2e, %.2f<|t|%.2f;#phi_{h}(Reco) [rad];",
							      q2edges[binq2],q2edges[binq2+1],
							      xBedges[binxB],xBedges[binxB+1],
							      tedges[bint],tedges[bint+1]), 
							      10, 0, 6.3);*/

      } // End of t bin loop
    }// End of xB bin loop
  } // End of Q2 bin loop

  //---------------------------------------------------------
  // Loop over files in list
  //---------------------------------------------------------
  int fileCounter{0};

  // 4-vectors for beam particles - need these defined outside of file loop
  P3EVector beame4(0,0,0,-1);     // Beam electron (generated)
  P3EVector beamp4(0,0,0,-1);     // Beam proton (generated)
  
  // Start file loop
  while(getline(fileListStream,fileName)){
    std::cout<<"Input file "<<fileCounter<<" : "<<fileName<<std::endl;
    
    // Open podio reader
    // New reader for each file
    auto reader = podio::ROOTReader();
    reader.openFile(fileName);
    Int_t nEntries = reader.getEntries("events");
    std::cout<<"File has "<<nEntries<<" events..."<<std::endl;
    
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
      else std::cout<<"Using beams from first file."<<std::endl;
    } // fi (!kUseEventBeams)

    // (Re)Run reader for main events
    for(size_t ev = 0; ev < reader.getEntries("events"); ev++){
      // Load next event
      //const auto event = podio::Frame(reader.readNextEntry("events"));
      const auto event = podio::Frame(reader.readEntry("events",ev));
      
      // 4-vectors for MC raw particles
      vector<P3EVector> scate4_gen;   // Scattered electron (generated)
      vector<P3EVector> scatp4_gen;   // Scattered proton (generated)
      vector<P3EVector> scatg4_gen;   // Scattered photon (generated)
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
      auto& mcparts = event.get<edm4hep::MCParticleCollection>("MCParticlesHeadOnFrameNoBeamFX");
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
	if(mcp.getGeneratorStatus() == 1){
	  P3EVector temp(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, 
			 calcE(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, mcp.getMass()));
	  
	  if(mcp.getPDG() == 11)   scate4_gen.push_back(temp);
	  if(mcp.getPDG() == 22)   scatg4_gen.push_back(temp);
	  if(mcp.getPDG() == 2212) scatp4_gen.push_back(temp);
	} // fi (mcp.getGeneratorStatus() == 1)
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

	      h_EoverP_elec->Fill(eoverp);

	      if(eoverp < 0.8 || eoverp > 1.2) continue;
	      
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

	      if(mcreco.getSim().getPDG() != 22) h_eta_FakePhot->Fill(temp_rec.Eta());
	    } // End of neutral finding
	  } // fi (!kUsePID)
	  
	} // fi (!kUseExplicitMatch)
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
      
      //cout<<"[DEBUG] Particles found - filling histos"<<endl;

      //---------------------------------------------------------
      // Fill histograms
      //---------------------------------------------------------
      // Single-particle kinematics
      // Eta and inclusive Q2/xB/y

      // MC particles
      // Need Q2 for electron cuts
      if(scate4_gen.size() == 0) fQ2 = 0;
      else fQ2 = calcQ2_Elec(beame4, scate4_gen[0]);
      
      if(applyCuts_Electron(beame4,scate4_gen)){
	h_eta_MCe->Fill(scate4_gen[0].Eta());
	h_Q2_MC->Fill(fQ2);
	h_xB_MC->Fill(calcX_Elec(beame4, beamp4, scate4_gen[0]));
	h_y_MC->Fill(calcY_Elec(beame4, beamp4, scate4_gen[0]));

	h_2D_EvEta_e->Fill(scate4_gen[0].Eta(), scate4_gen[0].E());
	h_2D_xVQ2_MC->Fill(calcX_Elec(beame4, beamp4, scate4_gen[0]), fQ2);
      }
      if(applyCuts_Photon(scatg4_gen)){
	h_eta_MCg->Fill(scatg4_gen[0].Eta());
	h_2D_EvEta_g->Fill(scatg4_gen[0].Eta(), scatg4_gen[0].E());
      }
      if(applyCuts_Proton(scatp4_gen, "all")){
	h_eta_MCp->Fill(scatp4_gen[0].Eta());
	h_2D_EvEta_p->Fill(scatp4_gen[0].Eta(), scatp4_gen[0].E());
      }

      // For cut histogram - don't apply proton theta cut
      if(scatp4_gen.size() == 1){
	h_theta_MCp->Fill(scatp4_gen[0].Theta()*1000);
      }
      if(applyCuts_Electron(beame4,scate4_gen) && applyCuts_Photon(scatg4_gen)){
	h_eta_ElecGamma_MC->Fill(scate4_gen[0].Eta(),scatg4_gen[0].Eta());
	h_EmPz2_MC->Fill((scate4_gen[0]+scatg4_gen[0]).E() - (scate4_gen[0]+scatg4_gen[0]).Pz());
      }
	
      if(applyCuts_Electron(beame4,scate4_aso)) h_eta_MCAe->Fill(scate4_aso[0].Eta());
      if(applyCuts_Photon(scatg4_aso))	        h_eta_MCAg->Fill(scatg4_gen[0].Eta());
      if(applyCuts_Proton(scatp4_aso, "all"))	h_eta_MCAp->Fill(scatp4_aso[0].Eta());
      
      // Reconstructed particles and associated MC
      // Need Q2 for electron cuts
      if(scate4_rec.size() == 0) fQ2 = 0;
      else fQ2 = calcQ2_Elec(beame4, scate4_rec[0]);

      if(scate4_rec.size() == 1) h_Q2_Reco->Fill(fQ2);
      if(applyCuts_Electron(beame4,scate4_rec) && applyCuts_Electron(beame4,scate4_aso)){
	h_eta_RPe->Fill(scate4_rec[0].Eta());
	h_Q2_Acc->Fill(calcQ2_Elec(beame4, scate4_aso[0]));
	h_xB_Acc->Fill(calcX_Elec(beame4, beamp4, scate4_aso[0]));
	h_xB_Reco->Fill(calcX_Elec(beame4, beamp4, scate4_rec[0]));
	h_y_Acc->Fill(calcY_Elec(beame4, beamp4, scate4_aso[0]));
	h_y_Reco->Fill(calcY_Elec(beame4, beamp4, scate4_rec[0]));

	h_2D_xVQ2_RP->Fill(calcX_Elec(beame4, beamp4, scate4_rec[0]), fQ2);
      }
      if(applyCuts_Photon(scatg4_rec))       h_eta_RPg->Fill(scatg4_rec[0].Eta());
      if(applyCuts_Proton(scatp4_rec, "B0")) h_eta_RPp->Fill(scatp4_rec[0].Eta());
      if(applyCuts_Proton(scatp4_rom, "RP")) h_eta_RPPp->Fill(scatp4_rom[0].Eta());
      // For cut histogram - don't apply proton theta cut
      if(scatp4_rec.size() == 1 && scatp4_rom.size() == 0) h_theta_B0p->Fill(scatp4_rec[0].Theta()*1000);
      if(scatp4_rom.size() == 1 && scatp4_rec.size() == 0) h_theta_RPp->Fill(scatp4_rom[0].Theta()*1000);
      if(applyCuts_Electron(beame4,scate4_rec) && applyCuts_Photon(scatg4_rec)) h_eta_ElecGamma_RP->Fill(scate4_rec[0].Eta(),scatg4_rec[0].Eta());
      
      // 2D inclusive distributions
      if(applyCuts_Electron(beame4,scate4_rec) && applyCuts_Electron(beame4,scate4_aso)){
	double q2_mc = calcQ2_Elec(beame4, scate4_aso[0]);
	double q2_rec = calcQ2_Elec(beame4, scate4_rec[0]);
	double q2_frac = (q2_rec-q2_mc)/q2_mc;
	h_Q2_2d->Fill(q2_mc,q2_rec);
	h_PctResQ2->Fill(q2_mc, q2_frac);

	double xB_mc = calcX_Elec(beame4, beamp4, scate4_aso[0]);
	double xB_rec = calcX_Elec(beame4, beamp4, scate4_rec[0]);
	double xB_frac = (xB_rec-xB_mc)/xB_mc;
	h_xB_2d->Fill(xB_mc,xB_rec);
	h_PctResxB->Fill(xB_mc, xB_frac);

	double y_mc = calcY_Elec(beame4, beamp4, scate4_aso[0]);
	double y_rec = calcY_Elec(beame4, beamp4, scate4_rec[0]);
	double y_frac = (y_rec-y_mc)/y_mc;
	h_y_2d->Fill(y_mc,y_rec);
	h_PctResy->Fill(y_mc, y_frac);
      }

      // Photon theta resolution
      Float_t th_rec{0}, th_gen{0};
      if(applyCuts_Photon(scatg4_rec) && applyCuts_Photon(scatg4_aso)){
	th_gen = scatg4_aso[0].Theta()*TMath::RadToDeg();
	th_rec = scatg4_rec[0].Theta()*TMath::RadToDeg();
	
	h_PhotRes_theta->Fill(th_rec-th_gen);
	h_PhotRes2D_theta->Fill(th_gen, th_rec-th_gen);
      }

      //cout<<"[DEBUG] Inclusive histos filled"<<endl;

      // Full DVCS event distributions - Mandelstam t and Trento phi
      // MC truth
      if(applyCuts_Electron(beame4,scate4_gen) && applyCuts_Photon(scatg4_gen) && applyCuts_Proton(scatp4_gen, "all") 
	 && TMath::Abs(calcM2Miss_3Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0])) < fMax_M2miss) 
	h_EmPz3_MC->Fill((scate4_gen[0]+scatp4_gen[0]+scatg4_gen[0]).E() - (scate4_gen[0]+scatp4_gen[0]+scatg4_gen[0]).Pz());
	 
      //cout<<"[DEBUG] MC E-Pz (full evt.) filled"<<endl;

      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all")){
	// Calculations
	Float_t t_gen = calcT_BABE(beamp4,scatp4_gen[0]);
	Float_t tphi_gen = calcTrentoPhi_qg(beame4, beamp4, scate4_gen[0], scatg4_gen[0]);
	// Distributions - Q2/xB integrated
	h_t_Truth->Fill(t_gen);
	h_TPhi_MC->Fill(tphi_gen);

	// Q2/xB diff.
	// Need Q2 and xB for events
	Float_t q2_gen = calcQ2_Elec(beame4, scate4_gen[0]);
	Float_t xB_gen = calcX_Elec(beame4, beamp4, scate4_gen[0]);
	
	// Find global bin numbers
	int binq2{-1}, binxB{-1}, bint{-1};
	for(int q{0}; q<nQ2bins; q++){
	  if(q2_gen >= q2edges[q] && q2_gen < q2edges[q+1]){
	    binq2 = q;
	    break;
	  }
	}
	for(int x{0}; x<nxBbins; x++){
	  if(xB_gen >= xBedges[x] && xB_gen < xBedges[x+1]){
	    binxB = x;
	    break;
	  }
	}
	for(int t{0}; t<ntbins; t++){
	  if(t_gen >= tedges[t] && t_gen < tedges[t+1]){
	    bint = t;
	    break;
	  }
	}
	if((binq2!=-1) && (binxB!=-1) && (bint!=-1)) h_TPhiDiff_MC[binq2][binxB][bint]->Fill(tphi_gen);
      }
      
      //cout<<"[DEBUG] MC TPhi filled"<<endl;
 
      // Reconstructed and MC accepted - B0 only
      if(applyCuts_Electron(beame4,scate4_rec) && applyCuts_Photon(scatg4_rec) && applyCuts_Proton(scatp4_rec, "B0") 
	 && scatp4_rom.size() == 0 
	 && TMath::Abs(calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0])) < fMax_M2miss)
	 h_EmPz3_RP->Fill((scate4_rec[0]+scatp4_rec[0]+scatg4_rec[0]).E() - (scate4_rec[0]+scatp4_rec[0]+scatg4_rec[0]).Pz());

      //cout<<"[DEBUG] Reco. E-Pz (full evt., B0 proton) filled"<<endl;

      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "B0") && scatp4_rom.size() == 0){	
	// Calculations
	Float_t t_acc = calcT_BABE(beamp4,scatp4_aso[0]);
	Float_t t_rec = calcT_BABE(beamp4,scatp4_rec[0]);
	Float_t tphi_acc = calcTrentoPhi_qg(beame4, beamp4, scate4_aso[0], scatg4_aso[0]);
	Float_t tphi_rec = calcTrentoPhi_qg(beame4, beamp4, scate4_rec[0], scatg4_rec[0]);
	// Distributions - Q2/xB integrated
	h_t_B0Acc->Fill(t_acc);
	h_t_B0Reco->Fill(t_rec);
	h_TPhi_B0Acc->Fill(tphi_acc);
	h_TPhi_B0Reco->Fill(tphi_rec);

	h_TPhiRes->Fill(tphi_rec-tphi_acc);

	// Q2/xB diff. - accepted
	// Need Q2 and xB for events
	Float_t q2_acc = calcQ2_Elec(beame4, scate4_aso[0]);
	Float_t xB_acc = calcX_Elec(beame4, beamp4, scate4_aso[0]);

	// Find global bin numbers
	int binq2{-1}, binxB{-1}, bint{-1};
	for(int q{0}; q<nQ2bins; q++){
	  if(q2_acc >= q2edges[q] && q2_acc < q2edges[q+1]){
	    binq2 = q;
	    break;
	  }
	}
	for(int x{0}; x<nxBbins; x++){
	  if(xB_acc >= xBedges[x] && xB_acc < xBedges[x+1]){
	    binxB = x;
	    break;
	  }
	}
	for(int t{0}; t<ntbins; t++){
	  if(t_acc >= tedges[t] && t_acc < tedges[t+1]){
	    bint = t;
	    break;
	  }
	}
	if((binq2!=-1) && (binxB!=-1) && (bint!=-1)) h_TPhiDiff_B0Acc[binq2][binxB][bint]->Fill(tphi_acc);
	
	// Q2/xB diff. - reconstructed
	// Need Q2 and xB for events
	Float_t q2_rec = calcQ2_Elec(beame4, scate4_rec[0]);
	Float_t xB_rec = calcX_Elec(beame4, beamp4, scate4_rec[0]);

	binq2 = -1;
	binxB = -1;
	bint = -1;
	for(int q{0}; q<nQ2bins; q++){
	  if(q2_rec >= q2edges[q] && q2_rec < q2edges[q+1]){
	    binq2 = q;
	    break;
	  }
	}
	for(int x{0}; x<nxBbins; x++){
	  if(xB_rec >= xBedges[x] && xB_rec < xBedges[x+1]){
	    binxB = x;
	    break;
	  }
	}
	for(int t{0}; t<ntbins; t++){
	  if(t_rec >= tedges[t] && t_rec < tedges[t+1]){
	    bint = t;
	    break;
	  }
	}
	if((binq2!=-1) && (binxB!=-1) && (bint!=-1)){
	  h_TPhiDiff_B0Reco[binq2][binxB][bint]->Fill(tphi_rec);
	  // Fill MCA histograms based on reco. kinematics
	  h_Q2Diff[binq2][binxB][bint]->Fill(q2_rec);
	  h_xBDiff[binq2][binxB][bint]->Fill(xB_rec);
	  h_tDiff[binq2][binxB][bint]->Fill(t_rec);
	}
      }
      
      //cout<<"[DEBUG] TPhi and differential histos. (B0 proton) filled"<<endl;

      // Reconstructed and accepted - RP only
      if(applyCuts_Electron(beame4,scate4_rec) && applyCuts_Photon(scatg4_rec) && applyCuts_Proton(scatp4_rom, "RP") 
	 && scatp4_rec.size() == 0 
	 && TMath::Abs(calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0])) < fMax_M2miss)
	 h_EmPz3_RP->Fill((scate4_rec[0]+scatp4_rom[0]+scatg4_rec[0]).E() - (scate4_rec[0]+scatp4_rom[0]+scatg4_rec[0]).Pz());

      //cout<<"[DEBUG] Reco. E-Pz (full evt., RP proton) filled"<<endl;

      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rom, scatg4_rec, "RP") && scatp4_rec.size() == 0){
	// Calculations
	Float_t t_acc = calcT_BABE(beamp4,scatp4_gen[0]);
	Float_t t_rec = calcT_BABE(beamp4,scatp4_rom[0]);
	Float_t tphi_acc = calcTrentoPhi_qg(beame4, beamp4, scate4_aso[0], scatg4_aso[0]);
	Float_t tphi_rec = calcTrentoPhi_qg(beame4, beamp4, scate4_rec[0], scatg4_rec[0]);
	// Distributions - Q2/xB integrated
       	h_t_RPAcc->Fill(t_acc);
	h_t_RPReco->Fill(t_rec);
	h_TPhi_RPAcc->Fill(tphi_acc);
	h_TPhi_RPReco->Fill(tphi_rec);	
	
	h_TPhiRes->Fill(tphi_rec-tphi_acc);

	// Q2/xB diff. - accepted
	// Need Q2 and xB for events
	Float_t q2_acc = calcQ2_Elec(beame4, scate4_aso[0]);
	Float_t xB_acc = calcX_Elec(beame4, beamp4, scate4_aso[0]);

	// Find global bin numbers
	int binq2{-1}, binxB{-1}, bint{-1};
	for(int q{0}; q<nQ2bins; q++){
	  if(q2_acc >= q2edges[q] && q2_acc < q2edges[q+1]){
	    binq2 = q;
	    break;
	  }
	}
	for(int x{0}; x<nxBbins; x++){
	  if(xB_acc >= xBedges[x] && xB_acc < xBedges[x+1]){
	    binxB = x;
	    break;
	  }
	}
	for(int t{0}; t<ntbins; t++){
	  if(t_acc >= tedges[t] && t_acc < tedges[t+1]){
	    bint = t;
	    break;
	  }
	}
	if((binq2!=-1) && (binxB!=-1) && (bint!=-1)) h_TPhiDiff_RPAcc[binq2][binxB][bint]->Fill(tphi_acc);

	// Q2/xB diff. - reconstructed
	Float_t q2_rec = calcQ2_Elec(beame4, scate4_rec[0]);
	Float_t xB_rec = calcX_Elec(beame4, beamp4, scate4_rec[0]);

	binq2 = -1;
	binxB = -1;
	bint = -1;
	for(int q{0}; q<nQ2bins; q++){
	  if(q2_rec >= q2edges[q] && q2_rec < q2edges[q+1]){
	    binq2 = q;
	    break;
	  }
	}
	for(int x{0}; x<nxBbins; x++){
	  if(xB_rec >= xBedges[x] && xB_rec < xBedges[x+1]){
	    binxB = x;
	    break;
	  }
	}
	for(int t{0}; t<ntbins; t++){
	  if(t_rec >= tedges[t] && t_rec < tedges[t+1]){
	    bint = t;
	    break;
	  }
	}
	
	//if(binxB==0) cout<<"Bin ["<<binq2<<"]["<<binxB<<"]["<<bint<<"]; x = "<<xB_rec<<endl;
	if((binq2!=-1) && (binxB!=-1) && (bint!=-1)){
	  h_TPhiDiff_RPReco[binq2][binxB][bint]->Fill(tphi_rec);
	  // Fill MCA histograms based on reco. kinematics
	  h_Q2Diff[binq2][binxB][bint]->Fill(q2_rec);
	  h_xBDiff[binq2][binxB][bint]->Fill(xB_rec);
	  h_tDiff[binq2][binxB][bint]->Fill(t_rec);
	}
      }
      
      //cout<<"[DEBUG] TPhi and differential histos. (B0 proton) filled"<<endl;

      // Semi-inclusive calculation - ignore if proton is detected or not
      if(applyCuts_Electron(beame4,scate4_rec) && applyCuts_Photon(scatg4_rec))	h_EmPz2_RP->Fill((scate4_rec[0]+scatg4_rec[0]).E() - (scate4_rec[0]+scatg4_rec[0]).Pz());

      //cout<<"[DEBUG] Reco. E-Pz (semi-inclusive) filled"<<endl;

      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "semi")){
	// Calculations
	Float_t t_acc = calcT_MethodL(beame4,beamp4,scate4_aso[0],fMass_proton,scatg4_aso[0]);
	Float_t t_rec = calcT_MethodL(beame4,beamp4,scate4_rec[0],fMass_proton,scatg4_rec[0]);
	
	// Distributions - Q2/xB integrated
	h_t_LCAcc->Fill(t_acc);
	h_t_LCReco->Fill(t_rec);
	
	// Q2/xB diff. - accepted
	Float_t q2_acc = calcQ2_Elec(beame4, scate4_aso[0]);
	Float_t xB_acc = calcX_Elec(beame4, beamp4, scate4_aso[0]);

	// Q2/xB diff. - reconstructed
	Float_t q2_rec = calcQ2_Elec(beame4, scate4_rec[0]);
	Float_t xB_rec = calcX_Elec(beame4, beamp4, scate4_rec[0]);
      }
      
      // Mandelstam t-resolution
      // ASSUME THAT GENERATED POSITIVE TRACK MATCHES RECONSTRUCTED POSITIVE TRACK
      Float_t t_rec{0}, t_gen{0};
      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all") 
	 && applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "B0") && scatp4_rom.size()==0){
	t_gen = calcT_BABE(beamp4, scatp4_gen[0]);
	t_rec = calcT_BABE(beamp4, scatp4_rec[0]);
	h_tResB0_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen));
	h_tResB0Pct_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen)/t_gen);
      }
      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all") 
	 && applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rom, scatg4_rec, "RP") && scatp4_rec.size()==0){
	t_gen = calcT_BABE(beamp4, scatp4_gen[0]);
	t_rec = calcT_BABE(beamp4, scatp4_rom[0]);
	h_tResRP_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen));
	h_tResRPPct_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen)/t_gen);
      }
      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "semi")
	 && applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "semi")){
	 // applyCuts_Electron(beame4,scate4_rec) && applyCuts_Photon(scatg4_rec)
	 //&& applyCuts_Electron(beame4,scate4_gen) && applyCuts_Photon(scatg4_gen)){
	t_gen = calcT_MethodL(beame4,beamp4,scate4_gen[0],fMass_proton,scatg4_gen[0]);
	t_rec = calcT_MethodL(beame4,beamp4,scate4_rec[0],fMass_proton,scatg4_rec[0]);
	h_tResLC_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen));
	h_tResLCPct_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen)/t_gen);
      }

      // Missing mass histos - only apply multiplicity and track quality cuts
      if(applyCuts_Electron(beame4,scate4_gen) && 
	 applyCuts_Photon(scatg4_gen)          &&
	 applyCuts_Proton(scatp4_gen, "all")     ) h_M2miss3_MC->Fill(calcM2Miss_3Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0]));
      if(applyCuts_Electron(beame4,scate4_aso) && 
	 applyCuts_Photon(scatg4_aso)          &&
	 applyCuts_Proton(scatp4_aso, "all")     ) h_M2miss3_MCA->Fill(calcM2Miss_3Body(beame4, beamp4, scate4_aso[0], scatp4_aso[0], scatg4_aso[0]));
      if(applyCuts_Electron(beame4,scate4_rec) && 
	 applyCuts_Photon(scatg4_rec)          &&
	 applyCuts_Proton(scatp4_rec, "B0")      ) h_M2miss3_RP->Fill(calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0]));
      if(applyCuts_Electron(beame4,scate4_rec) && 
	 applyCuts_Photon(scatg4_rec)          &&
	 applyCuts_Proton(scatp4_rom, "RP")      ) h_M2miss3_RP->Fill(calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0]));

    } // END OF EVENT LOOP - MAIN LOOP

    fileCounter++;
  } // END OF FILE LIST
  
  //------------------------------------------------------------
  // Write to output file
  //------------------------------------------------------------
  fOutFile->cd();
  
  // Eta - MC
  h_eta_MCp->Write();
  h_eta_MCe->Write();
  h_eta_MCg->Write();
  h_eta_MCAp->Write();
  h_eta_MCAe->Write();
  h_eta_MCAg->Write();
  // Eta - reco.
  h_eta_RPp->Write();
  h_eta_RPPp->Write();
  h_eta_RPe->Write();
  h_eta_RPg->Write();
  // Other single-particle plots
  h_EoverP_elec->Write();
  h_2D_EvEta_e->Write();
  h_2D_EvEta_g->Write();
  h_2D_EvEta_p->Write();
  // Photon theta resolution
  h_PhotRes_theta->Write();
  h_PhotRes2D_theta->Write();
  // t distributions
  h_t_Truth->Write();
  h_t_B0Acc->Write();
  h_t_RPAcc->Write();
  h_t_B0Reco->Write();
  h_t_RPReco->Write();
  h_t_LCAcc->Write();
  h_t_LCReco->Write();
  // 2D t resolution
  h_tResB0_2d->Write();
  h_tResRP_2d->Write();
  h_tResB0Pct_2d->Write();
  h_tResRPPct_2d->Write();
  h_tResLC_2d->Write();
  h_tResLCPct_2d->Write();
  // Trento phi
  h_TPhi_MC->Write();
  h_TPhi_B0Acc->Write();
  h_TPhi_RPAcc->Write();
  h_TPhi_B0Reco->Write();
  h_TPhi_RPReco->Write();
  h_TPhiRes->Write();
  // Inclusive kinematic distributions - 1D
  h_Q2_MC->Write();
  h_Q2_Acc->Write();
  h_Q2_Reco->Write();
  h_xB_MC->Write();
  h_xB_Acc->Write();
  h_xB_Reco->Write();
  h_y_MC->Write();
  h_y_Acc->Write();
  h_y_Reco->Write();
  // Inclusive kinematic distributions - 2D
  h_Q2_2d->Write();
  h_xB_2d->Write();
  h_y_2d->Write();
  // Inclusive kinematic resolutions
  h_PctResQ2->Write();
  h_PctResxB->Write();
  h_PctResy->Write();
  // Cross-variable coverages
  h_2D_xVQ2_MC->Write();
  h_2D_xVQ2_RP->Write();
  // Cuts plots
  h_theta_MCp->Write();
  h_theta_B0p->Write();
  h_theta_RPp->Write();
  h_M2miss3_MC->Write();
  h_M2miss3_MCA->Write();
  h_M2miss3_RP->Write();
  h_EmPz3_MC->Write();
  h_EmPz3_RP->Write();
  h_EmPz2_MC->Write();
  h_EmPz2_RP->Write();

  h_eta_FakePhot->Write();
  h_eta_ElecGamma_MC->Write();
  h_eta_ElecGamma_RP->Write();
  
  for(int q{0}; q<nQ2bins; q++){
    for(int x{0}; x<nxBbins; x++){
      for(int t{0}; t<ntbins; t++){
	h_Q2Diff[q][x][t]->Write();
	h_xBDiff[q][x][t]->Write();
	h_tDiff[q][x][t]->Write();
	h_TPhiDiff_MC[q][x][t]->Write();
	h_TPhiDiff_B0Acc[q][x][t]->Write();
	h_TPhiDiff_RPAcc[q][x][t]->Write();
	h_TPhiDiff_B0Reco[q][x][t]->Write();
	h_TPhiDiff_RPReco[q][x][t]->Write();
      }
    }// End of xB bin loop
  } // End of Q2 bin loop

  return;
}
