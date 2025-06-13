using namespace std;

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include "detectorResolution.h"

// Set histogram drawing options
void setHistOpts(TH1D* hist, TString type){
  // 1) Truth: kBlack, marker style 8 (filled circle), marker size 2
  // 2) Raw reco.: kBlue, marker style 71 (open circle, line thickness 2), marker size 2
  // 3) Corrected reco.: kGreen+3, marker style 33 (filled diamond), marker size 2

  if(type == "truth"){
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(8);
    hist->SetMarkerColor(kBlack);
  }
  else if(type == "reco"){
    hist->SetLineColor(kBlue);
    hist->SetMarkerStyle(24);
    hist->SetMarkerColor(kBlue);
  }
  else if(type == "corr"){
    hist->SetLineColor(kGreen+3);
    hist->SetMarkerStyle(33);
    hist->SetMarkerColor(kGreen+3);
  }
  else if(type == "scat"){
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(2);
    hist->SetMarkerColor(kBlack);
  }
  else{
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(8);
    hist->SetMarkerColor(kBlack);
  }

  return;
}

//---------------------------------------------------------------------
// MAIN
//---------------------------------------------------------------------
void QAPlots(TString campaign, TString energy, TString setting){
  //---------------------------------------------------------------------
  // Get data file from DVCSAnalysis script
  //---------------------------------------------------------------------
  TString inFileName = "$EIC_WORK_DIR/DVCS_Analysis/RootFiles/ePIC_DVCS_" + campaign + "_" + energy + "_QA.root";
  cout<<"Input file: "<<inFileName<<endl;
  TFile* inFile = new TFile(inFileName);

  //---------------------------------------------------------------------
  // Extract histograms
  //---------------------------------------------------------------------
  // Q2 
  TH1D* h_Q2_Truth = (TH1D*)inFile->Get("q2_truth");
  TH1D* h_Q2_Acc   = (TH1D*)inFile->Get("q2_acc");
  TH1D* h_Q2_Reco  = (TH1D*)inFile->Get("q2_reco");
  TH2D* h_Q2_Resp  = (TH2D*)inFile->Get("q2_resp");
  TH1D* h_Q2_Pur   = (TH1D*)inFile->Get("q2_pur");
  TH1D* h_Q2_Mig   = (TH1D*)inFile->Get("q2_mig");
  TH2D* h_dQ2vQ2   = (TH2D*)inFile->Get("dq2vq2");
  // x
  TH1D* h_xB_Truth = (TH1D*)inFile->Get("xb_truth");
  TH1D* h_xB_Acc   = (TH1D*)inFile->Get("xb_acc");
  TH1D* h_xB_Reco  = (TH1D*)inFile->Get("xb_reco");
  TH2D* h_xB_Resp  = (TH2D*)inFile->Get("xb_resp");
  TH1D* h_xB_Pur   = (TH1D*)inFile->Get("xb_pur");
  TH1D* h_xB_Mig   = (TH1D*)inFile->Get("xb_mig");
  TH2D* h_dxBvxB   = (TH2D*)inFile->Get("dxbvxb");
  // y
  TH1D* h_y_Truth = (TH1D*)inFile->Get("y_truth");
  TH1D* h_y_Acc   = (TH1D*)inFile->Get("y_acc");
  TH1D* h_y_Reco  = (TH1D*)inFile->Get("y_reco");
  TH2D* h_y_Resp  = (TH2D*)inFile->Get("y_resp");
  TH1D* h_y_Pur   = (TH1D*)inFile->Get("y_pur");
  TH1D* h_y_Mig   = (TH1D*)inFile->Get("y_mig");
  TH2D* h_dyvy    = (TH2D*)inFile->Get("dyvy");
  // Exclusive kinematic quantities
  TH1D* h_t_Truth  = (TH1D*)inFile->Get("t_truth");
  TH1D* h_t_B0Acc  = (TH1D*)inFile->Get("t_b0acc");
  TH1D* h_t_RPAcc  = (TH1D*)inFile->Get("t_rpacc");
  TH1D* h_t_B0Reco = (TH1D*)inFile->Get("t_b0reco");
  TH1D* h_t_RPReco = (TH1D*)inFile->Get("t_rpreco");
  TH2D* h_t_B0Resp = (TH2D*)inFile->Get("t_b0resp");
  TH2D* h_t_RPResp = (TH2D*)inFile->Get("t_rpresp");
  TH1D* h_t_B0Pur  = (TH1D*)inFile->Get("t_b0pur");
  TH1D* h_t_RPPur  = (TH1D*)inFile->Get("t_rppur");
  TH1D* h_t_B0Mig  = (TH1D*)inFile->Get("t_b0mig");
  TH1D* h_t_RPMig  = (TH1D*)inFile->Get("t_rpmig");
  TH2D* h_B0dtvt   = (TH2D*)inFile->Get("b0dtvt");
  TH2D* h_RPdtvt   = (TH2D*)inFile->Get("rpdtvt");
  //Single particle kinematics - protons
  TH1D* h_theta_p_Truth  = (TH1D*)inFile->Get("theta_p_truth");
  TH1D* h_theta_p_B0Acc  = (TH1D*)inFile->Get("theta_p_b0acc");
  TH1D* h_theta_p_RPAcc  = (TH1D*)inFile->Get("theta_p_rpacc");
  TH1D* h_theta_p_B0Reco = (TH1D*)inFile->Get("theta_p_b0reco");
  TH1D* h_theta_p_RPReco = (TH1D*)inFile->Get("theta_p_rpreco");
  TH2D* h_theta_p_B0Resp = (TH2D*)inFile->Get("theta_p_b0resp");
  TH2D* h_theta_p_RPResp = (TH2D*)inFile->Get("theta_p_rpresp");
  TH1D* h_theta_p_B0Pur  = (TH1D*)inFile->Get("theta_p_b0pur");
  TH1D* h_theta_p_RPPur  = (TH1D*)inFile->Get("theta_p_rppur");
  TH1D* h_theta_p_B0Mig  = (TH1D*)inFile->Get("theta_p_b0mig");
  TH1D* h_theta_p_RPMig  = (TH1D*)inFile->Get("theta_p_rpmig");
  TH1D* h_E_p_Truth  = (TH1D*)inFile->Get("E_p_truth");
  TH1D* h_E_p_B0Acc  = (TH1D*)inFile->Get("E_p_b0acc");
  TH1D* h_E_p_RPAcc  = (TH1D*)inFile->Get("E_p_rpacc");
  TH1D* h_E_p_B0Reco = (TH1D*)inFile->Get("E_p_b0reco");
  TH1D* h_E_p_RPReco = (TH1D*)inFile->Get("E_p_rpreco");
  TH2D* h_E_p_B0Resp = (TH2D*)inFile->Get("E_p_b0resp");
  TH2D* h_E_p_RPResp = (TH2D*)inFile->Get("E_p_rpresp");
  TH1D* h_E_p_B0Pur  = (TH1D*)inFile->Get("E_p_b0pur");
  TH1D* h_E_p_RPPur  = (TH1D*)inFile->Get("E_p_rppur");
  TH1D* h_E_p_B0Mig  = (TH1D*)inFile->Get("E_p_b0mig");
  TH1D* h_E_p_RPMig  = (TH1D*)inFile->Get("E_p_rpmig");
  //Single particle kinematics - electrons
  TH1D* h_theta_e_Truth = (TH1D*)inFile->Get("theta_e_truth");
  TH1D* h_theta_e_Acc   = (TH1D*)inFile->Get("theta_e_acc");  
  TH1D* h_theta_e_Reco  = (TH1D*)inFile->Get("theta_e_reco"); 
  TH2D* h_theta_e_Resp  = (TH2D*)inFile->Get("theta_e_resp"); 
  TH1D* h_theta_e_Pur   = (TH1D*)inFile->Get("theta_e_pur");
  TH1D* h_theta_e_Mig   = (TH1D*)inFile->Get("theta_e_mig");
  TH1D* h_E_e_Truth = (TH1D*)inFile->Get("E_e_truth");
  TH1D* h_E_e_Acc   = (TH1D*)inFile->Get("E_e_acc");
  TH1D* h_E_e_Reco  = (TH1D*)inFile->Get("E_e_reco");
  TH2D* h_E_e_Resp  = (TH2D*)inFile->Get("E_e_resp");
  TH1D* h_E_e_Pur   = (TH1D*)inFile->Get("E_e_pur");
  TH1D* h_E_e_Mig   = (TH1D*)inFile->Get("E_e_mig");
  //Single particle kinematics - photons
  TH1D* h_theta_g_Truth = (TH1D*)inFile->Get("theta_g_truth");
  TH1D* h_theta_g_Acc   = (TH1D*)inFile->Get("theta_g_acc");
  TH1D* h_theta_g_Reco  = (TH1D*)inFile->Get("theta_g_reco");
  TH2D* h_theta_g_Resp  = (TH2D*)inFile->Get("theta_g_resp");
  TH1D* h_theta_g_Pur   = (TH1D*)inFile->Get("theta_g_pur");
  TH1D* h_theta_g_Mig   = (TH1D*)inFile->Get("theta_g_mig");
  TH1D* h_E_g_Truth = (TH1D*)inFile->Get("E_g_truth");
  TH1D* h_E_g_Acc   = (TH1D*)inFile->Get("E_g_acc");
  TH1D* h_E_g_Reco  = (TH1D*)inFile->Get("E_g_reco");
  TH2D* h_E_g_Resp  = (TH2D*)inFile->Get("E_g_resp");
  TH1D* h_E_g_Pur   = (TH1D*)inFile->Get("E_g_pur");
  TH1D* h_E_g_Mig   = (TH1D*)inFile->Get("E_g_mig");
  // HFS quantities
  TH1D* h_HFSSigma_Truth = (TH1D*)inFile->Get("hfssigma_truth");
  TH1D* h_HFSSigma_Acc   = (TH1D*)inFile->Get("hfssigma_acc");
  TH1D* h_HFSSigma_Reco  = (TH1D*)inFile->Get("hfssigma_reco");
  TH2D* h_HFSSigma_Resp  = (TH2D*)inFile->Get("hfssigma_resp");
  TH1D* h_HFSSigma_Pur   = (TH1D*)inFile->Get("hfssigma_pur");
  TH1D* h_HFSSigma_Mig   = (TH1D*)inFile->Get("hfssigma_mig");
  TH1D* h_HFSpT2_Truth = (TH1D*)inFile->Get("hfspt2_truth");
  TH1D* h_HFSpT2_Acc   = (TH1D*)inFile->Get("hfspt2_acc");
  TH1D* h_HFSpT2_Reco  = (TH1D*)inFile->Get("hfspt2_reco");
  TH2D* h_HFSpT2_Resp  = (TH2D*)inFile->Get("hfspt2_resp");
  TH1D* h_HFSpT2_Pur   = (TH1D*)inFile->Get("hfspt2_pur");
  TH1D* h_HFSpT2_Mig   = (TH1D*)inFile->Get("hfspt2_mig");
  TH1D* h_FullSigma_Truth = (TH1D*)inFile->Get("fullsigma_truth");
  TH1D* h_FullSigma_Acc   = (TH1D*)inFile->Get("fullsigma_acc");
  TH1D* h_FullSigma_Reco  = (TH1D*)inFile->Get("fullsigma_reco");
  TH2D* h_FullSigma_Resp  = (TH2D*)inFile->Get("fullsigma_resp");
  TH1D* h_FullSigma_Pur   = (TH1D*)inFile->Get("fullsigma_pur");
  TH1D* h_FullSigma_Mig   = (TH1D*)inFile->Get("fullsigma_mig");
  // E/p for electron/photon
  TH1D* h_EOverp_e_Truth = (TH1D*)inFile->Get("eoverp_e_truth");
  TH1D* h_EOverp_e_Reco  = (TH1D*)inFile->Get("eoverp_e_reco");
  TH2D* h_EOverp_e_Resp  = (TH2D*)inFile->Get("eoverp_e_resp");
  TH1D* h_EOverp_g_Truth = (TH1D*)inFile->Get("eoverp_g_truth");
  TH1D* h_EOverp_g_Reco  = (TH1D*)inFile->Get("eoverp_g_reco");
  TH2D* h_EOverp_g_Resp  = (TH2D*)inFile->Get("eoverp_g_resp");
  // Inclusive t-resolutions
  TH2D* h_B0Inc_Res = (TH2D*)inFile->Get("b0incdt");
  TH2D* h_RPInc_Res = (TH2D*)inFile->Get("rpincdt");

  // Duplicate reco. histograms for detector corrected
  TH1D* h_Q2_Corr        = (TH1D*)h_Q2_Reco->Clone("q2_corr");
  TH1D* h_xB_Corr        = (TH1D*)h_xB_Reco->Clone("xb_corr");
  TH1D* h_y_Corr         = (TH1D*)h_y_Reco->Clone("y_corr");
  TH1D* h_t_B0Corr       = (TH1D*)h_t_B0Reco->Clone("t_b0corr");
  TH1D* h_t_RPCorr       = (TH1D*)h_t_RPReco->Clone("t_rpcorr");
  TH1D* h_theta_p_B0Corr = (TH1D*)h_theta_p_B0Reco->Clone("theta_p_b0corr");
  TH1D* h_theta_p_RPCorr = (TH1D*)h_theta_p_RPReco->Clone("theta_p_rpcorr");
  TH1D* h_E_p_B0Corr     = (TH1D*)h_E_p_B0Reco->Clone("E_p_b0corr");
  TH1D* h_E_p_RPCorr     = (TH1D*)h_E_p_RPReco->Clone("E_p_rpcorr");
  TH1D* h_theta_e_Corr   = (TH1D*)h_theta_e_Reco->Clone("theta_e_corr"); 
  TH1D* h_E_e_Corr       = (TH1D*)h_E_e_Reco->Clone("E_e_corr");
  TH1D* h_theta_g_Corr   = (TH1D*)h_theta_g_Reco->Clone("theta_g_corr");
  TH1D* h_E_g_Corr       = (TH1D*)h_E_g_Reco->Clone("E_g_corr");
  TH1D* h_HFSSigma_Corr  = (TH1D*)h_HFSSigma_Reco->Clone("hfssigma_corr");
  TH1D* h_HFSpT2_Corr    = (TH1D*)h_HFSpT2_Reco->Clone("hfspt2_corr");
  TH1D* h_FullSigma_Corr = (TH1D*)h_FullSigma_Reco->Clone("fullsigma_corr");
  
  //--------------------------------------------------------------------------------------
  // Extract detector corrected histograms
  //--------------------------------------------------------------------------------------
  // 1. Divide MC associated by MC generated - gets efficiency
  h_Q2_Acc->Divide(h_Q2_Truth);
  h_xB_Acc->Divide(h_xB_Truth);
  h_y_Acc->Divide(h_y_Truth);
  h_t_B0Acc->Divide(h_t_Truth);
  h_t_RPAcc->Divide(h_t_Truth);
  h_theta_p_B0Acc->Divide(h_theta_p_Truth);
  h_theta_p_RPAcc->Divide(h_theta_p_Truth);
  h_E_p_B0Acc->Divide(h_E_p_Truth);
  h_E_p_RPAcc->Divide(h_E_p_Truth);
  h_theta_e_Acc->Divide(h_theta_e_Truth);
  h_E_e_Acc->Divide(h_E_e_Truth);
  h_theta_g_Acc->Divide(h_theta_g_Truth);
  h_E_g_Acc->Divide(h_E_g_Truth);
  h_HFSSigma_Acc->Divide(h_HFSSigma_Truth);
  h_HFSpT2_Acc->Divide(h_HFSpT2_Truth);
  h_FullSigma_Acc->Divide(h_FullSigma_Truth);
  
  // 2. Divide (cloned) reconstructed by efficiency
  h_Q2_Corr->Divide(h_Q2_Acc);
  h_xB_Corr->Divide(h_xB_Acc);
  h_y_Corr->Divide(h_y_Acc);
  h_t_B0Corr->Divide(h_t_B0Acc);
  h_t_RPCorr->Divide(h_t_RPAcc);
  h_theta_p_B0Corr->Divide(h_theta_p_B0Acc);
  h_theta_p_RPCorr->Divide(h_theta_p_RPAcc);
  h_E_p_B0Corr->Divide(h_E_p_B0Acc);
  h_E_p_RPCorr->Divide(h_E_p_RPAcc);
  h_theta_e_Corr->Divide(h_theta_e_Acc);
  h_E_e_Corr->Divide(h_E_e_Acc);
  h_theta_g_Corr->Divide(h_theta_g_Acc);
  h_E_g_Corr->Divide(h_E_g_Acc);
  h_HFSSigma_Corr->Divide(h_HFSSigma_Acc);
  h_HFSpT2_Corr->Divide(h_HFSpT2_Acc);
  h_FullSigma_Corr->Divide(h_FullSigma_Acc);

  //--------------------------------------------------------------------------------------
  // Set histogram draw options
  //--------------------------------------------------------------------------------------
  // Q2
  setHistOpts(h_Q2_Truth, "truth");
  setHistOpts(h_Q2_Reco, "reco");
  setHistOpts(h_Q2_Corr, "corr");
  setHistOpts(h_Q2_Pur, "scat");
  setHistOpts(h_Q2_Mig, "scat");
  // x
  setHistOpts(h_xB_Truth, "truth");
  setHistOpts(h_xB_Reco, "reco");
  setHistOpts(h_xB_Corr, "corr");
  setHistOpts(h_xB_Pur, "scat");
  setHistOpts(h_xB_Mig, "scat");
  // y
  setHistOpts(h_y_Truth, "truth");
  setHistOpts(h_y_Reco, "reco");
  setHistOpts(h_y_Corr, "corr");
  setHistOpts(h_y_Pur, "scat");
  setHistOpts(h_y_Mig, "scat");
  // Exclusive kinematic quantities
  setHistOpts(h_t_Truth, "truth");
  setHistOpts(h_t_B0Reco, "reco");
  setHistOpts(h_t_RPReco, "reco");
  setHistOpts(h_t_B0Corr, "corr");
  setHistOpts(h_t_RPCorr, "corr");
  setHistOpts(h_t_B0Pur, "scat");
  setHistOpts(h_t_RPPur, "scat");
  setHistOpts(h_t_B0Mig, "scat");
  setHistOpts(h_t_RPMig, "scat");
  //Single particle kinematics - protons
  setHistOpts(h_theta_p_Truth, "truth");
  setHistOpts(h_theta_p_B0Reco, "reco");
  setHistOpts(h_theta_p_RPReco, "reco");
  setHistOpts(h_theta_p_B0Corr, "corr");
  setHistOpts(h_theta_p_RPCorr, "corr");
  setHistOpts(h_theta_p_B0Pur, "scat");
  setHistOpts(h_theta_p_RPPur, "scat");
  setHistOpts(h_theta_p_B0Mig, "scat");
  setHistOpts(h_theta_p_RPMig, "scat");
  setHistOpts(h_E_p_Truth, "truth");
  setHistOpts(h_E_p_B0Reco, "reco");
  setHistOpts(h_E_p_RPReco, "reco");
  setHistOpts(h_E_p_B0Corr, "corr");
  setHistOpts(h_E_p_RPCorr, "corr");
  setHistOpts(h_E_p_B0Pur, "scat");
  setHistOpts(h_E_p_RPPur, "scat");
  setHistOpts(h_E_p_B0Mig, "scat");
  setHistOpts(h_E_p_RPMig, "scat");
  //Single particle kinematics - electrons
  setHistOpts(h_theta_e_Truth, "truth");
  setHistOpts(h_theta_e_Reco, "reco");
  setHistOpts(h_theta_e_Corr, "corr");
  setHistOpts(h_theta_e_Pur, "scat");
  setHistOpts(h_theta_e_Mig, "scat");
  setHistOpts(h_E_e_Truth, "truth");
  setHistOpts(h_E_e_Reco, "reco");
  setHistOpts(h_E_e_Corr, "corr");
  setHistOpts(h_E_e_Pur, "scat");
  setHistOpts(h_E_e_Mig, "scat");
  //Single particle kinematics - photons
  setHistOpts(h_theta_g_Truth, "truth");
  setHistOpts(h_theta_g_Reco, "reco");
  setHistOpts(h_theta_g_Corr, "corr");
  setHistOpts(h_theta_g_Pur, "scat");
  setHistOpts(h_theta_g_Mig, "scat");
  setHistOpts(h_E_g_Truth, "truth");
  setHistOpts(h_E_g_Reco, "reco");
  setHistOpts(h_E_g_Corr, "corr");
  setHistOpts(h_E_g_Pur, "scat");
  setHistOpts(h_E_g_Mig, "scat");
  // HFS quantities
  setHistOpts(h_HFSSigma_Truth, "truth");
  setHistOpts(h_HFSSigma_Reco, "reco");
  setHistOpts(h_HFSSigma_Corr, "corr");
  setHistOpts(h_HFSSigma_Pur, "scat");
  setHistOpts(h_HFSSigma_Mig, "scat");
  setHistOpts(h_HFSpT2_Truth, "truth");
  setHistOpts(h_HFSpT2_Reco, "reco");
  setHistOpts(h_HFSpT2_Corr, "corr");
  setHistOpts(h_HFSpT2_Pur, "scat");
  setHistOpts(h_HFSpT2_Mig, "scat");
  setHistOpts(h_FullSigma_Truth, "truth");
  setHistOpts(h_FullSigma_Reco, "reco");
  setHistOpts(h_FullSigma_Corr, "corr");
  setHistOpts(h_FullSigma_Pur, "scat");
  setHistOpts(h_FullSigma_Mig, "scat");
  // E/p for electron/photon
  setHistOpts(h_EOverp_e_Truth, "truth");
  setHistOpts(h_EOverp_e_Reco, "reco");
  setHistOpts(h_EOverp_g_Truth, "truth");
  setHistOpts(h_EOverp_g_Reco, "reco");
    
  //--------------------------------------------------------------------------------------
  // Draw histograms
  //--------------------------------------------------------------------------------------
  // Draw histograms
  gStyle->SetOptStat(00000000);
  // Set header strings
  TString sHead1 = "#it{#bf{ePIC} Simulation}; "+campaign;
  TString sHead2 = "#it{ep " + energy + " GeV}";
  // Create header text objects
  TLatex* tHead1 = new TLatex(0.10, 0.91, sHead1);
  tHead1->SetNDC();
  tHead1->SetTextSize(25);
  tHead1->SetTextFont(43);
  tHead1->SetTextColor(kBlack);
  TLatex* tHead2 = new TLatex(0.63, 0.91, sHead2);
  tHead2->SetNDC();
  tHead2->SetTextSize(25);
  tHead2->SetTextFont(43);
  tHead2->SetTextColor(kBlack);
  
  // CANVAS 1: Q2
  TCanvas* cQ2 = new TCanvas("cQ2","",1200,800);
  cQ2->Divide(2,2);
  cQ2->cd(1);
  h_Q2_Truth->GetXaxis()->SetTitle("Q^{2} [(GeV/#it{c}^{2})^{2}]");
  h_Q2_Truth->GetXaxis()->SetTitleSize(0.05);
  h_Q2_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_Q2_Truth->GetXaxis()->SetLabelSize(0.05);
  h_Q2_Truth->Draw("ep");
  h_Q2_Reco->Draw("epsame");
  h_Q2_Corr->Draw("epsame");
  TLegend* lQ2 = new TLegend(0.55,0.49,0.87,0.83);
  lQ2->SetLineColorAlpha(kWhite,0);
  lQ2->SetLineWidth(0);
  lQ2->SetFillStyle(0);
  lQ2->SetFillColorAlpha(kWhite,0);
  lQ2->SetHeader("#splitline{Q^{2}; single e'}{Q^{2}#geq 1 GeV^{2}}");
  lQ2->AddEntry((TObject*)0, "", "");
  lQ2->AddEntry(h_Q2_Truth,"MC truth","pl");
  lQ2->AddEntry(h_Q2_Reco,"Reco.","pl");
  lQ2->AddEntry(h_Q2_Corr,"Reco., #varepsilon corr.","pl");
  lQ2->Draw();
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cQ2->cd(2);
  h_Q2_Resp->GetXaxis()->SetTitleSize(0.05);
  h_Q2_Resp->GetXaxis()->SetTitleOffset(0.9);
  h_Q2_Resp->GetXaxis()->SetLabelSize(0.05);
  h_Q2_Resp->GetYaxis()->SetTitleSize(0.05);
  h_Q2_Resp->GetYaxis()->SetTitleOffset(0.85);
  h_Q2_Resp->GetYaxis()->SetLabelSize(0.05);
  h_Q2_Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cQ2->cd(3);
  h_Q2_Pur->GetXaxis()->SetTitleSize(0.05);
  h_Q2_Pur->GetXaxis()->SetTitleOffset(0.9);
  h_Q2_Pur->GetXaxis()->SetLabelSize(0.05);
  h_Q2_Pur->GetYaxis()->SetTitleSize(0.05);
  h_Q2_Pur->GetYaxis()->SetTitleOffset(0.88);
  h_Q2_Pur->Draw("p");
  cQ2->cd(4);
  h_dQ2vQ2->GetXaxis()->SetTitleSize(0.05);
  h_dQ2vQ2->GetXaxis()->SetTitleOffset(0.9);
  h_dQ2vQ2->GetXaxis()->SetLabelSize(0.05);
  h_dQ2vQ2->GetYaxis()->SetTitleSize(0.05);
  h_dQ2vQ2->GetYaxis()->SetTitleOffset(0.8);
  h_dQ2vQ2->GetYaxis()->SetLabelSize(0.05);
  h_dQ2vQ2->Draw("colz");
  
  cQ2->Print("DVCSQA_temp001.pdf");
  cQ2->Close();

  TCanvas* cxB = new TCanvas("cxB","",1200,800);
  cxB->Divide(2,2);
  cxB->cd(1);
  gPad->SetLogx();
  h_xB_Truth->GetXaxis()->SetTitle("x_{B}");
  h_xB_Truth->GetXaxis()->SetTitleSize(0.05);
  h_xB_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_xB_Truth->GetXaxis()->SetLabelSize(0.05);
  h_xB_Truth->Draw("ep");
  h_xB_Reco->Draw("epsame");
  h_xB_Corr->Draw("epsame");
  TLegend* lxB = new TLegend(0.56,0.50,0.89,0.84);
  lxB->SetLineColorAlpha(kWhite,0);
  lxB->SetLineWidth(0);
  lxB->SetFillStyle(0);
  lxB->SetFillColorAlpha(kWhite,0);
  lxB->SetHeader("#splitline{x_{B}; single e'}{Q^{2}#geq 1 GeV^{2}}");
  lxB->AddEntry((TObject*)0, "", "");
  lxB->AddEntry(h_xB_Truth,"MC truth","pl");
  lxB->AddEntry(h_xB_Reco,"Reco.","pl");
  lxB->AddEntry(h_xB_Corr,"Reco., #varepsilon corr.","pl");
  lxB->Draw();
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cxB->cd(2);
  gPad->SetLogx();
  gPad->SetLogy();
  h_xB_Resp->GetXaxis()->SetTitleSize(0.05);
  h_xB_Resp->GetXaxis()->SetTitleOffset(0.9);
  h_xB_Resp->GetXaxis()->SetLabelSize(0.05);
  h_xB_Resp->GetYaxis()->SetTitleSize(0.05);
  h_xB_Resp->GetYaxis()->SetTitleOffset(0.9);
  h_xB_Resp->GetYaxis()->SetLabelSize(0.05);
  h_xB_Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cxB->cd(3);
  gPad->SetLogx();
  h_xB_Pur->GetXaxis()->SetTitleSize(0.05);
  h_xB_Pur->GetXaxis()->SetTitleOffset(0.9);
  h_xB_Pur->GetXaxis()->SetLabelSize(0.05);
  h_xB_Pur->GetYaxis()->SetTitleSize(0.05);
  h_xB_Pur->GetYaxis()->SetTitleOffset(0.88);
  h_xB_Pur->GetYaxis()->SetLabelSize(0.05);
  h_xB_Pur->Draw("p");
  cxB->cd(4);
  gPad->SetLogx();
  h_dxBvxB->GetXaxis()->SetTitleSize(0.05);
  h_dxBvxB->GetXaxis()->SetTitleOffset(0.9);
  h_dxBvxB->GetXaxis()->SetLabelSize(0.05);
  h_dxBvxB->GetYaxis()->SetTitleSize(0.05);
  h_dxBvxB->GetYaxis()->SetTitleOffset(0.8);
  h_dxBvxB->GetYaxis()->SetLabelSize(0.05);
  h_dxBvxB->Draw("colz");

  cxB->Print("DVCSQA_temp002.pdf");
  cxB->Close();

  TCanvas* cy = new TCanvas("cy","",1200,800);
  cy->Divide(2,2);
  cy->cd(1);
  h_y_Truth->GetXaxis()->SetTitle("y");
  h_y_Truth->GetXaxis()->SetTitleSize(0.05);
  h_y_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_y_Truth->GetXaxis()->SetLabelSize(0.05);
  h_y_Truth->Draw("ep");
  h_y_Reco->Draw("epsame");
  h_y_Corr->Draw("epsame");
  TLegend* ly = new TLegend(0.56,0.50,0.89,0.84);
  ly->SetLineColorAlpha(kWhite,0);
  ly->SetLineWidth(0);
  ly->SetFillStyle(0);
  ly->SetFillColorAlpha(kWhite,0);
  ly->SetHeader("#splitline{y; single e'}{Q^{2}#geq 1 GeV^{2}}");
  ly->AddEntry((TObject*)0, "", "");
  ly->AddEntry(h_y_Truth,"MC truth","pl");
  ly->AddEntry(h_y_Reco,"Reco.","pl");
  ly->AddEntry(h_y_Corr,"Reco., #varepsilon corr.","pl");
  ly->Draw();
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cy->cd(2);
  h_y_Resp->GetXaxis()->SetTitleSize(0.05);
  h_y_Resp->GetXaxis()->SetTitleOffset(0.9);
  h_y_Resp->GetXaxis()->SetLabelSize(0.05);
  h_y_Resp->GetYaxis()->SetTitleSize(0.05);
  h_y_Resp->GetYaxis()->SetTitleOffset(0.9);
  h_y_Resp->GetYaxis()->SetLabelSize(0.05);
  h_y_Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cy->cd(3);
  h_y_Pur->GetXaxis()->SetTitleSize(0.05);
  h_y_Pur->GetXaxis()->SetTitleOffset(0.9);
  h_y_Pur->GetXaxis()->SetLabelSize(0.05);
  h_y_Pur->GetYaxis()->SetTitleSize(0.05);
  h_y_Pur->GetYaxis()->SetTitleOffset(0.88);
  h_y_Pur->GetYaxis()->SetLabelSize(0.05);
  h_y_Pur->Draw("p");
  cy->cd(4);
  h_dyvy->GetXaxis()->SetTitleSize(0.05);
  h_dyvy->GetXaxis()->SetTitleOffset(0.9);
  h_dyvy->GetXaxis()->SetLabelSize(0.05);
  h_dyvy->GetYaxis()->SetTitleSize(0.05);
  h_dyvy->GetYaxis()->SetTitleOffset(0.8);
  h_dyvy->GetYaxis()->SetLabelSize(0.05);
  h_dyvy->Draw("colz");

  cy->Print("DVCSQA_temp003.pdf");
  cy->Close();

  TCanvas* cB0t = new TCanvas("cB0t","",1200,800);
  cB0t->Divide(2,2);
  cB0t->cd(1);
  gPad->SetLogy();
  h_t_Truth->GetXaxis()->SetTitle("|t| [(GeV/#it{c}^{2})^{2}]");
  h_t_Truth->GetXaxis()->SetTitleSize(0.05);
  h_t_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_t_Truth->GetXaxis()->SetLabelSize(0.05);
  h_t_Truth->GetYaxis()->SetLabelSize(0.05);
  h_t_Truth->Draw("ep");
  h_t_B0Reco->Draw("epsame");
  h_t_B0Corr->Draw("epsame");
  TLegend* lB0t = new TLegend(0.50,0.54,0.98,0.85);
  lB0t->SetLineColorAlpha(kWhite,0);
  lB0t->SetLineWidth(0);
  lB0t->SetFillStyle(0);
  lB0t->SetFillColorAlpha(kWhite,0);
  lB0t->SetHeader("#splitline{t; full e'p'#gamma final state}{Q^{2}#geq 1 GeV^{2}; MM^{2}#leq 1 GeV^{2}}");
  lB0t->AddEntry((TObject*)0, "", "");
  lB0t->AddEntry(h_t_Truth,"MC truth","pl");
  lB0t->AddEntry(h_t_B0Reco,"B0 Reco.","pl");
  lB0t->AddEntry(h_t_B0Corr,"B0 Reco., #varepsilon corr.","pl");
  lB0t->Draw();
  TLatex* tB0t = new TLatex(0.52, 0.48, "#theta_{p}(reco.)#in [5.5, 20] mrad");
  tB0t->SetNDC();
  tB0t->SetTextSize(20);
  tB0t->SetTextFont(43);
  tB0t->SetTextColor(kBlack);
  tB0t->Draw("same");
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cB0t->cd(2);
  h_t_B0Resp->GetXaxis()->SetTitleSize(0.05);
  h_t_B0Resp->GetXaxis()->SetTitleOffset(0.9);
  h_t_B0Resp->GetXaxis()->SetLabelSize(0.05);
  h_t_B0Resp->GetYaxis()->SetTitleSize(0.05);
  h_t_B0Resp->GetYaxis()->SetTitleOffset(0.9);
  h_t_B0Resp->GetYaxis()->SetLabelSize(0.04);
  h_t_B0Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cB0t->cd(3);
  h_t_B0Pur->GetXaxis()->SetTitleSize(0.05);
  h_t_B0Pur->GetXaxis()->SetTitleOffset(0.9);
  h_t_B0Pur->GetXaxis()->SetLabelSize(0.05);
  h_t_B0Pur->GetYaxis()->SetTitleSize(0.05);
  h_t_B0Pur->GetYaxis()->SetTitleOffset(0.88);
  h_t_B0Pur->GetYaxis()->SetLabelSize(0.05);
  h_t_B0Pur->Draw("p");
  cB0t->cd(4);
  h_B0dtvt->GetXaxis()->SetTitleSize(0.05);
  h_B0dtvt->GetXaxis()->SetTitleOffset(0.9);
  h_B0dtvt->GetXaxis()->SetLabelSize(0.05);
  h_B0dtvt->GetYaxis()->SetTitleSize(0.05);
  h_B0dtvt->GetYaxis()->SetTitleOffset(0.8);
  h_B0dtvt->GetYaxis()->SetLabelSize(0.05);
  h_B0dtvt->Draw("colz");

  cB0t->Print("DVCSQA_temp004.pdf");
  cB0t->Close();
    
  TCanvas* cRPt = new TCanvas("cRPt","",1200,800);
  cRPt->Divide(2,2);
  cRPt->cd(1);
  gPad->SetLogy();
  h_t_Truth->GetXaxis()->SetTitle("|t| [(GeV/#it{c}^{2})^{2}]");
  h_t_Truth->GetXaxis()->SetTitleSize(0.05);
  h_t_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_t_Truth->GetXaxis()->SetLabelSize(0.05);
  h_t_Truth->Draw("ep");
  h_t_RPReco->Draw("epsame");
  h_t_RPCorr->Draw("epsame");
  TLegend* lRPt = new TLegend(0.50,0.54,0.98,0.85);
  lRPt->SetLineColorAlpha(kWhite,0);
  lRPt->SetLineWidth(0);
  lRPt->SetFillStyle(0);
  lRPt->SetFillColorAlpha(kWhite,0);
  lRPt->SetHeader("#splitline{t; full e'p'#gamma final state}{Q^{2}#geq 1 GeV^{2}; MM^{2}#leq 1 GeV^{2}}");
  lRPt->AddEntry((TObject*)0, "", "");
  lRPt->AddEntry(h_t_Truth,"MC truth","pl");
  lRPt->AddEntry(h_t_RPReco,"RP Reco.","pl");
  lRPt->AddEntry(h_t_RPCorr,"RP Reco., #varepsilon corr.","pl");
  lRPt->Draw();
  TLatex* tRPt = new TLatex(0.52, 0.48, "#theta_{p}(reco.)#in [0, 5.5] mrad");
  tRPt->SetNDC();
  tRPt->SetTextSize(20);
  tRPt->SetTextFont(43);
  tRPt->SetTextColor(kBlack);
  tRPt->Draw("same");
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cRPt->cd(2);
  h_t_RPResp->GetXaxis()->SetTitleSize(0.05);
  h_t_RPResp->GetXaxis()->SetTitleOffset(0.9);
  h_t_RPResp->GetXaxis()->SetLabelSize(0.05);
  h_t_RPResp->GetYaxis()->SetTitleSize(0.05);
  h_t_RPResp->GetYaxis()->SetTitleOffset(0.9);
  h_t_RPResp->GetYaxis()->SetLabelSize(0.05);
  h_t_RPResp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cRPt->cd(3);
  h_t_RPPur->GetXaxis()->SetTitleSize(0.05);
  h_t_RPPur->GetXaxis()->SetTitleOffset(0.9);
  h_t_RPPur->GetXaxis()->SetLabelSize(0.05);
  h_t_RPPur->GetYaxis()->SetTitleSize(0.05);
  h_t_RPPur->GetYaxis()->SetTitleOffset(0.88);
  h_t_RPPur->GetYaxis()->SetLabelSize(0.05);
  h_t_RPPur->Draw("p");
  cRPt->cd(4);
  h_RPdtvt->GetXaxis()->SetTitleSize(0.05);
  h_RPdtvt->GetXaxis()->SetTitleOffset(0.9);
  h_RPdtvt->GetXaxis()->SetLabelSize(0.05);
  h_RPdtvt->GetYaxis()->SetTitleSize(0.05);
  h_RPdtvt->GetYaxis()->SetTitleOffset(0.8);
  h_RPdtvt->GetYaxis()->SetLabelSize(0.05);
  h_RPdtvt->Draw("colz");

  cRPt->Print("DVCSQA_temp005.pdf");
  cRPt->Close();

  TCanvas* cB0theta = new TCanvas("cB0theta","",1200,800);
  cB0theta->Divide(2,2);
  cB0theta->cd(1);
  gPad->SetLogy();
  h_theta_p_Truth->GetXaxis()->SetTitle("#theta_{p} [mrad]");
  h_theta_p_Truth->GetXaxis()->SetTitleSize(0.05);
  h_theta_p_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_theta_p_Truth->GetXaxis()->SetLabelSize(0.05);
  h_theta_p_Truth->GetYaxis()->SetLabelSize(0.05);
  h_theta_p_Truth->Draw("ep");
  h_theta_p_B0Reco->Draw("epsame");
  h_theta_p_B0Corr->Draw("epsame");
  TLegend* lB0theta = new TLegend(0.50,0.63,0.98,0.88);
  lB0theta->SetLineColorAlpha(kWhite,0);
  lB0theta->SetLineWidth(0);
  lB0theta->SetFillStyle(0);
  lB0theta->SetFillColorAlpha(kWhite,0);
  lB0theta->SetHeader("#theta_{p}; single B0 p'");
  lB0theta->AddEntry(h_theta_p_Truth,"MC truth","pl");
  lB0theta->AddEntry(h_theta_p_B0Reco,"B0 Reco.","pl");
  lB0theta->AddEntry(h_theta_p_B0Corr,"B0 Reco., #varepsilon corr.","pl");
  lB0theta->Draw();
  TLatex* tB0theta = new TLatex(0.52, 0.57, "#theta_{p}(reco.)#in [5.5, 20] mrad");
  tB0theta->SetNDC();
  tB0theta->SetTextSize(20);
  tB0theta->SetTextFont(43);
  tB0theta->SetTextColor(kBlack);
  tB0theta->Draw("same");
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cB0theta->cd(2);
  h_theta_p_B0Resp->GetXaxis()->SetTitleSize(0.05);
  h_theta_p_B0Resp->GetXaxis()->SetTitleOffset(0.9);
  h_theta_p_B0Resp->GetXaxis()->SetLabelSize(0.05);
  h_theta_p_B0Resp->GetYaxis()->SetTitleSize(0.05);
  h_theta_p_B0Resp->GetYaxis()->SetTitleOffset(0.9);
  h_theta_p_B0Resp->GetYaxis()->SetLabelSize(0.04);
  h_theta_p_B0Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cB0theta->cd(3);
  h_theta_p_B0Pur->GetXaxis()->SetTitleSize(0.05);
  h_theta_p_B0Pur->GetXaxis()->SetTitleOffset(0.9);
  h_theta_p_B0Pur->GetXaxis()->SetLabelSize(0.05);
  h_theta_p_B0Pur->GetYaxis()->SetTitleSize(0.05);
  h_theta_p_B0Pur->GetYaxis()->SetTitleOffset(0.88);
  h_theta_p_B0Pur->GetYaxis()->SetLabelSize(0.05);
  h_theta_p_B0Pur->Draw("p");

  cB0theta->Print("DVCSQA_temp006.pdf");
  cB0theta->Close();

  TCanvas* cRPtheta = new TCanvas("cRPtheta","",1200,800);
  cRPtheta->Divide(2,2);
  cRPtheta->cd(1);
  gPad->SetLogy();
  h_theta_p_Truth->GetXaxis()->SetTitle("#theta_{p} [mrad]");
  h_theta_p_Truth->GetXaxis()->SetTitleSize(0.05);
  h_theta_p_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_theta_p_Truth->GetXaxis()->SetLabelSize(0.05);
  h_theta_p_Truth->GetYaxis()->SetLabelSize(0.05);
  h_theta_p_Truth->Draw("ep");
  h_theta_p_RPReco->Draw("epsame");
  h_theta_p_RPCorr->Draw("epsame");
  TLegend* lRPtheta = new TLegend(0.48,0.63,0.98,0.88);
  lRPtheta->SetLineColorAlpha(kWhite,0);
  lRPtheta->SetLineWidth(0);
  lRPtheta->SetFillStyle(0);
  lRPtheta->SetFillColorAlpha(kWhite,0);
  lRPtheta->SetHeader("#theta_{p}; single RP p'");
  lRPtheta->AddEntry(h_theta_p_Truth,"MC truth","pl");
  lRPtheta->AddEntry(h_theta_p_RPReco,"RP Reco.","pl");
  lRPtheta->AddEntry(h_theta_p_RPCorr,"RP Reco., #varepsilon corr.","pl");
  lRPtheta->Draw();
  TLatex* tRPtheta = new TLatex(0.50, 0.57, "#theta_{p}(reco.)#in [0, 5.5] mrad");
  tRPtheta->SetNDC();
  tRPtheta->SetTextSize(20);
  tRPtheta->SetTextFont(43);
  tRPtheta->SetTextColor(kBlack);
  tRPtheta->Draw("same");
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cRPtheta->cd(2);
  h_theta_p_RPResp->GetXaxis()->SetTitleSize(0.05);
  h_theta_p_RPResp->GetXaxis()->SetTitleOffset(0.9);
  h_theta_p_RPResp->GetXaxis()->SetLabelSize(0.05);
  h_theta_p_RPResp->GetYaxis()->SetTitleSize(0.05);
  h_theta_p_RPResp->GetYaxis()->SetTitleOffset(0.9);
  h_theta_p_RPResp->GetYaxis()->SetLabelSize(0.04);
  h_theta_p_RPResp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cRPtheta->cd(3);
  h_theta_p_RPPur->GetXaxis()->SetTitleSize(0.05);
  h_theta_p_RPPur->GetXaxis()->SetTitleOffset(0.9);
  h_theta_p_RPPur->GetXaxis()->SetLabelSize(0.05);
  h_theta_p_RPPur->GetYaxis()->SetTitleSize(0.05);
  h_theta_p_RPPur->GetYaxis()->SetTitleOffset(0.88);
  h_theta_p_RPPur->GetYaxis()->SetLabelSize(0.05);
  h_theta_p_RPPur->Draw("p");

  cRPtheta->Print("DVCSQA_temp007.pdf");
  cRPtheta->Close();

  TCanvas* cB0E = new TCanvas("cB0E","",1200,800);
  cB0E->Divide(2,2);
  cB0E->cd(1);
  gPad->SetLogy();
  h_E_p_Truth->GetXaxis()->SetTitle("E_{p} [GeV]");
  h_E_p_Truth->GetXaxis()->SetTitleSize(0.05);
  h_E_p_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_E_p_Truth->GetXaxis()->SetLabelSize(0.05);
  h_E_p_Truth->GetYaxis()->SetLabelSize(0.05);
  h_E_p_Truth->Draw("ep");
  h_E_p_B0Reco->Draw("epsame");
  h_E_p_B0Corr->Draw("epsame");
  TLegend* lB0E = new TLegend(0.51,0.63,0.91,0.88);
  lB0E->SetLineColorAlpha(kWhite,0);
  lB0E->SetLineWidth(0);
  lB0E->SetFillStyle(0);
  lB0E->SetFillColorAlpha(kWhite,0);
  lB0E->SetHeader("E_{p}; single B0 p'");
  lB0E->AddEntry(h_E_p_Truth,"MC truth","pl");
  lB0E->AddEntry(h_E_p_B0Reco,"B0 Reco.","pl");
  lB0E->AddEntry(h_E_p_B0Corr,"B0 Reco., #varepsilon corr.","pl");
  lB0E->Draw();
  TLatex* tB0E = new TLatex(0.52, 0.57, "#theta_{p}(reco.)#in [5.5, 20] mrad");
  tB0E->SetNDC();
  tB0E->SetTextSize(20);
  tB0E->SetTextFont(43);
  tB0E->SetTextColor(kBlack);
  tB0E->Draw("same");
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cB0E->cd(2);
  h_E_p_B0Resp->GetXaxis()->SetTitleSize(0.05);
  h_E_p_B0Resp->GetXaxis()->SetTitleOffset(0.9);
  h_E_p_B0Resp->GetXaxis()->SetLabelSize(0.05);
  h_E_p_B0Resp->GetYaxis()->SetTitleSize(0.05);
  h_E_p_B0Resp->GetYaxis()->SetTitleOffset(0.9);
  h_E_p_B0Resp->GetYaxis()->SetLabelSize(0.04);
  h_E_p_B0Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cB0E->cd(3);
  h_E_p_B0Pur->GetXaxis()->SetTitleSize(0.05);
  h_E_p_B0Pur->GetXaxis()->SetTitleOffset(0.9);
  h_E_p_B0Pur->GetXaxis()->SetLabelSize(0.05);
  h_E_p_B0Pur->GetYaxis()->SetTitleSize(0.05);
  h_E_p_B0Pur->GetYaxis()->SetTitleOffset(0.88);
  h_E_p_B0Pur->GetYaxis()->SetLabelSize(0.05);
  h_E_p_B0Pur->Draw("p");

  cB0E->Print("DVCSQA_temp008.pdf");
  cB0E->Close();

  TCanvas* cRPE = new TCanvas("cRPE","",1200,800);
  cRPE->Divide(2,2);
  cRPE->cd(1);
  gPad->SetLogy();
  h_E_p_Truth->GetXaxis()->SetTitle("E_{p} [GeV]");
  h_E_p_Truth->GetXaxis()->SetTitleSize(0.05);
  h_E_p_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_E_p_Truth->GetXaxis()->SetLabelSize(0.05);
  h_E_p_Truth->GetYaxis()->SetLabelSize(0.05);
  h_E_p_Truth->Draw("ep");
  h_E_p_RPReco->Draw("epsame");
  h_E_p_RPCorr->Draw("epsame");
  TLegend* lRPE = new TLegend(0.51,0.63,0.91,0.88);
  lRPE->SetLineColorAlpha(kWhite,0);
  lRPE->SetLineWidth(0);
  lRPE->SetFillStyle(0);
  lRPE->SetFillColorAlpha(kWhite,0);
  lRPE->SetHeader("E_{p}; single RP p'");
  lRPE->AddEntry(h_E_p_Truth,"MC truth","pl");
  lRPE->AddEntry(h_E_p_RPReco,"RP Reco.","pl");
  lRPE->AddEntry(h_E_p_RPCorr,"RP Reco., #varepsilon corr.","pl");
  lRPE->Draw();
  TLatex* tRPE = new TLatex(0.52, 0.57, "#theta_{p}(reco.)#in [0, 5.5] mrad");
  tRPE->SetNDC();
  tRPE->SetTextSize(20);
  tRPE->SetTextFont(43);
  tRPE->SetTextColor(kBlack);
  tRPE->Draw("same");
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cRPE->cd(2);
  h_E_p_RPResp->GetXaxis()->SetTitleSize(0.05);
  h_E_p_RPResp->GetXaxis()->SetTitleOffset(0.9);
  h_E_p_RPResp->GetXaxis()->SetLabelSize(0.05);
  h_E_p_RPResp->GetYaxis()->SetTitleSize(0.05);
  h_E_p_RPResp->GetYaxis()->SetTitleOffset(0.9);
  h_E_p_RPResp->GetYaxis()->SetLabelSize(0.04);
  h_E_p_RPResp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cRPE->cd(3);
  h_E_p_RPPur->GetXaxis()->SetTitleSize(0.05);
  h_E_p_RPPur->GetXaxis()->SetTitleOffset(0.9);
  h_E_p_RPPur->GetXaxis()->SetLabelSize(0.05);
  h_E_p_RPPur->GetYaxis()->SetTitleSize(0.05);
  h_E_p_RPPur->GetYaxis()->SetTitleOffset(0.88);
  h_E_p_RPPur->GetYaxis()->SetLabelSize(0.05);
  h_E_p_RPPur->Draw("p");

  cRPE->Print("DVCSQA_temp009.pdf");
  cRPE->Close();

  TCanvas* cetheta = new TCanvas("cetheta","",1200,800);
  cetheta->Divide(2,2);
  cetheta->cd(1);
  gPad->SetLogy();
  h_theta_e_Truth->GetXaxis()->SetTitle("#theta_{e'} [rad]");
  h_theta_e_Truth->GetXaxis()->SetTitleSize(0.05);
  h_theta_e_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_theta_e_Truth->GetXaxis()->SetLabelSize(0.05);
  h_theta_e_Truth->GetYaxis()->SetLabelSize(0.05);
  h_theta_e_Truth->Draw("ep");
  h_theta_e_Reco->Draw("epsame");
  h_theta_e_Corr->Draw("epsame");
  TLegend* letheta = new TLegend(0.13,0.50,0.61,0.84);
  letheta->SetLineColorAlpha(kWhite,0);
  letheta->SetLineWidth(0);
  letheta->SetFillStyle(0);
  letheta->SetFillColorAlpha(kWhite,0);
  letheta->SetHeader("#splitline{#theta_{e}; single e'}{Q^{2}#geq 1 GeV^{2}}");
  letheta->AddEntry((TObject*)0, "", "");
  letheta->AddEntry(h_theta_e_Truth,"MC truth","pl");
  letheta->AddEntry(h_theta_e_Reco,"Reco.","pl");
  letheta->AddEntry(h_theta_e_Corr,"Reco., #varepsilon corr.","pl");
  letheta->Draw();
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cetheta->cd(2);
  h_theta_e_Resp->GetXaxis()->SetTitleSize(0.05);
  h_theta_e_Resp->GetXaxis()->SetTitleOffset(0.9);
  h_theta_e_Resp->GetXaxis()->SetLabelSize(0.05);
  h_theta_e_Resp->GetYaxis()->SetTitleSize(0.05);
  h_theta_e_Resp->GetYaxis()->SetTitleOffset(0.9);
  h_theta_e_Resp->GetYaxis()->SetLabelSize(0.04);
  h_theta_e_Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cetheta->cd(3);
  h_theta_e_Pur->GetXaxis()->SetTitleSize(0.05);
  h_theta_e_Pur->GetXaxis()->SetTitleOffset(0.9);
  h_theta_e_Pur->GetXaxis()->SetLabelSize(0.05);
  h_theta_e_Pur->GetYaxis()->SetTitleSize(0.05);
  h_theta_e_Pur->GetYaxis()->SetTitleOffset(0.88);
  h_theta_e_Pur->GetYaxis()->SetLabelSize(0.05);
  h_theta_e_Pur->Draw("p");

  cetheta->Print("DVCSQA_temp010.pdf");
  cetheta->Close();

  TCanvas* ceE = new TCanvas("ceE","",1200,800);
  ceE->Divide(2,2);
  ceE->cd(1);
  gPad->SetLogy();
  h_E_e_Truth->GetXaxis()->SetTitle("E_{e} [GeV]");
  h_E_e_Truth->GetXaxis()->SetTitleSize(0.05);
  h_E_e_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_E_e_Truth->GetXaxis()->SetLabelSize(0.05);
  h_E_e_Truth->GetYaxis()->SetLabelSize(0.05);
  h_E_e_Truth->Draw("ep");
  h_E_e_Reco->Draw("epsame");
  h_E_e_Corr->Draw("epsame");
  TLegend* leE = new TLegend(0.52,0.5,0.99,0.83);
  leE->SetLineColorAlpha(kWhite,0);
  leE->SetLineWidth(0);
  leE->SetFillStyle(0);
  leE->SetFillColorAlpha(kWhite,0);
  leE->SetHeader("#splitline{E_{e}; single e'}{Q^{2}#geq 1 GeV^{2}}");
  leE->AddEntry((TObject*)0, "", "");
  leE->AddEntry(h_E_e_Truth,"MC truth","pl");
  leE->AddEntry(h_E_e_Reco,"Reco.","pl");
  leE->AddEntry(h_E_e_Corr,"Reco., #varepsilon corr.","pl");
  leE->Draw();
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  ceE->cd(2);
  h_E_e_Resp->GetXaxis()->SetTitleSize(0.05);
  h_E_e_Resp->GetXaxis()->SetTitleOffset(0.9);
  h_E_e_Resp->GetXaxis()->SetLabelSize(0.05);
  h_E_e_Resp->GetYaxis()->SetTitleSize(0.05);
  h_E_e_Resp->GetYaxis()->SetTitleOffset(0.9);
  h_E_e_Resp->GetYaxis()->SetLabelSize(0.04);
  h_E_e_Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  ceE->cd(3);
  h_E_e_Pur->GetXaxis()->SetTitleSize(0.05);
  h_E_e_Pur->GetXaxis()->SetTitleOffset(0.9);
  h_E_e_Pur->GetXaxis()->SetLabelSize(0.05);
  h_E_e_Pur->GetYaxis()->SetTitleSize(0.05);
  h_E_e_Pur->GetYaxis()->SetTitleOffset(0.88);
  h_E_e_Pur->GetYaxis()->SetLabelSize(0.05);
  h_E_e_Pur->Draw("p");

  ceE->Print("DVCSQA_temp011.pdf");
  ceE->Close();

  TCanvas* cgtheta = new TCanvas("cgtheta","",1200,800);
  cgtheta->Divide(2,2);
  cgtheta->cd(1);
  h_theta_g_Truth->GetXaxis()->SetTitle("#theta_{#gamma} [rad]");
  h_theta_g_Truth->GetXaxis()->SetTitleSize(0.05);
  h_theta_g_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_theta_g_Truth->GetXaxis()->SetLabelSize(0.05);
  h_theta_g_Truth->GetYaxis()->SetLabelSize(0.05);
  h_theta_g_Truth->Draw("ep");
  h_theta_g_Reco->Draw("epsame");
  h_theta_g_Corr->Draw("epsame");
  TLegend* lgtheta = new TLegend(0.14,0.63,0.61,0.88);
  lgtheta->SetLineColorAlpha(kWhite,0);
  lgtheta->SetLineWidth(0);
  lgtheta->SetFillStyle(0);
  lgtheta->SetFillColorAlpha(kWhite,0);
  lgtheta->SetHeader("#theta_{#gamma}; single #gamma");
  lgtheta->AddEntry(h_theta_g_Truth,"MC truth","pl");
  lgtheta->AddEntry(h_theta_g_Reco,"Reco.","pl");
  lgtheta->AddEntry(h_theta_g_Corr,"Reco., #varepsilon corr.","pl");
  lgtheta->Draw();
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cgtheta->cd(2);
  h_theta_g_Resp->GetXaxis()->SetTitleSize(0.05);
  h_theta_g_Resp->GetXaxis()->SetTitleOffset(0.9);
  h_theta_g_Resp->GetXaxis()->SetLabelSize(0.05);
  h_theta_g_Resp->GetYaxis()->SetTitleSize(0.05);
  h_theta_g_Resp->GetYaxis()->SetTitleOffset(0.9);
  h_theta_g_Resp->GetYaxis()->SetLabelSize(0.04);
  h_theta_g_Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cgtheta->cd(3);
  h_theta_g_Pur->GetXaxis()->SetTitleSize(0.05);
  h_theta_g_Pur->GetXaxis()->SetTitleOffset(0.9);
  h_theta_g_Pur->GetXaxis()->SetLabelSize(0.05);
  h_theta_g_Pur->GetYaxis()->SetTitleSize(0.05);
  h_theta_g_Pur->GetYaxis()->SetTitleOffset(0.88);
  h_theta_g_Pur->GetYaxis()->SetLabelSize(0.05);
  h_theta_g_Pur->Draw("p");

  cgtheta->Print("DVCSQA_temp012.pdf");
  cgtheta->Close();

  TCanvas* cgE = new TCanvas("cgE","",1200,800);
  cgE->Divide(2,2);
  cgE->cd(1);
  gPad->SetLogy();
  h_E_g_Truth->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
  h_E_g_Truth->GetXaxis()->SetTitleSize(0.05);
  h_E_g_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_E_g_Truth->GetXaxis()->SetLabelSize(0.05);
  h_E_g_Truth->GetYaxis()->SetLabelSize(0.05);
  h_E_g_Truth->Draw("ep");
  h_E_g_Reco->Draw("epsame");
  h_E_g_Corr->Draw("epsame");
  TLegend* lgE = new TLegend(0.53,0.63,0.99,0.88);
  lgE->SetLineColorAlpha(kWhite,0);
  lgE->SetLineWidth(0);
  lgE->SetFillStyle(0);
  lgE->SetFillColorAlpha(kWhite,0);
  lgE->SetHeader("E_{#gamma}; single #gamma");
  lgE->AddEntry(h_E_g_Truth,"MC truth","pl");
  lgE->AddEntry(h_E_g_Reco,"Reco.","pl");
  lgE->AddEntry(h_E_g_Corr,"Reco., #varepsilon corr.","pl");
  lgE->Draw();
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cgE->cd(2);
  h_E_g_Resp->GetXaxis()->SetTitleSize(0.05);
  h_E_g_Resp->GetXaxis()->SetTitleOffset(0.9);
  h_E_g_Resp->GetXaxis()->SetLabelSize(0.05);
  h_E_g_Resp->GetYaxis()->SetTitleSize(0.05);
  h_E_g_Resp->GetYaxis()->SetTitleOffset(0.9);
  h_E_g_Resp->GetYaxis()->SetLabelSize(0.04);
  h_E_g_Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cgE->cd(3);
  h_E_g_Pur->GetXaxis()->SetTitleSize(0.05);
  h_E_g_Pur->GetXaxis()->SetTitleOffset(0.9);
  h_E_g_Pur->GetXaxis()->SetLabelSize(0.05);
  h_E_g_Pur->GetYaxis()->SetTitleSize(0.05);
  h_E_g_Pur->GetYaxis()->SetTitleOffset(0.88);
  h_E_g_Pur->GetYaxis()->SetLabelSize(0.05);
  h_E_g_Pur->Draw("p");

  cgE->Print("DVCSQA_temp013.pdf");
  cgE->Close();

  TCanvas* cHFSSigma = new TCanvas("cHFSSigma","",1200,800);
  cHFSSigma->Divide(2,2);
  cHFSSigma->cd(1);
  gPad->SetLogy();
  h_HFSSigma_Truth->GetXaxis()->SetTitle("#Sigma_{h} [GeV]");
  h_HFSSigma_Truth->GetXaxis()->SetTitleSize(0.05);
  h_HFSSigma_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_HFSSigma_Truth->GetXaxis()->SetLabelSize(0.05);
  h_HFSSigma_Truth->GetYaxis()->SetLabelSize(0.05);
  h_HFSSigma_Truth->GetYaxis()->SetRangeUser(1e0,1e7);
  h_HFSSigma_Truth->Draw("ep");
  h_HFSSigma_Reco->Draw("epsame");
  h_HFSSigma_Corr->Draw("epsame");
  TLegend* lHFSSigma = new TLegend(0.13,0.13,0.59,0.45);
  lHFSSigma->SetLineColorAlpha(kWhite,0);
  lHFSSigma->SetLineWidth(0);
  lHFSSigma->SetFillStyle(0);
  lHFSSigma->SetFillColorAlpha(kWhite,0);
  lHFSSigma->SetHeader("#Sigma_{h}; p'+#gamma reco.");
  lHFSSigma->AddEntry(h_HFSSigma_Truth,"MC truth","pl");
  lHFSSigma->AddEntry(h_HFSSigma_Reco,"Reco.","pl");
  lHFSSigma->AddEntry(h_HFSSigma_Corr,"Reco., #varepsilon corr.","pl");
  lHFSSigma->Draw();
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cHFSSigma->cd(2);
  h_HFSSigma_Resp->GetXaxis()->SetTitleSize(0.05);
  h_HFSSigma_Resp->GetXaxis()->SetTitleOffset(0.9);
  h_HFSSigma_Resp->GetXaxis()->SetLabelSize(0.05);
  h_HFSSigma_Resp->GetYaxis()->SetTitleSize(0.05);
  h_HFSSigma_Resp->GetYaxis()->SetTitleOffset(0.9);
  h_HFSSigma_Resp->GetYaxis()->SetLabelSize(0.04);
  h_HFSSigma_Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cHFSSigma->cd(3);
  h_HFSSigma_Pur->GetXaxis()->SetTitleSize(0.05);
  h_HFSSigma_Pur->GetXaxis()->SetTitleOffset(0.9);
  h_HFSSigma_Pur->GetXaxis()->SetLabelSize(0.05);
  h_HFSSigma_Pur->GetYaxis()->SetTitleSize(0.05);
  h_HFSSigma_Pur->GetYaxis()->SetTitleOffset(0.88);
  h_HFSSigma_Pur->GetYaxis()->SetLabelSize(0.04);
  h_HFSSigma_Pur->Draw("p");

  cHFSSigma->Print("DVCSQA_temp014.pdf");
  cHFSSigma->Close();

  TCanvas* cHFSpT2 = new TCanvas("cHFSpT2","",1200,800);
  cHFSpT2->Divide(2,2);
  cHFSpT2->cd(1);
  gPad->SetLogy();
  h_HFSpT2_Truth->GetXaxis()->SetTitle("p^{2}_{T,h} [(GeV/#it{c})^{2}]");
  h_HFSpT2_Truth->GetXaxis()->SetTitleSize(0.04);
  h_HFSpT2_Truth->GetXaxis()->SetTitleOffset(1.0);
  h_HFSpT2_Truth->GetXaxis()->SetLabelSize(0.04);
  h_HFSpT2_Truth->GetYaxis()->SetLabelSize(0.05);
  h_HFSpT2_Truth->GetYaxis()->SetRangeUser(1e0,1e7);
  h_HFSpT2_Truth->Draw("ep");
  h_HFSpT2_Reco->Draw("epsame");
  h_HFSpT2_Corr->Draw("epsame");
  TLegend* lHFSpT2 = new TLegend(0.12,0.12,0.58,0.41);
  lHFSpT2->SetLineColorAlpha(kWhite,0);
  lHFSpT2->SetLineWidth(0);
  lHFSpT2->SetFillStyle(0);
  lHFSpT2->SetFillColorAlpha(kWhite,0);
  lHFSpT2->SetHeader("p^{2}_{T,h}; p'+#gamma reco.");
  lHFSpT2->AddEntry(h_HFSpT2_Truth,"MC truth","pl");
  lHFSpT2->AddEntry(h_HFSpT2_Reco,"Reco.","pl");
  lHFSpT2->AddEntry(h_HFSpT2_Corr,"Reco., #varepsilon corr.","pl");
  lHFSpT2->Draw();
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cHFSpT2->cd(2);
  h_HFSpT2_Resp->GetXaxis()->SetTitleSize(0.04);
  h_HFSpT2_Resp->GetXaxis()->SetTitleOffset(1.0);
  h_HFSpT2_Resp->GetXaxis()->SetLabelSize(0.04);
  h_HFSpT2_Resp->GetYaxis()->SetTitleSize(0.05);
  h_HFSpT2_Resp->GetYaxis()->SetTitleOffset(0.9);
  h_HFSpT2_Resp->GetYaxis()->SetLabelSize(0.05);
  h_HFSpT2_Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cHFSpT2->cd(3);
  h_HFSpT2_Pur->GetXaxis()->SetTitleSize(0.04);
  h_HFSpT2_Pur->GetXaxis()->SetTitleOffset(1.0);
  h_HFSpT2_Pur->GetXaxis()->SetLabelSize(0.04);
  h_HFSpT2_Pur->GetYaxis()->SetTitleSize(0.05);
  h_HFSpT2_Pur->GetYaxis()->SetTitleOffset(0.88);
  h_HFSpT2_Pur->GetYaxis()->SetLabelSize(0.04);
  h_HFSpT2_Pur->Draw("p");

  cHFSpT2->Print("DVCSQA_temp015.pdf");
  cHFSpT2->Close();

  TCanvas* cFullSigma = new TCanvas("cFullSigma","",1200,800);
  cFullSigma->Divide(2,2);
  cFullSigma->cd(1);
  gPad->SetLogy();
  h_FullSigma_Truth->GetXaxis()->SetTitle("#Sigma [GeV]");
  h_FullSigma_Truth->GetXaxis()->SetTitleSize(0.05);
  h_FullSigma_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_FullSigma_Truth->GetXaxis()->SetLabelSize(0.05);
  h_FullSigma_Truth->GetYaxis()->SetLabelSize(0.05);
  h_FullSigma_Truth->GetYaxis()->SetRangeUser(1e0,1e7);
  h_FullSigma_Truth->Draw("ep");
  h_FullSigma_Reco->Draw("epsame");
  h_FullSigma_Corr->Draw("epsame");
  TLegend* lFullSigma = new TLegend(0.13,0.51,0.59,0.84);
  lFullSigma->SetLineColorAlpha(kWhite,0);
  lFullSigma->SetLineWidth(0);
  lFullSigma->SetFillStyle(0);
  lFullSigma->SetFillColorAlpha(kWhite,0);
  lFullSigma->SetHeader("#splitline{#Sigma; full e'p'#gamma final state}{Q^{2}#geq 1 GeV^{2}}");
  lFullSigma->AddEntry((TObject*)0, "", "");
  lFullSigma->AddEntry(h_FullSigma_Truth,"MC truth","pl");
  lFullSigma->AddEntry(h_FullSigma_Reco,"Reco.","pl");
  lFullSigma->AddEntry(h_FullSigma_Corr,"Reco., #varepsilon corr.","pl");
  lFullSigma->Draw();
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cFullSigma->cd(2);
  h_FullSigma_Resp->GetXaxis()->SetTitleSize(0.05);
  h_FullSigma_Resp->GetXaxis()->SetTitleOffset(0.9);
  h_FullSigma_Resp->GetXaxis()->SetLabelSize(0.05);
  h_FullSigma_Resp->GetYaxis()->SetTitleSize(0.05);
  h_FullSigma_Resp->GetYaxis()->SetTitleOffset(0.9);
  h_FullSigma_Resp->GetYaxis()->SetLabelSize(0.04);
  h_FullSigma_Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cFullSigma->cd(3);
  h_FullSigma_Pur->GetXaxis()->SetTitleSize(0.05);
  h_FullSigma_Pur->GetXaxis()->SetTitleOffset(0.9);
  h_FullSigma_Pur->GetXaxis()->SetLabelSize(0.05);
  h_FullSigma_Pur->GetYaxis()->SetTitleSize(0.05);
  h_FullSigma_Pur->GetYaxis()->SetTitleOffset(0.88);
  h_FullSigma_Pur->GetYaxis()->SetLabelSize(0.04);
  h_FullSigma_Pur->Draw("p");

  cFullSigma->Print("DVCSQA_temp016.pdf");
  cFullSigma->Close();

  TCanvas* ceEOverp = new TCanvas("ceEOverp","",1200,800);
  ceEOverp->Divide(2,1);
  ceEOverp->cd(1);
  h_EOverp_e_Truth->GetXaxis()->SetTitle("E_{e}/p_{e}");
  h_EOverp_e_Truth->GetXaxis()->SetTitleSize(0.05);
  h_EOverp_e_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_EOverp_e_Truth->GetXaxis()->SetLabelSize(0.05);
  h_EOverp_e_Truth->GetYaxis()->SetLabelSize(0.05);
  h_EOverp_e_Truth->Draw("ep");
  h_EOverp_e_Reco->Draw("epsame");
  TLegend* leEOverp = new TLegend(0.53,0.63,0.99,0.88);
  leEOverp->SetLineColorAlpha(kWhite,0);
  leEOverp->SetLineWidth(0);
  leEOverp->SetFillStyle(0);
  leEOverp->SetFillColorAlpha(kWhite,0);
  leEOverp->SetHeader("E_{e}/p_{e}; single e'");
  leEOverp->AddEntry(h_EOverp_e_Truth,"MC truth","pl");
  leEOverp->AddEntry(h_EOverp_e_Reco,"Reco.","pl");
  leEOverp->Draw();
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  ceEOverp->cd(2);
  h_EOverp_e_Resp->GetXaxis()->SetTitleSize(0.05);
  h_EOverp_e_Resp->GetXaxis()->SetTitleOffset(0.9);
  h_EOverp_e_Resp->GetXaxis()->SetLabelSize(0.05);
  h_EOverp_e_Resp->GetYaxis()->SetTitleSize(0.05);
  h_EOverp_e_Resp->GetYaxis()->SetTitleOffset(0.9);
  h_EOverp_e_Resp->GetYaxis()->SetLabelSize(0.04);
  h_EOverp_e_Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS

  ceEOverp->Print("DVCSQA_temp017.pdf");
  ceEOverp->Close();

  TCanvas* cgEOverp = new TCanvas("cgEOverp","",1200,800);
  cgEOverp->Divide(2,1);
  cgEOverp->cd(1);
  h_EOverp_g_Truth->GetXaxis()->SetTitle("E_{#gamma}/p_{#gamma}");
  h_EOverp_g_Truth->GetXaxis()->SetTitleSize(0.05);
  h_EOverp_g_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_EOverp_g_Truth->GetXaxis()->SetLabelSize(0.05);
  h_EOverp_g_Truth->GetYaxis()->SetLabelSize(0.05);
  h_EOverp_g_Truth->Draw("ep");
  h_EOverp_g_Reco->Draw("epsame");
  TLegend* lgEOverp = new TLegend(0.53,0.63,0.99,0.88);
  lgEOverp->SetLineColorAlpha(kWhite,0);
  lgEOverp->SetLineWidth(0);
  lgEOverp->SetFillStyle(0);
  lgEOverp->SetFillColorAlpha(kWhite,0);
  lgEOverp->SetHeader("E_{#gamma}/p_{#gamma}; single e'");
  lgEOverp->AddEntry(h_EOverp_g_Truth,"MC truth","pl");
  lgEOverp->AddEntry(h_EOverp_g_Reco,"Reco.","pl");
  lgEOverp->Draw();
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cgEOverp->cd(2);
  h_EOverp_g_Resp->GetXaxis()->SetTitleSize(0.05);
  h_EOverp_g_Resp->GetXaxis()->SetTitleOffset(0.9);
  h_EOverp_g_Resp->GetXaxis()->SetLabelSize(0.05);
  h_EOverp_g_Resp->GetYaxis()->SetTitleSize(0.05);
  h_EOverp_g_Resp->GetYaxis()->SetTitleOffset(0.9);
  h_EOverp_g_Resp->GetYaxis()->SetLabelSize(0.04);
  h_EOverp_g_Resp->Draw("colz");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS

  cgEOverp->Print("DVCSQA_temp018.pdf");
  cgEOverp->Close();

  TCanvas* cMig1 = new TCanvas("cMig1","",1200,800);
  cMig1->Divide(2,2);
  cMig1->cd(1);
  h_Q2_Mig->GetXaxis()->SetTitleSize(0.04);
  h_Q2_Mig->GetXaxis()->SetTitleOffset(1.0);
  h_Q2_Mig->GetXaxis()->SetLabelSize(0.04);
  h_Q2_Mig->GetYaxis()->SetTitleSize(0.05);
  h_Q2_Mig->GetYaxis()->SetTitleOffset(0.88);
  h_Q2_Mig->GetYaxis()->SetLabelSize(0.04);
  h_Q2_Mig->Draw("p");
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cMig1->cd(2);
  gPad->SetLogx();
  h_xB_Mig->GetXaxis()->SetTitleSize(0.04);
  h_xB_Mig->GetXaxis()->SetTitleOffset(1.0);
  h_xB_Mig->GetXaxis()->SetLabelSize(0.04);
  h_xB_Mig->GetYaxis()->SetTitleSize(0.05);
  h_xB_Mig->GetYaxis()->SetTitleOffset(0.88);
  h_xB_Mig->GetYaxis()->SetLabelSize(0.04);
  h_xB_Mig->Draw("p");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cMig1->cd(3);
  h_y_Mig->GetXaxis()->SetTitleSize(0.04);
  h_y_Mig->GetXaxis()->SetTitleOffset(1.0);
  h_y_Mig->GetXaxis()->SetLabelSize(0.04);
  h_y_Mig->GetYaxis()->SetTitleSize(0.05);
  h_y_Mig->GetYaxis()->SetTitleOffset(0.88);
  h_y_Mig->GetYaxis()->SetLabelSize(0.04);
  h_y_Mig->Draw("p");

  cMig1->Print("DVCSQA_temp019.pdf");
  cMig1->Close();

  TCanvas* cMig2 = new TCanvas("cMig2","",1200,800);
  cMig2->Divide(2,2);
  cMig2->cd(1);
  h_t_B0Mig->GetXaxis()->SetTitleSize(0.04);
  h_t_B0Mig->GetXaxis()->SetTitleOffset(1.0);
  h_t_B0Mig->GetXaxis()->SetLabelSize(0.04);
  h_t_B0Mig->GetYaxis()->SetTitleSize(0.05);
  h_t_B0Mig->GetYaxis()->SetTitleOffset(0.88);
  h_t_B0Mig->GetYaxis()->SetLabelSize(0.04);
  h_t_B0Mig->Draw("p");
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cMig2->cd(2);
  h_t_RPMig->GetXaxis()->SetTitleSize(0.04);
  h_t_RPMig->GetXaxis()->SetTitleOffset(1.0);
  h_t_RPMig->GetXaxis()->SetLabelSize(0.04);
  h_t_RPMig->GetYaxis()->SetTitleSize(0.05);
  h_t_RPMig->GetYaxis()->SetTitleOffset(0.88);
  h_t_RPMig->GetYaxis()->SetLabelSize(0.04);
  h_t_RPMig->Draw("p");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS

  cMig2->Print("DVCSQA_temp020.pdf");
  cMig2->Close();

  TCanvas* cMig3 = new TCanvas("cMig3","",1200,800);
  cMig3->Divide(2,2);
  cMig3->cd(1);
  h_theta_p_B0Mig->GetXaxis()->SetTitleSize(0.05);
  h_theta_p_B0Mig->GetXaxis()->SetTitleOffset(0.9);
  h_theta_p_B0Mig->GetXaxis()->SetLabelSize(0.05);
  h_theta_p_B0Mig->GetYaxis()->SetTitleSize(0.05);
  h_theta_p_B0Mig->GetYaxis()->SetTitleOffset(0.88);
  h_theta_p_B0Mig->GetYaxis()->SetLabelSize(0.05);
  h_theta_p_B0Mig->Draw("p");
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cMig3->cd(2);
  h_theta_p_RPMig->GetXaxis()->SetTitleSize(0.05);
  h_theta_p_RPMig->GetXaxis()->SetTitleOffset(0.9);
  h_theta_p_RPMig->GetXaxis()->SetLabelSize(0.05);
  h_theta_p_RPMig->GetYaxis()->SetTitleSize(0.05);
  h_theta_p_RPMig->GetYaxis()->SetTitleOffset(0.88);
  h_theta_p_RPMig->GetYaxis()->SetLabelSize(0.05);
  h_theta_p_RPMig->Draw("p");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cMig3->cd(3);
  h_E_p_B0Mig->GetXaxis()->SetTitleSize(0.05);
  h_E_p_B0Mig->GetXaxis()->SetTitleOffset(0.9);
  h_E_p_B0Mig->GetXaxis()->SetLabelSize(0.05);
  h_E_p_B0Mig->GetYaxis()->SetTitleSize(0.05);
  h_E_p_B0Mig->GetYaxis()->SetTitleOffset(0.88);
  h_E_p_B0Mig->GetYaxis()->SetLabelSize(0.05);
  h_E_p_B0Mig->Draw("p");
  cMig3->cd(4);
  h_E_p_RPMig->GetXaxis()->SetTitleSize(0.05);
  h_E_p_RPMig->GetXaxis()->SetTitleOffset(0.9);
  h_E_p_RPMig->GetXaxis()->SetLabelSize(0.05);
  h_E_p_RPMig->GetYaxis()->SetTitleSize(0.05);
  h_E_p_RPMig->GetYaxis()->SetTitleOffset(0.88);
  h_E_p_RPMig->GetYaxis()->SetLabelSize(0.05);
  h_E_p_RPMig->Draw("p");

  cMig3->Print("DVCSQA_temp021.pdf");
  cMig3->Close();

  TCanvas* cMig4 = new TCanvas("cMig4","",1200,800);
  cMig4->Divide(2,2);
  cMig4->cd(1);
  h_theta_e_Mig->GetXaxis()->SetTitleSize(0.05);
  h_theta_e_Mig->GetXaxis()->SetTitleOffset(0.9);
  h_theta_e_Mig->GetXaxis()->SetLabelSize(0.05);
  h_theta_e_Mig->GetYaxis()->SetTitleSize(0.05);
  h_theta_e_Mig->GetYaxis()->SetTitleOffset(0.88);
  h_theta_e_Mig->GetYaxis()->SetLabelSize(0.05);
  h_theta_e_Mig->Draw("p");
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cMig4->cd(2);
  h_theta_g_Mig->GetXaxis()->SetTitleSize(0.05);
  h_theta_g_Mig->GetXaxis()->SetTitleOffset(0.9);
  h_theta_g_Mig->GetXaxis()->SetLabelSize(0.05);
  h_theta_g_Mig->GetYaxis()->SetTitleSize(0.05);
  h_theta_g_Mig->GetYaxis()->SetTitleOffset(0.88);
  h_theta_g_Mig->GetYaxis()->SetLabelSize(0.05);
  h_theta_g_Mig->Draw("p");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cMig4->cd(3);
  h_E_e_Mig->GetXaxis()->SetTitleSize(0.05);
  h_E_e_Mig->GetXaxis()->SetTitleOffset(0.9);
  h_E_e_Mig->GetXaxis()->SetLabelSize(0.05);
  h_E_e_Mig->GetYaxis()->SetTitleSize(0.05);
  h_E_e_Mig->GetYaxis()->SetTitleOffset(0.88);
  h_E_e_Mig->GetYaxis()->SetLabelSize(0.05);
  h_E_e_Mig->Draw("p");
  cMig4->cd(4);
  h_E_g_Mig->GetXaxis()->SetTitleSize(0.05);
  h_E_g_Mig->GetXaxis()->SetTitleOffset(0.9);
  h_E_g_Mig->GetXaxis()->SetLabelSize(0.05);
  h_E_g_Mig->GetYaxis()->SetTitleSize(0.05);
  h_E_g_Mig->GetYaxis()->SetTitleOffset(0.88);
  h_E_g_Mig->GetYaxis()->SetLabelSize(0.05);
  h_E_g_Mig->Draw("p");

  cMig4->Print("DVCSQA_temp022.pdf");
  cMig4->Close();

  TCanvas* cMig5 = new TCanvas("cMig5","",1200,800);
  cMig5->Divide(2,2);
  cMig5->cd(1);
  h_HFSSigma_Mig->GetXaxis()->SetTitleSize(0.05);
  h_HFSSigma_Mig->GetXaxis()->SetTitleOffset(0.9);
  h_HFSSigma_Mig->GetXaxis()->SetLabelSize(0.05);
  h_HFSSigma_Mig->GetYaxis()->SetTitleSize(0.05);
  h_HFSSigma_Mig->GetYaxis()->SetTitleOffset(0.88);
  h_HFSSigma_Mig->GetYaxis()->SetLabelSize(0.04);
  h_HFSSigma_Mig->Draw("p");
  tHead1->Draw(); // DRAW THIS LAST ON THE CANVAS
  cMig5->cd(2);
  h_HFSpT2_Mig->GetXaxis()->SetTitleSize(0.04);
  h_HFSpT2_Mig->GetXaxis()->SetTitleOffset(1.0);
  h_HFSpT2_Mig->GetXaxis()->SetLabelSize(0.04);
  h_HFSpT2_Mig->GetYaxis()->SetTitleSize(0.05);
  h_HFSpT2_Mig->GetYaxis()->SetTitleOffset(0.88);
  h_HFSpT2_Mig->GetYaxis()->SetLabelSize(0.04);
  h_HFSpT2_Mig->Draw("p");
  tHead2->Draw(); // DRAW THIS LAST ON THE CANVAS
  cMig5->cd(3);
  h_FullSigma_Mig->GetXaxis()->SetTitleSize(0.05);
  h_FullSigma_Mig->GetXaxis()->SetTitleOffset(0.9);
  h_FullSigma_Mig->GetXaxis()->SetLabelSize(0.05);
  h_FullSigma_Mig->GetYaxis()->SetTitleSize(0.05);
  h_FullSigma_Mig->GetYaxis()->SetTitleOffset(0.88);
  h_FullSigma_Mig->GetYaxis()->SetLabelSize(0.04);
  h_FullSigma_Mig->Draw("p");

  cMig5->Print("DVCSQA_temp023.pdf");
  cMig5->Close();

  TCanvas* ctRes = new TCanvas("ctResB0","",1200,800);
  // Extract t-resolutions from 2D distributions
  TH1D* h_b0_extracted_t_resolution;
  TH1D* h_rp_extracted_t_resolution;
  h_b0_extracted_t_resolution = extractResolution("b0_extracted_t_resolution",h_B0Inc_Res);
  h_rp_extracted_t_resolution = extractResolution("rp_extracted_t_resolution",h_RPInc_Res);
  
  // Set drawing options
  h_b0_extracted_t_resolution->GetXaxis()->SetTitle("|t|_{MC} [GeV^{2}]");
  h_b0_extracted_t_resolution->GetXaxis()->SetTitleSize(0.05);
  h_b0_extracted_t_resolution->GetXaxis()->SetTitleOffset(0.80);
  h_b0_extracted_t_resolution->GetYaxis()->SetTitle("#delta t [GeV^{2}]");
  h_b0_extracted_t_resolution->GetYaxis()->SetTitleSize(0.05);
  h_b0_extracted_t_resolution->GetYaxis()->SetTitleOffset(0.85);
  h_b0_extracted_t_resolution->SetTitle("");
  h_b0_extracted_t_resolution->SetLineColor(kBlue);
  h_b0_extracted_t_resolution->SetMarkerColor(kBlue);
  h_b0_extracted_t_resolution->SetMarkerStyle(21);
  h_b0_extracted_t_resolution->SetMarkerSize(1.3);
  h_rp_extracted_t_resolution->SetLineColor(kCyan+2);
  h_rp_extracted_t_resolution->SetMarkerColor(kCyan+2);
  h_rp_extracted_t_resolution->SetMarkerStyle(33);
  h_rp_extracted_t_resolution->SetMarkerSize(2.);

  // Draw
  h_b0_extracted_t_resolution->GetYaxis()->SetRangeUser(0,0.15);
  h_b0_extracted_t_resolution->Draw("ep");
  h_rp_extracted_t_resolution->Draw("epsame");
  // Legend
  TLegend* lInctRes = new TLegend(0.17,0.63,0.50,0.87);
  lInctRes->SetLineColorAlpha(kWhite,0);
  lInctRes->SetLineWidth(0);
  lInctRes->SetFillStyle(0);
  lInctRes->SetFillColorAlpha(kWhite,0);
  lInctRes->SetHeader("Single p'; no other cuts");
  lInctRes->AddEntry(h_b0_extracted_t_resolution,"B0","pl");
  lInctRes->AddEntry(h_rp_extracted_t_resolution,"RP","pl");
  lInctRes->Draw();
  // New header text (needed for single frame canvas)
  TLatex* tHeadFull1 = new TLatex(0.10, 0.91, sHead1);
  tHeadFull1->SetNDC();
  tHeadFull1->SetTextSize(35);
  tHeadFull1->SetTextFont(43);
  tHeadFull1->SetTextColor(kBlack);
  tHeadFull1->Draw();
  TLatex* tHeadFull2 = new TLatex(0.71, 0.91, sHead2);
  tHeadFull2->SetNDC();
  tHeadFull2->SetTextSize(35);
  tHeadFull2->SetTextFont(43);
  tHeadFull2->SetTextColor(kBlack);
  tHeadFull2->Draw();

  ctRes->Print("DVCSQA_temp024.pdf");
  ctRes->Close();

  // Combine PDFs into one and clean up
  std::cout<<"...Cleaning up files..."<<std::endl;
  TString filePlots = "$EIC_WORK_DIR/DVCS_Analysis/Plots/QA/DVCSPlots_" + campaign + "_" + energy + "_QA.pdf";
  std::cout<<"Moving plots to "<<filePlots<<std::endl;
  TString pdfUniteCmd = "pdfunite DVCSQA*.pdf "+filePlots;
  gSystem->Exec(pdfUniteCmd);
  gSystem->Exec("rm DVCSQA_temp*.pdf");
  
  return;
}
