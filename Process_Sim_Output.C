// Stephen JD Kay - University of York
// Stephen.Kay@york.ac.uk
// 09/06/26
// Script to process EpIC generator output and plot real -t/x/Q2 distributions

using namespace ROOT::Math;

#include "TString.h"
#include "Process_Sim_Output.h"
#include "ePICStyle.C"
#include <vector>

void Process_Sim_Output(TString InFile = ""){

  if(CheckFile(InFile) == kFALSE){ // Check files exist, can be opened and contain tree with fn
    exit(1);
  }
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  
  TFile *File =  new TFile(InFile);
  
  //gDirectory->cd("Q2xB_Binned_Dists");
  // Just need to open other file in this case and get hists
  //  TFile *ofile = TFile::Open(Form("%s_Q2xB_Hists.root", (InFile.Remove(InFile.Length() - 5)).Data()),"RECREATE"); // 5 cuts the .root from the original input file name

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
  TCanvas* c_Q2xB_Results[9];
  TLatex *Q2_Range_Text[9];
  TLegend* Leg_tComp = new TLegend (0.1, 0.2, 0.9, 0.6);
  TString OutPdf = Form("%s_TrueQ2xB_Binned_tDists_Sim.pdf", (InFile.Remove(InFile.Length() - 5)).Data());
  for(int binq2{0}; binq2<nQ2bins; binq2++){
    c_Q2xB_Results[binq2] = new TCanvas(Form("c_Q2xB_Results_%i", binq2+1), Form("t Dists across xB bins, Q2 %i", binq2+1), 100, 0, 2560, 1920);
    c_Q2xB_Results[binq2]->Divide(4,3); 
    for(int binxB{0}; binxB<nxBbins; binxB++){
      c_Q2xB_Results[binq2]->cd(binxB+1);
      tmpHist1D = (TH1D*)(((TH1D*)File->Get(Form("Q2xB_Binned_Dists/h1_tMC_Q2xB_True[%i][%i]",binq2,binxB))));
      tmpHist1D->SetTitle(Form("%.2e<x_{B}<%.2e", xBedges[binxB],xBedges[binxB+1]));
      tmpHist1D->SetLineColor(kP6Gray);
      if(binq2 == 0 && binxB == 0){ //First bin only
	Leg_tComp->AddEntry(tmpHist1D, "t_{MC}");
      }
      tmpHist1D->Draw("HISTERR");
      tmpHist1D = (TH1D*)(((TH1D*)File->Get(Form("Q2xB_Binned_Dists/h1_tMCAcc_Q2xB_True[%i][%i]",binq2,binxB))));
      tmpHist1D->SetTitle(Form("%.2e<x_{B}<%.2e", xBedges[binxB],xBedges[binxB+1]));
      tmpHist1D->SetLineColor(kP6Violet);
      if(binq2 == 0 && binxB == 0){ //First bin only
	Leg_tComp->AddEntry(tmpHist1D, "t_{MC_Accepted}");
      }
      tmpHist1D->Draw("SAMEHISTERR");
      tmpHist1D = (TH1D*)(((TH1D*)File->Get(Form("Q2xB_Binned_Dists/h1_tMethL_Q2xB_True[%i][%i]",binq2,binxB))));
      tmpHist1D->SetTitle(Form("%.2e<x_{B}<%.2e", xBedges[binxB],xBedges[binxB+1]));
      tmpHist1D->SetLineColor(kP6Red);
      if(binq2 == 0 && binxB == 0){ //First bin only
	Leg_tComp->AddEntry(tmpHist1D, "t_{MethodL}");
      }
      tmpHist1D->Draw("SAMEHISTERR");
      tmpHist1D = (TH1D*)(((TH1D*)File->Get(Form("Q2xB_Binned_Dists/h1_tB0_Q2xB_True[%i][%i]",binq2,binxB))));
      tmpHist1D->SetTitle(Form("%.2e<x_{B}<%.2e", xBedges[binxB],xBedges[binxB+1]));
      tmpHist1D->SetLineColor(kP6Blue);
      if(binq2 == 0 && binxB == 0){ //First bin only
	Leg_tComp->AddEntry(tmpHist1D, "t_{B0}");
      }
      tmpHist1D->Draw("SAMEHISTERR");
      tmpHist1D = (TH1D*)(((TH1D*)File->Get(Form("Q2xB_Binned_Dists/h1_tRP_Q2xB_True[%i][%i]",binq2,binxB))));
      tmpHist1D->SetTitle(Form("%.2e<x_{B}<%.2e", xBedges[binxB],xBedges[binxB+1]));
      tmpHist1D->SetLineColor(kP6Yellow);
      if(binq2 == 0 && binxB == 0){ //First bin only
	Leg_tComp->AddEntry(tmpHist1D, "t_{RP}");
      }
      tmpHist1D->Draw("SAMEHISTERR");
      gPad->SetLogy(1);
    }
    c_Q2xB_Results[binq2]->cd(12);
    Q2_Range_Text[binq2] = new TLatex(0.2, 0.8, Form("%.1f<Q^{2}<%.1f GeV^{2}", q2edges[binq2],q2edges[binq2+1]));
    Q2_Range_Text[binq2]->Draw();
    Leg_tComp->Draw("SAME");
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

  TCanvas* c_Q2xB_Results2[9];
  TString OutPdf2 = Form("%s_RecQ2xB_Binned_tDists_Rec.pdf", InFile.Data());
  for(int binq2{0}; binq2<nQ2bins; binq2++){
    c_Q2xB_Results2[binq2] = new TCanvas(Form("c_Q2xB_Results2_%i", binq2+1), Form("t Dists across xB bins, Q2 %i", binq2+1), 100, 0, 2560, 1920);
    c_Q2xB_Results2[binq2]->Divide(4,3); 
    for(int binxB{0}; binxB<nxBbins; binxB++){
      c_Q2xB_Results2[binq2]->cd(binxB+1);
      tmpHist1D = (TH1D*)(((TH1D*)File->Get(Form("Q2xB_Binned_Dists/h1_tMC_Q2xB_True[%i][%i]",binq2,binxB))));
      tmpHist1D->SetTitle(Form("%.2e<x_{B}<%.2e", xBedges[binxB],xBedges[binxB+1]));
      tmpHist1D->SetLineColor(kP6Gray);
      tmpHist1D->Draw("HISTERR");
      tmpHist1D = (TH1D*)(((TH1D*)File->Get(Form("Q2xB_Binned_Dists/h1_tMCAcc_Q2xB_Rec[%i][%i]",binq2,binxB))));
      tmpHist1D->SetTitle(Form("%.2e<x_{B}<%.2e", xBedges[binxB],xBedges[binxB+1]));
      tmpHist1D->SetLineColor(kP6Violet);
      tmpHist1D->Draw("SAMEHISTERR");
      tmpHist1D = (TH1D*)(((TH1D*)File->Get(Form("Q2xB_Binned_Dists/h1_tMethL_Q2xB_Rec[%i][%i]",binq2,binxB))));
      tmpHist1D->SetTitle(Form("%.2e<x_{B}<%.2e", xBedges[binxB],xBedges[binxB+1]));
      tmpHist1D->SetLineColor(kP6Red);
      tmpHist1D->Draw("SAMEHISTERR");
      tmpHist1D = (TH1D*)(((TH1D*)File->Get(Form("Q2xB_Binned_Dists/h1_tB0_Q2xB_Rec[%i][%i]",binq2,binxB))));
      tmpHist1D->SetTitle(Form("%.2e<x_{B}<%.2e", xBedges[binxB],xBedges[binxB+1]));
      tmpHist1D->SetLineColor(kP6Blue);
      tmpHist1D->Draw("SAMEHISTERR");
      tmpHist1D = (TH1D*)(((TH1D*)File->Get(Form("Q2xB_Binned_Dists/h1_tRP_Q2xB_Rec[%i][%i]",binq2,binxB))));
      tmpHist1D->SetTitle(Form("%.2e<x_{B}<%.2e", xBedges[binxB],xBedges[binxB+1]));
      tmpHist1D->SetLineColor(kP6Yellow);
      tmpHist1D->Draw("SAMEHISTERR");
      gPad->SetLogy(1);
    }
    c_Q2xB_Results2[binq2]->cd(12);
    Q2_Range_Text[binq2] = new TLatex(0.2, 0.8, Form("%.1f<Q^{2}<%.1f GeV^{2}", q2edges[binq2],q2edges[binq2+1]));
    Q2_Range_Text[binq2]->Draw();
    Leg_tComp->Draw("SAME");
    if(binq2 == 0){
      c_Q2xB_Results2[binq2]->Print(OutPdf2 + "(");
    }
    else if(binq2 == 7){
      c_Q2xB_Results2[binq2]->Print(OutPdf2 + ")");
    }
    else{
      c_Q2xB_Results2[binq2]->Print(OutPdf2);
    }
  }

  TCanvas* c_tRes_Results =  new TCanvas("c_tResResults", "t Resolutions", 100, 0, 2560, 1920);
  TString OutPdf3 = Form("%s_tRes_Results.pdf", InFile.Data());
  TGraphErrors *tResGraphs[3];
  c_tRes_Results->Divide(3,2);
  c_tRes_Results->cd(1);
  tmpHist2D = (TH2D*)(((TH2D*)File->Get("t_Resolution_Plots/h2_tResB0Pct")));
  tmpHist2D->SetTitle("#Deltat/t_{MC}(t_{MC}) - p' in B0");
  tmpHist2D->Draw("COLZ");
  gPad->SetLogz(1);
  c_tRes_Results->cd(4);
  tResGraphs[0]=extractResolution(tmpHist2D);
  tResGraphs[0]->SetMarkerColor(kP6Red);
  tResGraphs[0]->SetMarkerStyle(21);
  tResGraphs[0]->GetXaxis()->SetRangeUser(0, 1.6);
  tResGraphs[0]->GetYaxis()->SetRangeUser(0, 20);
  tResGraphs[0]->SetTitle("");
  tResGraphs[0]->GetXaxis()->SetTitle("-t_{MC} (GeV^{2})");
  tResGraphs[0]->GetYaxis()->SetTitle("RMS(#Delta t) (%)");
  tResGraphs[0]->Draw("AP");
  c_tRes_Results->cd(2);
  tmpHist2D = (TH2D*)(((TH2D*)File->Get("t_Resolution_Plots/h2_tResRPPct")));
  tmpHist2D->SetTitle("#Deltat/t_{MC}(t_{MC}) - p' in RP");
  tmpHist2D->Draw("COLZ");
  gPad->SetLogz(1);
  c_tRes_Results->cd(5);
  tResGraphs[1]=extractResolution(tmpHist2D);
  tResGraphs[1]->SetMarkerColor(kP6Red);
  tResGraphs[1]->SetMarkerStyle(21);
  tResGraphs[1]->GetXaxis()->SetRangeUser(0, 0.8);
  tResGraphs[1]->GetYaxis()->SetRangeUser(0, 35);
  tResGraphs[1]->SetTitle("");
  tResGraphs[1]->GetXaxis()->SetTitle("-t_{MC} (GeV^{2})");
  tResGraphs[1]->GetYaxis()->SetTitle("RMS(#Delta t) (%)");
  tResGraphs[1]->Draw("AP");
  c_tRes_Results->cd(3);
  tmpHist2D = (TH2D*)(((TH2D*)File->Get("t_Resolution_Plots/h2_tResLCPct")));
  tmpHist2D->SetTitle("#Deltat/t_{MC}(t_{MC}) - Method L");
  tmpHist2D->Draw("COLZ");
  gPad->SetLogz(1);
  c_tRes_Results->cd(6);
  tResGraphs[2]=extractResolution(tmpHist2D);
  tResGraphs[2]->SetMarkerColor(kP6Red);
  tResGraphs[2]->SetMarkerStyle(21);
  tResGraphs[2]->GetXaxis()->SetRangeUser(0, 1.6);
  tResGraphs[2]->GetYaxis()->SetRangeUser(0, 65);
  tResGraphs[2]->SetTitle("");
  tResGraphs[2]->GetXaxis()->SetTitle("-t_{MC} (GeV^{2})");
  tResGraphs[2]->GetYaxis()->SetTitle("RMS(#Delta t) (%)");
  tResGraphs[2]->Draw("AP");
  c_tRes_Results->Print(OutPdf3);

}
