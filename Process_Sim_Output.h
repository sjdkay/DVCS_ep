// Stephen JD Kay - University of York
// Stephen.Kay@york.ac.uk
// 09/06/26
// Script to process EpIC generator output and plot real -t/x/Q2 distributions
// Header file

Bool_t CheckFile(TString File){
  Bool_t FileCheck = kTRUE;
  if (gSystem->AccessPathName(File) == kTRUE){
    FileCheck = kFALSE;
    cout << "Input file - " << File << " not found." << endl;
    cout << "Check pathing and edit input argument as needed!" << endl;
  }
  else{ // File exists, try to open it!
    TFile *tmpfile = TFile::Open(File);
    if (tmpfile && !tmpfile->IsZombie()) {
      tmpfile->Close();  // Close the file
      delete tmpfile;    // Delete the file pointer
    }
    else {
      cout << "Failed to open file: " << File <<endl;
      FileCheck = kFALSE;
      if (tmpfile) {
	tmpfile->Close();  // Close the file if it was partially opened
	delete tmpfile;    // Delete the file pointer
      }
    }
  }
  return FileCheck;
}

// Function to get a TGraphErrors of the resolution given an input 2D histogram of DeltaQuant as a fn of Quant
TGraphErrors* extractResolution(TH2D* twoDHisto){
  int num_bins  = twoDHisto->GetNbinsX();
  double xBinWidth = twoDHisto->GetXaxis()->GetBinWidth(1);
  double xMin = twoDHisto->GetXaxis()->GetBinCenter(1) - xBinWidth*0.5;
  double xMax = twoDHisto->GetXaxis()->GetBinCenter(num_bins) + xBinWidth*0.5;
  std::vector<double> x_vals;
  std::vector<double> xerr_vals;
  std::vector<double> y_vals;
  std::vector<double> yerr_vals;
  TH1D* tmp;
  double rmsReso = 0.0;
  double rmsErr = 0.0;
  for(int bin = 1; bin < num_bins+1; bin++){
    rmsReso = 0.0;
    tmp = (TH1D*)twoDHisto->ProjectionY("tmp_proj", bin, bin+1);
    rmsReso = tmp->GetRMS();
    rmsErr  = tmp->GetRMSError();
    x_vals.push_back(twoDHisto->GetXaxis()->GetBinCenter(bin));
    xerr_vals.push_back(xBinWidth*0.5);
    y_vals.push_back(rmsReso);
    yerr_vals.push_back(rmsErr);
    //delete func;
  }
  TGraphErrors* finalResoGraph = new TGraphErrors(num_bins, x_vals.data(), y_vals.data(), xerr_vals.data(), yerr_vals.data());
  return finalResoGraph;
}
