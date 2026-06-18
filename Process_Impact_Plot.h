// Stephen JD Kay - University of York
// Stephen.Kay@york.ac.uk
// 18/06/26
// Script to process EpIC generator output in combination with simulated output (efficiencies per x/Q2/t bin)
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
