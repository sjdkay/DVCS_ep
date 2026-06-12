// Stephen JD Kay - University of York
// Stephen.Kay@york.ac.uk
// 09/06/26
// Script to process EpIC generator output and plot real -t/x/Q2 distributions
// Header file

PxPyPzEVector Vec_e_beam;
PxPyPzEVector Vec_p_beam;
PxPyPzEVector Vec_eSc;
PxPyPzEVector Vec_pSc;
PxPyPzEVector Vec_gamma;
PxPyPzEVector Vec_Q2;
PxPyPzEVector Vec_t;

Double_t Q2, xB, t, y;
Bool_t eBeam, pBeam, eSc, pSc, Gamma, GoodEvent;

RotationY rY(0.025);

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
      TTree *tmptree = (TTree*)tmpfile->Get("events");
      // Ensure the file also contains a valid tree
      if (!tmptree) {
	cout << "Tree 'events' not found in file: " << File << endl;
	FileCheck = kFALSE;

      }
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
