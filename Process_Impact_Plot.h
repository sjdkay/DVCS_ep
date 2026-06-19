// Stephen JD Kay - University of York
// Stephen.Kay@york.ac.uk
// 18/06/26
// Script to process EpIC generator output in combination with simulated output (efficiencies per x/Q2/t bin)
// Header file

PxPyPzEVector Vec_e_beam;
PxPyPzEVector Vec_p_beam;
PxPyPzEVector Vec_eSc;
PxPyPzEVector Vec_pSc;
PxPyPzEVector Vec_gamma;
PxPyPzEVector Vec_Q2;
PxPyPzEVector Vec_t;
PxPyPzEVector Vec_CoM_Boost;
PxPyPzEVector Vec_HoF_Boost;
XYZVector Boost_CoM;
XYZVector Boost_HoF;
RotationX RotX;
RotationY RotY;

Double_t Q2, xB, t, y;
Bool_t eBeam, pBeam, eSc, pSc, Gamma, GoodEvent;

Bool_t CheckFile_Gen(TString File){
  Bool_t FileCheck = kTRUE;
  if (gSystem->AccessPathName(File) == kTRUE){
    FileCheck = kFALSE;
    cout << "Input file - " << File << " not found." << endl;
    cout << "Check pathing and edit input argument as needed!" << endl;
  }
  else{ // File exists, try to open it!
    TFile *tmpfile = TFile::Open(File);
    if (tmpfile && !tmpfile->IsZombie()) {
      TTree *tmptree = (TTree*)tmpfile->Get("hepmc3_tree");
      // Ensure the file also contains a valid tree
      if (!tmptree) {
	cout << "Tree 'hepmc3_tree' not found in generator file: " << File << endl;
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

Bool_t CheckFile_Sim(TString File){
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
