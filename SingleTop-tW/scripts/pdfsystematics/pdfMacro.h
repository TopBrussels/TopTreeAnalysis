//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul 13 03:28:14 2011 by ROOT version 5.27/06b
// from TChain myTree/myTree
//////////////////////////////////////////////////////////

#ifndef pdfMacro_h
#define pdfMacro_h

#include <TROOT.h>
#include <TChain.h>
#include <TString.h>
#include <TFile.h>
#include <iostream>
#include <string>

class pdfMacro {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   std::string     pdfSet_;
   void SetPdfSetName(TString value){pdfSet_=value;}

   // Declaration of leaf types
   Double_t        xlWeight;
   Double_t        luminosity;
   Double_t        metPx;
   Double_t        metPy;
   Int_t           id1;
   Int_t           id2;
   Float_t         x1;
   Float_t         x2;
   Float_t         q;
   vector<double>  *pxLepton;
   vector<double>  *pyLepton;
   vector<double>  *pzLepton;
   vector<double>  *eLepton;
   vector<double>  *qLepton;
   vector<double>  *pxJet;
   vector<double>  *pyJet;
   vector<double>  *pzJet;
   vector<double>  *eJet;
   vector<double>  *qJet;
   vector<double>  *btSSVHEJet;

   // List of branches
   TBranch        *b_xlWeight;   //!
   TBranch        *b_luminosity;   //!
   TBranch        *b_metPx;   //!
   TBranch        *b_metPy;   //!
   TBranch        *b_id1;   //!
   TBranch        *b_id2;   //!
   TBranch        *b_x1;   //!
   TBranch        *b_x2;   //!
   TBranch        *b_q;   //!
   TBranch        *b_pxLepton;   //!
   TBranch        *b_pyLepton;   //!
   TBranch        *b_pzLepton;   //!
   TBranch        *b_eLepton;   //!
   TBranch        *b_qLepton;   //!
   TBranch        *b_pxJet;   //!
   TBranch        *b_pyJet;   //!
   TBranch        *b_pzJet;   //!
   TBranch        *b_eJet;   //!
   TBranch        *b_qJet;   //!
   TBranch        *b_btSSVHEJet;   //!

   pdfMacro(TTree *tree=0);
   virtual ~pdfMacro();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual void ReadFromTextFile(TString textfile,float scaler);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef pdfMacro_cxx
pdfMacro::pdfMacro(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f) {
         f = new TFile("Memory Directory");
         f->cd("Rint:/");
      }
      tree = (TTree*)gDirectory->Get("myTree");

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("myTree","myTree");
      chain->Add("/afs/cern.ch/user/f/fblekman/pdfs/tW/micro_ee_tt.root/myTree");
      chain->Add("/afs/cern.ch/user/f/fblekman/pdfs/tW/micro_ee_tw.root/myTree");
      chain->Add("/afs/cern.ch/user/f/fblekman/pdfs/tW/micro_emu_tt.root/myTree");
      chain->Add("/afs/cern.ch/user/f/fblekman/pdfs/tW/micro_emu_tw.root/myTree");
      chain->Add("/afs/cern.ch/user/f/fblekman/pdfs/tW/micro_mumu_tt.root/myTree");
      chain->Add("/afs/cern.ch/user/f/fblekman/pdfs/tW/micro_mumu_tw.root/myTree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

pdfMacro::~pdfMacro()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pdfMacro::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pdfMacro::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void pdfMacro::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pxLepton = 0;
   pyLepton = 0;
   pzLepton = 0;
   eLepton = 0;
   qLepton = 0;
   pxJet = 0;
   pyJet = 0;
   pzJet = 0;
   eJet = 0;
   qJet = 0;
   btSSVHEJet = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("xlWeight", &xlWeight, &b_xlWeight);
   fChain->SetBranchAddress("luminosity", &luminosity, &b_luminosity);
   fChain->SetBranchAddress("metPx", &metPx, &b_metPx);
   fChain->SetBranchAddress("metPy", &metPy, &b_metPy);
   fChain->SetBranchAddress("id1", &id1, &b_id1);
   fChain->SetBranchAddress("id2", &id2, &b_id2);
   fChain->SetBranchAddress("x1", &x1, &b_x1);
   fChain->SetBranchAddress("x2", &x2, &b_x2);
   fChain->SetBranchAddress("q", &q, &b_q);
   fChain->SetBranchAddress("pxLepton", &pxLepton, &b_pxLepton);
   fChain->SetBranchAddress("pyLepton", &pyLepton, &b_pyLepton);
   fChain->SetBranchAddress("pzLepton", &pzLepton, &b_pzLepton);
   fChain->SetBranchAddress("eLepton", &eLepton, &b_eLepton);
   fChain->SetBranchAddress("qLepton", &qLepton, &b_qLepton);
   fChain->SetBranchAddress("pxJet", &pxJet, &b_pxJet);
   fChain->SetBranchAddress("pyJet", &pyJet, &b_pyJet);
   fChain->SetBranchAddress("pzJet", &pzJet, &b_pzJet);
   fChain->SetBranchAddress("eJet", &eJet, &b_eJet);
   fChain->SetBranchAddress("qJet", &qJet, &b_qJet);
   fChain->SetBranchAddress("btSSVHEJet", &btSSVHEJet, &b_btSSVHEJet);
   Notify();
}

Bool_t pdfMacro::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pdfMacro::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pdfMacro::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef pdfMacro_cxx
