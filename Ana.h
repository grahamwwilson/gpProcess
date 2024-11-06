//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct  8 11:12:40 2024 by ROOT version 6.32.04
// from TTree t/t
// found on file: Z-87.root
//////////////////////////////////////////////////////////

#ifndef Ana_h
#define Ana_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class Ana : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Long64_t> label1 = {fReader, "label1"};
   TTreeReaderValue<Double_t> E1 = {fReader, "E1"};
   TTreeReaderValue<Double_t> bx1 = {fReader, "bx1"};
   TTreeReaderValue<Double_t> by1 = {fReader, "by1"};
   TTreeReaderValue<Double_t> bz1 = {fReader, "bz1"};
   TTreeReaderValue<Double_t> x1 = {fReader, "x1"};
   TTreeReaderValue<Double_t> y1 = {fReader, "y1"};
   TTreeReaderValue<Double_t> z1 = {fReader, "z1"};
   TTreeReaderValue<Long64_t> f1 = {fReader, "f1"};
   TTreeReaderValue<Long64_t> label2 = {fReader, "label2"};
   TTreeReaderValue<Double_t> E2 = {fReader, "E2"};
   TTreeReaderValue<Double_t> bx2 = {fReader, "bx2"};
   TTreeReaderValue<Double_t> by2 = {fReader, "by2"};
   TTreeReaderValue<Double_t> bz2 = {fReader, "bz2"};
   TTreeReaderValue<Double_t> x2 = {fReader, "x2"};
   TTreeReaderValue<Double_t> y2 = {fReader, "y2"};
   TTreeReaderValue<Double_t> z2 = {fReader, "z2"};
   TTreeReaderValue<Long64_t> f2 = {fReader, "f2"};
   TTreeReaderValue<Long64_t> pabel1 = {fReader, "pabel1"};
   TTreeReaderValue<Double_t> E1p = {fReader, "E1p"};
   TTreeReaderValue<Double_t> bx1p = {fReader, "bx1p"};
   TTreeReaderValue<Double_t> by1p = {fReader, "by1p"};
   TTreeReaderValue<Double_t> bz1p = {fReader, "bz1p"};
   TTreeReaderValue<Double_t> x1p = {fReader, "x1p"};
   TTreeReaderValue<Double_t> y1p = {fReader, "y1p"};
   TTreeReaderValue<Double_t> z1p = {fReader, "z1p"};
   TTreeReaderValue<Long64_t> f1p = {fReader, "f1p"};
   TTreeReaderValue<Long64_t> pabel2 = {fReader, "pabel2"};
   TTreeReaderValue<Double_t> E2p = {fReader, "E2p"};
   TTreeReaderValue<Double_t> bx2p = {fReader, "bx2p"};
   TTreeReaderValue<Double_t> by2p = {fReader, "by2p"};
   TTreeReaderValue<Double_t> bz2p = {fReader, "bz2p"};
   TTreeReaderValue<Double_t> x2p = {fReader, "x2p"};
   TTreeReaderValue<Double_t> y2p = {fReader, "y2p"};
   TTreeReaderValue<Double_t> z2p = {fReader, "z2p"};
   TTreeReaderValue<Long64_t> f2p = {fReader, "f2p"};


   Ana(TTree * /*tree*/ =0) { }
   ~Ana() override { }
   Int_t   Version() const override { return 2; }
   void    Begin(TTree *tree) override;
   void    SlaveBegin(TTree *tree) override;
   void    Init(TTree *tree) override;
   bool    Notify() override;
   bool    Process(Long64_t entry) override;
   Int_t   GetEntry(Long64_t entry, Int_t getall = 0) override { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   void    SetOption(const char *option) override { fOption = option; }
   void    SetObject(TObject *obj) override { fObject = obj; }
   void    SetInputList(TList *input) override { fInput = input; }
   TList  *GetOutputList() const override { return fOutput; }
   void    SlaveTerminate() override;
   void    Terminate() override;

   ClassDefOverride(Ana,0);

};

#endif

#ifdef Ana_cxx
void Ana::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

bool Ana::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}


#endif // #ifdef Ana_cxx
