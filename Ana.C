#define Ana_cxx
// The class definition in Ana.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("Ana.C")
// root> T->Process("Ana.C","some options")
// root> T->Process("Ana.C+")
//

#include <cmath>
#include "Ana.h"
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <TStyle.h>

const double EBNOM = 45.6;
// Declare histograms globally here
TH1D *hE1 = new TH1D("hE1","; E1/Enominal; Events per bin ",1200,0.95,1.01); 
TH1D *hE2 = new TH1D("hE2","; E2/Enominal; Events per bin ",1200,0.95,1.01); 
TH1D *hECM = new TH1D("hECM","; ECM/ECMnominal; Events per bin ",1200,0.95,1.01);
TH1D *hECMp = new TH1D("hECMp","; ECMp/ECMnominal; Events per bin ",1200,0.95,1.01); 
TH1D *hdECM = new TH1D("hdECM","; (ECMp - ECM)/ECMnominal; Events per bin ",1200,-0.05,0.01); 

TH1D *heleDpx = new TH1D("heleDpx", "; Electron (pxp - px)/Ebnominal; Events per bin ",20000,-0.005,0.005);
TH1D *heleDpy = new TH1D("heleDpy", "; Electron (pyp - py)/Ebnominal; Events per bin ",20000,-0.005,0.005);
TH1D *heleDpz = new TH1D("heleDpz", "; Electron (pzp - pz)/Ebnominal; Events per bin ",20000,-0.005,0.005);
TH1D *heleDp  = new TH1D("heleDp",  "; Electron (pp - p)/Ebnominal; Events per bin ",20000,-0.005,0.005);
TH1D *hposDpx = new TH1D("hposDpx", "; Positron (pxp - px)/Ebnominal; Events per bin ",20000,-0.005,0.005);
TH1D *hposDpy = new TH1D("hposDpy", "; Positron (pyp - py)/Ebnominal; Events per bin ",20000,-0.005,0.005);
TH1D *hposDpz = new TH1D("hposDpz", "; Positron (pzp - pz)/Ebnominal; Events per bin ",20000,-0.005,0.005);
TH1D *hposDp  = new TH1D("hposDp",  "; Positron (pp - p)/Ebnominal; Events per bin ",20000,-0.005,0.005);

TH1D *heleDeflectionTheta =  new TH1D("heleDeflectionTheta",   "; Electron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *heleDeflectionTheta0 = new TH1D("heleDeflectionTheta0",  "; Electron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *heleDeflectionTheta1 = new TH1D("heleDeflectionTheta1",  "; Electron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *heleDeflectionTheta2 = new TH1D("heleDeflectionTheta2",  "; Electron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *heleDeflectionTheta3 = new TH1D("heleDeflectionTheta3",  "; Electron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *heleDeflectionTheta4 = new TH1D("heleDeflectionTheta4",  "; Electron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *heleDeflectionTheta5 = new TH1D("heleDeflectionTheta5",  "; Electron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *heleDeflectionTheta6 = new TH1D("heleDeflectionTheta6",  "; Electron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *heleDeflectionTheta7 = new TH1D("heleDeflectionTheta7",  "; Electron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *heleDeflectionTheta8 = new TH1D("heleDeflectionTheta8",  "; Electron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *heleDeflectionTheta9 = new TH1D("heleDeflectionTheta9",  "; Electron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);

TH1D *heleDeflectionTTheta0 = new TH1D("heleDeflectionTTheta0",  "; Electron EM theta deflection [urad]; Events per bin ",80,-40.0,40.0);
TH1D *heleDeflectionTTheta1 = new TH1D("heleDeflectionTTheta1",  "; Electron EM theta deflection [urad]; Events per bin ",80,-40.0,40.0);

TH1D *hposDeflectionTheta =  new TH1D("hposDeflectionTheta",   "; Positron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *hposDeflectionTheta0 = new TH1D("hposDeflectionTheta0",  "; Positron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *hposDeflectionTheta1 = new TH1D("hposDeflectionTheta1",  "; Positron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *hposDeflectionTheta2 = new TH1D("hposDeflectionTheta2",  "; Positron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *hposDeflectionTheta3 = new TH1D("hposDeflectionTheta3",  "; Positron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *hposDeflectionTheta4 = new TH1D("hposDeflectionTheta4",  "; Positron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *hposDeflectionTheta5 = new TH1D("hposDeflectionTheta5",  "; Positron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *hposDeflectionTheta6 = new TH1D("hposDeflectionTheta6",  "; Positron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *hposDeflectionTheta7 = new TH1D("hposDeflectionTheta7",  "; Positron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *hposDeflectionTheta8 = new TH1D("hposDeflectionTheta8",  "; Positron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);
TH1D *hposDeflectionTheta9 = new TH1D("hposDeflectionTheta9",  "; Positron EM theta deflection [urad]; Events per bin ",125,-30.0,370.0);

TH1D *hposDeflectionTTheta0 = new TH1D("hposDeflectionTTheta0",  "; Positron EM theta deflection [urad]; Events per bin ",80,-40.0,40.0);
TH1D *hposDeflectionTTheta1 = new TH1D("hposDeflectionTTheta1",  "; Positron EM theta deflection [urad]; Events per bin ",80,-40.0,40.0);

TH1D *heleDeflectionPhi = new TH1D("heleDeflectionPhi",  "; Electron EM phi deflection [mrad]; Events per bin ",100,-2.0,2.0);
TH1D *hposDeflectionPhi = new TH1D("hposDeflectionPhi",  "; Positron EM phi deflection [mrad]; Events per bin ",100,-2.0,2.0);

TProfile *heleDeflThetavsPhi = new TProfile("heleDeflThetavsPhi", "; Phi [rad]; Average EM Theta Deflection [urad]",40,-M_PI,M_PI,-20.0,300.0);
TProfile *hposDeflThetavsPhi = new TProfile("hposDeflThetavsPhi", "; Phi [rad]; Average EM Theta Deflection [urad]",40,-M_PI,M_PI,-20.0,300.0);
TProfile *heleDeflPhivsPhi = new TProfile("heleDeflPhivsPhi", "; Phi [rad]; Average EM Phi Deflection [mrad]",40,-M_PI,M_PI,-2.0,2.0);
TProfile *hposDeflPhivsPhi = new TProfile("hposDeflPhivsPhi", "; Phi [rad]; Average EM Phi Deflection [mrad]",40,-M_PI,M_PI,-2.0,2.0);

TProfile *heleDeflThetavsZ = new TProfile("heleDeflThetavsZ", "; Zvtx [um]; Average EM Theta Deflection [urad]",40,-1000.0,1000.0,-20.0,300.0);
TProfile *hposDeflThetavsZ = new TProfile("hposDeflThetavsZ", "; Zvtx [um]; Average EM Theta Deflection [urad]",40,-1000.0,1000.0,-20.0,300.0);

int NEVS = 1007105.0;
TProfile *heleDeflThetavsT = new TProfile("heleDeflThetavsT", "; t; Average EM Theta Deflection [urad]",40,-0.5,double(NEVS)+0.5,-20.0,300.0);
TProfile *hposDeflThetavsT = new TProfile("hposDeflThetavsT", "; t; Average EM Theta Deflection [urad]",40,-0.5,double(NEVS)+0.5,-20.0,300.0);


struct fvec{
//   double FVector[4];
     double px;
     double py;
     double pz;
     double E;
};

struct dilepton{
   fvec Lepton0;
   fvec Lepton1;
};

double momentum(struct fvec fv){
  
    double px = fv.px;
    double py = fv.py;
    double pz = fv.pz;
    double p = sqrt(px*px + py*py + pz*pz);
    return p;
    
}

double theta(struct fvec fv){
  
    double px = fv.px;
    double py = fv.py;
    double pz = fv.pz;
    double pt = sqrt(px*px + py*py);
    double p = sqrt(px*px + py*py + pz*pz);

//    double costh = pz/p;
//    double theta = acos(costh);
    
    double sinth = pt/p;
    double theta = asin(sinth);
    if(pz < 0.0)theta = M_PI - theta;
    
    return theta;
    
}

double phi(struct fvec fv){
  
    double px = fv.px;
    double py = fv.py;

    double phival = atan2(py,px);
    
    return phival;
    
}

double pt(struct fvec fv){
  
    double px = fv.px;
    double py = fv.py;
    double val = sqrt(px*px + py*py);
    return val;
    
}

double twobodymass(struct dilepton ll){
 
    fvec fv1 = ll.Lepton0;
    fvec fv2 = ll.Lepton1;
     
    double E12  = fv1.E  + fv2.E;
    double p12x = fv1.px + fv2.px;     
    double p12y = fv1.py + fv2.py; 
    double p12z = fv1.pz + fv2.pz;     
    
    double val = sqrt(E12*E12 - p12x*p12x - p12y*p12y - p12z*p12z);
    return val;
    
}    


void Ana::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();  
   TH1::SetDefaultSumw2(kTRUE);

}

void Ana::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

bool Ana::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);
   
// Based on pairs0 information close to collision point
   fvec fele;   fvec fpos;
   fele.px =  (*E1)*(*bx1);
   fele.py =  (*E1)*(*by1);
   fele.pz =  (*E1)*(*bz1);
   fele.E =    *E1;        
   
   fpos.px = -(*E2)*(*bx2);
   fpos.py = -(*E2)*(*by2);
   fpos.pz = -(*E2)*(*bz2);
   fpos.E =  -(*E2);    
   
   fvec felep; fvec fposp;
   felep.px =  (*E1p)*(*bx1p);
   felep.py =  (*E1p)*(*by1p);
   felep.pz =  (*E1p)*(*bz1p);
   felep.E =    *E1p;        
   
   fposp.px = -(*E2p)*(*bx2p);
   fposp.py = -(*E2p)*(*by2p);
   fposp.pz = -(*E2p)*(*bz2p);
   fposp.E =  -(*E2p);   
    
   dilepton bhabha0, bhabhap;
   bhabha0.Lepton0 = fele;
   bhabha0.Lepton1 = fpos;
   
// Apparent dilepton 4-vectors after propagation through field
   bhabhap.Lepton0 = felep;
   bhabhap.Lepton1 = fposp;   
   
   double mll = twobodymass(bhabha0);
   double mllp = twobodymass(bhabhap);   
   
   hE1->Fill((*E1)/EBNOM);
   hE2->Fill(-(*E2)/EBNOM); 
   hECM->Fill(mll/(2.0*EBNOM)); 
   hECMp->Fill(mllp/(2.0*EBNOM));  
   hdECM->Fill( (mllp - mll)/(2.0*EBNOM));  
   
   heleDpx->Fill( (felep.px - fele.px)/EBNOM ); 
   heleDpy->Fill( (felep.py - fele.py)/EBNOM );
   heleDpz->Fill( (felep.pz - fele.pz)/EBNOM ); 
   heleDp->Fill(  (momentum(felep) - momentum(fele))/EBNOM );      
   
   hposDpx->Fill( (fposp.px - fpos.px)/EBNOM ); 
   hposDpy->Fill( (fposp.py - fpos.py)/EBNOM );
   hposDpz->Fill( (fposp.pz - fpos.pz)/EBNOM );
   hposDp->Fill(  (momentum(fposp) - momentum(fpos))/EBNOM ); 
   
   double thele = theta(fele);
   double thelep = theta(felep);
   double thpos = theta(fpos);
   double thposp = theta(fposp);   
   
   heleDeflectionTheta->Fill(1.0e6*(thele - thelep));  
   hposDeflectionTheta->Fill(-1.0e6*(thpos - thposp)); 
   
   if ( (*label1)%10==0 ) heleDeflectionTheta0->Fill(1.0e6*(thele - thelep));
   if ( (*label1)%10==1 ) heleDeflectionTheta1->Fill(1.0e6*(thele - thelep));
   if ( (*label1)%10==2 ) heleDeflectionTheta2->Fill(1.0e6*(thele - thelep));
   if ( (*label1)%10==3 ) heleDeflectionTheta3->Fill(1.0e6*(thele - thelep));
   if ( (*label1)%10==4 ) heleDeflectionTheta4->Fill(1.0e6*(thele - thelep));  
   if ( (*label1)%10==5 ) heleDeflectionTheta5->Fill(1.0e6*(thele - thelep));
   if ( (*label1)%10==6 ) heleDeflectionTheta6->Fill(1.0e6*(thele - thelep));
   if ( (*label1)%10==7 ) heleDeflectionTheta7->Fill(1.0e6*(thele - thelep));
   if ( (*label1)%10==8 ) heleDeflectionTheta8->Fill(1.0e6*(thele - thelep));
   if ( (*label1)%10==9 ) heleDeflectionTheta9->Fill(1.0e6*(thele - thelep)); 
   
   if ( (*label1)%2==0 ) heleDeflectionTTheta0->Fill(1.0e6*(thele - thelep));
   if ( (*label1)%2==1 ) heleDeflectionTTheta1->Fill(1.0e6*(thele - thelep));     
   
   if ( (*label1)%10==0 ) hposDeflectionTheta0->Fill(-1.0e6*(thpos - thposp));
   if ( (*label1)%10==1 ) hposDeflectionTheta1->Fill(-1.0e6*(thpos - thposp));
   if ( (*label1)%10==2 ) hposDeflectionTheta2->Fill(-1.0e6*(thpos - thposp));
   if ( (*label1)%10==3 ) hposDeflectionTheta3->Fill(-1.0e6*(thpos - thposp));
   if ( (*label1)%10==4 ) hposDeflectionTheta4->Fill(-1.0e6*(thpos - thposp));  
   if ( (*label1)%10==5 ) hposDeflectionTheta5->Fill(-1.0e6*(thpos - thposp));
   if ( (*label1)%10==6 ) hposDeflectionTheta6->Fill(-1.0e6*(thpos - thposp));
   if ( (*label1)%10==7 ) hposDeflectionTheta7->Fill(-1.0e6*(thpos - thposp));
   if ( (*label1)%10==8 ) hposDeflectionTheta8->Fill(-1.0e6*(thpos - thposp));
   if ( (*label1)%10==9 ) hposDeflectionTheta9->Fill(-1.0e6*(thpos - thposp)); 
   
   if ( (*label1)%2==0 ) hposDeflectionTTheta0->Fill(-1.0e6*(thpos - thposp));
   if ( (*label1)%2==1 ) hposDeflectionTTheta1->Fill(-1.0e6*(thpos - thposp));                           
   
   if( 1.0e6*(thele - thelep) < -20.0 ){
      cout << "e- Underflow event (urad) " << (*label1) << " " << 1.0e6*thele << " " << 1.0e6*thelep << " " << 1.0e6*(thele - thelep) << endl;
   }
   if( -1.0e6*(thpos - thposp) < -20.0 ){
      cout << "e+ Underflow event (urad) " << (*label1) << " " << -1.0e6*thpos << " " << -1.0e6*thposp << " " << -1.0e6*(thpos - thposp) << endl;
   }   
   
   double phele = phi(fele);
   double phpos = phi(fpos);
   double phelep = phi(felep);
   double phposp = phi(fposp);   
   
   double zave= ((*z1) + (*z2))/2.0;
   double zvtx = 1.0e-3*zave;
   
   heleDeflectionPhi->Fill(1.0e3*(phele - phelep));  
   hposDeflectionPhi->Fill(1.0e3*(phpos - phposp));   
   
   heleDeflThetavsZ->Fill(zvtx, 1.0e6*(thele - thelep));  
   hposDeflThetavsZ->Fill(zvtx, -1.0e6*(thpos - thposp));

   heleDeflThetavsT->Fill(*label1,  1.0e6*(thele - thelep));  
   hposDeflThetavsT->Fill(*label1, -1.0e6*(thpos - thposp));

   
   heleDeflThetavsPhi->Fill(phele, 1.0e6*(thele - thelep));  
   hposDeflThetavsPhi->Fill(phpos,-1.0e6*(thpos - thposp));   
   
   heleDeflPhivsPhi->Fill(phele, 1.0e3*(phele - phelep));  
   hposDeflPhivsPhi->Fill(phpos, 1.0e3*(phpos - phposp));         

   return true;
}

void Ana::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Ana::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   
// Get a list of all objects in memory (for saving histograms)
   TList *list = gDirectory->GetList();
//   list->ls("-m");

// Save the histograms to a file - good for eg. fitting in the future
   TFile *fout = new TFile("Ana.root","RECREATE");

// We loop over all objects in memory and save any TH1 histograms. 
   TObject *obj;
   TIter iter(list);
   list->ls("-m");
// See https://root.cern.ch/root/roottalk/roottalk98/1015.html
   while( (obj = (TObject*) iter()) ){
        if( (obj->IsA()->InheritsFrom("TH1")) || (obj->IsA()->InheritsFrom("TH2")) || (obj->IsA()->InheritsFrom("TEfficiency")) ){
           obj->Draw();
           obj->Write();
        }
   }
   
   fout->Close();

}
