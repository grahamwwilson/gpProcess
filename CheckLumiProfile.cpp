//
// Read the part of the Guinea-PIG GPResults file that has the binned luminosity numbers like
// the following with a certain number of steps (computational parameter dependent)
// and the luminosity per step (total and peak). Peak is defined by the 
// ecm_min=0.99*2.0*energy.1; line of the acc.dat configuration file.
// Then put these in a couple of standard plots.
//
// For now I've edited the test input file from run Z-127 to be easily readable.
//
// lumi_ee=4.24398e+33;
// lumi_ee_high=4.18886e+33;
// step   lumi total   lumi peak 
// 1  1.24235e+25  1.24235e+25
// 2  6.04611e+25  6.04611e+25
// 3  9.44187e+25  9.44187e+25
// 4  1.96291e+26  1.96291e+26
// 5  3.64423e+26  3.64423e+26
// 6  5.79764e+26  5.79764e+26
#include "CLI11.hpp"   // See https://github.com/CLIUtils/CLI11    This is v2.4.2 of CLI11.
#include <iostream> 
#include <algorithm> //sort
#include <cmath>
#include <random>
#include <TRandom3.h>
#include <TMath.h>   //TMath::Prob
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <vector>    
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <array>


using namespace std;

// Make a vector of structs with luminosity steps from a Guinea-PIG output file

struct LumiStep {

    int step {};         // The step number (starts from 1)
    double lumiTot {};   // The total luminosity in the step (in more convenient units  TBD)
    double lumiPeak {};  // The total luminosity in the step (in more convenient units  TBD)

 // Specify a standard sorting criterion (here based on event id
    bool operator < (const LumiStep & st) const
    {
        return step < st.step;
    }
    
};

std::vector<LumiStep> v;

void Reader(std::string infilename, std::string rfilename){
    
    std::ifstream myfile;   
    myfile.open(infilename);    
    
    TFile *f = new TFile(rfilename.c_str(),"RECREATE");              

    const double LSCALEFACTOR=1.0e-31;          // multiplicative scale factor to inverse mb per step per BX

    int step;
    double lumiTot, lumiPeak;

    string line;

    if( myfile.is_open() ){
 
        while( getline( myfile, line )) 
        {
             stringstream ss(line);
             ss >> step >> lumiTot >> lumiPeak;
                  
             LumiStep st = { step, LSCALEFACTOR*lumiTot, LSCALEFACTOR*lumiPeak};
             v.push_back(st);
        }
        myfile.close();
    }
    else{
        cout << "Unable to open file" << endl;
    }
    
    int nbins = v.size();
    cout << "nbins = " << nbins << endl;
    TH1D* lTot  = new TH1D("lTot","One ILC Z Bunch Crossing; Computational Step Number; Total Lumi [/(mb*step*BX)]",nbins, 0.5, double(nbins)+0.5 );
    TH1D* lPeak = new TH1D("lPeak","One ILC Z Bunch Crossing; Computational Step Number; Peak Lumi [/(mb*step*BX)]",nbins, 0.5, double(nbins)+0.5 );  
    TH1D* lPeakTotRatio = new TH1D("lPeakTotRatio","One ILC Z Bunch Crossing; Computational Step Number; Peak Lumi / Total Lumi",nbins, 0.5, double(nbins)+0.5 );       
    
// Print a header
    cout << " " << endl;
    cout << " LumiStep struct. lumiTot and lumiPeak are in units of inverse milli-barn per computational step per BX" << endl;
    cout << " " << endl;   
    cout << "       step        lumiTot      lumiPeak  " << endl;
    cout << "------------------------------------------" << endl;
    double lumiSumTot = 0.0;
    double lumiSumPeak = 0.0;
    for (auto & el : v){
        cout << setw(10) << el.step << " " 
             << scientific << setprecision(5) << setw(12) << el.lumiTot << " "  
             << scientific << setprecision(5) << setw(12) << el.lumiPeak << endl;
        lTot->Fill(el.step, el.lumiTot); 
        lPeak->Fill(el.step, el.lumiPeak);
        lPeakTotRatio->Fill(el.step, el.lumiPeak/el.lumiTot);        
        lumiSumTot += el.lumiTot;
        lumiSumPeak += el.lumiPeak; 
    }
    cout << "Summed LumiTot = "  << scientific << setprecision(5) << setw(12) << lumiSumTot << endl;
    cout << "Summed LumiPeak = " << scientific << setprecision(5) << setw(12) << lumiSumPeak << endl;
    cout << "Ratio = " << fixed << setprecision(6) << setw(12) << lumiSumPeak/lumiSumTot << endl;
    double frep = 3.7;  // Repetition rate (Hz) - NB this is reduced from the canonical 5 Hz when running at the Z pole (see arXiv:1908.08212)
    double nb = 1312.0; // Number of bunches per bunch train
    cout << "Scaling up per BX lumi to an instantaneous luminosity of" << scientific << setprecision(5) << setw(12) << (lumiSumTot/LSCALEFACTOR)*1.0e-4*frep*nb << " cm^{-2} s^{-1} " << endl;
    cout << "The reference value (arXiv:1908.08212) is 2.05e33 cm^{-2} s^{-1} " << endl;
    
    f->Write();
    f->Close();   
}

int main(int argc, char** argv){

    CLI::App app{"Histogram binned luminosity information"};  
        
    std::string filename = "GPLumi-Summary.txt";
    app.add_option("-i,--ifile", filename, "Input data file (default: GPLumi-Summary.txt)"); 
    
    std::string rfilename = "GPLumi.root";
    app.add_option("-r,--rfile", rfilename, "Output ROOT file (default: GPLumi.root)");     
    
    CLI11_PARSE(app, argc, argv);

    Reader(filename, rfilename);
       
    return 0;
    
}

