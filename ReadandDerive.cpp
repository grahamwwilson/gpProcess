//
// Read lumi.ee output file from Guinea-PIG and 
// assemble information into a LumiEvent amd DerivedLumiEvent structs.
// Then superimpose a Bhabha event on the collision and create a corresponding BhabhaLumiEvent struct.
//
// Currently the super-imposed Bhabha event is a 2->2 event with the always forward scattered electron 
// at a polar angle of 31.3 mrad.
//
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

typedef std::mt19937 RandomNumberGenerator;

using namespace std;

// Make a vector of structs with luminosity events from a Guinea-PIG run

struct LumiEvent {

    double E1 {};   // Colliding electron energy [GeV]   (lab-frame)
    double E2 {};   // Colliding positron energy [GeV]   (lab-frame)
    double x {};    // Collision vertex in x [nm]        (lab-frame)
    double y {};    // Collision vertex in y [nm]        (lab-frame)
    double z {};    // Collision vertex in z [um]        (lab-frame)
    int t {};       // Discrete collision time. A non-negative integer related to the longitudinal motion of the grids. (lab-frame)
    double x1p {};  // Electron dx/ds [rad] (s is the path-length along the reference trajectory - in this case z)
    double y1p {};  // Electron dy/ds [rad]
    double x2p {};  // Positron dx/ds [rad]
    double y2p {};  // Positron dy/ds [rad]
    int s1x {};     // Electron spin-x      Note: All this spin info is not set
    int s1y {};     // Electron spin-y
    int s1z {};     // Electron spin-z
    int s2x {};     // Positron spin-x
    int s2y {};     // Positron spin-y
    int s2z {};     // Positron spin-z    
    int id {};      // The event ID which is tightly correlated with the interaction time

 // Specify a standard sorting criterion (here based on event id
    bool operator < (const LumiEvent & ev) const
    {
        return id < ev.id;
    }
    
};

struct DerivedLumiEvent : public LumiEvent{
 
    double ECM {};  // Center-of-mass energy [GeV]
    double betamag {}; // Speed of the lab-frame electron-positron colliding system (in units of c) - usually not very much
    std::array<double, 3> ubeta {};   // bx/b, by/b, bz/b (beta as a unit 3-vector).

// Constructor that initializes DerivedLumiEvent using an instance of LumiEvent
    DerivedLumiEvent(const LumiEvent& baseInstance, double ECMValue, double betamagValue, const std::array<double, 3>& betaValues ) 
        : LumiEvent(baseInstance), ECM(ECMValue), betamag(betamagValue), ubeta(betaValues) {}  // Initialize base with baseInstance
        
 // Specify a standard sorting criterion (here based on event id
    bool operator < (const DerivedLumiEvent & ev) const
    {
        return id < ev.id;
    }        

};


struct BhabhaLumiEvent : public LumiEvent{
 
    double ECM {};  // Center-of-mass energy [GeV]
    double betamag {}; // Speed of the lab-frame electron-positron colliding system (in units of c)
    std::array<double, 3> ubeta {}; // bx/b, by/b, bz/b (beta as a unit 3-vector).    
    std::array<double, 4> peCM {};  // 4-vector of the electron in the CoM frame (in units of GeV)
    std::array<double, 4> ppCM {};  // 4-vector of the positron in the CoM frame (in units of GeV)
    std::array<double, 4> pe {};    // 4-vector of the electron in the lab frame (in units of GeV)
    std::array<double, 4> pp {};    // 4-vector of the positron in the lab frame (in units of GeV)    

// Constructor that initializes BhabhaLumiEvent using an instance of LumiEvent
    BhabhaLumiEvent(const LumiEvent& baseInstance, double ECMValue, double betamagValue, 
                    const std::array<double, 3>& betaValues,
                    const std::array<double, 4>& peCMValues, const std::array<double, 4>& ppCMValues, 
                    const std::array<double, 4>& peValues, const std::array<double, 4>& ppValues ) 
        : LumiEvent(baseInstance), ECM(ECMValue), betamag(betamagValue), ubeta(betaValues), 
                                   peCM(peCMValues), ppCM(ppCMValues), pe(peValues), pp(ppValues) {}  // Initialize base with baseInstance
        
 // Specify a standard sorting criterion (here based on event id
    bool operator < (const BhabhaLumiEvent & ev) const
    {
        return id < ev.id;
    }        

};

std::vector<LumiEvent> vec;
std::vector<DerivedLumiEvent> vecd;
std::vector<BhabhaLumiEvent> vecb;

void Reader(int nevents, std::string infilename, unsigned long int seed){

    RandomNumberGenerator g(seed);
    std::uniform_real_distribution<double> uniform; 
    
    std::ifstream myfile;   
    
//    ifstream myfile("/home/graham/beamstats/lumi-LEPZ-105-100k.outfile");
//    ifstream myfile("/home/graham/beamstats/lumi-LEPZ-113.ee.out");
//    ifstream myfile("/home/graham/beamstats/lumi-Z-89.ee.out");    
    myfile.open(infilename);    
    
    TFile *f = new TFile("EMD-Analysis-Run113.root","RECREATE");
    TH1D* hECM = new TH1D("hECM","; Center-of-mass energy/Nominal; Events per 0.001 bin",750, 0.95, 1.025);    
    TH1D* hE1 = new TH1D("hE1","; Beam 1 E/Enominal; Events per 0.001 bin",750, 0.95, 1.025);
    TH1D* hE2 = new TH1D("hE2","; Beam 2 E/Enominal; Events per 0.001 bin",750, 0.95, 1.025);    
    TH1D* hx1p = new TH1D("hx1p","; Beam 1 x' [urad]; Events per 5 urad bin",240,-600.0,600.0);
    TH1D* hx2p = new TH1D("hx2p","; Beam 2 x' [urad]; Events per 5 urad bin",240,-600.0,600.0);
    TH1D* hy1p = new TH1D("hy1p","; Beam 1 y' [urad]; Events per 5 urad bin",240,-600.0,600.0);
    TH1D* hy2p = new TH1D("hy2p","; Beam 2 y' [urad]; Events per 5 urad bin",240,-600.0,600.0); 
    TH1D* hx = new TH1D("hx","; x [um]; Events per 8 um bin",200,-800.0,800.0);
    TH1D* hy = new TH1D("hy","; y [um]; Events per 0.2 um bin",200,-20.0,20.0);
    TH1D* hz = new TH1D("hz","; z [mm]; Events per 0.25 mm bin",400,-50.0,50.0);
    TH1D* hz2 = new TH1D("hz2","; z [um]; Events per 10 um bin",400,-2000.0,2000.0); 
    TH1D* ht0 = new TH1D("ht0","; t [discrete]; Events per bin",16400,-0.5,16399.5);       
    TH1D* ht1 = new TH1D("ht1","; t [discrete]; Events per bin",4100,-0.5,4099.5);
    TH1D* ht2 = new TH1D("ht2","; t [discrete]; Events per bin",410,-0.5,4099.5);         
    TH1D* hbeta = new TH1D("hbeta","; beta ; Events per 1e-5 bin",500,0.0,0.005);
    TH1D* htheta = new TH1D("htheta","; theta [mrad]; Events per 0.01 mrad bin",200,30.0,32.0);  
    TH1D* hphi = new TH1D("hphi","; phi [rad]; Events per 0.02",320,0.0,6.4);                            

    const double me = 0.51099895e-3;  // Electron mass in GeV
    const double Ebnominal = 45.6;    // Nominal beam energy [GeV]

    double E1, E2, x, y, z;
    int t;
    double x1p, y1p, x2p, y2p;
    int s1x, s1y, s1z, s2x, s2y, s2z;
    int id; 
    int iev = 0;   
    int NTOPRINT = 100;

    string line;

    if( myfile.is_open() ){
 
        while( getline( myfile, line ) && iev < nevents) 
        {

             stringstream ss(line);
             ss >> E1 >> E2 >> x >> y >> z >> t >> x1p >> y1p >> x2p >> y2p 
                >> s1x >> s1y >> s1z >> s2x >> s2y >> s2z >> id;
                
             iev +=1;
                  
             LumiEvent ev = { E1, E2, x, y, z, t, x1p, y1p, x2p, y2p, s1x, s1y, s1z, s2x, s2y, s2z, id };
             vec.push_back(ev);

// In practice we take E1 and E2 as the momenta of the colliding beams given that m << E.
             
             double p1[4];   // 4-vector of the electron (lab-frame z-axis aligned with nominal electron direction)
             double p2[4];   // 4-vector of the positron (lab-frame z-axis aligned with nominal electron direction)
             double p12[4];  // 4-vector of the electron-positron system (lab-frame z-axis aligned with nominal electron direction)

             p1[0] = sqrt(E1*E1 + me*me);
             p1[1] = E1*sin(atan2(x1p,1.0));
             p1[2] = E1*sin(atan2(y1p,1.0));
             p1[3] = sqrt(E1*E1 - p1[1]*p1[1] - p1[2]*p1[2]);
             
             p2[0] = sqrt(E2*E2 + me*me);
             p2[1] = E2*sin(atan2(x2p,1.0));
             p2[2] = E2*sin(atan2(y2p,1.0));
             p2[3] = -sqrt(E2*E2 - p2[1]*p2[1] - p2[2]*p2[2]);             
             
             std::array<double,3> beta;
             for (int i=0; i<=3; i++){
                 p12[i] = p1[i] + p2[i];
                 if ( i>=1){
                     beta[i-1] = p12[i]/p12[0]; 
                 }
             }
             double betamag = sqrt(beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2]);
             double ecm = sqrt(p12[0]*p12[0] - p12[1]*p12[1] - p12[2]*p12[2]  - p12[3]*p12[3]);            
             double gamma = 1.0/sqrt(1.0 - betamag*betamag);  // Lorentz gamma factor for the (small) boost from colliding beam (lab) frame to the CoM system 
             
             if ( iev <= NTOPRINT ){
                 cout << " " << endl;
                 cout << "Event " << iev << endl;
                 cout << "p1  " << p1[0] << " " << p1[1] << " " << p1[2] << " " << p1[3] << endl;
                 cout << "p2  " << p2[0] << " " << p2[1] << " " << p2[2] << " " << p2[3] << endl;             
                 cout << "p12 " << p12[0] << " " << p12[1] << " " << p12[2] << " " << p12[3] << endl;           
                 cout << "E1, E2, ECM, betamag, gamma, beta[0], beta[1], beta[2] " << E1 << " " << E2 << " " << ecm << " " <<  betamag << " " 
                      << gamma << " " <<  beta[0] << " " << beta[1] << " " << beta[2] << endl;
             }     
                  
             std::array<double, 3> b;
             b[0] = beta[0]/betamag;
             b[1] = beta[1]/betamag;
             b[2] = beta[2]/betamag; 
                                             
             double p1star[4]; double p2star[4]; // These general boosts are from the lab frame to the rest frame
             p1star[0] =         gamma*(p1[0]                       - beta[0]*p1[1]                       - beta[1]*p1[2]                        - beta[2]*p1[3]);
             p1star[1] = -gamma*beta[0]*p1[0] + (1.0 + (gamma-1.0)*b[0]*b[0])*p1[1] +         (gamma-1.0)*b[0]*b[1]*p1[2]  +         (gamma-1.0)*b[0]*b[2]*p1[3];
             p1star[2] = -gamma*beta[1]*p1[0] +         (gamma-1.0)*b[1]*b[0]*p1[1] + (1.0 + (gamma-1.0)*b[1]*b[1])*p1[2]  +         (gamma-1.0)*b[1]*b[2]*p1[3];
             p1star[3] = -gamma*beta[2]*p1[0] +         (gamma-1.0)*b[2]*b[0]*p1[1] +         (gamma-1.0)*b[2]*b[1]*p1[2]  + (1.0 + (gamma-1.0)*b[2]*b[2])*p1[3]; 
             
             p2star[0] =         gamma*(p2[0]                       - beta[0]*p2[1]                       - beta[1]*p2[2]                        - beta[2]*p2[3]);
             p2star[1] = -gamma*beta[0]*p2[0] + (1.0 + (gamma-1.0)*b[0]*b[0])*p2[1] +         (gamma-1.0)*b[0]*b[1]*p2[2]  +         (gamma-1.0)*b[0]*b[2]*p2[3];
             p2star[2] = -gamma*beta[1]*p2[0] +         (gamma-1.0)*b[1]*b[0]*p2[1] + (1.0 + (gamma-1.0)*b[1]*b[1])*p2[2]  +         (gamma-1.0)*b[1]*b[2]*p2[3];
             p2star[3] = -gamma*beta[2]*p2[0] +         (gamma-1.0)*b[2]*b[0]*p2[1] +         (gamma-1.0)*b[2]*b[1]*p2[2]  + (1.0 + (gamma-1.0)*b[2]*b[2])*p2[3];                            
             
// Now generate 4-vectors for Bhabha scattering in the rest frame that is moving with velocity beta with respect to the CoM system.
// Let's make the necessary unit 3-vectors for positioning our scattered electron on a cone of half-angle theta* with random phi* 
// with the cone centered on the boosted electron direction (p1star).

             // First form a unit vector along p1star
             double r[3];
             double p1cm = sqrt(p1star[1]*p1star[1] + p1star[2]*p1star[2] + p1star[3]*p1star[3]);
             r[0] = p1star[1]/p1cm;
             r[1] = p1star[2]/p1cm;
             r[2] = p1star[3]/p1cm;
             
             double smallest = abs(r[0]);
             int which = 0;
             if( abs(r[1]) < smallest ){
                 which = 1;
                 smallest = abs(r[1]);
             }
             if( abs(r[2]) < smallest ){
                 which = 2;
             }            
             double a[3] = {0.0, 0.0, 0.0};
             a[which] = 1.0;
             
             double u[3];
// Now define u = a x r which will be perpendicular to r
             u[0] = a[1]*r[2] - a[2]*r[1];
             u[1] = a[2]*r[0] - a[0]*r[2];
             u[2] = a[0]*r[1] - a[1]*r[0];   
             double umag = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
             u[0] = u[0]/umag;
             u[1] = u[1]/umag;
             u[2] = u[2]/umag;                                                    

             double v[3];
// Now define v = u x r which will be the axis perpendicular to the plane of u and r
             v[0] = u[1]*r[2] - u[2]*r[1];
             v[1] = u[2]*r[0] - u[0]*r[2];
             v[2] = u[0]*r[1] - u[1]*r[0];   
           
             if ( iev <= NTOPRINT ){
                 cout << "p1star  " << p1star[0] << " " << p1star[1] << " " << p1star[2] << " " << p1star[3] << endl;
                 cout << "p2star  " << p2star[0] << " " << p2star[1] << " " << p2star[2] << " " << p2star[3] << endl; 
                 cout << "b vector " << b[0] << " " << b[1] << " " << b[2] << " " << sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]) << endl;   
                 cout << "a vector " << a[0] << " " << a[1] << " " << a[2] << " " << sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) << endl;
                 cout << "u vector " << u[0] << " " << u[1] << " " << u[2] << " " << scientific << setprecision(12) << sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]) << endl;  
                 cout << "v vector magnitude " << scientific << setprecision(12) << sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) << endl;  
                 cout << "v vector " << v[0] << " " << v[1] << " " << v[2] << " " << scientific << setprecision(12) << sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) << endl; 
                 cout <<  "Summary " << endl;
                 cout << "r vector " << r[0] << " " << r[1] << " " << r[2] << endl;
                 cout << "u vector " << u[0] << " " << u[1] << " " << u[2] << endl;
                 cout << "v vector " << v[0] << " " << v[1] << " " << v[2] << endl;
             }
             
// So we now have our needed unit vectors.
             double thetaCM = 31.3e-3 ;      
// Choose random value of phiCM in [0, 2pi]
             double phiCM = 2.0*M_PI*uniform(g);
             hphi->Fill(phiCM);
             
             double evec[3];    // Unit vector on the cone
             for (int i=0; i<=2; i++){
                 evec[i] = sin(thetaCM)*(cos(phiCM)*u[i] + sin(phiCM)*v[i]) + cos(thetaCM)*r[i];
             }
             
// 4-vectors of the electron and positron in the CoM system       (Maybe use the sign of evec[2] to choose between forward and backward electron)
// For our purposes we impose that the electron is always forward scattered in the lab (ie. ignore the small s-channel contribution at these very forward angles)
             std::array<double,4> peCM; std::array<double,4> ppCM;
             peCM[0] = ecm/2.0;  ppCM[0] = peCM[0];
             double pcm = sqrt(peCM[0]*peCM[0] - me*me);
             if ( evec[2] > 0.0 ){
                 peCM[1] = pcm*evec[0];
                 peCM[2] = pcm*evec[1];
                 peCM[3] = pcm*evec[2];
             }
             else{
                 peCM[1] = -pcm*evec[0];
                 peCM[2] = -pcm*evec[1];
                 peCM[3] = -pcm*evec[2];                
             }
             ppCM[1] = -peCM[1];
             ppCM[2] = -peCM[2];
             ppCM[3] = -peCM[3];      
             
// And now we need to do a general Lorentz boost from the CM frame back to the lab (so in the -beta direction).
             std::array<double,4> peLab; std::array<double,4> ppLab; 
             
             peLab[0] =        gamma*(peCM[0] +                       beta[0]*peCM[1] +                       beta[1]*peCM[2]  +                       beta[2]*peCM[3]);
             peLab[1] = gamma*beta[0]*peCM[0] + (1.0 + (gamma-1.0)*b[0]*b[0])*peCM[1] +         (gamma-1.0)*b[0]*b[1]*peCM[2]  +         (gamma-1.0)*b[0]*b[2]*peCM[3];
             peLab[2] = gamma*beta[1]*peCM[0] +         (gamma-1.0)*b[1]*b[0]*peCM[1] + (1.0 + (gamma-1.0)*b[1]*b[1])*peCM[2]  +         (gamma-1.0)*b[1]*b[2]*peCM[3];
             peLab[3] = gamma*beta[2]*peCM[0] +         (gamma-1.0)*b[2]*b[0]*peCM[1] +         (gamma-1.0)*b[2]*b[1]*peCM[2]  + (1.0 + (gamma-1.0)*b[2]*b[2])*peCM[3]; 
             double peLabMag = sqrt(peLab[1]*peLab[1] + peLab[2]*peLab[2] + peLab[3]*peLab[3]);
             double theta = acos(peLab[3]/peLabMag);

             htheta->Fill(1.0e3*theta);
             
             ppLab[0] =        gamma*(ppCM[0] +                       beta[0]*ppCM[1] +                       beta[1]*ppCM[2]  +                       beta[2]*ppCM[3]);
             ppLab[1] = gamma*beta[0]*ppCM[0] + (1.0 + (gamma-1.0)*b[0]*b[0])*ppCM[1] +         (gamma-1.0)*b[0]*b[1]*ppCM[2]  +         (gamma-1.0)*b[0]*b[2]*ppCM[3];
             ppLab[2] = gamma*beta[1]*ppCM[0] +         (gamma-1.0)*b[1]*b[0]*ppCM[1] + (1.0 + (gamma-1.0)*b[1]*b[1])*ppCM[2]  +         (gamma-1.0)*b[1]*b[2]*ppCM[3];
             ppLab[3] = gamma*beta[2]*ppCM[0] +         (gamma-1.0)*b[2]*b[0]*ppCM[1] +         (gamma-1.0)*b[2]*b[1]*ppCM[2]  + (1.0 + (gamma-1.0)*b[2]*b[2])*ppCM[3];                    
             
             if ( iev <= NTOPRINT ){
                 cout << "peCM " << peCM[0] << " " << peCM[1] << " " << peCM[2] << " " << peCM[3] << endl;
                 cout << "ppCM " << ppCM[0] << " " << ppCM[1] << " " << ppCM[2] << " " << ppCM[3] << endl;     
                 cout << "peLab " << peLab[0] << " " << peLab[1] << " " << peLab[2] << " " << peLab[3] << " theta = " << fixed << 1.0e3*theta << " (mrad) " << endl;
                 cout << "ppLab " << scientific << ppLab[0] << " " << ppLab[1] << " " << ppLab[2] << " " << ppLab[3] << endl;             
                 cout << "psum  " << peLab[0]+ppLab[0] << " " << peLab[1]+ppLab[1] << " " << peLab[2]+ppLab[2] << " " << peLab[3]+ppLab[3] << endl; 
                 cout << "p12   " <<  p12[0] << " " << p12[1] << " " << p12[2] << " " << p12[3] << endl;
                 cout << "----------------------------------------------------------------------------------" << endl;
             }
                               
             hbeta->Fill(betamag);
             
             DerivedLumiEvent dev(ev, ecm, betamag, b);
             vecd.push_back(dev);
             
             BhabhaLumiEvent bha(ev, ecm, betamag, b, peCM, ppCM, peLab, ppLab);
             vecb.push_back(bha);             
             
             hECM->Fill(ecm/(2.0*Ebnominal)); hE1->Fill(E1/Ebnominal);  hE2->Fill(E2/Ebnominal);
             hx1p->Fill(1.0e6*x1p);   hx2p->Fill(1.0e6*x2p);  hy1p->Fill(1.0e6*y1p);  hy2p->Fill(1.0e6*y2p);
             hx->Fill(1.0e-3*x);  hy->Fill(1.0e-3*y);   hz->Fill(1.0e-3*z);  hz2->Fill(z);
             ht0->Fill(t); ht1->Fill(t); ht2->Fill(t);
             
        }
        myfile.close();
    }
    else{
        cout << "Unable to open file" << endl;
    }
    
//    std::sort(vecd.begin(), vecd.end() );        
    
// Print a header
    cout << " " << endl;
    cout << " DerivedLumiEvent structs potentially sorted by event id (smaller id is earlier in the collision) " << endl;
    cout << " " << endl;   
    cout << "        ID    ECM [GeV]     E1 [GeV]     E2 [GeV]     beta        bx/beta      by/beta      bz/beta     x [nm]       y [nm]       z [um]          t    x1p [rad]    y1p [rad]    x2p [rad]    y2p [rad] " << endl;
    cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    for (auto & el : vecd){
        cout << setw(10) << el.id << " " 
             << fixed << setprecision(6) << setw(12) << el.ECM << " "
             << fixed << setprecision(6) << setw(12) << el.E1 << " "  
             << fixed << setprecision(6) << setw(12) << el.E2 << " " 
             << scientific << setprecision(4) << setw(12) << el.betamag << " "  
             << scientific << setprecision(4) << setw(12) << el.ubeta[0] << " " 
             << scientific << setprecision(4) << setw(12) << el.ubeta[1] << " " 
             << scientific << setprecision(4) << setw(12) << el.ubeta[2] << " "             
             << scientific << setprecision(4) << setw(12) << el.x << " " 
             << scientific << setprecision(4) << setw(12) << el.y << " " 
             << scientific << setprecision(4) << setw(12) << el.z << " " 
             << scientific << setprecision(8) << setw(10) << el.t << " " 
             << scientific << setprecision(4) << setw(12) << el.x1p << " " 
             << scientific << setprecision(4) << setw(12) << el.y1p << " " 
             << scientific << setprecision(4) << setw(12) << el.x2p << " " 
             << scientific << setprecision(4) << setw(12) << el.y2p << endl;                                                               
    }
    
//    std::sort(vecb.begin(), vecb.end() );
    cout << " " << endl;
    cout << " BhabhaLumiEvent structs sorted by event id (smaller id is earlier in the collision).  " << endl;
    cout << " " << endl;   
    cout << "        ID    t           x [nm]          y [nm]         z [um]        E1          px1         py1         pz1        E2        px2        py2        pz2  " << endl;
    cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    for (auto & el : vecb){
        cout << setw(10) << el.id << " "
             << scientific << setprecision(8) << setw(10) << el.t << " "        
             << scientific << setprecision(6) << setw(14) << el.x << " " 
             << scientific << setprecision(6) << setw(14) << el.y << " " 
             << scientific << setprecision(6) << setw(14) << el.z << " " 
             << scientific << setprecision(6) << setw(14) << el.pe[0] << " " 
             << scientific << setprecision(6) << setw(14) << el.pe[1] << " " 
             << scientific << setprecision(6) << setw(14) << el.pe[2] << " " 
             << scientific << setprecision(6) << setw(14) << el.pe[3] << " " 
             << scientific << setprecision(6) << setw(14) << el.pp[0] << " " 
             << scientific << setprecision(6) << setw(14) << el.pp[1] << " " 
             << scientific << setprecision(6) << setw(14) << el.pp[2] << " " 
             << scientific << setprecision(6) << setw(14) << el.pp[3] << " "              
             << endl;                                                               
    }    
    
    f->Write();
    f->Close();   
}

int main(int argc, char** argv){

    CLI::App app{"Fill structs with GP lumi.ee.out luminosity events"};  
    
    int nevents = 1000;
    app.add_option("-n,--nevents", nevents, "Number of events to read for lumi file"); 
    
    std::string filename = "../lumi.ee_100k.outfile";
    app.add_option("-i,--ifile", filename, "Input data file"); 
    
    unsigned long int seed = 13579L;
    app.add_option("-s,--seed", seed, "Seed");     
    
    CLI11_PARSE(app, argc, argv);

    Reader(nevents, filename, seed);
       
    return 0;
    
}

