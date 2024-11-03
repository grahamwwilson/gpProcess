//using namespace ROOT;

void dfplot(std::string RUN){

   std::string fileName = "qcombpairs-"+RUN+".csv";
   cout << "fileName = " << fileName << endl;

/*
Parameters
    [in]	fileName	Path of the CSV file.
    [in]	readHeaders	true if the CSV file contains headers as first row, false otherwise (default true).
    [in]	delimiter	Delimiter character (default ',').
    [in]	linesChunkSize	bunch of lines to read, use -1 to read all
    [in]	colTypes	Allow user to specify custom column types, accepts an unordered map with keys being column type, values being type alias ('O' for boolean, 'D' for double, 'L' for Long64_t, 'T' for std::string)

Definition at line 558 of file RCsvDS.cxx.
*/

//   std::unordered_map<std::string, char> && colTypes = { {"z1", 'D' }, {"z2", 'D'}};

   auto df = ROOT::RDF::FromCSV(fileName, true, ' ', -1LL, {  

      {"x1", 'D'}, 
      {"y1", 'D'},
      {"z1", 'D' }, 
      {"x2", 'D'}, 
      {"y2", 'D'},
      {"z2", 'D'},
      {"x1p", 'D'}, 
      {"y1p", 'D'}, 
      {"z1p", 'D'}, 
      {"x2p", 'D'}, 
      {"y2p", 'D'}, 
      {"z2p", 'D'}       
        } );
   
// Define a new column for center-of-mass energy
//   auto funECM = [](double x, double y) { return 2.0*sqrt(x*y); };
//   auto ECMMean = df.Define("ECM", funECM, {"E1", "E2"}).Mean("ECM");
//   std::cout << *ECMMean << std::endl;
   
//   df.Define("ECM", "2.0*sqrt(E1*E2)");
    
//   auto df_all = df.Filter("E1 > -100.0", "All Events");

/*   
   auto hE = df.Histo1D( {"hE", "ILC250; E (GeV); N_{Events}", 1300, 0.0, 130.0}, "E"); 
   auto hx = df.Histo1D( {"hx", "ILC250; x (um); N_{Events}", 500, -2.0, 2.0}, "x");
   auto hy = df.Histo1D( {"hy", "ILC250; y (um); N_{Events}", 500, -0.05, 0.05}, "y");
   auto hxp = df.Histo1D( {"hxp", "ILC250; xp (urad); N_{Events}", 320, -160.0, 160.0}, "xp");
   auto hyp = df.Histo1D( {"hyp", "ILC250; yp (urad); N_{Events}", 300, -100.0, 100.0}, "yp");
   auto hz = df.Histo1D( {"hz", "ILC250; z (um); N_{Events}", 500, -1000.0, 1000.0}, "z");

   hE->DrawClone("hist"); 
   hz->DrawClone("hist");
   hx->DrawClone("hist");
   hy->DrawClone("hist");
   hxp->DrawClone("hist");
   hyp->DrawClone("hist");
      
// Can we do profiles too?
   auto hproftz = df.Profile1D({"hproftz", "Profile of z versus t", 65, -0.5, 5199.5}, "t", "z");
   auto hprofE1t = df.Profile1D({"hprofE1t", "Profile of E_{1} versus t", 65, -0.5, 5199.5}, "t", "E1");
   auto hprofE1z = df.Profile1D({"hprofE1z", "Profile of E_{1} versus z; z (microns); <E_{1}> (GeV)", 40, -1000.0, 1000.0}, "z", "E1");
   auto hprofE2z = df.Profile1D({"hprofE2z", "Profile of E_{2} versus z; z (microns); <E_{2}> (GeV)", 80, -1000.0, 1000.0}, "z", "E2");
   auto hprofECMz = df.Profile1D({"hprofECMz", "Profile of ECM versus z; z (microns); <ECM> (GeV)", 40, -1000.0, 1000.0}, "z", "ECM");
   auto hprofEdiffz = df.Profile1D({"hprofEdiffz", "Profile of Ediff versus z; z (microns); <Ediff> (GeV)", 40, -1000.0, 1000.0}, "z", "Ediff");
   auto hprofECMt = df.Profile1D({"hprofECMt", "Profile of ECM versus t; t (units?); <ECM> (GeV)", 65, -0.5, 5199.5}, "t", "ECM");
   auto hprofEdifft = df.Profile1D({"hprofEdifft", "Profile of Ediff versus t; t (units?); <Ediff> (GeV)", 65, -0.5, 5199.5}, "t", "Ediff");
   
// 2d plots   

// Statistics?
   auto allCutsReport = df.Report();
   allCutsReport->Print();   
   
   hproftz->DrawClone();
   hprofE1t->DrawClone();
   hprofE1z->DrawClone();
   hprofE2z->DrawClone();
   hprofECMz->DrawClone();
   hprofEdiffz->DrawClone();
   hprofECMt->DrawClone();
   hprofEdifft->DrawClone();
*/
 
   std::string rootfile = RUN+".root";
   
   cout << "rootfile " << rootfile << endl;
   
   df.Snapshot("t",rootfile.c_str());  

}
