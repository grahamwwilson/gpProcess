
// Invoked with postlumi.sh bash script from ~/beamstats directory

void dlplot(std::string RUN){

   std::string fileName = "lumiee-"+RUN+".csv";
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

      {"E1", 'D'}, 
      {"E2", 'D'},
      {"x", 'D' }, 
      {"y", 'D'}, 
      {"z", 'D'},
      {"time", 'L'},
      {"x1p", 'D'}, 
      {"y1p", 'D'}, 
      {"x2p", 'D'}, 
      {"y2p", 'D'}, 
      {"s1x", 'O'},    // try to minimize size
      {"s1y", 'O'},
      {"s1z", 'O'},            
      {"s2x", 'O'}, 
      {"s2y", 'O'},
      {"s2z", 'O'},          
      {"label", 'L'}       
        } );
   
 
   std::string rootfile = "lumi-"+RUN+".root";
   
   cout << "rootfile " << rootfile << endl;
   
   df.Snapshot("t",rootfile.c_str());  

}
