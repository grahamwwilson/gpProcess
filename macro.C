void macro(){
 
// The usual command from inside the interpreter is 
// t->Process("Ana.C+");
// but we need to escape the double quotes.
// 
    gROOT->ProcessLine("t->Process(\"Ana.C+\")");

}
