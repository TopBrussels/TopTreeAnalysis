#include "Riostream.h" // dedicated ROOT iostream
#include <TFile>
#include <TH1F>
#include <TString>

void read_and_merge(void){
    // based on http://root.cern.ch/root/html/tutorials/tree/basic.C.html
    // open text file. should have format: rootfile.root histonameinrootfile outfilename.root destinationhistname
    
    TString destname,sourcename,filename,outfilename;
    Int_t nlines=0;
    
    ifstream in;
    in.open("input.txt"); // file with format rootfile.root histonameinrootfile outfilename.root destinationhistname
    
    TFile *infile;
    TFile *outfile;
    
    while(1){
        in >> filename >> sourcename >> outfilename >> destname;
        if(!in.good()) break;
        
        infile = new TFile(filename,"read");
        //    infile->ls();
        TH1F *histoin = (TH1F*)infile->Get(sourcename);
        if(histoin==0) {
            cout << "tried to get " << sourcename << " from file " << filename << " and failed - check your input file! continueing without copy." << endl;
            continue;
        }
        if(nlines < 10)
            std::cout << "line " << nlines << " read: " << filename << ":" << sourcename  << " --> " << outfilename << ":" << destname << std::endl;
        
        
        histoin->SetName(destname);
        outfile = new TFile(outfilename,"update");
        outfile->cd();
        histoin->SetDirectory(outfile);
        
        outfile->Write(destname,TFile::kOverwrite);
//        outfile->ls();
        outfile->Close();
        infile->Close();
        nlines++;
    }
    
}
