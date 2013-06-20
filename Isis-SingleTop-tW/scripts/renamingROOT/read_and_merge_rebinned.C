#include "Riostream.h" // dedicated ROOT iostream
#include <TFile>
#include <TH1F>
#include <TString>

void read_and_merge_rebinned(void){
    // based on http://root.cern.ch/root/html/tutorials/tree/basic.C.html
    // open text file. should have format: rootfile.root histonameinrootfile outfilename.root destinationhistname
    
    TString destname,sourcename,filename,outfilename;
    Int_t nlines=0;
    
    ifstream in;
    in.open("input_stat.txt"); // file with format rootfile.root histonameinrootfile outfilename.root destinationhistname
    
    TFile *infile;
    TFile *outfile;
    Double_t array10[11]={1,8,11,14,17,20,23,27,32,40,200};
    Double_t array5[6]={1,11,17,23,32,200};

    
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
	
	
	
	
	
        if(nlines < 550)
            std::cout << "line " << nlines << " read: " << filename << ":" << sourcename  << " --> " << outfilename << ":" << destname << std::endl;
        
        
	/*
	// binning for: 
	//10
	Int_t ngoodbins=10;
	Float_t integral = 0;

	Int_t thresholds5[6];
	Int_t thresholds10[11];
	Int_t ithr5=0;
	Int_t ithr10=0;
	
	for(int ibin=1; ibin<histoin->GetNbinsX()+1; ibin++){

	  cout << ibin << " " << histoin->Integral(1,ibin)/histoin->Integral(1,histoin->GetNbinsX()+1) << endl;
	  if(histoin->Integral(1,ibin)/histoin->Integral(1,histoin->GetNbinsX()+1)>(float)ithr5/5.0){
	    thresholds5[ithr5]=ibin;
	    ithr5++;
	  }
	  if(histoin->Integral(1,ibin)/histoin->Integral(1,histoin->GetNbinsX()+1)>(float)ithr10/10.0){
            thresholds10[ithr10]=ibin;
            ithr10++;
          }


	}
	cout << "optimal binning for 5 bins: " << endl;
	for(int ii=0; ii<5; ii++){
	 cout << thresholds5[ii] << endl;
	}
       cout << "optimal binning for 10 bins: " << endl;
        for(int ii=0; ii<10; ii++){
          cout << thresholds10[ii] << endl;
        }

	
	*/
         //histoin->SetName(destname);
          
	//TH1F* historebin = histoin->Rebin(10, destname, array10); 
	TH1F* historebin = histoin->Rebin(5, destname, array5); 
        outfile = new TFile(outfilename,"update");
        outfile->cd();
        historebin->SetDirectory(outfile);
        
        outfile->Write(destname,TFile::kOverwrite);
//        outfile->ls();
        outfile->Close();
        infile->Close();
        nlines++;
    }
    
}
