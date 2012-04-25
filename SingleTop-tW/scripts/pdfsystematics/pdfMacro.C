#define pdfMacro_cxx
#include "pdfMacro.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
// the following was useable in CMSSW426
#include "LHAPDF/LHAPDF.h"
#include "Riostream.h"


void pdfMacro::ReadFromTextFile(TString textfile,float scaler){

  bool extravars=true;
//   In a ROOT session, you can do:
//      Root > .L pdfMacro.C
//      Root > pdfMacro t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   LHAPDF::initLHAPDF();
   std::string pdfname="MRST2006nnlo.LHgrid";
   std::vector<LHAPDF::PDFSetInfo> pdfinfo = LHAPDF::getAllPDFSetInfo();

   std::vector<int> interestingpdfs;
   std::cout << "loaded: " << pdfinfo.size() << " pdfs." << std::endl;
   for(size_t ii=0; ii<pdfinfo.size(); ii++){
     
     //     std::cout << "pdf " << ii << " name " << pdfinfo[ii].file << "  " << pdfinfo[ii].description << std::endl;
     //     if(pdfinfo[ii].file=="MRST2006nnlo.LHgrid")
     //       interestingpdfs.push_back(ii);
     if(pdfinfo[ii].file==pdfSet_){

       interestingpdfs.push_back(ii);
       //       break;
     }
   }
   //   for(size_t jj=0; jj<interestingpdfs.size(); jj++)
   LHAPDF::initPDFSet(pdfinfo[interestingpdfs[0]].file);
   Long64_t nbytes = 0, nb = 0;

   
   std::vector<float> baseweight;
   std::vector<float> myx1;
   std::vector<float> myx2;
   std::vector<float> myq;
   std::vector<int>   myid1;
   std::vector<int>   myid2;
   std::vector<float> myvar1;
   std::vector<float> myvar2;
   std::vector<std::vector<float> > weights;

 
   //open an asci file:
   ifstream instream;
   TString dirpath = "forPDF/txt/";
   TString filename = "nocommas_pdf_2j1t_0_atw_dr.txt";

   // Rebeca's default files need to be cleaned of commas to make for easy parsing. I executed the following command (once):
   // ls *.txt | awk '{print "cat "$1" | sed -e s/,//g > nocommas_"$1}' | sh

   std::cout << "SUMMARY: ++++++++++++++++++++++++++++++++++++++" << std::endl;
   if(textfile.Contains(".txt")){
     instream.open(textfile.Data());
     std::cout << "SUMMARY: opening file " << textfile << std::endl;
   }
   else{
     instream.open((dirpath+filename).Data());
     std::cout << "SUMMARY: opening file " << (dirpath+filename) << std::endl;
   }
   Float_t weight;
   TString name;
   TString mode;
   Float_t var1, var2;
   while(instream.good()){


     instream >> weight >>  x1 >> x2 >> q >> id1 >> id2 ;
     if(extravars)
       instream >> var1 >> var2 ;
     instream >> name >> mode ;

     //     std::cout << weight << ", " << x1 << ", " << x2 << ", " << q << ", " << id1 << ", " << id2 << ", " << name << ", " << mode << ", " << var1 << " " << var2 << std::endl;
     if(!instream.good())
       break;
     
     myx1.push_back(x1);
     myx2.push_back(x2);
     myq.push_back(q);
     myid1.push_back(id1);
     myid2.push_back(id2);
     myvar1.push_back(var1);
     myvar2.push_back(var2);
     baseweight.push_back(weight);
     weights.push_back(std::vector<float>(45,0));

   }
   // from appendix in AN-2009-048 (pdf note):
   //      LHAPDF::usePDFMember(interestingpdfs[0]);
   std::vector<float> totalsums;
   float minval=0;
   float maxval=0;
   float quadratmin=0;
   float quadratmax=0;
   float minmean=0;
   float maxmean=0;
   float nmin=0;
   float nmax=0;
   std::cout << "examining PDF: " <<interestingpdfs[0] << ", having " << LHAPDF::numberPDF() << " subsets." << std::endl;
   for(int imem=0; imem<LHAPDF::numberPDF(); imem++){
     float totsum=0;
     totalsums.push_back(totsum);
     LHAPDF::usePDFMember(imem);
     for(size_t iev=0; iev<weights.size(); iev++){
       weights[iev][imem]=baseweight[iev]*LHAPDF::xfx(myx1[iev],myq[iev],myid1[iev])* LHAPDF::xfx(myx2[iev],myq[iev],myid2[iev]);
       totalsums[imem]+=weights[iev][imem];
     }
     if(totalsums[imem]>0){
       //       std::cout << "imem: " << imem << " total sum events: " << totalsums[imem] << " " << 100*(totalsums[imem]-totalsums[0])/totalsums[0] << "%";
       if((totalsums[imem]-totalsums[0])/totalsums[0]<0){
	 minmean+=(totalsums[imem]-totalsums[0])/totalsums[0];
	 quadratmin+=pow((totalsums[imem]-totalsums[0])/totalsums[0],2);
	 nmin++;
       }
       else{
	 maxmean+=(totalsums[imem]-totalsums[0])/totalsums[0];
	 quadratmax+=pow((totalsums[imem]-totalsums[0])/totalsums[0],2);
	 nmax++;
       }
       if(minval>(totalsums[imem]-totalsums[0])/totalsums[0]){
	 //	 cout << "\t ---";
	 minval=(totalsums[imem]-totalsums[0])/totalsums[0];
	
       }
        if(maxval<(totalsums[imem]-totalsums[0])/totalsums[0]){
	  //	 cout << "\t +++";
	 maxval=(totalsums[imem]-totalsums[0])/totalsums[0];
       }
	//       std::cout<< std::endl;

     }
   }

   if(extravars){
     TH1F *histo1 = new TH1F("histo1","histo1",20,20,60);
     histo1->Sumw2();
     histo1->SetXTitle("var1");
     histo1->SetYTitle("a.u.");
     TH1F *histo1a = (TH1F*) histo1->Clone("histo1a");
     TH1F *histo1b = (TH1F*) histo1->Clone("histo1b");

     TH1F *histo2 = new TH1F("histo2","histo2",20,100,300);
     histo2->Sumw2();
     histo2->SetXTitle("var1");
     histo2->SetYTitle("a.u.");
     TH1F *histo2a = (TH1F*) histo2->Clone("histo2a");
     TH1F *histo2b = (TH1F*) histo2->Clone("histo2b");  
     
     float lowest, highest;
     for(size_t iev=0; iev<weights.size(); iev++){
       lowest=highest=weights[iev][0];
       for(int imem=1; imem<LHAPDF::numberPDF(); imem++){
	 if(lowest>weights[iev][imem])
	   lowest=weights[iev][imem];
	 if(highest<weights[iev][imem])
	   highest=weights[iev][imem];
       }
	
       histo1->Fill(myvar1[iev],weights[iev][0]);
       histo1a->Fill(myvar1[iev],lowest);
       histo1b->Fill(myvar1[iev],highest);

       histo2->Fill(myvar2[iev],weights[iev][0]);
       histo2a->Fill(myvar2[iev],lowest);
       histo2b->Fill(myvar2[iev],highest);       
     }
     std::cout << "Kolmogorov, lowest: " << histo1a->KolmogorovTest(histo1) << std::endl;
     std::cout << "Kolmogorov, highest: " << histo1b->KolmogorovTest(histo1) << std::endl;
     std::cout << "Kolmogorov, lowest: " << histo2a->KolmogorovTest(histo2) << std::endl;
     std::cout << "Kolmogorov, highest: " << histo2b->KolmogorovTest(histo2) << std::endl;


     TCanvas *canv = new TCanvas();
     canv->cd();
     TString workfile = textfile;
     workfile.ReplaceAll(".txt","_var1.jpg");
     histo1->SetTitle("var 1");
     histo2->SetTitle("var 2");
     histo1->Draw("hist");
     histo1a->SetLineStyle(2);
     histo1b->SetLineStyle(2);
     histo1a->Draw("histsame");
     histo1b->Draw("histsame");
     if(histo1a->GetMaximum()>histo1->GetMaximum())
       histo1->SetMaximum(histo1a->GetMaximum()*1.05);
     if(histo1b->GetMaximum()>histo1->GetMaximum())
       histo1->SetMaximum(histo1b->GetMaximum()*1.05);
     if(histo1a->GetMinimum()<histo1->GetMinimum())
       histo1->SetMinimum(histo1a->GetMinimum()*1.05);
     if(histo1b->GetMinimum()<histo1->GetMinimum())
       histo1->SetMinimum(histo1b->GetMinimum()*1.05);
     canv->Update();
     canv->Print(workfile);

     workfile.ReplaceAll("_var1.jpg","_var2.jpg");
     histo2->Draw("hist");
     histo2a->SetLineStyle(2);
     histo2b->SetLineStyle(2);
     histo2a->Draw("histsame");
     histo2b->Draw("histsame");
     if(histo2a->GetMaximum()>histo2->GetMaximum())
       histo2->SetMaximum(histo2a->GetMaximum()*1.05);
     if(histo2b->GetMaximum()>histo2->GetMaximum())
       histo2->SetMaximum(histo2b->GetMaximum()*1.05);
     if(histo2a->GetMinimum()<histo2->GetMinimum())
       histo2->SetMinimum(histo2a->GetMinimum()*1.05);
     if(histo2b->GetMinimum()<histo2->GetMinimum())
       histo2->SetMinimum(histo2b->GetMinimum()*1.05);
     canv->Update();
     canv->Print(workfile);
   }
   minval*=scaler;
   maxval*=scaler;
   minmean*=scaler;
   maxmean*=scaler;
   std::cout << std::setprecision(3) << "SUMMARY for pdf set: " << pdfSet_ << ":" << std::endl;
   std::cout << std::setprecision(3) << "SUMMARY (worse values): " << minval << " " << maxval << "" << std::endl;
   //   std::cout << std::setprecision(3) << "SUMMARY (mean values): " << minmean/nmin << " " << maxmean/nmax << "" << std::endl;
   std::cout << std::setprecision(3)<< "SUMMARY (all in quadrature): " << -sqrt(quadratmin)*scaler << " " << sqrt(quadratmax)*scaler <<"" << std::endl;
   //   std::cout << std::setprecision(3) << "SUMMARY (all in quadrature (mean)): " << -sqrt(quadratmin)/nmin*scaler << " " << sqrt(quadratmax)/nmax*scaler <<"" << std::endl;
   std::cout << "SUMMARY: ++++++++++++++++++++++++++++++++++++++" << std::endl;
}

void pdfMacro::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L pdfMacro.C
//      Root > pdfMacro t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   LHAPDF::initLHAPDF();
   std::string pdfname="MRST2006nnlo.LHgrid";
   std::vector<LHAPDF::PDFSetInfo> pdfinfo = LHAPDF::getAllPDFSetInfo();

   std::vector<int> interestingpdfs;
   std::cout << "loaded: " << pdfinfo.size() << " pdfs." << std::endl;
   for(size_t ii=0; ii<pdfinfo.size(); ii++){
     
     //     std::cout << "pdf " << ii << " name " << pdfinfo[ii].file << "  " << pdfinfo[ii].description << std::endl;
     //     if(pdfinfo[ii].file=="MRST2006nnlo.LHgrid")
     //       interestingpdfs.push_back(ii);
     if(pdfinfo[ii].file==pdfSet_){

       interestingpdfs.push_back(ii);
       //       break;
     }
   }
   //   for(size_t jj=0; jj<interestingpdfs.size(); jj++)
   LHAPDF::initPDFSet(pdfinfo[interestingpdfs[0]].file);
   Long64_t nbytes = 0, nb = 0;

   
      
   std::vector<float> myx1;
   std::vector<float> myx2;
   std::vector<float> myq;
   std::vector<int>   myid1;
   std::vector<int>   myid2;
   std::vector<std::vector<float> > weights;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      // don't care about doing full selection, just calculate pdfs


      myx1.push_back(x1);
      myx2.push_back(x2);
      myq.push_back(q);
      myid1.push_back(id1);
      myid2.push_back(id2);

      weights.push_back(std::vector<float>(45,0));

   }
   // from appendix in AN-2009-048 (pdf note):
   //      LHAPDF::usePDFMember(interestingpdfs[0]);
   std::vector<float> totalsums;
   float minval=0;
   float maxval=0;
   float quadratmin=0;
   float quadratmax=0;
   std::cout << "examining PDF: " <<interestingpdfs[0] << ", having " << LHAPDF::numberPDF() << " subsets." << std::endl;
   for(int imem=0; imem<LHAPDF::numberPDF(); imem++){
     float totsum=0;
     totalsums.push_back(totsum);
     LHAPDF::usePDFMember(imem);
     for(size_t iev=0; iev<weights.size(); iev++){
       weights[iev][imem]=LHAPDF::xfx(myx1[iev],myq[iev],myid1[iev])* LHAPDF::xfx(myx2[iev],myq[iev],myid2[iev]);
       totalsums[imem]+=weights[iev][imem];
     }
     if(totalsums[imem]>0){
       //       std::cout << "imem: " << imem << " total sum events: " << totalsums[imem] << " " << 100*(totalsums[imem]-totalsums[0])/totalsums[0] << "%";
       if(minval>100*(totalsums[imem]-totalsums[0])/totalsums[0]){
	 //	 cout << "\t ---";
	 minval=100*(totalsums[imem]-totalsums[0])/totalsums[0];
	 quadratmin+=pow((totalsums[imem]-totalsums[0])/totalsums[0],2);
       }
        if(maxval<100*(totalsums[imem]-totalsums[0])/totalsums[0]){
	  //	 cout << "\t +++";
	 maxval=100*(totalsums[imem]-totalsums[0])/totalsums[0];
	 quadratmax+=pow((totalsums[imem]-totalsums[0])/totalsums[0],2);
       }
	//       std::cout<< std::endl;

     }
   }
   std::cout << "SUMMARY for pdf set: " << pdfSet_ << ":" << std::endl;
   std::cout << "SUMMARY (worse values): " << minval << "% " << maxval << "%" << std::endl;
   std::cout << "SUMMARY (all in quadrature): " << -100*sqrt(quadratmin) << "% " << 100*sqrt(quadratmax) <<"%" << std::endl;
}


