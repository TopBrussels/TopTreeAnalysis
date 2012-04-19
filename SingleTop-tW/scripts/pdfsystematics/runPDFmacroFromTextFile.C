{
  // should be run on an lxplus-like system with a version of CMSSW set up!!

  // Make sure you link or obtain the LHAPDF functions. 
  // for the top group, also see: https://twiki.cern.ch/twiki/bin/view/CMS/TWikiTopRefSyst#PDF_uncertainties
 
//  ln -s /afs/cern.ch/cms/slc5_amd64_gcc434/external/lhapdf/5.6.0-cms4/lib/libLHAPDF*.so .                                                                                                   
// ln -s /afs/cern.ch/cms/slc5_amd64_gcc434/external/lhapdf/5.6.0-cms4/include/LHAPDF LHAPDF                                                                                                  
// ln -s /afs/cern.ch/cms/slc5_amd64_gcc434/external/lhapdf/5.6.0-cms4/share/lhapdf/PDFsets PDFsets        

// and EVERY TIME you log in (equivalent to cmsenv), set up the workspace by issueing the following command:
// scram setup lhapdffull

  // setting up fwlight, actually probably not necessary
  gSystem->Load("libFWCoreFWLite.so"); 
  AutoLibraryLoader::enable();
  gSystem->Load("libDataFormatsFWLite.so");
  gROOT->ProcessLine("namespace edm {typedef edm::Wrapper<vector<float> > Wrapper<vector<float,allocator<float> > >; }");
  gROOT->ProcessLine("namespace edm {typedef edm::Wrapper<vector<double> > Wrapper<vector<double,allocator<double> > >; }");


  // load the PDF libs
  //  gSystem->Load("libLHAPDFWrap.so"); // old version
  gSystem->Load("libLHAPDF.so"); // new version with namespace...


  // this is the macro where the magic happens, it is a MakeClass macro with some minor modifications that loads Rebeca's small trees and uses them to calculate pdf uncertainties.

  gSystem->CompileMacro("pdfMacro.C","k");

  // different files are read in ... NO SECURITIES for array sizes!!!!   
  TString filenames[45];
filenames[0]="forPDF/txt/nocommas_pdf_2j1t_0_atw_dr.txt";
filenames[1]="forPDF/txt/nocommas_pdf_2j1t_0_atw_ds.txt";
filenames[2]="forPDF/txt/nocommas_pdf_2j1t_0_tt.txt";
filenames[3]="forPDF/txt/nocommas_pdf_2j1t_0_tw_dr.txt";
filenames[4]="forPDF/txt/nocommas_pdf_2j1t_0_tw_ds.txt";
filenames[5]="forPDF/txt/nocommas_pdf_2j1t_1_atw_dr.txt";
filenames[6]="forPDF/txt/nocommas_pdf_2j1t_1_atw_ds.txt";
filenames[7]="forPDF/txt/nocommas_pdf_2j1t_1_tt.txt";
filenames[8]="forPDF/txt/nocommas_pdf_2j1t_1_tw_dr.txt";
filenames[9]="forPDF/txt/nocommas_pdf_2j1t_1_tw_ds.txt";
filenames[10]="forPDF/txt/nocommas_pdf_2j1t_2_atw_dr.txt";
filenames[11]="forPDF/txt/nocommas_pdf_2j1t_2_atw_ds.txt";
filenames[12]="forPDF/txt/nocommas_pdf_2j1t_2_tt.txt";
filenames[13]="forPDF/txt/nocommas_pdf_2j1t_2_tw_dr.txt";
filenames[14]="forPDF/txt/nocommas_pdf_2j1t_2_tw_ds.txt";
filenames[15]="forPDF/txt/nocommas_pdf_2j2t_0_atw_dr.txt";
filenames[16]="forPDF/txt/nocommas_pdf_2j2t_0_atw_ds.txt";
filenames[17]="forPDF/txt/nocommas_pdf_2j2t_0_tt.txt";
filenames[18]="forPDF/txt/nocommas_pdf_2j2t_0_tw_dr.txt";
filenames[19]="forPDF/txt/nocommas_pdf_2j2t_0_tw_ds.txt";
filenames[20]="forPDF/txt/nocommas_pdf_2j2t_1_atw_dr.txt";
filenames[21]="forPDF/txt/nocommas_pdf_2j2t_1_atw_ds.txt";
filenames[22]="forPDF/txt/nocommas_pdf_2j2t_1_tt.txt";
filenames[23]="forPDF/txt/nocommas_pdf_2j2t_1_tw_dr.txt";
filenames[24]="forPDF/txt/nocommas_pdf_2j2t_1_tw_ds.txt";
filenames[25]="forPDF/txt/nocommas_pdf_2j2t_2_atw_dr.txt";
filenames[26]="forPDF/txt/nocommas_pdf_2j2t_2_atw_ds.txt";
filenames[27]="forPDF/txt/nocommas_pdf_2j2t_2_tt.txt";
filenames[28]="forPDF/txt/nocommas_pdf_2j2t_2_tw_dr.txt";
filenames[29]="forPDF/txt/nocommas_pdf_2j2t_2_tw_ds.txt";
filenames[30]="forPDF/txt/nocommas_pdf_signal_0_atw_dr.txt";
filenames[31]="forPDF/txt/nocommas_pdf_signal_0_atw_ds.txt";
filenames[32]="forPDF/txt/nocommas_pdf_signal_0_tt.txt";
filenames[33]="forPDF/txt/nocommas_pdf_signal_0_tw_dr.txt";
filenames[34]="forPDF/txt/nocommas_pdf_signal_0_tw_ds.txt";
filenames[35]="forPDF/txt/nocommas_pdf_signal_1_atw_dr.txt";
filenames[36]="forPDF/txt/nocommas_pdf_signal_1_atw_ds.txt";
filenames[37]="forPDF/txt/nocommas_pdf_signal_1_tt.txt";
filenames[38]="forPDF/txt/nocommas_pdf_signal_1_tw_dr.txt";
filenames[39]="forPDF/txt/nocommas_pdf_signal_1_tw_ds.txt";
filenames[40]="forPDF/txt/nocommas_pdf_signal_2_atw_dr.txt";
filenames[41]="forPDF/txt/nocommas_pdf_signal_2_atw_ds.txt";
filenames[42]="forPDF/txt/nocommas_pdf_signal_2_tt.txt";
filenames[43]="forPDF/txt/nocommas_pdf_signal_2_tw_dr.txt";
filenames[44]="forPDF/txt/nocommas_pdf_signal_2_tw_ds.txt";

  // different pdfs are considered... NO SECURITIES for array sizes!!!!
  TString sets[4];
  sets[1]="MRST2006nnlo.LHgrid";
  sets[0]="cteq66.LHgrid";
  sets[2]="NNPDF10_100.LHgrid";
  sets[3]="MRST2007lomod.LHgrid";

  // loop over samples
  for(size_t isam=0; isam<45; isam++){
    TChain *ch = new TChain("myTree","myTree");
    //    ch->Add(filenames[isam]);

    // and over pdfss
    for(size_t iset=0; iset<1; iset++){
      // some linkin
      pdfMacro macro(ch);
      macro.SetPdfSetName(sets[iset]);
      //      macro.Loop();
      macro.ReadFromTextFile(filenames[isam]);
      std::cout << "looking at sample "<< filenames[isam] << std::endl;
    }
    std::cout << "DONE with sample "<< filenames[isam] << std::endl;
  }
  std::cout << "exiting..." << std::endl;
  gApplication->Terminate();
}
