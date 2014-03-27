#include "Headers.h"

bool inclusive=true;

using namespace std;
using namespace TopTree;
using namespace reweight;

typedef map<TString,TH1F*>      map_TH1F;  // <sample name , histo>
typedef map<TString,map_TH1F>   map2_TH1F; // <variable name , set of (1 histo/sample)>
typedef map<TString,map2_TH1F>  map3_TH1F; // <selection name , set of <var,set histos> >

///////////////////////
// DECLARE FUNCTIONS //
///////////////////////

int plot(TString, float, float);

//////////
// MAIN //
//////////

int main(int argc, char *argv[])
{

  TString path="test";

  float magnify=1.3;
  float magnifyLog=1.5;

  if(argc>1) path      = argv[1];

  plot(path, magnify, magnifyLog);

  return 0;

}

//////////////
// ANALYSIS //
//////////////

int plot(TString path="Results/forgot/", float magnify=1.3, float magnifyLog=1.5)
{

  /////////////////
  // ENVIRONMENT //
  /////////////////
  
  // Generic stuff
  clock_t start = clock();
  int verbose = 3;
  string pathPNG = path.Data();
  mkdir(pathPNG.c_str(),0777);

  // Outlog
  ofstream outlog(path+"/log.txt",ios::out);

  // TLegend coordinates
  float x1=0.55; float y1=0.66; float x2=0.88; float y2=0.88;

  // Input ROOT file
  TString finname = path+"/results.root";
  TFile *fin = new TFile(finname, "READ");
  if(fin->IsZombie()) {
    cout << "ERROR !!! FILE " << finname << " IS ZOMBIE !!!" << endl;
    return 3;
  }
  
  // Output ROOT file
  TString foutname = path+"/plots.root";
  cout << "foutname=" << foutname << endl;
  TFile *fout = new TFile (foutname, "RECREATE");

  ///////////////////////////////
  // GET PLOTS FROM INPUT FILE //
  ///////////////////////////////

  map3_TH1F AllPlots;

  const int nSel=22;
  const int nVar=1;
  const int nSamp=4;

  TString nomSel[nSel]={"AnaCut_none",  "AnaCut_PV",
			"AnaCut_2jets", "AnaCut_eta", "AnaCut_dphi", "AnaCut_pt1", "AnaCut_pt2", "AnaCut_kine",
			"AnaCut_jet1_emf0p2","AnaCut_jet2_emf0p2","AnaCut_pair_emf0p2","AnaCut_both_emf0p2",
			"AnaCut_jet1_chf0p1","AnaCut_jet2_chf0p1","AnaCut_pair_chf0p1","AnaCut_both_chf0p1",
			"AnaCut_jet1_chmult","AnaCut_jet2_chmult","AnaCut_pair_chmult","AnaCut_both_chmult",
			"AnaCut_both_emf0p2_chmult","AnaCut_both_emf0p2_chf0p1"};

  TString nomVar[nVar]={"NbOfVertices"};

  TString nomSamp[nSamp]={"QCD_HT_100_250", "QCD_HT_250_500", "QCD_HT_500_1000", "QCD_HT_1000_Inf"};
  TString nomplot="", nomdir="";
  TDirectory* myDir;

  for(    u_int iSel=0 ; iSel<nSel ; iSel++) {

    if(inclusive && iSel<=7) continue;

    for(  u_int iVar=0 ; iVar<nVar ; iVar++)  {
      for(u_int iSamp=0; iSamp<nSamp;iSamp++)  {
	
	nomdir  = "MultiSamplePlot_"+nomVar[iVar]+"_"+nomSel[iSel];
	nomplot = nomVar[iVar]+nomSel[iSel]+"_"+nomSamp[iSamp];

	myDir = fin->GetDirectory(nomdir);

	if(!myDir) {
	  if(verbose>1) cout << "ERROR : no such TDirectory named '" << nomdir << endl;
	  continue;
	}

	myDir->GetObject( nomplot , AllPlots[nomSel[iSel]][nomVar[iVar]][nomSamp[iSamp]] );

	if(!AllPlots[nomSel[iSel]][nomVar[iVar]][nomSamp[iSamp]]) {
	  if(verbose>1) cout << "ERROR : no such plot named '" << nomplot << endl;
	  continue;
	}

	if(iSamp==0)
	  AllPlots[nomSel[iSel]][nomVar[iVar]]["total"] = (TH1F*) AllPlots[nomSel[iSel]][nomVar[iVar]][nomSamp[iSamp]]->Clone();
	else
	  AllPlots[nomSel[iSel]][nomVar[iVar]]["total"]->Add( AllPlots[nomSel[iSel]][nomVar[iVar]][nomSamp[iSamp]] );
      }
    }
  }

  /////////////////////////
  // PLAY WITH THE PLOTS //
  /////////////////////////
  
  // plot invariant mass
  TLegend* leg_ = new TLegend(x1,y1,x2,y2);
  leg_->SetFillColor(0);
  leg_->SetTextFont(42);
  leg_->SetLineColor(1);
  leg_->SetLineWidth(1);
  leg_->SetLineStyle(0);
  leg_->SetBorderSize(1);	
  

  cout << "-- end loop over AllPlots" << endl
       << "-- closes files" << endl;
  
  fout->Close();
  fin->Close();
  
  delete fin;
  delete fout;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << "s to run the program" << endl;
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}




  

  //AllPlots[nomSel[iSel]]["RhoCorrection"]     = (MultiSamplePlot*) fin->Get("MultiSamplePlot_RhoCorrection"+nomSel[iSel]);  
//     //
//     AllPlots[nomSel[iSel]]["NbOfSelectedJets"]  =    

//     // Jet Variables
//     string nameJet="";
//     string nameJetH="";
//     string titleJet[m_nJet]={"Leading Jet", "Sub-Leading Jet"};

//     AllPlots[nomSel[iSel]]["em_E_frac_dijet"     ] = 
//     AllPlots[nomSel[iSel]]["charged_E_frac_dijet"] = 
//     AllPlots[nomSel[iSel]]["chargedMult_dijet"   ] = 
//     AllPlots[nomSel[iSel]]["mass_dijet"          ] = 
    
//     for(int iJ=0 ; iJ<m_nJet ; iJ++) {

//       nameJet  = Form("_Jet%d",iJ+1);
//       nameJetH = nameJet + "_" + nomSel[iSel];

//       cout << "--- nameJet=" << nameJet << endl;
    
//       AllPlots[nomSel[iSel]]["Eta"+nameJet] = 
//       AllPlots[nomSel[iSel]]["Phi"+nameJet] = 
//       AllPlots[nomSel[iSel]]["Pt" +nameJet] = 
//       AllPlots[nomSel[iSel]]["E" +nameJet]  = 

//       AllPlots[nomSel[iSel]]["chargedHad_E_frac"+nameJet] = 
//       AllPlots[nomSel[iSel]]["chargedEm_E_frac" +nameJet] = 
//       AllPlots[nomSel[iSel]]["chargedMu_E_frac" +nameJet] = 
//       AllPlots[nomSel[iSel]]["neutralEm_E_frac" +nameJet] = 
//       AllPlots[nomSel[iSel]]["neutralHad_E_frac"+nameJet] = 

//       AllPlots[nomSel[iSel]]["em_E_frac"     +nameJet]= 
//       AllPlots[nomSel[iSel]]["charged_E_frac"+nameJet]= 
//       AllPlots[nomSel[iSel]]["neutral_E_frac"+nameJet]= 

//       AllPlots[nomSel[iSel]]["chargedMult"      +nameJet] = 
//       AllPlots[nomSel[iSel]]["neutralMult"      +nameJet] = 
//       AllPlots[nomSel[iSel]]["muonMult"         +nameJet] = 

//       AllPlots[nomSel[iSel]]["chargedMultTot"  +nameJet] = 
//     }
//  }
