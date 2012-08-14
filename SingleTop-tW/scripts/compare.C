#include "setTDRStyle.C"

void compare(){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  setTDRStyle();
  //  gROOT->SetBatch(1);


  pt_leading_twdr->SetLineColor(kBlue);
  pt_leading_tt->SetLineColor(kRed);
  
  pt_leading_twdr->SetLineWidth(2);
  pt_leading_tt->SetLineWidth(2);
  
 /*
 // Rebin 
  pt_leading_twdr->Rebin(2);
  pt_leading_tt->Rebin(2);
  */
  
  TCanvas *c1 = new TCanvas();
  pt_leading_twdr->Draw("h");
  pt_leading_tt->Draw("h, sames");
  
  leg = new TLegend(0.7,0.7,0.94,0.94);
  leg ->SetFillStyle(1001);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(1);
  leg ->AddEntry(pt_leading_twdr, "tW signal", "l");
  leg ->AddEntry(pt_leading_tt, "tt bckgr.", "l");
  leg->Draw();
  
  c1->SaveAs("plots/compare_tw_tt_pt_leading.png");
  
  TCanvas *c2 = new TCanvas();
  // This is only for visualization
  pt_leading_twdr->SetNormFactor(1);
  pt_leading_tt->SetNormFactor(1);
  
  /*
  //Scale 
  pt_leading_twdr->Scale(1); -> Means that all the entries are multiplied by 1
  pt_leading_tt->Scale(1);
  */
  
  pt_leading_twdr->Draw("h");
  pt_leading_tt->Draw("h, sames");
  
  leg->Draw();
  
  c2->SaveAs("plots/compare_tw_tt_pt_leading.png");


}
