#include "setTDRStyle.C"

void compare(){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0);
  setTDRStyle();
  //  gROOT->SetBatch(1);


  pt_leading_twdr->SetLineColor("kBlue");
  pt_leading_tt->SetLineColor("kRed");
  
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
  pt_leading_twdr->SetNormFactor(1);
  pt_leading_tt->SetNormFactor(1);
  
  pt_leading_twdr->Draw("h");
  pt_leading_tt->Draw("h, sames");
  
  leg->Draw();
  
  c2->SaveAs("plots/compare_tw_tt_pt_leading.png");


}
