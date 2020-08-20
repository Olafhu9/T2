#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "Math/Vector4Dfwd.h"
#include "TStyle.h"
#include <vector>
#include <ROOT/RSnapshotOptions.hxx>

using namespace ROOT::VecOps;
using rvec_f = RVec<float>;
using rvec_b = RVec<bool>;
using rvec_i = RVec<int>;
using Impl_i = vector<int, ROOT::Detail::VecOps::RAdoptAllocator<int>>;

template <typename T>
void plot(T hist, TString name){
  TCanvas * c = new TCanvas("c",Form("c_%s", name.Data()));
  hist->Write();
 // c->Print(Form("%s.pdf",name.Data()));
}

//categorize events ( 0:lf, 1:b, 2:bb, 3:c, 4:, cc )
int categorize(int signal){
  int cat;
  if( signal == 12 ) cat = 0;
  else if( signal == 10) cat = 1;
  else if( signal == 9) cat = 2;
  else if( signal == 6) cat = 3;
  else if( signal == 5) cat = 4;
  else if( signal == 3) cat = 5;
  else cat = 6;
  return cat;
}

//sorting jets by btag
RVec<size_t> sortbtag(rvec_f jet_pt, rvec_f jet_eta, rvec_f btag){
    auto sort_btag = Reverse(Argsort(btag));
    RVec<size_t> sortjet_idx;
    for( size_t i = 0; i < jet_pt.size(); i++ ){
        if( jet_pt[sort_btag[i]] > 30 && abs(jet_eta[sort_btag[i]]) < 2.4 ){
            sortjet_idx.push_back(sort_btag[i]);
            }
    }
    for( int i = 0; i < sortjet_idx.size(); i++ ){
        cout << i << "eta: " << jet_eta[sortjet_idx[i]] << " pt: " << jet_pt[sortjet_idx[i]] << " btag: " << btag[sortjet_idx[i]] << endl;
    }
    return sortjet_idx;
}

//signal number
int signal_number(rvec_f jet_pt, rvec_f jet_eta, rvec_f jet_phi, rvec_f jet_e, rvec_f btag,               //jet
                     float addbjet1_pt, float addbjet1_eta, float addbjet1_phi, float addbjet1_e,         //addbjet1
                     float addbjet2_pt, float addbjet2_eta, float addbjet2_phi, float addbjet2_e){        //addbjet2  
  //sorting jets by btag
  auto sort_btag = Reverse(Argsort(btag));
  RVec<size_t> sortjet_idx;
  for( size_t i = 0; i < jet_pt.size(); i++ ){
      if( jet_pt[sort_btag[i]] > 30 && abs(jet_eta[sort_btag[i]]) < 2.4 ){
          sortjet_idx.push_back(sort_btag[i]);
          }
  }
//  for(int i = 0; i < jet_pt.size(); i++){
//      cout << i << " " << jet_pt[i] << " " << btag[i] << " " << jet_eta[i] << endl;
//  }
//
//  for( size_t i = 0; i < sortjet_idx.size(); i++){
//      cout << i << "btag: " << btag[sortjet_idx[i]] << endl;
//  }
  //LorentsVector
  TLorentzVector jet1, jet2, jet3, jet4, addbjet1, addbjet2;
  
  jet1.SetPtEtaPhiE(jet_pt[sortjet_idx[0]], jet_eta[sortjet_idx[0]], jet_phi[sortjet_idx[0]], jet_e[sortjet_idx[0]]);
  jet2.SetPtEtaPhiE(jet_pt[sortjet_idx[1]], jet_eta[sortjet_idx[1]], jet_phi[sortjet_idx[1]], jet_e[sortjet_idx[1]]);
  jet3.SetPtEtaPhiE(jet_pt[sortjet_idx[2]], jet_eta[sortjet_idx[2]], jet_phi[sortjet_idx[2]], jet_e[sortjet_idx[2]]);
  jet4.SetPtEtaPhiE(jet_pt[sortjet_idx[3]], jet_eta[sortjet_idx[3]], jet_phi[sortjet_idx[3]], jet_e[sortjet_idx[3]]);
  addbjet1.SetPtEtaPhiE(addbjet1_pt, addbjet1_eta, addbjet1_phi, addbjet1_e);
  addbjet2.SetPtEtaPhiE(addbjet2_pt, addbjet2_eta, addbjet2_phi, addbjet2_e);
  
  int signal = 0;

  if( jet1.DeltaR(addbjet1) < 0.4 || jet1.DeltaR(addbjet2) < 0.4){
      signal = signal + 8;
  }
  if( jet2.DeltaR(addbjet1) < 0.4 || jet2.DeltaR(addbjet2) < 0.4){
      signal = signal + 4;
  }
  if( jet3.DeltaR(addbjet1) < 0.4 || jet3.DeltaR(addbjet2) < 0.4){
      signal = signal + 2;
  }
  if( jet4.DeltaR(addbjet1) < 0.4 || jet4.DeltaR(addbjet2) < 0.4){
      signal = signal + 1;
  }

//  cout << signal << endl;

  return signal;
}


void ana_ttbb(){
  ROOT::RDataFrame df("ttbbLepJets/tree", "/cms/ldap_home/sarakm0704/public/ttbb/V10_3/sync/TTLJ_PowhegPythia_ttbb.root");

  //Lepton cut
  auto df_goodlepton = df.Filter("lepton_pt > 30 && abs(lepton_eta) < 2.4")
                         .Define("goodjets","jet_pt > 30 && abs(jet_eta) < 2.4")
                         .Define("ngoodjets", "Sum(goodjets)");
  //nJets >= 4
  auto df_goodjet = df_goodlepton.Filter("ngoodjets >=4 ", "Events with at least 4 goodjets")
                         //sort btag
                         .Define("sortjets", sortbtag, {"jet_pt", "jet_eta", "jet_deepJet"})
                         //signal
                         .Define("signal", signal_number, {"jet_pt", "jet_eta", "jet_phi", "jet_e", "jet_deepJet",
                                                           "addbjet1_pt", "addbjet1_eta", "addbjet1_phi", "addbjet1_e",
                                                           "addbjet2_pt", "addbjet2_eta", "addbjet2_phi", "addbjet2_e"})
                         //categorize 
                         .Define("category", categorize, {"signal"})
                         //jet1
                         .Define("jet1_pt", "jet_pt[sortjets[0]]")
                         .Define("jet1_eta", "jet_eta[sortjets[0]]")
                         .Define("jet1_phi", "jet_phi[sortjets[0]]")
                         .Define("jet1_e", "jet_e[sortjets[0]]")
                         //jet2
                         .Define("jet2_pt", "jet_pt[sortjets[1]]")
                         .Define("jet2_eta", "jet_eta[sortjets[1]]")
                         .Define("jet2_phi", "jet_phi[sortjets[1]]")
                         .Define("jet2_e", "jet_e[sortjets[1]]")
                         //jet3
                         .Define("jet3_pt", "jet_pt[sortjets[2]]")
                         .Define("jet3_eta", "jet_eta[sortjets[2]]")
                         .Define("jet3_phi", "jet_phi[sortjets[2]]")
                         .Define("jet3_e", "jet_e[sortjets[2]]")
                         //jet4
                         .Define("jet4_pt", "jet_pt[sortjets[3]]")
                         .Define("jet4_eta", "jet_eta[sortjets[3]]")
                         .Define("jet4_phi", "jet_phi[sortjets[3]]")
                         .Define("jet4_e", "jet_e[sortjets[3]]");
  //signal resions 
  auto df_cat0 = df_goodjet.Filter("category == 0");
  auto df_cat1 = df_goodjet.Filter("category == 1");
  auto df_cat2 = df_goodjet.Filter("category == 2");
  auto df_cat3 = df_goodjet.Filter("category == 3");
  auto df_cat4 = df_goodjet.Filter("category == 4");
  auto df_cat5 = df_goodjet.Filter("category == 5");
  auto df_cat6 = df_goodjet.Filter("category == 6");

  cout << "hi" << endl;

  //snapshot
  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "UPDATE";
  
  df_goodjet.Snapshot("dnn_input", "ttbb_merged.root", {"signal", "category", "jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e", "jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "jet4_pt", "jet4_eta", "jet4_phi", "jet4_e"});
  df_cat0.Snapshot("catetory_0", "ttbb_merged.root", {"category", "jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e", "jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "jet4_pt", "jet4_eta", "jet4_phi", "jet4_e"}, opts);
  df_cat1.Snapshot("catetory_1", "ttbb_merged.root", {"category", "jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e", "jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "jet4_pt", "jet4_eta", "jet4_phi", "jet4_e"}, opts);
  df_cat2.Snapshot("catetory_2", "ttbb_merged.root", {"category", "jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e", "jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "jet4_pt", "jet4_eta", "jet4_phi", "jet4_e"}, opts);
  df_cat3.Snapshot("catetory_3", "ttbb_merged.root", {"category", "jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e", "jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "jet4_pt", "jet4_eta", "jet4_phi", "jet4_e"}, opts);
  df_cat4.Snapshot("catetory_4", "ttbb_merged.root", {"category", "jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e", "jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "jet4_pt", "jet4_eta", "jet4_phi", "jet4_e"}, opts);
  df_cat5.Snapshot("catetory_5", "ttbb_merged.root", {"category", "jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e", "jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "jet4_pt", "jet4_eta", "jet4_phi", "jet4_e"}, opts);
  df_cat6.Snapshot("catetory_6", "ttbb_merged.root", {"category", "jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e", "jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "jet4_pt", "jet4_eta", "jet4_phi", "jet4_e"}, opts);
      
}
