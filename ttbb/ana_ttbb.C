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

//********** sorting jets by btag
RVec<int> sortbtag(rvec_f jet_pt, rvec_f jet_eta, rvec_f btag){
    auto sort_btag = Reverse(Argsort(btag));
    RVec<int> sortjet_idx;
    for( int i = 0; i < jet_pt.size(); i++ ){
        if( jet_pt[sort_btag[i]] > 20 && abs(jet_eta[sort_btag[i]]) < 2.4 ){
            sortjet_idx.push_back(sort_btag[i]);
            }
    }
    return sortjet_idx;
}
//********** set TLorentzVector
TLorentzVector jet_vector(float pt, float eta, float phi, float e){
    TLorentzVector jet;
    jet.SetPtEtaPhiE(pt, eta, phi, e);
    return jet;
}


//********** signal count
int count_matched(rvec_i Idx, rvec_f jet_pt, rvec_f jet_eta, rvec_f jet_phi, rvec_f jet_e, rvec_f jet_btag,               //jet
                  float addbjet1_pt, float addbjet1_eta, float addbjet1_phi, float addbjet1_e,         //addbjet1
                  float addbjet2_pt, float addbjet2_eta, float addbjet2_phi, float addbjet2_e){        //addbjet2  

    //*** LorentsVector
    TLorentzVector addbjet1, addbjet2;
    addbjet1.SetPtEtaPhiE(addbjet1_pt, addbjet1_eta, addbjet1_phi, addbjet1_e);
    addbjet2.SetPtEtaPhiE(addbjet2_pt, addbjet2_eta, addbjet2_phi, addbjet2_e);

    TLorentzVector jet;
    float dR1;
    float dR2;
    int cnt12 = 0;
    int cnt1 = 0;
    int cnt2 = 0;
    for (int i = 0; i < 4; i++){
        jet.SetPtEtaPhiE(jet_pt[Idx[i]], jet_eta[Idx[i]], jet_phi[Idx[i]], jet_e[Idx[i]]);

        dR1 = jet.DeltaR(addbjet1);
        dR2 = jet.DeltaR(addbjet2);

        if (dR1 < 0.4 && dR2 < 0.4){
            cnt12 = cnt12 + 100;

        }else if (dR1 < 0.4){
            cnt1 = cnt1 + 10;

        }else if (dR2 < 0.4){
            cnt2 = cnt2 + 1;
        }
    }
    int count = cnt12 + cnt1 + cnt2;

    return count;
}

//********** matched jet index
vector <int> matchedjet(rvec_i Idx, rvec_f jet_pt, rvec_f jet_eta, rvec_f jet_phi, rvec_f jet_e, TLorentzVector addbjet){
    
    TLorentzVector jet;
    vector <int> matchedjet_idx;
    float dR;

    for(int i = 0; i < 4; i++){
        jet.SetPtEtaPhiE(jet_pt[Idx[i]], jet_eta[Idx[i]], jet_phi[Idx[i]], jet_e[Idx[i]]);
        dR = jet.DeltaR(addbjet);

        if (dR < 0.4) matchedjet_idx.push_back(Idx[i]);
    }
    return matchedjet_idx;
}

//********** matched jet pt
vector <float> matchedjet_pt(vector <int> matchedjet_idx, rvec_f jet_pt, rvec_f jet_eta, rvec_f jet_phi, rvec_f jet_e){
    TLorentzVector jet;
    vector <float> Pt;
    float pt; 
    for (int i = 0; i < matchedjet_idx.size(); i++){
        jet.SetPtEtaPhiE(jet_pt[matchedjet_idx[i]], jet_eta[matchedjet_idx[i]], jet_phi[matchedjet_idx[i]], jet_e[matchedjet_idx[i]]);
        pt = jet.Pt();
        Pt.push_back(pt);
    }
    return Pt;
}

//*********** matched jet  deltaR
vector <float> matchedjet_dR(vector <int> matchedjet_idx, rvec_f jet_pt, rvec_f jet_eta, rvec_f jet_phi, rvec_f jet_e, TLorentzVector addbjet){
    TLorentzVector jet;
    vector <float> deltaR;
    float dR;
    for (int i = 0; i < matchedjet_idx.size(); i++){
        jet.SetPtEtaPhiE(jet_pt[matchedjet_idx[i]], jet_eta[matchedjet_idx[i]], jet_phi[matchedjet_idx[i]], jet_e[matchedjet_idx[i]]);
        dR = jet.DeltaR(addbjet);
        deltaR.push_back(dR);
    }
    return deltaR;    
}

int getsize(vector <int> vec){
    int size = vec.size();
//    cout << "*************************" << endl;
//    cout << "size: " << size << endl;
//    for (int i = 0; i < size; i++)
//        cout << i << " " << vec[i] << endl;
    return size;
}


void ana_ttbb(){
  ROOT::RDataFrame df("ttbbLepJets/tree", "/cms/ldap_home/sarakm0704/public/ttbb/V10_3/sync/TTLJ_PowhegPythia_ttbb.root");

    //*** Lepton cut
    auto df_goodlepton = df.Filter("lepton_pt > 30 && abs(lepton_eta) < 2.4")
                           .Define("goodjets","jet_pt > 20 && abs(jet_eta) < 2.4")
                           .Define("ngoodjets", "Sum(goodjets)")
                           .Define("bjets_m", "jet_deepJet > 0.277 && jet_pt > 20 && abs(jet_eta) < 2.4")
                           .Define("bjets_t", "jet_deepJet > 0.7264 && jet_pt > 20 && abs(jet_eta) < 2.4")
                           .Define("nbjets_t", "Sum(bjets_t)")
                           .Define("nbjets_m", "Sum(bjets_m)");
    //*** nJets >= 4
    auto df_goodjet = df_goodlepton.Filter("ngoodjets >=4 ", "Events with at least 4 goodjets");
    auto df_bjet = df_goodjet.Filter("nbjets_m >= 2", "Events with at least 2 medium b-jets")
                           //*** sort btag
                           .Define("sortjets", sortbtag, {"jet_pt", "jet_eta", "jet_deepJet"})
                           //*** additional b-jets
                           .Define("addbjet1", jet_vector, {"addbjet1_pt", "addbjet1_eta", "addbjet1_phi", "addbjet1_e"})
                           .Define("addbjet2", jet_vector, {"addbjet2_pt", "addbjet2_eta", "addbjet2_phi", "addbjet2_e"})
                           //*** count # of mathed jets
                           .Define("count", count_matched, {"sortjets", "jet_pt", "jet_eta", "jet_phi", "jet_e", "jet_deepJet",
                                                            "addbjet1_pt", "addbjet1_eta", "addbjet1_phi", "addbjet1_e",
                                                            "addbjet2_pt", "addbjet2_eta", "addbjet2_phi", "addbjet2_e"})
                           //*** matched jet index
                           .Define("matched_add1_idx", matchedjet, {"sortjets", "jet_pt", "jet_eta", "jet_phi", "jet_e", "addbjet1"})
                           .Define("matched_add2_idx", matchedjet, {"sortjets", "jet_pt", "jet_eta", "jet_phi", "jet_e", "addbjet2"})
                           .Define("nmatched1", getsize, {"matched_add1_idx"})
                           .Define("nmatched2", getsize, {"matched_add2_idx"})

                         
                           .Define("matchedjet_pt1", matchedjet_pt, {"matched_add1_idx", "jet_pt", "jet_eta", "jet_phi", "jet_e"})
                           .Define("matchedjet_pt2", matchedjet_pt, {"matched_add2_idx", "jet_pt", "jet_eta", "jet_phi", "jet_e"})
                           .Define("matchedjet_dR1", matchedjet_dR, {"matched_add1_idx", "jet_pt", "jet_eta", "jet_phi", "jet_e", "addbjet1"})
                           .Define("matchedjet_dR2", matchedjet_dR, {"matched_add2_idx", "jet_pt", "jet_eta", "jet_phi", "jet_e", "addbjet2"})


                           //*** jet1
                           .Define("jet1_pt", "jet_pt[sortjets[0]]")
                           .Define("jet1_eta", "jet_eta[sortjets[0]]")
                           .Define("jet1_phi", "jet_phi[sortjets[0]]")
                           .Define("jet1_e", "jet_e[sortjets[0]]")
                           .Define("jet1_btag", "jet_deepJet[sortjets[0]]")
                           .Define("jet1_CvsB", "jet_deepJetCvsB[sortjets[0]]")
                           .Define("jet1_CvsL", "jet_deepJetCvsL[sortjets[0]]")
                           //*** jet2
                           .Define("jet2_pt", "jet_pt[sortjets[1]]")
                           .Define("jet2_eta", "jet_eta[sortjets[1]]")
                           .Define("jet2_phi", "jet_phi[sortjets[1]]")
                           .Define("jet2_e", "jet_e[sortjets[1]]")
                           .Define("jet2_btag", "jet_deepJet[sortjets[1]]")
                           .Define("jet2_CvsB", "jet_deepJetCvsB[sortjets[1]]")
                           .Define("jet2_CvsL", "jet_deepJetCvsL[sortjets[1]]")
                           //*** jet3
                           .Define("jet3_pt", "jet_pt[sortjets[2]]")
                           .Define("jet3_eta", "jet_eta[sortjets[2]]")
                           .Define("jet3_phi", "jet_phi[sortjets[2]]")
                           .Define("jet3_e", "jet_e[sortjets[2]]")
                           .Define("jet3_btag", "jet_deepJet[sortjets[2]]")
                           .Define("jet3_CvsB", "jet_deepJetCvsB[sortjets[2]]")
                           .Define("jet3_CvsL", "jet_deepJetCvsL[sortjets[2]]")
                           //*** jet4
                           .Define("jet4_pt", "jet_pt[sortjets[3]]")
                           .Define("jet4_eta", "jet_eta[sortjets[3]]")
                           .Define("jet4_phi", "jet_phi[sortjets[3]]")
                           .Define("jet4_e", "jet_e[sortjets[3]]")
                           .Define("jet4_btag", "jet_deepJet[sortjets[3]]")
                           .Define("jet4_CvsB", "jet_deepJetCvsB[sortjets[3]]")
                           .Define("jet4_CvsL", "jet_deepJetCvsL[sortjets[3]]");

    auto df_matched1 = df_bjet.Filter("nmatched1 > 0", "Events with matched jet with add1");
    auto df_matched2 = df_bjet.Filter("nmatched2 > 0", "Events with matched jet with add2");

    cout << "hi" << endl;

    auto h_matched1 = df_matched1.Histo2D({"h_matched1", "h_matched1", 100, 0, 0.4, 180, 20, 200}, "matchedjet_dR1", "matchedjet_pt1");
    auto h_matched2 = df_matched2.Histo2D({"h_matched2", "h_matched2", 100, 0, 0.4, 180, 20, 200}, "matchedjet_dR2", "matchedjet_pt2"); 
  //snapshot
  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "UPDATE";
  
  //df_bjet.Snapshot("dnn_input", "ttbb_pt20_test5.root", {"nmatched1", "count", "matched_add1_idx", "matched_add2_idx"});
    df_bjet.Snapshot("dnn_input", "ttbb_pt20_test6.root", {"nmatched1", "nmatched2", "count", "matched_add1_idx", "matched_add2_idx", "matchedjet_pt1", "matchedjet_pt2", "matchedjet_dR1", "matchedjet_dR2"});
    df_matched1.Snapshot("add1", "ttbb_pt20_test6.root", {"count", "matched_add1_idx", "matched_add2_idx", "matchedjet_pt1", "matchedjet_pt2", "matchedjet_dR1", "matchedjet_dR2"}, opts);
    df_matched2.Snapshot("add2", "ttbb_pt20_test6.root", {"count", "matched_add1_idx", "matched_add2_idx", "matchedjet_pt1", "matchedjet_pt2", "matchedjet_dR1", "matchedjet_dR2"}, opts);


    TFile f("ttbb_pt20_test6.root", "UPDATE");

        plot(h_matched1, "h_matched1");
        plot(h_matched2, "h_matched2");

    f.Close();
}
