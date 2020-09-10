
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
#include "TMath.h"

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
RVec<size_t> sortbtag(rvec_f jet_pt, rvec_f jet_eta, rvec_f btag){
    auto sort_btag = Reverse(Argsort(btag));
    RVec<size_t> sortjet_idx;
    for( size_t i = 0; i < jet_pt.size(); i++ ){
        if( jet_pt[sort_btag[i]] > 30 && abs(jet_eta[sort_btag[i]]) < 2.4 ){
            sortjet_idx.push_back(sort_btag[i]);
            }
    }
    return sortjet_idx;
}

//********** signal
int categorize(rvec_f jet_pt, rvec_f jet_eta, rvec_f jet_phi, rvec_f jet_e, rvec_f btag,               //jet
	           float addbjet1_pt, float addbjet1_eta, float addbjet1_phi, float addbjet1_e,         //addbjet1
	           float addbjet2_pt, float addbjet2_eta, float addbjet2_phi, float addbjet2_e){        //addbjet2  

    //********** sorting jets by btag
    auto sort_btag = Reverse(Argsort(btag));
    RVec<size_t> sortjet_idx;
    for( size_t i = 0; i < jet_pt.size(); i++ ){
        if( jet_pt[sort_btag[i]] > 30 && abs(jet_eta[sort_btag[i]]) < 2.4 ){
            sortjet_idx.push_back(sort_btag[i]);
            }
    }

    //********** save additional b-jets in Lorentzvector
    TLorentzVector addbjet1, addbjet2;
    addbjet1.SetPtEtaPhiE(addbjet1_pt, addbjet1_eta, addbjet1_phi, addbjet1_e);
    addbjet2.SetPtEtaPhiE(addbjet2_pt, addbjet2_eta, addbjet2_phi, addbjet2_e);

    //********** save jets in Lorentzvector and make a vector of jets
    vector<TLorentzVector> Jets;
    TLorentzVector jet;
    for (int i = 0; i < sortjet_idx.size(); i++){
        jet.SetPtEtaPhiE(jet_pt[sortjet_idx[i]], jet_eta[sortjet_idx[i]], jet_phi[sortjet_idx[i]], jet_e[sortjet_idx[i]]);
        Jets.push_back(jet);
    }
    
    //********** match jets with additional b-jets
    int size = Jets.size();
    if (Jets.size() > 6){
        size = 6;
    }

    int b1cnt = 1;
    int b2cnt = 1;
    int cnt = 2;
    float dR1;
    float dR2;

    vector<TLorentzVector> out; //vector of mathced jets

    for (int i = 0; i < size && cnt > 0; i++){
        dR1 = Jets[i].DeltaR(addbjet1);
        dR2 = Jets[i].DeltaR(addbjet2);
        if (dR1 < 0.4 && dR2 < 0.4){
            cnt--;
            out.push_back(Jets[i]);
        } else if (b1cnt > 0 && dR1 < 0.4){
            b1cnt--;
            cnt--;
            out.push_back(Jets[i]);
        } else if (b2cnt > 0 && dR2 < 0.4){
            b2cnt--;
            cnt--;
            out.push_back(Jets[i]);
        }
    }

    //********** find the indice of mathced jets and save as a number 
    int signal_num = 0;
    int cnt2 = out.size();
    int num = 0;

    for (int i = 0; i < Jets.size() && cnt2 > 0; i++){
        for (int j = 0; j < out.size(); j++){
            if (Jets[i].Pt() == out[j].Pt() && Jets[i].Eta() == out[j].Eta() && Jets[i].Phi() == out[j].Phi()){
                num = TMath::Power(2, i);
                signal_num = signal_num + num;
                cnt2--;
           }
        }
    }

    //*********** categorize
    int cat;
    if (signal_num == 3) cat = 0;
    else if (signal_num == 5) cat = 1;
    else if (signal_num == 9) cat = 2;
    else if (signal_num == 6) cat = 3;
    else if (signal_num == 10) cat = 4;
    else if (signal_num == 12) cat = 5;
    else if (signal_num == 17) cat = 6;
    else if (signal_num == 18) cat = 7;
    else if (signal_num == 20) cat = 8;
    else if (signal_num == 24) cat = 9;
    else if (signal_num == 33) cat = 10;
    else if (signal_num == 34) cat = 11;
    else if (signal_num == 36) cat = 12;
    else if (signal_num == 40) cat = 13;
    else if (signal_num == 48) cat = 14;
    else cat = 15;

    return cat;
}

//********** delta R
float deltaR (float jet1_pt, float jet1_eta, float jet1_phi, float jet1_e,          //jet1
              float jet2_pt, float jet2_eta, float jet2_phi, float jet2_e){         //jet2
    
    TLorentzVector jet1, jet2;

    jet1.SetPtEtaPhiE(jet1_pt, jet1_eta, jet1_phi, jet1_e);
    jet2.SetPtEtaPhiE(jet2_pt, jet2_eta, jet2_phi, jet2_e);

    float delR = jet1.DeltaR(jet2);

    return delR;
}

//********** delta Phi
float deltaPhi (float jet1_pt, float jet1_eta, float jet1_phi, float jet1_e,          //jet1
              float jet2_pt, float jet2_eta, float jet2_phi, float jet2_e){         //jet2
    
    TLorentzVector jet1, jet2;

    jet1.SetPtEtaPhiE(jet1_pt, jet1_eta, jet1_phi, jet1_e);
    jet2.SetPtEtaPhiE(jet2_pt, jet2_eta, jet2_phi, jet2_e);

    float delPhi = jet1.DeltaPhi(jet2);

    return delPhi;
}

//********** delta Eta
float deltaEta (float jet1_pt, float jet1_eta, float jet1_phi, float jet1_e,          //jet1
              float jet2_pt, float jet2_eta, float jet2_phi, float jet2_e){         //jet2
    
    TLorentzVector jet1, jet2;

    jet1.SetPtEtaPhiE(jet1_pt, jet1_eta, jet1_phi, jet1_e);
    jet2.SetPtEtaPhiE(jet2_pt, jet2_eta, jet2_phi, jet2_e);

    float delE = jet1.Eta() - jet2.Eta();

    return delE;
}

//********** invariant mass
float invmass (float jet1_pt, float jet1_eta, float jet1_phi, float jet1_e,
               float jet2_pt, float jet2_eta, float jet2_phi, float jet2_e){

    TLorentzVector jet1, jet2;

    jet1.SetPtEtaPhiE(jet1_pt, jet1_eta, jet1_phi, jet1_e);
    jet2.SetPtEtaPhiE(jet2_pt, jet2_eta, jet2_phi, jet2_e);

    float invariantmass = (jet1+jet2).M();

    return invariantmass;
}



void ana_ttbb(){
    ROOT::RDataFrame df("ttbbLepJets/tree", "/cms/ldap_home/sarakm0704/public/ttbb/V10_3/sync/TTLJ_PowhegPythia_ttbb.root");

    //********** Lepton cut
    auto df_goodlepton = df.Filter("lepton_pt > 30 && abs(lepton_eta) < 2.4")
                           .Define("goodjets","jet_pt > 30 && abs(jet_eta) < 2.4")
                           .Define("ngoodjets", "Sum(goodjets)")
                           .Define("bjets_m", "jet_deepJet > 0.277 && jet_pt > 30 && abs(jet_eta) < 2.4")
                           .Define("bjets_t", "jet_deepJet > 0.7264 && jet_pt > 30 && abs(jet_eta) < 2.4")
                           .Define("nbjets_t", "Sum(bjets_t)")
                           .Define("nbjets_m", "Sum(bjets_m)");
    //********** nJets >= 4
    auto df_goodjet = df_goodlepton.Filter("ngoodjets >=4 ", "Events with at least 4 goodjets");
    auto df_bjet = df_goodjet.Filter("nbjets_m >= 2", "Events with at least 2 medium b-jets")
                           //*** sort btag
                           .Define("sortjets", sortbtag, {"jet_pt", "jet_eta", "jet_deepJet"})
                           //*** categorize
                           .Define("category", categorize, {"jet_pt", "jet_eta", "jet_phi", "jet_e", "jet_deepJet",
                                                        "addbjet1_pt", "addbjet1_eta", "addbjet1_phi", "addbjet1_e",
                                                        "addbjet2_pt", "addbjet2_eta", "addbjet2_phi", "addbjet2_e"})
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
                           .Define("jet4_CvsL", "jet_deepJetCvsL[sortjets[3]]")
                           //*** deltaR between jets
                           .Define("deltaR_12", deltaR, {"jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e"})
                           .Define("deltaR_13", deltaR, {"jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet3_pt", "jet3_eta", "jet3_phi", "jet3_e"})
                           .Define("deltaR_14", deltaR, {"jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet4_pt", "jet4_eta", "jet4_phi", "jet4_e"})
                           .Define("deltaR_23", deltaR, {"jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e"})
                           .Define("deltaR_24", deltaR, {"jet4_pt", "jet4_eta", "jet4_phi", "jet4_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e"})
                           .Define("deltaR_34", deltaR, {"jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "jet4_pt", "jet4_eta", "jet4_phi", "jet4_e"})
                           //*** deltaR between jet and lepton
                           .Define("deltaR_lepj1", deltaR, {"jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "lepton_pt", "lepton_eta", "lepton_phi", "lepton_e"})
                           .Define("deltaR_lepj2", deltaR, {"jet2_pt", "jet2_eta", "jet2_phi", "jet2_e", "lepton_pt", "lepton_eta", "lepton_phi", "lepton_e"})
                           .Define("deltaR_lepj3", deltaR, {"jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "lepton_pt", "lepton_eta", "lepton_phi", "lepton_e"})
                           .Define("deltaR_lepj4", deltaR, {"jet4_pt", "jet4_eta", "jet4_phi", "jet4_e", "lepton_pt", "lepton_eta", "lepton_phi", "lepton_e"})
                           //*** invariant mass of jets
                           .Define("invmass_12", invmass, {"jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e"})
                           .Define("invmass_13", invmass, {"jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet3_pt", "jet3_eta", "jet3_phi", "jet3_e"})
                           .Define("invmass_14", invmass, {"jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "jet4_pt", "jet4_eta", "jet4_phi", "jet4_e"})
                           .Define("invmass_23", invmass, {"jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e"})
                           .Define("invmass_24", invmass, {"jet4_pt", "jet4_eta", "jet4_phi", "jet4_e", "jet2_pt", "jet2_eta", "jet2_phi", "jet2_e"})
                           .Define("invmass_34", invmass, {"jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "jet4_pt", "jet4_eta", "jet4_phi", "jet4_e"})
                           //*** invariant mass of jet&lepton
                           .Define("invmass_lepj1", invmass, {"jet1_pt", "jet1_eta", "jet1_phi", "jet1_e", "lepton_pt", "lepton_eta", "lepton_phi", "lepton_e"})
                           .Define("invmass_lepj2", invmass, {"jet2_pt", "jet2_eta", "jet2_phi", "jet2_e", "lepton_pt", "lepton_eta", "lepton_phi", "lepton_e"})
                           .Define("invmass_lepj3", invmass, {"jet3_pt", "jet3_eta", "jet3_phi", "jet3_e", "lepton_pt", "lepton_eta", "lepton_phi", "lepton_e"})
                           .Define("invmass_lepj4", invmass, {"jet4_pt", "jet4_eta", "jet4_phi", "jet4_e", "lepton_pt", "lepton_eta", "lepton_phi", "lepton_e"});
    cout << "good" << endl;
    df_bjet.Snapshot("dnn_input", "ttbb_15cat.root", {"category", "jet1_pt", "jet1_eta", "jet1_btag", "jet1_e", "jet2_pt", "jet2_eta", "jet2_btag", "jet2_e", "jet3_pt", "jet3_eta", "jet3_btag", "jet3_e", "jet4_pt", "jet4_eta", "jet4_btag", "jet4_e", "deltaR_12", "deltaR_13", "deltaR_14", "deltaR_23", "deltaR_24", "deltaR_34", "lepton_pt", "lepton_eta", "lepton_e", "deltaR_lepj1", "deltaR_lepj2", "deltaR_lepj3", "deltaR_lepj4", "invmass_12", "invmass_13", "invmass_14", "invmass_23", "invmass_24", "invmass_34", "invmass_lepj1", "invmass_lepj2", "invmass_lepj3", "invmass_lepj4", "jet1_CvsB", "jet1_CvsL", "jet2_CvsB", "jet2_CvsL", "jet3_CvsB", "jet3_CvsL", "jet4_CvsB", "jet4_CvsL"});

}



















