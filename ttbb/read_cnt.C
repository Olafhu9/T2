void read_cnt(){


    TFile *f = new TFile("ttbb_pt25.root", "read");
    //TFile *f = new TFile("ttbb_cnt.root", "read");
    TTree *tree = (TTree*) f->Get("dnn_input");

    int count;
    tree->SetBranchAddress("count", &count);

    TH1F *h_cnt = new TH1F("h_cnt", "count", 113, 0, 113);
    //TH1F *h_cnt = new TH1F("h_cnt", "count", 115, 0, 115);

    for (int irow = 0; irow < tree->GetEntries(); ++irow){
        tree->GetEntry(irow);
        h_cnt->Fill(count);
    }

    TCanvas *c = new TCanvas;
    //c->SetCanvasSize(550, 500);
    c->SetCanvasSize(2400, 400);
    c->SetWindowSize(600, 600);
    c->SetLeftMargin(0.2);
    //c->SetTopMargin(0.2);

    gStyle->SetOptStat("e");
    gStyle->SetStatX(0.53);
    gStyle->SetStatY(0.7);
    gStyle->SetStatW(0.3);
    gStyle->SetStatBorderSize(0);
    gStyle->SetStatFontSize(0.06);
    //gStyle->SetStatTextColor(kBlue);

    h_cnt->SetXTitle("catetory");
    h_cnt->SetYTitle("Entries");
    h_cnt->SetLineWidth(2);
    h_cnt->SetLineColor(kBlue-6);
    h_cnt->SetMarkerSize(1.6);
    h_cnt->Draw("hist text");

    c->Print(Form("count.pdf"));

}


