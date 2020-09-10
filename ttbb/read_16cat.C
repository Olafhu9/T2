void read_16cat(){


    TFile *f = new TFile("ttbb_16cat.root", "read");
    TTree *tree = (TTree*) f->Get("dnn_input");

    int category;
    tree->SetBranchAddress("category", &category);

    TH1F *h_cat = new TH1F("h_cat", "Category", 16, 0, 16);

    for (int irow = 0; irow < tree->GetEntries(); ++irow){
        tree->GetEntry(irow);
        h_cat->Fill(category);
    }

    TCanvas *c = new TCanvas;
    c->SetCanvasSize(550, 500);
    c->SetWindowSize(580, 580);
    c->SetLeftMargin(0.2);

    gStyle->SetOptStat("e");
    gStyle->SetStatX(0.55);
    gStyle->SetStatY(0.7);
    gStyle->SetStatW(0.3);
    gStyle->SetStatBorderSize(0);
    gStyle->SetStatFontSize(0.06);
    //gStyle->SetStatTextColor(kBlue);

    h_cat->SetXTitle("catetory");
    h_cat->SetYTitle("Entries");
    h_cat->SetLineWidth(2);
    h_cat->SetLineColor(kBlue-6);
    h_cat->SetMarkerSize(1.6);
    h_cat->Draw("hist text");

    c->Print(Form("16categories.pdf"));

}


