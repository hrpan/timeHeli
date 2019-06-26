void plotIBDHists(){
    TFile *f = new TFile("data/hists.root","READ");
    TList *l = f->GetListOfKeys();
    TObject *obj;
    TKey *key;
    TIter next(l);
    char buf[500];
    /*
    while(key = (TKey*) next()){
        TH1 *h = (TH1*) f->Get(key->GetName());
        sprintf(buf, "./plots/ibd/%s.png", h->GetName());
        if(h->GetDimension() == 1){
            h->Rebin(5);
            h->Draw();
            c1->SetLogy(1);
        }else if(h->GetDimension() == 2){
            h->Draw("colz");
            c1->SetLogy(0);
            c1->SetLogz(1);
        }        
        c1->SaveAs(buf);
    }
    */
    char *ana_name[4] = {
        "Unified (w/o dist. cut)",
        "nGd",
        "nH",
        "unified (w/ dist. cut)",
    };

    

    for(int var=0; var<4; ++var){
        for(int site=0;site<3;++site){
            TH1 *_h[4];
            TLegend *leg = new TLegend(0.65, 0.7, 0.9, 0.9);
            double _min = 1e9;
            for(int ana=0;ana<4;++ana){
                sprintf(buf, "h_var1d_%d_%d_%d", var, ana, site);
                _h[ana] = (TH1*) f->Get(buf);
                _h[ana]->SetLineColor(ana+1);
                _h[ana]->SetLineWidth(2);
                _h[ana]->Rebin(5);
                if(_h[ana]->GetMinimum(0) < _min)
                    _min = _h[ana]->GetMinimum(0);
                leg->AddEntry(_h[ana], ana_name[ana]);
            }
            for(int ana=0;ana<4;++ana){
                if(ana==0){
                    _h[ana]->SetMinimum(_min);
                    _h[ana]->Draw();
                }else{
                    _h[ana]->Draw("same");
                }
            }
            leg->Draw("same");
            c1->SetLogy();
            sprintf(buf, "./plots/ibd/h_var1d_%d_%d.png", var, site);
            c1->SaveAs(buf);
            sprintf(buf, "./plots/ibd/h_var1d_%d_%d.pdf", var, site);
            c1->SaveAs(buf);
        }
    }
}
