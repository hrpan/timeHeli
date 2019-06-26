void plotHists(){
	string filename;
	cout << "Enter File Name: ";
	cin >> filename;
	TFile *f = new TFile(filename.c_str());
	
	TList *l = f->GetListOfKeys();

	for(size_t i=0;i<l->GetEntries();++i){
		char *buf = l->At(i)->GetName();
		cout << buf << endl;
		TH1 *h = f->Get(buf);
		if(h->GetDimension()==2)
			h->Draw("colz");
		else
			h->Draw();
		
		char out[255];
		sprintf(out,"./plots/%s.png",buf);
		gPad->SetLogz();
		c1->SaveAs(out);

	}

}
