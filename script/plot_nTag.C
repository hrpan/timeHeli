const int n_range = 7;

void plot_nTag(){
	TFile *f = new TFile("./data/proced.root");

	TTree *tr = f->Get("tr");

	for(int i=0;i<3;++i){

		char siteCut[255];
		int s = (i==2)? 4:i+1;
		sprintf(siteCut,"site==%d",s);

		for(int r=0;r<n_range;++r){

			char buf[255];
			
			sprintf(buf,"dtlSH_%d_tag_n>>h(80,0,80)",r);

			tr->Draw(buf,siteCut);
			
			sprintf(buf,"plots/EH%d_nTag_%d.png",i+1,r);
			gPad->SetLogy();
			c1->SaveAs(buf);	

		}

	}

}
