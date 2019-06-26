void printProgress(size_t current, size_t total){
	printf("\r%d/%d(%.3f%%)",current,total,float(current*100)/total);
	fflush(stdout);
}

void splitMuons(){
	TFile *f_mu = new TFile("data/muons.root");
	TTree *tr_mu = f_mu->Get("tr_mu");

	TFile *f_mu_out = new TFile("data/muons_split.root","RECREATE");
	TTree *tr = new TTree("tr","tr");

	int site, detector;
	double nPESum;
	int nTag_n, wpTag_n;
	vector<int> *wpTag_nHit;
	vector<float> *nTag_e, *nTag_dt, *wpTag_dt;
	
	tr_mu->SetBranchAddress("site",&site);
	tr_mu->SetBranchAddress("detector",&detector);
	tr_mu->SetBranchAddress("nPESum",&nPESum);
	tr_mu->SetBranchAddress("nTag_n",&nTag_n);
	tr_mu->SetBranchAddress("nTag_e",&nTag_e);
	tr_mu->SetBranchAddress("nTag_dt",&nTag_dt);
	tr_mu->SetBranchAddress("wpTag_n",&wpTag_n);
	tr_mu->SetBranchAddress("wpTag_nHit",&wpTag_nHit);
	tr_mu->SetBranchAddress("wpTag_dt",&wpTag_dt);


	tr->Branch("site",&site);
	tr->Branch("detector",&detector);
	tr->Branch("nPESum",&nPESum);
	tr->Branch("nTag_n",&nTag_n);
	tr->Branch("nTag_e",nTag_e);
	tr->Branch("nTag_dt",nTag_dt);
	tr->Branch("wpTag_n",&wpTag_n);
	tr->Branch("wpTag_nHit",wpTag_nHit);
	tr->Branch("wpTag_dt",wpTag_dt);


	size_t entries = tr_mu->GetEntries();

	size_t splitSize = 1000000;

	size_t filled[3] = {0, 0, 0};
		
	for(size_t i=0;i<entries;++i){

		tr_mu->GetEntry(i);
		int _s = site==4? 2:site-1;
		if(filled[_s] < splitSize){	
			tr->Fill();
			++filled[_s];
			if(filled[_s] % 1234 == 0){
				printf("\rFilled: %8d/%8d  %8d/%8d  %8d/%8d",filled[0],splitSize,filled[1],splitSize,filled[2],splitSize);
				fflush(stdout); 
			}
		}
		if(filled[2] == splitSize && filled[1] == splitSize && filled[0] == splitSize)
			break;
	}
	tr->Write("tr");
	f_mu_out->Close();
}
