void preprocess(){
	TChain *chain = new TChain("Heli");

	chain->Add("../p17b/data_heli/run1000.root");
	
	TFile *f_out = new TFile("./data/proced.root","RECREATE");

	size_t entries = chain->GetEntries();

	int site, detector;

	float ep, ed;
	float dt, dist;
	
	vector<double> *dtlSH;
	vector<double> *nPESum;
	
	vector<int> *nTag_n;

	vector<float> *nTag_e;
	vector<float> *nTag_dt;

	chain->SetBranchAddress("site",&site);
	chain->SetBranchAddress("detector",&detector);
	chain->SetBranchAddress("ep",&ep);
	chain->SetBranchAddress("ed",&ed);
	chain->SetBranchAddress("dt",&dt);
	chain->SetBranchAddress("dist",&dist);
	chain->SetBranchAddress("dtlSH",&dtlSH);
	chain->SetBranchAddress("nPESum",&nPESum);
	chain->SetBranchAddress("nTag_n",&nTag_n);
	chain->SetBranchAddress("nTag_e",&nTag_e);
	chain->SetBranchAddress("nTag_dt",&nTag_dt);

	double range[3] = {3e3,1.5e5,3e5};
	
	double _dtlSH[3];
	double _dtlSH_tag[3];

	f_out->cd();

	TTree *tr_out = new TTree("tr","tr");
	tr_out->Branch("site",&site);
	tr_out->Branch("detector",&detector);
	tr_out->Branch("ep",&ep);
	tr_out->Branch("ed",&ed);
	tr_out->Branch("dt",&dt);
	tr_out->Branch("dist",&dist);
	tr_out->Branch("dtlSH_0",&_dtlSH[0]);	
	tr_out->Branch("dtlSH_1",&_dtlSH[1]);	
	tr_out->Branch("dtlSH_2",&_dtlSH[2]);	
	tr_out->Branch("dtlSH_0_tag",&_dtlSH_tag[0]);	
	tr_out->Branch("dtlSH_1_tag",&_dtlSH_tag[1]);	
	tr_out->Branch("dtlSH_2_tag",&_dtlSH_tag[2]);	

	float distCut = 1000;

	for(size_t i=0;i<entries;++i){
		cout << "\r" << i << "/" << entries << flush;
		chain->GetEntry(i);
		
		if(dist>distCut) continue;

		size_t muons = dtlSH->size();

		for(int j=0;j<3;++j){
			_dtlSH[j] = -1;
			_dtlSH_tag[j] = -1;
		}	

		size_t n_idx_0 = 0;

		for(size_t mu = 0;mu < muons; ++mu){
		
			if( (*dtlSH)[mu] < 0 ){cout << " " << muons - mu << endl; break;}

			bool tagged = nTagged(nTag_e, nTag_dt, n_idx_0, n_idx_0 + (*nTag_n)[mu]);

			int r_idx = findRange(range,3,(*nPESum)[mu]);

			_dtlSH[r_idx] = (*dtlSH)[mu];
					
			_dtlSH_tag[r_idx] = (*dtlSH)[mu];

		}

		tr_out->Fill();
	
	}
	tr_out->Write("tr");	
	f_out->Close();
}

bool nTagged(vector<float> *nTag_e, vector<float> *nTag_dt, size_t start, size_t end){
	
	bool tagged = false;
	for(size_t i=start; i<end && !tagged && (*nTag_dt)[i] < 200.0 ;++i){

		bool nTag_e_cut = (*nTag_e)[i] > 1.8;
		bool nTag_dt_cut = (*nTag_dt)[i] > 20.0;

		tagged = nTag_e_cut && nTag_dt_cut;
	}
	
	return tagged;
}

int findRange(double *range, int n_range, double npe){

	int idx = n_range - 1;

	for(int i=0;i<n_range-1;++i){
		if( npe < range[i+1] && npe > range[i]){
			idx = i;
			break;
		}
	}
	return idx;

}
