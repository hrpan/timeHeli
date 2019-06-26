/*
const size_t n_range = 3;

double range[3] = {3e3,1.5e5,3e5};
*/
/*
const size_t n_range = 11;
double range[11] = {
			3.0e3,3.0e4,6.0e4,
			9.0e4,1.2e5,1.5e5,
			1.8e5,2.1e5,2.4e5,
			2.7e5,3.0e5
			};
*/
const size_t n_range = 7;
double range[n_range] = {
			3.0e3,0.5e5,1.0e5,
			1.5e5,2.0e5,2.5e5,
			3.0e5
			};


int findRange(double *range, int n_range, double npe);

void getMuons(){
	TFile *file = new TFile("../p17b/data/heli.root");

	TTree *tr = file->Get("Heli");

	size_t entries = tr->GetEntries();

	int site, detector;

	float dist;	
	vector<double> *nPESum;

	vector<int>	*nTag_n;

	vector<float> *nTag_e;
	vector<float> *nTag_dt;

	vector<int> *wpTag_n;
	vector<int> *wpTag_nHit;
	vector<float> *wpTag_dt;

	tr->SetBranchAddress("site",&site);
	tr->SetBranchAddress("detector",&detector);
	tr->SetBranchAddress("dist",&dist);
	tr->SetBranchAddress("nPESum",&nPESum);
	tr->SetBranchAddress("nTag_n",&nTag_n);
	tr->SetBranchAddress("nTag_e",&nTag_e);
	tr->SetBranchAddress("nTag_dt",&nTag_dt);
	tr->SetBranchAddress("wpTag_n",&wpTag_n);
	tr->SetBranchAddress("wpTag_nHit",&wpTag_nHit);
	tr->SetBranchAddress("wpTag_dt",&wpTag_dt);	

	TFile *f_out = new TFile("./data/muons.root","RECREATE");
	TTree *tr_mu = new TTree("tr_mu","tr_mu");

	double _nPESum;
	int _nTag_n;
	vector<float> _nTag_e;
	vector<float> _nTag_dt;
	int _wpTag_n;
	vector<int> _wpTag_nHit;
	vector<float> _wpTag_dt;	

	tr_mu->Branch("site",&site);
	tr_mu->Branch("detector",&detector);
	tr_mu->Branch("nPESum",&_nPESum);
	tr_mu->Branch("nTag_n",&_nTag_n);
	tr_mu->Branch("nTag_e",&_nTag_e);
	tr_mu->Branch("nTag_dt",&_nTag_dt);
	tr_mu->Branch("wpTag_n",&_wpTag_n);
	tr_mu->Branch("wpTag_nHit",&_wpTag_nHit);
	tr_mu->Branch("wpTag_dt",&_wpTag_dt);

	for(size_t i=0;i<entries;++i){
		if(i%(entries/10000)==0){
			printf("\r%d/%d(%.2f%%)",i,entries,float(i*100)/entries);
			fflush(stdout);
		}
		
		tr->GetEntry(i);
		if(dist>1000) continue;	
		int det = detector==4? 2:detector-1;
	
		size_t muons = nTag_n->size();
		size_t n_idx_0 = 0;
		size_t wp_idx_0 = 0;
		for(size_t mu = 0; mu < muons; ++mu){

			_nTag_e.clear();
			_nTag_dt.clear();

			_wpTag_nHit.clear();
			_wpTag_dt.clear();

			_nPESum = (*nPESum)[mu];
			
			_nTag_n = (*nTag_n)[mu];
			for(size_t nd = 0; nd < _nTag_n; ++nd){
				_nTag_e.push_back((*nTag_e)[ n_idx_0 + nd ]);
				_nTag_dt.push_back((*nTag_dt)[ n_idx_0 + nd ]);
			}
			n_idx_0 += _nTag_n;

			_wpTag_n = (*wpTag_n)[mu];
			for(size_t nd = 0; nd < _wpTag_n; ++nd){
				_wpTag_nHit.push_back((*wpTag_nHit)[ wp_idx_0 + nd ]);
				_wpTag_dt.push_back((*wpTag_dt)[ wp_idx_0 + nd ]);
			}
			wp_idx_0 += _wpTag_n;

			tr_mu->Fill();

		}	

	}
	tr_mu->Write("tr_mu");
	f_out->Close();
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
