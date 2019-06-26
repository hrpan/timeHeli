#include "range.h"

bool nTagged(vector<float> *nTag_e, vector<float> *nTag_dt, size_t start, size_t end);
int findRange(double *range, int n_range, double npe);

const double nTag_dt_cut = 800;

void checkPoisson(){
	TChain *chain = new TChain("Heli");

	chain->Add("../p17b/data/heli.root");

	TFile *f_out = new TFile("./data/poisson.root","RECREATE");

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

	double _last[n_range];
	double _last_tag[n_range];
	double _last_atag[n_range];
	double _t;
	double _nPESum;
	int _tagged;
	TTree *tr = new TTree("tr","tr");
	tr->Branch("site",&site);
	tr->Branch("detector",&detector);
	tr->Branch("nPESum",&_nPESum);
	tr->Branch("tagged",&_tagged);
	tr->Branch("t",&_t);
	char buf[255];
	for(size_t i=0;i<n_range;++i){
		sprintf(buf,"last_%d",i);
		tr->Branch(buf, &_last[i]); 
		sprintf(buf,"last_tag_%d",i);
		tr->Branch(buf, &_last_tag[i]); 
		sprintf(buf,"last_atag_%d",i);
		tr->Branch(buf, &_last_atag[i]); 
	}


	float distCut = 1000;

	for(size_t i=0;i<entries;++i){
		if(i%(entries/10000)==0){
			printf("\r%d/%d(%.3f%%)",i,entries,float(i*100)/entries);
			fflush(stdout);
		}
		if(i > entries / 100) break;
		chain->GetEntry(i);
		
		if(dist>distCut) continue;

		int s = site==4? 2:site-1;
		for(int j=0;j<n_range;++j){
			_last[j] = -1;
			_last_tag[j] = -1;
			_last_atag[j] = -1;
		}	

		size_t n_idx_0 = 0;
		size_t muons = dtlSH->size();
		
		for(size_t mu = 0;mu < muons; ++mu){
			bool tagged = nTagged(nTag_e, nTag_dt, n_idx_0, n_idx_0 + (*nTag_n)[mu]);
			_tagged = (int)tagged;

			int r_idx = findRange(range,n_range,(*nPESum)[mu]);
			_t = (*dtlSH)[mu];
			_nPESum = (*nPESum)[mu];
			tr->Fill();

			_last[r_idx] = (*dtlSH)[mu];
			
			if(tagged)
				_last_tag[r_idx] = (*dtlSH)[mu];
			else
				_last_atag[r_idx] = (*dtlSH)[mu];
			
			
			n_idx_0 += (*nTag_n)[mu];


		}

	}

	tr->Write("tr_out");

	f_out->Close();

}

bool nTagged(vector<float> *nTag_e, vector<float> *nTag_dt, size_t start, size_t end){
	
	bool tagged = false;

	for(size_t i=start; i<end && !tagged && (*nTag_dt)[i] < nTag_dt_cut ;++i){

		bool nTag_e_cut = (*nTag_e)[i] > 1.8;
		bool nTag_dt_cut = (*nTag_dt)[i] > 20.0;

		tagged = nTag_e_cut && nTag_dt_cut;
		//tagged = (*nTag_dt)[i] > 10;
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
