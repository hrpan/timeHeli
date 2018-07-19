#include "range.h"

bool nTagged(vector<float> *nTag_e, vector<float> *nTag_dt, size_t start, size_t end);
int findRange(double *range, int n_range, double npe);

const double nTag_dt_cut = 800;

void preprocess(){
	TChain *chain = new TChain("Heli");

	chain->Add("../p17b/data/heli.root");
	
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

	
	double _dtlSH[n_range];
	double _dtlSH_tag[n_range];
	double _dtlSH_atag[n_range];
	
	int _dtlSH_tag_n[n_range];


	f_out->cd();

	TTree *tr_out = new TTree("tr","tr");
	tr_out->Branch("site",&site);
	tr_out->Branch("detector",&detector);
	tr_out->Branch("ep",&ep);
	tr_out->Branch("ed",&ed);
	tr_out->Branch("dt",&dt);
	tr_out->Branch("dist",&dist);
	for(size_t i=0;i<n_range;++i){
		char buf[255];
		sprintf(buf,"dtlSH_%d",i);
		tr_out->Branch(buf,&_dtlSH[i]);
		sprintf(buf,"dtlSH_%d_tag",i);
		tr_out->Branch(buf,&_dtlSH_tag[i]);
		sprintf(buf,"dtlSH_%d_tag_n",i);
		tr_out->Branch(buf,&_dtlSH_tag_n[i]);
		sprintf(buf,"dtlSH_%d_atag",i);
		tr_out->Branch(buf,&_dtlSH_atag[i]);
	}
	float distCut = 1000;

	for(size_t i=0;i<entries;++i){
		if(i%(entries/10000)==0){
			printf("\r%d/%d(%.3f%%)",i,entries,float(i*100)/entries);
			fflush(stdout);
		}

		chain->GetEntry(i);
		
		if(dist>distCut) continue;

		size_t muons = dtlSH->size();

		for(int j=0;j<n_range;++j){
			_dtlSH[j] = -1;
			_dtlSH_tag[j] = -1;
			_dtlSH_tag_n[j] = -1;
			_dtlSH_atag[j] = -1;
		}	

		size_t n_idx_0 = 0;

		for(size_t mu = 0;mu < muons; ++mu){
		
			if( (*dtlSH)[mu] < 0 ) break;

			bool tagged = nTagged(nTag_e, nTag_dt, n_idx_0, n_idx_0 + (*nTag_n)[mu]);

			int r_idx = findRange(range,n_range,(*nPESum)[mu]);

			_dtlSH[r_idx] = (*dtlSH)[mu];
			
			if(tagged){
				_dtlSH_tag[r_idx] = (*dtlSH)[mu];
				_dtlSH_tag_n[r_idx] = (*nTag_n)[mu];
			}else{
				_dtlSH_atag[r_idx] = (*dtlSH)[mu];
			}
			
			n_idx_0 += (*nTag_n)[mu];

		}

		tr_out->Fill();
	
	}
	tr_out->Write("tr");	
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
