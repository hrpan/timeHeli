#include "range.h"

bool nTagged(vector<float> *nTag_e, vector<float> *nTag_dt, size_t start, size_t end);
int findRange(double *range, int n_range, double npe);

const double nTag_e_min = 1.8;
const double nTag_dt_min = 20;
const double nTag_dt_max = 400;
const double mmv_dt_cut = 500;
const double iso_dt_cut = 0.5;
const float distCut = 1000;

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

	
	double _dtlSH[3][n_range];
	double _dtlSH_sdt[3][n_range];
	double _dtlSH_dt[3][n_range];	
	double _dtlSH_nPESum[3][n_range];	
	int _dtlSH_tag_n[n_range];
	int _nSH[3][n_range];

	f_out->cd();

	TTree *tr_out = new TTree("tr","tr");
	tr_out->Branch("site",&site);
	tr_out->Branch("detector",&detector);
	tr_out->Branch("ep",&ep);
	tr_out->Branch("ed",&ed);
	tr_out->Branch("dt",&dt);
	tr_out->Branch("dist",&dist);

	char *suffix[3] = {
		"",
		"_tag",
		"_atag"
	};

	for(size_t i=0;i<n_range;++i){
		char buf[255];
		for(size_t j=0;j<3;++j){
			sprintf(buf,"dtlSH_%d%s",i,suffix[j]);
			tr_out->Branch(buf,&_dtlSH[j][i]);
			sprintf(buf,"dtlSH_sdt_%d%s",i,suffix[j]);
			tr_out->Branch(buf,&_dtlSH_sdt[j][i]);
			sprintf(buf,"dtlSH_dt_%d%s",i,suffix[j]);
			tr_out->Branch(buf,&_dtlSH_dt[j][i]);
			sprintf(buf,"dtlSH_nPESum_%d%s",i,suffix[j]);
			tr_out->Branch(buf,&_dtlSH_nPESum[j][i]);
			sprintf(buf,"nSH_%d%s",i,suffix[j]);
			tr_out->Branch(buf,&_nSH[j][i]);
		}
		sprintf(buf,"dtlSH_%d_tag_n",i);
		tr_out->Branch(buf,&_dtlSH_tag_n[i]);
	}

	for(size_t i=0;i<entries;++i){
		if(i%(entries/10000)==0){
			printf("\r%d/%d(%.3f%%)",i,entries,float(i*100)/entries);
			fflush(stdout);
		}

		chain->GetEntry(i);
		
		if(dist>distCut) continue;

		size_t muons = dtlSH->size();

		double _dtlSH_all[n_range];
		for(size_t j=0;j<n_range;++j){
			for(size_t k=0;k<3;++k){
				_dtlSH[k][j] = -1;
				_dtlSH_sdt[k][j] = -1;
				_dtlSH_dt[k][j] = -1;
				_dtlSH_nPESum[k][j] = -1;
				_nSH[k][j] = 0;
			}
			_dtlSH_all[j] = -1;
			_dtlSH_tag_n[j] = -1;
		}	

		size_t n_idx_0 = 0;
		
		for(size_t mu = 0; mu < muons; ++mu){
		
			if( (*dtlSH)[mu] < 0 ) break;

			bool tagged = nTagged(nTag_e, nTag_dt, n_idx_0, n_idx_0 + (*nTag_n)[mu]);
			n_idx_0 += (*nTag_n)[mu];

			int r_idx = findRange(range,n_range,(*nPESum)[mu]);

			_dtlSH_all[r_idx] = (*dtlSH)[mu];
			bool pass = true;
			/*
			for(int r_tmp = r_idx+1; r_tmp < n_range && pass; ++r_tmp)
				if(r_tmp == r_idx) 
					continue;
				else if(_dtlSH_all[r_tmp] > 0)
					pass = (_dtlSH_all[r_tmp] - (*dtlSH)[mu]) > mmv_dt_cut;
				else if(_dtlSH_all[r_tmp] < 0)
					pass = (*dtlSH)[mu] < 5000 - mmv_dt_cut;
			*/
			if( mu > 0 )
				pass = (*dtlSH)[mu-1] - (*dtlSH)[mu] > iso_dt_cut;
			else 
				pass = (*dtlSH)[mu] < 5000 - iso_dt_cut;

			if( mu + 1 < muons )
				pass = pass && ( (*dtlSH)[mu] - (*dtlSH)[mu+1] > iso_dt_cut );
			else
				pass = pass && ( (*dtlSH)[mu] > iso_dt_cut );

			if(!pass) continue;

			if(_dtlSH[0][r_idx] > 0)
				_dtlSH_sdt[0][r_idx] = _dtlSH[0][r_idx] - (*dtlSH)[mu];
			_dtlSH[0][r_idx] = (*dtlSH)[mu];
			++_nSH[0][r_idx];
			if(mu > 0){
				_dtlSH_dt[0][r_idx] = (*dtlSH)[mu-1] - (*dtlSH)[mu];
				_dtlSH_nPESum[0][r_idx] = (*nPESum)[mu-1]; 
			}
			if(tagged){
				if(_dtlSH[1][r_idx] > 0)
					_dtlSH_sdt[1][r_idx] = _dtlSH[1][r_idx] - (*dtlSH)[mu];
				_dtlSH[1][r_idx] = (*dtlSH)[mu];
				++_nSH[1][r_idx];
				_dtlSH_tag_n[r_idx] = (*nTag_n)[mu];
				if(mu > 0){
					_dtlSH_dt[1][r_idx] = (*dtlSH)[mu-1] - (*dtlSH)[mu];
					_dtlSH_nPESum[1][r_idx] = (*nPESum)[mu-1]; 
				}
			}else{
				if(_dtlSH[2][r_idx] > 0)
					_dtlSH_sdt[2][r_idx] = _dtlSH[2][r_idx] - (*dtlSH)[mu];
				_dtlSH[2][r_idx] = (*dtlSH)[mu];
				++_nSH[2][r_idx];
				if(mu > 0){
					_dtlSH_dt[2][r_idx] = (*dtlSH)[mu-1] - (*dtlSH)[mu];
					_dtlSH_nPESum[2][r_idx] = (*nPESum)[mu-1]; 
				}
			}

		}

		tr_out->Fill();
	
	}
	tr_out->Write("tr");

	TString meta = "";	
	for(size_t i=0; i<n_range;++i)
		meta += TString::Format("%.2f ", range[i]); 
	meta += TString::Format("nTag_e_min: %.2f ", nTag_e_min); 
	meta += TString::Format("nTag_dt_min: %.2f ", nTag_dt_min); 
	meta += TString::Format("nTag_dt_max: %.2f ", nTag_dt_max); 
	meta += TString::Format("iso_dt_cut: %.2f ", iso_dt_cut); 
	TNamed(meta, meta).Write("Metadata");
	
	f_out->Close();
}

bool nTagged(vector<float> *nTag_e, vector<float> *nTag_dt, size_t start, size_t end){
	
	bool tagged = false;

	for(size_t i=start; i<end && !tagged && (*nTag_dt)[i] < nTag_dt_max ;++i){

		bool _e_cut = (*nTag_e)[i] > nTag_e_min;
		bool _dt_cut = (*nTag_dt)[i] > nTag_dt_min;

		tagged = _e_cut && _dt_cut;
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
