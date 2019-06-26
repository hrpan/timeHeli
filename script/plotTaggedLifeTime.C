#include "range.h"

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

void plotTaggedLifeTime(){
	//TFile *f = new TFile("./data/muons_split.root");
	TFile *f = new TFile("./data/muons.root");
	TTree *tr = f->Get("tr_mu");

	int site, detector;
	double nPESum;
	int nTag_n;
	vector<float> *nTag_e;
	vector<float> *nTag_dt;

	tr->SetBranchAddress("site",&site);		
	tr->SetBranchAddress("detector",&detector);		
	tr->SetBranchAddress("nPESum",&nPESum);
	tr->SetBranchAddress("nTag_n",&nTag_n);		
	tr->SetBranchAddress("nTag_e",&nTag_e);		
	tr->SetBranchAddress("nTag_dt",&nTag_dt);		


	TFile *f_out = new TFile("./data/tagged_dt.root","RECREATE");
	TH1 *h[3][n_range];
	char buf[255];

	for(int i=0;i<3;++i){
		for(int j=0;j<n_range;++j){
			sprintf(buf,"tagged_dt_%d_%d",i,j);
			h[i][j] = new TH1D(buf,buf,800,0,800);
		}
	}	


	size_t events = tr->GetEntries();
	cout << events << endl;

	for(size_t evt=0;evt<events;++evt){
		if(evt%12345==0){
			printf("\r%d/%d",evt,events);
			fflush(stdout);
		}
		tr->GetEntry(evt);
		int s = site==4? 2:site-1;
		size_t r_idx = findRange(range,n_range,nPESum);
		
		for(int n_idx=0;n_idx<nTag_n;++n_idx){
			
			h[site==4? 2:site-1][r_idx]->Fill((*nTag_dt)[n_idx]);

		}

	}

	for(int i=0;i<3;++i){
		for(int j=0;j<n_range;++j){
			h[i][j]->Write();
		}
	}
	
	f_out->Close();


}
