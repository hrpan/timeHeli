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

bool nTagged(vector<float> *nTag_e, vector<float> *nTag_dt);

void plotMuons(){
	TFile *f = new TFile("./data/muons_split.root");
	
	TTree *tr_mu = f->Get("tr");

	int site, detector;
	double nPESum;
	int nTag_n;
	int wpTag_n;
	vector<float> *nTag_e;
	vector<float> *nTag_dt;
	vector<float> *wpTag_dt;
	vector<int> *wpTag_nHit;
	tr_mu->SetBranchAddress("site",&site);		
	tr_mu->SetBranchAddress("detector",&detector);		
	tr_mu->SetBranchAddress("nPESum",&nPESum);
	tr_mu->SetBranchAddress("nTag_n",&nTag_n);		
	tr_mu->SetBranchAddress("nTag_e",&nTag_e);		
	tr_mu->SetBranchAddress("nTag_dt",&nTag_dt);		
	tr_mu->SetBranchAddress("wpTag_n",&wpTag_n);		
	tr_mu->SetBranchAddress("wpTag_nHit",&wpTag_nHit);		
	tr_mu->SetBranchAddress("wpTag_dt",&wpTag_dt);		


	TFile *f_out = new TFile("./data/muon_plots.root","RECREATE");
	TH2D *h_npe_n[3][3];
	TH2D *h_npe_e[3][3];
	TH2D *h_npe_dt[3][3];
	TH2D *h_npe_wp_n[3][3];
	TH2D *h_npe_wp_nHit[3][3];
	TH2D *h_npe_wp_dt[3][3];
	TH2D *h_n_e[3][3];
	TH2D *h_n_dt[3][3];
	TH2D *h_e_dt[3][3];
	TH1D *h_dt[3][n_range][3];
	TH1D *h_e[3][n_range][3];
	TH1D *h_n[3][n_range][3];

	size_t bins_e = 100;
	size_t bins_dt = 100;
	size_t bins_n = 30;
	size_t bins_npe = 100;


	for(size_t i=0;i<3;++i){
		char buf[255];
		for(size_t j=0;j<3;++j){
			sprintf(buf,"h_npe_n_%d_%d",i,j);
			h_npe_n[i][j] = new TH2D(buf,buf,bins_npe,3e3,8e5,bins_n,0,bins_n);
			sprintf(buf,"h_npe_e_%d_%d",i,j);
			h_npe_e[i][j] = new TH2D(buf,buf,bins_npe,3e3,8e5,bins_e,1.5,12);
			sprintf(buf,"h_npe_dt_%d_%d",i,j);
			h_npe_dt[i][j] = new TH2D(buf,buf,bins_npe,3e3,8e5,bins_dt,0,800);

			sprintf(buf,"h_npe_wp_n_%d_%d",i,j);
			h_npe_wp_n[i][j] = new TH2D(buf,buf,bins_npe,3e3,8e5,20,0,20);
			sprintf(buf,"h_npe_wp_nHit_%d_%d",i,j);
			h_npe_wp_nHit[i][j] = new TH2D(buf,buf,bins_npe,3e3,8e5,88,12,100);
			sprintf(buf,"h_npe_wp_dt_%d_%d",i,j);
			h_npe_wp_dt[i][j] = new TH2D(buf,buf,bins_npe,3e3,8e5,100,-2,10);

			sprintf(buf,"h_n_e_%d_%d",i,j);
			h_n_e[i][j] = new TH2D(buf,buf,bins_n,0,bins_n,bins_e,1.5,12);
			sprintf(buf,"h_n_dt_%d_%d",i,j);
			h_n_dt[i][j] = new TH2D(buf,buf,bins_n,0,bins_n,bins_dt,0,800);
			sprintf(buf,"h_e_dt_%d_%d",i,j);
			h_e_dt[i][j] = new TH2D(buf,buf,bins_e,1.5,12,bins_dt,0,800);
			/*
			for(size_t k=0;k<n_range;++k){
				sprintf(buf,"h_dt_%d_%d_%d",i,k,j);
				h_dt[i][k][j] = new TH1D(buf,buf,bins_dt,0,800);
				sprintf(buf,"h_e_%d_%d_%d",i,k,j);
				h_e[i][k][j] = new TH1D(buf,buf,bins_e,1.5,12);
				sprintf(buf,"h_n_%d_%d_%d",i,k,j);
				h_n[i][k][j] = new TH1D(buf,buf,bins_n,0,80);
			}
			*/	
		}
	}


	size_t entries = tr_mu->GetEntries();

	for(size_t i=0;i<entries;++i){
		if(i%(entries/10000)==0){
			printf("\r%d/%d(%.2f%%)",i,entries,float(i*100)/entries);
			fflush(stdout);
		}

		tr_mu->GetEntry(i);

		int r_idx = findRange(range,n_range,nPESum);
		
		int idx = site==4? 2:site-1;	

		bool tagged = nTagged(nTag_e, nTag_dt);
		//bool tagged = wpTagged(wpTag_dt);

		h_npe_n[idx][0]->Fill(nPESum,nTag_n);
			
		for(size_t j=0;j<nTag_n;++j){
			float e = (*nTag_e)[j];
			float dt = (*nTag_dt)[j];
			h_npe_e[idx][0]->Fill(nPESum,e);
			h_npe_dt[idx][0]->Fill(nPESum,dt);
			h_n_e[idx][0]->Fill(nTag_n,e);
			h_n_dt[idx][0]->Fill(nTag_n,dt);
			h_e_dt[idx][0]->Fill(e,dt);
//			if( e > 6.0 && e < 12.0 && dt > 20 && dt < 200 ){
			if( tagged ){
				h_npe_e[idx][1]->Fill(nPESum,e);
				h_npe_dt[idx][1]->Fill(nPESum,dt);
				h_n_e[idx][1]->Fill(nTag_n,e);
				h_n_dt[idx][1]->Fill(nTag_n,dt);
				h_e_dt[idx][1]->Fill(e,dt);
				//h_e[idx][r_idx][1]->Fill(e);
				//h_dt[idx][r_idx][1]->Fill(dt);
			}else{
				h_npe_e[idx][2]->Fill(nPESum,e);
				h_npe_dt[idx][2]->Fill(nPESum,dt);
				h_n_e[idx][2]->Fill(nTag_n,e);
				h_n_dt[idx][2]->Fill(nTag_n,dt);
				h_e_dt[idx][2]->Fill(e,dt);
				
				//h_e[idx][r_idx][2]->Fill(e);
				//h_dt[idx][r_idx][2]->Fill(dt);
			}
		}


		if( tagged )
			h_npe_n[idx][1]->Fill(nPESum,nTag_n);
		else
			h_npe_n[idx][2]->Fill(nPESum,nTag_n);

		h_npe_wp_n[idx][0]->Fill(nPESum,wpTag_n);
			
		for(size_t j=0;j<wpTag_n;++j){
			int nHit = (*wpTag_nHit)[j];
			float dt = (*wpTag_dt)[j];
			h_npe_wp_nHit[idx][0]->Fill(nPESum,nHit);
			h_npe_wp_dt[idx][0]->Fill(nPESum,dt);
			if( tagged ){
				h_npe_wp_nHit[idx][1]->Fill(nPESum,nHit);
				h_npe_wp_dt[idx][1]->Fill(nPESum,dt);
			}else{
				h_npe_wp_nHit[idx][2]->Fill(nPESum,nHit);
				h_npe_wp_dt[idx][2]->Fill(nPESum,dt);
			}
		}

		if( tagged )
			h_npe_wp_n[idx][1]->Fill(nPESum,wpTag_n);
		else
			h_npe_wp_n[idx][2]->Fill(nPESum,wpTag_n);

		//h_n[idx][r_idx][0]->Fill(nTag_n);
		//h_n[idx][r_idx][1]->Fill(pass_n);
		//h_n[idx][r_idx][2]->Fill(nTag_n-pass_n);
	}
/*
	for(size_t i=0;i<3;++i){
		for(size_t j=0;
		h_npe_n[i]->Write(); 
		h_npe_e[i] ->Write();
		h_npe_dt[i]->Write();
		h_n_e[i]->Write(); 
		h_n_dt[i]->Write(); 
		h_e_dt[i]->Write(); 
	}
*/
	for(size_t i=0;i<3;++i){
		for(size_t j=0;j<3;++j){
			h_npe_n[i][j]->Write();
			h_npe_e[i][j]->Write();
			h_npe_dt[i][j]->Write();
			h_npe_wp_n[i][j]->Write();
			h_npe_wp_nHit[i][j]->Write();
			h_npe_wp_dt[i][j]->Write();

			h_n_e[i][j]->Write();
			h_n_dt[i][j]->Write();
			h_e_dt[i][j]->Write();
		}
	}


	f_out->Close();
}

bool nTagged(vector<float> *nTag_e, vector<float> *nTag_dt){
	
	bool tagged = false;

	size_t nd = nTag_e->size();

	for(size_t i=0; i<nd && !tagged; ++i){

		bool nTag_e_cut = (*nTag_e)[i] > 1.5;
		bool nTag_dt_cut = (*nTag_dt)[i] > 20.0 && (*nTag_dt)[i] < 800.0;

		//tagged = nTag_e_cut && nTag_dt_cut && ((end-start) > multCut);
		tagged = nTag_e_cut && nTag_dt_cut;
	}
	
	return tagged;
}

bool wpTagged(vector<float> *wpTag_dt){

	bool tagged = false;

	for(size_t i=0; i<wpTag_dt->size()-1;++i){

		tagged = ((*wpTag_dt)[i] - (*wpTag_dt)[i+1] == 0);

		if(tagged)
			break;

	}

	return tagged;
	
}

