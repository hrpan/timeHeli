#include "TChain.h"
#include "range.h"
#include<vector>
#include<iostream>
bool isIsolated(vector<double> *dtlSH, size_t mu, int r_idx, int r_shower);
int findRange(double npe);
using namespace std;
const double maxT = 5000; //ms

const bool gd_only = true;

const double nTag_e_min = 1.8; //mev
const double nTag_dt_min = 20; //us
const double nTag_dt_max = 400; //us
const double iso_dt_cut = 1.0; //ms
const float distCut = 500; //mm

void calMuRates(){
	TChain *chain = new TChain("Heli");

	chain->Add("../p17b/data_heli/*.root");

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
	
	double _dtlSH[3][n_range];
	double _dtlSH_sdt[3][n_range];
	double _dtlSH_dt[3][n_range];	
	double _dtlSH_nPESum[3][n_range];	
	int _dtlSH_tag_n[n_range];
	int _nSH[3][n_range];


    size_t n_mus[3][n_range];
    size_t cands[3];

    for(int s=0;s<3;++s){
        for(int r=0;r<n_range;++r)
            n_mus[s][r] = 0;
        cands[s] = 0;
    }

	int r_shower = findRange(4e5);
	for(size_t i=0;i<entries;++i){
		if(i%(entries/10000)==0){
			printf("\r%d/%d(%.3f%%)",i,entries,float(i*100)/entries);
			fflush(stdout);
		}

		chain->GetEntry(i);
        /*         
		for(size_t j=0;j<n_range;++j){
			for(size_t k=0;k<3;++k){
				_dtlSH[k][j] = -1;
				_dtlSH_nPESum[k][j] = -1;
				_nSH[k][j] = 0;
			}
			_dtlSH_tag_n[j] = -1;
		}	
        */
        if(site==4)
            site=3;

        --site;

		++cands[site];

		for(size_t mu = 0; mu < dtlSH->size(); ++mu){
		
			if( (*dtlSH)[mu] < 0 ) break;

			int r_idx = findRange((*nPESum)[mu]);
			
			if(!isIsolated(dtlSH, mu, r_idx, r_shower)) continue;

            ++n_mus[site][r_idx];

		}

	
	}
    for(int s=0;s<3;++s){
        for(int r=0;r<n_range;++r)
            cout << n_mus[s][r] << " ";
        cout << endl << cands[s] << endl;
    }

}

bool isIsolated(vector<double> *dtlSH, size_t mu, int r_idx, int r_shower){
	if(r_idx < r_shower){
		bool pass;
		
		if(mu > 0)
			pass = (*dtlSH)[mu-1] - (*dtlSH)[mu] > iso_dt_cut;
		else 
			pass = (*dtlSH)[mu] < maxT - iso_dt_cut;

		if(!pass) return false;

		if(mu + 1 < dtlSH->size())
			pass = (*dtlSH)[mu] - (*dtlSH)[mu+1] > iso_dt_cut;
		else
			pass = (*dtlSH)[mu] > iso_dt_cut;

		if(!pass) return false;
	}

	return true;
}


int findRange(double npe){

	for(int i=n_range-1;i>=0;--i){
		if( npe > range[i] )
			return i;
		
	}
	return 0;

}
