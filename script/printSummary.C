#include "range.h"

const int fix_rmu = 0;
const int fixB12 = 0;
const int fixHe8 = 0;
const int fitMin = 50;

const double livetime[3][4] = {
	{1544.17, 1737.38, 0, 0},
	{1741.12, 1554.04, 0, 0},
	{1739.26, 1739.26, 1739.26, 1551.78}
};

const double eff_muon[3][4] = {
	{0.8227, 0.8194, 0, 0},
	{0.8518, 0.8506, 0, 0},
	{0.9832, 0.9829, 0.9828, 0.9833}
};

const double eff_mult[3][4] = {
	{0.9495, 0.9502, 0, 0},
	{0.9520, 0.9518, 0, 0},
	{0.9525, 0.9525, 0.9523, 0.9524}
};

const double eff_sh = 0.21;

void printSummary(){
	string filename;
	cin >> filename;
	TFile *f1 = new TFile(filename.c_str());
	TTree *tr1 = f1->Get("tr");

	double result[3][n_range][4];

	double lt[3] = {0, 0, 0};
	for(int i=0;i<3;++i){
		for(int j=0;j<4;++j){
			lt[i] += livetime[i][j] * eff_muon[i][j] * eff_mult[i][j];
		}
	}
	
	parseTree(tr1, result, lt);	

	for(int i=0;i<3;++i){
		double _sum[3] = {0, 0, 0};
		double _sum_eff[3] = {0, 0, 0};
		for(int r=0;r<n_range;++r){
			cout << i << " " << r << " ";
			for(int j=0;j<4;++j)			
				printf("%6.3f ", result[i][r][j]);
			_sum[0] += result[i][r][0];
			_sum[1] += result[i][r][1] * result[i][r][1];
			_sum[2] += result[i][r][2] * result[i][r][2];
			if(r == n_range-1){
				_sum_eff[0] += eff_sh * result[i][r][0];
				double _tmp = eff_sh * result[i][r][1];
				_sum_eff[1] += _tmp * _tmp;
				_tmp = eff_sh * result[i][r][2];
				_sum_eff[2] += _tmp * _tmp;
			}else{
				_sum_eff[0] += result[i][r][0];
				_sum_eff[1] += result[i][r][1] * result[i][r][1];
				_sum_eff[2] += result[i][r][2] * result[i][r][2];
			}
			
			cout << endl;
		}
		
		printf("sum %6.3f/%6.3f/%6.3f/%6.3f (%6.3f/%6.3f/%6.3f/%6.3f)\n", _sum[0], sqrt(_sum[1]), sqrt(_sum[2]), sqrt(_sum[1] + _sum[2]), _sum_eff[0], sqrt(_sum_eff[1]), sqrt(_sum_eff[2]), sqrt(_sum_eff[1] + _sum_eff[2]));
		cout << "=============" << endl;
	}

}

void parseTree(TTree *tr, double result[3][n_range][4], double lt[3]){
	char buf[255], file[255], title[255];
	for(int i=0;i<3;++i){
		for(int r=0;r<n_range;++r){
//			sprintf(buf,"site==%d&&range==%d&&fixHe8==1&&eps_pull==1",i+1,r);
			sprintf(buf,"site==%d&&range==%d&&eps_pull==1",i+1,r);

			sprintf(title,"EH%d range%d N_Li",i+1,r); 
			TH1 *h;
			tr->Draw("n_lihe", buf);
			h = (TH1*) gPad->GetPrimitive("htemp");	
			result[i][r][0] = h->GetMean() / lt[i];
			result[i][r][2] = h->GetStdDev() / lt[i];
			sprintf(title,"EH%d range%d N_Li_err^2",i+1,r); 
			tr->Draw("n_lihe_err", buf);
			h = (TH1*) gPad->GetPrimitive("htemp");	
			cout << h->GetMean() << " " << h->GetRMS() << endl;
			result[i][r][1] = h->GetMean() / lt[i];
			result[i][r][3] = sqrt(result[i][r][1] * result[i][r][1] + result[i][r][2] * result[i][r][2]);			
		}	
	}	
	

}
