#include "range.h"

const int fix_rmu = 0;
const int fixB12 = 0;
const int fixHe8 = 0;
const int fitMin = 50;

const double livetime[3] = {
	3282.02942397,
	3295.28353276,
	6766.68404126
};

const double eff_sh = 0.21;

void printSummary(){
	string filename;
	cin >> filename;
	TFile *f1 = new TFile(filename.c_str());
	TTree *tr1 = f1->Get("tr");

	double result[3][n_range][3];
	parseTree(tr1, result);	
	
	for(int i=0;i<3;++i){
		double _sum[3] = {0, 0, 0};
		double _sum_eff[3] = {0, 0, 0};
		for(int r=0;r<n_range;++r){
			cout << i << " " << r << " ";
			for(int j=0;j<3;++j)			
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
		
		printf("sum %6.3f/%6.3f/%6.3f (%6.3f/%6.3f/%6.3f)\n", _sum[0], sqrt(_sum[1]), sqrt(_sum[2]), _sum_eff[0], sqrt(_sum_eff[1]), sqrt(_sum_eff[2]));
		cout << "=============" << endl;
	}

}

void parseTree(TTree *tr, double result[3][n_range][3]){
	char buf[255], file[255], title[255];
	for(int i=0;i<3;++i){
		for(int r=0;r<n_range;++r){
//			sprintf(buf,"site==%d&&range==%d&&fitMin==5&&fix_rmu==1",i+1,r);
//			sprintf(buf,"site==%d&&range==%d&&fix_rmu==%d&&fixB12==%d&&fixHe8==%d&&fitMin<%d",i+1,r,fix_rmu,fixB12,fixHe8,fitMin);
//			sprintf(buf,"site==%d&&range==%d&&fixHe8==1&&eps_pull==1",i+1,r);
			sprintf(buf,"site==%d&&range==%d&&eps_pull==1",i+1,r);

			sprintf(title,"EH%d range%d N_Li",i+1,r); 
			TH1D *h = new TH1D("h",title,100,0,5e4);
			tr->Draw("n_lihe>>h", buf);
			sprintf(file,"./plots/cfits_summary/EH%d_%d_n_li.png",i+1,r);
//			c1->SaveAs(file);
			result[i][r][0] = h->GetMean() / livetime[i];
			result[i][r][2] = h->GetStdDev() / livetime[i];
			sprintf(title,"EH%d range%d N_Li_err^2",i+1,r); 
			TH1D *h_err = new TH1D("h_err",title,100,0,1e13);
			tr->Draw("n_lihe_err * n_lihe_err>>h_err", buf);
			sprintf(file,"./plots/cfits_summary/EH%d_%d_n_li_err2.png",i+1,r);
//			c1->SaveAs(file);
			result[i][r][1] = sqrt(h_err->GetMean()) / livetime[i];

			h->Delete();
			h_err->Delete();
		}	
	}	
	

}
