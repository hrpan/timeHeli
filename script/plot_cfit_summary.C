#include "range.h"

int site, _range;
int fitMin, fitMax, fixB12, fixHe8, bound_eps, fix_lifetime, fix_rmu, eps_pull;
double r_mu_tag[2];
double r_mu_atag[2];
double n_dc[2], n_lihe[2], eps_lihe[2], r_lihe[2], n_bo[2], eps_bo[2];
const char *titles[3] = {
	"Low",
	"Mid",
	"High"
};
const char *prefix = "./plots/summary/";
void plot_cfit_summary(){
	gROOT->ProcessLine("#include <vector>");
	string input_file;
	cin >> input_file;
	cout << "INPUT FILE:" << input_file << endl;
	TFile *f = new TFile(input_file.c_str());
	TTree *tr = f->Get("tr");

	tr->SetBranchAddress("fitMin", &fitMin);
	tr->SetBranchAddress("fitMax", &fitMax);
	tr->SetBranchAddress("fixB12", &fixB12);
	tr->SetBranchAddress("fixHe8", &fixHe8);
	tr->SetBranchAddress("bound_eps", &bound_eps);
	tr->SetBranchAddress("fix_lifetime", &fix_lifetime);
	tr->SetBranchAddress("fix_rmu", &fix_rmu);
	tr->SetBranchAddress("eps_pull", &eps_pull);
	tr->SetBranchAddress("site", &site);
	tr->SetBranchAddress("range", &_range);
	tr->SetBranchAddress("r_mu_tag", &r_mu_tag[0]);
	tr->SetBranchAddress("r_mu_tag_err", &r_mu_tag[1]);
	tr->SetBranchAddress("r_mu_atag", &r_mu_atag[0]);
	tr->SetBranchAddress("r_mu_atag_err", &r_mu_atag[1]);
	tr->SetBranchAddress("n_dc", &n_dc[0]);
	tr->SetBranchAddress("n_dc_err", &n_dc[1]);
	tr->SetBranchAddress("n_lihe", &n_lihe[0]);
	tr->SetBranchAddress("n_lihe_err", &n_lihe[1]);
	tr->SetBranchAddress("eps_lihe", &eps_lihe[0]);
	tr->SetBranchAddress("eps_lihe_err", &eps_lihe[1]);
	tr->SetBranchAddress("r_lihe", &r_lihe[0]);
	tr->SetBranchAddress("r_lihe_err", &r_lihe[1]);
	tr->SetBranchAddress("n_bo", &n_bo[0]);
	tr->SetBranchAddress("n_bo_err", &n_bo[1]);
	tr->SetBranchAddress("eps_bo", &eps_bo[0]);
	tr->SetBranchAddress("eps_bo_err", &eps_bo[1]);
	
	char *prefix = "./plots/summary/";
	char buf[255], filename[255];
	vector<int> _site;
	vector<int> __range;
	vector<int> _fitMin;
	vector<int> _eps_pull;
	vector<double> _eps_lihe;
	vector<double> _eps_lihe_err;
	vector<double> _n_lihe;
	vector<double> _n_lihe_err;

	for(int i=0;i<tr->GetEntries();++i){
		tr->GetEntry(i);
		_site.push_back(site-1);
		__range.push_back(_range);
		_fitMin.push_back(fitMin);
		_eps_pull.push_back(eps_pull);	
		_eps_lihe.push_back(eps_lihe[0]);	
		_eps_lihe_err.push_back(eps_lihe[1]);	
		_n_lihe.push_back(n_lihe[0]);	
		_n_lihe_err.push_back(n_lihe[1]);	
	}
	plotMultiGraph(
		_site,
		__range,
		_fitMin, 
		_eps_lihe, 
		_eps_lihe_err, 
		_eps_pull, 
		"neutron tagging efficiency",	
		"minimum fitting range[ms]",
		"efficiency",
		"pull",
		"eps_lihe_fitMin");
	plotMultiGraph(
		_site,
		__range,
		_fitMin, 
		_n_lihe, 
		_n_lihe_err, 
		_eps_pull, 
		"number of Li9/He8",	
		"minimum fitting range[ms]",
		"N_Li9/He8",
		"pull",
		"n_lihe_fitMin");

	TGraphErrors *gr[3];
	TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle("Shower muon tagging efficiency");
	char buf[255];
	vector<double> _eff_tmp, _eff_err_tmp;
	for(int i=0;i<3;++i){
		vector<double> x, y, y_err;
		for(int j=0;j<_fitMin.size();++j){
			if(i != _site[j]) continue;
			if(__range[j] != 2) continue;
			if(_eps_pull == 1) continue;
			_eff_tmp.push_back(_eps_lihe[j]);
			_eff_err_tmp.push_back(_eps_lihe_err[j]);
			x.push_back(_fitMin[j]+i);		
			y.push_back(_eps_lihe[j]);
			y_err.push_back(_eps_lihe_err[j]);
		}
		gr[i] = new TGraphErrors(x.size(), &x[0], &y[0], 0, &y_err[0]);
		gr[i]->SetMarkerStyle(8);
		gr[i]->SetMarkerColor(i + 2);
		gr[i]->SetMarkerSize(1);
		sprintf(buf, "EH%d", i+1);
		gr[i]->SetTitle(buf);
		gr[i]->SetFillStyle(0);
		gr[i]->SetFillColor(0);
		mg->Add(gr[i]);
	}
	mg->Draw("AP");
	mg->GetXaxis()->SetTitle("minimum fitting range[ms]");
	mg->GetYaxis()->SetTitleOffset(1.5);
	mg->GetYaxis()->SetTitle("Tagging efficiency");
	c1->BuildLegend(0.9, 0.6, 0.97, 0.75);
	sprintf(buf,"%s/tag_eff.png",prefix);
	c1->SaveAs(buf);
/*
	for(int i=0;i<3;++i){
		for(int j=0;j<n_range;++j){
			sprintf(buf,"site==%d&&range==%d", i+1, j);
				
			//tr->Draw("n_lihe", buf);
			//sprintf(filename,"%sEH%d_%d_n_lihe.png", prefix, i+1, j);
			//c1->SaveAs(filename);

			tr->Draw("n_bo:fixHe8+2*fixB12+4*fix_rmu+8*eps_pull", buf, "colz");
			sprintf(filename,"%sEH%d_%d_n_bo_fixs.png", prefix, i+1, j);
			c1->SaveAs(filename);

			tr->Draw("n_lihe:fixHe8+2*fixB12+4*fix_rmu+8*eps_pull", buf, "colz");
			sprintf(filename,"%sEH%d_%d_n_lihe_fixs.png", prefix, i+1, j);
			c1->SaveAs(filename);

			tr->Draw("n_lihe_err:fixHe8+2*fixB12+4*fix_rmu+8*eps_pull", buf, "colz");
			sprintf(filename,"%sEH%d_%d_n_lihe_err_fixs.png", prefix, i+1, j);
			c1->SaveAs(filename);

			tr->Draw("n_lihe:n_bo", buf, "colz");
			sprintf(filename,"%sEH%d_%d_n_lihe_n_bo.png", prefix, i+1, j);
			c1->SaveAs(filename);

			tr->Draw("r_lihe:fixHe8+2*fixB12+4*fix_rmu+8*eps_pull", buf, "colz");
			sprintf(filename,"%sEH%d_%d_r_lihe_fixs.png", prefix, i+1, j);
			c1->SaveAs(filename);
			tr->Draw("eps_lihe:fixHe8+2*fixB12+4*fix_rmu+8*eps_pull", buf, "colz");
			sprintf(filename,"%sEH%d_%d_eps_lihe_fixs.png", prefix, i+1, j);
			c1->SaveAs(filename);
		}
	}
*/
}

void plotMultiGraph(vector<int> &site, vector<int> &range, vector<int> &x_axis, vector<double> &y_axis, vector<double> &y_axis_err, vector<int> &opt, char *graph_title, char *x_axis_title, char *y_axis_title, char *opt_title, char *output_title){
	TGraphErrors *gr[3][n_range][2];
	vector<double> x[2], y[2], y_err[2];
	char buf[255], filename[255];
	for(int s=0;s<3;++s){
		for(int r=0;r<n_range;++r){
			vector<double> x[2], y[2], y_err[2];
			for(int i=0;i<x_axis.size();++i){
				int _opt = opt[i];
				if(site[i] != s || range[i] != r) continue;	
				//cout << x_axis[i] << " " << y_axis[i] << " " << y_axis_err[i] << endl;
				x[_opt].push_back(x_axis[i]+r);		
				y[_opt].push_back(y_axis[i]);		
				y_err[_opt].push_back(y_axis_err[i]);		
			}
			for(int _opt=0;_opt<2;++_opt){
				gr[s][r][_opt] = new TGraphErrors(x[_opt].size(), &x[_opt][0], &y[_opt][0], 0, &y_err[_opt][0]);
				gr[s][r][_opt]->SetMarkerStyle(8);
				gr[s][r][_opt]->SetMarkerColor(r + 2);
				gr[s][r][_opt]->SetMarkerSize(1);
				sprintf(buf, "%s", titles[r]);
				gr[s][r][_opt]->SetTitle(buf);
				gr[s][r][_opt]->SetFillStyle(0);
				gr[s][r][_opt]->SetFillColor(0);
			}
		}
		for(int _opt=0;_opt<2;++_opt){
			TMultiGraph *mg = new TMultiGraph();
			for(int r=0;r<n_range;++r)
				mg->Add(gr[s][r][_opt]);
			
			if(_opt == 0)
				sprintf(buf, "EH%d %s (without %s)", s+1, graph_title, opt_title);
			else
				sprintf(buf, "EH%d %s (with %s)", s+1, graph_title, opt_title);
			mg->SetTitle(buf);
			mg->Draw("AP");
			mg->GetXaxis()->SetTitle(x_axis_title);
			mg->GetYaxis()->SetTitle(y_axis_title);
			mg->GetYaxis()->SetTitleOffset(1.5);
		//	gPad->Update();
			c1->BuildLegend(0.9, 0.6, 0.97, 0.75);
			sprintf(filename,"%sEH%d_%s_%d.png", prefix, s+1, output_title, _opt);
			c1->SaveAs(filename);
			c1->Clear();
			mg->Delete();
		}


	}
		
}

double max_vec(vector<double> &vec){
	double _max = -1000;
	for(int i=0;i<vec.size();++i)
		_max = _max < vec[i]? vec[i]:_max;
			
	return _max;
}

double min_vec(vector<double> &vec){
	double _min = 1000;
	for(int i=0;i<vec.size();++i)
		_min = _min < vec[i]? _min:vec[i];
			
	return _min;
}

double avg_vec(vector<double> &vec){
	double _sum = 0;
	for(int i=0;i<vec.size();++i)
		_sum += vec[i];

	return _sum / vec.size();
}
