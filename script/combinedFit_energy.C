#include "util.C"
#include "range.h"

const char *hists_prefix = "./hists_e_slices";

const double tau_li9 = 256.366;
const double tau_li9_err = 0.866;

const double tau_he8 = 171.68;

const double tau_b12 = 29.142;
const double tau_b12_err = 0.0288;

double fitMin = 1.5;
double fitMax = 5000;

const double simplex_args[2] = {1000000, 0.1};
const double migrad_args[2] = {1000000, 0.1};
const double minos_args[1] = {1000000};

const bool withB12 = true;

const int n_muon_rates = 3 * n_range;
const int e_slices = 20;
const int n_slice = e_slices + 1;

const int npars = 3 * n_range + 3 * e_slices + 1 + 1;
/*
		n_range //muon rates  
		n_range //n_lihe
		n_range //n_b12
		e_slices-1  //dc effs
		e_slices-1  //lihe effs
		e_slices-1  //b12 effs
		1 //n_tot
		1 //lihe ratio
*/
char par_names[npars][255];
const int npars_single = 8;
//[0]:mu rate
//[1]:n_dc
//[2]:n_li9he8
//[3]:t_li9
//[4]:t_he9
//[5]:ratio_lihe
//[6]:n_b12
//[7]:t_b12
const int n_pdfs = 3;
const char *_pdfs[n_pdfs] = {
	"[1] * [0] * exp(-[0] * x )",
	"[2] * ([5] * ([0] + 1 / [3]) * exp(-([0] + 1 / [3]) * x ) + (1 - [5]) * ([0] + 1 / [4]) * exp(-([0] + 1 / [4]) * x) )",
	"[6] * ([0] + 2 / [7]) * exp(-([0] + 2 / [7]) * x )",
};


void printResults(double results[3][npars][2]){
	char buf[512];
	printf("%15s %21s %21s %21s\n", "site", "EH1", "EH2", "EH3");
	for(int i=0;i<npars;++i){
		printf("%15s", par_names[i]);
		for(int j=0;j<3;++j)
			printf(" %10.3e/%10.3e", results[j][i][0], results[j][i][1]);
		cout << endl;
	}

}


ROOT::Math::WrappedMultiTF1 *wf[n_range][n_slice];
ROOT::Math::IMultiGenFunction *fPNLL[n_range][n_slice];

TF1 *func[n_range][n_slice];

TH1 *h[n_range][n_slice];

void parTrans(const double *par, double _par[n_range][n_slice][npars_single]){
/*
		n_range //muon rates  
		n_range //n_lihe
		n_range //n_b12
		( e_slices - 1 ) //dc effs
		( e_slices - 1 ) //lihe effs
		( e_slices - 1 ); //b12 effs
		1 //n_tot
		1 //lihe ratio
*/
	double scale = h[0][0]->GetBinWidth(1);

	double r_mu[n_range];
	double n_lihe[n_range];
	double n_b12[n_range];
	double eps_dc[e_slices];
	double eps_lihe[e_slices];
	double eps_b12[e_slices];
	double r_lihe;
	double n_dc;
	//double r_lihe;

	int p_idx = 0;
	for(int i=0;i<n_range;++i){
		r_mu[i] = par[p_idx]; p_idx++;
		n_lihe[i] = par[p_idx] * scale; p_idx++;
		n_b12[i] = par[p_idx] * scale; p_idx++;
	}

	for(int i=0;i<e_slices;++i){
		eps_dc[i] = par[p_idx]; p_idx++;
		eps_lihe[i] = par[p_idx]; p_idx++;
		eps_b12[i] = par[p_idx]; p_idx++;
	}
	r_lihe = par[p_idx]; p_idx++;

	n_dc = par[p_idx] * scale; p_idx++;

	for(int i = 0; i < n_range; ++i){
		double n_dc_fake_lihe = 0;
		double n_dc_fake_b12 = 0;
		for(int j=0;j<n_range;++j){
			if(i==j) continue;
			n_dc_fake_lihe += n_lihe[j];
			n_dc_fake_b12 += n_b12[j];
		}
		for(int j = 0; j < n_slice; ++j){
			_par[i][j][0] = r_mu[i];
			if(j==0){	
				_par[i][j][1] = n_dc + (n_dc_fake_lihe + n_dc_fake_b12);
				_par[i][j][2] = n_lihe[i];
				_par[i][j][6] = n_b12[i];
			}else{
				_par[i][j][1] = n_dc * eps_dc[j-1] + (n_dc_fake_lihe * eps_lihe[j-1] + n_dc_fake_b12 * eps_b12[j-1]);
				_par[i][j][2] = eps_lihe[j-1] * n_lihe[i];
				_par[i][j][6] = eps_b12[j-1] * n_b12[i];
			}
			_par[i][j][3] = tau_li9;
			_par[i][j][4] = tau_he8;
			_par[i][j][5] = r_lihe;
			_par[i][j][7] = tau_b12;
			//for(int k=0;k<8;++k)
			//	cout << _par[i][j][k] << " ";
			//cout << endl;
		}
	}
//	for(int i=0;i<npars;++i)
//		cout << par[i] << endl;
}

double max_diff = 0;

double _likelihood(const double *par){
	double _par[n_range][n_slice][npars_single];

	parTrans(par, _par);
	double result = 0;
	vector<double> _ll(n_slice * n_range);	
	for(int i=0;i<n_range;++i){
		for(int j=0;j<n_slice;++j){
			double _tmp = (*fPNLL[i][j])(_par[i][j]);
			_ll[i * n_slice + j] = _tmp;
		}
	}
	result = neumaierSum(&_ll[0], n_slice * n_range);
	//cout << result << endl;
	//if(max_diff < 1e-10)
	//	getchar();
	//printf("result: %20.15e %20.15e %20.15e %20.15e\n", result, safe_sum, max_diff, fabs(result-safe_sum));
	return result;
}

void wrap(int &npar, double *g, double &result, double *par, int flag){
	result = _likelihood(par);
}

void fillPars(double *par, TF1 *f[n_range][n_slice]){

	double _par[n_range][n_slice][npars_single];

	parTrans(par, _par);

	for(int i=0;i<n_range;++i)
		for(int j=0;j<n_slice;++j)
			f[i][j]->SetParameters(_par[i][j]);

}

void plotHists(int site, TH1 *h[n_range][n_slice], TF1 *f[n_range][n_slice], double *par){
	fillPars(par, func);
	char *dir = "./plots/cfits_energy";
	char buf[255];
	double _par[n_range][n_slice][npars_single];
	parTrans(par, _par);	
	for(int i=0;i<n_range;++i){
		for(int j=0;j<n_slice;++j){
			TF1 ibd("ibd","[1] * [0] * exp(-[0] * x)", 0, 5000);
			ibd.SetParameters(_par[i][j]);
			ibd.SetLineColor(kGreen);
			sprintf(buf, "%s/cfit_%d_%d_%d.png", dir, site, i, j);
			h[i][j]->Draw("E!");
			f[i][j]->Draw("same");
			ibd.Draw("same");
			gPad->SetLogx();
			c1->SaveAs(buf);
		}
	}

}

void initialize_minimizer(int site, ROOT::Math::Minimizer *minim){
	minim->SetErrorDef(0.5);
	minim->SetTolerance(0.01);
	minim->SetPrintLevel(2);
	minim->SetStrategy(2);
	minim->SetMaxFunctionCalls(1000000);
/*
		n_range //muon rates  
		n_range //n_lihe
		n_range //n_b12
		e_slices //dc effs
		e_slices //lihe effs
		e_slices //b12 effs
		1 //n_tot
		1 //lihe ratio
*/

	double step_ratio = 1e-2;
	double mu_rate_init = 1e-4;
	double n_lihe_init = 1e3;
	double n_b12_init = 1e1;
	double eps_init = 1.0 / e_slices;
	double r_lihe_init = 0.5;	
	char buf[255];
	int p_idx = 0;

	for(int i=0;i<n_range;++i){

		h[i][0]->Fit("expo","LQN0","", 300, 5000);
		mu_rate_init = fabs(expo->GetParameter(1));
		sprintf(buf, "mu_rate_%d", i);
		sprintf(par_names[p_idx], "%s", buf);
		minim->SetVariable(p_idx, buf, mu_rate_init, mu_rate_init * step_ratio);
		p_idx++;

		sprintf(buf, "n_lihe_%d", i);
		sprintf(par_names[p_idx], "%s", buf);
		minim->SetLowerLimitedVariable(p_idx, buf, n_lihe_init, n_lihe_init * step_ratio, 0);
		p_idx++;

		sprintf(buf, "n_b12_%d", i);
		sprintf(par_names[p_idx], "%s", buf);
		if(withB12)
			minim->SetLowerLimitedVariable(p_idx, buf, n_b12_init, n_b12_init * step_ratio, 0);
		else
			minim->SetFixedVariable(p_idx, buf, 0);
		
		p_idx++;
	}
	for(int i=0;i<e_slices;++i){
		eps_init = h[0][i]->GetEntries() / h[0][0]->GetEntries();
		sprintf(buf, "eps_dc_%d", i);
		sprintf(par_names[p_idx], "%s", buf);
		minim->SetLimitedVariable(p_idx, buf, eps_init, eps_init * step_ratio, 0, 1);
		p_idx++;

		sprintf(buf, "eps_lihe_%d", i);
		sprintf(par_names[p_idx], "%s", buf);
		minim->SetLimitedVariable(p_idx, buf, eps_init, eps_init * step_ratio, 0, 1);
		p_idx++;

		sprintf(buf, "eps_b12_%d", i);
		sprintf(par_names[p_idx], "%s", buf);
		if(withB12)
			minim->SetLimitedVariable(p_idx, buf, eps_init, eps_init * step_ratio, 0, 1);
		else
			minim->SetFixedVariable(p_idx, buf, 0);
		
		p_idx++;
	}

	sprintf(buf,"r_lihe",i);
	sprintf(par_names[p_idx], "%s", buf);
	minim->SetLimitedVariable(p_idx, buf, r_lihe_init, r_lihe_init * step_ratio, 0, 1);
	p_idx++;

	double n_dc_init = h[0][0]->GetEntries();
	sprintf(par_names[p_idx], "n_dc");
	minim->SetVariable(p_idx, "n_dc", n_dc_init, n_dc_init * step_ratio);
	p_idx++;
	
}

void combinedFit_energy(){

	char _f_sum[512];
	for(int i=0;i<n_pdfs;++i)
		if(i==0)
			sprintf(_f_sum,"%s",_pdfs[i]);
		else
			sprintf(_f_sum,"%s+%s",_f_sum,_pdfs[i]);

	cout << "INITIALIZING FUNCTION WRAPPERS" << endl;
	char buf[255];	
	for(int r=0;r<n_range;++r){
		for(int s=0;s<n_slice;++s){
			sprintf(buf,"func%d_%d",r,s);
			func[r][s] = new TF1(buf, _f_sum, 0, 5000);
			wf[r][s] = new ROOT::Math::WrappedMultiTF1(*func[r][s], 1);
		}
	}
	
	double results[3][npars][2];
	
	ROOT::Fit::DataOptions opt;
	opt.fUseEmpty = true;
	
	for(int site=0;site<3;++site){
		do_fit(site, opt, results);
	}
	//printLatexTable(results);
	printResults(results);

}

void do_fit(int site, ROOT::Fit::DataOptions &opt, double results[3][npars][2]){
	char buf[255];
	ROOT::Fit::BinData *data[n_range][n_slice];
	ROOT::Fit::DataRange rangeD(fitMin, fitMax);
	for(int r=0;r<n_range;++r){
		for(int s=0;s<n_slice;++s){
			if(s==0)
				sprintf(buf, "%s/EH%d_dtlSH_%d.root", hists_prefix, site+1, r);
			else
				sprintf(buf, "%s/EH%d_dtlSH_%d_%d.root", hists_prefix, site+1, r, s-1);
			TFile *f_tmp = new TFile(buf, "READ");
			h[r][s] = (TH1*)f_tmp->Get("h");			
			data[r][s] = new ROOT::Fit::BinData(opt, rangeD);
			ROOT::Fit::FillData(*data[r][s], h[r][s]);	
			fPNLL[r][s] = new ROOT::Fit::PoissonLikelihoodFCN<ROOT::Math::IBaseFunctionMultiDim>(*data[r][s], *wf[r][s]);
		}
	}

	ROOT::Math::Minimizer *minim = ROOT::Math::Factory::CreateMinimizer("Minuit2");
	ROOT::Math::Functor _func(_likelihood, npars);
	minim->SetFunction(_func);
	initialize_minimizer(site, minim);
	minim->Hesse();
	minim->Minimize();
	minim->Hesse();
	minim->PrintResults();

	TH2D *h_corr = new TH2D("h_corr","Correlations",npars, 0, npars, npars, 0, npars);
	for(int i=0;i<npars;++i)
		for(int j=0;j<npars;++j)
			h_corr->SetBinContent(i,j,minim->Correlation(i,j));
	h_corr->Draw("colz");
	gPad->SetLogx(0);
	gPad->SetLogy(0);
	sprintf(buf,"./plots/cfits_energy/corr%d.png",site+1);
	c1->SaveAs(buf);

	plotHists(site, h, func, minim->X());
	plotEnergy(site, minim->X(), minim->Errors());

	double *_x = minim->X();
	double *_x_err = minim->Errors();
	for(int i=0;i<npars;++i){
		results[site][i][0] = _x[i];
		results[site][i][1] = _x_err[i];
	}
}

void plotEnergy(int site, double *par, double *par_err){
	double x[e_slices], x_err[e_slices], y[3][e_slices], y_err[3][e_slices];

	int offset = 3 * n_range;
	for(int i=0;i<e_slices;++i){
		x[i] = 0.7 + (12-0.7)/e_slices * (i+1);
		x_err[i] = (12-0.7) / ( 2 * e_slices );
		y[0][i] = par[offset + 3*i];
		y_err[0][i] = par_err[offset + 3*i];
		y[1][i] = par[offset + 3*i+1];
		y_err[1][i] = par_err[offset + 3*i+1];
		y[2][i] = par[offset + 3*i+2];
		y_err[2][i] = par_err[offset + 3*i+2];
	}
	TMultiGraph *mg = new TMultiGraph();
	TGraph *gr[3];
	for(int i=0;i<2;++i){
		gr[i] = new TGraphErrors(e_slices, x, y[i], x_err, y_err[i]);
		gr[i]->SetMarkerStyle(20);
		gr[i]->SetMarkerColor(i+2);
		mg->Add(gr[i]);
	}
	mg->SetMinimum(0);
	mg->Draw("APL");
	char buf[255];
	char *dir = "./plots/cfits_energy";
	sprintf(buf,"%s/cfit_e_%d.png",dir,site+1);
	gPad->SetLogx(0);
	gPad->SetLogy(0);
	c1->SaveAs(buf);
}
