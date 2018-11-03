#include "util.C"
#include "range.h"
#include "slice.h"

const char *hists_prefix = "./hists";

const double tau_li9 = 256.366;
const double tau_li9_err = 0.866;

const double tau_he8 = 171.68;

const double tau_b12 = 29.142;
const double tau_b12_err = 0.0288;

double fitMin = 1.5;
double fitMax = 5000;

const int npars_max = 200;
const bool use_tagging = false;

char par_names[npars_max][255]; 

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

void printResults(double results[3][npars_max][2]){
	char buf[512];
	printf("%15s %21s %21s %21s\n", "site", "EH1", "EH2", "EH3");
	for(int i=0;i<npars_max;++i){
		printf("%15s", par_names[i]);
		for(int j=0;j<3;++j)
			printf(" %10.3e/%10.3e", results[j][i][0], results[j][i][1]);
		cout << endl;
		if(par_names[i+1][0] == 0) break;
	}

}

ROOT::Math::IMultiGenFunction *fPNLL[n_range][3][slice_types][slices_max];
ROOT::Math::WrappedMultiTF1 *wf[n_range][3][slice_types][slices_max];

TF1 *func[n_range][3][slice_types][slices_max];
TH1 *h[n_range][3][slice_types][slices_max];

void parTrans(const double *par, double _par[n_range][3][slice_types][slices_max][npars_single]){

	double scale = h[0][0][0]->GetBinWidth(1);

	double r_mu[n_range];
	double r_mu_tag[n_range];
	double r_mu_atag[n_range];
	double n_lihe[n_range];
	double n_bo[n_range];
	double eps_t_lihe[n_range];
	//double _eps_lihe[eps_slices];
	//double eps_bo[n_range];
	double eps_t_bo;
	double eps_e_dc[e_slices];
	double eps_e_lihe[e_slices];
	double eps_e_bo[e_slices];
	//double eps_bo;
	double n_dc;
	double r_lihe;

	int p_idx = 0;
	for(int r = 0; r < n_range; ++r){
		if(use_tagging){
			r_mu_tag[r] = par[p_idx++];
			r_mu_atag[r] = par[p_idx++];
			n_lihe[r] = par[p_idx++] * scale;
			eps_t_lihe[r] = par[p_idx++];
			n_bo[r] = par[p_idx++] * scale;
		}else{
			r_mu[r] = par[p_idx++];
			n_lihe[r] = par[p_idx++] * scale;
			n_bo[r] = par[p_idx++] * scale;
		}
	}
	if(use_tagging)
		eps_t_bo = par[p_idx++];
	eps_e_dc[0] = 1;
	eps_e_lihe[0] = 1;
	eps_e_bo[0] = 1;

	double _tmp_sum[3] = {1, 1, 1};
	for(int e=0;e<e_slices-1;++e){
		eps_e_dc[e+1] = par[p_idx++];
		_tmp_sum[0] += eps_e_dc[e+1];
		eps_e_lihe[e+1] = par[p_idx++];
		_tmp_sum[1] += eps_e_lihe[e+1];
		eps_e_bo[e+1] = par[p_idx++];
		_tmp_sum[2] += eps_e_bo[e+1];
	}
	for(int e=0;e<e_slices;++e){
		eps_e_dc[e] /= _tmp_sum[0];
		eps_e_lihe[e] /= _tmp_sum[1];
		eps_e_bo[e] /= _tmp_sum[2];
	}
	
	r_lihe = par[p_idx++];
	n_dc = par[p_idx++] * scale;

	for(int r=0; r<n_range; ++r){
		double n_dc_fake_lihe = 0;
		double n_dc_fake_bo = 0;
		for(int j=0;j<n_range;++j){
			if(r==j) continue;
			n_dc_fake_lihe += n_lihe[j];
			n_dc_fake_bo += n_bo[j];
		}
		for(int t=0; t<3; ++t){
			for(int e=0; e<e_slices; ++e){
				if(t==0){
					if(use_tagging)
						_par[r][t][e][0] = r_mu_tag[r] + r_mu_atag[r];
					else
						_par[r][t][e][0] = r_mu[r];
					_par[r][t][e][1] = n_dc * eps_e_dc[e] 
						+ (n_dc_fake_lihe * eps_e_lihe[e] 
						+ n_dc_fake_bo * eps_e_bo[e]);
					_par[r][t][e][2] = n_lihe[r] * eps_e_lihe[e];			
					_par[r][t][e][6] = n_bo[r] * eps_e_bo[e];
				}else if(j==1){
					_par[r][t][e][0] = r_mu_tag[r];
					_par[r][t][e][1] = n_dc * eps_e_dc[e] 
						+ ( (n_dc_fake_lihe + ( 1 - eps_t_lihe[r] ) * n_lihe[i] ) 
						* eps_e_lihe[e] 
						+ ( n_dc_fake_bo + ( 1 - eps_t_bo ) * n_bo[i]) 
						* eps_e_bo[e]);
					_par[r][t][e][2] = eps_t_lihe[r] * eps_e_lihe[e] * n_lihe[r];
					_par[r][t][e][6] = eps_t_bo * eps_e_bo[e] * n_bo[r];
				}else{
					_par[r][t][e][0] = r_mu_atag[r];
					_par[r][t][e][1] = n_dc * eps_e_dc[e] 
						+  ( (n_dc_fake_lihe + eps_t_lihe[r] * n_lihe[r]) 
						* eps_e_lihe[e] 
						+ (n_dc_fake_bo + eps_t_bo * n_bo[r]) 
						* eps_e_bo[e]);
					_par[r][t][e][2] = ( 1 - eps_t_lihe[r] ) * n_lihe[r] * eps_e_lihe[e];
					_par[r][t][e][6] = ( 1 - eps_t_bo ) * n_bo[r] * eps_e_bo[e];
				}
				_par[r][t][e][3] = tau_li9;
				_par[r][t][e][4] = tau_he8;
				_par[r][t][e][5] = r_lihe[e];
				_par[r][t][e][7] = tau_b12;
/*	
				cout << i << " " << j << " " << k << " ";
				for(int l=0;l<npars_single;++l)
					printf("%6e ", _par[i][j][k][l]);
				cout << endl;
*/
			} //e
			if(!use_tagging) break;
		} //t
		
	}

}

double likelihood(const double *par){
	double _par[n_range][3][e_slices][npars_single];

	parTrans(par, _par);
	
	double result = 0;
	vector<double> _nlls(n_range * 3 * e_slices);
	int idx = 0;
	for(int r = 0; r < n_range; ++r){
		for(int t = 0; t < 3; ++t){
			for(int e = 0; e < e_slices; ++e){
				double _tmp = (*fPNLL[r][t][e])(_par[r][t][e]);
				//if(r==3)
				//	cout << r << " " << s << " " << e << " " << _tmp << endl;
				_nlls[ idx++ ] = _tmp;
				result += _tmp;
			}
			if(!use_tagging) break;
		}
	}
	
	double safesum = neumaierSum(&_nlls[0], _nlls.size());
	//cout << result << " " << safesum << " " << result - safesum << endl;
	return safesum;
}


void fillPars(double *par, TF1 *f[n_range][3][e_slices]){

	double _par[n_range][3][e_slices][npars_single];

	parTrans(par, _par);

	for(int r=0;r<n_range;++r)
		for(int t=0;t<3;++t){
			for(int e=0;e<e_slices;++e)
				f[r][t][e]->SetParameters(_par[r][t][e]);
			if(!use_tagging) break;
		}

}

void plotHists(int site, double *par){
	fillPars(par, func);
	char *dir = "./plots/cfits";
	char buf[255];
	double _par[n_range][3][e_slices][npars_single];
	parTrans(par, _par);
	for(int r=0;r<n_range;++r)
		for(int t=0;t<3;++t){
			for(int e=0;e<e_slices;++e){
				TF1 ibd("ibd","[1] * [0] * exp(-[0] * x)",0, 5000);
				ibd.SetParameters(_par[r][t][e]);
				ibd.SetLineColor(kGreen);
				sprintf(buf, "%s/cfit_%d_%d_%d_%d.png", dir, site, r, t, e);
				h[r][t][e]->Draw("E1");
				func[r][t][e]->Draw("same");
				ibd.Draw("same");
				gPad->SetLogx();
				c1->SaveAs(buf);
			}
			if(!use_tagging) break;
		}
	
}

void initialize_minimizer(int site, ROOT::Math::Minimizer *minim){
	cout << "INITIALIZING MINIMIZER" << endl;
	minim->SetErrorDef(0.5);
	minim->SetTolerance(0.01);
	minim->SetPrintLevel(2);
	minim->SetPrecision(1e-16);
	minim->SetStrategy(2);
	minim->SetMaxFunctionCalls(1000000);
	minim->SetMaxIterations(1000000);

	double rmu_init;
	double n_lihe_init = 1e3;
	double n_bo_init = 1e2;
	double eps_init = 0.5;
	double r_lihe_init = 0.8;

	double step_ratio = 1e-1;
	char buf[255];
	int p_idx = 0;
	char *init_str = "Initializing parameter %2d %20s %10.2e %10.2e\n";
	for(int r=0;r<n_range;++r){
		if(use_tagging){		
			h[r][1][0]->Fit("expo","LQN0","", 300, 5000);
			rmu_init = fabs(expo->GetParameter(1));
			sprintf(buf, "rmu_tag_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, rmu_init, rmu_init * step_ratio);
			printf(init_str, p_idx-1, par_names[p_idx-1], rmu_init, rmu_init * step_ratio);

			h[r][2][0]->Fit("expo","LQN0","", 300, 5000);
			rmu_init = fabs(expo->GetParameter(1));
			sprintf(buf, "rmu_atag_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, rmu_init, rmu_init * step_ratio);
			printf(init_str, p_idx-1, par_names[p_idx-1], rmu_init, rmu_init * step_ratio);

			sprintf(buf, "n_lihe_%d",i);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLowerLimitedVariable(p_idx++, buf, n_lihe_init, n_lihe_init * step_ratio , 0);
			printf(init_str, p_idx-1, par_names[p_idx-1], n_lihe_init, n_lihe_init * step_ratio);
			//
			sprintf(buf, "eps_tag_lihe_%d",i);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0, 1);
			printf(init_str, p_idx-1, par_names[p_idx-1], eps_init, eps_init * step_ratio);

			sprintf(buf, "n_bo_%d",i);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLowerLimitedVariable(p_idx++, buf, n_bo_init, n_bo_init * step_ratio , 0);
			printf(init_str, p_idx-1, par_names[p_idx-1], n_bo_init, n_bo_init * step_ratio);
			//minim->SetVariable(p_idx++, buf, n_lihe_init, n_lihe_init * step_ratio);

		}else{
			h[r][0][0]->Fit("expo","LQN0","", 300, 5000);
			rmu_init = fabs(expo->GetParameter(1));
			sprintf(buf, "rmu_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetVariable(p_idx++, buf, rmu_init, rmu_init * step_ratio);
			printf(init_str, p_idx-1, par_names[p_idx-1], rmu_init, rmu_init * step_ratio);

			sprintf(buf, "n_lihe_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLowerLimitedVariable(p_idx++, buf, n_lihe_init, n_lihe_init * step_ratio , 0);
			printf(init_str, p_idx-1, par_names[p_idx-1], n_lihe_init, n_lihe_init * step_ratio);
			
			sprintf(buf, "n_bo_%d",r);
			sprintf(par_names[p_idx], "%s", buf);
			minim->SetLowerLimitedVariable(p_idx++, buf, n_bo_init, n_bo_init * step_ratio , 0);
			printf(init_str, p_idx-1, par_names[p_idx-1], n_bo_init, n_bo_init * step_ratio);

		}
	}
	if(use_tagging){
		sprintf(buf, "eps_tag_bo", i);
		sprintf(par_names[p_idx], "%s", buf);
		minim->SetLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0, 1);
		printf(init_str, p_idx-1, par_names[p_idx-1], eps_init, eps_init * step_ratio);
	}

	for(int e=0;e<e_slices-1;++e){
		eps_init = h[0][0][e+1]->GetEntries() / h[0][0][0]->GetEntries();
		sprintf(buf, "eps_e_dc_%d", e);
		sprintf(par_names[p_idx], "%s", buf);
		minim->SetLowerLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0);
		printf(init_str, p_idx-1, par_names[p_idx-1], eps_init, eps_init * step_ratio);

		sprintf(buf, "eps_e_lihe_%d", e);
		sprintf(par_names[p_idx], "%s", buf);
		minim->SetLowerLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0);
		//minim->SetLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0, 1);
		printf(init_str, p_idx-1, par_names[p_idx-1], eps_init, eps_init * step_ratio);

		sprintf(buf, "eps_e_bo_%d", e);
		sprintf(par_names[p_idx], "%s", buf);
		minim->SetLowerLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0);
		//minim->SetLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0, 1);
		printf(init_str, p_idx-1, par_names[p_idx-1], eps_init, eps_init * step_ratio);
	}
	
	sprintf(par_names[p_idx], "r_lihe");
	minim->SetLimitedVariable(p_idx++, "r_lihe", r_lihe_init, r_lihe_init * step_ratio, 0, 1);
	printf(init_str, p_idx-1, par_names[p_idx-1], r_lihe_init, r_lihe_init * step_ratio);
	
	double n_dc_init = 0;
	for(int e=0;e<e_slices;++e){
		n_dc_init += h[0][0][e]->GetEntries();
	}
	sprintf(par_names[p_idx], "n_dc");
	minim->SetVariable(p_idx++, "n_dc", n_dc_init, n_dc_init * step_ratio);
	//minim->SetLowerLimitedVariable(p_idx++, "n_dc", n_dc_init, n_dc_init * step_ratio, 0);
	printf(init_str, p_idx-1, par_names[p_idx-1], n_dc_init, n_dc_init * step_ratio);

	par_names[p_idx][0] = 0;
}

void combinedFit(){


	gStyle->SetOptFit(1);
	char _f_sum[512];
	for(int i=0;i<n_pdfs;++i)
		if(i==0)
			sprintf(_f_sum,"%s",_pdfs[i]);
		else
			sprintf(_f_sum,"%s+%s",_f_sum,_pdfs[i]);

	cout << "INITIALIZING FUNCTION WRAPPERS" << endl;
	char buf[255];	
	for(int r=0;r<n_range;++r){
		for(int i=0;i<3;++i){
			for(int e=0;e<e_slices;++e){
				sprintf(buf,"func%d_%d_%d",r,i,e);
				func[r][i][e] = new TF1(buf, _f_sum, 0, 5000);
				wf[r][i][e] = new ROOT::Math::WrappedMultiTF1(*func[r][i][e], 1);
			}
		}
	}
	double results[3][npars_max][2];
	
	ROOT::Fit::DataOptions opt;
	opt.fUseEmpty = true;
	for(int site=0;site<3;++site){
		do_fit(site, opt, results);
	}
	printResults(results);
}

void do_fit(int site, ROOT::Fit::DataOptions &opt, double results[3][npars_max][2]){

	char buf[255];

	char *_suffix[3] = { "", "_tag", "_atag" };

	ROOT::Fit::BinData *data[n_range][3][e_slices];
	ROOT::Fit::DataRange rangeD(fitMin, fitMax);
	for(int r=0;r<n_range;++r){
		for(int s=0;s<3;++s){
			for(int e=0;e<e_slices;++e){
				sprintf(buf,"%s/EH%d_dtlSH_%d%s_%d.root", hists_prefix, site+1, r, _suffix[s], e);
				cout << "READING FILE:" << buf << endl;
				TFile *f = new TFile(buf, "READ");
				h[r][s][e] = (TH1*) f->Get("h");
				data[r][s][e] = new ROOT::Fit::BinData(opt, rangeD);
				ROOT::Fit::FillData(*data[r][s][e], h[r][s][e]);
				fPNLL[r][s][e] = new ROOT::Fit::PoissonLikelihoodFCN<ROOT::Math::IBaseFunctionMultiDim>(*data[r][s][e], *wf[r][s][e]);	
			}
			if(!use_tagging) break;
		}
	}

	int npars;
	if(use_tagging)
		npars = 5 * n_range + 1 + 3 * e_slices + 2;
	else
	//	npars = 3 * n_range + 3 * ( e_slices - 1 ) + 2;
		npars = 3 * n_range + 3 * ( e_slices - 1 ) + 1 + 1;

	ROOT::Math::Functor f_tmp(likelihood, npars);
	ROOT::Minuit2::Minuit2Minimizer *minim = ROOT::Math::Factory::CreateMinimizer("Minuit2");
	minim->SetFunction(f_tmp);
	initialize_minimizer(site, minim);
	
	minim->Hesse();
	minim->Minimize();
	minim->Hesse();

	cout << "COV STATUS: " << minim->CovMatrixStatus() << endl;
	minim->PrintResults();

/*	
	double minos_errs[npars][2];	
	int minos_status[npars];
	for(int i=0;i<npars;++i){
		double low, up;
		minim->GetMinosError(i, low, up);
		minos_errs[i][0] = low;
		minos_errs[i][1] = up;	
		minos_status[i] = minim->Status();	
	}
	minim->Hesse();
	minim->PrintResults();
	for(int i=0;i<npars;++i)
		cout << par_names[i] << " " << minos_errs[i][0] << " " << minos_errs[i][1] << " " << minos_status[i] << endl;
	plotHists(site, minim->X());
*/
	for(int i=0;i<npars;++i){
		results[site][i][0]	= minim->X()[i];
		results[site][i][1] = minim->Errors()[i];
	}

	plotHists(site, minim->X());
	plotEnergy(site, minim->X(), minim->Errors());

}

void plotEnergy(int site, double *par, double *par_err){
	
	int offset=0;
	for(;offset<npars_max;++offset)
		if(strstr(par_names[offset],"eps_e")!=NULL)	break;
	cout << offset << " " << par_names[offset] << endl;
	/*	
	double _tmp[3][e_slices][2];
	double _tmp_sum[3] = {1, 1, 1};
	for(int i=0;i<3;++i){
		_tmp[i][0][0] = 1;	
		_tmp[i][0][1] = 0;
		for(int j=0;j<e_slices;++j){
			_tmp[i][j][0] = par[offset + 3 * j + i];
			_tum_sum[i] += _tmp[i][j][0];
		}
	}
	*/
	double x[e_slices], x_err[e_slices], y[3][e_slices], y_err[3][e_slices], y_sum[3];

	for(int i=0;i<3;++i)
		y_sum[i] = 0;

	for(int i=0;i<e_slices;++i){
		x[i] = 0.7 + (12-0.7)/e_slices * (i+0.5);
		x_err[i] = (12-0.7) / ( 2 * e_slices );
		if(i==0){
			for(int j=0;j<3;++j){
				y[j][i] = 1;
				y_err[j][i] = 0;	
			}
		}else{
			y[0][i] = par[offset + 3*(i-1)];
			y_err[0][i] = par_err[offset + 3*(i-1)];
			y[1][i] = par[offset + 3*(i-1)+1];
			y_err[1][i] = par_err[offset + 3*(i-1)+1];
			y[2][i] = par[offset + 3*(i-1)+2];
			y_err[2][i] = par_err[offset + 3*(i-1)+2];
		}
		for(int j=0;j<3;++j)
			y_sum[j] += y[j][i];
	}
	for(int i=0;i<3;++i)
		for(int j=0;j<e_slices;++j){
			y[i][j] /= y_sum[i];
			y_err[i][j] /= y_sum[i];
		}
	TMultiGraph *mg = new TMultiGraph();
	TGraph *gr[3];
	for(int i=0;i<3;++i){
		gr[i] = new TGraphErrors(e_slices, x, y[i], x_err, 0);
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
