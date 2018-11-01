#include "util.C"
#include "range.h"

const char *hists_prefix = "./hists";

const double tau_li9 = 256.366;
const double tau_li9_err = 0.866;

const double tau_he8 = 171.68;

const double tau_b12 = 29.142;
const double tau_b12_err = 0.0288;

double fitMin = 1.5;
double fitMax = 5000;
const int eps_slices = 3;
const int npars = 4 * n_range + eps_slices + 1 + 2 ;
/*
 * r_mu_tag, r_mu_atag : 2 * n_range
 * n_lihe, n_bo : 2 * n_range
 * eps_bo, eps_lihe : 2 * n_range 
 * ratio_lihe : 1 
 * n_dc : 1
 */
char par_names[npars][255]; 

double eps_pulls[3][eps_slices][2] = {
	{
		{},
		{},
		{}
	},
	{
		{},
		{},
		{}
	},
	{
		{ 0.99, 0.53 },
		{ 0.73, 0.12},
		{ 0.97, 0}
	}
};

const bool latex_print[npars] = {
	false, false, false,
	true, true, true,
	false, false,
	true, true, false
};

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

void printLatexTable(double results[3][n_range][npars][2]){
	for(int i=0;i<3;++i){
		for(int j=0;j<n_range;++j){
			if(j==n_range-1)
				printf("%4d & [%10.1e,%10s] ",i+1,range[j],"inf");
			else
				printf("%4d & [%10.1e,%10.1e] ",i+1,range[j],range[j+1]);
			for(int k=0;k<npars;++k)
				if(latex_print[k])
					printf("& $%10.2f\\pm%10.2f$ ",fabs(results[i][j][k][0]),results[i][j][k][1]);
			
			printf("\\\\ \n");
		}
		printf("\\hline \n");
	}

}

ROOT::Math::IMultiGenFunction *fPNLL[n_range][3];
ROOT::Math::WrappedMultiTF1 *wf[n_range][3];

TF1 *func[n_range][3];
TH1 *h[n_range][3];

void parTrans(const double *par, double _par[n_range][3][npars_single]){

	double scale = h[0][0]->GetBinWidth(1);

	double r_mu_tag[n_range];
	double r_mu_atag[n_range];
	double n_lihe[n_range];
	double n_bo[n_range];
	double eps_lihe[n_range];
	double _eps_lihe[eps_slices];
	//double eps_bo[n_range];
	double eps_bo;
	double n_dc;
	double r_lihe;

	int p_idx = 0;
	for(int i = 0; i < n_range; ++i){
		r_mu_tag[i] = par[p_idx++];
		r_mu_atag[i] = par[p_idx++];
		n_lihe[i] = par[p_idx++] * scale;
		//eps_lihe[i] = par[p_idx++];
		n_bo[i] = par[p_idx++] * scale;
	//	eps_bo[i] = par[p_idx++];
	}
	for(int i=0;i<eps_slices;++i)
		_eps_lihe[i] = par[p_idx++];
	eps_lihe[0] = _eps_lihe[0];
	eps_lihe[1] = _eps_lihe[0];
	eps_lihe[2] = _eps_lihe[0];
	eps_lihe[3] = _eps_lihe[1];
	eps_lihe[4] = _eps_lihe[2];
	eps_bo[i] = par[p_idx++];
	n_dc = par[p_idx++] * scale;
	r_lihe = par[p_idx++];

	for(int i=0;i<n_range;++i){
		double n_dc_fake = 0;
		for(int j=0;j<n_range;++j){
			if(i==j) continue;
			n_dc_fake += (n_lihe[j] + n_bo[j]);
		}
		for(int j=0;j<3;++j){
			if(j==0){
				_par[i][j][0] = r_mu_tag[i] + r_mu_atag[i];
				_par[i][j][1] = n_dc + n_dc_fake;
				_par[i][j][2] = n_lihe[i];
				_par[i][j][6] = n_bo[i];
			}else if(j==1){
				_par[i][j][0] = r_mu_tag[i];
				_par[i][j][1] = n_dc + (n_dc_fake + (( 1 - eps_lihe[i] ) * n_lihe[i] + ( 1 - eps_bo[i] ) * n_bo[i]));
				_par[i][j][2] = eps_lihe[i] * n_lihe[i];
				_par[i][j][6] = eps_bo[i] * n_bo[i];
			}else{
				_par[i][j][0] = r_mu_atag[i];
				_par[i][j][1] = n_dc + (n_dc_fake + (eps_lihe[i] * n_lihe[i] + eps_bo[i] * n_bo[i]));
				_par[i][j][2] = ( 1 - eps_lihe[i] ) * n_lihe[i];
				_par[i][j][6] = ( 1 - eps_bo[i] ) * n_bo[i];
			}
			_par[i][j][3] = tau_li9;
			_par[i][j][4] = tau_he8;
			_par[i][j][5] = r_lihe;
			_par[i][j][7] = tau_b12;
/*
			cout << i << " " << j << " ";
			for(int k=0;k<npars_single;++k)
				printf("%6e ", _par[i][j][k]);
			cout << endl;
*/
		}

	}

}

double likelihood(const double *par){
	double _par[n_range][3][npars_single];

	parTrans(par, _par);
	
	double result = 0;
	vector<double> _nlls(n_range * 3);
	for(int r = 0; r < n_range; ++r){
		for(int s = 0; s < 3; ++s){
			double _tmp = (*fPNLL[r][s])(_par[r][s]);
			_nlls[ r * 3 + s ] = _tmp;
			result += _tmp;
		}
	}
	double safesum = neumaierSum(&_nlls[0], n_range * 3);
	//cout << result << " " << safesum << " " << result - safesum << endl;
	return safesum;
}


void fillPars(double *par, TF1 *f[n_range][3]){

	double _par[n_range][3][npars_single];

	parTrans(par, _par);

	for(int i=0;i<n_range;++i)
		for(int j=0;j<3;++j)
			f[i][j]->SetParameters(_par[i][j]);

}

void plotHists(int site, double *par){
	fillPars(par, func);
	char *dir = "./plots/cfits";
	char buf[255];
	double _par[n_range][3][npars_single];
	parTrans(par, _par);
	for(int r=0;r<n_range;++r)
		for(int i=0;i<3;++i){
			TF1 ibd("ibd","[1] * [0] * exp(-[0] * x)",0, 5000);
			ibd.SetParameters(_par[r][i]);
			ibd.SetLineColor(kGreen);
			sprintf(buf, "%s/cfit_%d_%d_%d.png", dir, site, r, i);
			h[r][i]->Draw("E1");
			func[r][i]->Draw("same");
			ibd.Draw("same");
			gPad->SetLogx();
			c1->SaveAs(buf);
		}
	
}

void initialize_minimizer(int site, ROOT::Math::Minimizer *minim){
	minim->SetErrorDef(0.5);
	minim->SetTolerance(0.01);
	minim->SetPrintLevel(1);
	//minim->SetPrecision(1e-16);
	minim->SetStrategy(2);
	minim->SetMaxFunctionCalls(1000000);

	double rmu_init;
	double n_lihe_init = 1e3;
	double n_bo_init = 1e1;
	double eps_init = 0.7;
	double r_lihe_init = 0.8;

	double step_ratio = 1e-2;
	char buf[255];
	int p_idx = 0;
	for(int i=0;i<n_range;++i){
		h[i][1]->Fit("expo","LQN0","", 300, 5000);
		rmu_init = fabs(expo->GetParameter(1));
		sprintf(buf, "rmu_tag_%d",i);
		sprintf(par_names[p_idx], "%s", buf);
		minim->SetVariable(p_idx++, buf, rmu_init, rmu_init * step_ratio);

		h[i][2]->Fit("expo","LQN0","", 300, 5000);
		rmu_init = fabs(expo->GetParameter(1));
		sprintf(buf, "rmu_atag_%d",i);
		sprintf(par_names[p_idx], "%s", buf);
		minim->SetVariable(p_idx++, buf, rmu_init, rmu_init * step_ratio);

		sprintf(buf, "n_lihe_%d",i);
		sprintf(par_names[p_idx], "%s", buf);
		minim->SetLowerLimitedVariable(p_idx++, buf, n_lihe_init, n_lihe_init * step_ratio , 0);
		//minim->SetVariable(p_idx++, buf, n_lihe_init, n_lihe_init * step_ratio);
	/*
		sprintf(buf, "eps_lihe_%d",i);
		sprintf(par_names[p_idx], "%s", buf);
		//minim->SetLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0, 1);
		minim->SetVariable(p_idx++, buf, eps_init, eps_init * step_ratio);
	*/
		sprintf(buf, "n_bo_%d",i);
		sprintf(par_names[p_idx], "%s", buf);
		minim->SetLowerLimitedVariable(p_idx++, buf, n_bo_init, n_bo_init * step_ratio , 0);
		//minim->SetVariable(p_idx++, buf, n_lihe_init, n_lihe_init * step_ratio);

//		sprintf(buf, "eps_bo_%d",i);
//		sprintf(par_names[p_idx], "%s", buf);
//		minim->SetLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0, 1);
	}
	for(int i=0;i<eps_slices;++i){
		sprintf(buf, "eps_lihe_%d",i);
		sprintf(par_names[p_idx], "%s", buf);
		minim->SetLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0, 1);
	}
	sprintf(buf, "eps_bo",i);
	sprintf(par_names[p_idx], "%s", buf);
	minim->SetLimitedVariable(p_idx++, buf, eps_init, eps_init * step_ratio, 0, 1);


	double n_dc_init = h[0][0]->GetEntries();
	sprintf(par_names[p_idx], "n_dc");
	minim->SetVariable(p_idx++, "n_dc", n_dc_init, n_dc_init * step_ratio);

	sprintf(par_names[p_idx], "r_lihe");
	minim->SetLimitedVariable(p_idx++, "r_lihe", r_lihe_init, r_lihe_init * step_ratio, 0, 1);
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
			sprintf(buf,"func%d_%d",r,i);
			func[r][i] = new TF1(buf, _f_sum, 0, 5000);
			wf[r][i] = new ROOT::Math::WrappedMultiTF1(*func[r][i], 1);
		}
	}
	double results[3][npars][2];
	
	ROOT::Fit::DataOptions opt;
	opt.fUseEmpty = true;
	for(int site=0;site<3;++site){
		do_fit(site, opt, results);
	}
	printResults(results);
}

void do_fit(int site, ROOT::Fit::DataOptions &opt, double results[3][npars][2]){

	char buf[255];

	char *_suffix[3] = { "", "_tag", "_atag" };

	ROOT::Fit::BinData *data[n_range][3];
	ROOT::Fit::DataRange rangeD(fitMin, fitMax);
	for(int r=0;r<n_range;++r){
		for(int s=0;s<3;++s){
			sprintf(buf,"%s/EH%d_dtlSH_%d%s.root", hists_prefix, site+1, r, _suffix[s]);
			cout << "READING FILE:" << buf << endl;
			TFile *f = new TFile(buf, "READ");
			h[r][s] = (TH1*) f->Get("h");
			data[r][s] = new ROOT::Fit::BinData(opt, rangeD);
			ROOT::Fit::FillData(*data[r][s], h[r][s]);
			fPNLL[r][s] = new ROOT::Fit::PoissonLikelihoodFCN<ROOT::Math::IBaseFunctionMultiDim>(*data[r][s], *wf[r][s]);	
		}
	}
	
	ROOT::Math::Functor f_tmp(likelihood, npars);
	ROOT::Minuit2::Minuit2Minimizer *minim = ROOT::Math::Factory::CreateMinimizer("Minuit2");
	minim->SetFunction(f_tmp);
	initialize_minimizer(site, minim);
	minim->Hesse();
	minim->Minimize();
	minim->Hesse();
	cout << "COV STATUS: " << minim->CovMatrixStatus() << endl;
	minim->PrintResults();
	
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
	
}


