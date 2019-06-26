double eff[3][5][2] = {
	{
		{  1.17,      7.16 },
		{  0.11,      0.08 },
		{  0.81,      0.19 },
		{  0.98,      0.17 },
		{  0.90,      0.03 }
	},
	{
		{  0.49,      0.87 },
		{ -0.16,      0.14 },
		{  0.56,      0.14 },
		{  0.64,      0.09 },
		{  0.95,      0.04 }
	},
	{
		{  1.92,      1.20 },
		{  0.22,      0.08 },
		{  0.51,      0.09 },
		{  1.01,      0.16 },
		{  0.95,      0.03 }
	}

};

double range[5] = {
	0.03,1.00,2.00,
	3.00,4.00
};

double upperbound = 8.0;
int start = 0;
double axis_ext = 0.3;

void fitWeibull(){

	gStyle->SetOptFit(1);

	TF1 *func = new TF1("func","1-exp(-(x/[0])**[1])",0,10);

	func->SetParameter(0,3);
	func->SetParameter(1,1);
	func->SetParLimits(1,0,5);
	//func->FixParameter(1,1);

	char buf[255];
	for(int s=0;s<3;++s){
		TGraphErrors *g = new TGraphErrors(4);
		TAxis *axis = g->GetXaxis();
		axis->SetLimits(0,upperbound+axis_ext);
		for(int i=start;i<5;++i){
			if(eff[s][i][0] < 0) continue;
			if(i<4){
				g->SetPoint(i,(range[i]+range[i+1])/2,eff[s][i][0]);
				g->SetPointError(i,(range[i+1]-range[i])/2,eff[s][i][1]);
			}else{
				g->SetPoint(i,(range[i]+upperbound)/2,eff[s][i][0]);
				g->SetPointError(i,(upperbound-range[i])/2,eff[s][i][1]);
			}
		}
		g->GetHistogram()->SetMinimum(-axis_ext);
		g->GetHistogram()->SetMaximum(1 + axis_ext);
		g->Fit(func);
		g->Draw("A P");
		func->Draw("same");
		sprintf(buf,"./plots/eff_fit_%d.png",s);
		c1->SaveAs(buf);
		g->Delete();
	}
}
