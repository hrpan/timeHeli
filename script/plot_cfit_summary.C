void plot_cfit_summary(){
	TFile *f = new TFile("./cfit_summary.root");
	TTree *tr = f->Get("tr");
	char *prefix = "./plots/summary/";
	char buf[255], filename[255];
	for(int i=0;i<3;++i){
		for(int j=0;j<3;++j){
			sprintf(buf,"site==%d&&range==%d&&fixB12==0&&eps_pull==1", i+1, j);
		
			tr->Draw("n_lihe", buf);
			sprintf(filename,"%sEH%d_%d_n_lihe.png", prefix, i+1, j);
			c1->SaveAs(filename);

			tr->Draw("n_lihe:fitMin", buf, "colz");
			sprintf(filename,"%sEH%d_%d_n_lihe_fitMin.png", prefix, i+1, j);
			c1->SaveAs(filename);	

			tr->Draw("n_lihe_err:fitMin", buf, "colz");
			sprintf(filename,"%sEH%d_%d_n_lihe_err_fitMin.png", prefix, i+1, j);
			c1->SaveAs(filename);
	
			tr->Draw("n_bo:fitMin", buf, "colz");
			sprintf(filename,"%sEH%d_%d_n_bo_fitMin.png", prefix, i+1, j);
			c1->SaveAs(filename);

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

			tr->Draw("r_lihe:fitMin", buf, "colz");
			sprintf(filename,"%sEH%d_%d_r_lihe_fitMin.png", prefix, i+1, j);
			c1->SaveAs(filename);

			tr->Draw("r_lihe:fixHe8+2*fixB12+4*fix_rmu+8*eps_pull", buf, "colz");
			sprintf(filename,"%sEH%d_%d_r_lihe_fixs.png", prefix, i+1, j);
			c1->SaveAs(filename);

			tr->Draw("eps_lihe:fitMin", buf, "colz");
			sprintf(filename,"%sEH%d_%d_eps_lihe_fitMin.png", prefix, i+1, j);
			c1->SaveAs(filename);

			tr->Draw("eps_lihe:fixHe8+2*fixB12+4*fix_rmu+8*eps_pull", buf, "colz");
			sprintf(filename,"%sEH%d_%d_eps_lihe_fixs.png", prefix, i+1, j);
			c1->SaveAs(filename);

		}
	}
}
