#include <string>

const int nArgs = 8;
const char *arg_names[nArgs] = {
	"fitMin",
	"fitMax",
	"fixB12",
	"fixHe8",
	"bound_eps",
	"fix_lifetime",
	"fix_rmu",
	"eps_pull"
};

int site, range;
double r_mu_tag[2];
double r_mu_atag[2];
double n_dc[2], n_lihe[2], eps_lihe[2], r_lihe[2], n_bo[2], eps_bo[2];

void print_Args(int *args){
	for( size_t i = 0; i < nArgs; ++i)
		cout << arg_names[i] << " " << args[i] << endl;
}

void parse_file_name(const char *fn, int *args){
	int idx = 0;
	for(size_t i = 0; i < strlen(fn); ++i){
		if( fn[i] == '_' ){
			char buf[255];
			memset(buf, 0, 255);
			for(size_t j = 1; i + j < strlen(fn); ++j){
				if( fn[i+j] == '_' ){
					memcpy( buf, &fn[i+1], j-1);
					buf[j] = '\0';
					i += j-1;	
					break;
				}else if( i + j == strlen(fn) - 1 ){
					memcpy( buf, &fn[i+1], j);
					buf[j] = '\0';
					i += j;
				}
			}
			args[idx++] = atoi(buf);
		}
	}
}

void parse_word(string str, double *val){
	int slash = str.find("/");
	val[0] = atof(str.substr(0, slash).c_str());
	val[1] = atof(str.substr(slash+1, -1).c_str());
}

void parse_and_fill_output(const char *fn, TTree *tr){
	const int n_range = 3;
	cout << fn << endl;
	ifstream ifile(fn);
	string line;
	bool ignore = true;
	while(getline(ifile, line))
		if( line.find("range") != std::string::npos ) break;

	for(int i=0;i<3;++i){
		for(int r=0;r<n_range;++r){
			getline(ifile, line);
			vector<string> tmp;
			int start = 0;
			int end = line.find(" ");
			while( end != std::string::npos ){
				if(end-start>0)
					tmp.push_back(line.substr(start, end-start));
				start = end + 1;
				end = line.find(" ", start);
			}		
			tmp.push_back(line.substr(start, end-start));	
			site = atoi(tmp[0].c_str());
			range = r;
			parse_word(tmp[3], r_mu_tag);
			parse_word(tmp[4], r_mu_atag);
			parse_word(tmp[5], n_dc);
			parse_word(tmp[6], n_lihe);
			parse_word(tmp[7], eps_lihe);
			parse_word(tmp[8], r_lihe);
			parse_word(tmp[9], n_bo);
			parse_word(tmp[10], eps_bo);
			tr->Fill();	
		}		
		getline(ifile, line);
	}
		
	

}

void parse_cfit_outputs(){
	const char *output_dir = "./cfit_outputs/";
	TSystemDirectory *dir = new TSystemDirectory(output_dir, output_dir);
	TList *files = dir->GetListOfFiles();

	int bArgs[nArgs];
	TTree *tr = new TTree("tr","tr");
	for(size_t i = 0; i < nArgs; ++i)
		tr->Branch(arg_names[i], &bArgs[i]);
	tr->Branch("site", &site);
	tr->Branch("range", &range);
	tr->Branch("r_mu_tag", &r_mu_tag[0]);
	tr->Branch("r_mu_tag_err", &r_mu_tag[1]);
	tr->Branch("r_mu_atag", &r_mu_atag[0]);
	tr->Branch("r_mu_atag_err", &r_mu_atag[1]);
	tr->Branch("n_dc", &n_dc[0]);
	tr->Branch("n_dc_err", &n_dc[1]);
	tr->Branch("n_lihe", &n_lihe[0]);
	tr->Branch("n_lihe_err", &n_lihe[1]);
	tr->Branch("eps_lihe", &eps_lihe[0]);
	tr->Branch("eps_lihe_err", &eps_lihe[1]);
	tr->Branch("r_lihe", &r_lihe[0]);
	tr->Branch("r_lihe_err", &r_lihe[1]);
	tr->Branch("n_bo", &n_bo[0]);
	tr->Branch("n_bo_err", &n_bo[1]);
	tr->Branch("eps_bo", &eps_bo[0]);
	tr->Branch("eps_bo_err", &eps_bo[1]);

	for(size_t i = 2; i < files->GetSize(); ++i){
		cout << files->At(i)->GetName() << endl;		
		char *fn = files->At(i)->GetName();
		
		parse_file_name(fn, bArgs);
		
		char buf[255];
		sprintf(buf, "%s%s", output_dir, fn);
		parse_and_fill_output(buf, tr);
	}	
	tr->SaveAs("cfit_summary.root");
}
