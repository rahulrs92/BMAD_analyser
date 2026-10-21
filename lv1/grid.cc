#include<lv1.hh>

TTree* lv1::grid_RF(double iQx, double iQy, TFile* ifile, bool pdf)
{
	auto grid_tree = ifile->Get<TTree>("grid");
	if(grid_tree) return grid_tree;
	else grid_tree = new TTree("grid","grid");

	double Cx, Cy, epx, epy, dta, oms, eta0, nus;
	grid_tree->Branch("#xi_{x}",&Cx);
	grid_tree->Branch("#xi_{y}",&Cy);
	grid_tree->Branch("#epsilon_{x}",&epx);
	grid_tree->Branch("#epsilon_{y}",&epy);
	grid_tree->Branch("#delta_{a}",&dta);
	grid_tree->Branch("#omega_{s}",&oms);
	grid_tree->Branch("#eta_{0}",&eta0);
	grid_tree->Branch("#nu_{s}",&nus);

	auto itree = ifile->Get<TTree>("RF_on");
	auto xitree = ifile->Get<TTree>("xi");
	int c, r;
	double cx, cy;
	vector<double> *n, *x, *y, *z, *px, *py, *pz, *sx, *sy, *sz;
	itree->SetBranchAddress("r",&r);
	itree->SetBranchAddress("c",&c);
	itree->SetBranchAddress("l",&n);
	itree->SetBranchAddress("x",&x);
	itree->SetBranchAddress("y",&y);
	itree->SetBranchAddress("z",&z);
	itree->SetBranchAddress("px",&px);
	itree->SetBranchAddress("py",&py);
	itree->SetBranchAddress("pz",&pz);
	itree->SetBranchAddress("sx",&sx);
	itree->SetBranchAddress("sy",&sy);
	itree->SetBranchAddress("sz",&sz);
	xitree->SetBranchAddress("Cx",&cx);
	xitree->SetBranchAddress("Cy",&cy);

	int Nc = itree->GetEntries("r==1");
	int Nr = itree->GetEntries("c==1");
	for(size_t i=0; i<Nr; i++)
	{
		xitree->GetEntry(i);
		for(size_t j=0; j<Nc; j++)
		{
			itree->GetEntry(i*Nc + j);
			proton_RF P(iQx, iQy, cx, cy, *n, *x, *y, *z, *px, *py, *pz, *sx, *sy, *sz);
			auto CSx = P.Courant_Snyder(0,pdf);
			auto CSy = P.Courant_Snyder(1,pdf);
			auto syn = P.longitudinal(pdf);
			auto eta = P.phase_slip(pdf);
			nus = P.spin_tune(pdf);

			Cx = cx; Cy = cy;
			epx = CSx[0]; epy = CSy[0];
			dta = syn[3]; oms = syn[4]; eta0 = eta[1];

			grid_tree->Fill();
		}
	}
	grid_tree->Write();
	return grid_tree;
}

TTree* lv1::grid(double iQx, double iQy, TFile* ifile, bool pdf)
{

	auto grid_tree = ifile->Get<TTree>("grid");
	if(grid_tree) return grid_tree;
	else grid_tree = new TTree("grid","grid");
	double Cx, Cy, epx, epy, dt0, ps, nus;
	grid_tree->Branch("#xi_{x}",&Cx);
	grid_tree->Branch("#xi_{y}",&Cy);
	grid_tree->Branch("#epsilon_{x}",&epx);
	grid_tree->Branch("#epsilon_{y}",&epy);
	grid_tree->Branch("#delta_{0}",&dt0);
	grid_tree->Branch("T",&ps);
	grid_tree->Branch("#nu_{s}",&nus);

	//TMultiGraph *x_px, *y_py;
	//grid_tree->Branch("plot_x-px",&x_px);
	//grid_tree->Branch("plot_y-py",&y_py);

	auto itree = ifile->Get<TTree>("RF_off");
	auto xitree = ifile->Get<TTree>("xi");
	int c, r;
	double cx, cy;
	vector<vector<double>> C(10,vector<double>());
	/*
	itree->SetBranchAddress("l",&C[0]);
	itree->SetBranchAddress("x",&C[1]);
	itree->SetBranchAddress("y",&C[2]);
	itree->SetBranchAddress("z",&C[3]);
	itree->SetBranchAddress("px",&C[4]);
	itree->SetBranchAddress("py",&C[5]);
	itree->SetBranchAddress("pz",&C[6]);
	itree->SetBranchAddress("sx",&C[7]);
	itree->SetBranchAddress("sy",&C[8]);
	itree->SetBranchAddress("sz",&C[9]);
	*/
	vector<double> *n = &C[0], *x = &C[1], *y = &C[2], *z = &C[3],
									*px = &C[4], *py = &C[5], *pz = &C[6],
									*sx = &C[7], *sy = &C[8], *sz = &C[9];
	itree->SetBranchAddress("r",&r);
	itree->SetBranchAddress("c",&c);
	itree->SetBranchAddress("l",&n);
	itree->SetBranchAddress("x",&x);
	itree->SetBranchAddress("y",&y);
	itree->SetBranchAddress("z",&z);
	itree->SetBranchAddress("px",&px);
	itree->SetBranchAddress("py",&py);
	itree->SetBranchAddress("pz",&pz);
	itree->SetBranchAddress("sx",&sx);
	itree->SetBranchAddress("sy",&sy);
	itree->SetBranchAddress("sz",&sz);
	xitree->SetBranchAddress("Cx",&cx);
	xitree->SetBranchAddress("Cy",&cy);

	int Nc = itree->GetEntries("r==1");
	int Nr = itree->GetEntries("c==1");
	//Nr = 1; Nc = 1;
	//for(size_t i=9994; i<9995; i++)
	for(size_t i=0; i<Nr; i++)
	{
		xitree->GetEntry(i);
		cout<<"Analysing Simulation "<<i;//<<": spin tune: ";//
		for(size_t j=0; j<Nc; j++)
		{
			itree->GetEntry(i*Nc + j);
			proton P(iQx, iQy, cx, cy, C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8], C[9]);
			//P.setPath("/home/Rahul/JEDI/PTEDM/1.0/cluster/cluster_data/eS/j" + to_string(j));

			auto CSx = P.Courant_Snyder(0,pdf);
			auto CSy = P.Courant_Snyder(1,pdf);
			auto eta = P.phase_slip(pdf);
			nus = P.spin_tune(pdf);
			//x_px = &P.mg["CSx"];
			//y_py = &P.mg["CSy"];

			Cx = cx; Cy = cy;
			epx = CSx[0]; epy = CSy[0];
			dt0 = C[6][0]; // j*0.0002;//
			ps = eta[1]+eta[0];
			//cout<<" "<<nus;
			/*
			nus = 0;
			Cx = cx; Cy = cy;
			epx = 0; epy = 0;
			dt0 = C[6][0]; ps = 0;//eta[1]+eta[0];
			*/

			grid_tree->Fill();
			//delete P;
		}
		cout<<endl;
	}
	grid_tree->Write();
	return grid_tree;
}

void lv1::grid(acquisition& ia,string path,bool out)
{
	fstream f;
	string gridfile = path + (ia.zero() ? "/grid0" : "/grid");
	f.open(gridfile,fstream::in);
	if(f.is_open())
	{
		f.close();
		return;
	}
	f.close();
	f.open(path + "/grid_info",fstream::out | fstream::trunc);
	f<<"t_{f}\nt_{d}\nt_{s}\n#xi_{x}\n#xi_{y}\n";
	f<<"#epsilon_{x}\n#epsilon_{y}\n";
	f<<"#epsilon_{x0}\n#epsilon_{y0}\n";
	f<<"#delta_{0}\n#nu_{s}\nT"<<endl;
	f.close();
	int nth = 5, q = 0;
	int Nr = ia.Nr, Nc = ia.Nc;
	double Qx = ia.Qx, Qy = ia.Qy, L = ia.L;
	vector<vector<double>> xi(5);
	ia.xi(xi[0],xi[1],xi[2],xi[3],xi[4]);
	vector<array<double,7>> D;
	vector<acquisition> A;
	for(int k=0; k<nth; k++) A.push_back(acquisition({Nc,Nr,Nc*Nr},{L,Qx,Qy},k));
	for(int k=0; k<nth; k++) A[k].set_zero_flag(ia.zero());
	vector<thread> vth;
	f.open(gridfile.data(),fstream::out | fstream::trunc);
	int crap = 0;
	for(int i=0; i<Nr; i++)
	{
		for(int j=0; j<Nc; j++)
		{
			// This part is to override the inexplicable missing data in the first thread set
			if(crap == 1)
			{
				i=0; j=0; crap++;
				f.close(); f.open(gridfile.data(),fstream::out | fstream::trunc);
			}
			//

			D.push_back({0,0,0,0,0});
			array<double,7> &ar = D[q];
			double &cx = xi[3][i], &cy = xi[4][i];
			acquisition &aq = A[q];
			vth.push_back
			(
				thread
				(
					[&ar,&aq,i,j,cx,cy,out,this]()
					{
						ar = point(aq,i+1,j+1,cx,cy,out);
						//cout<<"Contents of D: ";
						//for(int l=0; l<4; l++) cout<<ar[l]<<" ";
						//cout<<ar[4]<<endl;
					}
				)
			);
			//vth.push_back(&tt);
			q++;
			if(q == nth || (i == (Nr-1) && j == (Nc-1)))
			{
				crap++;
				cout<<"Size of D: "<<D.size()<<endl;
				cout<<"Size of vth: "<<vth.size()<<endl;
				for(int k=0; k<vth.size(); k++) vth[k].join();
				for(int k=0; k<D.size(); k++)
				{
					for(int l=0; l<5; l++) f<<xi[l][i]<<"\t";
					for(int l=0; l<6; l++) f<<D[k][l]<<"\t";
					f<<D[k][6]<<endl;
				}
				q = 0;
				D.clear();
				vth.clear();
			}
		}
	}
	f.close();
}

void lv1::grid(string path,bool out)
{
	fstream f;
	string gridfile = path + "/grid";
	f.open(gridfile,fstream::in);
	if(f.is_open())
	{
		f.close();
		return;
	}
	f.close();
	f.open(path + "/grid_info",fstream::out | fstream::trunc);
	f<<"t_{f}\nt_{d}\nt_{s}\n#xi_{x}\n#xi_{y}\n";
	f<<"#epsilon_{x}\n#epsilon_{y}\n";
	f<<"#epsilon_{x0}\n#epsilon_{y0}\n";
	f<<"#delta_{0}\nT\n#nu_{s}"<<endl;
	f.close();
	int NL=0, Nc=0, Nr=0, Ncr=0;
	string s_NL = "cat " + path + "/r1/c1/spin.dat | wc -l";
	string s_Nc = "find " + path + " -type d -name 'c*' | wc -l";
	string s_Nr = "find " + path + " -type d -name 'r*' | wc -l";
	string s_Ncr = "find " + path + " -type f -name 'spin.dat' | wc -l";
	NL = bash::to_int(s_NL);
	Nc = bash::to_int(s_Nc);
	Nr = bash::to_int(s_Nr);
	Ncr = bash::to_int(s_Ncr);
	cout<<"NL = "<<NL<<endl;
	cout<<"Nc = "<<Nc<<endl;
	cout<<"Nr = "<<Nr<<endl;
	cout<<"Ncr = "<<Ncr<<endl;

	if(Ncr == Nc * Nr) cout<<"local repository consistent!"<<endl;
	else cout<<"local repository NOT consistent!"<<endl;
	int nth = 1, q = 0;
	vector<array<double,7>> D;
	vector<thread> vth;
	f.open(gridfile.data(),fstream::out | fstream::trunc);
	int crap = 0;
	for(int i=0; i<Nr; i++)
	{
		for(int j=0; j<Nc; j++)
		{
			// This part is to override the inexplicable missing data in the first thread set
			if(crap == 1)
			{
				i=0; j=0; crap++;
				f.close(); f.open(gridfile.data(),fstream::out | fstream::trunc);
			}
			//

			D.push_back({0,0,0,0,0});
			array<double,7> &ar = D[q];
			array<double,5> xi;
			vth.push_back(thread(
				[&ar,&xi,path,i,j,out,this]()
				{
					ar = point(path,i+1,j+1,out);
					xi = lv1::xi(path,i+1,j+1,out);
					//cout<<"Contents of D: ";
					//for(int l=0; l<4; l++) cout<<ar[l]<<" ";
					//cout<<ar[4]<<endl;
				}
			));
			//vth.push_back(&tt);
			q++;
			if(q == nth || (i == (Nr-1) && j == (Nc-1)))
			{
				crap++;
				cout<<"Size of D: "<<D.size()<<endl;
				cout<<"Size of vth: "<<vth.size()<<endl;
				for(int k=0; k<vth.size(); k++) vth[k].join();
				for(int k=0; k<D.size(); k++)
				{
					for(int l=0; l<5; l++) f<<xi[l]<<"\t";
					for(int l=0; l<6; l++) f<<D[k][l]<<"\t";
					f<<D[k][6]<<endl;
				}
				q = 0;
				D.clear();
				vth.clear();
			}
		}
	}
	f.close();
}

void lv1::grid_xi(string path,bool out)
{
	fstream f;
	string gridfile = path + "/grid";
	f.open(gridfile,fstream::in);
	if(f.is_open())
	{
		f.close();
		return;
	}
	f.close();
	f.open(path + "/grid_info",fstream::out | fstream::trunc);
	f<<"t_{f}\nt_{d}\nt_{s}\n#xi_{x}\n#xi_{y}\n";
	f<<"#epsilon_{x}\n#epsilon_{y}\n";
	f<<"#epsilon_{x0}\n#epsilon_{y0}\n";
	f<<"#delta_{0}\nT\n#nu_{s}"<<endl;
	f.close();
	int NL=0, Nc=0, Nr=0, Ncr=0;
	string s_NL = "cat " + path + "/r1/c1/spin.dat | wc -l";
	string s_Nc = "find " + path + "/r1 -type f -name 'spin.dat' | wc -l";
	string s_Nr = "find " + path + "/r[0-9]*/c1 -type f -name 'spin.dat' | wc -l";
	string s_Ncr = "find " + path + " -type f -name 'spin.dat' | wc -l";
	NL = bash::to_int(s_NL);
	Nc = bash::to_int(s_Nc);
	Nr = bash::to_int(s_Nr);
	Ncr = bash::to_int(s_Ncr);
	cout<<"NL = "<<NL<<endl;
	cout<<"Nc = "<<Nc<<endl;
	cout<<"Nr = "<<Nr<<endl;
	cout<<"Ncr = "<<Ncr<<endl;

	if(Ncr == Nc * Nr) cout<<"local repository consistent!"<<endl;
	else cout<<"local repository NOT consistent!"<<endl;
	Data<double> xi(path + "/xi");
	int nth = 5, q = 0;
	vector<array<double,7>> D;
	vector<thread> vth;
	f.open(gridfile.data(),fstream::out | fstream::trunc);
	int crap = 0;
	for(int i=0; i<Nr; i++)
	{
		for(int j=0; j<Nc; j++)
		{
			// This part is to override the inexplicable missing data in the first thread set
			if(crap == 1)
			{
				i=0; j=0; crap++;
				f.close(); f.open(gridfile.data(),fstream::out | fstream::trunc);
			}
			//

			D.push_back({0,0,0,0,0});
			array<double,7> &ar = D[q];
			vth.push_back(thread(
				[&ar,path,i,j,out,this]()
				{
					ar = point(path,i+1,j+1,out);
					//cout<<"Contents of D: ";
					//for(int l=0; l<4; l++) cout<<ar[l]<<" ";
					//cout<<ar[4]<<endl;
				}
			));
			//vth.push_back(&tt);
			q++;
			if(q == nth || (i == (Nr-1) && j == (Nc-1)))
			{
				crap++;
				cout<<"Size of D: "<<D.size()<<endl;
				cout<<"Size of vth: "<<vth.size()<<endl;
				for(int k=0; k<vth.size(); k++) vth[k].join();
				for(int k=0; k<D.size(); k++)
				{
					for(int l=0; l<5; l++) f<<xi[l+1][i+1]<<"\t";
					for(int l=0; l<6; l++) f<<D[k][l]<<"\t";
					f<<D[k][6]<<endl;
				}
				q = 0;
				D.clear();
				vth.clear();
			}
		}
	}
	f.close();
}

