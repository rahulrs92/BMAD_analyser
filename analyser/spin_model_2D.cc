#include<analyser.hh>

//void lv2::spin_model(string path,bool out)
tuple<matrix<1,3>,matrix<3,1>,matrix<3,3>,matrix<3,1>> analyser::spin_model_2D(bool out)
{
	fstream f;
	string path1 = path + "/eS";
	auto De = proton::line(path1 + "/re",81);
	auto Dc = proton::line(path1 + "/rc",81);

	auto e = lv1::e(path1 + "/re",De);
	auto c = lv1::c(path1 + "/rc",Dc,e);
	auto pt = lv1::pt(path1 + "/rc",Dc,out);

	f.open((path1 + "/eta").data(),fstream::out | fstream::trunc);
	for(int i=0; i<81; i++)
	{
		int N = 21; string path2 = path1 + "/r" + to_string(i+1);
		auto D = proton::line(path2,N,out);
		auto eta = lv1::eta(path2,D,pt,out);
		auto sig = lv1::sigma(path2,D,e,c);
		Data<string> info(path2 + "/c1/info.dat");
		f<<D.col("#xi_{x}")[0]<<"\t"<<D.col("#xi_{y}")[0]<<"\t"
			<<D.col("#epsilon_{x}")[0]<<"\t"<<D.col("#epsilon_{y}")[0]<<"\t"
			<<eta[0]<<"\t"<<eta[1]<<"\t"<<eta[2]<<"\t"
			<<sig[0]<<"\t"<<sig[1]<<"\t"<<sig[2]<<endl;
	}
	f.close();
	f.open((path + "/eta_info").data(),fstream::out | fstream::trunc);
	f<<"#xi_{x}\n#xi_{y}\n#epsilon_{x}\n#epsilon_{y}\n"
		<<"T\n#eta_{0}\n#eta_{1}\nS\n#sigma_{0}\n#sigma_{1}"<<endl;
	f.close();

	Data<double> D(path1 + "/eta");
	D.label(1) = "#xi_{x}";
	D.label(2) = "#xi_{y}";
	D.label(3) = "#epsilon_{x}";
	D.label(4) = "#epsilon_{y}";
	D.label(5) = "T";
	D.label(6) = "#eta_{0}";
	D.label(7) = "#eta_{1}";
	D.label(8) = "S";
	D.label(9) = "#sigma_{0}";
	D.label(10) = "#sigma_{1}";

	auto eta1 = lv2::eta1(path1,D);
	auto sigma1 = lv2::sigma1(path1,D);
	double eta0 = D.col("#eta_{0}").sum() / D.order()[0];
	double sigma0 = D.col("#sigma_{0}").sum() / D.order()[0];
	/*
	matrix<3,3> M1 = {c[0] - (sigma0/eta0)*pt[0],0,e[0],
									 0,c[1] - (sigma0/eta0)*pt[1],e[1],
									 sigma1[0] - (sigma0/eta0)*eta1[0],
									 sigma1[1] - (sigma0/eta0)*eta1[1],
									 sigma1[2] - (sigma0/eta0)*eta1[2]};
	*/
	matrix<3,3> M1 = {c[0] - (sigma0/eta0)*pt[0],0,0,
									 0,c[1] - (sigma0/eta0)*pt[1],0,
									 sigma1[0] - (sigma0/eta0)*eta1[0],
									 sigma1[1] - (sigma0/eta0)*eta1[1],
									 - (sigma0/eta0)*eta1[2]};
	M = {c[0] - (sigma0/eta0)*pt[0],0,0,
			 0,c[1] - (sigma0/eta0)*pt[1],0,
			 sigma1[0], sigma1[1], -(sigma0/eta0)};

	cout<<"M = "<<M<<endl;

	auto Dt = proton_RF::line(path1 + "/rt",81,0);

	double ex = Dt.col("#epsilon_{x}").sum() / Dt.order()[0];
	double ey = Dt.col("#epsilon_{y}").sum() / Dt.order()[0];
	double da = Dt.col("#delta_{a}").sum() / Dt.order()[0];

	ep = {ex,ey,pow(da,2)/2};
	a = {e[0],e[1],sigma1[2]};
	xi = -(M^-1)*a;
	cout<<"expected optimized point:\n"<<xi<<endl;

	cout<<"ep = "<<ep<<endl;
	cout<<"ep*a = "<<ep*a<<endl;
	matrix<1,3> fcr = ep*M1;
	//cout<<"ep*M = "<<fcr<<endl;
	fcr[2] += ep*a;
	cout<<"spin tune as a function of chromaticity: "<<fcr<<endl;

	cout<<"Properties at optimized point: "<<endl;
	cout<<"second-order phase slip: "<<eta1[0]*xi[0] + eta1[1]*xi[1] + eta1[2]<<endl;
	cout<<"xi = "<<xi<<endl;
	cout<<"M*xi = "<<M*xi<<endl;
	cout<<"ep*M*xi = "<<ep*M*xi<<endl;
	cout<<"ep*a + ep*M*xi = "<<ep*a + ep*M*xi<<endl;

	f.open((path1 + "/matrix").data(),fstream::out | fstream::trunc);
	f<<"M = "<<M<<endl;
	f<<"\nep*M = "<<ep*M<<endl;
	f<<"\nM*xi = "<<M*xi<<endl;
	f<<"\nep*M*xi = "<<ep*M*xi<<endl;
	f<<"ep*a + ep*M*xi = "<<ep*a + ep*M*xi<<endl;
	f<<"spin tune as a function of chromaticity: "<<fcr<<endl;
	f.close();



	return make_tuple(ep,a,M,xi);
	/*
	TApplication app("app", &argc, argv);
	//
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

	cout << "done!" << endl;
	app.Run();
	*/
}


