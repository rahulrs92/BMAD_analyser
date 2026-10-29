#include<thread>

#include "TF1.h"
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"

#include<proton.hh>
#include<lv1.hh>
#include<lv2.hh>
#include<acquisition.hh>

#include<matrix3.hh>
#include<bash_tools.hh>

using namespace std;

string proton_RF::ID = "proton"; // "deuteron"; // 
string proton_RF::mode = "COSY"; // "COSY_mis"; // "PSR"; // "PSR-electric"; // "PSR-magnetic"; // 

string acquisition::source = "PTEDM/source_COSY/sct_test";
//string acquisition::source = "PTEDM/prototype_particle/eMP";
string acquisition::destination = "/Users/utente/JEDI/PTEDM/1.0/source_prototype2/build/data/sct2";

int main_spin_tune_analysis(int argc, char **argv);

int main_beam_sct_analysis(int argc, char **argv);

int main_multiplane_analysis(int argc, char **argv);

int main_multiplane_RFoff(int argc, char **argv);

int main_single_plane_analysis(int argc, char **argv);

int main_travel_time_analysis(int argc, char **argv);

int main_spin_tune_analysis_RFoff(int argc, char **argv);

int main_tracking_chromaticity(int argc, char **argv);

int main_calculate_SCT(int argc, char **argv)
{
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/sct1/c1";
	lv2 az;

	proton_RF p(path);
		
	TApplication app("app", &argc, argv);

	p.SCT();
	p.spin_tune();
	
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	return 0;
}

int main_reference(int argc, char **argv)
{
	proton P("/Users/utente/JEDI/PTEDM/1.0/source_COSY_mis/build/data/Deuteron/reference/c1");
	//auto CSx = P.Courant_Snyder(0);
	//auto CSy = P.Courant_Snyder(1);
	//auto lon = P.longitudinal();
	double nu = P.spin_tune();
	double nu01, nu02, nu03;
	nu01 = -0.161656;
	nu02 = 3.23282e-08;
	nu03 = -6.15587e-15;
	nu = nu - nu01 - nu02 - nu03;
	cout<<"spin tune of deuteron in COSY: "
	//		<<setprecision(30)<<nu<<endl;
			<<nu<<endl;
	return 0;
}

int main_SCT_merge(int argc, char **argv)
{
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/sct_test";
	lv1 l1;
	acquisition aq;
	l1.SCT_merge(aq,path);
	return 0;
}

int main_c_vs_emittance(int argc, char **argv);

int main_h_vs_chromaticity(int argc, char **argv);

int main_w_vs_chromaticity(int argc, char **argv);

int main_cluster_particle(int argc, char **argv)
{
	int c0 = 0;
	if(argc<2) cout<<"No argument provided, assuming 0"<<endl;
	else c0 = stod(argv[1]);

	string path = "/home/shankar/PTEDM/COSY_particle/wp_Max/eMP_proton/r1";

	lv1 l1;

	l1.line(path,c0,0);

	return 0;
}

int main_cluster_beam(int argc, char **argv)
{
	string path = "/home/shankar/PTEDM/source_COSY/wp_Max/sct_proton_test";

	auto nc = bash::to_int("find " + path + " -type d -name \"c*\" | wc -l");
	cout<<"nc = "<<nc<<endl;
	proton_RF::beam(path,nc,0,0);

	return 0;
}


int main(int argc, char **argv)
{
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/sct_beam_based/test1";
	// string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/wp_Marcel/Deuteron/sphere";
	// string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/r_tunes_Marcel_Michael";
	// string wmpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/rrnd";
	//string path = "/home/shankar/PTEDM/hybrid_particle/eMP/r1";
	// string sample = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/r_tunes_Marcel_Michael/sample_SCT_beam";
	
	// auto nc = bash::to_int("find " + path + " -type d -name \"c*\" | wc -l");
	// cout<<"nc = "<<nc<<endl;
	// proton_RF::beam(path,nc,0,1);

	// lv1 l1;
	//acquisition aq;
	// lv2 az;
	
	//l1.point_RF(aq,1,164,1);
	
	// l1.multiplane(aq,path);
	// main_travel_time_analysis(argc, argv);
	// main_spin_tune_analysis_RFoff(argc, argv);
	// main_beam_sct_analysis(argc, argv);
	main_spin_tune_analysis(argc, argv);
	// main_c_vs_emittance(argc, argv);
	// main_h_vs_chromaticity(argc, argv);
	// main_w_vs_chromaticity(argc, argv);
	// main_SCT_merge(argc, argv);
	// main_multiplane_analysis(argc, argv);
	// main_multiplane_RFoff(argc, argv);
	// main_single_plane_analysis(argc, argv);
	// main_calculate_SCT(argc, argv);
	// main_reference(argc, argv);
	// main_tracking_chromaticity(argc, argv);
	
	// l1.line("",0,1);
	
	//l1.SCT_merge(aq,path);
	//bash::to_int("cp " + path + "/info.dat " + path + "/c1/info.dat");
	
	
	
	// auto d = 
	// proton_RF::line("/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/wp_Marcel_Michael/Deuteron/sphere/sample",9,0);
	// d = d + Data<double>(wmpath + "/line_error");

	// Data<double> d(path + "/data.txt");

	// Deuteron in COSY, Max Vitz
	// matrix<3,1> Q0 = {3.6333,3.719004,0.06089445};	
	// matrix<3,1> k0 = {-0.291997,0.363592,0.446509}; 

	// Deuteron in COSY, Michael Hartmann
	// matrix<3,1> Q0 = {3.609650,3.692445,5.07583153};	
	// matrix<3,1> k0 = {-0.28268,0.380387,0.380387};

	// Deuteron in COSY, Marcel Rosenthal (from Max's data)
	// matrix<3,1> Q0 = {3.623062,3.586645,-0.06198729};	
	// matrix<3,1> k0 = {-0.253797342919,0.344620755769,0.444622012417};

	// Deuteron in COSY, Marcel Rosenthal (from Michael's data)
	// matrix<3,1> Q0 = {3.620310,3.585020,-0.03205275};	
	// matrix<3,1> k0 = {-0.258838943366,0.345993764181,0.441952504399};

	// cout<<k0<<endl;
	
	// TApplication app("app", &argc, argv);

	// proton::load_offset("/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/Deuteron/rrndT");

	// proton_RF p("/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/wp_Marcel_Michael/Deuteron/opt1/c1");
	// p.betatron_tune(0,1);
	// p.betatron_tune(1,1);
	// p.phase_slip(0);
	// p.longitudinal(1);
	// p.spin_tune(1);
	// p.SCT(1);
	
	// az.COSY_tune_map(d,Q0,k0,path,1);
	// az.windmap(d,wmpath,1);
	// az.sampler_sphere(50,0.01,{0.672608520183,2.88510340087,-0.457974089331},sample,1);
	
	// TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	// rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	// app.Run();
	
	
	return 0;
}


/*
int main(int argc, char **argv)
{
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/sct_test/r1/c5";
	lv2 az;

	proton_RF p(path);
		
	TApplication app("app", &argc, argv);

	p.SCT();
	p.spin_tune();
	
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	return 0;
}

int main(int argc, char **argv)
{
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/eMP";
	string sample = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/sample";

	lv1 l1;
	acquisition aq;

	//l1.point_RF(aq,1,164,1);

	//l1.multiplane(aq,sample,path,9716);

	
	
	TApplication app("app", &argc, argv);
	
	proton_RF p(path);
	p.spin_tune(1);
	
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	

	return 0;
}
int main(int argc, char **argv)
{
	lv2 az;
	string rfpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/eMP";
	//string fpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/opt1";
	
	auto ssrc = bash::to_int("ls -d " + rfpath + "/r* | wc -l");
	string path;
	
	fstream f(rfpath + "/nu_opt_info", fstream::out | fstream::trunc);
	f<<"t_{s}"<<"\n"<<"#xi_{x}"<<"\n"<<"#xi_{y}"<<"\n"<<"#nu_{s}"<<"\n"<<"area"<<endl;
	f.close();
	f.open(rfpath + "/matrix_info", fstream::out | fstream::trunc);
	f<<"b[0]"<<"\n"<<"b[0]"<<"\n"<<"b[0]"<<"\n";
	f<<"M[0][0]"<<"\n"<<"M[0][1]"<<"\n";
	f<<"M[1][0]"<<"\n"<<"M[1][1]"<<"\n";
	f<<"M[2][0]"<<"\n"<<"M[2][1]"<<endl;
	f.close();
	f.open(rfpath + "/nu_opt", fstream::out | fstream::trunc);
	fstream f2(rfpath + "/matrix", fstream::out | fstream::trunc);
	for(int i=1; i<=ssrc; i++)
	{
		string path1 = rfpath + "/r" + to_string(i);
		cout<<"Currently working on "<<path1<<endl;
		
		auto nc = bash::to_int("ls -d " + path1 + "/c* | wc -l");
		cout<<"nc = "<<nc<<endl;
		auto ln = proton_RF::line(path1,nc,1);
		matrix<1,3> ep = {1.00164e-07,1.77605e-06,0.000990515};
		matrix<3,1> b;
		matrix<3,2> M;
		matrix<2,1> xi;
		tie(b,M,xi) = az.spin_model_2R(ln,path1);
		f<<ln.col("t_{s}")[0]<<"\t"<<xi[0]<<"\t"<<xi[1]<<"\t"
		<<ep*b + ep*M*xi<<"\t";
		f2<<b[0]<<"\t"<<b[1]<<"\t"<<b[2]<<"\t"
		<<M[0][0]<<"\t"<<M[0][1]<<"\t"
		<<M[1][0]<<"\t"<<M[1][1]<<"\t"
		<<M[2][0]<<"\t"<<M[2][1]<<endl;
		
		matrix<3,3> Q;
		matrix<2,2> m1, m2, m3;
		Q.sub<3,2>(0,0) = M; Q(2) = -b;
		m1 = M.sub<2,2>(1,0);
		m2[0] = M[2]; m2[1] = M[0];
		m3 = M.sub<2,2>(0,0);
		cout<<"Q = "<<Q<<endl;
		f<<pow(det(Q),2)/(2*det(m1)*det(m2)*det(m3))<<endl;
	}
	f.close();
	f2.close();
	
	Data<double> oD(rfpath + "/nu_opt");
	
	auto a1 = proton_RF::line(fpath,121,1);
	
	TApplication app("app", &argc, argv);
	az.optimize(oD,rfpath,1);
	
	az.spin_model_2R(a1,fpath);
	
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	
return 0;
}


int main(int argc, char **argv)
{
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data";
	lv2 az;
		
	TApplication app("app", &argc, argv);

	az.sampler(10000,17,1.0,path,1);
	
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();

	return 0;
}

int main(int argc, char **argv)
{
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/r0/c1";

	TApplication app("app", &argc, argv);
	
	proton_RF p(path);
	p.spin_tune(1);

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();

	return 0;
}
int main(int argc, char **argv)
{
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/rE";
	int nc = bash::to_int("ls -d " + path + "/c* | wc -l");
	proton::line(path,nc,1);
	lv2 az;

	proton_RF::line("/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/eSM/r1",121,1);
		
	TApplication app("app", &argc, argv);

	az.tune_vs_emittance(path + "/line");
	
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	return 0;
}

int main(int argc, char **argv)
{
	acquisition aq;
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/sct3";

	//lv1 l1;

	//cout << bash::to_string("echo \"$PATH\"") << endl;
	//l1.SCT_merge(aq,path);

	proton_RF p(path + "/c1");
	p.setChromaticity(-0.14148808,19.906467);

	TApplication app("app", &argc, argv);
	
	p.SCT(1);
	p.spin_tune(1);
	
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	return 0;
}


int main(int argc, char **argv)
{
	proton P("/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/reference/c0");
	//auto CSx = P.Courant_Snyder(0);
	//auto CSy = P.Courant_Snyder(1);
	//auto lon = P.longitudinal();
	double nu = P.spin_tune();
	double nu01, nu02, nu03;
	nu01 = 0.0786942;
	nu02 = -4.70913e-08;
	nu03 = -5.97969e-15;
	nu = nu - nu01 - nu02 - nu03;
	cout<<"spin tune of proton in magnetic prototype: "
	//		<<setprecision(30)<<nu<<endl;
			<<nu<<endl;
}

int main(int argc, char **argv)
{
	string spinfile = "/Users/utente/JEDI/PTEDM/1.0/cluster/cluster_data/beam3/spin0.673461_1.03817_0.0";
	
	Data<double> D(spinfile);
	
	TApplication app("app", &argc, argv);
	
	proton_RF::simple_SCT(spinfile,801.6,700.74037);
    
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	return 0;
}

int main(int argc, char **argv)
{
	string wmpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/r1";

	TApplication app("app", &argc, argv);

	lv1 l1;
	l1.windmap(wmpath,1);

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();

	return 0;
}
int main(int argc, char **argv)
{
	lv2 az;
	string rfpath = "/home/shankar/JEDI/PTEDM/1.0/source_COSY/build/data/eSM3";
	string fpath = "/home/shankar/JEDI/PTEDM/1.0/source_COSY/build/data/opt1";

	auto ssr = bash::to_ss("ls -d " + rfpath + "/r*");
	string path;

	fstream f(rfpath + "/nu_opt_info", fstream::out | fstream::trunc);
	f<<"t_{s}"<<"\n"<<"#xi_{x}"<<"\n"<<"#xi_{y}"<<"\n"<<"#nu_{s}"<<"\n"<<"area"<<endl;
	f.close();
	f.open(rfpath + "/nu_opt", fstream::out | fstream::trunc);
	while(ssr>>path)
	{
		cout<<"Currently working on "<<path<<endl;

		auto nc = bash::to_int("ls -d " + path + "/c* | wc -l");
		cout<<"nc = "<<nc<<endl;
		auto ln = proton_RF::line(path,nc,0);
		matrix<1,3> ep = {1.00164e-07,1.77605e-06,0.000990515};
		matrix<3,1> b;
		matrix<3,2> M;
		matrix<2,1> xi;
		tie(b,M,xi) = az.spin_model_2R(ln,path);
		f<<ln.col("t_{s}")[0]<<"\t"<<xi[0]<<"\t"<<xi[1]<<"\t"
			<<ep*b + ep*M*xi<<"\t";

		matrix<3,3> Q;
		matrix<2,2> m1, m2, m3;
		Q.sub<3,2>(0,0) = M; Q(2) = -b;
		m1 = M.sub<2,2>(1,0);
		m2[0] = M[2]; m2[1] = M[0];
		m3 = M.sub<2,2>(0,0);
		cout<<"Q = "<<Q<<endl;
		f<<pow(det(Q),2)/(2*det(m1)*det(m2)*det(m3))<<endl;
	}
	f.close();

	Data<double> oD(rfpath + "/nu_opt");

	auto a1 = proton_RF::line(fpath,121,1);
	az.spin_model_2R(a1,fpath);

	TApplication app("app", &argc, argv);
	az.optimize(oD,rfpath,1);

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();

	return 0;
}
int main(int argc, char **argv)
{
	lv2 az;
	string path = "/home/shankar/JEDI/PTEDM/1.0/source_COSY/build/data/e2";
	string rfpath = "/home/shankar/JEDI/PTEDM/1.0/source_COSY/build/data/RF_multiplane";

	//lv1 l1;
	//l1.grid(path,1);


	Data<double> D(path + "/grid");//, D0(path + "/grid0");
	//proton::line("/home/shankar/JEDI/PTEDM/1.0/cluster/simulator/source_prototype/eS3/r1",20,1);
	auto rfD = proton_RF::line(rfpath,2000,0);

	az.spin_model_2(D,path);
	TApplication app("app", &argc, argv);
	az.spin_model_RF(rfD,rfpath,1);
	//az.SM_RF_test2(rfD,rfpath,1);

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();


	return 0;
}
int main(int argc, char **argv)
{
	//acquisition aq;
	//acquisition aq0 = aq;
	//aq0.set_zero_flag(1);
	//cout<<"zero flags: "<<aq.zero()<<" "<<aq0.zero()<<endl;
	lv2 az;
	//string path = "/home/shankar/JEDI/PTEDM/1.0/cluster/cluster_data/eS3";
	string path = "/home/shankar/JEDI/PTEDM/1.0/cluster/simulator/source_prototype/eS3";
	string wmpath = "/home/shankar/JEDI/PTEDM/1.0/source_prototype2/build/data/ew";
	string wmepath = "/home/shankar/JEDI/PTEDM/1.0/source_prototype2/build/data/ewe";
	string epath = "/home/shankar/JEDI/PTEDM/1.0/source_prototype2/build/data/re";
	string exp_path = "/home/shankar/JEDI/PTEDM/1.0/source_prototype2/build/data/ew/r1/c11";
	string rfpath = "/home/shankar/JEDI/PTEDM/1.0/source_prototype2/build/data/rtt";
	matrix<3,1> xi0;
	matrix<3,3> Xi;

	//az.tune_vs_emittance(epath);
	tie(xi0,Xi) = az.windmap(wmpath);

	//lv1 l1;
	//l1.grid(path,0);
	//l1.grid(aq0,path,0);

	Data<double> D(path + "/grid");//, D0(path + "/grid0");
	az.spin_model(D,D,path);
	auto rfD = proton_RF::line(rfpath,200,0);
	//proton::line("/home/shankar/JEDI/PTEDM/1.0/cluster/simulator/source_prototype/eS3/r1",20,1);

	TApplication app("app", &argc, argv);
	//az.windmap_e(wmepath);
	az.SM_RF_test2(rfD,rfpath,1);

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();

	return 0;
}

int main(int argc, char **argv)
{
	//acquisition aq;
	//acquisition aq0 = aq;
	//aq0.set_zero_flag(1);
	//cout<<"zero flags: "<<aq.zero()<<" "<<aq0.zero()<<endl;
	lv2 az;
	//string path = "/home/shankar/JEDI/PTEDM/1.0/cluster/cluster_data/eS3";
	string path = "/home/shankar/JEDI/PTEDM/1.0/cluster/simulator/source_prototype/eS3";
	string wmpath = "/home/shankar/JEDI/PTEDM/1.0/source_prototype2/build/data/ew";
	string wmepath = "/home/shankar/JEDI/PTEDM/1.0/source_prototype2/build/data/ewe";
	string epath = "/home/shankar/JEDI/PTEDM/1.0/source_prototype2/build/data/re";
	string exp_path = "/home/shankar/JEDI/PTEDM/1.0/source_prototype2/build/data/ew/r1/c11";
	string rfpath = "/home/shankar/JEDI/PTEDM/1.0/source_prototype2/build/data/rtt";
	matrix<3,1> xi0;
	matrix<3,3> Xi;

	//az.tune_vs_emittance(epath);
	tie(xi0,Xi) = az.windmap(wmpath);

	//lv1 l1;
	//l1.grid(path,0);
	//l1.grid(aq0,path,0);

	Data<double> D(path + "/grid");//, D0(path + "/grid0");
	az.spin_model(D,D,path);
	auto rfD = proton_RF::line(rfpath,200,0);
	//proton::line("/home/shankar/JEDI/PTEDM/1.0/cluster/simulator/source_prototype/eS3/r1",20,1);

	TApplication app("app", &argc, argv);
	//az.windmap_e(wmepath);
	az.SM_RF_test2(rfD,rfpath,1);

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();



	return 0;
}
int main(int argc, char **argv)
{
	string rfpath = "/home/shankar/JEDI/PTEDM/1.0/source_magnetic/build/data/rtt";
	auto rfD = proton_RF::line(rfpath,200,1,"magnetic");
	return 0;
}
int main(int argc, char **argv)
{
	proton P("/home/shankar/JEDI/PTEDM/1.0/cluster/cluster_data/beam1/r1/c1");
	P.set_mode("electric");
	double nu = P.spin_tune();
	cout<<"spin tune of proton in electric prototype: "
	//		<<setprecision(30)<<nu<<endl;
			<<nu<<endl;
}
int main()
{
	acquisition aq;
	vector<vector<double>> C(3,vector<double>());
	aq.mom(3,3,C[0],C[1],C[2]);
	for(int i=0; i<C[0].size(); i++) cout<<C[0][i]<<"\t"<<C[1][i]<<"\t"<<C[2][i]<<endl;
	return 0;
}
int main(int argc, char **argv)
{
	string path = "/home/Rahul/JEDI/PTEDM/1.0/cluster/cluster_data/eS";
	//string path = "/run/media/Rahul/Hiruzen/Rahul's HDD/JEDI/prototype_singles/Spin Model Simulations/Qx1.855Qy1.023/eS2";
	//int NS = 1000;
	TFile* file = TFile::Open((path + "/eS.root").data(),"UPDATE");//);//
	lv2 az(path,file);
	//lv1 l1;
	//cout<<"YAY!!!"<<endl;

	TApplication app("app", &argc, argv);
	//l1.grid(0,0,file,1);
	az.spin_model_2();
	//az.sampler();//
	//az.windmap();
	file->Close();
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();


	return 0;
}

int main(int argc, char **argv)
{
	string path = "/run/media/Rahul/Hiruzen/Rahul's HDD/JEDI/prototype_singles/Spin Model Simulations/Qx1.855Qy1.023/eS2";
	//int NS = 1000;
	TFile* file = TFile::Open((path + "/eS.root").data());//,"UPDATE");//
	lv2 az(path,file);

	az.spin_model_2();
	TApplication app("app", &argc, argv);
	az.sampler();
	file->Close();
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();


	return 0;
}
*/


