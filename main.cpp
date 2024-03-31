#include "TF1.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"

#include<proton.hh>
#include<lv1.hh>
#include<lv2.hh>
#include<analyser.hh>


using namespace std;
int main(int argc, char **argv)
{
	int NS = 343; //1000; //
	string path = "/run/media/Rahul/Hiruzen/Rahul's HDD/JEDI/prototype_singles/Spin Model Simulations/Qx1.855Qy1.023/1";
	//string path = "/home/Rahul/JEDI/PTEDM/1.0/source_prototype2/build/data/Qx1.855Qy1.023"; //eS"; //
	//string path = "/home/Rahul/JEDI/PTEDM/1.0/cluster/cluster_data/eS";
	analyser az(path,NS);

	matrix<1,3> ep;
	matrix<3,3> Xi, M;
	matrix<3,1> a, xi, xi0;

	fstream f;
	f.open(path + "/spin_model",fstream::in);
	bool sme_flag = 0;
	if(f.is_open())
	{
		f>>az.sig0;
		for(int i=0; i<3; i++) f>>az.a[i];
		for(int i=0; i<3; i++) for(int j=0; j<3; j++) f>>az.V[i][j];
		f>>az.eta0;
		for(int i=0; i<3; i++) for(int j=0; j<3; j++) f>>az.T[i][j];
		az.M = az.V - (az.sig0/az.eta0)*az.T;
		f.close();
		M = az.M;
		a = az.a;
		sme_flag = 1;
	}
	else tie(a,M,xi) = az.spin_model();
	cout<<"a = "<<a<<endl;
	cout<<"M = "<<M<<endl;

	f.open(path + "/matrix",fstream::in);
	if(f.is_open())
	{
		for(int i=0; i<3; i++) for(int j=0; j<3; j++) f>>Xi[i][j];
		for(int i=0; i<3; i++) f>>xi0[i];
		az.Xi = Xi; az.xi0 = xi0;
		f.close();
	}
	else tie(xi0,Xi) = az.windmap();
	cout<<"xi0 = "<<xi0<<endl;
	cout<<"Xi = "<<Xi<<endl;
	auto op = (Xi^-1) * (xi - xi0);

	f.open(path + "/op.txt",fstream::out | fstream::trunc);
	f<<"k.sexf      =   "<<op[0]
	<<"\nk.sexd      =   "<<op[1]
	<<"\nk.sexs      =   "<<op[2]<<endl;
	f<<endl;
	f<<op[0]<<"\t"<<op[1]<<"\t"<<op[2]<<"\t1000000"<<endl;
	f.close();

	matrix<2,1> xi_o = {-1.35123,-4.527};
	double ts_o = -0.976576;

	//Data<double> Dte = proton_RF::line(path + "/rte",64);
	Data<double> Dtc = proton_RF::line(path + "/rtc",NS); //497); //

	auto t2_o = (Xi.sub<2,2>(0,0)^-1) * (xi_o - xi0.sub<2,1>() - ts_o*Xi.sub<2,1>(0,2));
	matrix<3,1> t_o = {t2_o[0],t2_o[1],ts_o};
	cout<<"Optimized sextupole fields expected at: \n"<<op<<endl;
	cout<<"Calculated optimized optical settings: \n"<<xi0 + Xi*op<<endl;
	cout<<"Measured optimized sextupole fields: \n"<<t_o<<endl;
	cout<<"Measured optimized optical settings: \n"<<xi0 + Xi*t_o<<endl;

	TApplication app("app", &argc, argv);
	//az.sampler();//windmap();//
	//if(sme_flag) az.spin_model_errors();
	az.SM_RF_test2(Dtc);
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	/*
	*/
	return 0;
}

/*
int main(int argc, char **argv)
{
	string path = "/home/Rahul/JEDI/PTEDM/1.0/source_prototype_v4/build/data/eS";
	analyser az(path);

	Data<double> Dte = proton_RF::line(path + "/rtep",121);

	TApplication app("app", &argc, argv);
	az.E_RF_test(Dte);
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();

	return 0;
}
*/


/*
int main(int argc, char **argv)
{
	if(argc < 2) { cerr<<"Analyser requires a path!"<<endl; return 1; }
	string path = argv[1];

	analyser az(path);

	matrix<3,3> Xi;
	matrix<3,1> xi0;

	tie(xi0,Xi) = az.windmap();

	auto oo = (Xi^-1) * (-xi0);
	cout<<oo[0]<<" "<<oo[1]<<" "<<oo[2]<<endl;

	return 0;
}
*/

/*
int main(int argc, char **argv)
{
	TApplication app("app", &argc, argv);
	string path = "/home/Rahul/JEDI/PTEDM/1.0/source_prototype2/build/data/e3/c10";

	proton P(path);
	P.Courant_Snyder();
	P.Courant_Snyder(1);

	cout << "done!" << endl;
	app.Run();
	return 0;
}
int main(int argc, char **argv)
{
	string path = "/run/media/Rahul/Windows/LINUX_EXCHANGE/PTEDM/prototype_singles/e6";
	int N = 51;
	auto D = proton::line(path,N);

	cout << "done!" << endl;

	return 0;
}

int main(int argc, char **argv)
{
	string path = "/home/Rahul/JEDI/PTEDM/1.0/source_prototype2/build/data/e8";

	int N = 81;
	auto D = proton::line(path,N);

	fstream f;
	f.open((path + "/stv").data(),fstream::out | fstream::trunc);
	for(int i=1; i<=N; i++)
	{
		Data<string> info(path + "/c" + to_string(i) + "/info.dat");
		f<<info[2][2]<<"\t"<<info[2][3]<<"\t"<<D[3][i]<<"\t"<<D[4][i]<<"\t"<<D[7][i]<<endl;
	}
	f.close();
	f.open((path + "/stv_info").data(),fstream::out | fstream::trunc);
	f<<"Qx\nQy\nep_x\nep_y\nnu_s"<<endl;
	f.close();

	cout << "done!" << endl;

	return 0;
}
int main(int argc, char **argv)
{
	TApplication app("app", &argc, argv);
	string path = "/home/Rahul/JEDI/PTEDM/1.0/source_prototype2/build/data/e1/c0";

	proton_RF P(path);

	P.Courant_Snyder();
	P.phase_slip();
	P.spin_tune();
	app.Run();
	return 0;
}

int main(int argc, char **argv)
{
	string path = "/home/Rahul/JEDI/PTEDM/1.0/source_prototype2/build/data";

	fstream f;
	f.open((path + "/spin_tune_variation").data(),fstream::out | fstream::trunc);
	for(int i=0; i<7; i++)
	{
		auto D = proton::line(path + "/e" + to_string(i+1),64);
		auto stv = proton::spin_tune_variation(path + "/e" + to_string(i+1),D);

		Data<string> info(path + "/e" + to_string(i+1) + "/c1/info.dat");

		f<<info[2][2]<<"\t"<<info[2][3]<<"\t"<<stv[0]<<"\t"<<stv[1]<<"\t"<<stv[2]<<endl;
	}
	f.close();

	cout << "done!" << endl;

	return 0;
}
*/
	/*
int main(int argc, char **argv)
{
	TApplication app("app", &argc, argv);
	TCanvas *c = new TCanvas("c", "Something", 0, 0, 800, 600);
	TF1 *f1 = new TF1("f1", "sin(x)", -5, 5);
	f1->SetLineColor(kBlue + 1);
	f1->SetTitle("My graph;x; sin(x)");
	f1->Draw();
	c->Modified(); c->Update();
	TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	return 0;
}
*/
