#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"

#include<lv2.hh>

int main_multiplane_analysis(int argc, char **argv)
{
	lv2 az;
	string rfpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/eMP";
	//string fpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/opt1";

	auto ssrc = bash::to_int("find " + rfpath + " -type d -name \"r*\" | wc -l");
	cout<<"ssrc = "<<ssrc<<endl;
	string path;
	
	fstream f(rfpath + "/nu_opt_info", fstream::out | fstream::trunc);
	f<<"t_{s}"<<"\n"<<"#xi_{x}"<<"\n"<<"#xi_{y}"<<"\n"<<"#nu_{s}"<<"\n"<<"area"<<endl;
	f.close();
	f.open(rfpath + "/matrix_info", fstream::out | fstream::trunc);
	f<<"w[0]"<<"\n"<<"w[0]"<<"\n"<<"w[0]"<<"\n";
	f<<"W[0][0]"<<"\n"<<"W[0][1]"<<"\n";
	f<<"W[1][0]"<<"\n"<<"W[1][1]"<<"\n";
	f<<"W[2][0]"<<"\n"<<"W[2][1]"<<endl;
	f.close();
	f.open(rfpath + "/nu_opt", fstream::out | fstream::trunc);
	fstream f2(rfpath + "/matrix", fstream::out | fstream::trunc);

	matrix<2,1> xi;

	for(int i=1; i<=ssrc; i++)
	{
		string path1 = rfpath + "/r" + to_string(i);
		cout<<"Currently working on "<<path1<<endl;
		
		auto nc = bash::to_int("ls -d " + path1 + "/c* | wc -l");
		cout<<"nc = "<<nc<<endl;
		auto ln = proton_RF::line(path1,nc,1);
		matrix<1,3> ep = {1.00164e-07,1.77605e-06,0.000990515};
		matrix<3,1> w;
		matrix<3,2> W;
		tie(w,W,xi) = az.spin_model_2R(ln,path1);
		f<<ln.col("t_{s}")[0]<<"\t"<<xi[0]<<"\t"<<xi[1]<<"\t"
		<<ep*w + ep*W*xi<<"\t";
		f2<<w[0]<<"\t"<<w[1]<<"\t"<<w[2]<<"\t"
		<<W[0][0]<<"\t"<<W[0][1]<<"\t"
		<<W[1][0]<<"\t"<<W[1][1]<<"\t"
		<<W[2][0]<<"\t"<<W[2][1]<<endl;
		
		matrix<3,3> Q;
		matrix<2,2> w1, w2, w3;
		Q.sub<3,2>(0,0) = W; Q(2) = -w;
		w1 = W.sub<2,2>(1,0);
		w2[0] = W[2]; w2[1] = W[0];
		w3 = W.sub<2,2>(0,0);
		//cout<<"Q = "<<Q<<endl;
		f<<pow(det(Q),2)/(2*det(w1)*det(w2)*det(w3))<<endl;
	}
	f.close();
	f2.close();
	
	Data<double> oMP(rfpath + "/multiplane");
	matrix<2,1> xi0;

	lv1 l1;
	TApplication app("app", &argc, argv);

	if (proton_RF::mode == "hybrid")
	{
		matrix<2,2> Xi;
		tie(xi0,Xi) = l1.windmap_2D(oMP,rfpath,1);
		
		xi = {0.673461, 1.03817};
		cout<<"xi = "<<xi<<endl;
		cout<<"xi0 = "<<xi0<<endl;
		cout<<"Xi = "<<Xi<<endl;

		auto t = (Xi^-1) * (xi - xi0);

		cout<<"t = "<<t.tp()<<endl;
	}
	else
	{
		matrix<2,3> Xi;
		tie(xi0,Xi) = l1.windmap_plane(oMP,rfpath,1);
		matrix<2,2> Xi2 = Xi.sub<2,2>(); matrix<2,1> Xi3 = Xi(2);
		
		Data<double> oD(rfpath + "/nu_opt");
		Data<double> oM(rfpath + "/matrix");
		//auto a1 = proton_RF::line(fpath,121,1);

		auto a = az.optimize(oD,oM,rfpath,1);
		matrix<2,1> xi = {a[0],a[1]};
		auto t12 = (Xi2^-1) * (xi - xi0 - Xi3*a[2]);
		matrix<3,1> t = {t12[0],t12[1],a[2]};
		
		cout<<"t = "<<t.tp()<<endl;
	}
	//matrix<2,3> Xi;

	//proton_RF p("/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/sct4/c5");
	//p.SCT();
	//p.spin_tune();
	
	//az.spin_model_2R(a1,fpath);
	
	
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	return 0;
}
