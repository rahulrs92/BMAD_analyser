#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"

#include<lv2.hh>

int main_multiplane_RFoff(int argc, char **argv)
{
	string rfpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/eMP_T";
	//string fpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/opt1";
	
	auto ssrc = bash::to_int("find " + rfpath + " -type d -name \"r*\" | wc -l");
	cout<<"ssrc = "<<ssrc<<endl;
	string path;
	
	fstream f(rfpath + "/sm_matrix_info", fstream::out | fstream::trunc);
	f<<"#sigma_{0}"<<"\n"<<"h[0]"<<"\n"<<"h[1]"<<"\n"<<"h[2]"<<"\n";
	f<<"H[0][0]"<<"\n"<<"H[0][1]"<<"\n";
	f<<"H[1][0]"<<"\n"<<"H[1][1]"<<"\n";
	f<<"H[2][0]"<<"\n"<<"H[2][1]"<<endl;
	f.close();
	f.open(rfpath + "/tt_matrix_info", fstream::out | fstream::trunc);
	f<<"#eta_{0}"<<"\n"<<"u[0]"<<"\n"<<"u[1]"<<"\n"<<"u[2]"<<"\n";
	f<<"U[0][0]"<<"\n"<<"U[0][1]"<<"\n";
	f<<"U[1][0]"<<"\n"<<"U[1][1]"<<"\n";
	f<<"U[2][0]"<<"\n"<<"U[2][1]"<<endl;
	f.close();
	f.open(rfpath + "/xi_opt_info", fstream::out | fstream::trunc);
	f<<"#xi_{x}"<<"\n"<<"#xi_{y}"<<endl;
	f.close();
	f.open(rfpath + "/sm_matrix", fstream::out | fstream::trunc);
	fstream f2(rfpath + "/tt_matrix", fstream::out | fstream::trunc);
	fstream f3(rfpath + "/xi_opt", fstream::out | fstream::trunc);
	
	
	for(int i=1; i<=ssrc; i++)
	{
		lv2 az;
		string path1 = rfpath + "/r" + to_string(i);
		cout<<"Currently working on "<<path1<<endl;
		
		auto nc = bash::to_int("ls -d " + path1 + "/c* | wc -l");
		cout<<"nc = "<<nc<<endl;
		auto ln = proton::line(path1,nc,1);
		matrix<1,3> ep = {1.00164e-07,1.77605e-06,0.000990515};
		double sig0;
		matrix<3,1> h;
		matrix<3,2> H;
		matrix<2,1> xi;
		tie(sig0,h,H,xi) = az.spin_model_2(ln,path1);
		double eta0 = az.eta0;
		matrix<3,1> u = az.u, w = az.w;
		matrix<3,2> U = az.U, W = az.W;
		xi = ((W.tp()*W)^-1) * W.tp() * (-w);
		f<<sig0<<"\t"<<h[0]<<"\t"<<h[1]<<"\t"<<h[2]<<"\t"
		<<H[0][0]<<"\t"<<H[0][1]<<"\t"
		<<H[1][0]<<"\t"<<H[1][1]<<"\t"
		<<H[2][0]<<"\t"<<H[2][1]<<endl;
		f2<<eta0<<"\t"<<u[0]<<"\t"<<u[1]<<"\t"<<u[2]<<"\t"
		<<U[0][0]<<"\t"<<U[0][1]<<"\t"
		<<U[1][0]<<"\t"<<U[1][1]<<"\t"
		<<U[2][0]<<"\t"<<U[2][1]<<endl;
		f3<<xi[0]<<"\t"<<xi[1]<<endl;
	}
	f.close();
	f2.close();
	f3.close();
	
	/*
	TApplication app("app", &argc, argv);
	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	*/
	
	return 0;
}
