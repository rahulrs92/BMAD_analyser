#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"

#include<lv2.hh>

int main_single_plane_analysis(int argc, char **argv)
{
	lv2 az;
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_hybrid/build/data/eMP";
	//string fpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/opt1";

	fstream f(path + "/matrix", fstream::out | fstream::trunc);
	
	Data<double> D(path + "/plane");
	
	matrix<1,3> ep = {1.00164e-07,1.77605e-06,0.000990515};
	matrix<3,1> w;
	matrix<3,2> W;
	matrix<2,1> xi, xi0;
	tie(w,W,xi) = az.spin_model_2R(D,path);
	f<<w[0]<<"\t"<<w[1]<<"\t"<<w[2]<<"\t"
	<<W[0][0]<<"\t"<<W[0][1]<<"\t"
	<<W[1][0]<<"\t"<<W[1][1]<<"\t"
	<<W[2][0]<<"\t"<<W[2][1]<<endl;
	f.close();
	
	TApplication app("app", &argc, argv);

	lv1 l1;
	matrix<2,2> Xi;
	tie(xi0,Xi) = l1.windmap_2D(D,path,1);
	
	cout<<"xi = "<<xi<<endl;
	xi = {0.673461, 1.03817};
	cout<<"xi = "<<xi<<endl;
	cout<<"xi0 = "<<xi0<<endl;
	cout<<"Xi = "<<Xi<<endl;

	auto t = (Xi^-1) * (xi - xi0);

	cout<<"t = "<<t.tp()<<endl;

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	return 0;
}
