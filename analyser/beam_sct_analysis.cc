#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"

#include<lv2.hh>

int main_beam_sct_analysis(int argc, char **argv)
{
	fstream f;

	lv2 az;
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/sct_beam_based";
	
	auto nc = bash::to_int("find " + path + " -type d -name \"c*\" | wc -l");
	cout<<"nc = "<<nc<<endl;
	auto ln = proton_RF::beam(path,nc,0,1) + Data<double>(path + "/beam_error");
	
	ln.clean_nans();

	f.open(path + "/beam_matrix_info", fstream::out | fstream::trunc);
	f<<"m0"<<endl;
	f<<"m[0]"<<"\n"<<"m[1]"<<"\n"<<"m[2]"<<endl;
	f<<"M[0][0]"<<"\n"<<"M[0][1]"<<"\n"<<"M[0][2]"<<endl;
	f<<"M[1][0]"<<"\n"<<"M[1][1]"<<"\n"<<"M[1][2]"<<endl;
	f<<"M[2][0]"<<"\n"<<"M[2][1]"<<"\n"<<"M[2][2]"<<endl;
	f.close();
	f.open(path + "/beam_matrix", fstream::out | fstream::trunc);

	TApplication app("app", &argc, argv);

	az.optimize_beam(ln,path);
	// cout<<"YAY!!"<<endl;

	f<<az.m0<<endl;
	f<<az.m[0]<<"\n"<<az.m[1]<<"\n"<<az.m[2]<<endl;
	f<<az.M[0][0]<<"\t"<<az.M[0][1]<<"\t"<<az.M[0][2]<<endl;
	f<<az.M[1][0]<<"\t"<<az.M[1][1]<<"\t"<<az.M[1][2]<<endl;
	f<<az.M[2][0]<<"\t"<<az.M[2][1]<<"\t"<<az.M[2][2]<<endl;
	f.close();

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	return 0;
}
