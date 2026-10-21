#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"

#include<lv2.hh>

int main_spin_tune_analysis_RFoff(int argc, char **argv)
{
	lv2 az;
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/Deuteron/rrndT";
	//string fpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/opt1";

	auto nc = bash::to_int("find " + path + " -type d -name \"c*\" | wc -l");
	cout<<"nc = "<<nc<<endl;
	auto ln = proton_RF::line(path,nc,1);
	
	TApplication app("app", &argc, argv);

	//az.spin_model_RFoff(ln,path);
	az.windmap(ln,ln.col("#eta_{1}"),path,1);

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	return 0;
}
