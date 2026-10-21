#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"

#include<lv2.hh>

int main_travel_time_analysis(int argc, char **argv)
{
	lv2 az;
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/Deuteron/rrnd_mis";
	//string fpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/opt1";
	
	auto nc = bash::to_int("ls -d " + path + "/c* | wc -l");
	cout<<"nc = "<<nc<<endl;
	auto ln = proton::line(path,nc,1);
	
	TApplication app("app", &argc, argv);

	az.travel_time_model(ln,path);
	proton p(path + "/c1");
	cout<<"p.Tp_model() = "<<p.Tp_model()<<endl;
	p.phase_slip(1);

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	return 0;
}
