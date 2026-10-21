#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"

#include<lv2.hh>

int main_eta_1_RF(int argc, char **argv)
{
    lv2 az;
    string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/rE";
    //string fpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/opt1";

    auto nc = bash::to_int("find " + path + " -type d -name \"c*\" | wc -l");
    cout<<"nc = "<<nc<<endl;

    auto ln = proton_RF::line(path, nc);

    TApplication app("app", &argc, argv);

    //az.windmap(ln,path,1);

    TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
    rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    app.Run();
    
    return 0;
}