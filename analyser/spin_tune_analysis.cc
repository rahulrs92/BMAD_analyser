#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"

#include<lv2.hh>

int main_spin_tune_analysis(int argc, char **argv)
{
	lv2 az;
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/Deuteron/sphere";
	string rpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY_mis/build/data/Deuteron/ref_sphere";
	string wmpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/Deuteron/rrnd";
	// string sample = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/Deuteron/sample_SCT_beam";
	// string fpath = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/opt1";
	
	auto nc = bash::to_int("find " + path + " -type d -name \"c*\" | wc -l");
	cout<<"nc = "<<nc<<endl;
	Data<double> ln = proton_RF::line(path,nc,0,0) + Data<double>(path + "/line_error");
	// ln = ln ;
	// cout<<"YAY!!!"<<endl;

	if(proton_RF::mode == "COSY_mis")
	{
		auto r_nc = bash::to_int("find " + rpath + " -type d -name \"c*\" | wc -l");
		cout<<"r_nc = "<<r_nc<<endl;
		auto r_ln = proton_RF::line(rpath,r_nc,0,0);
		r_ln = r_ln + Data<double>(rpath + "/line_error");
		
		// zeroing unused columns
		for(int j=0; j<r_ln.order()[0]; ++j) r_ln[3][j+1] = 0.0;
		for(int i=6; i<11; ++i) for(int j=0; j<r_ln.order()[0]; ++j) r_ln[i+1][j+1] = 0.0;
		for(int j=0; j<r_ln.order()[0]; ++j) r_ln[15][j+1] = 0.0;
		for(int i=6; i<11; ++i) for(int j=0; j<r_ln.order()[0]; ++j) r_ln[i+13][j+1] = 0.0;
		// cout<<"r_ln after zeroing unused columns: "<<endl;
		// cout<<r_ln<<endl;
		r_ln.clean_nans();
		// cout<<ln.col("#epsilon_{x}")<<endl;

		az.tensor_orbit = az.spin_orbit_correction(r_ln,rpath,1);
	}
	//Data<double> ln(path + "/multiplane");
	//Data<double> ln(path + "/line");


	// clean nans
	ln.clean_nans();

	fstream f;
	// fstream f(wmpath + "/wm_matrix", fstream::in);
	// f>>az.xi03[0]>>az.xi03[1]>>az.xi03[2]
	// >>az.Xi3[0][0]>>az.Xi3[0][1]>>az.Xi3[0][2]
	// >>az.Xi3[1][0]>>az.Xi3[1][1]>>az.Xi3[1][2]
	// >>az.Xi3[2][0]>>az.Xi3[2][1]>>az.Xi3[2][2];
	// f.close();

	// f.open(wmpath + "/wm_matrix_error", fstream::in);
	// f>>az.xi03_e[0]>>az.xi03_e[1]>>az.xi03_e[2]
	// >>az.Xi3_e[0][0]>>az.Xi3_e[0][1]>>az.Xi3_e[0][2]
	// >>az.Xi3_e[1][0]>>az.Xi3_e[1][1]>>az.Xi3_e[1][2]
	// >>az.Xi3_e[2][0]>>az.Xi3_e[2][1]>>az.Xi3_e[2][2];
	// f.close();

	if(proton_RF::mode == "COSY_mis")
	{
		f.open(wmpath + "/wm_tensor", fstream::in);
		for(int i=0; i<3; ++i) for(int j=0; j<10; ++j) f>>az.tensor_Xi[i][j];
		f.close();
		
		f.open(wmpath + "/wm_tensor_error", fstream::in);
		for(int i=0; i<3; ++i) for(int j=0; j<10; ++j) f>>az.e_tensor_Xi[i][j];
		f.close();
	}

	f.open(path + "/sm3_matrix_info", fstream::out | fstream::trunc);
	f<<"w[0]"<<"\n"<<"w[1]"<<"\n"<<"w[2]"<<"\n";
	f<<"W[0][0]"<<"\n"<<"W[0][1]"<<"\n"<<"W[0][2]"<<"\n";
	f<<"W[1][0]"<<"\n"<<"W[1][1]"<<"\n"<<"W[1][2]"<<"\n";
	f<<"W[2][0]"<<"\n"<<"W[2][1]"<<"\n"<<"W[2][2]"<<endl;
	f.close();
	f.open(path + "/sm3_matrix", fstream::out | fstream::trunc);

	TApplication app("app", &argc, argv);

	az.windmap(ln,path,1);
	
	if(proton_RF::mode == "COSY_mis") az.spin_model_tensor(ln,path);
	else az.spin_model_3(ln,path);
	
	f<<az.w[0]<<"\t"<<az.w[1]<<"\t"<<az.w[2]<<"\t"
	<<az.W3[0][0]<<"\t"<<az.W3[0][1]<<"\t"<<az.W3[0][2]<<"\t"
	<<az.W3[1][0]<<"\t"<<az.W3[1][1]<<"\t"<<az.W3[1][2]<<"\t"
	<<az.W3[2][0]<<"\t"<<az.W3[2][1]<<"\t"<<az.W3[2][2]<<endl;
	f.close();
	
	/*
	az.sampler_sphere(50,0.05,{0.578567121117,3.64403694752,-0.636129393325},sample,1); // COSY deuteron
	az.sampler_sphere(50,0.05,{1.73615868366,-2.35158847301,0.709609199553},sample,1);	// COSY proton
	az.sampler_sphere(50,1,path + "/sample_SCT_beam",1);
	TCanvas *c = new TCanvas("ceta","eta1 Comparison",800,600);
	Data<double> oD(path + "/RFxi_out");
	Data<double> oD_off(path + "/../rrnd/xi_out");
	plot<double> pD(oD_off.col("#eta_{1}"),oD.col("#eta_{1}"),"#eta_{1} using different algorithms",
	"#eta_{1} RF off","#eta_{1} RF on");
	pD->SetMarkerStyle(20);
	pD->Draw("AP");
	c->SetGrid();
	*/

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	return 0;
}
