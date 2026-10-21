#include<chrono>

//#include<TFile.h>
#include<TTree.h>

#include<proton.hh>

Data<double> proton_RF::beam(string path, int N, int n0, bool pdf)
{
	fstream f, fe;
	f.open(path + "/beam",fstream::in);
	if(f.is_open() && n0==0)
	{
		f.close();
		Data<double> oD(path + "/beam");
		return oD;
	}
	f.close();
	f.open((path + "/beam").data(),fstream::out | fstream::app);
	fe.open((path + "/beam_error").data(),fstream::out | fstream::app);
	proton_RF *P;
	for(int i=n0; i<N; i++)
	{
		if(i>100) pdf = 0;
		if(proton_RF::mode == "COSY_mis") proton::load_offset("/Users/utente/JEDI/PTEDM/1.0/source_COSY_mis/build/data/Deuteron/ref_sphere/c" + to_string(i+1));
		//proton::load_spin_offset("/Users/utente/JEDI/PTEDM/1.0/source_COSY_mis/build/data/Deuteron/ref_sphere/c" + to_string(i+1));
		// cout<<"YAY!! "<<path + "/c" + to_string(i+1)<<endl;
		P = new proton_RF(path + "/c" + to_string(i+1));
		if(proton_RF::mode == "COSY_mis") P->apply_orbit_offset();
		//P.set_mode(mode);
		double nu_s = P->spin_tune(pdf);
		double nu_s_e = P->e_nu_s;
		double SCT = P->SCT(pdf);
		double SCT_e = P->e_sct;
		//cout<<"T = "<<eta[0]<<endl;
		
		f<<P->Kt1()<<'\t'<<P->Kt2()<<'\t'<<P->Kt3()<<'\t'
			<<nu_s<<'\t'<<SCT<<endl;
		
		fe<<"0.0\t0.0\t0.0\t"
			<<nu_s_e<<'\t'<<SCT_e<<endl;
		
		delete P;
	}
	f.close();
	fe.close();
	Data<double> oD(path + "/beam");
	oD.label(1) = "t_{f}";
	oD.label(2) = "t_{d}";
	oD.label(3) = "t_{s}";
	oD.label(4) = "#nu_{s}";
	oD.label(5) = "#tau";

	f.open((path + "/beam_info").data(),fstream::out | fstream::trunc);
	fe.open((path + "/beam_error_info").data(),fstream::out | fstream::trunc);
	for(int i=0; i<oD.order()[1]; i++) f<<oD.label(i+1)<<endl;
	for(int i=0; i<oD.order()[1]; i++) fe<<"#Delta"<<oD.label(i+1)<<endl;
	f.close();
	fe.close();

	return oD;
}

