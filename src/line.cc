#include<chrono>

//#include<TFile.h>
#include<TTree.h>

#include<proton.hh>

Data<double> proton_RF::line(string path, int N, int n0, bool pdf)
{
	fstream f, fe;
	f.open(path + "/line",fstream::in);
	if(f.is_open() && n0==0)
	{
		f.close();
		Data<double> oD(path + "/line");
		return oD;
	}
	f.close();
	f.open((path + "/line").data(),fstream::out | fstream::app);
	fe.open((path + "/line_error").data(),fstream::out | fstream::app);
	proton_RF *P;
	for(int i=n0; i<N; i++)
	{
		if(i>100) pdf = 0;
		proton::load_offset("/Users/utente/JEDI/PTEDM/1.0/source_COSY_mis/build/data/Deuteron/ref_sphere/c" + to_string(i+1));
		//proton::load_spin_offset("/Users/utente/JEDI/PTEDM/1.0/source_COSY_mis/build/data/Deuteron/ref_sphere/c" + to_string(i+1));
		P = new proton_RF(path + "/c" + to_string(i+1));
		P->apply_orbit_offset();
		//P.set_mode(mode);
		auto CSx = P->Courant_Snyder(0,pdf);
		auto CSx_e = P->e_CSx;
		auto CSy = P->Courant_Snyder(1,pdf);
		auto CSy_e = P->e_CSy;
		auto Qxm = P->betatron_tune(0,pdf);
		auto Qxm_e = P->e_Qxm;
		auto Qym = P->betatron_tune(1,pdf);
		auto Qym_e = P->e_Qym;
		auto syn = P->longitudinal(pdf);
		auto syn_e = P->e_syn;
		auto eta = P->phase_slip(pdf);
		auto eta_e = P->e_PSF;
		double nu_s = P->spin_tune(pdf);
		double nu_s_e = P->e_nu_s;
		//cout<<"T = "<<eta[0]<<endl;
		
		f<<P->Kt1()<<'\t'<<P->Kt2()<<'\t'<<P->Kt3()<<'\t'
			<<P->Cx()<<'\t'<<P->Cy()<<'\t'
			<<Qxm[1]<<'\t'<<Qym[1]<<'\t'<<eta[2]<<'\t'
			<<CSx[0]<<'\t'<<CSy[0]<<'\t'
			<<syn[3]<<'\t'<<syn[4]<<'\t'<<eta[1]<<'\t'
			<<nu_s<<endl;
		
		fe<<"0.0\t0.0\t0.0\t0.0\t0.0\t"
			<<Qxm_e[1]<<'\t'<<Qym_e[1]<<'\t'<<eta_e[2]<<'\t'
			<<CSx_e[0]<<'\t'<<CSy_e[0]<<'\t'
			<<syn_e[3]<<'\t'<<syn_e[4]<<'\t'<<eta_e[1]<<'\t'
			<<nu_s_e<<endl;
		
		delete P;
	}
	f.close();
	fe.close();
	Data<double> oD(path + "/line");
	oD.label(1) = "t_{f}";
	oD.label(2) = "t_{d}";
	oD.label(3) = "t_{s}";
	oD.label(4) = "#xi_{x}";
	oD.label(5) = "#xi_{y}";
	oD.label(6) = "#xi_{xm}";
	oD.label(7) = "#xi_{ym}";
	oD.label(8) = "#eta_{1}";
	oD.label(9) = "#epsilon_{x}";
	oD.label(10) = "#epsilon_{y}";
	oD.label(11) = "#delta_{a}";
	oD.label(12) = "#omega_{s}";
	oD.label(13) = "#eta_{0}";
	oD.label(14) = "#nu_{s}";

	f.open((path + "/line_info").data(),fstream::out | fstream::trunc);
	fe.open((path + "/line_error_info").data(),fstream::out | fstream::trunc);
	for(int i=0; i<oD.order()[1]; i++) f<<oD.label(i+1)<<endl;
	for(int i=0; i<oD.order()[1]; i++) fe<<"#Delta"<<oD.label(i+1)<<endl;
	f.close();
	fe.close();

	return oD;
}

Data<double> proton::line(string path, int N, bool pdf)
{
	fstream f;
	f.open(path + "/line",fstream::in);
	if(f.is_open())
	{
		f.close();
		Data<double> oD(path + "/line");
		return oD;
	}
	f.close();
	f.open((path + "/line").data(),fstream::out | fstream::trunc);
	for(int i=0; i<N; i++)
	{
		proton P(path + "/c" + to_string(i+1));
		//P.set_mode(mode);
		auto CSx = P.Courant_Snyder(0,pdf);
		auto CSy = P.Courant_Snyder(1,pdf);

		auto eta = P.phase_slip(pdf);
		//cout<<"Phase slip per turn = "<<eta[1]+eta[0]<<endl;
		
		double nu_s = P.spin_tune(pdf);
		
		f<<P.Cx()<<'\t'<<P.Cy()<<'\t'
		<<P.Kt1()<<'\t'<<P.Kt2()<<'\t'<<P.Kt3()<<'\t'
		<<CSx[0]<<'\t'<<CSy[0]<<'\t'<<P.pz[0]<<'\t'
		//<<eta[0]<<'\t'<<eta[1]<<'\t'<<nu_s<<endl;
		<<eta[1]<<'\t'<<nu_s<<'\t'
		<<CSx[6]<<'\t'<<CSy[6]<<endl;

		cout<<"T at line = "<<eta[1]<<endl;
	}
	f.close();
	Data<double> oD(path + "/line");
	oD.label(1) = "#xi_{x}";
	oD.label(2) = "#xi_{y}";
	oD.label(3) = "t_{f}";
	oD.label(4) = "t_{d}";
	oD.label(5) = "t_{s}";
	oD.label(6) = "#epsilon_{x}";
	oD.label(7) = "#epsilon_{y}";
	oD.label(8) = "#delta_{0}";
	//oD.label(6) = "T_{model}";
	oD.label(9) = "T";
	oD.label(10) = "#nu_{s}";
	oD.label(11) = "#epsilon_{x0}";
	oD.label(12) = "#epsilon_{y0}";

	f.open((path + "/line_info").data(),fstream::out | fstream::trunc);
	for(int i=0; i<oD.order()[1]; i++) f<<oD.label(i+1)<<endl;
	f.close();

	return oD;
}


