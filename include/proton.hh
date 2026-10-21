#pragma once

#include<proton_RF.hh>

using namespace ROOT;
using namespace std;

class proton : public proton_RF
{
protected:
	array<double,2> PSF;

public:
	proton(string ipath) : proton_RF(ipath) {}

	proton(const char* ipath) : proton(string(ipath)) {}

	proton(double iqx,double iqy,double icx,double icy,vector<double> &in,
						vector<double> &ix,vector<double> &iy,vector<double> &iz,
						vector<double> &ipx,vector<double> &ipy,vector<double> &ipz,
						vector<double> &isx,vector<double> &isy,vector<double> &isz) : proton_RF(iqx,iqy,icx,icy,in,
						ix,iy,iz,
						ipx,ipy,ipz,
						isx,isy,isz) {}

	~proton() { mg.clear(); }

	static void load_offset(string filepath)
	{
		if(mode != "COSY_mis") return; //only for misaligned COSY simulations
		try
		{
			fstream f(filepath + "/vec.dat", fstream::in);
			if(!f.is_open()) throw runtime_error("File not found");
			f.close();
		}
		catch(...)
		{
			cerr<<"Offset files not found in " + filepath + ", orbit offset not applied!"<<endl;
			return;
		}
		Data<double> offset_data = Data<double>(filepath + "/vec.dat") + Data<double>(filepath + "/mom.dat");
		for(int i=0; i<6; ++i)  orbit_offset[i] = valarray<double>(offset_data.col(i+1)).sum() / offset_data.order()[0];
	}
	static void load_spin_offset(string filepath)
	{
		if(mode != "COSY_mis") return; //only for misaligned COSY simulations
		try
		{
			fstream f(filepath + "/spin.dat", fstream::in);
			if(!f.is_open()) throw runtime_error("File not found");
			f.close();
		}
		catch(...)
		{
			cerr<<"Offset files not found in " + filepath + ", orbit offset not applied!"<<endl;
			return;
		}

		proton P(filepath);
		nu0[ID][mode] = {0.0, 0.0, 0.0}; // reset spin tune corrections due to misalignment
		nu0[ID][mode][0] = P.spin_tune(1); P.f_ST = 0;// set only the zero-order correction
		nu0[ID][mode][1] = P.spin_tune(1); P.f_ST = 0; // set only the first-order correction
		nu0[ID][mode][2] = P.spin_tune(1); P.f_ST = 0; // set only the second-order correction
	}

	array<double,7> Courant_Snyder(bool = 0, bool = 1);

	array<double,2> phase_slip(bool Out = 1);

	double spin_tune(bool = 1);

	static void phase_slip_comparison(string, Data<double>&);
	static array<double,3> spin_tune_variation(string, Data<double>&);
	static array<double,3> get_etas(string, Data<double>&);
	static array<double,3> get_sigmas(string, Data<double>&);

	static Data<double> line(string, int, bool = 1);
	TTree* grid(double, double, TFile*, bool = 1);
};




