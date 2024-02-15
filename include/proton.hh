#pragma once

#include<matrix2.hh>
#include<data_tools.hh>
#include<root/plot_tools.hh>

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

	array<double,6> Courant_Snyder(bool = 0, bool = 1);

	array<double,2> phase_slip(bool Out = 1);

	double spin_tune(bool = 1);

	static Data<double> line(string, int, bool = 1);
	static void phase_slip_comparison(string, Data<double>&);
	static array<double,3> spin_tune_variation(string, Data<double>&);
	static array<double,3> get_etas(string, Data<double>&);
	static array<double,3> get_sigmas(string, Data<double>&);

};




