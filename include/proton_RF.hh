#pragma once

#include<matrix2.hh>
#include<data_tools.hh>
#include<root/plot_tools.hh>

using namespace ROOT;
using namespace std;


class proton_RF
{
protected:
	string path;
	static double constexpr pi = TMath::Pi(), mp = 938.27208816, G = (5.5856946893 - 2) / 2;

	double z0 = 1.19376;
	double z1 = -1.17916;
	double r = 0.7147, p = 294.057, L, E, gm, bt;

	Data<double> D;
	Data<string> info;

	array<double,6> CSx, CSy;
	array<double,5> syn;
	array<double,3> PSF;
	array<double,2> sig, zet;
	double Tps, nu_s;
	//pair<array<double,2>,valarray<double>> wm_pz;			//wave mean, amplitude and zeroes, useful for timestamping
	bool f_CSx, f_CSy, f_syn, f_PSF, f_ST, f_sig, f_zet;

	vector<double> trash;
	stringstream ss;

	proton_RF(Data<double> iD) : D(iD), n(D[0]), x(D[1]), y(D[2]), z(D[3]),
															px(D[4]), py(D[5]), pz(D[6]),
															sx(D[7]), sy(D[8]), sz(D[9])
	{
		E = pow(pow(p,2) + pow(mp,2),0.5);
		gm = E/ mp;
		bt = p/E;

		N = n.size();
		CSx = {0,0,0,0,0,0}; f_CSx = 0;
		CSy = {0,0,0,0,0,0}; f_CSy = 0;
		syn = {0,0,0,0,0}; f_syn = 0;
		PSF = {0,0,0}; f_PSF = 0;
		sig = {0,0}; f_sig = 0;
		zet = {0,0}; f_zet = 0;
		Tps = 0; nu_s = 0; f_ST = 0;

	}
public:
	long N;
	double md;
	vector<double> &n, &x, &y, &z, &px, &py, &pz, &sx, &sy, &sz;

	map<string,TMultiGraph*> mg;

	double Lr() { return stod(info[2][1]); }
	double Qx() { return stod(info[2][2]); }
	double Qy() { return stod(info[2][3]); }
	double Cx() { return stod(info[2][4]); }
	double Cy() { return stod(info[2][5]); }

	proton_RF(string ipath) : proton_RF(Data<double>(ipath + "/vec.dat") +
																Data<double>(ipath + "/mom.dat") +
																Data<double>(ipath + "/spin.dat"))
	{
		path = ipath;
		info = Data<string>(ipath + "/info.dat");
		D.title() = "Phase Space and spin coordinates @ Q_{X} = " + info[2][2] + ", Q_{Y} = " + info[2][3];
		D.labels() = {"turn no.","x","y","z","p_{x}","p_{y}","#delta","s_{x}","s_{y}","s_{z}"};
		L = Lr();
	}

	proton_RF(const char* ipath) : proton_RF(string(ipath)) {}

	double Tp_model()
	{
		if(!f_CSx) Courant_Snyder(0,0);
		if(!f_CSy) Courant_Snyder(1,0);

		Tps = - (pi/L) * (CSx[0]*Cx() + CSy[0]*Cy());
		return Tps;
	}

	static pair<array<double,2>,valarray<double>> wave_mean(vector<double>&, vector<double>&);

	array<double,6> Courant_Snyder(bool = 0, bool = 1);

	array<double,5> longitudinal(bool = 1);

	array<double,3> phase_slip(bool = 1);

	double spin_tune(bool = 1);
	array<double,2> spin_sigma(bool = 1);
	array<double,2> spin_zeta(bool = 1);

	static Data<double> line(string, int, bool = 1);
};

