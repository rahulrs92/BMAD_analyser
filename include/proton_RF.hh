#pragma once

#include<TFile.h>
#include<TTree.h>

#include<matrix3.hh>
#include<data_tools.hh>
#include<root/plot_tools.hh>

using namespace ROOT;
using namespace std;

class proton;

class proton_RF
{
protected:
	string path;

	static map<string,map<string,array<double,3>>> nu0;		//reference spin tune up to three levels of precision for each (mode,ID) pair
	static map<string,double> L;							//reference ring circumference for each (mode,ID) pair
	static map<string,double> p;							//reference momentum for each (mode,ID) pair
	//static array<vector<double>,9> orbit_offset;
	static array<double,6> orbit_offset;
	
	double E, gm, bt;										//total energy, relativistic gamma, beta
	double qx, qy, cx, cy;									//theoretical betatron tunes, chromaticities
	double kt1 = 0, kt2 = 0, kt3 = 0;		 				//sextupole strengths
	
	Data<double> D;
	Data<string> info;
	
	array<double,7> CSx, e_CSx, CSy, e_CSy; 	//Emittance, Alpha, Beta, Gamma, x0, px0, zero-emittance
	array<double,5> syn, e_syn;					//z-mean, z-amp, pz-mean, pz-amp, omega
	array<double,3> PSF, e_PSF;					//Phase Slip Factor: Theoretical (dT/T)_T, eta0, eta1
	array<double,2> Qxm, Qym, e_Qxm, e_Qym;		//measured betatron tunes, chromaticities and errors
	
	array<double,2> sig, zet;
	double Tps, nu_s, e_nu_s, sct, e_sct;
	//pair<array<double,2>,valarray<double>> wm_pz;			//wave mean, amplitude and zeroes, useful for timestamping
	bool f_CSx, f_CSy, f_syn, f_PSF, f_ST, f_BTx, f_BTy, f_SCT, f_sig, f_zet, f_tree;
	
	stringstream ss;
	
	proton_RF(Data<double> iD) : D(iD), n(D[0]), x(D[1]), y(D[2]), z(D[3]),
	px(D[4]), py(D[5]), pz(D[6]),
	sx(D[7]), sy(D[8]), sz(D[9])
	{
		E = pow(pow(p[mode],2) + pow(m0[ID],2),0.5);
		gm = E/ m0[ID];
		bt = p[mode]/E;
		
		N = n.size();
		CSx = {0,0,0,0,0,0,0}; f_CSx = 0;
		CSy = {0,0,0,0,0,0,0}; f_CSy = 0;
		syn = {0,0,0,0,0}; f_syn = 0;
		PSF = {0,0,0}; f_PSF = 0;
		sig = {0,0}; f_sig = 0;
		zet = {0,0}; f_zet = 0;
		Qxm = {0,0}; f_BTx = 0;
		Qym = {0,0}; f_BTy = 0;
		Tps = 0; nu_s = 0; f_ST = 0;
		sct = 0; f_SCT = 0;
		
		f_tree = 0;

		//apply_orbit_offset();
	}

	/*
	void apply_orbit_offset()
	{
		if(orbit_offset[0].size() != N) { cerr<<"Orbit offset size mismatch! Offset not applied."<<endl; return; }
		
		for(int i=0; i<orbit_offset[0].size(); i++)
		{
			x[i] -= orbit_offset[0][i];
			y[i] -= orbit_offset[1][i];
			z[i] -= orbit_offset[2][i];
			px[i] -= orbit_offset[3][i];
			py[i] -= orbit_offset[4][i];
			pz[i] -= orbit_offset[5][i];
		}
		vector<double> sx1(N,0.0), sy1(N,0.0), sz1(N,1.0);
		double ct = 1, st = 0;
		matrix<3,1> u;
		// Adjust spin vectors: absolute change in spin direction due to emittance
		for(int i=1; i<orbit_offset[6].size(); i++)
		{
			matrix<3,1> r = {orbit_offset[6][i], orbit_offset[7][i], orbit_offset[8][i]};
			matrix<3,1> s = {sx[i], sy[i], sz[i]};
			matrix<3,1> s0 = {sx[0], sy[0], sz[0]};
			
			u = r * s;
			u = u / (u.norm());
			// if(i<4) cout<<"r = "<<r<<endl;
			// if(i<4) cout<<"s = "<<s<<endl;
			// if(i<4) cout<<"u = "<<u<<endl;
			ct = r.tp()*s; st = pow(1 - pow(ct>1 ? 1 : ct<-1 ? -1 : ct,2),0.5);
			// if(i<4) cout<<"ct = "<<ct<<", st = "<<st<<endl;
			auto s1 = ct*s0 + st*(u * s0) + (1 - ct)*(u.tp()*s0)*u;
			sx1[i] = s1[0]; sy1[i] = s1[1]; sz1[i] = s1[2];
		}
		
		// Adjust spin vectors: only horizontal change in spin direction due to emittance
		
		// for(int i=1; i<orbit_offset[6].size(); i++)
		// {
		// 	matrix<3,1> r1 = {orbit_offset[6][i], orbit_offset[7][i], orbit_offset[8][i]};
		// 	// matrix<3,1> r0 = {orbit_offset[6][i-1], orbit_offset[7][i-1], orbit_offset[8][i-1]};
		// 	// double tt1 = (r1[0] - r0[0]) / (r1[2] - r0[2]), tt2 = (sx[i] - sx[i-1]) / (sz[i] - sz[i-1]);
		// 	double tt1 = r1[0] / r1[2], tt2 = sx[i] / sz[i];
		// 	double dt = atan((tt2 - tt1) / (1 + tt1*tt2));
		// 	double ct = cos(dt), st = sin(dt);
		// 	sx1[i] = ct*sx1[i-1] + st*sz1[i-1];
		// 	sy1[i] = sy1[i-1];
		// 	sz1[i] = - st*sx1[i-1] + ct*sz1[i-1];
		// }
		
		sx = sx1; sy = sy1; sz = sz1;
		// for(int i=0; i<50; i++) cout<<sx[i]<<"\t"<<sy[i]<<"\t"<<sz[i]<<endl;
	}
	*/
public:
	static string ID;
	static string mode;

	static map<string,double> m0;

	static double constexpr pi = TMath::Pi(), G = (5.5856946893 - 2) / 2, cel = 299792458.0;

	long N;
	double md;
	vector<double> &n, &x, &y, &z, &px, &py, &pz, &sx, &sy, &sz;

	map<string,TMultiGraph*> mg;

	/*
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
			//cerr<<"Offset files not found in " + filepath+", orbit offset not applied!"<<endl;
			return;
		}
		Data<double> offset_data = Data<double>(filepath + "/vec.dat") + Data<double>(filepath + "/mom.dat") + Data<double>(filepath + "/spin.dat");
		for(int i=0; i<9; ++i)  orbit_offset[i] = offset_data[i+1];
	}
	*/

	double Lr() { return L[mode]; }
	double Qx() { return qx; }
	double Qy() { return qy; }
	double Cx() { return cx; }
	double Cy() { return cy; }
	double Kt1() { return kt1; }
	double Kt2() { return kt2; }
	double Kt3() { return kt3; }
	
	double Cx_m() { if(!f_BTx) betatron_tune(0,0); return Qxm[1]; }
	double Cy_m() { if(!f_BTy) betatron_tune(1,0); return Qym[1]; }

	proton_RF(double iqx,double iqy,double icx,double icy,vector<double> &in,
						vector<double> &ix,vector<double> &iy,vector<double> &iz,
						vector<double> &ipx,vector<double> &ipy,vector<double> &ipz,
						vector<double> &isx,vector<double> &isy,vector<double> &isz) :
															n(in), x(ix), y(iy), z(iz),
															px(ipx), py(ipy), pz(ipz),
															sx(isx), sy(isy), sz(isz)
	{
		E = pow(pow(p[mode],2) + pow(m0[ID],2),0.5);
		gm = E/ m0[ID];
		bt = p[mode]/E;
		
		N = n.size();
		CSx = {0,0,0,0,0,0,0}; f_CSx = 0;
		CSy = {0,0,0,0,0,0,0}; f_CSy = 0;
		syn = {0,0,0,0,0}; f_syn = 0;
		PSF = {0,0,0}; f_PSF = 0;
		sig = {0,0}; f_sig = 0;
		zet = {0,0}; f_zet = 0;
		Qxm = {0,0}; f_BTx = 0;
		Qym = {0,0}; f_BTy = 0;
		Tps = 0; nu_s = 0; f_ST = 0;
		sct = 0; f_SCT = 0;

		qx = iqx; qy = iqy; cx = icx; cy = icy;
		f_tree = 1;

		//apply_orbit_offset();
	}

	proton_RF(string ipath) : proton_RF(Data<double>(ipath + "/vec.dat") +
																Data<double>(ipath + "/mom.dat") +
																Data<double>(ipath + "/spin.dat"))
	{
		path = ipath;
		info = Data<string>(ipath + "/info.dat");
		// cout<<"YAY!! "<<endl;
		D.title() = "Phase Space and spin coordinates @ Q_{X} = " + info[2][2] + ", Q_{Y} = " + info[2][3];
		D.labels() = {"turn no.","x","y","z","p_{x}","p_{y}","#delta","s_{x}","s_{y}","s_{z}"};
		L[mode] = stod(info[2][1]);
		qx = stod(info[2][2]);
		qy = stod(info[2][3]);
		cx = stod(info[2][4]);
		cy = stod(info[2][5]);
		if(info.order()[0] > 5) kt1 = stod(info[2][6]); else cerr<<"kt1 not found; set to zero!"<<endl;
		if(info.order()[0] > 6) kt2 = stod(info[2][7]); else cerr<<"kt2 not found; set to zero!"<<endl;
		if(info.order()[0] > 7) kt3 = stod(info[2][8]); else cerr<<"kt3 not found; set to zero!"<<endl;
	}

	void apply_orbit_offset()
	{
		for(int i=0; i<N; i++)
		{
			x[i] -= orbit_offset[0];
			y[i] -= orbit_offset[1];
			z[i] -= orbit_offset[2];
			px[i] -= orbit_offset[3];
			py[i] -= orbit_offset[4];
			pz[i] -= orbit_offset[5];
		}
	}
	
	proton_RF(const char* ipath) : proton_RF(string(ipath)) {}

	~proton_RF() { mg.clear(); }

	array<double,5> xi() { return {cx,cy,kt1,kt2,kt3}; }

	static void set_mode(string imode) { mode = imode; }

	void reset() { f_CSx=0; f_CSy=0; f_syn=0; f_PSF=0; f_ST=0; f_sig=0; f_zet=0; }

	void setPath(string ipath) { path = ipath; }

	void setChromaticity(double icx,double icy) { cx = icx; cy = icy; }

	double Tp_model()
	{
		if(!f_CSx) Courant_Snyder(0,0);
		if(!f_CSy) Courant_Snyder(1,0);

		Tps = - (pi/L[mode]) * ((isnan(CSx[0]) ? 0 : CSx[0])*Cx() + (isnan(CSy[0]) ? 0 : CSy[0])*Cy());
		return Tps;
	}

	static pair<array<double,2>,valarray<double>> wave_mean(vector<double>&, vector<double>&);
	static array<double,2> simple_SCT(string path);

	array<double,7> Courant_Snyder(bool = 0, bool = 1);
	matrix<6,5> Courant_Snyder_Jacobian(matrix<5,1> M);
	array<double,7> Courant_Snyder_errors(bool u = 0, bool out = 0)
	{
		if(!(u ? f_CSy : f_CSx)) Courant_Snyder(u,out);
		return (u ? e_CSy : e_CSx);
	}
	
	array<double,5> longitudinal(bool = 1);
	matrix<5,5> longitudinal_Jacobian(matrix<5,1> M);
	array<double,5> longitudinal_errors(bool out = 0)
	{
		if(!f_syn) longitudinal(out);
		return e_syn;
	}

	array<double,3> phase_slip(bool = 1);
	array<double,3> phase_slip_errors(bool out = 0)
	{
		if(!f_PSF) phase_slip(out);
		return e_PSF;
	}

	double spin_tune(bool = 1);
	double spin_tune_error(bool out = 0)
	{
		if(!f_ST) spin_tune(out);
		return e_nu_s;
	}

	array<double,2> betatron_tune(bool,bool = 1);
	array<double,2> betatron_tune_error(bool iy,bool out = 0)
	{
		if(!(iy ? f_BTy : f_BTx)) betatron_tune(iy,out);
		return iy ? e_Qym : e_Qxm;
	}

	array<double,2> spin_sigma(bool = 1);

	double SCT(bool = 1);

	static Data<double> line(string, int, int = 0, bool = 1);
	static Data<double> beam(string, int, int = 0, bool = 1);
};

