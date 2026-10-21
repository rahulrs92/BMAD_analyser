#pragma once

#include<lv1.hh>

using namespace ROOT;
using namespace std;

class lv2 : public lv1
{
public:

	int NS;
	bool f_tree;

	double nu_s0, sig0, eta0, m0;

	matrix<1,3> ep;
	matrix<3,2> U, H, W;
	matrix<3,3> U3, H3, W3, e_W3, M;
	matrix<3,1> u, h, w, e_w, psa, xi3, m;
	matrix<2,1> xi;
	//matrix<2,1> xi0;
	matrix<3,10> tensor_Xi, e_tensor_Xi;
	matrix<1,10> tensor_orbit;
	matrix<3,1,valarray<double>> tensor_xi3;

	matrix<1,2> c;

	string path;
	TFile* file;
	map<string,TObject*> mg;

	lv2() : NS(0), ep(0), W(0), u(0), h(0), xi(0), w(0), path("") {}
	lv2(int iNS) : lv2() { NS = iNS; }
	lv2(string ip,int iNS) : lv2(iNS) { path = ip; f_tree = 0; }

	lv2(string ip,TFile* ifile) : lv2()
	{
		path = ip; file = ifile; f_tree = 1;
		auto xitree = file->Get<TTree>("xi");
		NS = xitree->GetEntries();
	}

	tuple<double,matrix<3,1>> tune_vs_emittance(string,bool = 1);
	matrix<3,1> tune_vs_emittance_RF(string,bool = 1);
	tuple<double,matrix<1,2>> tune_vs_chromaticity(string,bool = 1);
	void tune_vs_sextupoles(string,bool = 1);

	matrix<3,1,valarray<double>> sampler(int,int,array<double,3>,double,string,bool = 1);
	matrix<3,1,valarray<double>> sampler_sphere(int,double,matrix<3,1>,string,bool = 1);
	matrix<3,1,valarray<double>> sampler_collapse(int,array<double,3>,double,string,bool = 1);
	//static void windmap(string,bool = 1);
	//static void spin_model(string,bool = 1);
	tuple<double,matrix<3,1>,matrix<3,2>,matrix<2,1>> spin_model_2(Data<double>&,string,bool = 1);
	tuple<double,matrix<3,1>,matrix<3,2>> travel_time_model(Data<double>&,string,bool = 1);
	tuple<double,matrix<3,1>,matrix<3,2>> spin_model_RFoff(Data<double>&,string,bool = 1);
	tuple<matrix<3,1>,matrix<3,2>,matrix<2,1>> spin_model_2R(Data<double>&,string,bool = 1);
	tuple<matrix<3,1>,matrix<3,3>,matrix<3,1>> spin_model_tensor(Data<double>&,string,bool = 1);
	tuple<matrix<3,1>,matrix<3,3>,matrix<3,1>> spin_model_3(Data<double>&,string,bool = 1);
	tuple<double,matrix<3,1>,matrix<3,3>> optimize_beam(Data<double>&,string,bool = 1);
	//tuple<matrix<3,1>,matrix<3,3>,matrix<3,1>> spin_model_3R_old(Data<double>&,string,bool = 1);
	matrix<1,10> spin_orbit_correction(Data<double>&,string,bool = 1);
	array<double,3> optimize(Data<double>&,Data<double>&,string,bool = 1);
	
	void correlation(string,bool = 1);
};


