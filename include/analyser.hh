#pragma once

#include<lv2.hh>

using namespace ROOT;
using namespace std;

class analyser : public lv2
{
public:

	int NS;

	double sig0, eta0;

	matrix<1,3> ep;
	matrix<3,3> Xi, V, T, M;
	matrix<3,1> a, psa, xi, xi0;

	string path;
	map<string,TObject*> mg;

	analyser() : NS(0), ep(0), Xi(0), M(0), a(0), xi(0), xi0(0), path("") {}
	analyser(int iNS) : analyser() { NS = iNS; }
	analyser(string ip,int iNS) : analyser(iNS) { path = ip; }

	tuple<matrix<3,1>,matrix<3,3>> windmap(bool = 1);
	matrix<3,1,valarray<double>> sampler(bool = 1);
	//static void windmap(string,bool = 1);
	//static void spin_model(string,bool = 1);
	tuple<matrix<1,3>,matrix<3,1>,matrix<3,3>,matrix<3,1>> spin_model_2D(bool = 1);
	tuple<matrix<3,1>,matrix<3,3>,matrix<3,1>> spin_model(bool = 1);
	void spin_model_errors(bool = 1);
	void SM_RF_test(Data<double>&,bool = 1);
	void E_RF_test(Data<double>&,bool = 1);
	void SM_RF_test2(Data<double>&,bool = 1);
};


