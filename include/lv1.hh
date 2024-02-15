#pragma once

#include<proton.hh>

using namespace ROOT;
using namespace std;

class lv1
{
public:
	map<string,TObject*> mg;
	array<double,2> v_e, v_pt, v_c;
	array<double,3> v_eta, v_sigma;
	bool f_e, f_pt, f_c, f_eta, f_sigma;

	lv1() : f_e(0), f_pt(0), f_c(0), f_eta(0), f_sigma(0) {}

	array<double,2> e(string,Data<double>&);
	array<double,2> pt(string,Data<double>&,bool = 1);
	array<double,2> c(string,Data<double>&,array<double,2> = {0,0});
	array<double,3> eta(string,Data<double>&,array<double,2> = {0,0},bool = 1);
	array<double,3> sigma(string,Data<double>&,array<double,2> = {0,0},array<double,2> = {0,0});
};
