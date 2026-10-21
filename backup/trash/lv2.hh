#pragma once

#include<lv1.hh>

using namespace ROOT;
using namespace std;

class lv2 : public lv1
{
public:
	map<string,TMultiGraph*> mg;
	array<double,3> v_eta1, v_sigma1;
	bool f_eta1, f_sigma1;

	lv2() : f_eta1(0), f_sigma1(0) {}

	array<double,3> eta1(string,Data<double>&);
	array<double,3> sigma1(string,Data<double>&);
};

