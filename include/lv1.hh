#pragma once

#include<thread>

#include<acquisition.hh>
#include<proton.hh>

using namespace ROOT;
using namespace std;

class lv1
{
protected:
	vector<vector<vector<double>>> tI;
	vector<vector<double>> tO;

public:
	map<string,TObject*> mg;
	array<double,2> v_e, v_pt, v_c;
	array<double,3> v_eta, v_sigma;
	bool f_e, f_pt, f_c, f_eta, f_sigma;
	matrix<2,1> xi0;
	matrix<2,3> Xi;
	matrix<2,2> Xi_2D;
	matrix<3,1> xi03, xi03_e;
	matrix<3,3> Xi3, Xi3_e;
	matrix<12,12> Cov_Xi;

	matrix<3, 1, std::valarray<double>> tensor_solutions;

	lv1() : f_e(0), f_pt(0), f_c(0), f_eta(0), f_sigma(0) {}

	array<double,2> e(string,Data<double>&);
	array<double,2> pt(string,Data<double>&,bool = 1);
	array<double,2> c(string,Data<double>&,array<double,2> = {0,0});
	array<double,5> eta(string,Data<double>&,array<double,2> = {0,0},bool = 1);
	array<double,3> sigma(string,Data<double>&,array<double,2> = {0,0},array<double,2> = {0,0});

	tuple<matrix<2,1>,matrix<2,3>> windmap_plane(Data<double>&,string,bool = 1);
	tuple<matrix<3,1>,matrix<3,3>> windmap_tensor(Data<double>&,string,bool = 1);
	tuple<matrix<3,1>,matrix<3,3>> windmap(Data<double>&,string,bool = 1);
	tuple<matrix<3,1>,matrix<3,3>> windmap(Data<double>&,valarray<double>,string,bool = 1);
	tuple<matrix<2,1>,matrix<2,2>> windmap_2D(Data<double>&,string,bool = 1);

	matrix<3,3> COSY_tune_map(Data<double>&,matrix<3,1>,matrix<3,1>,string,bool = 1);

	TTree* grid_RF(double, double, TFile*,bool = 1);
	TTree* grid(double, double, TFile*,bool = 1);

	array<double,7> point(acquisition&,uint,uint,double = 0.0,double = 0.0,bool = 1);
	array<double,11> point_RF(acquisition&,uint,uint,bool = 1);
	array<double,11> point_RF(string,uint,uint,bool = 1);
	array<double,7> point(string,uint,uint,bool = 1);
	array<double,5> xi(string,uint,uint,bool = 1);
	void grid(acquisition&,string,bool = 1);
	void grid(string,bool = 1);
	void grid_xi(string,bool = 1);
	void multiplane(acquisition&,string,uint = 0,bool = 0);
	void line(string,uint = 0,bool = 0);

	matrix<3,1,valarray<double>>& solve_tensor_system(matrix<3,10>,matrix<3,1>);

	void SCT_merge(acquisition&,string);
};
