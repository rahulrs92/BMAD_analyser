#include<lv1.hh>
#include<proton.hh>

using namespace std;

/*
array<vector<double>,9> proton_RF::orbit_offset =
{
	vector<double>(),vector<double>(),vector<double>(),		//x,y,z
	vector<double>(),vector<double>(),vector<double>(),		//px,py,pz
	vector<double>(),vector<double>(),vector<double>()		//sx,sy,sz
};
*/
array<double,6> proton_RF::orbit_offset =
{
	0.0,0.0,0.0,		//x,y,z
	0.0,0.0,0.0,		//px,py,pz
};

map<string,double> proton_RF::m0 =
{
	{"proton",938.27208943},
	{"deuteron",1875.612928}
};

map<string,double> proton_RF::p =
{
	{"PSR",294.057},
	{"PSR-electric",294.057},
	{"PSR-magnetic",294.057},
	{"COSY",970.0},
	{"COSY_mis",970.0},
	{"Hybrid",700.74037}
};

map<string,double> proton_RF::L =
{
	{"PSR",123.358},
	{"PSR-electric",123.358},
	{"PSR-magnetic",129.758},
	{"COSY",183.473},
	{"COSY_mis",183.473},
	{"Hybrid",801.6}
};

map<string,map<string,array<double,3>>> proton_RF::nu0 = 
{
	{
		"proton",
		{
			{"PSR",{0.0, 0.0, 0.0}},
			{"PSR-electric",{0.143855, -2.30725e-7, -3.46717e-13}},
			{"PSR-magnetic",{-0.121167, 2.33025e-7, 2.79126e-14}},
			{"COSY",{0.0786942, -4.70913e-8, -2.92657e-15}},
			{"COSY_mis",{0.0786942, -4.70913e-8, -2.92657e-15}},
			{"Hybrid",{0.0, 0.0, 0.0}}
		}
	},
	{
		"deuteron",
		{
			{"PSR",{0.0, 0.0, 0.0}},
			{"PSR-electric",{0.143855, -2.30725e-7, -3.46717e-13}},
			{"PSR-magnetic",{-0.121167, 2.33025e-7, 2.79126e-14}},
			{"COSY",{-0.160977, -1.78822e-07, -4.93645e-13}},
			{"COSY_mis",{-0.161656, 3.23282e-08, -6.15587e-15}},
			//{"COSY_mis",{0.0, 0.0, 0.0}}, // already applied orbit offset
			{"Hybrid",{0.0, 0.0, 0.0}}
		}
	}
};

/*
map<string,TObject*> lv1::mg =
{
	{"pt",new TList()},
	{"eta_1",new TMultiGraph()},
	{"eta_2",new TMultiGraph()}
};
*/
