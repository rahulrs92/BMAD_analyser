#pragma once

#include<array>
#include<vector>
#include<bash_tools.hh>

using namespace std;

class acquisition
{
protected:
	static uint constexpr Vec = 0;	//vec
	static uint constexpr Mom = 1;	//mom
	static uint constexpr Spin = 2;	//spin

	
	template<uint> void get(uint,uint,vector<double>&,vector<double>&,vector<double>&);
	
public:
	static string source, destination;
	string login = "/opt/homebrew/bin/sshpass -p1qa12wqs ssh -i ~/key_shankar shankar@jedi01.ikp.kfa-juelich.de";
	string dlogin = "/opt/homebrew/bin/sshpass -p1qa12wqs scp -i ~/key_shankar -r shankar@jedi01.ikp.kfa-juelich.de";
	
	uint th, NL, Nc, Nc0, Nr, Ncr;
	bool f_zero;
	double L, Qx, Qy;					//Ring-length, Betatron Tunes
	acquisition(uint ith = 0) : th(ith), f_zero(0)
	{
		string s_NL = "cat \\$(find " + source + " -type f -name 'spin.dat' | sed -n \"1p\") | wc -l"; // "cat " + path + "/r1/c1/spin.dat | wc -l";
		string s_Nc = "find ./" + source + "/r1 -type f -name 'spin.dat' | wc -l";
		string s_Nc0 = "find ./" + source + "/r1/c*_0 -type f -name 'spin.dat' | wc -l";
		string s_Nr = "find " + source + " -type d -name 'r*' | wc -l";
		string s_Ncr = "find ./" + source + " -type f -name 'spin.dat' | wc -l";

		auto ss = bash::to_ss(login + " \"" + s_NL + "; " + s_Nc + "; " + s_Nc0 + "; " + s_Nr + "; " + s_Ncr + "\"");
		ss>>NL>>Nc>>Nc0>>Nr>>Ncr;
		if(Ncr == Nc * Nr) cout<<"repository consistent!"<<endl;
		else cerr<<"Warning: repository NOT consistent!"<<endl;
		Nc = Nc - Nc0;

		cout<<"NL = "<<NL<<endl;
		cout<<"Nc = "<<Nc<<endl;
		cout<<"Nc0 = "<<Nc<<endl;
		cout<<"Nr = "<<Nr<<endl;
		cout<<"Ncr = "<<Ncr<<endl;

		auto inf = info();
		L = inf[0]; Qx = inf[1]; Qy = inf[2];
		cout<<"L = "<<L<<endl;
		cout<<"Qx = "<<Qx<<endl;
		cout<<"Qy = "<<Qy<<endl;
	}

	acquisition(array<int,3> iNs,array<double,3> iinfo,uint ith = 0) : th(ith), f_zero(0)
	{
		NL = iNs[0]; Nc = iNs[1]; Nr = iNs[2]; Ncr = Nc * Nr;
		L = iinfo[0]; Qx = iinfo[1]; Qy = iinfo[2];
	}

	acquisition(const acquisition& ia, uint ith) : th(ith), f_zero(0)
	{
		NL = ia.NL; Nc = ia.Nc; Nc0 = ia.Nc0; Nr = ia.Nr; Ncr = ia.Ncr;
		L = ia.L; Qx = ia.Qx; Qy = ia.Qy;
	}

	//acquisition(const acquisition& ia) : acquisition(ia, ia.th) {}

	void set_zero_flag(bool ifz) { f_zero = ifz; }
	bool zero() { return f_zero; }

	void set_max_r(uint imr) { Nr = imr > Ncr/(Nc+Nc0) ? Nr : imr; }

	void n(vector<double>&);
	void vec(uint,uint,vector<double>&,vector<double>&,vector<double>&);
	void mom(uint,uint,vector<double>&,vector<double>&,vector<double>&);
	void spin(uint,uint,vector<double>&,vector<double>&,vector<double>&);
	void xi(vector<double>&,vector<double>&,vector<double>&,vector<double>&,vector<double>&);
	array<double,8> info(uint,uint);
	array<double,8> info();
};
