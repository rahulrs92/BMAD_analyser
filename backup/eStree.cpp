#include<thread>

#include<matrix2.hh>
#include<bash_tools.hh>

#include<root/plot_tools.hh>
#include<root/TApplication.h>
#include<root/TCanvas.h>
#include<root/TRootCanvas.h>
#include<root/TTree.h>
#include<root/TFile.h>

using namespace std;

string ssh_login = "sshpass -p1qa12wqs ssh -i ~/key_shankar shankar@jedi01.ikp.kfa-juelich.de";
int nth = 8;
int NL, Nc, Nr, Ncr, r, c;
vector<uint> L, *l;
vector<double> *x, *y, *z, *px, *py, *pz, *sx, *sy, *sz;
vector<vector<vector<double>>> vx, vy, vz,
															vpx, vpy, vpz,
															vsx, vsy, vsz;
TTree* T;

void fill_r(int ir)
{
	stringstream loop;
	loop<<"for i in \\$(seq "<<Nc<<");"<<endl;
	loop<<"do"<<endl;
	loop<<"cat PTEDM/source_prototype/eS/r"<<ir<<"/c\\$i/vec.dat";
	loop<<" PTEDM/source_prototype/eS/r"<<ir<<"/c\\$i/mom.dat";
	loop<<" PTEDM/source_prototype/eS/r"<<ir<<"/c\\$i/spin.dat"<<endl;
	loop<<"done"<<endl;
	cout<<"Accessing the contents of PTEDM/source_prototype/eS/r"<<ir<<endl;
	stringstream ss = bash::to_ss(ssh_login + " \"" + loop.str() + "\"",to_string(ir));
	for(int i=0; i<Nc; i++)
	{
		for(int j=0; j<NL; j++) ss>>vx[(ir - 1)%nth][i][j]>>vy[(ir - 1)%nth][i][j]>>vz[(ir - 1)%nth][i][j];
		for(int j=0; j<NL; j++) ss>>vpx[(ir - 1)%nth][i][j]>>vpy[(ir - 1)%nth][i][j]>>vpz[(ir - 1)%nth][i][j];
		for(int j=0; j<NL; j++) ss>>vsx[(ir - 1)%nth][i][j]>>vsy[(ir - 1)%nth][i][j]>>vsz[(ir - 1)%nth][i][j];
	}
}

void eStree()
{
	string s_NL = "cat PTEDM/source_prototype/eS/r1/c1/spin.dat | wc -l";
	string s_Nc = "find ./PTEDM/source_prototype/eS/r1 -type f -name 'spin.dat' | wc -l";
	string s_Nr = "find ./PTEDM/source_prototype/eS/r[0-9]*/c1 -type f -name 'spin.dat' | wc -l";
	string s_Ncr = "find ./PTEDM/source_prototype/eS -type f -name 'spin.dat' | wc -l";

	auto ss = bash::to_ss(ssh_login + " \"" + s_NL + "; " + s_Nc + "; " + s_Nr + "; " + s_Ncr + "\"");
	ss>>NL>>Nc>>Nr>>Ncr;
	if(Ncr == Nc * Nr) cout<<"repository consistent!"<<endl;
	else cout<<"repository NOT consistent!"<<endl;

	cout<<"NL = "<<NL<<endl;
	cout<<"Nc = "<<Nc<<endl;
	cout<<"Nr = "<<Nr<<endl;
	cout<<"Ncr = "<<Ncr<<endl;

	TFile *hfile = TFile::Open("/home/Rahul/JEDI/PTEDM/1.0/cluster/cluster_data/eS/eS.root","RECREATE");//,"UPDATE");//

	T = new TTree("RF_off","RF_off");

	L = vector<uint>(NL);
	vx = vector<vector<vector<double>>>(nth,vector<vector<double>>(Nc,vector<double>(NL,0)));
	vy = vector<vector<vector<double>>>(nth,vector<vector<double>>(Nc,vector<double>(NL,0)));
	vz = vector<vector<vector<double>>>(nth,vector<vector<double>>(Nc,vector<double>(NL,0)));
	vpx = vector<vector<vector<double>>>(nth,vector<vector<double>>(Nc,vector<double>(NL,0)));
	vpy = vector<vector<vector<double>>>(nth,vector<vector<double>>(Nc,vector<double>(NL,0)));
	vpz = vector<vector<vector<double>>>(nth,vector<vector<double>>(Nc,vector<double>(NL,0)));
	vsx = vector<vector<vector<double>>>(nth,vector<vector<double>>(Nc,vector<double>(NL,0)));
	vsy = vector<vector<vector<double>>>(nth,vector<vector<double>>(Nc,vector<double>(NL,0)));
	vsz = vector<vector<vector<double>>>(nth,vector<vector<double>>(Nc,vector<double>(NL,0)));

	l = &L;
	//x = &vx[0]; y = &vy[0]; z = &vz[0];
	//px = &vpx[0]; py = &vpy[0]; pz = &vpy[0];
	//sx = &vsx[0]; sy = &vsy[0]; sz = &vsz[0];

	for(int j=0; j<NL; j++) L[j] = j;

	T->Branch("r",&r);
	T->Branch("c",&c);
	T->Branch("l",&l);
	T->Branch("x",&x);
	T->Branch("y",&y);
	T->Branch("z",&z);
	T->Branch("px",&px);
	T->Branch("py",&py);
	T->Branch("pz",&pz);
	T->Branch("sx",&sx);
	T->Branch("sy",&sy);
	T->Branch("sz",&sz);

	int nth = 8; r = 0;
	while(r+1 < Nr)
	{
		vector<thread> th;
		for(int i=0; i<nth; i++) { th.push_back(thread(fill_r,i+r+1)); }
		for(int i=0; i<nth; i++) th[i].join();
		for(int i=0; i<nth; i++)
		{
			r++;
			for(int j=0; j<Nc; j++)
			{
				c = j + 1;
				x = &vx[i][j]; y = &vy[i][j]; z = &vz[i][j];
				px = &vpx[i][j]; py = &vpy[i][j]; pz = &vpz[i][j];
				sx = &vsx[i][j]; sy = &vsy[i][j]; sz = &vsz[i][j];
				T->Fill();
			}
		}
	}

	//int oN = T->GetEntries();
	//for(int i=0; i<oN; i++) T->Show(i);
	T->Write();

	TTree* Tx = new TTree("xi","xi");
	int n = 0;
	double tf, td, ts, cx, cy;
	cout<<"yay!!!"<<endl;
	Tx->Branch("n",&n);
	Tx->Branch("tf",&tf);
	Tx->Branch("td",&td);
	Tx->Branch("ts",&ts);
	Tx->Branch("Cx",&cx);
	Tx->Branch("Cy",&cy);

	stringstream loop, ss1;
	loop<<"cat PTEDM/source_prototype/eS/xi";

	ss1 = bash::to_ss(ssh_login + " \"" + loop.str() + "\"");

	for(int i=0; i<Nr; i++)
	{
		n++;
		ss1>>tf>>td>>ts>>cx>>cy;
		cout<<n<<"\t"<<tf<<"\t"<<td<<"\t"<<ts<<"\t"<<cx<<"\t"<<cy<<endl;
		Tx->Fill();
	}
	Tx->Write();
	hfile->Close();
	/*
	int oN = T->GetEntries();
	cout<<"Number of entries in tree: "<<oN;

	vector<double> ov;
	vector<double> *ovr = &ov;
	int o_r, o_c;
	T->SetBranchAddress("r",&o_r);
	T->SetBranchAddress("c",&o_c);
	T->SetBranchAddress("sz",&ovr);

	for(int i=0; i<oN; i++)
	{
		T->GetEntry(i);
		cout<<"\nr = "<<o_r<<"; c = "<<o_c<<"\n"<<endl;
		for(int j=0; j<ov.size()-1; j++) cout<<ov[j]<<", ";
		cout<<ov[ov.size()-1]<<endl;
	}

	stringstream loop;
	ss<<"for i in \\$(seq "<<Nr<<");"<<endl;
	ss<<"do"<<endl;
	ss<<"for i in \\$(seq "<<Nc<<");"<<endl;
	ss<<"do"<<endl;
	*/
}

