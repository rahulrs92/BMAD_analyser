#include<acquisition.hh>

/*
template<uint u> void acquisition::get(uint ir,uint ic,
														vector<double>& x,vector<double>& y,vector<double>& z)
{
	string space = u==Vec ? "vec" : (u==Mom ? "mom" : (u==Spin ? "spin" : ""));
	cout<<"Acquiring: r"<<ir<<"/c"<<ic<<"/"<<space<<".dat"<<endl;
	if(space == "") { cerr<<"template argument not valid!"<<endl; return; }
	stringstream ss;
	ss<<login<<" \"cat "<<path<<"/r"<<ir<<"/c"<<ic<<(f_zero ? "_0/" : "/")<<space<<".dat\"";
	string com = ss.str();
	stringstream ss2 = bash::to_ss(com,to_string(th));
	while(1)
	{
		double ix, iy, iz;
		ss2>>ix;
		x.push_back(ix); if(ss2.eof()) break;
		ss2>>iy;
		y.push_back(iy); if(ss2.eof()) break;
		ss2>>iz;
		z.push_back(iz); if(ss2.eof()) break;
	}
}
*/

template<uint u> void acquisition::get(uint ir,uint ic,
														vector<double>& x,vector<double>& y,vector<double>& z)
{
	string space = u==Vec ? "vec" : (u==Mom ? "mom" : (u==Spin ? "spin" : ""));
	cout<<"Acquiring: r"<<ir<<"/c"<<ic<<"/"<<space<<".dat"<<endl;
	if(space == "") { cerr<<"template argument not valid!"<<endl; return; }
	stringstream ss;
	string tfn = "acq_" + space + "_" + to_string(th) + ".dat";
 	ss<<dlogin<<":~/"<<source<<"/r"<<ir<<"/c"<<ic<<(f_zero ? "_0/" : "/")<<space<<".dat " + tfn + " && echo \"0\"";
	string com = ss.str();
	//cout<<"Running command: "<<com<<endl;
	if(bash::to_int(com,to_string(th)))
	{
		cerr<<"Error running command: "<<com<<endl;
		return;
	}
	fstream f(tfn,fstream::in);
	if(!f.is_open())
	{
		cerr<<"Error opening file: "<<tfn<<endl;
		return;
	}
	while(1)
	{
		double ix, iy, iz;
		f>>ix;
		x.push_back(ix); if(f.eof()) break;
		f>>iy;
		y.push_back(iy); if(f.eof()) break;
		f>>iz;
		z.push_back(iz); if(f.eof()) break;
	}
	f.close();
	bash::to_int("rm " + tfn,to_string(th));
}

void acquisition::n(vector<double>& n) { for(int i=0; i<NL; i++) n.push_back(i+1); }

void acquisition::vec(uint ir, uint ic,
											vector<double>& x,vector<double>& y,vector<double>& z)
{ get<Vec>(ir,ic,x,y,z); }

void acquisition::mom(uint ir, uint ic,
											vector<double>& x,vector<double>& y,vector<double>& z)
{ get<Mom>(ir,ic,x,y,z); }

void acquisition::spin(uint ir, uint ic,
											vector<double>& x,vector<double>& y,vector<double>& z)
{ get<Spin>(ir,ic,x,y,z); }

void acquisition::xi(vector<double>& tf,vector<double>& td,vector<double>& ts,
										 vector<double>& cx,vector<double>& cy)
{
	stringstream ss;
	ss<<login<<" \"cat "<<source<<"/xi\"";
	string com = ss.str();
	stringstream ss2 = bash::to_ss(com,to_string(th));
	while(1)
	{
		double itf, itd, its, icx, icy;
		ss2>>itf;
		tf.push_back(itf); if(ss2.eof()) break;
		ss2>>itd;
		td.push_back(itd); if(ss2.eof()) break;
		ss2>>its;
		ts.push_back(its); if(ss2.eof()) break;
		ss2>>icx;
		cx.push_back(icx); if(ss2.eof()) break;
		ss2>>icy;
		cy.push_back(icy); if(ss2.eof()) break;
	}

}

array<double,8> acquisition::info(uint ir,uint ic)
{
	stringstream ss;
	ss<<login<<" \"cat "<<source<<"/r"<<ir<<"/c"<<ic<<"/info.dat\"";
	string com = ss.str();
	stringstream ss2 = bash::to_ss(com,to_string(th));
	array<double,8> oa;
	for(int i=0; i<8; i++) ss2>>com>>oa[i];
	return oa;
}
array<double,8> acquisition::info()
{
	stringstream ss;
	ss<<login<<" \"cat \\$(find "<<source<<" -type f -name \"info.dat\" | sed -n \"1p\")\"";
	string com = ss.str();
	stringstream ss2 = bash::to_ss(com,to_string(th));
	array<double,8> oa;
	for(int i=0; i<8; i++) ss2>>com>>oa[i];
	return oa;
}