#include<lv1.hh>
#include<mutex>

mutex mtx;

array<double,7> lv1::point(acquisition& ia,uint ir,uint ic,double icx,double icy,bool out)
{
	vector<vector<double>> C(10,vector<double>());
	ia.n(C[0]);
	ia.vec(ir,ic,C[1],C[2],C[3]);
	ia.mom(ir,ic,C[4],C[5],C[6]);
	ia.spin(ir,ic,C[7],C[8],C[9]);
	auto info = ia.info(ir,ic);
	proton P(info[1], info[2], info[3], info[4], C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8], C[9]);
	auto CSx = P.Courant_Snyder(0,out);
	auto CSy = P.Courant_Snyder(1,out);
	auto eta = P.phase_slip(out);
	double nus = P.spin_tune(out);
	//x_px = &P.mg["CSx"];
	//y_py = &P.mg["CSy"];

	double epx = CSx[0], epy = CSy[0];
	double epx0 = CSx[6], epy0 = CSy[6];
	double dt0 = C[6][0]; // j*0.0002;//
	double ps = eta[1];//+eta[0];

	return {epx,epy,epx0,epy0,dt0,ps,nus};
}
/*
array<double,11> lv1::point_RF(acquisition& ia,uint ir,uint ic,bool out)
{
	vector<vector<double>> C(10,vector<double>());
	ia.n(C[0]);
	ia.vec(ir,ic,C[1],C[2],C[3]);
	ia.mom(ir,ic,C[4],C[5],C[6]);
	ia.spin(ir,ic,C[7],C[8],C[9]);
	auto inf = ia.info(ir,ic);
	proton_RF P(inf[1], inf[2], inf[3], inf[4], C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8], C[9]);
	//cout<<"size of C: "<<C.size()<<" : "<<C[1].size()<<endl;
	auto CSx = P.Courant_Snyder(0,out);
	auto CSy = P.Courant_Snyder(1,out);
	auto syn = P.longitudinal(out);
	//cout<<"Thread = "<<ia.th<<endl;
	auto eta = P.phase_slip(out);
	
	double kt1 = inf[5], kt2 = inf[6], kt3 = inf[7];
	double Cx = inf[3], Cy = inf[4];
	double epx = CSx[0], epy = CSy[0];
	double dta = syn[3], om = syn[4];
	double ps = eta[1];
	double nus = P.spin_tune(out);
	
	//cout<<"Spin tune: "<<nus<<endl;
	
	return {kt1,kt2,kt3,Cx,Cy,epx,epy,dta,om,ps,nus};
}
*/

array<double,11> lv1::point_RF(acquisition& ia,uint ir,uint ic,bool out)
{
	stringstream ss;
 	ss<<ia.dlogin<<":~/"<<acquisition::source<<"/r"<<ir<<"/c"<<ic<<" "<<ia.destination<<"/c"<<ia.th<<" || echo \"1\" && echo \"0\"";
	string com = ss.str();
	cout<<"Currently downloading: "<<acquisition::source<<"/r"<<ir<<"/c"<<ic<<endl;

	while(1)
	{
		//mtx.lock();
		int r = bash::to_int(com,to_string(ia.th));
		//cout<<"r = "<<r<<endl;
		if(r == 0)
		{
			//cout<<"if condition entered!!!"<<endl;
			fstream f(ia.destination + "/c" + to_string(ia.th) + "/vec.dat",fstream::in);
			if(f.is_open())
			{
				//cout<<"File opened: "<<ia.destination + "/c" + to_string(ia.th) + "/vec.dat"<<endl;
				f.close();
				break;
			}
			else
			{
				cerr<<"Error opening file: "<<ia.destination + "/c" + to_string(ia.th) + "/vec.dat"<<endl;
				cout<<"\nRetrying Download..."<<endl;
				continue;
			}
		}
		else
		{
			cerr<<"Error running command: "<<com<<endl;
			continue;
		}
		//mtx.unlock();
		
	}
	proton_RF P(ia.destination + "/c" + to_string(ia.th));
	auto CSx = P.Courant_Snyder(0,out);
	auto CSy = P.Courant_Snyder(1,out);
	auto syn = P.longitudinal(out);
	//cout<<"Thread = "<<ia.th<<endl;
	auto eta = P.phase_slip(out);
	
	double kt1 = P.Kt1(), kt2 = P.Kt2(), kt3 = P.Kt3();
	double Cx = P.Cx(), Cy = P.Cy();
	double epx = CSx[0], epy = CSy[0];
	double dta = syn[3], om = syn[4];
	double ps = eta[1];
	double nus = P.spin_tune(out);

	bash::to_int("rm -r " + ia.destination + "/c" + to_string(ia.th),to_string(ia.th));
	
	//cout<<"Spin tune: "<<nus<<endl;
	
	return {kt1,kt2,kt3,Cx,Cy,epx,epy,dta,om,ps,nus};
}

array<double,7> lv1::point(string path,uint ir,uint ic,bool out)
{
	proton P(path + "/r" + to_string(ir) + "/c" + to_string(ic));
	auto CSx = P.Courant_Snyder(0,out);
	auto CSy = P.Courant_Snyder(1,out);
	auto eta = P.phase_slip(out);
	double nus = P.spin_tune(out);

	double epx = CSx[0], epy = CSy[0];
	double epx0 = CSx[6], epy0 = CSy[6];
	double dt0 = P.pz[0]; // j*0.0002;//
	double ps = eta[1];//+eta[0];

	return {epx,epy,epx0,epy0,dt0,ps,nus};
}

array<double,11> lv1::point_RF(string path,uint ir,uint ic,bool out)
{
	proton_RF P(path + "/r" + to_string(ir) + "/c" + to_string(ic));
	auto CSx = P.Courant_Snyder(0,out);
	auto CSy = P.Courant_Snyder(1,out);
	auto syn = P.longitudinal(out);
	auto eta = P.phase_slip(out);
	double nus = P.spin_tune(out);

	double kt1 = P.Kt1(), kt2 = P.Kt2(), kt3 = P.Kt3();
	double Cx = P.Cx(), Cy = P.Cy();
	double epx = CSx[0], epy = CSy[0];
	double dta = syn[3], om = syn[4];
	double ps = eta[1];
	
	return {kt1,kt2,kt3,Cx,Cy,epx,epy,dta,om,ps,nus};
}

array<double,5> lv1::xi(string path,uint ir,uint ic,bool out)
{
	proton P(path + "/r" + to_string(ir) + "/c" + to_string(ic));
	//cout<<ir<<", "<<ic<<endl;
	return P.xi();
}
