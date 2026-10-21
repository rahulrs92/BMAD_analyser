#include<lv1.hh>

void lv1::SCT_merge(acquisition& aq,string path)
{
	int nr = aq.Nr, nc = aq.Nc, nl = aq.NL;
	cout<<"nr: "<<nr<<endl;
	cout<<"nc: "<<nc<<endl;
	cout<<"nl: "<<nl<<endl;
	for(int i=1; i<=nr; i++)
	{
		bash::to_string("mkdir " + path + "/c" + to_string(i));
		for(int q=0; q<3; q++)
		{
			vector<double> x(nl,0.0), y(nl,0.0), z(nl,0.0);
			for(int j=1; j<=nc; j++)
			{
				vector<double> jx, jy, jz;
				if(q==0) aq.vec(i,j,jx,jy,jz);
				if(q==1) aq.mom(i,j,jx,jy,jz);
				if(q==2) aq.spin(i,j,jx,jy,jz);
				for(int k=0; k<nl; k++)
				{
					x[k] += jx[k]/nc;
					y[k] += jy[k]/nc;
					z[k] += jz[k]/nc;
				}
			}
			stringstream ss;
			ss<<path<<"/c"<<i;
			if(q==0) ss<<"/vec.dat";
			if(q==1) ss<<"/mom.dat";
			if(q==2) ss<<"/spin.dat";
			cout<<"Writing file: "<<ss.str()<<endl;
			fstream f(ss.str().data(),fstream::out | fstream::trunc);
			for(int j=0; j<nl; j++) f<<x[j]<<"\t"<<y[j]<<"\t"<<z[j]<<endl;
			f.close();
		}
		auto info = aq.info(i,1);
		fstream f(path + "/c" + to_string(i) + "/info.dat",fstream::out | fstream::trunc);
		f<<"L:\t"<<info[0]<<"\n";
		f<<"Qx:\t"<<info[1]<<"\n";
		f<<"Qy:\t"<<info[2]<<"\n";
		f<<"Cx:\t"<<info[3]<<"\n";
		f<<"Cy:\t"<<info[4]<<"\n";
		f<<"tF:\t"<<info[5]<<"\n";
		f<<"tD:\t"<<info[6]<<"\n";
		f<<"tS:\t"<<info[7]<<"\n";
		f.close();
		cout<<"Merged data for r"<<i<<endl;
	}
}
