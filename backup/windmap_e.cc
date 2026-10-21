#include<analyser.hh>

//void lv2::windmap(string path,bool out)
tuple<matrix<3,1>,matrix<3,3>> lv2::windmap_e(string wpath,bool out)
{
	ulong NS = 125; //343; //497; //
	//auto Dc = proton::line(wpath + "/rc", NS);
	//auto pt = lv1::pt(wpath + "/rc",Dc,0);

	Data<double> Dt(wpath + "/xi");
	valarray<double> Cx = Dt.col(4), Cy = Dt.col(5);
	matrix<3,1,valarray<double>> tx = {Dt.col(1),Dt.col(2),Dt.col(3)};
	valarray<double> eta1 = xi0[2] + Xi[2]*tx;
	fstream f;

	f.open((wpath + "/eps").data(),fstream::in);
	if(f.is_open()) { f.close(); }
	else
	{
		f.open((wpath + "/eps").data(),fstream::out | fstream::trunc);
		for(int i=0; i<NS; i++)
		{
			ulong N1 = 9; string path2 = wpath + "/r" + to_string(i+1);
			auto D = proton::line(path2,N1);
			matrix<1,3,valarray<double>> cr = {D.col("#xi_{x}"),D.col("#xi_{y}"),pow(D.col("#xi_{y}"),0)};
			cr[2] *= eta1[i];
			valarray<double> epx = D.col("#epsilon_{x}");
			valarray<double> epy = D.col("#epsilon_{y}");
			valarray<double> x0(N1), y0(N1);
			for(int j=0; j<N1; j++)
			{
				x0[j] = 0.01*((j/3) + 1);
				y0[j] = 0.01*((j%3) + 1);
				f<<x0[j]<<"\t"<<y0[j]<<"\t"
					<<cr[0][j]<<"\t"<<cr[1][j]<<"\t"<<cr[2][j]<<"\t"
					<<epx[j]<<"\t"<<epy[j]<<endl;//endl;//
			}

		}
		f.close();
		f.open((wpath + "/eps_info").data(),fstream::out | fstream::trunc);
		f<<"x0\ny0\n#xi_{x}\n#xi_{y}\n#eta_{1}\n"
			<<"#epsilon_{x}\n#epsilon_{y}"<<endl;//"<<endl;//
		f.close();
	}
	Data<double> Dw(wpath + "/eps");
	matrix<1,3,valarray<double>> cr = {Dw.col("#xi_{x}"),Dw.col("#xi_{y}"),Dw.col("#eta_{1}")};
	valarray<double> epx = Dw.col("#epsilon_{x}");
	valarray<double> epy = Dw.col("#epsilon_{y}");
	valarray<double> x0 = Dw.col("x0");
	valarray<double> y0 = Dw.col("y0");
	if(out)
	{
		new TCanvas("clatplot","clatplot",1900,1000);
		plot<double> lat_plot(cr[1][x0==0.01],epx[x0==0.01],"Sample Lattice",
		//plot<double> lat_plot(cr[2][x0==0.01 && y0==0.01],epx[x0==0.01 && y0==0.01],"Sample Lattice",
								"#eta_{1}","#epsilon_{x} (m rad)");
		plot2<double> lat_plot2(cr[0][x0==0.01 && y0==0.01],cr[1][x0==0.01 && y0==0.01],epx[x0==0.01 && y0==0.01],"Sample Lattice",
								"#xi_{x}","#xi_{y}","#epsilon_{x} (m rad)");
		//lat_plot.axis_equal();
		lat_plot2->Draw("AP0");
	}
	return make_tuple(xi0,Xi);
}

tuple<matrix<3,1>,matrix<3,3>> lv2::windmap_e(bool out) { return lv2::windmap_e(path,out); }



