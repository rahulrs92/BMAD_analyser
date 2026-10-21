#include<analyser.hh>
#include<TH3.h>
#include<TRandom3.h>

//void lv2::windmap(string path,bool out)
matrix<3,1,valarray<double>> lv2::sampler(bool out)
{
	/*
	Xi = { 50.679239553,  4.05289115646, 24.106319242,
				-38.7701652089,-40.7524635569,-17.4288629738,
				-18.8894366569,-2.38869014091,-9.49567097425};

	xi0 = {-0.0945,-0.0036,0.254957};
	Ci = Ci + C0 - xi0.tp();
	*/

	int sn = 10000;
	double rad = 0.01;//, dU = 1;
	matrix<3,1,valarray<double>> sample = valarray<double>(0.0,sn);

	valarray<double> evX;
	matrix<3,3> EVX;
	Xi.Eigens(evX,EVX);
	cout<<"Xi = \n"<<Xi<<endl;
	cout<<"EVX = "<<EVX<<endl;
	cout<<"evX = "<<evX<<endl;

	/*
	matrix<3,3> U, U1;//,Ui2;
	for(int l=0; l<3; l++)
	{
		U[l][l] = evX[l];
		dU *= evX[l];
	}
	cout<<dU<<endl;
	dU = dU/abs(dU) * pow(abs(dU),1.0/3.0);
	cout<<dU<<endl;
	cout<<"U = "<<U/dU<<endl;
	//cout<<"U1 = "<<U1<<endl;

	int n = 10;
	matrix<3,1,valarray<double>> sample;
	matrix<3,1> tm;
	double r = rad / n;
	matrix<3,1,valarray<double>> tmv_s = {{1,1,1,1,-1,-1,-1,-1},
																				{1,1,-1,-1,1,1,-1,-1},
																				{1,-1,1,-1,1,-1,1,-1}}, tmv;

	for(double i=0; i<=rad; i+=r)// *pow(U1[0][0],-0.25))
		for(double j=0; j<=rad; j+=r)// *pow(U1[1][1],-0.25))
			for(double k=0; k<=rad; k+=r)// *pow(U1[2][2],-0.25))
			{
				tm = {i,j,k};
				bool ell = tm.tp() * tm < pow(rad,2);
				for(int l=0; l<3; l++) tmv[l] = tmv_s[l] * tm[l];
				tmv = (EVX^-1) * U * tmv;
				if(ell)
					for(int m=0; m<8; m++)
					{
						if((i==0 && tmv_s[0][m]<0) || (j==0 && tmv_s[1][m]<0) || (k==0 && tmv_s[2][m]<0)) continue;
						for(int l=0; l<3; l++) v_sam[l].push_back(tmv[l][m]);
					}
			}
	for(int l=0; l<3; l++) sample[l] = valarray<double>(v_sam[l].data(),v_sam[l].size());
	*/
	auto trn3 = new TRandom3();
	double spr;
	for(int i=0; i<sn; i++)
	{
		spr = trn3->Rndm();
		spr *= 4*TMath::Pi()*pow(rad,2); spr = 0.5 * pow(spr/TMath::Pi(),0.5);
		trn3->Sphere(sample[0][i],sample[1][i],sample[2][i],spr);
	}

	//sample = (Xi^-1) * EVX * (U/dU) * (EVX^-1) * sample;
	//sample = (Xi^-1) * (sample - xi0);
	for(int i=0; i<3; i++) sample[i][sn-1] = 0.0;

	matrix<3,1> opt = {-1.35123,-4.527,1.1364379};
	//opt = EVX * (U/dU) * (EVX^-1) * opt;
	auto opt_mag = pow(opt.tp() * opt,0.5);
	cout<<"distance of optimized point from center: "<<opt_mag<<endl;
	//sample = sample - (Xi^-1)*(-xi0);
	auto sam_xi = xi0 + Xi*sample;
	cout<<"size of sample: "<<sample[0].size()<<endl;
	cout<<endl;
	cout<<"sample x_max: "<<sample[0].max()<<endl;
	cout<<"sample y_max: "<<sample[1].max()<<endl;
	cout<<"sample z_max: "<<sample[2].max()<<endl;
	cout<<endl;
	cout<<"sample x_min: "<<sample[0].min()<<endl;
	cout<<"sample y_min: "<<sample[1].min()<<endl;
	cout<<"sample z_min: "<<sample[2].min()<<endl;
	cout<<endl;
	cout<<endl;
	cout<<"sample xi_x_max: "<<sam_xi[0].max()<<endl;
	cout<<"sample xi_y_max: "<<sam_xi[1].max()<<endl;
	cout<<"sample eta_1_max: "<<sam_xi[2].max()<<endl;
	cout<<endl;
	cout<<"sample xi_x_min: "<<sam_xi[0].min()<<endl;
	cout<<"sample xi_y_min: "<<sam_xi[1].min()<<endl;
	cout<<"sample eta_1_min: "<<sam_xi[2].min()<<endl;

	cout<<"sam_xi:\n"<<sam_xi<<endl;

	if(out)
	{
		fstream f;
		f.open((path + "/sample").data(), fstream::out | fstream::trunc);
		for(int i=0; i<sample[0].size(); i++)
			f<<sample[0][i]<<"\t"<<sample[1][i]<<"\t"<<sample[2][i]<<endl;
		f.close();
		f.open((path + "/sample_stats").data(), fstream::out | fstream::trunc);
		f<<"size of sample: "<<sample[0].size()<<endl;
		f<<endl;
		f<<"sample x_max: "<<sample[0].max()<<endl;
		f<<"sample y_max: "<<sample[1].max()<<endl;
		f<<"sample z_max: "<<sample[2].max()<<endl;
		f<<endl;
		f<<"sample x_min: "<<sample[0].min()<<endl;
		f<<"sample y_min: "<<sample[1].min()<<endl;
		f<<"sample z_min: "<<sample[2].min()<<endl;
		f<<endl;
		f<<endl;
		f<<"sample xi_x_max: "<<sam_xi[0].max()<<endl;
		f<<"sample xi_y_max: "<<sam_xi[1].max()<<endl;
		f<<"sample eta_1_max: "<<sam_xi[2].max()<<endl;
		f<<endl;
		f<<"sample xi_x_min: "<<sam_xi[0].min()<<endl;
		f<<"sample xi_y_min: "<<sam_xi[1].min()<<endl;
		f<<"sample eta_1_min: "<<sam_xi[2].min()<<endl;
		f.close();

		new TCanvas("clatplot","clatplot",1900,1000);
		plot2<double> lat_plot(sample[0],sample[1],sample[2],"Sample Lattice",
										"t_{f} (Tesla/m^{2})","t_{d} (Tesla/m^{2})","t_{s} (Tesla/m^{2})");
		//lat_plot.axis_equal();
		lat_plot->Draw("AP0");

		new TCanvas("clatplot2","clatplot2",1900,1000);
		plot2<double> lat_plot2(sam_xi[0],sam_xi[1],sam_xi[2],"Sample Lattice",
										"#xi_{x}","#xi_{y}","#eta_{1}");
		//lat_plot.axis_equal();
		lat_plot2->Draw("surf1");
	}
	return sample;
}



