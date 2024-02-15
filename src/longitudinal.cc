#include<proton.hh>


array<double,5> proton_RF::longitudinal(bool Out)
{
	if(!f_syn)
	{
		/*
		vector<double> v_rz, v_rpz, v_rn;
		valarray<double> vx(x.data(),N), vpx(px.data(),N);
		valarray<double> vmag = pow(pow(vx,2) + pow(vpx,2),0.5);
		double min_vmag = vmag.max();
		for(int i=1; i<N; i++)
		{
			if(abs((px[i-1] - vpx.min())/(px[i] - vpx.min())) > 100)
			{
				v_rz.push_back(z[i]);
				v_rpz.push_back(pz[i]);
				v_rn.push_back(n[i]);
				//cout<<"yay!!!"<<endl;
			}
		}
		cout<<"Size of v_rz: "<<v_rz.size()<<endl;
		valarray<double> vz(v_rz.data(),v_rz.size()), vpz(v_rpz.data(),v_rpz.size()), vn(v_rn.data(),v_rn.size());
		*/
		valarray<double> vz(z.data(),N), vpz(pz.data(),N), vn(n.data(),N);

		matrix<4,1,valarray<Double_t>> vA({pow(vz,2),pow(vpz,2),vz,vpz});
		auto A = sum(vA*vA.tp());
		auto B = (A^(-1))*sum(vA);

		double eta0 = phase_slip(Out)[1];
		Double_t k = -L *eta0;
		Double_t om = k * sqrt(B[0]/B[1]);
		Double_t zm = - B[2]/(2*B[0]);
		Double_t dm = - B[3]/(2*B[1]);
		Double_t da = sqrt(1/B[1] + pow(B[2],2)/(4*B[0]*B[1]) + pow(B[3],2)/(4*pow(B[1],2)));
		Double_t za = da*k/om;

		syn = {zm,za,dm,da,om};

		valarray<Double_t> elx = zm + za*sin(om*vn);
		valarray<Double_t> ely = dm + da*cos(om*vn);

		mg["lon"] = new TMultiGraph("lon","Longitudinal Phase Space;z;#delta");
		plot<double> lon(z,pz);
		lon->SetMarkerStyle(0);
		lon->SetLineColor(10);
		mg["lon"]->Add(lon.get(),"P");
		/*
		plot<double> lon_r(vz,vpz);
		lon_r->SetMarkerStyle(2);
		lon_r->SetLineColor(10);
		mg["lon"]->Add(lon_r.get(),"P");
		*/
		ss.clear();
		ss<<"#frac{#(){z - "<<zm<<"}^{2}}{"<<k*da/om<<"^{2}} + "
			<<"#frac{#(){#delta - "<<dm<<"}^{2}}{"<<da<<"^{2}} = 1";
		string ts; getline(ss,ts);
		plot<double> lon_fit(elx,ely,ts);
		lon_fit->SetMarkerStyle(0);
		lon_fit->SetMarkerColor(2);
		lon_fit->SetLineColor(2);
		mg["lon"]->Add(lon_fit.get(),"L");
	}
	if(Out)
	{
		cout<<"Synchrotron Frequency = "<<syn[4]<<" radians per turn at eta0 = "<<PSF[1]<<endl;
		cout<<"Synchrotron period = "<<2*pi/syn[4]<<" turns"<<endl;
		cout<<"Mean shift in z = "<<syn[0]<<" metres"<<endl;
		cout<<"Mean shift in momentum offset = "<<syn[2]<<endl;
		cout<<"Amplitude of z oscillations = "<<syn[1]<<" metres"<<endl;
		cout<<"Amplitude of momentum offset oscillations = "<<syn[3]<<endl;

		auto clon = new TCanvas("clon","Phase Space Ellipse");
		mg["lon"]->Draw("AP");
		mg["lon"]->GetYaxis()->SetTitleOffset(0.9);
		auto leg = clon->BuildLegend(0.11,0.25,0.64,0.11);//(0.7,0.9,0.9,0.8);
		clon->SetGrid();
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		clon->SaveAs((path + "/longitudinal.pdf").data(),"pdf");
	}
	return syn;
}

