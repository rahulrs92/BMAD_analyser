#include<proton.hh>

array<double,3> proton_RF::phase_slip(bool Out)
{
	if(!f_PSF)
	{
		Tp_model();

		//Tps = isnan(Tps) ? 0 : Tps;

		vector<double> vmdt, vps, r_mdt, r_ps;
		for(int i=0; i<N-1; i++)
		{
			vmdt.push_back((pz[i]+pz[i+1])/2);
			vps.push_back(-(z[i+1]-z[i]) / L[mode]);
		}
		/*
		map<double,array<double,2>> map_dp;
		for(int i=0; i<N-1; i++) map_dp[abs(atan2(x[i],px[i]))] = {vmdt[i], vps[i]};
		
		//for(auto& [x_val, dp] : map_dp) cout<<x_val<<"\t"<<dp[0]<<"\t"<<dp[1]<<endl;
		
		for(auto& [x_val, dp] : map_dp)
		{
			r_mdt.push_back(dp[0]);
			r_ps.push_back(dp[1]);
			if(r_mdt.size() > N/10 || (N < 1000 && r_mdt.size() > 10)) break;
		}
		int N1 = r_mdt.size();
		*/

		valarray<double> mdt(vmdt.data(),N-1);//, rmdt(r_mdt.data(),N1), rps(r_ps.data(),N1);
		valarray<double> ps = valarray<double>(vps.data(),N-1);
		valarray<double> n(0.0,N-1);
		for(int i=0; i<N-1; i++) n[i] = i;

		matrix<4,1,valarray<double>> vA = {pow(mdt,0),mdt,pow(mdt,2),pow(mdt,3)};
		auto A = sum(vA * vA.tp())^-1;
		auto B = A * sum(vA * valarray<double>(ps - Tps));

		// if(Out) cout<<"Theoretical Transverse Contribution: "<<Tps<<endl;
		// if(Out) cout<<"Travel-Time Intercept: "<<B[0]<<endl;

		valarray<double> res_ps = ps - (Tps + B[0] + B[1]*mdt + B[2]*pow(mdt,2) + B[3]*pow(mdt,3));
		auto E = A * sum((vA * res_ps * res_ps * vA.tp())) * A;
		// if(Out) cout<<"Coefficient variances of Phase Slip Fit: "<<endl;
		// if(Out) cout<<E<<endl;

		mg["eta"] = new TMultiGraph("eta",
			("Change in Travel Time at each turn @ " + mode + " (" + ID + ");Turn no.;#Sigma#Deltat/t").data());

		plot<double> ps_data(vps,"#DeltaT/T at given turn");
		ps_data->SetMarkerStyle(1);
		ps_data->SetLineColor(10);
		mg["eta"]->Add(ps_data.get(),"P");

		ss.clear();
		ss<<"#frac{#DeltaT}{T} = "<<(B[0]<0 ? "- " : "")<<abs(B[0])
		<<(B[1]<0 ? " - " : " + ")<<abs(B[1])<<" #delta"
		<<(B[2]<0 ? " - " : " + ")<<abs(B[2])<<" #delta^{2}"
		<<(B[3]<0 ? " - " : " + ")<<abs(B[3])<<" #delta^{3}"<<endl;
		string ts; getline(ss,ts);

		plot<double> ps_fit(n,Tps + B[0] + B[1]*mdt + B[2]*pow(mdt,2) + B[3]*pow(mdt,3),ts);
		ps_fit->SetMarkerStyle(1);
		ps_fit->SetLineColor(2);
		mg["eta"]->Add(ps_fit.get(),"L");
		
		mg["eta-fit"] = new TMultiGraph("eta-fit",
			("Change in Travel Time vs. #delta @ " + mode + " (" + ID + ");#delta;#Sigma#Deltat/t").data());

		plot<double> ef_data(mdt, ps,"#DeltaT/T at given turn");
		ef_data->SetMarkerStyle(1);
		ef_data->SetLineColor(10);
		mg["eta-fit"]->Add(ef_data.get(),"P");
		/*
		plot<double> rps_data(r_mdt, r_ps,"#DeltaT/T at given turn");
		rps_data->SetMarkerStyle(4);
		rps_data->SetLineColor(10);
		mg["eta-fit"]->Add(rps_data.get(),"P");
		*/
		plot<double> ef_fit(mdt,Tps + B[0] + B[1]*mdt + B[2]*pow(mdt,2) + B[3]*pow(mdt,3),ts);
		ef_fit->SetMarkerStyle(1);
		ef_fit->SetLineColor(2);
		mg["eta-fit"]->Add(ef_fit.get(),"L");

		PSF = array<double,3>({Tps,B[1],B[2]});
		e_PSF = array<double,3>({0,pow(E[1][1],0.5),pow(E[2][2],0.5)});

		if(Out) cout<<"Phase Slip Factors: "<<PSF[0]<<", "<<PSF[1]<<", "<<PSF[2]<<endl;
		if(Out) cout<<"Phase Slip Factor Errors: "<<e_PSF[0]<<", "<<e_PSF[1]<<", "<<e_PSF[2]<<endl;
		f_PSF = 1;
	}
	if(Out)
	{
		TCanvas* cps = new TCanvas("cmc","eta");
		mg["eta"]->Draw("AP");
		cps->SetGrid();
		mg["eta"]->GetYaxis()->SetTitleOffset(0.9);
		auto leg = cps->BuildLegend(0.11,0.25,0.64,0.11);
		leg->SetFillColorAlpha(10,0.2);
		cps->SaveAs((path + "/eta.pdf(").data(),"pdf");
		mg["eta"]->GetXaxis()->SetRange(N>1e3 ? N/2 - 5e2 : 0,N>1e3 ? N/2 + 5e2 : N);
		cps->SaveAs((path + "/eta.pdf)").data(),"pdf");
		mg["eta"]->GetXaxis()->SetRange(0,N);

		TCanvas* cef = new TCanvas("cef","eta-fit");
		mg["eta-fit"]->Draw("AP");
		cef->SetGrid();
		mg["eta-fit"]->GetYaxis()->SetTitleOffset(0.9);
		auto leg2 = cef->BuildLegend(0.11,0.25,0.64,0.11);
		leg2->SetFillColorAlpha(10,0.2);
		cef->SaveAs((path + "/eta-fit.pdf").data(),"pdf");
	}
	return PSF;
}

array<double,2> proton::phase_slip(bool Out)
{
	if(!f_PSF)
	{
		Tp_model();

		/*
		valarray<double> vn(n.data(),N), ps(z.data(),N);
		ps = -ps / L[mode];
		
		double eta = (ps *vn).sum() / pow(vn,2).sum();
		PSF = array<double,2>({Tps,eta});
		*/
		double dL = 4.10810986E-02; //manual path length correction for COSY misaligned lattice

		valarray<double> vn(n.data(),N), ps(z.data(),N);
		ps = -(ps) / (L[mode]);

		matrix<2,1,valarray<double>> vA = {pow(vn,0),pow(vn,1)};
		auto eta = (sum(vA * vA.tp())^-1) * sum(vA * ps);
		PSF = array<double,2>({Tps,eta[1] - dL/L[mode]});

		if(!f_tree || Out)
		{
			mg["eta"] = new TMultiGraph("eta",("Cumulative Phase Slip at each turn" +
																				string(pz[0]==0 ? "@ #delta = 0;" : ";" ) +
																				"turn no.;#Sigma#Deltat/t").data());
			plot<double> ps_data(ps,"Total Phase Slip");
			ps_data->SetMarkerStyle(1);
			ps_data->SetLineColor(10);
			mg["eta"]->Add(ps_data.get(),"P");
			plot<double> ps_t(Tps * vn,"Theoretical Transverse Contribution");
			ps_t->SetMarkerStyle(1);
			ps_t->SetMarkerColor(4);
			ps_t->SetLineColor(4);
			mg["eta"]->Add(ps_t.get(),"L");
			ss.str() = "";
			ss<<"#frac{#DeltaT}{T} = "<<(eta[0]<0 ? "- " : "")<<abs(eta[0])<<(eta[1]<0 ? " - " : " + ")<<abs(eta[1])<<" n"<<endl;
			string ts; getline(ss,ts);
			plot<double> ps_fit(eta.tp()*vA,ts);
			ps_fit->SetLineColor(2);
			ps_fit->SetMarkerStyle(1);
			mg["eta"]->Add(ps_fit.get(),"L");
		}
		f_PSF = 1;
	}
	if(Out)
	{
		cout<<"eta = "<<PSF[1]<<endl;
		TCanvas* cps = new TCanvas("cps","phase_slip");
		mg["eta"]->Draw("AP");
		mg["eta"]->GetYaxis()->SetTitleOffset(0.9);
		auto leg = cps->BuildLegend(0.36,0.25,0.89,0.11);
		leg->SetFillColorAlpha(10,0.2);
		cps->SetGrid();
		cps->SaveAs((path + "/eta-fit.pdf").data(),"pdf");
	}
	return PSF;
}
