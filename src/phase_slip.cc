#include<proton.hh>

array<double,3> proton_RF::phase_slip(bool Out)
{
	if(!f_PSF)
	{
		Tp_model();

		vector<double> vmdt, vps;
		for(int i=0; i<N-1; i++)
		{
			vmdt.push_back((pz[i]+pz[i+1])/2);
			vps.push_back(-(z[i+1]-z[i])/ L);
		}

		valarray<double> mdt = valarray<double>(vmdt.data(),N-1);
		valarray<double> ps = valarray<double>(vps.data(),N-1);
		valarray<double> n(0.0,N-1);
		for(int i=0; i<N-1; i++) n[i] = i;

		matrix<2,1,valarray<double>> vA = {mdt,pow(mdt,2)};
		auto A = sum(vA * vA.tp());
		auto B = (A^-1) * sum(vA * (ps - Tps));

		mg["eta"] = new TMultiGraph("eta","Phase Slip;turn no.;#DeltaT/T");

		plot<double> ps_data(n, ps,"#DeltaT/T per turn");
		ps_data->SetMarkerStyle(1);
		ps_data->SetLineColor(10);
		mg["eta"]->Add(ps_data.get(),"P");
		ss.clear();
		ss<<"#frac{#DeltaT}{T} = "<<(Tps<0 ? "- " : "")<<abs(Tps)
															<<(B[0]<0 ? " - " : " + ")<<abs(B[0])<<" #delta"
															<<(B[1]<0 ? " - " : " + ")<<abs(B[1])<<" #delta^{2}"<<endl;
		string ts; getline(ss,ts);
		plot<double> ps_fit(n,Tps + B[0]*mdt + B[1]*pow(mdt,2),ts);
		ps_fit->SetMarkerStyle(1);
		ps_fit->SetLineColor(2);
		mg["eta"]->Add(ps_fit.get(),"L");

		PSF = array<double,3>({Tps,B[0],B[1]});
		f_PSF = 1;
	}
	if(Out)
	{
		TCanvas* cps = new TCanvas("cmc","eta-fit");
		mg["eta"]->Draw("AP");
		cps->SetGrid();
		mg["eta"]->GetYaxis()->SetTitleOffset(0.9);
		auto leg = cps->BuildLegend(0.11,0.25,0.64,0.11);
		leg->SetFillColorAlpha(10,0.2);
		cps->SaveAs((path + "/eta-fit.pdf(").data(),"pdf");
		mg["eta"]->GetXaxis()->SetRange(N>1e3 ? N/2 - 5e2 : 0,N>1e3 ? N/2 + 5e2 : N);
		cps->SaveAs((path + "/eta-fit.pdf)").data(),"pdf");
		mg["eta"]->GetXaxis()->SetRange(0,N);
	}
	return PSF;
}

array<double,2> proton::phase_slip(bool Out)
{
	if(!f_PSF)
	{
		Tp_model();

		valarray<double> vn(n.data(),N), ps(z.data(),N);
		ps = -ps / L;

		double eta = (ps *vn).sum() / pow(vn,2).sum();
		PSF = array<double,2>({Tps,eta- Tps});

		mg["eta"] = new TMultiGraph("eta",("Cumulative Phase Slip at each turn" +
																			string(pz[0]==0 ? "@ #delta = 0;" : ";" ) +
																			"turn no.;#Sigma#DeltaT/T").data());
		plot<double> ps_data(ps,"Total Phase Slip");
		ps_data->SetMarkerStyle(1);
		ps_data->SetLineColor(10);
		mg["eta"]->Add(ps_data.get(),"P");
		plot<double> ps_t(Tps *vn,"Theoretical Transverse Contribution");
		ps_t->SetMarkerStyle(1);
		ps_t->SetMarkerColor(4);
		ps_t->SetLineColor(4);
		mg["eta"]->Add(ps_t.get(),"L");
		ss.str() = "";
		ss<<"#frac{#DeltaT}{T} = "<<(eta<0 ? "- " : "")<<abs(eta)<<" n"<<endl;
		string ts; getline(ss,ts);
		plot<double> ps_fit(eta*vn,ts);
		ps_fit->SetLineColor(2);
		ps_fit->SetMarkerStyle(1);
		mg["eta"]->Add(ps_fit.get(),"L");
		f_PSF = 1;
	}
	if(Out)
	{
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
