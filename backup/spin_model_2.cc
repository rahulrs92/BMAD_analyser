#include<lv2.hh>

tuple<double,matrix<3,1>,matrix<3,2>,matrix<2,1>> lv2::spin_model_2(Data<double>& D,string path,bool out)
{
	valarray<double> dt = D.col("#delta_{0}");
	valarray<double> nu = D.col("#nu_{s}");
	valarray<double> ps = D.col("T");

	matrix<1,3,valarray<double>> eps =
						{D.col("#epsilon_{x}"),D.col("#epsilon_{y}"),pow(dt,2)};
	matrix<3,1,valarray<double>> tx =
						{D.col("t_{f}"),D.col("t_{d}"),D.col("t_{s}")};
	matrix<2,1,valarray<double>> cr =
						{D.col("#xi_{x}"),D.col("#xi_{y}")};

	matrix<10,1,valarray<double>> regs_ps;
	regs_ps[0] = dt;
	for(int i=0; i<3; i++) regs_ps[i+1] = eps[i];
	for(int i=0; i<6; i++) regs_ps[i+4] = eps[i/2] * cr[i%2];
	auto arr_U = (sum(regs_ps * regs_ps.tp())^-1) * sum(regs_ps * (ps));
	eta0 = arr_U[0];
	//matrix<3,1> s;
	for(int i=0; i<3; i++) u[i] = arr_U[i+1];
	//matrix<3,2> T;
	for(int i=0; i<6; i++) U[i/2][i%2] = arr_U[i+4];

	cout<<"eta0 = "<<eta0<<endl;
	cout<<"u = "<<u<<endl;
	cout<<"U = "<<U<<endl;

	valarray<double> res_ps = ps - (eta0*dt + eps*u + eps*U*cr);
	double std_res_ps = pow(pow(res_ps,2).sum() / res_ps.size(),0.5);

	matrix<10,1,valarray<double>> regs_nu;
	regs_nu[0] = dt;
	for(int i=0; i<3; i++) regs_nu[i+1] = eps[i];
	for(int i=0; i<6; i++) regs_nu[i+4] = eps[i/2] * cr[i%2];
	auto arr_H = (sum(regs_nu * regs_nu.tp())^-1) * sum(regs_nu * nu);
	double sig0 = arr_H[0];
	//matrix<3,1> a;
	for(int i=0; i<3; i++) h[i] = arr_H[i+1];
	//matrix<3,2> V;
	for(int i=0; i<6; i++) H[i/2][i%2] = arr_H[i+4];

	cout<<"sig0 = "<<sig0<<endl;
	cout<<"h = "<<h<<endl;
	cout<<"H = "<<H<<endl;

	valarray<double> res_nu = nu - (sig0*dt + eps*h + eps*H*cr);
	double std_res_nu = pow(pow(res_nu,2).sum() / res_nu.size(),0.5);

	W = H - (sig0/eta0)*U;
	w = h - (sig0/eta0)*u;

	if(out)
	{
		TCanvas* cpts_ps = new TCanvas("cpts2","cpts2",2800,1000);
		plot<double> pts_ps(ps - eta0*dt,//[dt == dt.max()/2],//
											"Travel Time: Actual vs Linear Map","Data Point Index",
										"#frac{#DeltaT}{T} - #eta_{0}#delta");
		//plot<double> psm((arr_T2.tp() * regs_ps2) + eps[2]*cr[2]);//[dt == dt.max()/2]);//
		plot<double> psm(eps*u + eps*U*cr);//[dt == dt.max()/2]);//
		pts_ps->Draw("AP");
		pts_ps->SetMarkerStyle(3);
		pts_ps->SetLineColor(0);
		psm->Draw("P same");
		psm->SetMarkerColor(2);
		psm->SetMarkerStyle(2);
		psm->SetLineColor(0);
		auto leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts_ps.get(),"Data Point (Measured Values)");
		leg->AddEntry(psm.get(),"Phase Slip Model (Calculated Values)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cpts_ps->SetGrid();
		cpts_ps->SaveAs((path + "/error_plot_2.pdf(").data(),"pdf");

		TCanvas* cpts_nu = new TCanvas("cpts","cpts",2800,1000);

		plot<double> pts_nu(nu - sig0*dt,//[dt == dt.max()/2],//
										 "Spin tunes: Actual vs Linear Map","Data Point Index",
									 "#nu_{s} - #sigma_{0}#delta");
		plot<double> stm(eps*h + eps*H*cr);//arr_H.tp()*regs_nu);//(sig0*dt + (sig0*dt);//
		pts_nu->Draw("AP");
		pts_nu->SetMarkerStyle(3);
		pts_nu->SetLineColor(0);
		stm->Draw("P same");
		stm->SetMarkerColor(2);
		stm->SetMarkerStyle(2);
		stm->SetLineColor(0);
		leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts_nu.get(),"Data Point (Measured Values)");
		leg->AddEntry(stm.get(),"Spin Tune Model (Calculated Values)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cpts_nu->SetGrid();
		cpts_nu->SaveAs((path + "/error_plot_2.pdf)").data(),"pdf");

		TCanvas* ceh_ps = new TCanvas("ceh2","ceh2",1024,768);
		TH1F* hps = new TH1F("hps",
			"Distribution of residuals of travel-time fit;#frac{#DeltaT}{T}_{simulated} - #frac{#DeltaT}{T}_{model};Frequency",
			51,-5* std_res_ps,5* std_res_ps);
		for(int i=0; i<res_ps.size(); i++) hps->Fill(res_ps[i]);
		//hps->Fit("gaus");
		hps->Draw();
		ceh_ps->SaveAs((path + "/error_dist_2.pdf(").data(),"pdf");

		TCanvas* ceh_nu = new TCanvas("ceh","ceh",1024,768);
		TH1F* hnu = new TH1F("hnu",
			"Distribution of residuals of spin-tune fit;#nu_{simulated} - #nu_{model};Frequency",
			51,-5*std_res_nu,5*std_res_nu);
		for(int i=0; i<res_nu.size(); i++) hnu->Fill(res_nu[i]);
		//hnu->Fit("gaus");
		hnu->Draw();
		ceh_nu->SaveAs((path + "/error_dist_2.pdf)").data(),"pdf");
	}
	return make_tuple(sig0,u,U,matrix<2,1>({0,0}));
}




