#include<analyser.hh>

void lv2::tune_vs_emittance(string epath,bool out)
{
	auto D = proton::line(epath,891);
	D.clean_nans();
	valarray<double> dt = D.col("#delta_{0}");
	valarray<double> nu = D.col("#nu_{s}");
	valarray<double> ps = D.col("T");
	valarray<double> epx0 = D.col("#epsilon_{x}");
	valarray<double> epy0 = D.col("#epsilon_{y}");

	for(int i=0; i<epx0.size(); i++) epx0[i] = epx0[11 * (i/11)];
	for(int i=0; i<epy0.size(); i++) epy0[i] = epy0[11 * (i/11)];

	matrix<1,3,valarray<double>> eps =
						{epx0,epy0,pow(dt,2)};
	matrix<3,1> a0;

	matrix<3,1,valarray<double>> regs_nu;
	//for(int i=0; i<2; i++) regs_nu[i] = pow(eps[i],0.5);
	regs_nu[0] = dt;
	for(int i=0; i<2; i++) regs_nu[i+1] = eps[i];
	auto arr_V = (sum(regs_nu * regs_nu.tp())^-1) * sum(regs_nu * (nu));
	sig0 = arr_V[0];
	//for(int i=0; i<3; i++) a0[i] = arr_V[i];
	for(int i=0; i<2; i++) a[i] = arr_V[i+1];
	cout<<"sig0 = "<<sig0<<endl;
	cout<<"a0 = "<<a0<<endl;
	cout<<"a = "<<a<<endl;

	fstream f;
	f.open((epath + "/spin_model_oo").data(), fstream::out | fstream::trunc);
	f<<sig0<<"\n"<<endl;
	f<<a<<endl;
	f.close();

	valarray<double> res_nu = nu - arr_V.tp()*regs_nu;
	double std_res_nu = pow(pow(res_nu,2).sum() / res_nu.size(),0.5);

	if(out)
	{
		TCanvas* cpts = new TCanvas("cpts","cpts",2800,1000);

		plot<double> pts(nu,//[dt == dt.max()/2],//
										 "Spin tunes: Actual vs Model at Optical Origin","Data Point Index","#nu_{s}");
		plot<double> stm(arr_V.tp()*regs_nu);//(sig0*dt + eps*a + eps*V*cr));//
		pts->Draw("AP");
		pts->SetMarkerStyle(3);
		pts->SetLineColor(0);
		stm->Draw("P same");
		stm->SetMarkerColor(2);
		stm->SetMarkerStyle(2);
		stm->SetLineColor(0);
		auto leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts.get(),"Data Point (Measured Values)");
		leg->AddEntry(stm.get(),"Spin Tune Model (Calculated Values)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cpts->SetGrid();
		cpts->SaveAs((epath + "/error_plot_2.pdf").data(),"pdf");

		TCanvas* ceh = new TCanvas("ceh","ceh",1024,768);

		TH1F* hnu = new TH1F("hnu",
			"Distribution of residuals of spin-tune fit;#nu_{simulated} - #nu_{model};Frequency",
			51,-5* std_res_nu,5* std_res_nu);
		for(int i=0; i<res_nu.size(); i++) hnu->Fill(res_nu[i]);
		//hnu->Fit("gaus");
		hnu->Draw();
		ceh->SaveAs((epath + "/error_dist_2.pdf").data(),"pdf");

	}

}
