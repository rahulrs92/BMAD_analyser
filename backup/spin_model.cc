#include<analyser.hh>

tuple<matrix<3,1>,matrix<3,3>,matrix<3,1>> lv2::spin_model(Data<double>& D,Data<double>& D0,string path,bool out)
{
	//lv1::grid(ia,path,0);
	//D.clean_nans();
	//D0.clean_nans();
	valarray<double> dt = D.col("#delta_{0}");
	valarray<double> nu = D.col("#nu_{s}");
	valarray<double> ps = D.col("T");

	matrix<1,3,valarray<double>> eps =
						{D.col("#epsilon_{x}"),D.col("#epsilon_{y}"),pow(dt,2)};
	matrix<3,1,valarray<double>> tx =
						{D.col("t_{f}"),D.col("t_{d}"),D.col("t_{s}")};
	matrix<3,1,valarray<double>> cr =
						{D.col("#xi_{x}"),D.col("#xi_{y}"),pow(dt,2)};

	/*
	matrix<4,1,valarray<double>> regs_wm = {valarray<double>(1.0,tx[0].size()),tx[0],tx[1],tx[2]};
	auto X1 = ((sum(regs_wm * regs_wm.tp())^-1) * sum(regs_wm * cr.sub<2,1>().tp())).tp();
	xi0.sub<2,1>() = X1.sub<2,1>();
	Xi.sub<2,3>() = X1.sub<2,3>(0,1);

	//calculate the null space of the 2x3 transformation matrix:
	Xi[2][2] = 1;
	Xi.sub<1,2>(2,0) = ((X1.sub<2,2>(0,1)^-1) * (-X1).sub<2,1>(0,3)).tp();
	cr[2] = Xi[2] * tx;

	matrix<3,3> TX;
	auto xt = Xi.sub<2,3>() * tx;
	matrix<13,1,valarray<double>> regs_ps;
	regs_ps[0] = dt;
	for(int i=0; i<3; i++) { regs_ps[i+1] = eps[i]; }
	for(int i=0; i<9; i++) { regs_ps[i+4] = eps[i/3] * tx[i%3]; }
	auto arr_T = (sum(regs_ps * regs_ps.tp())^-1) * sum(regs_ps * ps);
	eta0 = arr_T[0];
	xi0[2] = arr_T[3];
	for(int i=0; i<9; i++) TX[i/3][i%3] = arr_T[i+4];
	Xi[2] = TX[2];

	cout<<"First localization...\nxi0 = "<<xi0<<endl;
	cout<<"Xi = "<<Xi<<endl;
	cout<<"T = "<<TX * (Xi^-1)<<endl;
	cout<<"eta0 = "<<eta0<<endl;
	*/


	/*
	T.sub<2,3>() = TX.sub<2,3>() * (Xi^-1);
	cout<<"Confirmation of xi0..."<<endl;
	cout<<"Tpj * xi0 = "<<T.sub<2,3>() * xi0<<endl;
	cout<<"new Tpj * xi0 = "<<xi0_1<<endl;
	TX[2][2] = 1;



	T = TX * (Xi^-1);
	cout<<"new xi0 = "<<(T^-1) * b<<endl;
	Xi[2] *= T[2][2];
	T(2) /= T[2][2];
	cout<<"T = "<<T<<endl;

	cout<<"current xi0 = "<<xi0<<endl;
	xi0[2] = (T^-1)[2] * b;
	cout<<"new xi0 = "<<(T^-1) * b<<endl;

	Xi[2] = TX[2];
	T.sub<2,2>() = TX.sub<2,2>() * (Xi.sub<2,2>()^-1);
	cout<<"Confirmation of Tpq...\nCp3 = "<<TX.sub<2,1>(0,2)<<endl;
	cout<<"Tpq*Xq3 = "<<T.sub<2,2>()*Xi.sub<2,1>(0,2)<<endl;

	cout<<"current xi0 = "<<xi0<<endl;
	//xi0.sub<2,1>() = (T.sub<2,2>()^-1) * b.sub<2,1>();
	xi0[2] = b[2];
	cout<<"new xi0 = "<<xi0<<endl;
	T[2][2] = 1;
	cout<<"eta0 = "<<eta0<<endl;
	*/

	//cr[2] = valarray<double>(v_cy.data(),v_cy.size());
	cr[2] = xi0[2] + Xi[2]*tx;

	/*
	matrix<5,1,valarray<double>> regs_ps2;
	//regs_nu[0] = nu*0 + 1;
	regs_ps2[0] = dt;
	for(int i=0; i<4; i++) { regs_ps2[i+1] = eps[i/2] * cr[i%2]; }
	auto arr_T2 = (sum(regs_ps2 * regs_ps2.tp())^-1) * sum(regs_ps2 * (ps - cr[2]*eps[2]));
	eta0 = arr_T2[0];
	for(int i=0; i<9; i++) { T[i/3][i%3] = arr_T2[i+1]; }
	*/

	for(int i=0; i<4; i++)
	{
		matrix<10,1,valarray<double>> regs_ps2;
		//regs_nu[0] = nu*0 + 1;
		regs_ps2[0] = dt;
		//for(int i=0; i<2; i++) { regs_ps2[i+1] = eps[i]; }
		for(int i=0; i<9; i++) { regs_ps2[i+1] = eps[i/3] * cr[i%3]; }
		//for(int i=0; i<2; i++) { regs_ps2[i+1] = eps[i] * cr[i]; }
		//for(int i=0; i<2; i++) { regs_ps2[i+3] = eps[i] * cr[2]; }
		auto arr_T2 = (sum(regs_ps2 * regs_ps2.tp())^-1) * sum(regs_ps2 * (ps));// - eps[2]*cr[2]));//
		eta0 = arr_T2[0];
		//matrix<2,1> T0;
		//for(int i=0; i<2; i++) T0[i] = arr_T2[i+1];
		//cout<<"T0 = "<<T0<<endl;
		for(int i=0; i<9; i++) { T[i/3][i%3] = arr_T2[i+1]; }
		//for(int i=0; i<2; i++) { T[i][i] = arr_T2[i+1]; }
		//for(int i=0; i<2; i++) { T[i][2] = arr_T2[i+3]; }
		//T[2][2] = 1;
		cout<<"Localization "<<i + 1<<"...\nT = "<<T<<endl;
		cout<<"eta0 = "<<eta0<<endl;
		if(i==3) break;

		auto xi_mod = eye<3>();
		xi_mod[2] = T[2];
		cout<<"xi_mod = "<<xi_mod<<endl;
		Xi = xi_mod * Xi;
		xi0 = xi_mod * xi0;
		cr[2] = xi0[2] + Xi[2]*tx;

		T = T * (xi_mod^-1);
		cout<<"Post-adjustment "<<i + 1<<"...\nT = "<<T<<endl;
	}

	double sig1 = 0;
	sig0 = -0.117945;
	matrix<13,1,valarray<double>> regs_nu;
	//regs_nu[0] = nu*0 + 1;
	regs_nu[0] = dt;
	//for(int i=0; i<9; i++) { regs_nu[i+1] = eps[i/3] * cr[i%3]; }
	for(int i=0; i<3; i++) regs_nu[i+1] = eps[i];
	for(int i=0; i<9; i++) { regs_nu[i+4] = eps[i/3] * cr[i%3]; }
	//for(int i=0; i<2; i++) { regs_nu[i+4] = eps[i] * cr[i]; }
	//for(int i=0; i<2; i++) { regs_nu[i+6] = eps[i] * cr[2]; }
	auto arr_V = (sum(regs_nu * regs_nu.tp())^-1) * sum(regs_nu * (nu));
	sig0 = arr_V[0];
	//V = arr_V[0] * T;
	for(int i=0; i<3; i++) a[i] = arr_V[i+1];
	for(int i=0; i<9; i++) { V[i/3][i%3] = arr_V[i+4]; }
	//for(int i=0; i<9; i++) { V[i/3][i%3] = arr_V[i+1]; }
	//for(int i=0; i<2; i++) { V[i][i] = arr_V[i+4]; }
	//for(int i=0; i<2; i++) { V[i][2] = arr_V[i+6]; }

	cout<<"sig0 = "<<sig0<<endl;
	cout<<"sig1 = "<<sig1<<endl;
	//cout<<"arr_V[0] = "<<arr_V[0]<<endl;
	cout<<"a = "<<a<<endl;
	cout<<"V = "<<V<<endl;

	fstream f;
	f.open((path + "/spin_model_2").data(), fstream::out | fstream::trunc);
	f<<sig0<<"\n"<<endl;
	f<<a<<endl;
	f<<V<<endl;
	f<<eta0<<"\n"<<endl;
	f<<T<<endl;
	f.close();
	f.open((path + "/optics_model").data(), fstream::out | fstream::trunc);
	f<<xi0<<endl;
	f<<Xi<<endl;
	f.close();

	//M = V - (sig0/eta0)*T;

	valarray<double> res_nu = nu - (sig0*dt + eps*a + eps*V*cr);//(sig0/eta0);//arr_V.tp()*regs_nu;//
	double std_res_nu = pow(pow(res_nu,2).sum() / res_nu.size(),0.5);
	//valarray<double> res_ps = ps - (eta0*dt + eps*T*cr);
	valarray<double> res_ps = ps - (eta0*dt + eps*T*cr);
	double std_res_ps = pow(pow(res_ps,2).sum() / res_ps.size(),0.5);

	valarray<double> mag_O = pow(cr.tp() * cr,0.5);
	valarray<double> mag_M = pow((xi0 + Xi*tx).tp() * (xi0 + Xi*tx),0.5);
	valarray<double> res_mag = mag_O - mag_M;
	double std_res_mag = pow(pow(res_mag,2).sum()/ res_mag.size(),0.5);

	if(out)
	{
		TCanvas* cpts = new TCanvas("cpts","cpts",2800,1000);

		plot<double> pts(nu - sig0*dt,//[dt == dt.max()/2],//
										 "Spin tunes: Actual vs Linear Map","Data Point Index",
									 "#nu_{s} - #sigma_{0}#delta");
		plot<double> stm(eps*a + eps*V*cr);//arr_V.tp()*regs_nu);//(sig0*dt + (sig0*dt);//
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
		cpts->SaveAs((path + "/error_plot_2.pdf(").data(),"pdf");

		TCanvas* cpts2 = new TCanvas("cpts2","cpts2",2800,1000);
		plot<double> pts2(ps - eta0*dt,//[dt == dt.max()/2],//
											"Phase slips: Actual vs Linear Map","Data Point Index",
										"#frac{#DeltaT}{T} - #eta_{0}#delta");
		//plot<double> psm((arr_T2.tp() * regs_ps2) + eps[2]*cr[2]);//[dt == dt.max()/2]);//
		plot<double> psm(eps*T*cr);//[dt == dt.max()/2]);//
		pts2->Draw("AP");
		pts2->SetMarkerStyle(3);
		pts2->SetLineColor(0);
		psm->Draw("P same");
		psm->SetMarkerColor(2);
		psm->SetMarkerStyle(2);
		psm->SetLineColor(0);
		leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts2.get(),"Data Point (Measured Values)");
		leg->AddEntry(psm.get(),"Phase Slip Model (Calculated Values)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cpts2->SetGrid();
		cpts2->SaveAs((path + "/error_plot_2.pdf)").data(),"pdf");

		TCanvas* ceh = new TCanvas("ceh","ceh",1024,768);

		TH1F* hnu = new TH1F("hnu",
			"Distribution of residuals of spin-tune fit;#nu_{simulated} - #nu_{model};Frequency",
			51,-5* std_res_nu,5* std_res_nu);
		for(int i=0; i<res_nu.size(); i++) hnu->Fill(res_nu[i]);
		//hnu->Fit("gaus");
		hnu->Draw();
		ceh->SaveAs((path + "/error_dist_2.pdf(").data(),"pdf");

		TCanvas* ceh2 = new TCanvas("ceh2","ceh2",1024,768);
		TH1F* hps = new TH1F("hps",
			"Distribution of residuals of travel-time fit;#frac{#DeltaT}{T}_{simulated} - #frac{#DeltaT}{T}_{model};Frequency",
			51,-5* std_res_ps,5* std_res_ps);
		for(int i=0; i<res_ps.size(); i++) hps->Fill(res_ps[i]);
		//hps->Fit("gaus");
		hps->Draw();
		ceh2->SaveAs((path + "/error_dist_2.pdf)").data(),"pdf");

		TCanvas* cwmpts = new TCanvas("cwmpts","cwmpts",2800,1000);
		plot<double> wmpts(mag_O,//[abs(PE)<10],//
											 "Vector Magnitudes of Second Order Parameters: Actual vs Linear Map",
										"Data Point Index","#xi_{x}^{2} + #xi_{y}^{2} + #eta_{1}^{2}");
		plot<double> wmp(mag_M);//[abs(PE)<10]);//
		wmpts->Draw("AP");
		wmpts->SetMarkerStyle(3);
		wmpts->SetLineColor(0);
		wmp->Draw("P same");
		wmp->SetMarkerColor(2);
		wmp->SetMarkerStyle(2);
		wmp->SetLineColor(0);
		leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(wmpts.get(),"Data Point (Measured Values)");
		leg->AddEntry(wmp.get(),"Linear Map (Calculated Values)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cwmpts->SetGrid();
		cwmpts->SaveAs((path + "/wm_error_plot_2.pdf").data(),"pdf");

		TCanvas* ceh_wm = new TCanvas("ceh_wm","ceh_wm",1024,768);

		TH1F* hwme = new TH1F("h1","Distribution of residuals of Linear Model fit;Error;Frequency",51,-5* std_res_mag,5* std_res_mag);
		for(int i=0; i<res_mag.size(); i++) hwme->Fill(res_mag[i]);
		//hnu->Fit("gaus");
		hwme->Draw();
		ceh_wm->SaveAs((path + "/wm_error_dist_2.pdf").data(),"pdf");
	}

	M = V - (sig0/eta0)*T;
	//a = a - (sig0/eta0)*psa;
	cout<<"M = "<<M<<endl;
	cout<<"(M^-1) = "<<(M^-1)<<endl;
	xi = (M^-1)*(-a);
	cout<<"expected optical coordinates of optimized point:\n"<< xi <<endl;

	matrix<2,1> xi_o = {-1.35123,-4.527};
	double ts_o = -0.976576;

	auto op = (Xi^-1) * (xi - xi0);
	auto t2_o = (Xi.sub<2,2>(0,0)^-1) * (xi_o - xi0.sub<2,1>() - ts_o*Xi.sub<2,1>(0,2));
	matrix<3,1> t_o = {t2_o[0],t2_o[1],ts_o};
	cout<<"Optimized sextupole fields expected at: \n"<<op<<endl;
	cout<<"Calculated optimized optical settings: \n"<<xi0 + Xi*op<<endl;
	cout<<"Measured optimized sextupole fields: \n"<<t_o<<endl;
	cout<<"Measured optimized optical settings: \n"<<xi0 + Xi*t_o<<endl;

	cout<<"sig1 = "<<sig1<<endl;

	return make_tuple(a,M,xi);
}



