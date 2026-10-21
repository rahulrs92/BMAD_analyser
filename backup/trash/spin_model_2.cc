#include<analyser.hh>
#include<TH1F.h>

tuple<matrix<3,1>,matrix<3,3>,matrix<3,1>> lv2::spin_model_2(double iQx,double iQy,bool out)
{
	int N1 = 10;
	vector<vector<double>> C(10,vector<double>(NS * N1,0.0));
	vector<double> Ci(10,0.0);
	vector<uint> nans;
	//vector<double> v_ex, v_ey, v_dt, v_cx, v_cy, v_nu, v_ps, v_tf, v_td, v_ts;
	if(f_tree)
	{
		TTree *xitree = file->Get<TTree>("xi");
		xitree->SetBranchAddress("tf",&Ci[0]);
		xitree->SetBranchAddress("td",&Ci[1]);
		xitree->SetBranchAddress("ts",&Ci[2]);
		TTree *ctree = lv1::grid(iQx,iQy,file,0);
		ctree->SetBranchAddress("#epsilon_{x}",&Ci[3]);
		ctree->SetBranchAddress("#epsilon_{y}",&Ci[4]);
		ctree->SetBranchAddress("#delta_{0}",&Ci[5]);
		ctree->SetBranchAddress("#xi_{x}",&Ci[6]);
		ctree->SetBranchAddress("#xi_{y}",&Ci[7]);
		ctree->SetBranchAddress("#nu_{s}",&Ci[8]);
		ctree->SetBranchAddress("T",&Ci[9]);
		for(int i=0; i<NS; i++)
		{
			xitree->GetEntry(i);
			for(int j=0; j<N1; j++)
			{
				ctree->GetEntry(i*N1 + j);
				if(isnan(Ci[8]) || isnan(Ci[9])) nans.push_back(i*N1 + j);
				else
				{
					for(int k=0; k<10; k++) C[k][i*N1 + j] = Ci[k];
				}
			}
		}
	}
	else
	{
		Data<double> Dw(path + "/xi");
		vector<vector<double>> v_tr;
		for(int i=0; i<3; i++) v_tr.push_back(Dw.col<vector<double>>(i+1));
		for(int i=0; i<7; i++) v_tr.push_back(Dw.col<vector<double>>(1));
		for(int i=0; i< NS; i++)
		{
			string path2 = path + "/r" + to_string(i+1);
			auto D = proton::line(path2,N1,out);
			v_tr[3] = D.col<vector<double>>("#epsilon_{x}");
			v_tr[4] = D.col<vector<double>>("#epsilon_{y}");
			v_tr[5] = D.col<vector<double>>("#delta_{0}");
			v_tr[6] = D.col<vector<double>>("#xi_{x}");
			v_tr[7] = D.col<vector<double>>("#xi_{y}");
			v_tr[8] = D.col<vector<double>>("#nu_{s}");
			v_tr[9] = D.col<vector<double>>("T");
			for(int j=0; j<N1; j++)
			{
				if(isnan(v_tr[8][j]) || isnan(v_tr[9][j])) nans.push_back(i*N1 + j);
				else
				{
					for(int k=0; k<3; k++) C[k][i*N1 + j] = v_tr[k][i];
					for(int k=3; k<10; k++) C[k][i*N1 + j] = v_tr[k][j];
				}
			}
		}
	}
	for(int i=0; i<nans.size(); i++) for(int k=0; k<10; k++) C[k].erase(C[k].begin() + nans[i]);

	matrix<3,1,valarray<double>> cr, tx;
	matrix<1,3,valarray<double>> eps;

	for(int i=0; i<3; i++) tx[i] = valarray<double>(C[i].data(),C[i].size());
	for(int i=0; i<2; i++) eps[i] = valarray<double>(C[i+3].data(),C[i+3].size());
	valarray<double> dt(C[5].data(),C[5].size());
	eps[2] = pow(dt,2);
	for(int i=0; i<2; i++) cr[i] = valarray<double>(C[i+6].data(),C[i+6].size());
	valarray<double> nu(C[8].data(),C[8].size());
	valarray<double> ps(C[9].data(),C[9].size());

	cout<<"dt0 = "; for(int i=0; i<100; i++) cout<<C[5][i]<<" "; cout<<endl;

	matrix<3,3> T1;
	matrix<3,1> b;
	matrix<13,1,valarray<double>> regs_ps;
	regs_ps[0] = dt;
	for(int i=0; i<3; i++) regs_ps[i+1] = eps[i];
	for(int i=0; i<9; i++) { regs_ps[i+4] = eps[i/3] * tx[i%3]; }
	auto arr_T = (sum(regs_ps * regs_ps.tp())^-1) * sum(regs_ps * ps);
	eta0 = arr_T[0];
	for(int i=0; i<3; i++) b[i] = arr_T[i+1];
	for(int i=0; i<9; i++) T1[i/3][i%3] = arr_T[i+4];

	matrix<4,1,valarray<double>> regs_wm = {valarray<double>(1.0,tx[0].size()),tx[0],tx[1],tx[2]};
	auto X1 = ((sum(regs_wm * regs_wm.tp())^-1) * sum(regs_wm * cr.sub<2,1>().tp())).tp();
	xi0.sub<2,1>() = X1.sub<2,1>();
	Xi.sub<2,3>() = X1.sub<2,3>(0,1);

	//calculate the null space of the 2x3 transformation matrix:
	Xi[2][2] = 1;
	Xi.sub<1,2>(2,0) = ((X1.sub<2,2>(0,1)^-1) * (-X1).sub<2,1>(0,3)).tp();

	T = T1 * (Xi^-1);
	Xi[2] *= T[2][2];
	T(2) /= T[2][2];
	cout<<"current xi0 = "<<xi0<<endl;
	xi0[2] = (T^-1)[2] * b;
	cout<<"new xi0 = "<<(T^-1) * b<<endl;

	//cr[2] = valarray<double>(v_cy.data(),v_cy.size());
	cr[2] = xi0[2] + Xi[2]*tx;

	matrix<10,1,valarray<double>> regs_ps2;
	//regs_nu[0] = nu*0 + 1;
	regs_ps2[0] = dt;
	for(int i=0; i<9; i++) { regs_ps2[i+1] = eps[i/3] * cr[i%3]; }
	auto arr_T2 = (sum(regs_ps2 * regs_ps2.tp())^-1) * sum(regs_ps2 * ps);
	eta0 = arr_T2[0];
	for(int i=0; i<9; i++) { T[i/3][i%3] = arr_T2[i+1]; }

	double nu0 = 0;
	matrix<13,1,valarray<double>> regs_nu;
	//regs_nu[0] = nu*0 + 1;
	regs_nu[0] = dt;
	for(int i=0; i<3; i++) regs_nu[i+1] = eps[i];
	for(int i=0; i<9; i++) { regs_nu[i+4] = eps[i/3] * cr[i%3]; }
	auto arr_V = (sum(regs_nu * regs_nu.tp())^-1) * sum(regs_nu * nu);
	sig0 = arr_V[0];
	for(int i=0; i<3; i++) a[i] = arr_V[i+1];
	for(int i=0; i<9; i++) { V[i/3][i%3] = arr_V[i+4]; }

	cout<<"sig0 = "<<sig0<<endl;
	cout<<"eta0 = "<<eta0<<endl;
	cout<<"a = "<<a<<endl;
	cout<<"V = "<<V<<endl;
	cout<<"T = "<<T<<endl;

	fstream f;
	f.open((path + "/spin_model_2").data(), fstream::out | fstream::trunc);
	f<<sig0<<"\n"<<endl;
	f<<a<<endl;
	f<<V<<endl;
	f<<eta0<<"\n"<<endl;
	f<<T<<endl;
	f.close();

	M = V - (sig0/eta0)*T;

	valarray<double> res_nu = nu - (sig0*dt + eps*a + eps*V*cr);
	double std_res_nu = pow(pow(res_nu,2).sum()/ res_nu.size(),0.5);
	valarray<double> res_ps = ps - (eta0*dt + eps*T*cr);
	double std_res_ps = pow(pow(res_ps,2).sum()/ res_ps.size(),0.5);
	/*
	ch2 /= pow(nu,2).sum() / (nu.size() - 1);
	ch2 /= nu.size() - 13;
	*/
	valarray<double> mag_O = pow(cr.tp() * cr,0.5);
	valarray<double> mag_M = pow((xi0 + Xi*tx).tp() * (xi0 + Xi*tx),0.5);
	//auto PE = (Ci - (Xi*Fi).tp()) * (Ci.tp() - (Xi*Fi));
	//PE = pow(PE,0.5);

	valarray<double> PE = (cr - (xi0 + Xi*tx)).tp() * ((xi0 + Xi*tx) / mag_M);
	PE[mag_M == 0] = 0;
	cout<<"PE = "<<valarray<double>(PE[slice(0,100,1)])<<endl;

	double std_PE = pow(pow(PE,2).sum()/ PE.size(),0.5);

	if(out)
	{
		TCanvas* cpts = new TCanvas("cpts","cpts",1900,500);

		plot<double> pts(nu,//[dt == dt.max()/2],//
										 "Spin tunes: Actual vs Linear Map","Data Point Index","#nu_{s}");
		plot<double> stm((sig0*dt + eps*a + eps*V*cr));//[dt == dt.max()/2]);//
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

		plot<double> pts2(ps,//[dt == dt.max()/2],//
											"Phase slips: Actual vs Linear Map","Data Point Index","#frac{#DeltaT}{T}");
		//plot<double> psm((eta0*dt + eps*b + eps*T1*tx)[dt == dt.max()/2]);//);//
		plot<double> psm((eta0*dt + eps*T*cr));//[dt == dt.max()/2]);//
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
		cpts->SetGrid();
		cpts->SaveAs((path + "/error_plot_2.pdf)").data(),"pdf");

		TCanvas* ceh = new TCanvas("ceh","ceh",600,400);

		TH1F* hnu = new TH1F("hnu",
			"Distribution of residuals of spin-tune fit;#nu_{simulated} - #nu_{model};Frequency",
			51,-5* std_res_nu,5* std_res_nu);
		for(int i=0; i<nu.size(); i++) hnu->Fill(res_nu[i]);
		//hnu->Fit("gaus");
		hnu->Draw();
		ceh->SaveAs((path + "/error_dist_2.pdf(").data(),"pdf");

		TH1F* hps = new TH1F("hps",
			"Distribution of residuals of travel-time fit;#frac{#DeltaT}{T}_{simulated} - #frac{#DeltaT}{T}_{model};Frequency",
			51,-5* std_res_ps,5* std_res_ps);
		for(int i=0; i<nu.size(); i++) hps->Fill(res_ps[i]);
		//hps->Fit("gaus");
		hps->Draw();
		ceh->SaveAs((path + "/error_dist_2.pdf)").data(),"pdf");

		TCanvas* cwmpts = new TCanvas("cwmpts","cwmpts",1900,500);
		plot<double> wmpts(mag_O[abs(PE)<10],"Vector Magnitudes of Second Order Parameters: Actual vs Linear Map",
										"Data Point Index","#xi_{x}^{2} + #xi_{y}^{2} + #eta_{1}^{2}");
		plot<double> wmp(mag_M[abs(PE)<10]);
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

		TCanvas* ceh_wm = new TCanvas("ceh_wm","ceh_wm",600,400);

		TH1F* hwme = new TH1F("h1","Distribution of residuals of Linear Model fit;Error;Frequency",51,-5* std_PE,5* std_PE);
		for(int i=0; i<PE.size(); i++) hwme->Fill(PE[i]);
		//hnu->Fit("gaus");
		hwme->Draw();
		ceh_wm->SaveAs((path + "/wm_error_dist_2.pdf").data(),"pdf");
	}

	M = V - (sig0/eta0)*T;
	//a = a - (sig0/eta0)*psa;
	cout<<"M = "<<M<<endl;
	xi = -(M^-1)*a;
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

	cout<<"nu0 = "<<nu0<<endl;

	return make_tuple(a,M,xi);
}




