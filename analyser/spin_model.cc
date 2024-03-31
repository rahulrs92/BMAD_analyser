#include<analyser.hh>
#include<TH1F.h>

tuple<matrix<3,1>,matrix<3,3>,matrix<3,1>> analyser::spin_model(bool out)
{
	int N1 = 11;
	fstream f;
	f.open((path + "/eta").data(),fstream::in);
	if(f.is_open()) { f.close(); }
	else
	{
		auto Dc = proton::line(path + "/rc", NS);
		auto pt = lv1::pt(path + "/rc",Dc);
		f.open((path + "/eta").data(),fstream::out | fstream::trunc);
		for(int i=1; i<= NS; i++)
		{
			string path2 = path + "/r" + to_string(i);
			auto D = proton::line(path2,N1,out);
			auto eta = lv1::eta(path2,D,pt,out);
			f<<D.col("#xi_{x}")[0]<<"\t"<<D.col("#xi_{y}")[0]<<"\t"
				<<D.col("#epsilon_{x}")[0]<<"\t"<<D.col("#epsilon_{y}")[0]<<"\t"
				<<eta[0]<<"\t"<<eta[1]<<"\t"<<eta[2]<<endl;
		}
		f.close();
		f.open((path + "/eta_info").data(),fstream::out | fstream::trunc);
		f<<"#xi_{x}\n#xi_{y}\n#epsilon_{x}\n#epsilon_{y}\n"
			<<"T\n#eta_{0}\n#eta_{1}"<<endl;
		f.close();
	}
	Data<double> Dw(path + "/eta");
	Dw.clean_nans();

	vector<double> v_ex, v_ey, v_dt, v_cx, v_cy, v_eta1, v_nu, v_ps;
	auto tr_eta1 = Dw.col("#eta_{1}");
	for(int i=0; i< NS; i++)
	{
		string path2 = path + "/r" + to_string(i+1);
		auto D = proton::line(path2,N1,out);
		auto tr_ex = D.col("#epsilon_{x}");
		auto tr_ey = D.col("#epsilon_{y}");
		auto tr_dt = D.col("#delta_{0}");
		auto tr_cx = D.col("#xi_{x}");
		auto tr_cy = D.col("#xi_{y}");
		auto tr_nu = D.col("#nu_{s}");
		auto tr_ps = D.col("T");
		for(int j=0; j<N1; j++)
		{
			v_ex.push_back(tr_ex[j]);
			v_ey.push_back(tr_ey[j]);
			v_dt.push_back(tr_dt[j]);
			v_cx.push_back(tr_cx[j]);
			v_cy.push_back(tr_cy[j]);
			v_nu.push_back(tr_nu[j]);
			v_ps.push_back(tr_ps[j]);
			v_eta1.push_back(tr_eta1[i]);
			if(isnan(tr_ps[j]) || isnan(tr_nu[j]))
			{
				v_ex.pop_back();
				v_ey.pop_back();
				v_dt.pop_back();
				v_cx.pop_back();
				v_cy.pop_back();
				v_nu.pop_back();
				v_ps.pop_back();
				v_eta1.pop_back();
			}
		}
	}
	matrix<3,1,valarray<double>> cr;
	matrix<1,3,valarray<double>> eps;

	eps[0] = valarray<double>(v_ex.data(),v_ex.size());
	eps[1] = valarray<double>(v_ey.data(),v_ey.size());
	eps[2] = pow(valarray<double>(v_dt.data(),v_dt.size()),2);
	cr[0] = valarray<double>(v_cx.data(),v_cx.size());
	cr[1] = valarray<double>(v_cy.data(),v_cx.size());
	cr[2] = valarray<double>(v_eta1.data(),v_eta1.size());
	valarray<double> nu(v_nu.data(),v_nu.size());
	valarray<double> ps(v_ps.data(),v_ps.size());
	valarray<double> dt(v_dt.data(),v_dt.size());

	double nu0 = 0;
	matrix<13,1,valarray<double>> regs_nu;
	matrix<10,1,valarray<double>> regs_ps;
	//regs_nu[0] = nu*0 + 1;
	regs_nu[0] = dt; regs_ps[0] = dt;
	for(int i=0; i<3; i++) regs_nu[i+1] = eps[i];
	for(int i=0; i<9; i++) { regs_nu[i+4] = eps[i/3] * cr[i%3]; regs_ps[i+1] = eps[i/3] * cr[i%3]; }
	auto arr_V = (sum(regs_nu * regs_nu.tp())^-1) * sum(regs_nu * nu);
	auto arr_T = (sum(regs_ps * regs_ps.tp())^-1) * sum(regs_ps * ps);
	//nu0 = arr_V[0];
	sig0 = arr_V[0]; eta0 = arr_T[0];
	for(int i=0; i<3; i++) a[i] = arr_V[i+1];
	for(int i=0; i<9; i++) { V[i/3][i%3] = arr_V[i+4]; T[i/3][i%3] = arr_T[i+1]; }

	cout<<"sig0 = "<<sig0<<endl;
	cout<<"a = "<<a<<endl;
	cout<<"V = "<<V<<endl;
	cout<<"T = "<<T<<endl;

	f.open((path + "/spin_model").data(), fstream::out | fstream::trunc);
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
	if(out)
	{
		TCanvas* cpts = new TCanvas("cpts","cpts",1900,500);

		plot<double> pts(nu[dt == dt.max()/2],//,//
										 "Spin tunes: Actual vs Linear Map","Data Point Index","#nu_{s}");
		plot<double> wmp((sig0*dt + eps*a + eps*V*cr)[dt == dt.max()/2]);//);//
		pts->Draw("AP");
		pts->SetMarkerStyle(3);
		pts->SetLineColor(0);
		wmp->Draw("P same");
		wmp->SetMarkerColor(2);
		wmp->SetMarkerStyle(2);
		wmp->SetLineColor(0);
		auto leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts.get(),"Data Point (Measured Values)");
		leg->AddEntry(wmp.get(),"Spin Tune Model (Calculated Values)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cpts->SetGrid();
		cpts->SaveAs((path + "/error_plot.pdf(").data(),"pdf");

		plot<double> pts2(ps[dt == dt.max()/2],//,//
											"Phase slips: Actual vs Linear Map","Data Point Index","#frac{#DeltaT}{T}");
		plot<double> wmp2((eta0*dt + eps*T*cr)[dt == dt.max()/2]);//);//
		pts2->Draw("AP");
		pts2->SetMarkerStyle(3);
		pts2->SetLineColor(0);
		wmp2->Draw("P same");
		wmp2->SetMarkerColor(2);
		wmp2->SetMarkerStyle(2);
		wmp2->SetLineColor(0);
		leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts2.get(),"Data Point (Measured Values)");
		leg->AddEntry(wmp2.get(),"Phase Slip Model (Calculated Values)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cpts->SetGrid();
		cpts->SaveAs((path + "/error_plot.pdf)").data(),"pdf");

		TCanvas* ceh = new TCanvas("ceh","ceh",600,400);

		TH1F* hnu = new TH1F("hnu",
			"Distribution of residuals of spin-tune fit;#nu_{simulated} - #nu_{model};Frequency",
			51,-5* std_res_nu,5* std_res_nu);
		for(int i=0; i<nu.size(); i++) hnu->Fill(res_nu[i]);
		//hnu->Fit("gaus");
		hnu->Draw();
		ceh->SaveAs((path + "/error_dist.pdf(").data(),"pdf");

		TH1F* hps = new TH1F("hps",
			"Distribution of residuals of travel-time fit;#frac{#DeltaT}{T}_{simulated} - #frac{#DeltaT}{T}_{model};Frequency",
			51,-5* std_res_ps,5* std_res_ps);
		for(int i=0; i<nu.size(); i++) hps->Fill(res_ps[i]);
		//hps->Fit("gaus");
		hps->Draw();
		ceh->SaveAs((path + "/error_dist.pdf)").data(),"pdf");
	}

	M = V - (sig0/eta0)*T;
	//a = a - (sig0/eta0)*psa;
	cout<<"M = "<<M<<endl;
	xi = -(M^-1)*a;
	cout<<"expected optical coordinates of optimized point:\n"<< xi <<endl;

	cout<<"nu0 = "<<nu0<<endl;

	return make_tuple(a,M,xi);
}



