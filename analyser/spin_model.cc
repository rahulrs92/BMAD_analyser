#include<analyser.hh>

//void lv2::spin_model(string path,bool out)
tuple<matrix<1,3>,matrix<3,1>,matrix<3,3>,matrix<3,1>> analyser::spin_model(bool out)
{
	int N = 343, N1 = 11;
	auto De = proton::line(path + "/re",275);
	auto Ds = proton::line(path + "/re/0e",11);

	valarray<double> dt = De.col("#delta_{0}");
	valarray<double> nu = De.col("#nu_{s}");
	valarray<double> ex = De.col("#epsilon_{x}");
	valarray<double> ey = De.col("#epsilon_{y}");
	valarray<double> ps = De.col("T");

	matrix<1,3,valarray<double>> ep = {ex,ey,pow(dt,2)};


	auto sig = sigma(path + "/re/0e",Ds);
	auto eta = lv1::eta(path + "/re/0e",Ds);
	a = (sum(ep.tp()*ep)^-1) * sum(ep.tp() * (nu - sig[1]*dt));
	psa = (sum(ep.tp()*ep)^-1) * sum(ep.tp() * (ps - eta[1]*dt));
	cout<<"a = "<<a<<endl;

	for(int i=0; i<3; i++)
	{
		sig = lv1::sigma(path + "/re/0e",Ds,{a[0],a[1]});
		eta = lv1::eta(path + "/re/0e",Ds,{psa[0],psa[1]});
		a = (sum(ep.tp()*ep)^-1) * sum(ep.tp() * (nu - sig[1]*dt));
		psa = (sum(ep.tp()*ep)^-1) * sum(ep.tp() * (ps - eta[1]*dt));
	}
	cout<<"a = "<<a<<endl;

	if(out)
	{
		TCanvas* cpts = new TCanvas("cpts","cpts",1900,500);
		plot<double> pts(nu,"Spin tunes: Actual vs Linear Map","Data Point Index","#nu_{s}");
		plot<double> wmp(sig[1]*dt + ep*a);
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
		cpts->SaveAs((path + "/re/ai_error_plot.pdf(").data(),"pdf");

		plot<double> pts2(ps,"Phase slip: Actual vs Linear Map","Data Point Index","#frac{#DeltaT}{T}");
		plot<double> wmp2(eta[1]*dt + ep*psa);
		pts2->Draw("AP");
		pts2->SetMarkerStyle(3);
		pts2->SetLineColor(0);
		wmp2->Draw("P same");
		wmp2->SetMarkerColor(2);
		wmp2->SetMarkerStyle(2);
		wmp2->SetLineColor(0);
		auto leg2 = new TLegend(0.3,0.9,0.7,0.8);
		leg2->AddEntry(pts.get(),"Data Point (Measured Values)");
		leg2->AddEntry(wmp.get(),"Phase Slip Model (Calculated Values)");
		leg2->SetFillColorAlpha(10,0.2);
		leg2->Draw();
		cpts->SetGrid();
		cpts->SaveAs((path + "/re/ai_error_plot.pdf)").data(),"pdf");
	}

	fstream f;
	f.open((path + "/eta").data(),fstream::in);
	if(f.is_open()) { f.close(); }
	else
	{
		auto Dc = proton::line(path + "/rc",49);
		auto pt = lv1::pt(path + "/rc",Dc);
		for(int i=1; i<=N; i++)
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

	vector<double> v_ex, v_ey, v_dt, v_cx, v_cy, v_eta1, v_nu, v_ps;
	auto tr_eta1 = Dw.col("#eta_{1}");
	for(int i=0; i<N; i++)
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
	ep[0] = valarray<double>(v_ex.data(),v_ex.size());
	ep[1] = valarray<double>(v_ey.data(),v_ey.size());
	ep[2] = pow(valarray<double>(v_dt.data(),v_dt.size()),2);
	cr[0] = valarray<double>(v_cx.data(),v_cx.size());
	cr[1] = valarray<double>(v_cy.data(),v_cx.size());
	cr[2] = valarray<double>(v_eta1.data(),v_eta1.size());
	nu = valarray<double>(v_nu.data(),v_nu.size());
	ps = valarray<double>(v_ps.data(),v_ps.size());
	dt = valarray<double>(v_dt.data(),v_dt.size());

	cout<<"dt = "<<valarray<double>(ps).sum()<<endl;


	matrix<9,1,valarray<double>> ep_xi;
	for(int i=0; i<9; i++) ep_xi[i] = ep[i/3] * cr[i%3];
	auto arr_V = (sum(ep_xi*ep_xi.tp())^-1) * sum(ep_xi * (nu - sig[1]*dt - ep*a));
	auto arr_T = (sum(ep_xi*ep_xi.tp())^-1) * sum(ep_xi * (ps - eta[1]*dt - ep*psa));
	for(int i=0; i<9; i++) { V[i/3][i%3] = arr_V[i]; T[i/3][i%3] = arr_T[i]; }

	cout<<"a = "<<a<<endl;
	cout<<"V = "<<V<<endl;
	cout<<"T = "<<T<<endl;

	if(out)
	{
		TCanvas* cpts = new TCanvas("cpts","cpts",1900,500);

		plot<double> pts(nu[dt == dt.max()/2],
										 "Spin tunes: Actual vs Linear Map","Data Point Index","#nu_{s}");
		plot<double> wmp((sig[1]*dt + ep*a + ep*(V* cr))[dt == dt.max()/2]);
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

		plot<double> pts2(ps[dt == dt.max()/2],
											"Phase slips: Actual vs Linear Map","Data Point Index","#frac{#DeltaT}{T}");
		plot<double> wmp2((eta[1]*dt + ep*psa + ep*(T* cr))[dt == dt.max()/2]);
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
	}

	sig0 = sig[1]; eta0 = eta[1];
	M = V - (sig0/eta0)*T;
	a = a - (sig0/eta0)*psa;
	cout<<"M = "<<M<<endl;
	xi = -(M^-1)*a;
	cout<<"expected optical coordinates of optimized point:\n"<< xi <<endl;

	return make_tuple(sum(ep),a,M,xi);
}



