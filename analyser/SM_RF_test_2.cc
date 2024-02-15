#include<analyser.hh>
#include <TMatrixD.h>

#include<TH3D.h>

void analyser::SM_RF_test2(Data<double>& Dc,Data<double>& De,bool out)
{
	valarray<double> nu_e = De.col("#nu_{s}");
	matrix<1,3,valarray<double>> ep_e = {De.col("#epsilon_{x}"), De.col("#epsilon_{y}"),pow(De.col("#delta_{a}"),2)/2};
	//for(int i=0; i<ep_e[0].size(); i++) ep_e[0][i] = ep_e[0][i/4];
	matrix<1,3,valarray<double>> ep_f = {ep_e[0],ep_e[1],ep_e[2]};//,pow(ep_e[0]*ep_e[2],0.5)};

	auto a1 = (sum(ep_f.tp()*ep_f)^-1) * sum(ep_f.tp() * nu_e);

	Data<double> Dt(path + "/rtc/xi");

	valarray<double> nu = Dc.col("#nu_{s}");
	matrix<1,3,valarray<double>> ep = {Dc.col("#epsilon_{x}"), Dc.col("#epsilon_{y}"),pow(Dc.col("#delta_{a}"),2)/2};
	matrix<1,3,valarray<double>> epf = {ep[0],ep[1],ep[2]};//,pow(ep[0]*ep[2],0.5)};
	matrix<3,1,valarray<double>> t = {Dt.col(1),Dt.col(2),Dt.col(3)};

	matrix<3,1,valarray<double>> xi1 = xi0 + Xi*t;

	/*
	matrix<12,1,valarray<double>> ep_xi;
	for(int i=0; i<12; i++) ep_xi[i] = epf[i/3] * xi1[i%3];
	//auto arr_M = (sum(ep_xi*ep_xi.tp())^-1) * sum(ep_xi * nu);
	auto big_M = sum(ep_xi*ep_xi.tp());
	cout<<big_M<<endl;
	cout<<"calculation starts!"<<endl;
	TArrayD t_h(144);
	for(int i=0; i<144; i++) t_h[i] = big_M[i/12][i%12];
	TMatrixD t_ep_xi; t_ep_xi.Use(12,12,t_h.GetArray());
	t_ep_xi.Invert();
	matrix<12,12> inv_big_M;
	for(int i=0; i<144; i++) inv_big_M[i/12][i%12] = t_h[i];
	cout<<"calculation ends!"<<endl;
	cout<<inv_big_M<<endl;
	auto arr_M = (inv_big_M) * sum(ep_xi * (nu - epf*a1));
	matrix<4,3> M1;
	for(int i=0; i<12; i++) M1[i/3][i%3] = arr_M[i];
	*/

	matrix<9,1,valarray<double>> ep_xi;
	for(int i=0; i<9; i++) ep_xi[i] = ep[i/3] * xi1[i%3];
	auto arr_M = (sum(ep_xi*ep_xi.tp())^-1) * sum(ep_xi * (nu - epf*a1));
	matrix<3,3> M1;
	for(int i=0; i<9; i++) M1[i/3][i%3] = arr_M[i];

	cout<<"mean RMS error of RF-off estimates: "
			<<pow(pow(nu - ep*a - ep*((V - (sig0/eta0)*T)*xi1),2).sum(),0.5) / nu.size()<<endl;
	cout<<"mean RMS error of RF-on estimates: "
			<<pow(pow(nu - epf*a1 - ep*(M1* xi1),2).sum(),0.5) / nu.size()<<endl;
	cout<<"mean RMS error of RF-on estimates at optical origin: "
			<<pow(pow(nu_e - ep_f*a1,2).sum(),0.5) / nu_e.size()<<endl;

	auto op_off = ((V - (sig0/eta0)*T)^-1) * (-a);
	auto op_on = (M1^-1) * (-a1);
	cout<<"optimized point (RF-off Estimate): \n"<<op_off<<endl;
	cout<<"optimized point (RF-on Estimate): \n"<<op_on<<endl;
	cout<<"optimized fields (RF-off Estimate): \n"<<(Xi^-1) * (op_off - xi0)<<endl;
	cout<<"optimized fields (RF-on Estimate): \n"<<(Xi^-1) * (op_on - xi0)<<endl;



	if(out)
	{
		TCanvas* cpts = new TCanvas("cpts","cpts",1900,500);

		plot<double> pts(nu,"Spin tunes: Actual vs RF-off predictions","Data Point Index","#nu_{s}");
		plot<double> wmp(ep*a + ep*((V - (sig0/eta0)*T)*xi1));
		pts->Draw("AP");
		pts->SetMarkerStyle(3);
		pts->SetLineColor(0);
		wmp->Draw("P same");
		wmp->SetMarkerColor(2);
		wmp->SetMarkerStyle(2);
		wmp->SetLineColor(0);
		auto leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts.get(),"Data Point (Measured Values)");
		leg->AddEntry(wmp.get(),"Spin Tune Model (Calculated Values) (RF off)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cpts->SetGrid();
		cpts->SaveAs((path + "/rf_error_plot.pdf(").data(),"pdf");

		plot<double> ptse(nu_e,"Spin tunes: Actual vs RF-on predictions (optical origin)","Data Point Index","#nu_{s}");
		plot<double> wmpe(ep_f*a1);
		ptse->Draw("AP");
		ptse->SetMarkerStyle(3);
		ptse->SetLineColor(0);
		wmpe->Draw("P same");
		wmpe->SetMarkerColor(2);
		wmpe->SetMarkerStyle(2);
		wmpe->SetLineColor(0);
		auto lege = new TLegend(0.3,0.9,0.7,0.8);
		lege->AddEntry(pts.get(),"Data Point (Measured Values)");
		lege->AddEntry(wmp.get(),"Spin Tune Model (Calculated Values) (RF on)");
		lege->SetFillColorAlpha(10,0.2);
		lege->Draw();
		cpts->SetGrid();
		cpts->SaveAs((path + "/rf_error_plot.pdf").data(),"pdf");

		plot<double> pts2(nu,"Spin tunes: Actual vs RF-on predictions","Data Point Index","#nu_{s}");
		plot<double> wmp2(epf*a1 + ep*(M1* xi1));
		pts2->Draw("AP");
		pts2->SetMarkerStyle(3);
		pts2->SetLineColor(0);
		wmp2->Draw("P same");
		wmp2->SetMarkerColor(2);
		wmp2->SetMarkerStyle(2);
		wmp2->SetLineColor(0);
		auto leg2 = new TLegend(0.3,0.9,0.7,0.8);
		leg2->AddEntry(pts.get(),"Data Point (Measured Values)");
		leg2->AddEntry(wmp.get(),"Spin Tune Model (Calculated Values) (RF on)");
		leg2->SetFillColorAlpha(10,0.2);
		leg2->Draw();
		cpts->SetGrid();
		cpts->SaveAs((path + "/rf_error_plot.pdf").data(),"pdf");

		plot<double> pts3(xi1[2],"2° Phase-slip factors: RF-off vs Wave measurements","Data Point Index","#eta_{1}");
		plot<double> wmp3(Dc.col("#eta_{1}"));
		pts3->Draw("AP");
		pts3->SetMarkerStyle(3);
		pts3->SetLineColor(0);
		wmp3->Draw("P same");
		wmp3->SetMarkerColor(2);
		wmp3->SetMarkerStyle(2);
		wmp3->SetLineColor(0);
		auto leg3 = new TLegend(0.3,0.9,0.7,0.8);
		leg3->AddEntry(pts.get(),"RF-off measurements");
		leg3->AddEntry(wmp.get(),"Wave-based measurements");
		leg3->SetFillColorAlpha(10,0.2);
		leg3->Draw();
		cpts->SetGrid();
		cpts->SaveAs((path + "/rf_error_plot.pdf)").data(),"pdf");

		plot<double> pts4(xi1[2],pow(nu - ep*a - ep*((V - (sig0/eta0)*T)*xi1),2),
											"Correlation of errors","#eta_{1}","#chi^{2} error on #nu_{s}");
		pts4->Draw("AP");
		pts4->SetMarkerStyle(3);
		pts4->SetLineColor(0);
		cpts->SetGrid();
		cpts->SaveAs((path + "/error_correlations.pdf(").data(),"pdf");

		plot<double> pts5(xi1[0],pow(nu - ep*a - ep*((V - (sig0/eta0)*T)*xi1),2),
											"Correlation of errors","#xi_{x}","#chi^{2} error on #nu_{s}");
		pts5->Draw("AP");
		pts5->SetMarkerStyle(3);
		pts5->SetLineColor(0);
		cpts->SetGrid();
		cpts->SaveAs((path + "/error_correlations.pdf").data(),"pdf");

		plot<double> pts6(xi1[1],pow(nu - ep*a - ep*((V - (sig0/eta0)*T)*xi1),2),
											"Correlation of errors","#xi_{y}","#chi^{2} error on #nu_{s}");
		pts6->Draw("AP");
		pts6->SetMarkerStyle(3);
		pts6->SetLineColor(0);
		cpts->SetGrid();
		cpts->SaveAs((path + "/error_correlations.pdf").data(),"pdf");

		plot<double> pts7(pow(t.tp()*t,0.5),pow(nu - ep*a - ep*((V - (sig0/eta0)*T)*xi1),2),
											"Correlation of errors","(#xi_{x}^{2} + #xi_{y}^{2} + #eta_{1}^{2})^{#frac{1}{2}}",
											"#chi^{2} error on #nu_{s}");
		pts7->Draw("AP");
		pts7->SetMarkerStyle(3);
		pts7->SetLineColor(0);
		cpts->SetGrid();
		cpts->SaveAs((path + "/error_correlations.pdf)").data(),"pdf");

		TCanvas* ch1 = new TCanvas("ch1","ch1",1900,1000);
		//TH3D* H = new TH3D("H","Distribution of spin-tune errors;#xi_{x};#xi_{y};#eta_{1}",
											//20,-40,40,20,-30,30,20,-15,15);
		TH3D* H = new TH3D("H","Distribution of spin-tune errors;t_{r};t_{g};t_{b}",
												7,-0.09,0.09,7,-0.09,0.09,7,-0.09,0.09);
		for(int i=0; i<xi1[0].size(); i++)
			H->Fill(t[0][i],t[1][i],t[2][i],pow(nu - ep*a - ep*((V - (sig0/eta0)*T)*xi1),2)[i]);
		H->Draw("ISO");

		TCanvas* ch2 = new TCanvas("ch2","ch2",1900,1000);
		//TH3D* H = new TH3D("H","Distribution of spin-tune errors;#xi_{x};#xi_{y};#eta_{1}",
											//20,-40,40,20,-30,30,20,-15,15);
		TH3D* H2 = new TH3D("H","Distribution of spin-tunes;t_{r};t_{g};t_{b}",
												7,-0.09,0.09,7,-0.09,0.09,7,-0.09,0.09);
		for(int i=0; i<xi1[0].size(); i++)
			H2->Fill(t[0][i],t[1][i],t[2][i],pow(nu,4)[i]);
		H2->Draw("ISO");
	}
}


