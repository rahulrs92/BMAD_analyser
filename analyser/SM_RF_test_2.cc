#include<analyser.hh>
#include <TMatrixD.h>

#include<TH3D.h>

void analyser::SM_RF_test2(Data<double>& Dc,bool out)
{
	Data<double> Dt(path + "/rtc/xi");

	valarray<double> nu = Dc.col("#nu_{s}");
	matrix<1,3,valarray<double>> eps = {Dc.col("#epsilon_{x}"), Dc.col("#epsilon_{y}"),pow(Dc.col("#delta_{a}"),2)/2};
	matrix<1,3,valarray<double>> epf = {eps[0], eps[1], eps[2]};//,pow(ep[0]*ep[2],0.5)};
	matrix<3,1,valarray<double>> t = {Dt.col(1),Dt.col(2),Dt.col(3)};

	matrix<3,1,valarray<double>> xi1 = xi0 + Xi*t;

	matrix<12,1,valarray<double>> ep_xi;
	for(int i=0; i<3; i++) ep_xi[i] = eps[i];
	for(int i=0; i<9; i++) ep_xi[i+3] = eps[i/3] * xi1[i%3];
	auto arr_M = (sum(ep_xi*ep_xi.tp())^-1) * sum(ep_xi * nu);
	matrix<3,3> M1;
	matrix<3,1> a1;
	for(int i=0; i<3; i++) a1[i] = arr_M[i];
	for(int i=0; i<9; i++) M1[i/3][i%3] = arr_M[i+3];

	cout<<"mean RMS error of RF-off estimates: "
			<<pow(pow(nu - eps *a - eps *((V - (sig0/eta0)*T)*xi1),2).sum(),0.5) / nu.size()<<endl;
	cout<<"mean RMS error of RF-on estimates: "
			<<pow(pow(nu - epf*a1 - eps *(M1* xi1),2).sum(),0.5) / nu.size()<<endl;

	auto op_off = ((V - (sig0/eta0)*T)^-1) * (-a);
	auto op_on = (M1^-1) * (-a1);
	//cout<<"optimized point (RF-off Estimate): \n"<<op_off<<endl;
	cout<<"optimized point (RF-on Estimate): \n"<<op_on<<endl;
	//cout<<"optimized fields (RF-off Estimate): \n"<<(Xi^-1) * (op_off - xi0)<<endl;
	cout<<"optimized fields (RF-on Estimate): \n"<<(Xi^-1) * (op_on - xi0)<<endl;

	valarray<double> res_nu_off = nu - (eps*a + eps*M*xi1);
	double std_res_nu_off = pow(pow(res_nu_off,2).sum()/ res_nu_off.size(),0.5);
	valarray<double> res_nu_on = nu - (eps*a1 + eps*M1*xi1);
	double std_res_nu_on = pow(pow(res_nu_on,2).sum()/ res_nu_on.size(),0.5);

	if(out)
	{
		TCanvas* cpts = new TCanvas("cpts","cpts",1900,500);

		plot<double> pts(nu,"Spin tunes: Actual vs RF-off predictions","Data Point Index","#nu_{s}");
		plot<double> wmp(eps *a + eps *((V - (sig0/eta0)*T)*xi1));
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

		plot<double> pts2(nu,"Spin tunes: Actual vs RF-on predictions","Data Point Index","#nu_{s}");
		plot<double> wmp2(epf*a1 + eps *(M1* xi1));
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
		cpts->SaveAs((path + "/rf_error_plot.pdf)").data(),"pdf");

		/*
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
		*/

		auto ceh = new TCanvas("creh","creh",600,400);
		TH1F* hnu_off = new TH1F("hnu_off",
			"Distribution of residuals of spin-tune fit;#nu_{simulated} - #nu_{model_rf_off};Frequency",
			51,-5* std_res_nu_off,5* std_res_nu_off);
		for(int i=0; i<nu.size(); i++) hnu_off->Fill(res_nu_off[i]);
		//hnu->Fit("gaus");
		hnu_off->Draw();
		ceh->SaveAs((path + "/rf_error_dist.pdf(").data(),"pdf");
		TH1F* hnu_on = new TH1F("hnu_on",
			"Distribution of residuals of spin-tune fit;#nu_{simulated} - #nu_{model_rf_on};Frequency",
			51,-5* std_res_nu_on,5* std_res_nu_on);
		for(int i=0; i<nu.size(); i++) hnu_on->Fill(res_nu_on[i]);
		//hnu->Fit("gaus");
		hnu_on->Draw();
		ceh->SaveAs((path + "/rf_error_dist.pdf)").data(),"pdf");
	}
}


