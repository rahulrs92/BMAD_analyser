#include<analyser.hh>

#include<TH3D.h>
/*
void lv2::spin_model_RF(Data<double>& D,string path,bool out)
{
	valarray<bool> sel = D.col("t_{s}") == 0.0;

	matrix<3,1,valarray<double>> t = {D.col("t_{f}")[sel],D.col("t_{d}")[sel],D.col("t_{s}")[sel]};
	valarray<double> nu = D.col("#nu_{s}")[sel];
	valarray<double> dta = D.col("#delta_{a}")[sel];

	matrix<1,3,valarray<double>> eps = {D.col("#epsilon_{x}")[sel], D.col("#epsilon_{y}")[sel],pow(dta,2)/2};

	matrix<2,1,valarray<double>> xi = {D.col("#xi_{x}")[sel], D.col("#xi_{y}")[sel]};

	matrix<9,1,valarray<double>> regs_M;
	for(int i=0; i<3; i++) regs_M[i] = (eps[i]);
	cout<<"YAY!"<<endl;
	for(int i=0; i<6; i++) regs_M[i+3] = (eps[i/2] * xi[i%2]);
	auto arr_M = (sum(regs_M * regs_M.tp())^-1) * sum(regs_M * nu);
	matrix<3,2> M1;
	matrix<3,1> a1;
	for(int i=0; i<3; i++) a1[i] = arr_M[i];
	for(int i=0; i<6; i++) M1[i/2][i%2] = arr_M[i+3];

	cout<<"a1 = "<<a1<<endl;
	cout<<"M1 = "<<M1<<endl;

	cout<<"b = "<<b<<endl;
	cout<<"M = "<<M<<endl;

	cout<<"mean RMS error of RF-on estimates: "
			<<pow(pow(nu - eps*a1 - eps*M1*xi,2).sum(),0.5) / nu.size()<<endl;

	//auto op_on = (M1^-1) * (-a1);
	auto op_on = ((M1.tp() * M1)^-1) * M1.tp() * (-a1);
	//cout<<"optimized point (RF-off Estimate): \n"<<op_off<<endl;
	cout<<"optimized point (RF-on Estimate): \n"<<op_on<<endl;
	matrix<1,3> eps0 = {eps[0][0], eps[1][0], eps[2][0]};
	cout<<"spin tune at optimized point: "<<(eps0*a1 + eps0*M1*op_on)<<endl;
	//cout<<"optimized fields (RF-off Estimate): \n"<<(Xi^-1) * (op_off - xi0)<<endl;
	//cout<<"optimized fields (RF-on Estimate): \n"<<(Xi^-1) * (op_on - xi0)<<endl;
	//cout<<"optimized fields (RF-on Estimate) (new calc): \n"<<(A.tp()^-1) * (op_on - C0.tp())<<endl;

	valarray<double> res_nu_on = nu - (eps*a1 + eps*M1*xi);
	double std_res_nu_on = pow(pow(res_nu_on,2).sum()/ res_nu_on.size(),0.5);

	valarray<double> res_nu_off = nu - (eps*a + eps*M*xi);
	double std_res_nu_off = pow(pow(res_nu_off,2).sum()/ res_nu_off.size(),0.5);
	cout<<"res_nu_off = "<<res_nu_off.max()<<endl;
	cout<<"std_res_nu_off = "<<std_res_nu_off<<endl;

	if(out)
	{
		TCanvas* cpts = new TCanvas("cpts","cpts",1900,500);
		plot<double> pts(nu,"Spin tunes: Actual vs RF-off predictions","Data Point Index","#nu_{s}");
		plot<double> wmp(eps*b + eps*M*xi);
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
		plot<double> wmp2(eps*a1 + eps*M1*xi);
		pts2->Draw("AP");
		pts2->SetMarkerStyle(3);
		pts2->SetLineColor(0);
		wmp2->Draw("P same");
		wmp2->SetMarkerColor(2);
		wmp2->SetMarkerStyle(2);
		wmp2->SetLineColor(0);
		auto leg2 = new TLegend(0.3,0.9,0.7,0.8);
		leg2->AddEntry(pts2.get(),"Data Point (Measured Values)");
		leg2->AddEntry(wmp2.get(),"Spin Tune Model (Calculated Values) (RF on)");
		leg2->SetFillColorAlpha(10,0.2);
		leg2->Draw();
		cpts->SetGrid();
		cpts->SaveAs((path + "/rf_error_plot.pdf)").data(),"pdf");

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
*/
