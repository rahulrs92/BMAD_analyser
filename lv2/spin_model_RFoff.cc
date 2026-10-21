#include <lv2.hh>
#include<TH1F.h>

tuple<double,matrix<3,1>,matrix<3,2>> lv2::spin_model_RFoff(Data<double>& D,string path,bool out)
{
	fstream f(path + "/tt_matrix_info", fstream::out | fstream::trunc);
	f<<"#eta_{0}"<<"\n";
    f<<"u3[0]"<<"\n"<<"u3[1]"<<"\n";
	f<<"U3[0][0]"<<"\n"<<"U3[0][1]"<<"\n"<<"U3[0][2]"<<"\n";
	f<<"U3[1][0]"<<"\n"<<"U3[1][1]"<<"\n"<<"U3[1][2]"<<"\n";
	f.close();
    f.open(path + "/sm_matrix_info", fstream::out | fstream::trunc);
	f<<"#sigma_{0}"<<"\n";
    f<<"h3[0]"<<"\n"<<"h3[1]"<<"\n"<<"h3[2]"<<"\n";
	f<<"H3[0][0]"<<"\n"<<"H3[0][1]"<<"\n"<<"H3[0][2]"<<"\n";
	f<<"H3[1][0]"<<"\n"<<"H3[1][1]"<<"\n"<<"H3[1][2]"<<"\n";
	f<<"H3[2][0]"<<"\n"<<"H3[2][1]"<<"\n"<<"H3[2][2]"<<"\n";
	f.close();
    f.open(path + "/xi_out_info", fstream::out | fstream::trunc);
    f<<"t_{f}"<<"\n"<<"t_{d}"<<"\n"<<"t_{s}"<<"\n";
    f<<"#xi_{x}"<<"\n"<<"#xi_{y}"<<"\n"<<"#eta_{1}"<<"\n";
    f.close();

	//tie(h,H,ignore) = az.spin_model_2(ln,path);
    valarray<double> dt = D.col("#delta_{0}");
	matrix<1,3,valarray<double>> eps = {D.col("#epsilon_{x}"),D.col("#epsilon_{y}"),pow(dt,2)};
	matrix<3,1,valarray<double>> cr = {D.col("#xi_{x}"),D.col("#xi_{y}"),D.col("t_{s}")};
	matrix<3,1,valarray<double>> t = {D.col("t_{f}"),D.col("t_{d}"),D.col("t_{s}")};
	valarray<double> ps = D.col("T");
    valarray<double> nu = D.col("#nu_{s}");
    
	matrix<1,3> to_eta2;
	matrix<3,3> To_eta2;
	matrix<31,1,valarray<double>> regs_ps;
    regs_ps[0] = dt;
    for(int i=0; i<3; i++) regs_ps[i+1] = eps[i];
	for(int i=0; i<9; i++) regs_ps[i+4] = eps[i/3] * t[i%3];
	for(int i=0; i<9; i++) regs_ps[i+13] = eps[i/3] * t[i%3] * t[i%3];
	for(int i=0; i<9; i++) regs_ps[i+22] = eps[i/3] * t[i%3] * t[(i+1)%3];
	auto arr_U = (sum(regs_ps * regs_ps.tp())^-1) * sum(regs_ps * (ps));
    eta0 = arr_U[0];
    for(int i=0; i<3; i++) u[i] = arr_U[i+1];
    for(int i=0; i<3; i++) to_eta2[i] = arr_U[i+4+6];
	for(int i=0; i<3; i++) To_eta2[i][i] = arr_U[i+13+6];

	for(int i=0; i<3; i++) To_eta2[i][(i+1)%3] = arr_U[i+22+6];
	for(int i=0; i<3; i++) To_eta2[i][(i+2)%3] = To_eta2[(i+1)%3][i];
    
    cr[2] = u[2] + to_eta2*t + t.tp()*To_eta2*t;
    //cr[2] = u[2] + U3[2][0]*cr[0] + U3[2][1]*cr[1] + U3[2][2]*cr[2];
    //u[2] = 0; U3[2][0] = 0; U3[2][1] = 0; U3[2][2] = 1;
    
	
    matrix<10,1,valarray<double>> regs_ps2;
    regs_ps2[0] = dt;
    for(int i=0; i<3; i++) regs_ps2[i+1] = eps[i];
	for(int i=0; i<6; i++) regs_ps2[i+4] = eps[i/3] * cr[i%3];
	auto arr_U2 = (sum(regs_ps2 * regs_ps2.tp())^-1) * sum(regs_ps2 * (ps));
    eta0 = arr_U2[0];
    for(int i=0; i<3; i++) u[i] = arr_U2[i+1];
	for(int i=0; i<9; i++) to_eta2[i] = arr_U2[i+18+4];

	//valarray<double> res_ps = ps - (eta0*dt + eps*u + eps*U3*cr);
	valarray<double> res_ps = ps - arr_U.tp()*regs_ps;
	double std_res_ps = pow(pow(res_ps,2).sum() / res_ps.size(),0.5);
    
	f.open(path + "/tt_matrix", fstream::out | fstream::trunc);
	f<<eta0<<"\t"
    <<u[0]<<"\t"<<u[1]<<"\t"<<u[2]<<"\t"
    <<U3[0][0]<<"\t"<<U3[0][1]<<"\t"<<U3[0][2]<<"\t"
	<<U3[1][0]<<"\t"<<U3[1][1]<<"\t"<<U3[1][2]<<endl;
	f.close();

	f.open(path + "/xi_out", fstream::out | fstream::trunc);
	for(int i=0; i<t[0].size(); i++)
    {
        for(int j=0; j<3; j++) f<<t[j][i]<<"\t";
        for(int j=0; j<2; j++) f<<cr[j][i]<<"\t"; f<<cr[2][i]<<endl;
	}
	f.close();
    
	cout<<"eta0 = "<<eta0<<endl;
	cout<<"h = "<<h<<endl;
	cout<<"H = "<<H<<endl;
    
	matrix<13,1,valarray<double>> regs_nu;
    regs_nu[0] = dt;
    for(int i=0; i<3; i++) regs_nu[i+1] = eps[i];
	for(int i=0; i<9; i++) regs_nu[i+4] = eps[i/3] * cr[i%3];
	auto arr_H = (sum(regs_nu * regs_nu.tp())^-1) * sum(regs_nu * (nu));
    sig0 = arr_H[0];
    for(int i=0; i<3; i++) h[i] = arr_H[i+1];
	for(int i=0; i<9; i++) H3[i/3][i%3] = arr_H[i+4];
    
    f.open(path + "/sm_matrix", fstream::out | fstream::trunc);
	f<<sig0<<"\t"
    <<h[0]<<"\t"<<h[1]<<"\t"<<h[2]<<"\t"
    <<H3[0][0]<<"\t"<<H3[0][1]<<"\t"<<H3[0][2]<<"\t"
	<<H3[1][0]<<"\t"<<H3[1][1]<<"\t"<<H3[1][2]<<"\t"
	<<H3[2][0]<<"\t"<<H3[2][1]<<"\t"<<H3[2][2]<<endl;
	f.close();
    
    valarray<double> res_nu = nu - (sig0*dt + eps*h + eps*H3*cr);
	double std_res_nu = pow(pow(res_nu,2).sum() / res_nu.size(),0.5);
    
    
    matrix<3,3> Xi3;
    matrix<3,1> xi03;
	//windmap(D,path,out);
    tie(xi03,Xi3) = windmap(D,cr[2],path,out);

    W3 = H3 - (sig0/eta0)*U3;
	w = h - (sig0/eta0)*u;

    auto xi_opt = (W3^-1) * (-w);
    cout<<"xi_opt = "<<xi_opt<<endl;
    cout<<"t_opt = "<<(Xi3^-1) * (xi_opt-xi03)<<endl;

	if(out)
	{
		TCanvas* cpts_ps = new TCanvas("cpts_ps","cpts_ps",2800,1000);
		plot<double> pts_ps(ps,
								"Travel Time: Actual vs Linear Map","Data Point Index",
								"#frac{#DeltaT}{T}");// - #eta_{0}#delta");//
		//plot<double> psm((arr_T2.tp() * regs_ps2) + eps[2]*cr[2]);//[dt == dt.max()/2]);//
		plot<double> psm(ps - res_ps);
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
		cpts_ps->SaveAs((path + "/error_plot.pdf(").data(),"pdf");

		TCanvas* ceh_ps = new TCanvas("ceh_ps","ceh_ps",1024,768);
		TH1F* hps = new TH1F("hps",
			"Distribution of residuals of travel-time fit;#frac{#DeltaT}{T}_{simulated} - #frac{#DeltaT}{T}_{model};Frequency",
			51,-5* std_res_ps,5* std_res_ps);
		for(int i=0; i<res_ps.size(); i++) hps->Fill(res_ps[i]);
		//hps->Fit("gaus");
		hps->Draw();
		ceh_ps->SaveAs((path + "/error_dist.pdf(").data(),"pdf");

		TCanvas* cpts_nu = new TCanvas("cpts_nu","cpts_nu",2800,1000);
		plot<double> pts_nu(nu - sig0*dt,
								"Spin Tune: Simulation vs Model","Data Point Index",
								"#nu_{s} - #sigma_{0}#delta");
		//plot<double> psm((arr_T2.tp() * regs_ps2) + eps[2]*cr[2]);//[dt == dt.max()/2]);//
		plot<double> num(eps*h + eps*H3*cr);
		pts_nu->Draw("AP");
		pts_nu->SetMarkerStyle(3);
		pts_nu->SetLineColor(0);
		num->Draw("P same");
		num->SetMarkerColor(2);
		num->SetMarkerStyle(2);
		num->SetLineColor(0);
		auto leg2 = new TLegend(0.3,0.9,0.7,0.8);
		leg2->AddEntry(pts_nu.get(),"Data Point (Simulated Values)");
		leg2->AddEntry(num.get(),"Spin Tune Model (Calculated Values)");
		leg2->SetFillColorAlpha(10,0.2);
		leg2->Draw();
		cpts_nu->SetGrid();
		cpts_nu->SaveAs((path + "/error_plot.pdf)").data(),"pdf");

		TCanvas* ceh_nu = new TCanvas("ceh_nu","ceh_nu",1024,768);
		TH1F* hnu = new TH1F("hnu",
			"Distribution of residuals of spin tune fit;(#nu_{s})_{simulated} - (#nu_{s})_{model};Frequency",
			51,-5* std_res_nu,5* std_res_nu);
		for(int i=0; i<res_nu.size(); i++) hnu->Fill(res_nu[i]);
		//hnu->Fit("gaus");
		hnu->Draw();
		ceh_nu->SaveAs((path + "/error_dist.pdf)").data(),"pdf");
    }

    return make_tuple(sig0,h,H);
}