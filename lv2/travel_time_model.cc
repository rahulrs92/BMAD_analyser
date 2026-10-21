#include <lv2.hh>
#include<TH1F.h>


tuple<double,matrix<3,1>,matrix<3,2>> lv2::travel_time_model(Data<double>& D,string path,bool out)
{
	fstream f(path + "/tt_matrix_info", fstream::out | fstream::trunc);
	f<<"#eta_{0}"<<"\n";
    f<<"u[0]"<<"\n"<<"u[1]"<<"\n"<<"u[2]"<<"\n";
	f<<"U[0][0]"<<"\n"<<"U[0][1]"<<"\n";
	f<<"U[1][0]"<<"\n"<<"U[1][1]"<<"\n";
	f<<"U[2][0]"<<"\n"<<"U[2][1]"<<"\n";
	f.close();
	f.open(path + "/tt_matrix", fstream::out | fstream::trunc);
		
	//tie(h,H,ignore) = az.spin_model_2(ln,path);
    valarray<double> dt = D.col("#delta_{0}");
	matrix<1,3,valarray<double>> eps = {D.col("#epsilon_{x}"),D.col("#epsilon_{y}"),pow(dt,2)};
	matrix<2,1,valarray<double>> cr = {D.col("#xi_{x}"),D.col("#xi_{y}")};
	valarray<double> ps = D.col("T");

	matrix<10,1,valarray<double>> regs_ps;
    regs_ps[0] = dt;
    for(int i=0; i<3; i++) regs_ps[i+1] = eps[i];
	for(int i=0; i<6; i++) regs_ps[i+4] = eps[i/2] * cr[i%2];
	auto arr_U = (sum(regs_ps * regs_ps.tp())^-1) * sum(regs_ps * (ps));
    eta0 = arr_U[0];
    for(int i=0; i<3; i++) u[i] = arr_U[i+1];
	for(int i=0; i<6; i++) U[i/2][i%2] = arr_U[i+4];

	f<<eta0<<"\t"
    <<u[0]<<"\t"<<u[1]<<"\t"<<u[2]<<"\t"
    <<U[0][0]<<"\t"<<U[0][1]<<"\t"
	<<U[1][0]<<"\t"<<U[1][1]<<"\t"
	<<U[2][0]<<"\t"<<U[2][1]<<"\t"<<endl;
	f.close();

	cout<<"eta0 = "<<eta0<<endl;
	cout<<"u = "<<u<<endl;
	cout<<"U = "<<U<<endl;

	valarray<double> res_ps = ps - (eta0*dt + eps*u + eps*U*cr);
	double std_res_ps = pow(pow(res_ps,2).sum() / res_ps.size(),0.5);

	if(out)
	{
		TCanvas* cpts_ps = new TCanvas("cpts2","cpts2",2800,1000);
		plot<double> pts_ps(ps - eta0*dt,
								"Travel Time: Actual vs Linear Map","Data Point Index",
								"#frac{#DeltaT}{T} - #eta_{0}#delta");
		//plot<double> psm((arr_T2.tp() * regs_ps2) + eps[2]*cr[2]);//[dt == dt.max()/2]);//
		plot<double> psm(eps*u + eps*U*cr);
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
		cpts_ps->SaveAs((path + "/error_plot_2.pdf").data(),"pdf");

		TCanvas* ceh_ps = new TCanvas("ceh2","ceh2",1024,768);
		TH1F* hps = new TH1F("hps",
			"Distribution of residuals of travel-time fit;#frac{#DeltaT}{T}_{simulated} - #frac{#DeltaT}{T}_{model};Frequency",
			51,-5* std_res_ps,5* std_res_ps);
		for(int i=0; i<res_ps.size(); i++) hps->Fill(res_ps[i]);
		//hps->Fit("gaus");
		hps->Draw();
		ceh_ps->SaveAs((path + "/error_dist_2.pdf").data(),"pdf");
    }

    return make_tuple(eta0,u,U);
}