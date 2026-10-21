#include<analyser.hh>

tuple<matrix<3,1>,matrix<3,2>,matrix<2,1>> lv2::spin_model_2R(Data<double>& D,string path,bool out)
{
	valarray<double> dt = D.col("#delta_{a}");
	valarray<double> nu = D.col("#nu_{s}");
	valarray<double> ps = D.col("T");

	matrix<1,3,valarray<double>> eps =
						{D.col("#epsilon_{x}"),D.col("#epsilon_{y}"),pow(dt,2)};
	matrix<3,1,valarray<double>> tx =
						{D.col("t_{f}"),D.col("t_{d}"),D.col("t_{s}")};
	matrix<2,1,valarray<double>> cr =
						{D.col("#xi_{x}"),D.col("#xi_{y}")};

	//cout<<"b = "<<b<<endl;
	//cout<<"M = "<<M<<endl;

	matrix<9,1,valarray<double>> regs_nu;
	for(int i=0; i<3; i++) regs_nu[i] = eps[i];
	for(int i=0; i<6; i++) regs_nu[i+3] = eps[i/2] * cr[i%2];
	auto arr_V = (sum(regs_nu * regs_nu.tp())^-1) * sum(regs_nu * nu);
	for(int i=0; i<3; i++) w[i] = arr_V[i];
	for(int i=0; i<6; i++) W[i/2][i%2] = arr_V[i+3];

	cout<<"b = "<<w<<endl;
	cout<<"M = "<<W<<endl;

	valarray<double> res_nu = nu - (eps*w + eps*W*cr);
	double std_res_nu = pow(pow(res_nu,2).sum() / res_nu.size(),0.5);

	string path1 = path.substr(0,path.find_last_of("/"));
	fstream f(path1 + "/model_info",fstream::out | fstream::trunc);
	f<<"#nu_{s}\t#nu_{s}_model"<<endl;
	f.close();
	f.open(path1 + "/model",fstream::out | fstream::app);
	valarray<double> model = (eps*w + eps*W*cr);
	for(int i=0; i<nu.size(); i++) f<<nu[i]<<"\t"<<model[i]<<endl;
	f.close();

	if(out)
	{
		TCanvas* cpts_nu = new TCanvas("cpts","cpts",2800,1000);

		plot<double> pts_nu(nu,//[dt == dt.max()/2],//
										 "Spin tunes: Actual vs Model","Data Point Index",
									 "#nu_{s}");
		plot<double> stm(eps*w + eps*W*cr);//arr_V.tp()*regs_nu);//(sig0*dt + (sig0*dt);//
		pts_nu->Draw("AP");
		pts_nu->SetMarkerStyle(3);
		pts_nu->SetLineColor(0);
		stm->Draw("P same");
		stm->SetMarkerColor(2);
		stm->SetMarkerStyle(2);
		stm->SetLineColor(0);
		auto leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts_nu.get(),"Data Point (Measured Values)");
		leg->AddEntry(stm.get(),"Spin Tune Model (Calculated Values)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cpts_nu->SetGrid();
		cpts_nu->SaveAs((path + "/error_plot_2.pdf").data(),"pdf");

		TCanvas* ceh_nu = new TCanvas("ceh","ceh",1024,768);
		TH1F* hnu = new TH1F("hnu",
			"Distribution of residuals of spin-tune fit;#nu_{simulated} - #nu_{model};Frequency",
			51,-5*std_res_nu,5*std_res_nu);
		for(int i=0; i<res_nu.size(); i++) hnu->Fill(res_nu[i]);
		//hnu->Fit("gaus");
		hnu->Draw();
		ceh_nu->SaveAs((path + "/error_dist_2.pdf").data(),"pdf");
	}
	//cout<<"unique elements of tx[2]: "<<unique(D.col("t_{s}"))<<endl;
	xi = ((W.tp() * W)^-1) * W.tp() * (-w);
	cout<<"xi = "<<xi<<endl;
	return make_tuple(w,W,xi);
}





