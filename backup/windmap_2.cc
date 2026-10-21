#include<analyser.hh>

//void lv2::windmap(string path,bool out)
tuple<matrix<3,1>,matrix<3,3>> lv2::windmap(string wpath,bool out)
{
	int NS = 125; //343; //497; //
	//auto Dc = proton::line(wpath + "/rc", NS);
	//auto pt = lv1::pt(wpath + "/rc",Dc,0);

	Data<double> Dt(wpath + "/xi");
	valarray<double> Cx = Dt.col(4), Cy = Dt.col(5);
	matrix<3,1,valarray<double>> tx = {Dt.col(1),Dt.col(2),Dt.col(3)};
	fstream f;

	f.open((wpath + "/eta").data(),fstream::in);
	if(f.is_open()) { f.close(); }
	else
	{
		f.open((wpath + "/eta").data(),fstream::out | fstream::trunc);
		for(int i=1; i<= NS; i++)
		{
			int N1 = 21; string path2 = wpath + "/r" + to_string(i);
			auto D = proton::line(path2,N1);
			valarray<double> dt = D.col("#delta_{0}");
			valarray<double> ps = D.col("T");
			valarray<double> epx = D.col("#epsilon_{x}");
			valarray<double> epy = D.col("#epsilon_{y}");
			matrix<3,1,valarray<double>> regs_ps = {pow(dt,0),dt,pow(dt,2)};//};//
			auto arr_ps = (sum(regs_ps * regs_ps.tp())^-1) * sum(regs_ps * ps);
			auto arr_epx = (sum(regs_ps * regs_ps.tp())^-1) * sum(regs_ps * epx);
			auto arr_epy = (sum(regs_ps * regs_ps.tp())^-1) * sum(regs_ps * epy);

			valarray<double> res_ps = ps - arr_ps.tp()*regs_ps;
			double std_res_ps = pow(pow(res_ps,2).sum()/ res_ps.size(),0.5);

			if(out)
			{
				stringstream ss;
				string ts;

				TCanvas* cpspts = new TCanvas("cpspts","cpspts",1000,500);
				plot<double> pspts(dt,ps,"Travel Time: Variation with Momentum Offset",
												"#delta","#frac{#DeltaT}{T}");
				plot<double> psm(dt,arr_ps.tp()*regs_ps);
				pspts->Draw("AP");
				pspts->SetMarkerStyle(3);
				pspts->SetLineColor(0);
				psm->Draw("L same");
				psm->SetMarkerColor(2);
				psm->SetMarkerStyle(0);
				psm->SetLineColor(2);
				auto leg = new TLegend(0.3,0.9,0.7,0.8);
				leg->AddEntry(pspts.get(),"Data Point (Simulated Values)");
				ss<<"#frac{#DeltaT}{T} = "<<(arr_ps[0]<0 ? "- " : "")<<abs(arr_ps[0])
																	<<(arr_ps[1]<0 ? " - " : " + ")<<abs(arr_ps[1])<<" #delta"
																	<<(arr_ps[2]<0 ? " - " : " + ")<<abs(arr_ps[2])<<" #delta^{2}"<<endl;
				getline(ss,ts);
				leg->AddEntry(psm.get(),ts.data());
				leg->SetFillColorAlpha(10,0.2);
				leg->Draw();
				cpspts->SetGrid();
				cpspts->SaveAs((path2 + "/ps_error_plot.pdf(").data(),"pdf");

				TCanvas* cepxpts = new TCanvas("cepxpts","cepxpts",1000,500);
				plot<double> epxpts(dt,epx,"Horizontal Emittance: Variation with Momentum Offset",
												"#delta","#epsilon_{x}");
				plot<double> epxm(dt,arr_epx.tp()*regs_ps);
				epxpts->Draw("AP");
				epxpts->SetMarkerStyle(3);
				epxpts->SetLineColor(0);
				epxm->Draw("L same");
				epxm->SetMarkerColor(2);
				epxm->SetMarkerStyle(0);
				epxm->SetLineColor(2);
				auto leg2 = new TLegend(0.3,0.9,0.7,0.8);
				leg2->AddEntry(epxpts.get(),"Data Point (Simulated Values)");
				ss.clear();
				ss<<"#epsilon_{x} = "<<(arr_epx[0]<0 ? "- " : "")<<abs(arr_epx[0])
															<<(arr_epx[1]<0 ? " - " : " + ")<<abs(arr_epx[1])<<" #delta"
															<<(arr_epx[2]<0 ? " - " : " + ")<<abs(arr_epx[2])<<" #delta^{2}"<<endl;
				getline(ss,ts);
				leg2->AddEntry(epxm.get(),ts.data());
				leg2->SetFillColorAlpha(10,0.2);
				leg2->Draw();
				cepxpts->SetGrid();
				cepxpts->SaveAs((path2 + "/ps_error_plot.pdf").data(),"pdf");

				TCanvas* cepypts = new TCanvas("cepypts","cepypts",1000,500);
				plot<double> epypts(dt,epy,"Vertical Emittance: Variation with Momentum Offset",
												"#delta","#epsilon_{y}");
				plot<double> epym(dt,arr_epy.tp()*regs_ps);
				epypts->Draw("AP");
				epypts->SetMarkerStyle(3);
				epypts->SetLineColor(0);
				epym->Draw("L same");
				epym->SetMarkerColor(2);
				epym->SetMarkerStyle(0);
				epym->SetLineColor(2);
				auto leg3 = new TLegend(0.3,0.9,0.7,0.8);
				leg3->AddEntry(epypts.get(),"Data Point (Simulated Values)");
				ss<<"#epsilon_{y} = "<<(arr_epy[0]<0 ? "- " : "")<<abs(arr_epy[0])
															<<(arr_epy[1]<0 ? " - " : " + ")<<abs(arr_epy[1])<<" #delta"
															<<(arr_epy[2]<0 ? " - " : " + ")<<abs(arr_epy[2])<<" #delta^{2}"<<endl;
				getline(ss,ts);
				leg3->AddEntry(epym.get(),ts.data());
				leg3->SetFillColorAlpha(10,0.2);
				leg3->Draw();
				cepypts->SetGrid();
				cepypts->SaveAs((path2 + "/ps_error_plot.pdf)").data(),"pdf");

				TCanvas* cep0 = new TCanvas("cep0","cep0",1900,1000);
				plot<double> pep0(dt,D.col("#epsilon_{x0}"),"Correlation",
												"#delta","ellipse offset");
				//lat_plot.axis_equal();
				pep0->Draw("AP*");
				cep0->SetGrid();
				cep0->SaveAs((path2 + "/ep0_correlation.pdf").data(),"pdf");
				/*
				TCanvas* cpseh = new TCanvas("cpseh","cpseh",1024,768);
				TH1F* hnu = new TH1F("h1","Distribution of residuals of Linear Model fit;Error;Frequency",51,-5* std_res_ps,5* std_res_ps);
				for(int i=0; i<res_ps.size(); i++) hnu->Fill(res_ps[i]);
				//hnu->Fit("gaus");
				hnu->Draw();
				cpseh->SaveAs((path2 + "/ps_error_dist.pdf").data(),"pdf");
				*/
			}

			f<<D.col("#xi_{x}")[0]<<"\t"<<D.col("#xi_{y}")[0]<<"\t"
				<<arr_ps[0]<<"\t"<<arr_ps[1]<<"\t"<<arr_ps[2]<<"\t"//endl;//
				<<arr_epx[0]<<"\t"<<arr_epx[1]<<"\t"<<arr_epx[2]<<"\t"//endl;//
				<<arr_epy[0]<<"\t"<<arr_epy[1]<<"\t"<<arr_epy[2]<<endl;//endl;//
		}
		f.close();
		f.open((wpath + "/eta_info").data(),fstream::out | fstream::trunc);
		f<<"#xi_{x}\n#xi_{y}\nT\n#eta_{0}\n#eta_{1}\n"
			<<"#epsilon_{x0}\n#epsilon_{x1}\n#epsilon_{x2}\n"//"<<endl;//
			<<"#epsilon_{y0}\n#epsilon_{y1}\n#epsilon_{y2}"<<endl;//"<<endl;//
		f.close();
	}
	Data<double> Dw(wpath + "/eta");

	matrix<1,3,valarray<double>> cr({Cx,Cy,Dw.col("#eta_{1}")});
	matrix<1,2,valarray<double>> eps({Dw.col("#epsilon_{x0}"),Dw.col("#epsilon_{y0}")});

	matrix<4,1,valarray<double>> regs_wm = {pow(tx[0],0),tx[0],tx[1],tx[2]};
	auto X1 = ((sum(regs_wm * regs_wm.tp())^-1) * sum(regs_wm * cr)).tp();
	cout<<"YAY!!"<<endl;
	xi0 = X1.sub<3,1>();
	Xi = X1.sub<3,3>(0,1);
	cout<<"xi0 = "<<xi0<<endl;
	cout<<"Xi = "<<Xi<<endl;
	cout<<"t0 = "<<(Xi^-1) * (-xi0)<<endl;

	//cout<<"measured eta1 at first point: "<<cr[2]<<endl;
	//cout<<"calculated eta1 at first point: "<<(X1[2] * regs_wm)[1]<<endl;

	f.open((wpath + "/matrix_info").data(), fstream::out | fstream::trunc);
	f<<"                                                 --->  ---->    -->"<<endl;
	f<<"This matrix can be used to perform the operation: xi  = xi0 + D  f"<<endl;
	f<<endl;
	f<<"         --->    [  Cx  ]  -->   [ TF  ]    ---->  [  Cx0  ]"<<endl;
	f<<"...where  xi  =  [  Cy  ] , f  = [ TD  ] and xi0 = [  Cy0  ]"<<endl;
	f<<"                 [ eta1 ]      = [ TSS ]           [ eta10 ]"<<endl;
	f<<endl;
	f<<"Lines 1-3:"<<endl;
	f<<"D = "<<Xi<<endl;
	f<<"Lines 5-7:"<<endl;
	f<<"xi0 = "<<xi0<<endl;
	/*
	f<<endl;
	f<<"These are the errors in D:"<<endl;
	f<<endl;
	f<<"De = "<<~De<<endl;
	*/
	f.close();
	f.open((wpath + "/matrix").data(), fstream::out | fstream::trunc);
	f<<Xi<<endl;
	f<<xi0<<endl;
	f.close();

	matrix<4,1,valarray<double>> regs_T;
	//regs_T[0] = pow(eps[0],0);
	//for(int i=0; i<2; i++) { regs_T[i+1] = eps[i]; }
	for(int i=0; i<4; i++) { regs_T[i] = eps[i/2] * cr[i%2]; }
	auto arr_T = (sum(regs_T * regs_T.tp())^-1) * sum(regs_T * Dw.col("T"));
	//double T00 = arr_T[0];
	//matrix<2,1> T0;
	//for(int i=0; i<2; i++) T0[i] = arr_T[i+1];
	matrix<2,2> T;
	for(int i=0; i<4; i++) T[i/2][i%2] = arr_T[i];
	//cout<<"T00 = "<<T00<<endl;
	//cout<<"T0 = "<<T0<<endl;
	cout<<"T = "<<T<<endl;



	valarray<double> mag_O = cr * cr.tp();
	valarray<double> mag_M = (xi0 + Xi*tx).tp() * (xi0 + Xi*tx);

	valarray<double> res_mag = mag_O - mag_M;
	double std_res_mag = pow(pow(res_mag,2).sum()/ res_mag.size(),0.5);

	valarray<double> res_tc = Dw.col("T") - (arr_T.tp() * regs_T);
	double std_res_tc = pow(pow(res_tc,2).sum()/ res_tc.size(),0.5);

	if(out)
	{
		TCanvas* cpts = new TCanvas("cpts","cpts",1900,500);
		plot<double> pts(mag_O,"Vector Magnitudes of Second Order Parameters: Actual vs Linear Map",
										"Data Point Index","#xi_{x}^{2} + #xi_{y}^{2} + #eta_{1}^{2}");
		plot<double> wmp(mag_M);
		pts->Draw("AP");
		pts->SetMarkerStyle(3);
		pts->SetLineColor(0);
		wmp->Draw("P same");
		wmp->SetMarkerColor(2);
		wmp->SetMarkerStyle(2);
		wmp->SetLineColor(0);
		auto leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts.get(),"Data Point (Measured Values)");
		leg->AddEntry(wmp.get(),"Linear Map (Calculated Values)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cpts->SetGrid();
		cpts->SaveAs((wpath + "/wm_error_plot.pdf").data(),"pdf");

		TCanvas* cRpts = new TCanvas("cRpts","cRpts",1900,500);
		plot<double> Rpts(Dw.col("T"),
										"Transverse contribution to travel time: Actual vs T Model",
										"Data Point Index","#frac{#DeltaT}{T}_{transverse}");
		plot<double> Rmp(arr_T.tp()*regs_T);
		Rpts->Draw("AP");
		Rpts->SetMarkerStyle(3);
		Rpts->SetLineColor(0);
		Rmp->Draw("P same");
		Rmp->SetMarkerColor(2);
		Rmp->SetMarkerStyle(2);
		Rmp->SetLineColor(0);
		leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts.get(),"Data Point (Measured Values)");
		leg->AddEntry(wmp.get(),"T Model (Calculated Values)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cRpts->SetGrid();
		cRpts->SaveAs((wpath + "/TC_plot.pdf").data(),"pdf");

		TCanvas* ceh = new TCanvas("ceh","ceh",1024,768);

		TH1F* hnu = new TH1F("hnu","Distribution of residuals of Linear Model fit;Error;Frequency",
												 51,-5* std_res_mag,5* std_res_mag);
		for(int i=0; i<res_mag.size(); i++) hnu->Fill(res_mag[i]);
		//hnu->Fit("gaus");
		hnu->Draw();
		ceh->SaveAs((wpath + "/wm_error_dist.pdf").data(),"pdf");

		TCanvas* cteh = new TCanvas("cteh","cteh",1024,768);

		TH1F* htc = new TH1F("htc","Distribution of residuals of Transverse Travel-Time fit;Error;Frequency",
												 51,-5* std_res_tc,5* std_res_tc);
		for(int i=0; i<res_tc.size(); i++) htc->Fill(res_tc[i]);
		//hnu->Fit("gaus");
		htc->Draw();
		cteh->SaveAs((wpath + "/TC_error_dist.pdf").data(),"pdf");


		new TCanvas("clatplot","clatplot",1900,1000);
		plot<double> lat_plot(cr[2],Dw.col("#epsilon_{x1}"),"Correlation",
										"#eta_{1}","#epsilon_{x2}");
		//lat_plot.axis_equal();
		lat_plot->Draw("AP0");
	}
	return make_tuple(xi0,Xi);
}

tuple<matrix<3,1>,matrix<3,3>> lv2::windmap(bool out) { return lv2::windmap(path,out); }


