#include<analyser.hh>
#include<TH3.h>
#include<TRandom3.h>

//void lv2::windmap(string path,bool out)
tuple<matrix<3,1>,matrix<3,3>> analyser::windmap(bool out)
{
	//int NS = 343; //497; //
	auto Dc = proton::line(path + "/rc", NS);
	auto pt = lv1::pt(path + "/rc",Dc,0);

	Data<double> Dt(path + "/xi");
	valarray<double> tf = Dt.col(1), td = Dt.col(2), ts = Dt.col(3), Cx = Dt.col(4), Cy = Dt.col(5);
	fstream f;

	f.open((path + "/eta").data(),fstream::in);
	if(f.is_open()) { f.close(); }
	else
	{
		f.open((path + "/eta").data(),fstream::out | fstream::trunc);
		for(int i=1; i<= NS; i++)
		{
			int N1 = 11; string path2 = path + "/r" + to_string(i);
			auto D = proton::line(path2,N1,!out);
			auto eta = lv1::eta(path2,D,pt,!out);
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

	valarray<double> n1 = Dw.col("#eta_{1}");
	matrix<1,3> C0;//({Dw.col("#xi_{x}")[0],Dw.col("#xi_{y}")[0],n1[0]});
	matrix<3,1,valarray<double>> Fi({tf,td,ts});
	matrix<1,3,valarray<double>> Ci({Cx,Cy,n1});
	if(xi0[0] == 0 && xi0[1] == 0 && xi0[2] == 0)
	{
		C0[0] = valarray<double>(Dw.col("#xi_{x}")[tf==0.0 && td==0.0 && ts==0.0])[0];
		C0[1] = valarray<double>(Dw.col("#xi_{y}")[tf==0.0 && td==0.0 && ts==0.0])[0];
		C0[2] = valarray<double>(n1[tf==0.0 && td==0.0 && ts==0.0])[0];
		Ci = Ci - C0;
		cout<<"C0 = "<<C0<<endl;
		xi0 = C0.tp();
		matrix<3,3> A = (sum(Fi*Fi.tp())^-1) * sum(Fi*Ci);
		Xi = A.tp();
	}
	else Ci = Ci - xi0.tp();


	f.open((path + "/matrix_info").data(), fstream::out | fstream::trunc);
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
	f.open((path + "/matrix").data(), fstream::out | fstream::trunc);
	f<<Xi<<endl;
	f<<xi0<<endl;
	f.close();

	valarray<double> mag_O = pow(Ci * Ci.tp(),0.5);
	valarray<double> mag_M = pow((Xi*Fi).tp() * (Xi*Fi),0.5);
	//auto PE = (Ci - (Xi*Fi).tp()) * (Ci.tp() - (Xi*Fi));
	//PE = pow(PE,0.5);

	valarray<double> PE = (Ci - (Xi*Fi).tp()) * ((Xi*Fi) / mag_M);
	PE[mag_M == 0] = 0;
	cout<<"PE = "<<PE<<endl;

	double std_PE = pow(pow(PE,2).sum()/ PE.size(),0.5);

	if(out)
	{
		TCanvas* cpts = new TCanvas("cpts","cpts",1900,500);
		plot<double> pts(mag_O[abs(PE)<10],"Vector Magnitudes of Second Order Parameters: Actual vs Linear Map",
										"Data Point Index","#xi_{x}^{2} + #xi_{y}^{2} + #eta_{1}^{2}");
		plot<double> wmp(mag_M[abs(PE)<10]);
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
		cpts->SaveAs((path + "/wm_error_plot.pdf").data(),"pdf");

		TCanvas* ceh = new TCanvas("ceh","ceh",600,400);

		TH1F* hnu = new TH1F("h1","Distribution of residuals of Linear Model fit;Error;Frequency",51,-5* std_PE,5* std_PE);
		for(int i=0; i<PE.size(); i++) hnu->Fill(PE[i]);
		//hnu->Fit("gaus");
		hnu->Draw();
		ceh->SaveAs((path + "/wm_error_dist.pdf").data(),"pdf");
	}
	return make_tuple(xi0,Xi);
}


