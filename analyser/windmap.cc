#include<analyser.hh>

//void lv2::windmap(string path,bool out)
tuple<matrix<3,1>,matrix<3,3>> analyser::windmap(bool out)
{

	auto Dc = proton::line(path + "/rc",49);
	auto pt = lv1::pt(path + "/rc",Dc,0);

	Data<double> Dt(path + "/xi");
	valarray<double> tf = Dt.col(1), td = Dt.col(2), ts = Dt.col(3), Cx = Dt.col(4), Cy = Dt.col(5);
	fstream f;

	f.open((path + "/eta").data(),fstream::out | fstream::trunc);
	for(int i=1; i<=343; i++)
	{
		int N = 11; string path2 = path + "/r" + to_string(i);
		auto D = proton::line(path2,N,out);
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

	Data<double> Dw(path + "/eta");
	/*
	Dw.label(1) = "#xi_{x}";
	Dw.label(2) = "#xi_{y}";
	Dw.label(3) = "#epsilon_{x}";
	Dw.label(4) = "#epsilon_{y}";
	Dw.label(5) = "T";
	Dw.label(6) = "#eta_{0}";
	Dw.label(7) = "#eta_{1}";
	*/

	valarray<double> n1 = Dw.col("#eta_{1}");
	matrix<1,3> C0;//({Dw.col("#xi_{x}")[0],Dw.col("#xi_{y}")[0],n1[0]});
	C0[0] = valarray<double>(Dw.col("#xi_{x}")[tf==0.0 && td==0.0 && ts==0.0])[0];
	C0[1] = valarray<double>(Dw.col("#xi_{y}")[tf==0.0 && td==0.0 && ts==0.0])[0];
	C0[2] = valarray<double>(n1[tf==0.0 && td==0.0 && ts==0.0])[0];
	//cout<<C0<<endl;
	//n1 = valarray<double>(n1[slice(1,125,1)]);
	//cout<<"yay "<<n1.size()<<endl;
	matrix<3,1,valarray<double>> Fi({tf,td,ts});
	matrix<1,3,valarray<double>> Ci({Cx,Cy,n1});
	Ci = Ci - C0;
	matrix<3,3> A = (sum(Fi*Fi.tp())^-1) * sum(Fi*Ci);
	Xi = A.tp();
	xi0 = C0.tp();

	//cout<<"xi0 = "<<xi0<<endl;
	//cout<<"Xi = "<<Xi<<endl;
	//matrix<3,3> Xi_ns = {55.603416,3.658112,26.693864,-35.76376,-41.153072,-15.665432,-21.513543,-3.842864,-10.340451};
	//cout<<"Xi^-1 = "<<(Xi^-1)<<endl;
	//cout<<"Xi^-1 = "<<(Xi^-1)<<endl;

	f.open((path + "/matrix").data(), fstream::out | fstream::trunc);
	f<<"                                                 --->  ---->    -->"<<endl;
	f<<"This matrix can be used to perform the operation: xi  = xi0 + D  f"<<endl;
	f<<endl;
	f<<"         --->    [  Cx  ]  -->   [ TF  ]    ---->  [  Cx0  ]"<<endl;
	f<<"...where  xi  =  [  Cy  ] , f  = [ TD  ] and xi0 = [  Cy0  ]"<<endl;
	f<<"                 [ eta1 ]      = [ TSS ]           [ eta10 ]"<<endl;
	f<<endl;
	f<<"D = "<<Xi<<endl;
	f<<"xi0 = "<<xi0<<endl;
	/*
	f<<endl;
	f<<"These are the errors in D:"<<endl;
	f<<endl;
	f<<"De = "<<~De<<endl;
	*/
	f.close();

	if(out)
	{
		auto mag = Ci * Ci.tp();
		TCanvas* cpts = new TCanvas("cpts","cpts",1900,500);
		plot<double> pts(mag[mag<5000],"Vector Magnitudes of Second Order Parameters: Actual vs Linear Map",
										"Data Point Index","#xi_{x}^{2} + #xi_{y}^{2} + #eta_{1}^{2}");
		plot<double> wmp(((Xi*Fi).tp() * (Xi*Fi))[mag<5000]);
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
	}
	return make_tuple(xi0,Xi);
}


