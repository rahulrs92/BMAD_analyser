#include<lv1.hh>

array<double,2> lv1::pt(string path,Data<double>& D,bool out)
{
	valarray<double> ps = D.col("T");
	valarray<double> ex = D.col("#epsilon_{x}");
	valarray<double> ey = D.col("#epsilon_{y}");
	valarray<double> cx = D.col("#xi_{x}");
	valarray<double> cy = D.col("#xi_{y}");

	auto vx = cx*ex;
	auto vy = cy*ey;
	mg["pt"] = new TList();
	plot2<double> pl1(vx,vy,ps);
	pl1->SetLineColor(0);
	valarray<double> svx(0.0,100), svy(0.0,100);
	for(int i=0; i<10; i++)
	{
		svx[slice(i*10,10,1)] = vx.min() + i * (vx.max() - vx.min())/9;
		svy[slice(i,10,10)] = vy.min() + i * (vy.max() - vy.min())/9;
	}

	matrix<3,1,valarray<double>> vA = {vx,vy,0*vy + 1};
	matrix<3,1,valarray<double>> vS = {svx,svy,0*svy + 1};
	auto A = (sum(vA * vA.tp())^-1) * sum(vA * ps);
	auto B = A.tp() * vS;
	plot2<double> pl1_fit(svx,svy,B,"Change in Phase Slip with Chromaticity",
													"#xi_{x} #epsilon_{x}","#xi_{y} #epsilon_{x}","#DeltaT/T");

	((TList*)mg["pt"])->Add(pl1_fit.get());
	((TList*)mg["pt"])->Add(pl1.get());

	auto leg = new TLegend(0.3,0.9,0.9,0.85);
	stringstream ss;
	ss<<"#frac{#DeltaT}{T} = "<<(A[0]<0 ? "- " : "")<<abs(A[0])<<"#xi_{x} #epsilon_{x}"
														<<(A[1]<0 ? " - " : " + ")<<abs(A[1])<<"#xi_{y} #epsilon_{y}"
														<<(A[2]<0 ? " - " : " + ")<<abs(A[2])<<endl;
	string ts; getline(ss,ts);
	leg->AddEntry(pl1_fit.get(),ts.data());
	leg->SetFillColorAlpha(10,0.2);
	((TList*)mg["pt"])->Add(leg);
	if(out) cout<<ts<<endl;

	if(out)
	{
		TCanvas* cpt = new TCanvas("cpt","cpt");

		((TGraph2D*)((TList*)mg["pt"])->At(0))->Draw("surf1");
		((TGraph2D*)((TList*)mg["pt"])->At(1))->Draw("P0 same");
		((TLegend*)((TList*)mg["pt"])->At(2))->Draw();

		//((TList*)mg["pt"])->Draw("surf1");
		cpt->SetGrid();
		cpt->SaveAs((path + "/pt_plot2.pdf(").data(),"pdf");
		cpt->SetTheta(15);
		cpt->SetPhi(137);
		cpt->SaveAs((path + "/pt_plot2.pdf)").data(),"pdf");
	}
	return {A[0],A[1]};
}

