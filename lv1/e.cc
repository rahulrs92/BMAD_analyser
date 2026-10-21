#include<lv1.hh>

array<double,2> lv1::e(string path,Data<double>& iD)
{
	valarray<double> ex = iD.col("#epsilon_{x}");
	valarray<double> ey = iD.col("#epsilon_{y}");
	valarray<double> nu = iD.col("#nu_{s}");

	plot2<double> pl1(ex,ey,nu);
	pl1->SetMarkerStyle(2);

	valarray<double> sex(0.0,100), sey(0.0,100);
	for(int i=0; i<10; i++)
	{
		sex[slice(i*10,10,1)] = ex.min() + i * (ex.max() - ex.min())/9;
		sey[slice(i,10,10)] = ey.min() + i * (ey.max() - ey.min())/9;
	}

	matrix<3,1,valarray<double>> vA = {ex,ey,0*ey + 1};
	matrix<3,1,valarray<double>> vS = {sex,sey,0*sey + 1};
	auto A = (sum(vA * vA.tp())^-1) * sum(vA * nu);
	auto B = A.tp() * vS;
	plot2<double> pl1_fit(sex, sey,B,"Change in Spin Tune with Emittance @ #xi_{x},#xi_{y} = 0",
													"#epsilon_{x}","#epsilon_{y}","#nu_{s}");

	TCanvas* ce = new TCanvas("ce","ce");
	pl1_fit->Draw("surf1");
	pl1->Draw("p0 same");

	auto leg = new TLegend(0.3,0.9,0.9,0.85);
	stringstream ss;
	ss<<"#nu_{s} = "<<(A[0]<0 ? "- " : "")<<abs(A[0])<<" #epsilon_{x}"
									<<(A[1]<0 ? " - " : " + ")<<abs(A[1])<<" #epsilon_{y}"
									<<(A[2]<0 ? " - " : " + ")<<abs(A[2])<<endl;
	string ts; getline(ss,ts);
	leg->AddEntry(pl1_fit.get(),ts.data());
	leg->SetFillColorAlpha(10,0.2);
	leg->Draw();
	cout<<ts<<endl;

	ce->SetGrid();
	ce->SaveAs((path + "/e_plot2.pdf(").data(),"pdf");
	ce->SetTheta(15);
	ce->SetPhi(137);
	ce->SaveAs((path + "/e_plot2.pdf)").data(),"pdf");

	return {A[0],A[1]};
}


