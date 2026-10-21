#include<lv1.hh>

array<double,2> lv1::c(string path,Data<double>& D,array<double,2> ie)
{
	valarray<double> nu = D.col("#nu_{s}");
	valarray<double> ex = D.col("#epsilon_{x}");
	valarray<double> ey = D.col("#epsilon_{y}");
	valarray<double> cx = D.col("#xi_{x}");
	valarray<double> cy = D.col("#xi_{y}");

	auto vx = cx*ex;
	auto vy = cy*ey;
	plot2<double> pl1(vx,vy,nu);

	valarray<double> svx(0.0,100), svy(0.0,100);
	for(int i=0; i<10; i++)
	{
		svx[slice(i*10,10,1)] = vx.min() + i * (vx.max() - vx.min())/9;
		svy[slice(i,10,10)] = vy.min() + i * (vy.max() - vy.min())/9;
	}

	matrix<3,1,valarray<double>> vA = {vx,vy,0*vy + 1};
	matrix<3,1,valarray<double>> vS = {svx,svy,0*svy + 1};
	auto A = (sum(vA * vA.tp())^-1) * sum(vA * nu);
	auto B = A.tp() * vS;
	plot2<double> pl1_fit(svx,svy,B,"Change in Spin Tune with Chromaticity",
													"#xi_{x} #epsilon_{x}","#xi_{y} #epsilon_{x}","#nu_{s}");

	TCanvas* cc = new TCanvas("cc","cc");
	pl1_fit->Draw("surf1");
	pl1->Draw("p0 same");

	auto leg = new TLegend(0.3,0.9,0.9,0.85);
	stringstream ss;
	ss<<"#nu_{s} = "<<(A[0]<0 ? "- " : "")<<abs(A[0])<<"#xi_{x} #epsilon_{x}"
									<<(A[1]<0 ? " - " : " + ")<<abs(A[1])<<"#xi_{y} #epsilon_{y}"
									<<(A[2]<0 ? " - " : " + ")<<abs(A[2])<<endl;
	string ts; getline(ss,ts);
	leg->AddEntry(pl1_fit.get(),ts.data());
	leg->SetFillColorAlpha(10,0.2);
	leg->Draw();
	cout<<ts<<endl;

	cc->SetGrid();
	cc->SetTheta(30);
	cc->SetPhi(-43);
	cc->SaveAs((path + "/c_plot2.pdf(").data(),"pdf");
	cc->SetPhi(47);
	cc->SaveAs((path + "/c_plot2.pdf").data(),"pdf");

	cc = new TCanvas("cc2","cc2");
	nu = nu - (ie[0]*ex + ie[1]*ey);
	plot2<double> pl2(vx,vy,nu);
	A = (sum(vA * vA.tp())^-1) * sum(vA * nu);
	B = A.tp() * vS;
	plot2<double> pl2_fit(svx,svy,B,"Change in Spin Tune with Chromaticity (adjusted for emittance)",
													"#xi_{x} #epsilon_{x}","#xi_{y} #epsilon_{y}",
											 "#nu_{s} - #nu_{s}(#epsilon_{x},#epsilon_{y},#xi_{x}=0,#xi_{y}=0)");

	pl2_fit->Draw("surf1");
	pl2->Draw("p0 same");

	auto leg2 = new TLegend(0.3,0.9,0.9,0.85);
	ss.clear();
	ss<<"#nu_{s} = "<<(A[0]<0 ? "- " : "")<<abs(A[0])<<"#xi_{x} #epsilon_{x}"
									<<(A[1]<0 ? " - " : " + ")<<abs(A[1])<<"#xi_{y} #epsilon_{y}"
									<<(A[2]<0 ? " - " : " + ")<<abs(A[2])<<endl;
	getline(ss,ts);
	leg2->AddEntry(pl2_fit.get(),ts.data());
	leg2->SetFillColorAlpha(10,0.2);
	leg2->Draw();
	cout<<ts<<endl;

	cc->SetGrid();
	cc->SetTheta(30);
	cc->SetPhi(-43);
	cc->SaveAs((path + "/c_plot2.pdf").data(),"pdf");
	cc->SetPhi(47);
	cc->SaveAs((path + "/c_plot2.pdf)").data(),"pdf");

	return {A[0],A[1]};
}
