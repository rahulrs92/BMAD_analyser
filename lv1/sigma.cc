#include<lv1.hh>

using namespace std;

array<double,3> lv1::sigma(string path,Data<double>& D,array<double,2> ie,array<double,2> ic)
{
	valarray<double> dt = D.col("#delta_{0}");
	valarray<double> nu = D.col("#nu_{s}");

	TCanvas* cnuv = new TCanvas("cnuv","cnuv");
	plot<double> pl1(dt,nu,"Spin Tune Variation with #delta","#delta","#nu_{s}");
	pl1->Draw("AP");

	matrix<3,1,valarray<double>> vA = {pow(dt,0),dt,pow(dt,2)};
	auto A = (sum(vA * vA.tp())^-1) * sum(vA * nu);
	auto B = A.tp() * vA;
	//cout<<"C = "<<C<<endl;

	plot<double>pl1_fit(dt,B);
	pl1_fit->SetLineColor(2);
	pl1_fit->Draw("L same");

	auto leg = new TLegend(0.3,0.9,0.9,0.85);
	stringstream ss;
	ss<<"#nu_{s} = "<<(A[0] < 0 ? "- " : "")<<abs(A[0])
									<<(A[1] < 0 ? " - " : " + ")<<abs(A[1])<<" #delta"
									<<(A[2] < 0 ? " - " : " + ")<<abs(A[2])<<" #delta^{2}"<<endl;
	string ts; getline(ss,ts);
	leg->AddEntry(pl1_fit.get(),ts.data());
	leg->SetFillColorAlpha(10,0.2);
	leg->Draw();
	cnuv->SetGrid();
	cnuv->SaveAs((path + "/sigma_plot.pdf(").data(),"pdf");

	valarray<double> ex = D.col("#epsilon_{x}");
	valarray<double> ey = D.col("#epsilon_{y}");
	valarray<double> cx = D.col("#xi_{x}");
	valarray<double> cy = D.col("#xi_{y}");
	nu = nu - (ie[0]*ex + ie[1]*ey) - (ic[0]*ex*cx + ic[1]*ey*cy);
	A = (sum(vA * vA.tp())^-1) * sum(vA * nu);
	B = A.tp() * vA;
	//cout<<"C = "<<C<<endl;

	plot<double> pl2(dt,nu,"Spin Tune Factor Variation with #delta (adjusted for emittance and chromaticity)",
									 "#delta","#nu_{s} - #nu_{s}(#epsilon_{x},#epsilon_{y},#xi_{x},#xi_{y},#delta=0)");
	pl2->Draw("AP");
	plot<double> pl2_fit(dt,B);
	pl2_fit->SetLineColor(2);
	pl2_fit->Draw("L same");

	auto leg2 = new TLegend(0.3,0.9,0.9,0.85);
	ss.clear();
	ss<<"#nu_{s} = "<<(A[0] < 0 ? "- " : "")<<abs(A[0])
									<<(A[1] < 0 ? " - " : " + ")<<abs(A[1])<<" #delta"
									<<(A[2] < 0 ? " - " : " + ")<<abs(A[2])<<" #delta^{2}"<<endl;
	getline(ss,ts);
	leg2->AddEntry(pl1_fit.get(),ts.data());
	leg2->SetFillColorAlpha(10,0.2);
	leg2->Draw();
	cnuv->SetGrid();
	cnuv->SaveAs((path + "/sigma_plot.pdf)").data(),"pdf");

	return {A[0],A[1],A[2]};
}


