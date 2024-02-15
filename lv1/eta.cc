#include<lv1.hh>

using namespace std;

array<double,3> lv1::eta(string path,Data<double>& D,array<double,2> iptc,bool out)
{
	D.clean_nans();
	valarray<double> dt = D.col("#delta_{0}");
	valarray<double> ps = D.col("T");

	mg["eta_1"] = new TMultiGraph("eta_1","Phase Slip Variation with #delta;#delta;#frac{#DeltaT}{T}");
	plot<double> pl1(dt,ps,"Measured Values");
	pl1->SetLineColor(0);
	((TMultiGraph*)mg["eta_1"])->Add(pl1.get(),"P");

	matrix<3,1,valarray<double>> vA = {pow(dt,0),dt,pow(dt,2)};
	auto A = (sum(vA * vA.tp())^-1) * sum(vA * ps);
	auto B = A.tp() * vA;
	//cout<<"C = "<<C<<endl;

	stringstream ss;
	ss<<"#frac{#DeltaT}{T} = "<<(A[0] < 0 ? "- " : "")<<abs(A[0])
									<<(A[1] < 0 ? " - " : " + ")<<abs(A[1])<<" #delta"
									<<(A[2] < 0 ? " - " : " + ")<<abs(A[2])<<" #delta^{2}"<<endl;
	string ts; getline(ss,ts);
	plot<double>pl1_fit(dt,B,ts);
	pl1_fit->SetLineColor(2);
	pl1_fit->SetMarkerStyle(1);
	pl1_fit->SetMarkerColor(2);
	((TMultiGraph*)mg["eta_1"])->Add(pl1_fit.get(),"L");
	((TMultiGraph*)mg["eta_1"])->GetYaxis()->SetTitleOffset(1.0);

	valarray<double> cx = D.col("#xi_{x}");
	valarray<double> cy = D.col("#xi_{y}");
	valarray<double> ex = D.col("#epsilon_{x}");
	valarray<double> ey = D.col("#epsilon_{y}");
	ps = ps - (iptc[0]*ex*cx + iptc[1]*ey*cy);
	plot<double> pl2(dt,ps,"Measured Values");
	ss.clear();
	ss<<"Phase Slip Variation with #delta (adjusted for chromaticity);";
	ss<<"#delta;#frac{#DeltaT}{T} - #frac{#DeltaT}{T}";
	ss<<"(#epsilon_{x},#epsilon_{y},#xi_{x},#xi_{y},#delta=0)";
	getline(ss,ts);
	mg["eta_2"] = new TMultiGraph("eta_2",ts.data());
	pl2->SetLineColor(0);
	((TMultiGraph*)mg["eta_2"])->Add(pl2.get(),"P");
	A = (sum(vA * vA.tp())^-1) * sum(vA * ps);
	B = A.tp() * vA;
	//cout<<"C = "<<C<<endl;

	ss.clear();
	ss<<"#frac{#DeltaT}{T} = "<<(A[0] < 0 ? "- " : "")<<abs(A[0])
									<<(A[1] < 0 ? " - " : " + ")<<abs(A[1])<<" #delta"
									<<(A[2] < 0 ? " - " : " + ")<<abs(A[2])<<" #delta^{2}"<<endl;
	getline(ss,ts);
	plot<double>pl2_fit(dt,B,ts);
	pl2_fit->SetLineColor(2);
	pl2_fit->SetMarkerStyle(1);
	pl2_fit->SetMarkerColor(2);
	((TMultiGraph*)mg["eta_2"])->Add(pl2_fit.get(),"L");
	((TMultiGraph*)mg["eta_2"])->GetYaxis()->SetTitleOffset(1.0);

	if(out)
	{
		TCanvas* ceta = new TCanvas("ceta","ceta");
		mg["eta_1"]->Draw("AP");
		ceta->BuildLegend(0.3,0.9,0.9,0.8)->SetFillColorAlpha(10,0.2);
		ceta->SetGrid();
		ceta->SaveAs((path + "/eta_plot.pdf(").data(),"pdf");

		ceta = new TCanvas("ceta2","ceta2");
		mg["eta_2"]->Draw("AP");
		ceta->BuildLegend(0.3,0.9,0.9,0.8)->SetFillColorAlpha(10,0.2);
		ceta->SetGrid();
		ceta->SaveAs((path + "/eta_plot.pdf)").data(),"pdf");
	}

	return {A[0],A[1],A[2]};
}

