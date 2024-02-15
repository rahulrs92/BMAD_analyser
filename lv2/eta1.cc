#include<lv2.hh>

array<double,3> lv2::eta1(string path,Data<double>& D)
{
	if(f_eta1) { return v_eta1; } else { f_eta1 = true; }
	valarray<double> cx = D.col("#xi_{x}");
	valarray<double> cy = D.col("#xi_{y}");
	valarray<double> eta1 = D.col("#eta_{1}");

	valarray<double> scx(0.0,100), scy(0.0,100);
	for(int i=0; i<10; i++)
	{
		scx[slice(i*10,10,1)] = cx.min() + i * (cx.max() - cx.min())/9;
		scy[slice(i,10,10)] = cy.min() + i * (cy.max() - cy.min())/9;
	}

	plot2<double> pl1(cx,cy,eta1);

	matrix<3,1,valarray<double>> vA = {cx,cy,0*cx + 1}, vS = {scx, scy,0*scx + 1};
	auto A = (sum(vA * vA.tp())^-1) * sum(vA * eta1);
	auto B = A.tp() * vS;
	plot2<double> pl1_fit(scx, scy, B,"Variation of Second Order Phase Slip Factor with Chromaticity",
								"#xi_{x}","#xi_{y}","#eta_{1}");

	TCanvas* c1 = new TCanvas("ceta1","ceta1");
	pl1_fit->Draw("surf1");
	pl1->Draw("p0 same");

	auto leg = new TLegend(0.3,0.9,0.9,0.85);
	stringstream ss; string ts;
	ss<<"#eta_{1} = "<<(A[0]<0 ? "- " : "")<<abs(A[0])<<" #xi_{x}"
														<<(A[1]<0 ? " - " : " + ")<<abs(A[1])<<" #xi_{y}"
														<<(A[2]<0 ? " - " : " + ")<<abs(A[2])<<endl;
	getline(ss,ts);
	leg->AddEntry(pl1_fit.get(),ts.data());
	leg->SetFillColorAlpha(10,0.2);
	leg->Draw();
	cout<<ts<<endl;

	c1->SetGrid();
	c1->SaveAs((path + "/eta1_plot2.pdf(").data(),"pdf");
	c1->SetTheta(15);
	c1->SetPhi(137);
	c1->SaveAs((path + "/eta1_plot2.pdf)").data(),"pdf");

	v_eta1 = {A[0],A[1],A[2]};
	return v_eta1;
}

