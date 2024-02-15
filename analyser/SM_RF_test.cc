#include<analyser.hh>

void analyser::SM_RF_test(Data<double>& D,bool out)
{
	valarray<double> nu = D.col("#nu_{s}");
	valarray<double> cx = D.col("#xi_{x}");
	valarray<double> cy = D.col("#xi_{y}");

	plot2<double> pl1(cx,cy,nu);

	valarray<double> svx(0.0,100), svy(0.0,100);
	for(int i=0; i<10; i++)
	{
		svx[slice(i*10,10,1)] = cx.min() + i * (cx.max() - cx.min())/9;
		svy[slice(i,10,10)] = cy.min() + i * (cy.max() - cy.min())/9;
	}

	matrix<3,1,valarray<double>> vA = {cx,cy,0* cy + 1};
	matrix<3,1,valarray<double>> vS = {svx,svy,0*svy + 1};
	auto A = (sum(vA * vA.tp())^-1) * sum(vA * nu);
	auto B = A.tp() * vS;
	plot2<double> pl1_fit(svx,svy,B,"Change in Spin Tune with Chromaticity",
													"#xi_{x}","#xi_{y}","#nu_{s}");

	TCanvas* crft = new TCanvas("rftest","SM_RF_test");
	pl1_fit->Draw("surf1");
	pl1->Draw("p0 same");

	auto leg = new TLegend(0.3,0.9,0.9,0.85);
	stringstream ss;
	ss<<"#nu_{s} = "<<(A[0]<0 ? "- " : "")<<abs(A[0])<<"#xi_{x}"
									<<(A[1]<0 ? " - " : " + ")<<abs(A[1])<<"#xi_{y}"
									<<(A[2]<0 ? " - " : " + ")<<abs(A[2])<<endl;
	string ts; getline(ss,ts);
	leg->AddEntry(pl1_fit.get(),ts.data());
	leg->SetFillColorAlpha(10,0.2);
	leg->Draw();
	cout<<ts<<endl;

	crft->SetGrid();
	crft->SetTheta(30);
	crft->SetPhi(-43);
	crft->SaveAs((path + "/rft_plot.pdf(").data(),"pdf");
	crft->SetPhi(47);
	crft->SaveAs((path + "/rft_plot.pdf)").data(),"pdf");

	auto Xi_1 = Xi.sub<1,2>(2,0) * (Xi.sub<2,2>(0,0)^-1);
	//ep[0] = D.col("#epsilon_{x}").sum() / D.order()[0];
	//ep[1] = D.col("#epsilon_{y}").sum() / D.order()[0];
	//ep[2] = pow(D.col("#delta_{a}"),2).sum() / (D.order()[0]*2);
	double cx0 = abs(cx).min();
	ep[0] = valarray<double>(D.col("#epsilon_{x}")[abs(cx*cy)==abs(cx*cy).min()])[0];
	ep[1] = valarray<double>(D.col("#epsilon_{y}")[abs(cx*cy)==abs(cx*cy).min()])[0];
	ep[2] = valarray<double>(pow(D.col("#delta_{a}"),2)[abs(cx*cy)==abs(cx*cy).min()])[0] / 2;
	auto cc = ep*(M.sub<3,2>()) + (ep*matrix<3,1>(M(2)))*Xi_1;
	auto ck = ep*a + ep*matrix<3,1>(M(2))*(xi0[2] - (Xi_1*xi0.sub<2,1>()));

	cout<<"Spin tune as a function of chromaticity: "<<cc[0]<<"\t"<<cc[1]<<"\t"<<ck<<endl;
}

