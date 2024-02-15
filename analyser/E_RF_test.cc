#include<analyser.hh>

void analyser::E_RF_test(Data<double>& De,bool out)
{
	valarray<double> nu = De.col("#nu_{s}");
	matrix<1,3,valarray<double>> ep = {De.col("#epsilon_{x}"), De.col("#epsilon_{y}"),pow(De.col("#delta_{a}"),2)/2};
	matrix<1,3,valarray<double>> ep1 = {De.col("#epsilon_{x}"), De.col("#epsilon_{y}"),pow(De.col("#delta_{a}"),2)/2};
	for(int i=0; i<ep[0].size(); i++) ep1[0][i] = ep[0][i/11];

	Data<double> Dp(path + "/rtep/pos");
	valarray<double> x = Dp.col(1), y = Dp.col(2), pz = (Dp.col(3)/2);
	matrix<3,1,valarray<double>> pos_2 = {pow(ep1[0],1),pow(ep1[0]*ep1[2],0.5),pow(ep1[2],1)};//,
																	//pow(x,6),pow(pz,6),pow(pow(x,2)*pz,2),pow(x*pow(pz,2),2)};
	matrix<2,1,valarray<double>> ep2 = {pow(ep[0],1),pow(ep[2],1)};
	auto a1 = (sum(pos_2 * pos_2.tp())^-1) * sum(pos_2 * ep2.tp());

	auto s_ep = a1.tp() * pos_2;

	cout<<a1<<endl;
	if(out)
	{
		TCanvas* cpts = new TCanvas("cpts","cpts",1900,500);

		plot<double> pts(ep[2],ep[0],"Horizontal emittance variation with #frac{#delta_{a}}{2}",
										 "#frac{#delta_{a}}{2}","#epsilon_{x}");
		plot<double> wmp(s_ep[1],s_ep[0]);
		pts->Draw("AP");
		pts->SetMarkerStyle(3);
		pts->SetLineColor(0);
		wmp->Draw("P same");
		wmp->SetMarkerColor(2);
		wmp->SetMarkerStyle(2);
		wmp->SetLineColor(0);
		auto leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts.get(),"Data Point (Measured Values)");
		leg->AddEntry(wmp.get(),"Linear Map (Calculated Values) (RF off)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cpts->SetGrid();
		cpts->SaveAs((path + "/emittance_plot.pdf(").data(),"pdf");

		TCanvas* cmg = new TCanvas("cmg","cmg",1900,500);
		TMultiGraph* mg = new TMultiGraph("mg",
						"Horizontal emittance vs x variation with #frac{#delta_{a}}{2};x;#epsilon_{x}");
		//plot<double> wmp(s_ep[1][slice(110,11,1)],s_ep[0][slice(110,11,1)]);
		for(int i=0; i<11; i++)
		{
			plot<double> mpts(pow(x,1)[slice(i,11,11)],pow(ep[0],1)[slice(i,11,11)],
											 "#frac{#delta_{a}}{2} = " + to_string(ep[2][i]));
			mpts->SetMarkerStyle(3);
			mpts->SetLineColor(2);
			mpts->Fit("pol2");
			mg->Add(mpts.get());
		}
		mg->Draw("AP");
		//cpts->BuildLegend(0.3,0.9,0.7,0.8);
		cmg->SetGrid();
		cmg->SaveAs((path + "/emittance_plot.pdf)").data(),"pdf");
	}
}



