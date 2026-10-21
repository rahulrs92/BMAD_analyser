#include<analyser.hh>
#include<TH1F.h>
#include<TH2F.h>

void lv2::spin_model_errors(bool out)
{
	int N1 = 11;
	fstream f;
	f.open((path + "/eta").data(),fstream::in);
	if(f.is_open()) { f.close(); }
	else
	{
		auto Dc = proton::line(path + "/rc", NS);
		auto pt = lv1::pt(path + "/rc",Dc);
		f.open((path + "/eta").data(),fstream::out | fstream::trunc);
		for(int i=1; i<= NS; i++)
		{
			string path2 = path + "/r" + to_string(i);
			auto D = proton::line(path2,N1,out);
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
	}
	Data<double> Dw(path + "/eta");

	vector<double> v_ex, v_ey, v_dt, v_cx, v_cy, v_eta1, v_nu, v_ps;
	auto tr_eta1 = Dw.col("#eta_{1}");
	for(int i=0; i< NS; i++)
	{
		string path2 = path + "/r" + to_string(i+1);
		auto D = proton::line(path2,N1,out);
		auto tr_ex = D.col("#epsilon_{x}");
		auto tr_ey = D.col("#epsilon_{y}");
		auto tr_dt = D.col("#delta_{0}");
		auto tr_cx = D.col("#xi_{x}");
		auto tr_cy = D.col("#xi_{y}");
		auto tr_nu = D.col("#nu_{s}");
		auto tr_ps = D.col("T");
		for(int j=0; j<N1; j++)
		{
			v_ex.push_back(tr_ex[j]);
			v_ey.push_back(tr_ey[j]);
			v_dt.push_back(tr_dt[j]);
			v_cx.push_back(tr_cx[j]);
			v_cy.push_back(tr_cy[j]);
			v_nu.push_back(tr_nu[j]);
			v_ps.push_back(tr_ps[j]);
			v_eta1.push_back(tr_eta1[i]);
			if(isnan(tr_ps[j]) || isnan(tr_nu[j]))
			{
				v_ex.pop_back();
				v_ey.pop_back();
				v_dt.pop_back();
				v_cx.pop_back();
				v_cy.pop_back();
				v_nu.pop_back();
				v_ps.pop_back();
				v_eta1.pop_back();
			}
		}
	}

	matrix<3,1,valarray<double>> cr;
	matrix<1,3,valarray<double>> eps;

	eps[0] = valarray<double>(v_ex.data(),v_ex.size());
	eps[1] = valarray<double>(v_ey.data(),v_ey.size());
	eps[2] = pow(valarray<double>(v_dt.data(),v_dt.size()),2);
	cr[0] = valarray<double>(v_cx.data(),v_cx.size());
	cr[1] = valarray<double>(v_cy.data(),v_cx.size());
	cr[2] = valarray<double>(v_eta1.data(),v_eta1.size());
	valarray<double> nu(v_nu.data(),v_nu.size());
	valarray<double> ps(v_ps.data(),v_ps.size());
	valarray<double> dt(v_dt.data(),v_dt.size());

	cr = cr - xi0;

	valarray<double> res_nu = abs(nu - (sig0*dt + eps*a + eps*V*cr));
	double std_res_nu = pow(pow(res_nu,2).sum()/ res_nu.size(),0.5);

	valarray<double> th = TMath::Pi()/2 - acos(cr[2]/pow(cr.tp()*cr,0.5)), phi = atan2(cr[1],cr[0]);
	cout<<"theta: "<<th<<endl;
	cout<<"point radii: "<<valarray<double>(pow(cr.tp() * cr,0.5))<<endl;
	//cout<<"yay!"<<endl;
	/*
	th[cr[2]>0] = atan(pow(pow(cr[0],2) + pow(cr[1],2),0.5) / cr[2])[cr[2]>0];
	th[cr[2]<0] = atan(TMath::Pi() + pow(pow(cr[0],2) + pow(cr[1],2),0.5) / cr[2])[cr[2]<0];
	th[cr[2]==0] = (TMath::Pi() / 2);
	*/

	phi[cr[1]==0 && cr[0]==0] = 0;
	/*
	ch2 /= pow(nu,2).sum() / (nu.size() - 1);
	ch2 /= nu.size() - 13;
	*/
	if(out)
	{
		TCanvas* cpts = new TCanvas("cpts","cpts",1900,500);

		plot<double> pts(nu[dt == dt.max()/2],//,//
										 "Spin tunes: Actual vs Linear Map","Data Point Index","#nu_{s}");
		plot<double> wmp((sig0*dt + eps*a + eps*V*(xi0 + cr))[dt == dt.max()/2]);//);//
		pts->Draw("AP");
		pts->SetMarkerStyle(3);
		pts->SetLineColor(0);
		wmp->Draw("P same");
		wmp->SetMarkerColor(2);
		wmp->SetMarkerStyle(2);
		wmp->SetLineColor(0);
		auto leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts.get(),"Data Point (Measured Values)");
		leg->AddEntry(wmp.get(),"Spin Tune Model (Calculated Values)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cpts->SetGrid();
		cpts->SaveAs((path + "/error_plot2.pdf").data(),"pdf");

		TCanvas* csme = new TCanvas("csme","csme",900,600);

		TH2F* hnue = new TH2F("hnue",
			"Distribution of residuals of spin-tune fit;#phi;#theta;#nu_{simulated} - #nu_{model}",
			42,-TMath::Pi(),TMath::Pi(),21,-TMath::Pi()/2,TMath::Pi()/2);
		for(int i=0; i<nu.size(); i++) hnue->Fill(phi[i],th[i],res_nu[i]);
		//hnu->Fit("gaus");
		hnue->Draw("surf1 sph");
		//csme->SaveAs((path + "/error_dist2.pdf").data(),"pdf");
		new TCanvas("clatplot2","clatplot2",1900,1000);
		plot2<double> lat_plot2(cr[0],cr[1],cr[2],"Sample Lattice",
										"#xi_{x}","#xi_{y}","#eta_{1}");
		//lat_plot.axis_equal();
		lat_plot2->Draw("AP0");
	}
}



