#include<proton.hh>

double proton_RF::spin_tune(bool Out)
{
	if(!f_ST)
	{
		valarray<double> th(0.0,N), m(0.0,N), mdt(0.0,N-1), nu(0.0,N-1);
		for(int i=0; i<N; i++) m[i] = sx[i] / sz[i];
		m = (m.shift(1) - m) / (1 + m * m.shift(1));
		m = atan(m);
		th[0] = 0;
		for(int i=0; i<N-1; i++)
		{
			th[i+1] = th[i] + m[i];
			nu[i] = m[i]/(2*pi);
			mdt[i] = (pz[i] + pz[i+1]);
		}

		vector<double> rfi, rth;
		double pz_mean = valarray<double>(pz.data(),pz.size()).sum() / N;
		for(int i=1; i<N-1; i++)
			if((pz[i-1] - pz_mean)*(pz[i] - pz_mean) < 0)
			{
				double tt1 = pz[i]-pz[i-1];
				double tt2 = th[i]-th[i-1];
				double tq = - (pz[i-1] - pz_mean) / tt1;
				rfi.push_back(2*pi*(i-1+tq));
				rth.push_back(th[i-1] + tt2*tq);
			}
		cout<<"YOOOOO! "<<rfi.size()<<endl;
		for(int i=0; i<rfi.size()-1; i++)
		{
			rfi[i] = (rfi[i+1] + rfi[i])/2;
			rth[i] = (rth[i+1] + rth[i])/2;
		}
		rfi.pop_back(); rth.pop_back();

		if(rfi.size()==0) for(int i=0; i<N; i = i + N/100) { rfi.push_back(2*pi*n[i]) ; rth.push_back(th[i]); }

		matrix<2,1,valarray<double>> vA = {valarray<double>(rfi.data(),rfi.size()),valarray<double>(1.0,rfi.size())};
		auto A = sum(vA * vA.tp());
		auto B = (A^-1) * sum(vA * valarray<double>(rth.data(),rth.size()));
		auto C = B.tp() * vA;

		mg["nu_s"] = new TMultiGraph("nu_s",
										"Angle of spin vector at #delta=0;#phi (radians);#theta_{x} (radians)");

		plot<double> th_data(valarray<double>(n.data(),N-1)*2*pi,th,"#theta_{x} at given turn");
		th_data->SetMarkerStyle(1);
		th_data->SetLineColor(10);
		mg["nu_s"]->Add(th_data.get(),"P");
		plot<double> rth_data(rfi,rth,"#theta_{x} when #delta = 0");
		rth_data->SetMarkerStyle(2);
		rth_data->SetLineColor(10);
		mg["nu_s"]->Add(rth_data.get(),"P");
		ss.str() = "";
		ss<<"#theta_{x} = "<<(B[0]<0 ? "- " : "")<<abs(B[0])<<" #phi"
												<<(B[1]<0 ? " - " : " + ")<<abs(B[1])<<endl;
		string ts; getline(ss,ts);
		plot<double> fi_th(valarray<double>(rfi.data(),rfi.size()),C,ts);
		fi_th->SetMarkerSize(1);
		fi_th->SetLineColor(2);
		mg["nu_s"]->Add(fi_th.get(),"L");
		f_ST = 1;
		nu_s = -B[0];
	}
	if(Out)
	{
		TCanvas* cth = new TCanvas("cth","fi_th");
		mg["nu_s"]->Draw("AP");
		mg["nu_s"]->GetYaxis()->SetTitleOffset(0.9);
		auto leg = cth->BuildLegend(0.11,0.25,0.64,0.11);
		leg->SetFillColorAlpha(10,0.2);
		cth->SetGrid();
		cth->SaveAs((path + "/theta_fit.pdf").data(),"pdf");
	}
	return nu_s;
}

double proton::spin_tune(bool Out)
{
	if(!f_ST)
	{
		valarray<double> th(0.0,N), m(0.0,N), mdt(0.0,N-1), nu(0.0,N-1);
		for(int i=0; i<N; i++) m[i] = sx[i] / sz[i];
		m = (m.shift(1) - m) / (1 + m * m.shift(1));
		m = atan(m);
		th[0] = 0;
		for(int i=0; i<N-1; i++)
		{
			th[i+1] = th[i] + m[i];
			nu[i] = m[i]/(2*pi);
			mdt[i] = (pz[i] + pz[i+1]);
		}

		valarray<double> fi = 2*pi*valarray<double>(n.data(),n.size());

		matrix<2,1,valarray<double>> vA = {fi,valarray<double>(1.0,fi.size())};
		auto A = sum(vA * vA.tp());
		auto B = (A^-1) * sum(vA * th);
		auto C = B.tp() * vA;

		//B[0] -= 0.121167 - 2.33025e-07 - 2.96846e-14; // Prototype Magnetic

		mg["nu_s"] = new TMultiGraph("nu_s",
										"Angle of spin vector;#phi (radians);#theta_{x} (radians)");

		plot<double> th_data(fi,th,"#theta_{x} at given turn");
		th_data->SetMarkerStyle(1);
		th_data->SetLineColor(10);
		mg["nu_s"]->Add(th_data.get(),"P");
		ss.str() = "";
		ss<<"#theta_{x} = "<<(B[0]<0 ? "- " : "")<<abs(B[0])<<" #phi"
												<<(B[1]<0 ? " - " : " + ")<<abs(B[1])<<endl;
		string ts; getline(ss,ts);
		plot<double> fi_th(fi,C,ts);
		fi_th->SetMarkerSize(1);
		fi_th->SetLineColor(2);
		mg["nu_s"]->Add(fi_th.get(),"L");
		f_ST = 1;
		nu_s = -B[0];
	}
	if(Out)
	{
		TCanvas* cth = new TCanvas("cth","fi_th");
		mg["nu_s"]->Draw("AP");
		mg["nu_s"]->GetYaxis()->SetTitleOffset(0.9);
		auto leg = cth->BuildLegend(0.36,0.25,0.89,0.11);
		leg->SetFillColorAlpha(10,0.2);
		cth->SetGrid();
		cth->SaveAs((path + "/theta_fit.pdf").data(),"pdf");
	}
	return nu_s;
}

