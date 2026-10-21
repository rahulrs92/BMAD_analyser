#include<proton.hh>

double proton_RF::SCT(bool Out)
{
	if(!f_SCT)
	{
		stringstream ss;
		valarray<double> mag(0.0,N);
		for(int i=0; i<N; i++) mag[i] = pow(pow(sx[i],2) + pow(sy[i],2) + pow(sz[i],2),0.5);
		
		double to_t = L[mode] / (bt * cel);
		cout<<"bt = "<<bt<<endl;
		cout<<"to_t = "<<to_t<<endl;
		
		vector<double> rt, rmag;
		/*
		double pz_mean = valarray<double>(pz.data(),pz.size()).sum() / pz.size();
		
		for(int i=1; i<N-1; i++)
			if((pz[i-1] - pz_mean)*(pz[i] - pz_mean) < 0)
			{
				//cout<<"i: "<<i<<endl;
				double tt1 = pz[i]-pz[i-1];
				double tt2 = mag2[i]-mag2[i-1];
				double tq = - (pz[i-1] - pz_mean) / tt1;
				rt.push_back((i-1+tq) * to_t);
				rmag.push_back(mag2[i-1] + tt2*tq);
			}
		//N=N/10;
		for(int i=1; i<N-1; i++)
			if(z[i-1]*z[i] < 0)
			{
				double tt1 = z[i]-z[i-1];
				double tt2 = mag2[i]-mag2[i-1];
				double tq = - z[i-1] / tt1;
				rt.push_back((i-1+tq) * to_t);
				rmag.push_back(mag2[i-1] + tt2*tq);
			}
		cout<<"Size of magnitude subset with zero momentum offset: "<<rt.size()<<endl;
			
		for(int i=0; i<rt.size()-1; i++)
		{
			rt[i] = (rt[i+1] + rt[i])/2;
			rmag[i] = (rmag[i+1] + rmag[i])/2;
		}
		rt.pop_back(); rmag.pop_back();
		
		if(rt.size()==0) for(int i=0; i<N; i = i + N/100) { rt.push_back(n[i] * to_t); rmag.push_back(mag2[i]); }
		
		auto vrt = valarray<double>(rt.data(),rt.size()), vrth = valarray<double>(rmag.data(),rmag.size());
		*/
		valarray<double> vt = valarray<double>(n.data(),n.size()) * 2*pi; // * to_t; //

		matrix<2,1,valarray<double>> vA = {pow(vt,0),pow(vt,2)};
		auto A = (sum(vA * vA.tp())^-1);
		auto B = A * sum(vA * mag);
		auto C = B.tp() * vA;

		sct = pow((exp(-1.0) - B[0]) / B[1],0.5);

		valarray<double> res_mag = mag - C;
		auto E = A * sum((vA * res_mag * res_mag * vA.tp())) * A;
		e_sct = - ((exp(-1.0) - B[0]) / pow(B[1],2) * pow(E[1][1],0.5) + 1.0 / B[1] * pow(E[0][0],0.5)) / (2*sct);

		sct *= to_t / (2*pi);
		e_sct *= to_t / (2*pi);

		// Uncomment to flip the decoherence curve
		// mag = 1 - mag;
		// A[0] = 1 - A[0];
		// A[1] = -A[1];

		mg["SCT"] = new TMultiGraph("SCT",
			("Polarization Magnitude vs Phase Advance @ " + mode + 
				" (" + ID + ");#phi (radians);#left|#vec{P}(#phi)#right|^{2}").data());

		plot<double> mag_data(vt,mag,"#left|#vec{P}(#phi)#right| at given turn");
		mag_data->SetMarkerStyle(1);
		mag_data->SetMarkerColor(12);
		mag_data->SetLineColor(10);
		mg["SCT"]->Add(mag_data.get(),"P");
		/*
		plot<double> rth_data(rt,rmag,"#left|P(x)#right| when #delta = 0");
		rth_data->SetMarkerStyle(2);
		rth_data->SetLineColor(10);
		//mg["SCT"]->Add(rth_data.get(),"P");
		*/
		ss.str() = "";
		ss<<"#left|#vec{P}(#phi)#right| = "<<(B[1]<0 ? "- " : "")<<abs(B[1])<<" #phi^{2}"
			<<(B[0]<0 ? " - " : " + ")<<abs(B[0])<<endl;
		string ts; getline(ss,ts);
		cout<<"ts: "<<ts<<endl;
		plot<double> fi_th(vt,B.tp()*vA,ts);
		fi_th->SetMarkerColorAlpha(0,0);
		fi_th->SetLineColor(2);
		mg["SCT"]->Add(fi_th.get(),"L");
		f_SCT = 1;
	}
	if(Out)
	{
		cout<<"Spin Coherence Time: "<<sct<<" seconds"<<endl;
		TCanvas* csct = new TCanvas("csct","t_mag");
		mg["SCT"]->Draw("AP");
		mg["SCT"]->GetYaxis()->SetTitleOffset(0.9);
		auto leg = csct->BuildLegend(0.5,0.93,0.96,0.7);
		leg->AddEntry((TObject*)0,("Spin Coherence Time: " + to_string(sct) + " seconds").data(),"");
		leg->AddEntry((TObject*)0,("@ #vec{#xi} = (" + to_string(Cx()) + ", " + to_string(Cy()) + ")").data(),"");
		if(mode != "Hybrid") leg->AddEntry((TObject*)0,("@ t_{s} = " + to_string(Kt3())).data(),"");
		leg->SetFillColorAlpha(10,1.0);
		csct->SetGrid();
		csct->SaveAs((path + "/SCT_fit.pdf").data(),"pdf");
	}
	return sct;
}


