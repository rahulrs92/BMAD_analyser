#include<proton.hh>

array<double,2> proton_RF::spin_sigma(bool Out)
{
	if(!f_sig)
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

		matrix<2,1,valarray<double>> vA = {mdt,pow(mdt,2)};
		auto A = sum(vA * vA.tp());
		auto B = (A^-1) * sum(vA * nu);
		auto C = B.tp() * vA;

		mg["sig"] = new TMultiGraph("sig","Spin Tune Map;turn no.;#nu_{s}");

		plot<double> nu_data(valarray<double>(n.data(),N-1),nu,"#nu_{s} at given turn");
		nu_data->SetMarkerStyle(1);
		nu_data->SetLineColor(10);
		mg["sig"]->Add(nu_data.get(),"P");
		ss.str() = "";
		ss<<"#nu_{s} = "<<(B[0]<0 ? "- " : "")<<abs(B[0])<<" #delta"
												<<(B[1]<0 ? " - " : " + ")<<abs(B[1])<<" #delta^{2}"<<endl;
		string ts; getline(ss,ts);
		plot<double> nu_fit(valarray<double>(n.data(),N-1),C,ts);
		nu_fit->SetMarkerSize(1);
		nu_fit->SetLineColor(2);
		mg["sig"]->Add(nu_fit.get(),"L");
		f_sig = 1;
		sig = {B[0],B[1]};
	}
	if(Out)
	{
		TCanvas* csig = new TCanvas("csig","sigma");
		mg["sig"]->Draw("AP");
		mg["sig"]->GetYaxis()->SetTitleOffset(0.9);
		auto leg = csig->BuildLegend(0.11,0.25,0.64,0.11);
		leg->SetFillColorAlpha(10,0.2);
		csig->SetGrid();
		csig->SaveAs((path + "/sigma_fit.pdf(").data(),"pdf");
		mg["sig"]->GetXaxis()->SetRange(N>1e3 ? N/2 - 5e2 : 0,N>1e3 ? N/2 + 5e2 : N);
		csig->SaveAs((path + "/sigma_fit.pdf)").data(),"pdf");
		mg["sig"]->GetXaxis()->SetRange(0,N);
	}
	return sig;
}
