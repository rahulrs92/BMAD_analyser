#include<proton.hh>

array<double,2> proton_RF::betatron_tune(bool iy,bool Out)
{
	string key = (iy ? "y" : "x");
	array<double,2> &Qm = (iy ? Qym : Qxm), &e_Qm = (iy ? e_Qym : e_Qxm);
	bool &f_BT = (iy ? f_BTy : f_BTx);
	if(!f_BT)
	{
		stringstream ss;
		bool &f_CS = (iy ? f_CSy : f_CSx);
		auto CS = Courant_Snyder(iy,Out);
		vector<double> &u = (iy ? y : x), &pu = (iy ? py : px), un(N), pun(N);
		// Normalize position and momentum to Courant-Snyder phase space coordinates
		for(int i=0; i<N; i++)
		{
			un[i] = u[i] / pow(CS[2],0.5);
			pun[i] = (CS[1]*u[i] + CS[2]*pu[i]) / pow(CS[2],0.5);
		}
		// Initialize arrays for momentum offset and betatron tune samples
		valarray<double> mdt(0.0,N-2), Q(0.0,N-2);
		// Loop over each turn (excluding first and last) to compute local tune from trajectory curvature
		for(int i=1; i<N-1; i++)
		{
			// Extract three consecutive phase space points (previous, current, next)
			matrix<3,1> a = {un[i-1],pun[i-1],0}, b = {un[i],pun[i],0}, c = {un[i+1],pun[i+1],0};
			// Build coefficient matrix from differences between consecutive points
			matrix<2,2> A = {a[0] - b[0], a[1] - b[1], b[0] - c[0], b[1] - c[1]};
			// Build right-hand side from squared distances of points from origin
			matrix<2,1> d = {a.tp()*a - b.tp()*b, b.tp()*b - c.tp()*c};
			// Solve for circle center 'o' passing through the three points
			matrix<3,1> o; o.sub<2,1>() = 0.5 * d/A;

			// Lambda function to compute angle between two vectors
			auto dphi = [](matrix<3,1> p, matrix<3,1> q)
			{
				double ct = (p.tp()*q) / (p.norm() * q.norm());
				double st = (p*q).norm() / (p.norm() * q.norm());
				// Return clockwise angle accounting for full [0, 2π] range
				return (st>0 ? 2*pi - acos(ct) : acos(ct));
			};
			
			// Compute betatron tune as average rotation angle divided by 2π (normalized)
			Q[i-1] = (dphi(a-o, b-o) + dphi(b-o, c-o)) / (4*pi);
			// Store momentum offset for this turn
			mdt[i-1] = pz[i];
		}
		string ytt = (iy ? "Q_{y}" : "Q_{x}");

		matrix<3,1,valarray<double>> vA = {pow(mdt,0),mdt,pow(mdt,2)};
		auto A = sum(vA * vA.tp())^-1;
		auto B = A * sum(vA * Q);

		valarray<double> res_Q = Q - (Tps + B[0] + B[1]*mdt + B[2]*pow(mdt,2));
		auto E = A * sum((vA * res_Q * res_Q * vA.tp())) * A;

		if(Out) cout<<"Estimated Chromaticity: "<<B[1]<<" +/- "<<E[1][1]<<endl;

		// mg[("Q" + key)] = new TMultiGraph(("Q" + key).data(),
		// 								("Normalized x vs p_{x} @ " + mode + " (" + ID + ");x (metres);" +
		// 								"p_{x}").data());

		// plot<double> Q_data(vu / pow(CS[2],0.5),(CS[1]*vu + CS[2]*vpu) / pow(CS[2],0.5),"phase space points");
		// Q_data->SetMarkerStyle(1);
		// Q_data->SetLineColor(10);
		// mg[("Q" + key)]->Add(Q_data.get(),"P");

		mg[("Q" + key)] = new TMultiGraph(("Q" + key).data(),
										("Betatron Tune vs. Phase advance @ " + mode + " (" + ID + ");Turn no.;" +
											ytt).data());

		plot<double> Q_data(valarray<double>(n.data(),N-2)+1,Q,ytt + " at given turn");
		Q_data->SetMarkerStyle(1);
		Q_data->SetLineColor(10);
		mg[("Q" + key)]->Add(Q_data.get(),"P");

		ss.clear();
		ss<<"Q_{"<<key<<"} ="
		<<(B[0]<0 ? " - " : " ")<<abs(B[0])
		<<(B[1]<0 ? " - " : " + ")<<abs(B[1])<<" #delta"
		<<(B[2]<0 ? " - " : " + ")<<abs(B[2])<<" #delta^{2}"<<endl;
		string ts; getline(ss,ts);

		plot<double> Q_fit(valarray<double>(n.data(),N-2)+1,Tps + B[0] + B[1]*mdt + B[2]*pow(mdt,2),ts);
		Q_fit->SetMarkerStyle(1);
		Q_fit->SetLineColor(2);
		mg[("Q" + key)]->Add(Q_fit.get(),"L");

		mg[("Q" + key + "_fit")] = new TMultiGraph(("Q" + key + "_fit").data(),
										("Betatron Tune vs. Momentum Offset @ " + mode + " (" + ID + ");#delta;" +
											ytt).data());

		plot<double> Qd_data(mdt,Q,ytt + " at given turn");
		Qd_data->SetMarkerStyle(1);
		Qd_data->SetLineColor(10);
		mg[("Q" + key + "_fit")]->Add(Qd_data.get(),"P");

		plot<double> Qd_fit(mdt,B[0] + B[1]*mdt + B[2]*pow(mdt,2),ts);
		Qd_fit->SetMarkerStyle(1);
		Qd_fit->SetLineColor(2);
		mg[("Q" + key + "_fit")]->Add(Qd_fit.get(),"L");

		f_BT = 1;
		Qm[0] = B[0];
		e_Qm[0] = sqrt(E[0][0]);
		Qm[1] = B[1];
		e_Qm[1] = sqrt(E[1][1]);
	}
	if(Out)
	{
		TCanvas* cQ = new TCanvas(("cQ" + key).data(),"fi_Q");
		mg[("Q" + key)]->Draw("AP");
		mg[("Q" + key)]->GetYaxis()->SetTitleOffset(1.2);
		auto leg = cQ->BuildLegend(0.11,0.25,0.64,0.11);
		leg->SetFillColorAlpha(10,0.2);
		cQ->SetGrid();
		cQ->SaveAs((path + "/Q" + key + ".pdf(").data(),"pdf");
		mg[("Q" + key)]->GetXaxis()->SetRange(N>1e3 ? N/2 - 5e2 : 0,N>1e3 ? N/2 + 5e2 : N);
		cQ->SaveAs((path + "/Q" + key + ".pdf)").data(),"pdf");
		mg[("Q" + key)]->GetXaxis()->SetRange(0,N);

		TCanvas* cQ_fit = new TCanvas(("cQ" + key + "_fit").data(),"fi_Q_fit");
		mg[("Q" + key + "_fit")]->Draw("AP");
		cQ_fit->SetGrid();
		mg[("Q" + key + "_fit")]->GetYaxis()->SetTitleOffset(0.9);
		auto leg2 = cQ_fit->BuildLegend(0.11,0.25,0.64,0.11);
		leg2->SetFillColorAlpha(10,0.2);
		cQ_fit->SaveAs((path + "/Q" + key + "_fit.pdf").data(),"pdf");
	}
	return Qm;
}

