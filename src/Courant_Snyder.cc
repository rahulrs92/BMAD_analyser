#include<proton.hh>

array<double,7> proton_RF::Courant_Snyder(bool iy, bool Out)
{
	bool &f_CS = (iy ? f_CSy : f_CSx);
	string key = (iy ? "CSy" : "CSx");
	array<double,7> &CS = (iy ? CSy : CSx), &e_CS = (iy ? e_CSy : e_CSx);
	double m_pz = 0; for(auto a : pz) m_pz += a / pz.size();
	//cout<<"m_pz = "<<m_pz<<endl;

	if(!f_CS)
	{
		vector<double> &u = (iy ? y : x),&pu = (iy ? py : px),v_rx,v_rpx;
		for(int i=1; i<N; i++)
		{
			if((pz[i-1] - m_pz)*(pz[i] - m_pz) < 0)
			{
				v_rx.push_back(u[i-1]); v_rx.push_back(u[i]);
				v_rpx.push_back(pu[i-1]); v_rpx.push_back(pu[i]);
			}
		}
		valarray<double> rx(v_rx.data(),v_rx.size()), rpx(v_rpx.data(),v_rpx.size());

		if(v_rx.size()<6)
		{
			rx = valarray<double>(u.data(),u.size());
			rpx = valarray<double>(pu.data(),pu.size());
		}

		matrix<5,1,valarray<double>> A({pow(rx,2),rx*rpx,pow(rpx,2),rx,rpx});
		auto sA = sum(A);
		auto B = A*A.tp();
		auto sB = sum(B);
		
		//cout<<sBi<<endl;
		if(det(sB)==0)
		{
			cerr<<string(iy ? "Vertical" : "Horizontal") + " emittance is incalculable!"<<endl;
			f_CS = 1;
			return CS;
		}
		else
		{
			auto sBi = sB^(-1);
			auto C = sBi*sA;
			double ell = C[0]*C[2] - 0.25*C[1]*C[1];
			if(ell<0)
			{
				cerr<<"Fitting error: Hyperbola detected! (1/em1^2 = "<<ell<<")"<<endl;
			}
			double em1 = C[2]/abs(C[2]) * pow(abs(ell),-0.5);
			//cout<<"Wrong emittance: "<<em1<<endl;
			CS[1] = C[1]*em1/2;		//alpha
			CS[2] = C[2]*em1;		//beta
			CS[3] = C[0]*em1;		//gamma
			auto tC = - ((matrix<2,2>{CS[3],CS[1],CS[1],CS[2]}^-1) * matrix<2,1>{C[3],C[4]}) * (em1/2);
			CS[4] = tC[0];			//x offset
			CS[5] = tC[1];			//px offset
			CS[0] = em1 + CS[3]*pow(CS[4],2) + 2*CS[1]*CS[4]*CS[5] + CS[2]*pow(CS[5],2);
			CS[6] = pow(u[0],2) * CS[3];
			if(CS[0]<0)
			{
				cerr<<"Fitting error: Negative Emittance detected! ("<<CS[0]<<")"<<endl;
				CS = {-CS[0],-CS[1],-CS[2],-CS[3],CS[4],CS[5],-CS[6]};
				return CS;
			}

			if(Out) cout<<string(iy ? "Vertical" : "Horizontal") + " emittance: "<<CS[0]<<endl;
			if(Out) cout<<"CS_"<<(iy ? "y" : "x")<<": "<<CS[1]<<", "<<CS[2]<<", "<<CS[3]<<endl;
			if(Out) cout<<(iy ? "y" : "x")<<"-orbit offsets: "<<CS[4]<<", "<<CS[5]<<endl;

			valarray<double> res_u = 1 - A.tp()*C;
			auto E = sBi * sum((A * res_u * res_u * A.tp())) * sBi;
			auto J = Courant_Snyder_Jacobian(C);
			auto E_J = J * E * J.tp();
			valarray<double> sig = E_J.diag();
			e_CS[0] = pow(sig[5],0.5);
			for(int i=0; i<5; i++) e_CS[i+1] = pow(sig[i],0.5);
			e_CS[6] = (pow(u[0],2) * e_CS[3]);
			/*
			if(Out)
			{
				cout<<"Coefficients: \n"<<C<<endl;
				cout<<"ell = "<<ell<<", wrong emittance = "<<em1<<endl;
				cout<<"Coefficient variances of Ellipse Fit: "<<endl;
				cout<<E<<endl;
				cout<<"Coefficient variances of Courant-Snyder parameters: "<<endl;
				cout<<E_J<<endl;
				cout<<"Estimated uncertainties (std. dev.): "<<endl;
				cout<<string(iy ? "Vertical" : "Horizontal") + " emittance: "<<e_CS[0]<<endl;
				cout<<"CS_"<<(iy ? "y" : "x")<<": "<<e_CS[1]<<", "<<e_CS[2]<<", "<<e_CS[3]<<endl;
				cout<<(iy ? "y" : "x")<<"-orbit offsets: "<<e_CS[4]<<", "<<e_CS[5]<<endl<<endl;
			}
			*/
		}
		if(!f_tree || Out)
		{
			ss.str() = "";
			ss<<(iy ? "Vertical" : "Horizontal")<<" Phase-Space;"
				<<(iy ? "y;" : "x;")<<(iy ? "p_{y}" : "p_{x}");
			string ts; getline(ss,ts);
			mg[key] = new TMultiGraph(key.data(),ts.data());

			matrix<2,2> mcs = {CS[2],0,-CS[1],-1};
			valarray<double> fi(100);
			for(int i=0; i<100; i++) fi[i] = i*2*pi/100;
			matrix<2,1,valarray<double>> ph = {cos(fi),sin(fi)};
			matrix<2,1> u_off = {CS[4],CS[5]};
			auto uu = u_off - pow(CS[0]/CS[2],0.5) * (mcs * ph);

			plot<double> u_pu(u,pu,"Phase-space points");
			u_pu->SetMarkerStyle(1);
			u_pu->SetMarkerColor(1);
			u_pu->SetLineColor(10);
			mg[key]->Add(u_pu.get(),"P");

			if(v_rx.size()>=6)
			{
				plot<double> z_pz(rx,rpx,"Points at #delta = 0");
				z_pz->SetMarkerStyle(3);
				z_pz->SetLineColor(10);
				mg[key]->Add(z_pz.get(),"P");
			}

			ss.clear();
			ss.clear();
			ss<<CS[0]<<" = "<<(CS[3]<0 ? "- " : "")<<abs(CS[3])<<(iy ? " y'^{2}" : " x'^{2}")
											<<(CS[1]<0 ? " - " : " + ")<<abs(CS[1])<<(iy ? " y'p_{y}'" : " x'p_{x}'")
											<<(CS[2]<0 ? " - " : " + ")<<abs(CS[2])<<(iy ? " p_{y}'^{2}" : " p_{x}'^{2}")
											<<endl;
			ts = ""; getline(ss,ts);
			plot<double> z_pz_fit(uu[0],uu[1],ts);
			z_pz_fit->SetMarkerStyle(0);
			z_pz_fit->SetMarkerColor(2);
			z_pz_fit->SetLineColor(2);
			mg[key]->Add(z_pz_fit.get(),"L");
		}
		f_CS = 1;
	}
	if(Out)
	{
		TCanvas* c1 = new TCanvas(iy ? "y_py" : "x_px",iy ? "y_py" : "x_px");
		mg[key]->Draw("AP");
		mg[key]->GetYaxis()->SetTitleOffset(0.9);
		c1->SetGrid();
		ss<<(iy ? "y'" : "x'")<<" = "<<(iy ? "y" : "x")
			<<(CS[4]<0 ? " - " : " + ")<<abs(CS[4])
			<<", "<<(iy ? "p_{y}'" : "p_{x}'")<<" = "<<(iy ? "p_{y}" : "p_{x}")
			<<(CS[5]<0 ? " - " : " + ")<<abs(CS[5])<<endl;
		string ts; getline(ss,ts);
		auto leg = c1->BuildLegend(CS[1]<0 ? 0.36 : 0.11,0.25,CS[1]<0 ? 0.89 : 0.64,0.11);//(0.7,0.9,0.9,0.8);
		leg->AddEntry((TObject*)0,ts.data(),"");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		c1->SaveAs((path + (iy ? "/y_py.pdf" : "/x_px.pdf")).data(),"pdf");
		//c1->Close();
	}
	return CS;
}

array<double,7> proton::Courant_Snyder(bool iy, bool Out)
{
	bool &f_CS = (iy ? f_CSy : f_CSx);
	string key = (iy ? "CSy" : "CSx");
	array<double,7> &CS = (iy ? CSy : CSx);

	if(!f_CS)
	{
		valarray<double> u((iy ? y : x).data(),(iy ? y : x).size());
		valarray<double> pu((iy ? py : px).data(),(iy ? py : px).size());

		matrix<5,1,valarray<double>> A = {pow(u,2),u*pu,pow(pu,2),u,pu};
		auto sA = sum(A);
		auto B = A*A.tp();
		auto sBi = sum(B)^(-1);

		if(det(sBi)==0)
		{
			cerr<<string(iy ? "Vertical" : "Horizontal") + " emittance is incalculable!"<<endl;
			f_CS = 1;
			return CS;
		}

		auto C = sBi*sA;
		double ell = C[0]*C[2] - 0.25*C[1]*C[1];
		if(ell<0)
		{
			cerr<<"Fitting error: Hyperbola detected! (1/em1^2 = "<<ell<<")"<<endl;
		}
		double em1 = C[2]/abs(C[2]) * pow(abs(ell),-0.5);
		//cout<<"Wrong emittance: "<<em1<<endl;
		CS[1] = C[1]*em1/2;
		CS[2] = C[2]*em1;
		CS[3] = C[0]*em1;
		auto tC = -((matrix<2,2>{CS[3],CS[1],CS[1],CS[2]}^-1) * matrix<2,1>{C[3],C[4]})*(em1/2);
		CS[4] = tC[0];
		CS[5] = tC[1];
		CS[0] = em1 + CS[3]*pow(CS[4],2) + 2*CS[1]*CS[4]*CS[5] + CS[2]*pow(CS[5],2);
		CS[6] = pow(u[0],2) * CS[3];
		if(CS[0]<0)
		{
			cerr<<"Fitting error: Negative Emittance detected! ("<<CS[0]<<")"<<endl;
			CS = {-CS[0],-CS[1],-CS[2],-CS[3],CS[4],CS[5],-CS[6]};
			return CS;
		}

		if(Out) cout<<string(iy ? "Vertical" : "Horizontal") + " emittance: "<<CS[0]<<endl;
		if(Out) cout<<"CS_"<<(iy ? "y" : "x")<<": "<<CS[1]<<", "<<CS[2]<<", "<<CS[3]<<endl;
		if(Out) cout<<(iy ? "y" : "x")<<"-orbit offsets: "<<CS[4]<<", "<<CS[5]<<endl;

		if(!f_tree || Out)
		{
			//cout<<"YAY!!!"<<endl;
			ss.str() = "";
			ss<<(iy ? "Vertical" : "Horizontal")<<" Phase-Space;"
				<<(iy ? "y;" : "x;")<<(iy ? "p_{y}" : "p_{x}")<<endl;
			string ts; getline(ss,ts);
			mg[key] = new TMultiGraph(key.data(),ts.data());

			matrix<2,2> mcs = {CS[2],0,-CS[1],-1};
			valarray<double> fi(100);
			for(int i=0; i<100; i++) fi[i] = i*2*pi/100;
			matrix<2,1,valarray<double>> ph = {cos(fi),sin(fi)};
			matrix<2,1> u_off = {CS[4],CS[5]};
			auto uu = u_off - (pow(CS[0]/CS[2],0.5) * (mcs * ph));
			//cout<<u_off<<endl;
			//cout<<uu<<endl;

			plot<double> u_pu(u,pu,"Phase-space points");
			u_pu->SetMarkerStyle(0);
			u_pu->SetLineColor(10);
			mg[key]->Add(u_pu.get(),"P");

			ss.clear();
			ss<<CS[0]<<" = "<<(CS[3]<0 ? "- " : "")<<abs(CS[3])<<(iy ? " y'^{2}" : " x'^{2}")
											<<(CS[1]<0 ? " - " : " + ")<<abs(CS[1])<<(iy ? " y'p_{y}'" : " x'p_{x}'")
											<<(CS[2]<0 ? " - " : " + ")<<abs(CS[2])<<(iy ? " p_{y}'^{2}" : " p_{x}'^{2}")
											<<endl;
			ts = ""; getline(ss,ts);
			plot<double> z_pz_fit(uu[0],uu[1],ts);
			z_pz_fit->SetMarkerStyle(0);
			z_pz_fit->SetMarkerColor(2);
			z_pz_fit->SetLineColor(2);
			mg[key]->Add(z_pz_fit.get(),"L");
		}
		f_CS = 1;
	}
	if(Out)
	{
		TCanvas* c1 = new TCanvas(iy ? "y_py" : "x_px",iy ? "y_py" : "x_px");
		mg[key]->Draw("AP");
		mg[key]->GetYaxis()->SetTitleOffset(0.9);
		c1->SetGrid();
		ss.str() = "";
		ss<<(iy ? "y'" : "x'")<<" = "<<(iy ? "y" : "x")
			<<(CS[4]<0 ? " - " : " + ")<<abs(CS[4])
			<<", "<<(iy ? "p_{y}'" : "p_{x}'")<<" = "<<(iy ? "p_{y}" : "p_{x}")
			<<(CS[5]<0 ? " - " : " + ")<<abs(CS[5])<<endl;
		string ts; getline(ss,ts);
		auto leg = c1->BuildLegend(CS[1]<0 ? 0.36 : 0.11,0.25,CS[1]<0 ? 0.89 : 0.64,0.11);//(0.7,0.9,0.9,0.8);
		leg->AddEntry((TObject*)0,ts.data(),"");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		c1->SaveAs((path + (iy ? "/y_py.pdf" : "/x_px.pdf")).data(),"pdf");
		//c1->Close();
	}
	return isnan(CS[0]) ? array<double,7>{0,0,0,0,0,0,0} : CS;
}

