#include<lv2.hh>
#include <numeric>
#include <random>

#include<TH1F.h>

tuple<matrix<3,1>,matrix<3,3>,matrix<3,1>> lv2::spin_model_3(Data<double>& D,string path,bool out)
{
	fstream f;
	valarray<double> dt = D.col("#delta_{a}");
	valarray<double> dt_e = D.col("#Delta#delta_{a}");
	
	valarray<double> nu = D.col("#nu_{s}");
	valarray<double> nu_e = D.col("#Delta#nu_{s}");
	cout<<"YAY"<<endl;
	int N = nu.size();

	double eta0 = D.col("T").sum() / nu.size();   // ALWAYS CHECK LABELS

	matrix<1,3,valarray<double>> eps = {D.col("#epsilon_{x}"),D.col("#epsilon_{y}"),pow(dt,2)/2};
	matrix<1,3,valarray<double>> ep_e = {D.col("#Delta#epsilon_{x}"),D.col("#Delta#epsilon_{y}"),dt*dt_e};
	matrix<3,1,valarray<double>> tx = {D.col("t_{f}"),D.col("t_{d}"),D.col("t_{s}")};

	matrix<3,1,valarray<double>> cr = xi03 + Xi3*tx;
	matrix<(cr.M()+1)*cr.M(),cr.M(),valarray<double>> J_cr = valarray<double>(0.0,nu.size());
	for(int i=0; i<cr.M(); i++) { J_cr(i*(cr.M()+1),i) = pow(nu,0); J_cr.sub<3,1>(i*(cr.M()+1)+1,i) = tx; }
	cout << "J_cr = \n" << pull(J_cr,1) << endl;
	auto cov_cr = J_cr.tp() * Cov_Xi * J_cr;
	
	valarray<double> res_nu(N); // To store residuals for all data points
	valarray<double> model(N), model_no_off(N), nu_no_off(N);  // To store model predictions for all data points
	//valarray<double> model_long(N);  // To store longitudinal model predictions for all data points
	
	cout << "N = " << N << endl;
		
	matrix<12,1,valarray<double>> regs_nu;
	for(int i=0; i<3; i++) regs_nu[i] = eps[i];
	for(int i=0; i<9; i++) regs_nu[i+3] = eps[i/3] * cr[i%3];

	auto cov2err = [] (auto C)
	{
		matrix<C.M(),1> o1;
		for(int i=0; i<C.M(); i++) o1[i] = sqrt(C(i,i));
		return o1;
	};
	auto par2mat = [] (auto C)
	{
		matrix<C.M()/4,C.M()/4> res;
		for(int i=0; i<C.M()/4; i++) for(int j=0; j<C.M()/4; j++) res[i][j] = C[(i+1)*(C.M()/4) + j];
		return res;
	};
	    
	matrix<regs_nu.M(),1> C;
	matrix<regs_nu.M(),regs_nu.M()> Cov_C;

    // K-fold cross-validation parameters
	int K = 5; // You can make this a function parameter if desired
	vector<size_t> indices(N);
	std::iota(indices.begin(), indices.end(), 0);
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(indices.begin(), indices.end(), g);
	
	// Perform K-fold cross-validation
	for (int k = 0; k < K; ++k)
	{
		// Determine test indices for this fold
	    int fold_size = N / K + (k < (N % K) ? 1 : 0);
	    int fold_start = (N / K) * k + min(k, N % K);
	    int fold_end = fold_start + fold_size;
	    vector<size_t> test_idx(indices.begin() + fold_start, indices.begin() + fold_end);
	    vector<size_t> train_idx;
	    train_idx.reserve(N - fold_size);
	    for (int i = 0; i < N; ++i)  if (i < fold_start || i >= fold_end) train_idx.push_back(indices[i]);
		
		auto train_indices = valarray<size_t>(train_idx.data(), train_idx.size());
		
		cout << "Fold " << k + 1 << ": Training on " << train_idx.size() << " points, Testing on " << test_idx.size() << " points." << endl;
		
	    // Build regression matrix for training
        auto regs_nu_train = pull(regs_nu, train_indices);
		auto cr_train = pull(cr, train_indices);
		auto ep_train = pull(eps,train_indices);
		auto cov_ep_train = pull(eye(ep_e.v())^2,train_indices);
		// cout << "cov_ep_train =\n" << pull(cov_ep_train,1) << endl;
		auto cov_cr_train = pull(cov_cr,train_indices);
		// cout << "cov_cr_train =\n" << pull(cov_cr_train,fold_start) << endl;
		valarray<double> nu_train = nu[train_indices];
		valarray<double> nu_e_train = nu_e[train_indices];

	    auto A_train = sum(regs_nu_train * regs_nu_train.tp())^-1;
		auto C_train = A_train * sum(regs_nu_train * nu_train);

		auto J_ep = C_train.sub<3,1>() + par2mat(C_train)*cr_train;
		auto J_xi = ep_train * par2mat(C_train);
		
		valarray<double> w_train = pow(nu_e_train,2) + J_ep.tp()*cov_ep_train*J_ep + J_xi*cov_cr_train*J_xi.tp();
		w_train = 1.0 / w_train;
		
		auto Cov_C_train = sum(regs_nu_train * w_train * regs_nu_train.tp())^-1;
		C_train = Cov_C_train * sum(regs_nu_train * w_train * nu_train);
		
		// C += C_train * (fold_size / double(N));
		// Cov_C += Cov_C_train * pow(fold_size / double(N),2);
		
	    // Predict on test set
		valarray<double> nu_test = nu[valarray<size_t>(test_idx.data(), test_idx.size())];
		valarray<double> model_test = C_train.tp() * pull(regs_nu, valarray<size_t>(test_idx.data(), test_idx.size()));
		model[valarray<size_t>(test_idx.data(), test_idx.size())] = model_test;
		res_nu[valarray<size_t>(test_idx.data(), test_idx.size())] = nu_test - model_test;
	}
	double std_res_nu = pow(pow(res_nu,2).sum() / res_nu.size(),0.5);

	auto A = sum(regs_nu * regs_nu.tp())^-1;
	C = A * sum(regs_nu * nu);

	auto J_ep = C.sub<3,1>() + par2mat(C)*cr;
	auto J_xi = ep * par2mat(C);
	
	valarray<double> W = pow(nu_e,2) + J_ep.tp()*(eye(ep_e.v())^2)*J_ep + J_xi*cov_cr*J_xi.tp();
	W = 1.0 / W;
	
	Cov_C = sum(regs_nu * W * regs_nu.tp())^-1;
	C = Cov_C * sum(regs_nu * W * nu);

	cout << "Performance metrics: " << endl;
	auto RMSE = sqrt(pow(res_nu,2).sum() / res_nu.size());
	cout << "Root Mean Squared Error: " << RMSE << endl;
	auto R2 = pow(res_nu,2).sum() / pow(nu - nu.sum()/nu.size(),2).sum();
	cout << "1 - Coefficient of Determination: " << R2 << endl;
	auto rChi2 = pow(res_nu/nu_e,2).sum() / (res_nu.size() - C.M());
	cout << "Reduced Chi-Squared: " << rChi2 << endl;
	cout << endl;
	
	f.open(path + "/model_info",fstream::out | fstream::trunc);
	f<<"#nu_{s}\n#nu_{s}_model\n#nu_{s}_model_long"<<endl;
	f.close();
	f.open(path + "/model",fstream::out | fstream::trunc);
	for(int i=0; i<nu.size(); i++) f<<nu[i]<<"\t"<<model[i]<<endl;
	f.close();

	if(out)
	{
		TCanvas* cpts_nu = new TCanvas("cpts","cpts",2800,1000);

		plot<double> pts_nu(nu,//[dt == dt.max()/2],//
										 "Spin tunes: Actual vs Model","Data Point Index",
									 "#nu_{s}");
		plot<double> stm(model);//arr_V.tp()*regs_nu);//(sig0*dt + (sig0*dt);//
		pts_nu->Draw("AP");
		pts_nu->SetMarkerStyle(3);
		pts_nu->SetLineColor(0);
		stm->Draw("P same");
		stm->SetMarkerColor(2);
		stm->SetMarkerStyle(2);
		stm->SetLineColor(0);
		auto leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts_nu.get(),"Data Point (Measured Values)");
		leg->AddEntry(stm.get(),"Spin Tune Model (Calculated Values)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cpts_nu->SetGrid();
		cpts_nu->SaveAs((path + "/error_plot_2.pdf").data(),"pdf");

		TCanvas* ceh_nu = new TCanvas("ceh","ceh",1024,768);
		TH1F* hnu = new TH1F("hnu",
			"Distribution of normalized residuals of spin-tune fit;#frac{#nu_{simulated} - #nu_{model}}{#sigma_{#nu}};Frequency",
			101,-2*sqrt(rChi2),2*sqrt(rChi2));
		for(int i=0; i<res_nu.size(); i++) hnu->Fill(res_nu[i] / nu_e[i]);
		// hnu->Fit("gaus");
		hnu->Draw();
		ceh_nu->SaveAs((path + "/error_dist_2.pdf").data(),"pdf");
	}
    w = C.sub<3,1>();
	W3 = par2mat(C);
	
	auto xi_opt = - w / W3;

	matrix<3,1> e_xi_opt;
	auto sig_xi_w = Cov_C.sub<3,3>().divide_left(W3) / W3.tp();
	matrix<3,12> J_sig_xi_wW;
	J_sig_xi_wW.sub<3,3>() = (W3^-1);
	for(int i=0; i<W3.M(); i++) J_sig_xi_wW.sub<3,3>(0,(i+1)*3) = xi_opt[i] * J_sig_xi_wW.sub<3,3>();
	auto sig_xi_W3 = J_sig_xi_wW * Cov_C * J_sig_xi_wW.tp();
    auto cov_xi_opt = sig_xi_w + sig_xi_W3;
	for(int i=0; i<w.M(); i++) e_xi_opt[i] = sqrt(cov_xi_opt(i,i));
	cout<<"xi_opt +/- e_xi_opt =\n"<<cat(xi_opt,e_xi_opt)<<endl;
	
    auto t_o = (xi_opt - xi03) / Xi3;

	auto sig_t_xo = cov_xi_opt.divide_left(Xi3) / Xi3.tp();
	matrix<3,12> J_sig_t_xX;
	J_sig_t_xX.sub<3,3>() = (Xi3^-1);
	for(int i=0; i<Xi3.M(); i++) J_sig_t_xX.sub<3,3>(0,(i+1)*3) = t_o[i] * J_sig_xi_wW.sub<3,3>();
	auto sig_t_xX = J_sig_t_xX * Cov_Xi * J_sig_t_xX.tp();
	auto cov_t_o = sig_t_xo + sig_t_xX;
    auto e_t_o = cov2err(cov_t_o);
	cout<<"t_opt +/- e_t_opt =\n"<<cat(t_o,e_t_o)<<endl;

	cout<<"Spin chromaticities at optimal point: \n"<<(w + W3*xi_opt)<<endl;
	
	cout<<"Spin tune shift coefficients at origin: \n"<<w<<endl;
	cout<<"W3 (optimized) =\n"<<W3<<endl;

	matrix<3,3> W3_t;
	W3_t[0][0] = W3[0][0] * (eta0*183.473/M_PI);
	W3_t[1][1] = W3[1][1] * (eta0*183.473/M_PI);
	W3_t[2][2] = - W3[2][2] * eta0;
	cout<<"sigma_0 estimations =\n"<<W3_t<<endl;

	sampler_sphere(100,1.0,t_o,"bbo_sample");

	return make_tuple(w,W3,xi_opt);
}





