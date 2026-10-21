#include<lv2.hh>
#include <numeric>
#include <random>
#include<TH1F.h>


tuple<matrix<3,1>,matrix<3,3>,matrix<3,1>> lv2::spin_model_tensor(Data<double>& D,string path,bool out)
{
	fstream f;
	valarray<double> dt = D.col("#delta_{a}");
	valarray<double> nu = D.col("#nu_{s}");
	//valarray<double> ps = D.col("#eta_{0}");   // ALWAYS CHECK LABELS

	matrix<1,3,valarray<double>> eps = {D.col("#epsilon_{x}"),D.col("#epsilon_{y}"),pow(dt,2)};
	matrix<3,1,valarray<double>> tx = {D.col("t_{f}"),D.col("t_{d}"),D.col("t_{s}")};

	matrix<10,1,valarray<double>> tx2;
	tx2[0] = pow(tx[0],0);
	for(int i=0; i<3; i++) tx2[i+1] = tx[i];
	for(int i=0; i<3; i++) tx2[i+4] = pow(tx[i],2);
	for(int i=0; i<3; i++) tx2[i+7] = tx[i]*tx[(i+1)%3];
	
	valarray<double> eta1 = tensor_Xi[2]*tx2;
	nu -= tensor_orbit * tx2;

	matrix<3,1,valarray<double>> cr = {D.col("#xi_{x}"),D.col("#xi_{y}"), eta1};

	matrix<4,1,valarray<double>> ers = 
				{D.col("#Delta#epsilon_{x}"),D.col("#Delta#epsilon_{y}"),2*dt*D.col("#Delta#delta_{a}"),
								D.col("#Delta#nu_{s}")};

	matrix<4,4,valarray<double>> Q;
	for(int i=0; i<4; i++) Q[i][i] = ers[i] * ers[i];
	
	
	// K-fold cross-validation parameters
	int K = 4; // You can make this a function parameter if desired
	int N = nu.size();
	vector<size_t> indices(N);
	std::iota(indices.begin(), indices.end(), 0);
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(indices.begin(), indices.end(), g);
	
	valarray<double> res_nu(N); // To store residuals for all data points
	valarray<double> model(N), model_no_off(N), nu_no_off(N);  // To store model predictions for all data points
	//valarray<double> model_long(N);  // To store longitudinal model predictions for all data points
	
	cout << "N = " << N << endl;
	
	w = 0;
	W3 = 0;
	matrix<3,1,double> tw;
	
	matrix<30,1,valarray<double>> regs_nu;
	//for(int i=0; i<9; i++) regs_nu[i] = tx2[i+1];
	for(int i=0; i<3; i++) for(int j=0; j<10; j++) regs_nu[i*10 + j] = eps[i] * tx2[j];
	/*
	matrix<21,1,valarray<double>> regs_nu;
	for(int i=0; i<3; i++) regs_nu[i] = eps[i];
	for(int i=0; i<9; i++) regs_nu[i+3] = eps[i/3] * cr[i%3];
	for(int i=0; i<3; i++) regs_nu[i+12] = tx[i];
	for(int i=0; i<3; i++) regs_nu[i+15] = pow(tx[i],2);
	for(int i=0; i<3; i++) regs_nu[i+18] = tx[i]*tx[(i+1)%3];
	
	matrix<9,1,valarray<double>> regs_nu;
	for(int i=0; i<3; i++) regs_nu[i] = tx[i];
	for(int i=0; i<3; i++) regs_nu[i+3] = pow(tx[i],2);
	for(int i=0; i<3; i++) regs_nu[i+6] = tx[i]*tx[(i+1)%3];
	*/

	matrix<regs_nu.M(),1> V, e_V;
	    
	// Perform K-fold cross-validation
	for (int k = 0; k < K; ++k)
	{
		// Determine test indices for this fold
	    int fold_size = N / K + (k < N % K ? 1 : 0);
	    int fold_start = (N / K) * k + min(k, N % K);
	    int fold_end = fold_start + fold_size;
	    vector<size_t> test_idx(indices.begin() + fold_start, indices.begin() + fold_end);
	    vector<size_t> train_idx;
	    train_idx.reserve(N - fold_size);
	    for (int i = 0; i < N; ++i)  if (i < fold_start || i >= fold_end) train_idx.push_back(indices[i]);
		
		cout << "Fold " << k + 1 << ": Training on " << train_idx.size() << " points, Testing on " << test_idx.size() << " points." << endl;
		//cout << "Test indices: ";
		//for (const auto& idx : test_idx) cout << idx << " ";
		//cout << endl;
		
		auto train_indices = valarray<size_t>(train_idx.data(), train_idx.size());

	    // Build regression matrix for training
	    auto regs_nu_train = pull(regs_nu, train_indices);
		valarray<double> nu_train = nu[train_indices];
	    auto arr_V = (sum(regs_nu_train * regs_nu_train.tp())^-1) * sum(regs_nu_train * nu_train);
		
		matrix<arr_V.M(),1> d_arr_V = 1.0;
		matrix<arr_V.M(),1,valarray<double>> A;
		valarray<double> w_train, res_nu_train;
		int iter = 0;

		// Iteratively Reweighted Least Squares (IRLS)
		while(d_arr_V.norm() > 1e-8)
		{
			matrix<4,1,valarray<double>> ers_train = pull(ers, train_indices), J;
			
			matrix<4,4,valarray<double>> Q = valarray<double>(0.0,train_idx.size());
			for(int i=0; i<4; i++) Q[i][i] = ers_train[i] * ers_train[i];
			
			for(int i=0; i<3; i++) J[i] = - arr_V.sub<10,1>(10*i,0).tp() * pull(tx2,train_indices);
			J[3] = pow(nu_train,0);
			
			A = - regs_nu_train;
			res_nu_train = nu_train - (arr_V.tp() * regs_nu_train);
			// for(int i=0; i<4; i++) for(int j=0; j<4; j++) cout<<"sizes of Q: \n"<<Q[i][j].size()<<endl;
			w_train = 1.0 / (J.tp() * Q * J);
			// w_train = w_train / w_train.sum() * train_idx.size();
			// cout<<"Sum of weights: "<<w_train.sum()<<endl;
			d_arr_V = (sum(A * w_train * A.tp())^-1) * sum(A * w_train * res_nu_train);
			arr_V -= d_arr_V;
			iter++;
			cout<<"IRLS iteration "<<iter<<", delta V norm = "<<d_arr_V.norm()<<endl;
			if(iter > 100) break;
		}
		if(iter > 100) cout<<"IRLS did not converge in 100 iterations!"<<endl;
		else cout<<"IRLS converged in "<<iter<<" iterations."<<endl;

		auto cor_arr_V = sum(A * w_train * A.tp())^-1;
		matrix<arr_V.M(),1> e_arr_V = cor_arr_V.diag();
		
		res_nu_train = nu_train - arr_V.tp() * regs_nu_train;
		double chi = (pow(res_nu_train,2) * w_train).sum() / (res_nu_train.size() - arr_V.M());

		cout<<"Chi^2/NDF = "<<chi<<endl;
		cout<<"arr_V = "<<arr_V.tp()<<endl;
		cout<<"var_arr_V = "<<e_arr_V.tp()<<endl;
		
	    // matrix<3,1,double> w_fold, tw_fold;
	    // matrix<3,3,double> W_fold;
	    // for(int i=0; i<3; i++) w_fold[i] = arr_V[i];
	    // for(int i=0; i<9; i++) W_fold[i/3][i%3] = arr_V[i+3];
		// for(int i=0; i<3; i++) tw_fold[i] = arr_V[i+12];
		
		// cout<<"w_fold = "<<w_fold<<endl;
		// cout<<"W_fold = "<<W_fold<<endl;
		// cout<<"tw_fold = "<<tw_fold<<endl;
		V += arr_V * (fold_size / double(N));
		e_V += e_arr_V * pow(fold_size / double(N),2);
		
		// w = w + w_fold * (fold_size / double(N));	
		// W3 = W3 + W_fold * (fold_size / double(N));
		// tw = tw + tw_fold * (fold_size / double(N));

		
	    // Predict on test set
		valarray<double> nu_test = nu[valarray<size_t>(test_idx.data(), test_idx.size())];
		valarray<double> model_test = arr_V.tp() * pull(regs_nu, valarray<size_t>(test_idx.data(), test_idx.size()));
		model[valarray<size_t>(test_idx.data(), test_idx.size())] = model_test;
		res_nu[valarray<size_t>(test_idx.data(), test_idx.size())] = nu_test - model_test;
	}
	double std_res_nu = pow(pow(res_nu,2).sum() / res_nu.size(),0.5);
	for(int i=0; i<regs_nu.M(); i++) e_V[i] = pow(e_V[i],0.5);
	//cout<<"e_V = "<<e_V<<endl;
	//cout<<"res_nu = "<<res_nu<<endl;
	//cout<<"Standard deviation of residuals (nu): "<<std_res_nu<<endl;

	//string path1 = path.substr(0,path.find_last_of("/"));
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
			"Distribution of residuals of spin-tune fit;#nu_{simulated} - #nu_{model};Frequency",
			51,-5*std_res_nu,5*std_res_nu);
		for(int i=0; i<res_nu.size(); i++) hnu->Fill(res_nu[i]);
		//hnu->Fit("gaus");
		hnu->Draw();
		ceh_nu->SaveAs((path + "/error_dist_2.pdf").data(),"pdf");

		// cout<<"w = "<<w<<endl;
		// cout<<"W3 = "<<W3<<endl;
	}
	//cout<<"unique elements of tx[2]: "<<unique(D.col("t_{s}"))<<endl;
	//auto xi_opt = (W3^-1) * (-w);
	matrix<3,10> tensor_V;
	matrix<10,1,valarray<double>> t_o2;
	matrix<10,1> t_o_l2;
	for(int i=0; i<3; i++) for(int j=0; j<10; j++) tensor_V[i][j] = V[i*10 + j];
	auto t_o = solve_tensor_system(tensor_V, {0,0,0});
	cout<<"tensor_V = \n"<<tensor_V<<endl;
	cout<<"t_opt = "<<t_o<<endl;
	t_o2 = {pow(t_o[0],0), t_o[0], t_o[1], t_o[2],
			pow(t_o[0],2), pow(t_o[1],2), pow(t_o[2],2),
			t_o[0]*t_o[1], t_o[1]*t_o[2], t_o[2]*t_o[0]};
	cout<<"tensor_Xi = "<<tensor_Xi<<endl;
	auto xi_opt = tensor_Xi * t_o2;
	cout<<"xi_opt = "<<xi_opt<<endl;

	cout<<"Spin chromaticities at optimal point: \n"<<(tensor_V * t_o2)<<endl;
	cout<<"Spin chromaticities at origin: \n"<<(tensor_V * eye<10>()(0)).norm()<<endl;
	auto t_o_l = - (tensor_V.sub<3,3>(0,1)^-1) * (tensor_V(0));
	t_o_l2 = {pow(t_o_l[0],0), t_o_l[0], t_o_l[1], t_o_l[2],
			pow(t_o_l[0],2), pow(t_o_l[1],2), pow(t_o_l[2],2),
			t_o_l[0]*t_o_l[1], t_o_l[1]*t_o_l[2], t_o_l[2]*t_o_l[0]};
	cout<<"Spin chromaticities at linear solution: \n"<<(tensor_V * t_o_l2).norm()<<endl;
	//cout<<"t_opt = "<<(Xi3^-1) * (xi_opt-xi03)<<endl;
	//cout<<"t_opt = "<<xi03 + Xi3*xi_opt<<endl;
	return make_tuple(w,W3,pull(xi_opt,0));
}





