#include<lv2.hh>
#include <numeric>
#include <random>

tuple<matrix<3,1>,matrix<3,2>,matrix<2,1>> lv2::spin_model_2R(Data<double>& D,string path,bool out)
{
	valarray<double> dt = D.col("#delta_{a}");
	valarray<double> nu = D.col("#nu_{s}");
	valarray<double> ps = D.col("T");

	matrix<1,3,valarray<double>> eps =
						{D.col("#epsilon_{x}"),D.col("#epsilon_{y}"),pow(dt,2)};
	matrix<3,1,valarray<double>> tx =
						{D.col("t_{f}"),D.col("t_{d}"),D.col("t_{s}")};
	matrix<2,1,valarray<double>> cr =
						{D.col("#xi_{x}"),D.col("#xi_{y}")};

	//cout<<"b = "<<b<<endl;
	//cout<<"M = "<<M<<endl;

	/*
	matrix<9,1,valarray<double>> regs_nu;
	for(int i=0; i<3; i++) regs_nu[i] = eps[i];
	for(int i=0; i<6; i++) regs_nu[i+3] = eps[i/2] * cr[i%2];
	auto arr_V = (sum(regs_nu * regs_nu.tp())^-1) * sum(regs_nu * nu);
	for(int i=0; i<3; i++) w[i] = arr_V[i];
	for(int i=0; i<6; i++) W[i/2][i%2] = arr_V[i+3];
	
	cout<<"b = "<<w<<endl;
	cout<<"M = "<<W<<endl;
	*/

	// K-fold cross-validation parameters
	int K = 5; // You can make this a function parameter if desired
	int N = nu.size();
	vector<int> indices(N);
	std::iota(indices.begin(), indices.end(), 0);
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(indices.begin(), indices.end(), g);

	valarray<double> res_nu(N); // To store residuals for all data points
	valarray<double> model(N);  // To store model predictions for all data points
	valarray<double> model_long(N);  // To store longitudinal model predictions for all data points

	cout << "N = " << N << endl;

	w = 0;
	W = 0;

	// Perform K-fold cross-validation
	for (int k = 0; k < K; ++k)
	{
	    // Determine test indices for this fold
	    int fold_size = N / K + (k < N % K ? 1 : 0);
	    int fold_start = (N / K) * k + min(k, N % K);
	    int fold_end = fold_start + fold_size;
	    vector<int> test_idx(indices.begin() + fold_start, indices.begin() + fold_end);
	    vector<int> train_idx;
	    train_idx.reserve(N - fold_size);
	    for (int i = 0; i < N; ++i)  if (i < fold_start || i >= fold_end) train_idx.push_back(indices[i]);
		
		cout << "Fold " << k + 1 << ": Training on " << train_idx.size() << " points, Testing on " << test_idx.size() << " points." << endl;
		//cout << "Test indices: ";
		//for (const auto& idx : test_idx) cout << idx << " ";
		//cout << endl;

	    // Extract training data
	    valarray<double> dt_train(train_idx.size()), nu_train(train_idx.size());
	    matrix<1,3,valarray<double>> eps_train = valarray<double>(0.0, train_idx.size());
	    matrix<2,1,valarray<double>> cr_train = valarray<double>(0.0, train_idx.size());
	    for (size_t j = 0; j < train_idx.size(); ++j)
		{
			dt_train[j] = dt[train_idx[j]];
	        nu_train[j] = nu[train_idx[j]];
			eps_train[0][j] = eps[0][train_idx[j]];
			eps_train[1][j] = eps[1][train_idx[j]];
			cr_train[0][j] = cr[0][train_idx[j]];
			cr_train[1][j] = cr[1][train_idx[j]];
	    }
		eps_train[2] = pow(dt_train,2);
		
	    // Build regression matrix for training
	    matrix<9,1,valarray<double>> regs_nu_train;
	    for(int i=0; i<3; i++) regs_nu_train[i] = eps_train[i];
	    for(int i=0; i<6; i++) regs_nu_train[i+3] = eps_train[i/2] * cr_train[i%2];
	    auto arr_V = (sum(regs_nu_train * regs_nu_train.tp())^-1) * sum(regs_nu_train * nu_train);
	    matrix<3,1,double> w_fold;
	    matrix<3,2,double> W_fold;
	    for(int i=0; i<3; i++) w_fold[i] = arr_V[i];
	    for(int i=0; i<6; i++) W_fold[i/2][i%2] = arr_V[i+3];

		//cout<<"w_fold = "<<w_fold<<endl;
		//cout<<"W_fold = "<<W_fold<<endl;

		w = w + w_fold / K;	
		W = W + W_fold / K;

	    // Predict on test set
	    for (size_t j = 0; j < test_idx.size(); ++j)
		{
	        int idx = test_idx[j];
	        double dt_val = dt[idx];
	        double eps_x = eps[0][idx];
	        double eps_y = eps[1][idx];
	        double xi_x = cr[0][idx];
	        double xi_y = cr[1][idx];
	        matrix<1,3,double> eps_vec = {eps_x, eps_y, pow(dt_val,2)};
	        matrix<2,1,double> cr_vec = {xi_x, xi_y};
	        double pred = eps_vec*w_fold + eps_vec*W_fold*cr_vec;
	        model[idx] = pred;
			model_long[idx] = eps_vec[2]*w_fold[2] + eps_vec[2]*W_fold[2]*cr_vec;
	        res_nu[idx] = nu[idx] - pred;
	    }
	}
	double std_res_nu = pow(pow(res_nu,2).sum() / res_nu.size(),0.5);

	string path1 = path.substr(0,path.find_last_of("/"));
	fstream f(path1 + "/model_info",fstream::out | fstream::trunc);
	f<<"#nu_{s}\n#nu_{s}_model\n#nu_{s}_longitudinal_model"<<endl;
	f.close();
	f.open(path1 + "/model",fstream::out | fstream::app);
	for(int i=0; i<nu.size(); i++) f<<nu[i]<<"\t"<<model[i]<<"\t"<<model_long[i]<<endl;
	f.close();

	//cout<<"model: "<<model<<endl;

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

		cout<<"w = "<<w<<endl;
		cout<<"W = "<<W<<endl;
	}
	//cout<<"unique elements of tx[2]: "<<unique(D.col("t_{s}"))<<endl;
	xi = ((W.tp() * W)^-1) * W.tp() * (-w);
	cout<<"xi = "<<xi<<endl;
	return make_tuple(w,W,xi);
}





