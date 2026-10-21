#include<lv2.hh>
#include <numeric>
#include <random>

tuple<matrix<3,1>,matrix<3,3>,matrix<3,1>> lv2::spin_model_tensor(Data<double>& D,string path,bool out)
{
	valarray<double> dt = D.col("#delta_{a}");
	valarray<double> nu = D.col("#nu_{s}");
	//valarray<double> ps = D.col("#eta_{0}");   // ALWAYS CHECK LABELS

	matrix<1,3,valarray<double>> eps =
						{D.col("#epsilon_{x}"),D.col("#epsilon_{y}"),pow(dt,2)};
	matrix<3,1,valarray<double>> tx =
						{D.col("t_{f}"),D.col("t_{d}"),D.col("t_{s}")};
	valarray<double> eta1 = xi03[2] + Xi3[2]*tx;
	matrix<3,1,valarray<double>> cr =
						{D.col("#xi_{x}"),D.col("#xi_{y}"), eta1};

	//cout<<"b = "<<b<<endl;
	//cout<<"M = "<<M<<endl;

	fstream f(path + "/RFxi_out_info", fstream::out | fstream::trunc);
    f<<"t_{f}"<<"\n"<<"t_{d}"<<"\n"<<"t_{s}"<<"\n";
    f<<"#xi_{x}"<<"\n"<<"#xi_{y}"<<"\n"<<"#eta_{1}"<<"\n";
    f.close();
	f.open(path + "/RFxi_out", fstream::out | fstream::trunc);
	for(int i=0; i<tx[0].size(); i++)
    {
        for(int j=0; j<3; j++) f<<tx[j][i]<<"\t";
        for(int j=0; j<2; j++) f<<cr[j][i]<<"\t"; f<<cr[2][i]<<endl;
	}
	f.close();
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

	matrix<21,1,valarray<double>> regs_nu;
		for(int i=0; i<3; i++) regs_nu[i] = eps[i];
		for(int i=0; i<9; i++) regs_nu[i+3] = eps[i/3] * cr[i%3];
		for(int i=0; i<3; i++) regs_nu[i+12] = tx[i];
		for(int i=0; i<3; i++) regs_nu[i+15] = pow(tx[i],2);
		for(int i=0; i<3; i++) regs_nu[i+18] = tx[i]*tx[(i+1)%3];
	/*
		matrix<9,1,valarray<double>> regs_nu;
		for(int i=0; i<3; i++) regs_nu[i] = tx[i];
		for(int i=0; i<3; i++) regs_nu[i+3] = pow(tx[i],2);
		for(int i=0; i<3; i++) regs_nu[i+6] = tx[i]*tx[(i+1)%3];
	*/
	    
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
		
	    // Build regression matrix for training
	    auto regs_nu_train = pull(regs_nu, valarray<size_t>(train_idx.data(), train_idx.size()));
		valarray<double> nu_train = nu[valarray<size_t>(train_idx.data(), train_idx.size())];
	    auto arr_V = (sum(regs_nu_train * regs_nu_train.tp())^-1) * sum(regs_nu_train * nu_train);
	    matrix<3,1,double> w_fold, tw_fold;
	    matrix<3,3,double> W_fold;
	    for(int i=0; i<3; i++) w_fold[i] = arr_V[i];
	    for(int i=0; i<9; i++) W_fold[i/3][i%3] = arr_V[i+3];
		//for(int i=0; i<3; i++) tw_fold[i] = arr_V[i+12];

		//cout<<"w_fold = "<<w_fold<<endl;
		//cout<<"W_fold = "<<W_fold<<endl;
		//cout<<"tw_fold = "<<tw_fold<<endl;

		w = w + w_fold * (fold_size / double(N));	
		W3 = W3 + W_fold * (fold_size / double(N));
		tw = tw + tw_fold * (fold_size / double(N));

	    // Predict on test set
		valarray<double> nu_test = nu[valarray<size_t>(test_idx.data(), test_idx.size())];
		valarray<double> model_test = arr_V.tp() * pull(regs_nu, valarray<size_t>(test_idx.data(), test_idx.size()));
		model[valarray<size_t>(test_idx.data(), test_idx.size())] = model_test;
		res_nu[valarray<size_t>(test_idx.data(), test_idx.size())] = nu_test - model_test;
	}
	double std_res_nu = pow(pow(res_nu,2).sum() / res_nu.size(),0.5);
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

		cout<<"w = "<<w<<endl;
		cout<<"W3 = "<<W3<<endl;
	}
	//cout<<"unique elements of tx[2]: "<<unique(D.col("t_{s}"))<<endl;
	auto xi_opt = (W3^-1) * (-w);
	cout<<"xi_opt = "<<xi_opt<<endl;
	cout<<"t_opt = "<<(Xi3^-1) * (xi_opt-xi03)<<endl;
	//cout<<"t_opt = "<<xi03 + Xi3*xi_opt<<endl;
	return make_tuple(w,W3,xi_opt);
}





