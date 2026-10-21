#include<lv2.hh>
#include <numeric>
#include <random>
#include<TH1F.h>


matrix<1,10> lv2::spin_orbit_correction(Data<double>& D,string path,bool out)
{
	fstream f;
	valarray<double> nu = D.col("#nu_{s}");
	//valarray<double> ps = D.col("#eta_{0}");   // ALWAYS CHECK LABELS

	matrix<3,1,valarray<double>> tx = {D.col("t_{f}"),D.col("t_{d}"),D.col("t_{s}")};

	matrix<10,1,valarray<double>> tx2;
	tx2[0] = pow(tx[0],0);
	for(int i=0; i<3; i++) tx2[i+1] = tx[i];
	for(int i=0; i<3; i++) tx2[i+4] = pow(tx[i],2);
	for(int i=0; i<3; i++) tx2[i+7] = tx[i]*tx[(i+1)%3];
	
	valarray<double> e_nu_s = D.col("#Delta#nu_{s}");
	valarray<double> w = pow(e_nu_s,-2);

	// K-fold cross-validation parameters
	int K = 5; // You can make this a function parameter if desired
	int N = nu.size();
	vector<size_t> indices(N);
	std::iota(indices.begin(), indices.end(), 0);
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(indices.begin(), indices.end(), g);

	valarray<double> res_nu(N); // To store residuals for all data points
	valarray<double> model(N);  // To store model predictions for all data points

	cout << "N = " << N << endl;

	matrix<10,1,valarray<double>> regs_nu;
	for(int i=0; i<10; i++) regs_nu[i] = tx2[i];

	matrix<regs_nu.M(),1> R, e_R;
	    
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
		
	    // Build regression matrix for training
	    auto regs_nu_train = pull(regs_nu, valarray<size_t>(train_idx.data(), train_idx.size()));
		valarray<double> nu_train = nu[valarray<size_t>(train_idx.data(), train_idx.size())];
		valarray<double> w_train = w[valarray<size_t>(train_idx.data(), train_idx.size())];

		// Weighted linear regression to find coefficients
		auto arr_A = sum(regs_nu_train * w_train * regs_nu_train.tp())^-1;
	    auto arr_R = arr_A * sum(regs_nu_train * w_train * nu_train);
		
		R += arr_R * (fold_size / double(N));
		valarray<double> res_nu_train = nu_train - arr_R.tp() * regs_nu_train;
		e_R += (arr_A * sum((regs_nu_train * w_train * res_nu_train * res_nu_train * w_train * regs_nu_train.tp())) * arr_A).diag();
		e_R *= pow(fold_size / double(N),2);
		
	    // Predict on test set
		valarray<double> nu_test = nu[valarray<size_t>(test_idx.data(), test_idx.size())];
		valarray<double> model_test = arr_R.tp() * pull(regs_nu, valarray<size_t>(test_idx.data(), test_idx.size()));
		model[valarray<size_t>(test_idx.data(), test_idx.size())] = model_test;
		res_nu[valarray<size_t>(test_idx.data(), test_idx.size())] = nu_test - model_test;
	}
	double std_res_nu = pow(pow(res_nu,2).sum() / res_nu.size(),0.5);

	//string path1 = path.substr(0,path.find_last_of("/"));
	f.open(path + "/model_info",fstream::out | fstream::trunc);
	f<<"#nu_{s}\n#nu_{s}_model"<<endl;
	f.close();
	f.open(path + "/model",fstream::out | fstream::trunc);
	for(int i=0; i<nu.size(); i++) f<<nu[i]<<"\t"<<model[i]<<endl;
	f.close();

	f.open(path + "/orbit_tensor_info",fstream::out | fstream::trunc);
	f<<"#nu_{s0}\nr[0]\nr[1]\nr[2]\nR[0][0]\nR[1][1]\nR[2][2]\nR[0][1]\nR[1][2]\nR[2][0]"<<endl;
	f.close();
	f.open(path + "/orbit_tensor",fstream::out | fstream::trunc);
	f<<R<<endl;
	f.close();

	f.open(path + "/orbit_tensor_error_info",fstream::out | fstream::trunc);
	f<<"#Delta#nu_{s0}\ne_r[0]\ne_r[1]\ne_r[2]\ne_R[0][0]\ne_R[1][1]\ne_R[2][2]\ne_R[0][1]\ne_R[1][2]\ne_R[2][0]"<<endl;
	f.close();
	f.open(path + "/orbit_tensor_error",fstream::out | fstream::trunc);
	f<<e_R<<endl;
	f.close();

	if(out)
	{
		TCanvas* cpts_nu = new TCanvas("cpts","cpts",2800,1000);

		plot<double> pts_nu(nu,"Spin tunes: Actual vs Model","Data Point Index", "#nu_{s}");
		plot<double> stm(model);
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
	}
	return R.tp();
}





