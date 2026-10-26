#include<lv2.hh>
#include <numeric>
#include <random>

#include<TH1F.h>

tuple<double,matrix<3,1>,matrix<3,3>> lv2::optimize_beam(Data<double>& D,string path,bool out)
{
	fstream f;
	valarray<double> sct = D.col("#tau");
	valarray<double> R = pow(D.col("#tau"),-2);
	valarray<double> nu = D.col("#nu_{s}");
	//valarray<double> ps = D.col("#eta_{0}");   // ALWAYS CHECK LABELS

	matrix<3,1,valarray<double>> tx = {D.col("t_{f}"),D.col("t_{d}"),D.col("t_{s}")};

    // For testing with misaligned COSY: must be commented out for real data
    // matrix<10,1,valarray<double>> tx2;
	// tx2[0] = pow(tx[0],0);
	// for(int i=0; i<3; i++) tx2[i+1] = tx[i];
	// for(int i=0; i<3; i++) tx2[i+4] = pow(tx[i],2);
	// for(int i=0; i<3; i++) tx2[i+7] = tx[i]*tx[(i+1)%3];
	// valarray<double> eta1 = tensor_Xi[2]*tx2;
	// nu -= tensor_orbit * tx2;

	valarray<double> e_sct = D.col("#Delta#tau");
	valarray<double> e_R = abs(-2 * pow(D.col("#tau"),-3) * D.col("#Delta#tau"));
    valarray<double> e_nu = D.col("#Delta#nu_{s}");
	
	int N = nu.size();
	
	valarray<double> res_nu(N), res_R(N), res_sct(N); // To store residuals for all data points
	valarray<double> model_R(N), model_sct(N);  // To store model predictions for all data points
	
	cout << "N = " << N << endl;
		
	matrix<10,1,valarray<double>> regs_R, J_k;
    regs_R[0] = pow(tx[0],0);
	for(int i=0; i<3; i++) regs_R[i+1] = tx[i];
	for(int i=0; i<3; i++) regs_R[i+4] = pow(tx[i],2);
	for(int i=0; i<3; i++) regs_R[i+7] = tx[i] * tx[(i+1)%3];

    valarray<double> W_sct = pow(e_sct / sct, -2);
	cout << "W_sct sum: " << W_sct.sum() << endl;
	// W_sct = W_sct / W_sct.sum() * W_sct.size();
    valarray<double> W_R = pow(e_R, -2);
    // W_R = W_R / W_R.sum() * W_R.size();
    valarray<double> W_nu = pow(e_nu, -2);

    auto A = sum(regs_R * W_R * regs_R.tp())^(-1);
    auto B = A * sum(regs_R * W_R * R);

	res_sct = sct - 1.0 / sqrt(regs_R.tp() * B);

	auto B2M = [] (const matrix<10,1>& iB)
	{
		return matrix<3,3>{
			iB[4],   iB[7]/2, iB[9]/2,
			iB[7]/2, iB[5],   iB[8]/2,
			iB[9]/2, iB[8]/2, iB[6]
		};
	};
	
	m0 = B[0];
    m = {B[1], B[2], B[3]};
    M = B2M(B);
	
	B[0] = sqrt(max(1e-8,m0 - m.tp()*(m/M)));
	B.sub<3,1>(1,0) = - 0.5 * m / M;
	
	auto vertex_R = [&] (const matrix<10,1>& iB) -> valarray<double>
	{
		return iB[0]*iB[0] + (tx - iB.sub<3,1>(1,0)).tp() * B2M(iB) * (tx - iB.sub<3,1>(1,0));
	};

	auto regs_vertex = [&tx, &B2M] (const matrix<10,1>& iB)
	{
		matrix<10,1,valarray<double>> regs_V(pow(tx[0],0) * 0);
		auto dtx = tx - iB.sub<3,1>(1,0);
		regs_V[0] = 2 * iB[0];
		regs_V.sub<3,1>(1,0) = -2.0 * B2M(iB) * dtx;
		for(int i=0; i<3; i++) regs_V[i+4] = pow(dtx[i],2);
		for(int i=0; i<3; i++) regs_V[i+7] = (dtx[i]) * (dtx[(i+1)%3]);
		return regs_V;
	};

	int iter = 0;
	double tol_B = 1.0;
	double tol_cost = 1.0;

	while (tol_B > 1e-6 && tol_cost > 1e-8 && iter++ < 100)
	{
		auto B_k = B;
		auto R_p = vertex_R(B_k);
		// R_p[R_p < 1e-8] = 1e-8;
		valarray<double> D_k = -0.5 * pow(R_p,-1.5);
		J_k = D_k * regs_vertex(B_k);
		auto dB_k = sum(J_k * W_sct * res_sct) / sum(J_k * W_sct * J_k.tp());
		// cout<<"res_sct:\n"<<res_sct<<endl;
		// char pausech; cin>>pausech;

		// 3. Backtracking line search (Step Halving)
		double alpha = 1.0; // Step size fraction
		auto B_trial = B + alpha * dB_k;
		auto R_trial = vertex_R(B_trial);
		valarray<double> res_sct_trial = sct - 1.0/sqrt(R_trial);

		// If step pushes any R <= 0, halve alpha until all R > 0
		double cost_curr = (res_sct * W_sct * res_sct).sum();
		double cost_trial = (res_sct_trial * W_sct * res_sct_trial).sum();
		while ((B_trial[0] <= 1e-8 || R_trial.min() <= 1e-8 || std::isnan(cost_trial) || cost_trial >= cost_curr) && alpha > 1e-4)
		{
			alpha *= 0.5;
			B_trial = B + alpha * dB_k;
			R_trial = vertex_R(B_trial);
			res_sct_trial = sct - 1.0/sqrt(R_trial);
			cost_trial = (res_sct_trial * W_sct * res_sct_trial).sum(); // Re-evaluate cost
		}

		// Apply the safe update step
		B = B_trial;
		R_p = R_trial;
		// R_p[R_p < 1e-8] = 1e-8;
		res_sct = res_sct_trial;

		// tol = (alpha * dB_k).norm();
		// Relative parameter update tolerance: ||alpha * dB|| / (||B|| + 1e-12)
		tol_B = (alpha * dB_k).norm() / (B.norm() + 1e-12);
		
		// Relative cost function improvement tolerance: |cost_curr - cost_trial| / cost_curr
		tol_cost = std::abs(cost_curr - cost_trial) / cost_curr;

		cout << "loop " << iter << ", tol_B: " << tol_B << ", tol_cost: " << tol_cost << endl;
		// cout<<"loop "<<iter<<", tolerance: "<<tol<<endl;
	}
	cout<<"WNLLS converged in "<<iter<<" iterations..."<<endl;

	model_R = vertex_R(B);
    res_R = R - model_R;
	model_sct = 1.0 / sqrt(model_R);
	res_sct = sct - model_sct;
    
	m0 = B[0];
    m = {B[1], B[2], B[3]};
    M = B2M(B);
    cout<<"M = "<<M<<endl;
    cout<<"M^1 = "<<(M^-1)<<endl;
    // auto t_o = - 0.5 * (M^-1) * m;
	
    // auto cor_B = (A * sum((regs_R * W_R * res_R * res_R * W_R * regs_R.tp())) * A);
    auto cor_B = sum(J_k * W_sct * J_k.tp())^-1;
	matrix<B.M(),1> e_B = cor_B.diag();
		
    // double chi = (pow(res_R,2) * W_R).sum() / (res_R.size() - B.M());
	double chi = (pow(res_sct, 2) * W_sct).sum() / (res_sct.size() - B.M());
	double std_res_sct = sqrt((pow(res_sct, 2)).sum() / res_sct.size());

    cout<<"Chi^2/NDF = "<<chi<<endl;
    cout<<"B +/- e_B=\n"<<cat(B,e_B)<<endl;
    // cout<<"var_B = "<<e_B<<endl;
    
    // double std_res_sct = pow(pow(res_R,2).sum() / res_R.size(),0.5);
	for(int i=0; i<regs_R.M(); i++) e_B[i] = sqrt(e_B[i]);
	//cout<<"e_B = "<<e_B<<endl;
	//cout<<"res_R = "<<res_R<<endl;
	//cout<<"Standard deviation of residuals (R): "<<std_res_sct<<endl;

	//string path1 = path.substr(0,path.find_last_of("/"));
	f.open(path + "/model_info",fstream::out | fstream::trunc);
	f<<"#tau\n#tau_model"<<endl;
	f.close();
	f.open(path + "/model",fstream::out | fstream::trunc);
	for(int i=0; i<R.size(); i++) f<<sct[i]<<"\t"<<model_sct[i]<<endl;
	f.close();

	if(out)
	{
		TCanvas* cpts_sct = new TCanvas("cpts","cpts",2800,1000);

		plot<double> pts_sct(sct,//[dt == dt.max()/2],//
										 "Spin Coherence Times: Actual vs Model","Data Point Index",
									 "#tau (s)");
		plot<double> stm(model_sct);//arr_V.tp()*regs_R);//(sig0*dt + (sig0*dt);//
		pts_sct->Draw("AP");
		pts_sct->SetMarkerStyle(3);
		pts_sct->SetLineColor(0);
		stm->Draw("P same");
		stm->SetMarkerColor(2);
		stm->SetMarkerStyle(2);
		stm->SetLineColor(0);
		auto leg = new TLegend(0.3,0.9,0.7,0.8);
		leg->AddEntry(pts_sct.get(),"Data Point (Measured Values)");
		leg->AddEntry(stm.get(),"Spin Coherence Time Model (Calculated Values)");
		leg->SetFillColorAlpha(10,0.2);
		leg->Draw();
		cpts_sct->SetGrid();
		cpts_sct->SaveAs((path + "/error_plot_2.pdf").data(),"pdf");

		TCanvas* ceh_sct = new TCanvas("ceh","ceh",1024,768);
		TH1F* hsct = new TH1F("hsct",
			"Distribution of residuals of spin-coherence-time fit;#tau_{simulated} - #tau_{model};Frequency",
			51,-5*std_res_sct,5*std_res_sct);
		for(int i=0; i<res_R.size(); i++) hsct->Fill(res_sct[i]);
		//hnu->Fit("gaus");
		hsct->Draw();
		ceh_sct->SetGrid();
		ceh_sct->SaveAs((path + "/error_dist_2.pdf").data(),"pdf");

		// cout<<"w = "<<w<<endl;
		// cout<<"W3 = "<<W3<<endl;
	}
	//cout<<"unique elements of tx[2]: "<<unique(D.col("t_{s}"))<<endl;
	//auto xi_opt = (W3^-1) * (-w);
    
	// cout<<"Vertex: \n"<<matrix<3,2>{t_o[0],t_o[1],t_o[2],e_t_o[0][0],e_t_o[1][1],e_t_o[2][2]}<<endl;
	// cout<<"Vertex: \n"<<t_o<<endl;
    cout<<"Vertex +/- uncertainty:\n"<<cat(m,e_B.sub<3,1>(1,0))<<endl;
    // double ist_o = m0 + m.tp() * t_o + t_o.tp() * M * t_o;
    // cout<<"Spin Coherence Time at vertex: "<<(ist_o<0 ? -pow(-ist_o,-0.5) : pow(ist_o,-0.5))<<endl;
	cout<<"Spin Coherence Time at vertex: "<<1.0 / m0<<" +/- "<<e_B[0] / (m0 * m0)<<endl;
	cout<<"SCT's at segment passing through vertex (step size 0.001):"<<endl;
	for(double i=-0.01; i<=0.01; i+=0.001) cout<<1.0/sqrt(m0*m0 + (m - m + i).tp() * M * (m - m + i))<<" ";
	cout<<endl;
    return make_tuple(m0,m,M);
}





