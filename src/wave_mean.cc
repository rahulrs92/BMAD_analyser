#include<proton.hh>

pair<array<double,2>,valarray<double>> proton_RF::wave_mean(vector<double> &ix, vector<double> &iw)
{
	long N = iw.size();
	valarray<double> w(iw.data(),N), x(ix.data(),N), cw(N);
	double md1 = w.sum() / N;
	cw[0] = w[0];
	for(int i=1; i<N; i++) cw[i] = cw[i-1] + (x[i]-x[i-1])*w[i];
	valarray<double> w1 = w - md1;

	vector<double> rx, rw;
	for(int i=2; i<=N; i++)
		if(w1[i-1]*w1[i] < 0)
		{
			double tt1 = w1[i]-w1[i-1];
			double tt2 = x[i]-x[i-1];
			double tt3 = cw[i]-cw[i-1];
			double tq = - w1[i-1] / tt1;
			rx.push_back(x[i-1] + tt2*tq);
			rw.push_back(cw[i-1] + tt3*tq);
		}
	for(int i=0; i<rx.size()-1; i++)
	{
		rx[i] = (rx[i+1] + rx[i])/2;
		rw[i] = (rw[i+1] + rw[i])/2;
	}
	rx.pop_back();
	rw.pop_back();

	TCanvas* c1 = new TCanvas("t_cw","t_cw");
	plot<double>(x,cw,"Cumulative Wavefunction at zero displacement","x","#Sigmaw")->Draw();
	plot<double> t_cw(rx,rw);
	t_cw->Fit("pol1");
	t_cw->SetMarkerStyle(3);
	t_cw->Draw("P same");
	c1->SetGrid();

	double md = t_cw->GetFunction("pol1")->GetParameters()[1];

	valarray<double> w2 = w - md;
	double da = 0;
	for(int i=0; i<N; i++) da += pow(w2[i],2);
	da = pow(2 * da / N,0.5);
	//cout<<"da = "<<da<<endl;

	return make_pair(array<double,2>{md,da},valarray<double>(rx.data(),rx.size()));
}
