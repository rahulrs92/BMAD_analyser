#include<lv2.hh>

void lv2::tune_vs_sextupoles(string path,bool out)
{
    Data<double> D(path + "/line");
    valarray<double> nu = D.col("#nu_{s}");
    //valarray<double> ps = D.col("T");
    valarray<double> dt = D.col("#delta_{0}");
    valarray<double> dn = D.col(0);
    matrix<1,3,valarray<double>> eps = {D.col("#epsilon_{x}"), D.col("#epsilon_{y}"), pow(dt,2)};
    matrix<2,1,valarray<double>> cr = {D.col("#xi_{x}"), D.col("#xi_{y}")};
    matrix<2,1,valarray<double>> t = {D.col("t_{f}"), D.col("t_{d}")};

    //matrix<2,1> xi0;
    //matrix<2,2> Xi;
    tie(xi0,Xi_2D) = windmap_2D(D,path,out);

    //matrix<4,1,valarray<double>> regs_ep = {pow(ps,0),ps,ep_x,ep_y};
    //auto A = (sum(regs_ep * regs_ep.tp())^-1) * sum(regs_ep * nu);
    //auto B = A.tp() * regs_ep;

    matrix<10,1,valarray<double>> regs_nu;
    regs_nu[0] = dt;
    for(int i=0; i<3; i++) regs_nu[i+1] = eps[i];
	for(int i=0; i<6; i++) regs_nu[i+4] = eps[i/2] * cr[i%2];
	auto arr_H = (sum(regs_nu * regs_nu.tp())^-1) * sum(regs_nu * (nu));
    sig0 = arr_H[0];
    for(int i=0; i<3; i++) h[i] = arr_H[i+1];
	for(int i=0; i<6; i++) H[i/2][i%2] = arr_H[i+4];

    matrix<3,1> A = {h[0], H[0][0], H[0][1]};

    int n = pow(nu.size(),0.5);
    valarray<double> sm_tf(pow(n,2)), sm_td(pow(n,2));

    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
        {
            sm_tf[i*n + j] = t[0][i];
            sm_td[i*n + j] = t[1][j];
        }

    //cout<<"sm_tf = "<<sm_tf<<endl;
    //cout<<"sm_td = "<<sm_td<<endl;

    valarray<double> sm_cx = xi0[0] + Xi_2D[0][0]*sm_tf + Xi_2D[0][1]*sm_td;
    valarray<double> sm_cy = xi0[1] + Xi_2D[1][0]*sm_tf + Xi_2D[1][1]*sm_td;

    valarray<double> res_nu = nu - (sig0*dt + eps*h + eps*H*cr);
    double std_res_nu = pow(pow(res_nu,2).sum() / res_nu.size(),0.5);
    cout<<"standard deviation of residuals = "<<std_res_nu<<endl;

    if(out)
    {
        
        valarray<double> ydata = nu - (h[1] + H[1][0]*cr[0] + H[1][1]*cr[1])*eps[1]
                                    - (h[2] + H[2][0]*cr[0] + H[2][1]*cr[1])*eps[2];
        ydata -= sig0*dt;
        //ydata /= eps[0];
        
        //valarray<double> ydata = A[0] + A[1]*cr[0] + A[2]*cr[1] + res_nu;
        
        plot2<double> fit2(sm_tf,sm_td,A[0] + A[1]*sm_cx + A[2]*sm_cy,
            "h_{1} vs. Sextupole field strengths @ " + 
            proton_RF::mode + " (" + proton_RF::ID + ");t_{f};t_{d};h_{1}");
        fit2->SetMarkerColorAlpha(0,0);
        fit2->SetLineColor(2);

        plot2<double> graph2(t[0],t[1],ydata,"#nu_{s}");
        graph2->SetMarkerStyle(20);
        graph2->SetMarkerColor(4);
        graph2->SetLineColor(4);
        
        TCanvas *cnuep = new TCanvas("canvas2", "h1 vs. Sextupole field strengths", 800, 600);
        TPad* pad1 = new TPad("pad1", "Main plot", 0, 0.25, 1, 1.0);
        pad1->SetBottomMargin(0.07);
        pad1->Draw();
        TPad* pad2 = new TPad("pad2", "Residuals", 0, 0.0, 1, 0.25);
        pad2->SetTopMargin(0.1);
        pad2->SetBottomMargin(0.3);
        pad2->Draw();

        pad1->cd();
        fit2->Draw("surf1");
        graph2->Draw("P0 same");
        fit2->GetXaxis()->SetTitleOffset(2);
        fit2->GetYaxis()->SetTitleOffset(2);

        pad2->cd();
        valarray<double> fit_vals = A[0] + A[1]*cr[0] + A[2]*cr[1];
        valarray<double> residuals = ydata - fit_vals;
        valarray<double> zero(0.0,residuals.size());
        plot<double> gres(dn+1, residuals, "");
        plot<double> gzero(dn+1, zero, "Zero Line");
        gres->Sort();
        gres->SetMarkerStyle(7);
        gres->SetLineColor(2);
        gres->GetXaxis()->SetTitle("Data Point Index");
        gres->GetYaxis()->SetTitle("Residuals");
        gres->GetYaxis()->SetTitleSize(0.12);
        gres->GetYaxis()->SetLabelSize(0.10);
        gres->GetYaxis()->SetTitleOffset(0.3);
        gres->GetXaxis()->SetTitleSize(0.12);
        gres->GetXaxis()->SetLabelSize(0.10);
        gres->GetXaxis()->SetLimits(0 - 0.01*dn.size(), dn.size() + 0.01*dn.size());
        //gStyle->SetTitleFontSize(0.1);
        gres->Draw("APL");
        gzero->SetLineColor(1);
        gzero->SetLineWidth(2);
        gzero->Draw("L");
        pad2->SetGrid();

        cnuep->cd();
        cnuep->SaveAs((path + "/cr_correlation.pdf").data(), "pdf");
    }
}