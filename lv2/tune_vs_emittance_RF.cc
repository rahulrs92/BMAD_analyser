#include<lv2.hh>

matrix<3,1> lv2::tune_vs_emittance_RF(string path,bool out)
{
    Data<double> D(path + "/line");
    valarray<double> nu = D.col("#nu_{s}");
    valarray<double> ep_x = D.col("#epsilon_{x}");
    valarray<double> ep_y = D.col("#epsilon_{y}");
    valarray<double> dt = D.col("#delta_{0}");
    valarray<double> ps = D.col("T");
    valarray<double> dn = D.col(0);

    matrix<2,1,valarray<double>> regs_ep = {ep_x,ep_y};//,ep_y*pow(dt,2),pow(dt,2)*ep_x,pow(dt,4)};
    auto A = (sum(regs_ep * regs_ep.tp())^-1) * sum(regs_ep * nu);
    //auto B = A.tp() * regs_ep;
    cout<<"A = "<<A<<endl;
    //cout<<"regs_ep = "<<regs_ep<<endl;
    
    if(out)
    {
        int n = nu.size();
        cout<<"n = "<<n<<endl;
        double d_x = (ep_x.max() - ep_x.min()) / (n);
        double d_y = (ep_y.max() - ep_y.min()) / (n);
        valarray<double> sm_epx(pow(n,2)), sm_epy(pow(n,2));

        for(int i=0; i<n; i++)
        {
            for(int j=0; j<n; j++)
            {
                sm_epx[i*n + j] = ep_x.min() + i * d_x;
                sm_epy[i*n + j] = ep_y.min() + j * d_y;
            }
        }
        matrix<2,1,valarray<double>> sm_ep = {sm_epx,sm_epy};

        valarray<double> fit_vals = A.tp()*regs_ep;
        valarray<double> ydata = nu;
        //cout<<"fit_vals: "<<fit_vals<<endl;
        //cout<<"ydata: "<<ydata<<endl;

        plot2<double> fit2(sm_epx,sm_epy,A.tp()*sm_ep,
            "Spin tune vs. Transverse Emittances @ " + 
            proton_RF::mode + " (" + proton_RF::ID + ");#epsilon_{x} (m#upointrad);#epsilon_{y} (m#upointrad);#nu_{s}");
        fit2->SetMarkerColorAlpha(0,0);
        fit2->SetLineColor(2);

        plot2<double> graph2(ep_x,ep_y,ydata,"#nu_{s}");
        graph2->SetMarkerStyle(20);
        graph2->SetMarkerColor(4);
        graph2->SetLineColor(4);
        
        TCanvas *cnuep = new TCanvas("canvas2", "Spin tune vs. Emittance", 800, 600);
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

        //cout<<"fit_vals: "<<fit_vals<<endl;
        //cout<<"ydata: "<<ydata<<endl;

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
        cnuep->SaveAs((path + "/ep_correlation.pdf").data(), "pdf");
    }
    w = {A[0], A[1], 0};
    return w;
}