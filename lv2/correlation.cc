#include <lv2.hh>
#include <TRatioPlot.h>
#include <TStyle.h>

void lv2::correlation(string path, bool pdf)
{
    Data<double> D(path + "/line");
    valarray<double> ps = D.col("T");
    valarray<double> nu = D.col("#nu_{s}");
    valarray<double> ep_x = D.col("#epsilon_{x}");
    valarray<double> ep_y = D.col("#epsilon_{y}");
    valarray<double> dt = D.col("#delta_{0}");
    matrix<1,3,valarray<double>> eps = {ep_x, ep_y, pow(dt,2)};
    matrix<2,1,valarray<double>> cr = {D.col("#xi_{x}"), D.col("#xi_{y}")};

    //cout<<"T = "<<x<<endl;
    //cout<<"#nu_{s} = "<<y<<endl;

    matrix<2,1,valarray<double>> vA = {pow(ps,0), pow(ps,1)};// , ep_x, ep_y};
    auto A = (sum(vA * vA.tp())^-1) * sum(vA * nu);

    /*
    matrix<10,1,valarray<double>> regs_ps;
    regs_ps[0] = dt;
    for(int i=0; i<3; i++) regs_ps[i+1] = eps[i];
	for(int i=0; i<6; i++) regs_ps[i+4] = eps[i/2] * cr[i%2];
	auto arr_U = (sum(regs_ps * regs_ps.tp())^-1) * sum(regs_ps * (nu));
    eta0 = arr_U[0];
    for(int i=0; i<3; i++) h[i] = arr_U[i+1];
	for(int i=0; i<6; i++) H[i/2][i%2] = arr_U[i+4];
    */

    //matrix<2,1> A;

    if(pdf)
    {   
        TCanvas* c = new TCanvas("c", "Scatter Plot with Fit", 800, 600);
        // Create two pads: top (main), bottom (residuals)
        TPad* pad1 = new TPad("pad1", "Main plot", 0, 0.25, 1, 1.0);
        pad1->SetBottomMargin(0.07);
        pad1->Draw();
        TPad* pad2 = new TPad("pad2", "Residuals", 0, 0.0, 1, 0.25);
        pad2->SetTopMargin(0.1);
        pad2->SetBottomMargin(0.3);
        pad2->Draw();

        // Draw main plot in pad1
        pad1->cd();
        TMultiGraph* mg = new TMultiGraph("mg",
            ("Spin tune vs change in travel time @ " + proton_RF::mode + " (" + proton_RF::ID + ");#Deltat/t;#nu_{s}").data());
        plot<double> gr(ps,nu,// - A[2]*ep_x - A[3]*ep_y,
            "Single Particle Spin Tune");
        gr->SetMarkerStyle(20);
        gr->SetLineColorAlpha(0,0);
        mg->Add(gr.get(),"P");

        stringstream ss;
        ss<<"#nu_{s} = "<<(A[0]<0 ? "- " : "")<<abs(A[0])<<(A[1]<0 ? " - " : " + ")<<abs(A[1])<<" #frac{#Deltat}{t}"<<endl;
        string ts; getline(ss,ts);
        plot<double> fit(ps,A[0] + A[1]*ps,ts);
        fit->SetLineColor(kRed);
        fit->SetMarkerColorAlpha(0,0);
        mg->Add(fit.get(),"L");

        mg->Draw("AP");

        auto leg = pad1->BuildLegend(0.12,0.25,0.5,0.11);
        leg->SetFillColorAlpha(10,0.2);
        pad1->SetGrid();

        // Draw residuals in pad2
        pad2->cd();
        valarray<double> fit_vals = A[0] + A[1]*ps;
        valarray<double> ydata = nu;// - A[2]*ep_x - A[3]*ep_y;
        valarray<double> residuals = ydata - fit_vals;
        valarray<double> zero(0.0,residuals.size());
        plot<double> gres(ps, residuals, "");
        plot<double> gzero(ps, zero, "Zero Line");
        gres->Sort();
        gres->SetMarkerStyle(7);
        gres->SetLineColor(2);
        gres->GetXaxis()->SetTitle("#Deltat/t");
        gres->GetYaxis()->SetTitle("Residuals");
        gres->GetYaxis()->SetTitleSize(0.12);
        gres->GetYaxis()->SetLabelSize(0.10);
        gres->GetYaxis()->SetTitleOffset(0.3);
        gres->GetXaxis()->SetTitleSize(0.12);
        gres->GetXaxis()->SetLabelSize(0.10);
        //gStyle->SetTitleFontSize(0.1);
        gres->Draw("APL");
        gzero->SetLineColor(1);
        gzero->SetLineWidth(2);
        gzero->Draw("L");
        pad2->SetGrid();

        c->cd();
        c->SaveAs((path + "/correlation.pdf").data(),"pdf");
    }
}
