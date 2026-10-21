#include<lv1.hh>
#include <TH2F.h>

tuple<matrix<2,1>,matrix<2,3>> lv1::windmap_plane(Data<double>& D,string path,bool out)
{
    valarray<double> tf = D.col("t_{f}");
    valarray<double> td = D.col("t_{d}");
    valarray<double> ts = D.col("t_{s}");
    valarray<double> x = D.col("#xi_{x}");
    valarray<double> y = D.col("#xi_{y}");

    matrix<10,1,valarray<double>> regs_t = {pow(tf,0),tf,td,ts,tf*tf,td*td,ts*ts,tf*td,td*ts,ts*tf};
    matrix<2,1,valarray<double>> xi = {x,y};

    auto A = sum(xi * regs_t.tp()) * (sum(regs_t * regs_t.tp())^-1);
    cout<<"regs_t * regs_t.tp() = "<<det(sum(regs_t * regs_t.tp()))<<endl;
    cout<<"A = "<<A<<endl;
    xi0 = A(0);
    //cout<<"xi0 = "<<xi0<<endl;
    Xi = A.sub<2,3>(0,1);
    //cout<<"Xi = "<<Xi<<endl;

    auto sm_xi = A * regs_t;
    //cout<<"sm_xi = "<<sm_xi<<endl;
    auto res_sm_xi = xi - sm_xi;
    
    auto cor_sm_xi = sum(res_sm_xi * res_sm_xi.tp()) / sm_xi[0].size();

    if(out)
    {
        TCanvas* cpts = new TCanvas("cpts","cpts",1000,1000);
        plot<double> pts(x,y,"Chromaticities: Actual vs Linear Map;#xi_{x};#xi_{y}");
        plot<double> wmp(sm_xi[0],sm_xi[1]);
        pts->Draw("AP");
        pts->SetMarkerStyle(3);
        pts->SetLineColor(0);
        wmp->Draw("P same");
        wmp->SetMarkerColor(2);
        wmp->SetMarkerStyle(2);
        wmp->SetLineColor(0);
        auto leg = new TLegend(0.3,0.9,0.7,0.8);
        leg->AddEntry(pts.get(),"Data Point (Measured Values)");
        leg->AddEntry(wmp.get(),"Linear Map (Calculated Values)");
        leg->SetFillColorAlpha(10,0.2);
        leg->Draw();
        cpts->SetGrid();
        cpts->SaveAs((path + "/wm_error_plot2.pdf").data(),"pdf");

        TCanvas* ceh = new TCanvas("ceh","ceh",1024,768);

        TH2F* hnu = new TH2F("hnu","Distribution of residuals of Linear Model fit;#xi_{x};#xi_{y};Frequency",
                                    101,-5*pow(cor_sm_xi[0][0],0.5),5*pow(cor_sm_xi[0][0],0.5),
                                    101,-5*pow(cor_sm_xi[1][1],0.5),5*pow(cor_sm_xi[1][1],0.5));
        for(int i=0; i<res_sm_xi[0].size(); i++) hnu->Fill(res_sm_xi[0][i],res_sm_xi[1][i]);
        //hnu->Fit("gaus");
        hnu->Draw("colz");
        ceh->SaveAs((path + "/wm_error_dist2.pdf").data(),"pdf");
    }
    return make_tuple(xi0,Xi);
}