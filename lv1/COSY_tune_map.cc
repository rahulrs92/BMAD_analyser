#include<lv1.hh>
#include <TH2F.h>

matrix<3,3> lv1::COSY_tune_map(Data<double>& D,matrix<3,1> Q0,matrix<3,1> k0,string path,bool out)
{
    D.clean_nans();
    size_t constexpr ne = Q0.M();

    valarray<double> k1 = D.col(1);
    valarray<double> k2 = D.col(2);
    valarray<double> k3 = D.col(3);

    matrix<ne,1,valarray<double>> y;
    
    // valarray<double> eta = etas.tp() * etas / etas.M();
    
    // matrix<3,1,valarray<double>> y;
    for (size_t i=0; i<ne; i++) y[i] = D.col(i+3+1);
    y = y - Q0;
    
    k1 = k1 - k0[0];
    k2 = k2 - k0[1];
    k3 = k3 - k0[2];
    matrix<9,1,valarray<double>> regs_k = {k1,k2,k3,pow(k1,2),pow(k2,2),pow(k3,2),k1*k2,k2*k3,k3*k1};

    // for(size_t i=0; i<regs_k.M(); i++) regs_k[i] /= k0[i];
    
    fstream f(path + "/tunes_matrix_info", fstream::out | fstream::trunc);
    f<<"Q0"<<"\n"<<"A(0)"<<"\n"<<"A(1)"<<"\n"<<"A(2)"<<"\n";
	f.close();
    
    auto A = sum(regs_k * regs_k.tp())^-1;
    matrix<y.M(),regs_k.M()> C;
    C = sum(y * regs_k.tp()) * A;
    cout<<"C = "<<C<<endl;
    
    // matrix<C.M(),C.N()> U;
    // valarray<double> s;
    // matrix<C.N(),C.N()> V;
    // C.SVD(U,s,V);
    // cout<<"Singular values = "<<s<<endl;

    matrix<ne,1> Q_out;
    Q_out[0] = 3.62;
    Q_out[1] = 3.585;
    Q_out[2] = 0.0;

    matrix<y.M(),regs_k.M()+1> C_out;
    C_out.sub<y.M(),regs_k.M()>(0,1) = C;
    auto k_out = k0 + solve_tensor_system(C_out, Q_out - Q0);
    
    // auto k_out = k0 + ((C.tp() * C)^-1) * C.tp() * (Q_out - Q0);
    // for(size_t i=0; i<regs_k.M(); i++) k_out[i] = k_out[i] * k0[i] + k0[i];

    cout<<"k_out = \n"<<k_out<<endl;
    
    f.open(path + "/tune_sensitivity_matrix", fstream::out | fstream::trunc);
    f<<C<<endl;
    f.close();

    auto sm_y = C * regs_k;// + B * regs_t2;
    auto res_sm_y = y - sm_y;

    // sm_y = sm_y + Q0;
    // y = y + Q0;

    auto cor_sm_y = sum(res_sm_y * res_sm_y.tp()) / sm_y[0].size();

    if(out)
    {
        TCanvas* cpts = new TCanvas("cpts","cpts",1000,1000);

        plot<double> pts(y[0],y[1],"Betatron Tunes: Actual vs Local Linear Map;#DeltaQ_{x};#DeltaQ_{y}");
        plot<double> wmp(sm_y[0],sm_y[1]);
        pts->Draw("AP");
        pts->SetMarkerStyle(3);
        pts->SetLineColor(0);
        wmp->Draw("P same");
        wmp->SetMarkerColor(2);
        wmp->SetMarkerStyle(2);
        wmp->SetLineColor(0);
        auto leg = new TLegend(0.3,0.9,0.7,0.8);
        leg->AddEntry(pts.get(),"Data Point (Simulated Values)");
        leg->AddEntry(wmp.get(),"Linear Map (Calculated Values)");
        leg->SetFillColorAlpha(10,0.2);
        leg->Draw();
        cpts->SetGrid();
        cpts->SaveAs((path + "/wm_error_plot.pdf(").data(),"pdf");

        plot<double> pts2(y[1],y[2],"Betatron Tunes: Actual vs Local Linear Map;#DeltaQ_{y};#DeltaD_{BPMx1}");
        plot<double> wmp2(sm_y[1],sm_y[2]);
        pts2->Draw("AP");
        pts2->SetMarkerStyle(3);
        pts2->SetLineColor(0);
        wmp2->Draw("P same");
        wmp2->SetMarkerColor(2);
        wmp2->SetMarkerStyle(2);
        wmp2->SetLineColor(0);
        auto leg2 = new TLegend(0.3,0.9,0.7,0.8);
        leg2->AddEntry(pts2.get(),"Data Point (Simulated Values)");
        leg2->AddEntry(wmp2.get(),"Linear Map (Calculated Values)");
        leg2->SetFillColorAlpha(10,0.2);
        leg2->Draw();
        cpts->SetGrid();
        cpts->SaveAs((path + "/wm_error_plot.pdf").data(),"pdf");

        plot<double> pts3(y[0],y[2],"Betatron Tunes: Actual vs Local Linear Map;#DeltaQ_{x};#DeltaD_{BPMx1}");
        plot<double> wmp3(sm_y[0],sm_y[2]);
        pts3->Draw("AP");
        pts3->SetMarkerStyle(3);
        pts3->SetLineColor(0);
        wmp3->Draw("P same");
        wmp3->SetMarkerColor(2);
        wmp3->SetMarkerStyle(2);
        wmp3->SetLineColor(0);
        auto leg3 = new TLegend(0.3,0.9,0.7,0.8);
        leg3->AddEntry(pts3.get(),"Data Point (Simulated Values)");
        leg3->AddEntry(wmp3.get(),"Linear Map (Calculated Values)");
        leg3->SetFillColorAlpha(10,0.2);
        leg3->Draw();
        cpts->SetGrid();
        cpts->SaveAs((path + "/wm_error_plot.pdf)").data(),"pdf");

        TCanvas* ceh = new TCanvas("ceh","ceh",1024,768);

        TH2F* hxy = new TH2F("hxy","Distribution of residuals of Local Linear Model fit;#DeltaQ_{x} residual;#DeltaQ_{y} residual;Frequency",
                                    101,-5*pow(cor_sm_y[0][0],0.5),5*pow(cor_sm_y[0][0],0.5),
                                    101,-5*pow(cor_sm_y[1][1],0.5),5*pow(cor_sm_y[1][1],0.5));
        for(int i=0; i<res_sm_y[0].size(); i++) hxy->Fill(res_sm_y[0][i],res_sm_y[1][i]);
        //hxy->Fit("gaus");
        hxy->Draw("colz");
        ceh->SaveAs((path + "/wm_error_dist.pdf(").data(),"pdf");

        TH1F* hx = new TH1F("hx",
            "Distribution of residuals of Local Linear Model fit;#DeltaQ_{x} residual;Frequency",
            51,-5*pow(cor_sm_y[0][0],0.5),5*pow(cor_sm_y[0][0],0.5));
        for(int i=0; i<res_sm_y[0].size(); i++) hx->Fill(res_sm_y[0][i]);
        hx->Fit("gaus");
        hx->Draw();
        ceh->SaveAs((path + "/wm_error_dist.pdf").data(),"pdf");

        TH1F* hy = new TH1F("hy",
            "Distribution of residuals of Local Linear Model fit;#DeltaQ_{y} residual;Frequency",
            51,-5*pow(cor_sm_y[1][1],0.5),5*pow(cor_sm_y[1][1],0.5));
        for(int i=0; i<res_sm_y[1].size(); i++) hy->Fill(res_sm_y[1][i]);
        hy->Fit("gaus");
        hy->Draw();
        ceh->SaveAs((path + "/wm_error_dist.pdf").data(),"pdf");

        TH1F* hD = new TH1F("hD",
            "Distribution of residuals of Local Linear Model fit;#DeltaD_{BPMx1};Frequency",
            51,-5*pow(cor_sm_y[2][2],0.5),5*pow(cor_sm_y[2][2],0.5));
        for(int i=0; i<res_sm_y[2].size(); i++) hD->Fill(res_sm_y[2][i]);
        hD->Fit("gaus");
        hD->Draw();
        ceh->SaveAs((path + "/wm_error_dist.pdf)").data(),"pdf");

    }
    return C.sub<3,3>(0,0);
}