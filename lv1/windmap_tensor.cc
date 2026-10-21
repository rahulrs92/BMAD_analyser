#include<lv1.hh>
#include <TH2F.h>

tuple<matrix<3,1>,matrix<3,3>> lv1::windmap_tensor(Data<double>& D,string path,bool out)
{
    D.clean_nans();

    valarray<double> tf = D.col("t_{f}");
    valarray<double> td = D.col("t_{d}");
    valarray<double> ts = D.col("t_{s}");
    valarray<double> x = D.col("#xi_{x}");
    valarray<double> y = D.col("#xi_{y}");
    valarray<double> z = D.col("#eta_{1}");
    valarray<double> z_e = D.col("#Delta#eta_{1}");

    valarray<double> w = pow(z_e,-2);

    matrix<10,1,valarray<double>> regs_t = {pow(tf,0),tf,td,ts,tf*tf,td*td,ts*ts,tf*td,td*ts,ts*tf};
    //matrix<4,1,valarray<double>> regs_t = {pow(tf,0),tf,td,ts};
    matrix<3,1,valarray<double>> xi = {x,y,z}, t = {tf,td,ts};

    fstream f(path + "/wm_matrix_info", fstream::out | fstream::trunc);
    f<<"xi0[0]"<<"\n"<<"xi0[1]"<<"\n"<<"xi0[2]"<<"\n";
	f<<"Xi[0][0]"<<"\n"<<"Xi[0][1]"<<"\n"<<"Xi[0][2]"<<"\n";
	f<<"Xi[1][0]"<<"\n"<<"Xi[1][1]"<<"\n"<<"Xi[1][2]"<<"\n";
	f<<"Xi[2][0]"<<"\n"<<"Xi[2][1]"<<"\n"<<"Xi[2][2]"<<"\n";
	f.close();

    f.open(path + "/wm_matrix_error_info", fstream::out | fstream::trunc);
    f<<"e_xi0[0]"<<"\n"<<"e_xi0[1]"<<"\n"<<"e_xi0[2]"<<"\n";
	f<<"e_Xi[0][0]"<<"\n"<<"e_Xi[0][1]"<<"\n"<<"e_Xi[0][2]"<<"\n";
	f<<"e_Xi[1][0]"<<"\n"<<"e_Xi[1][1]"<<"\n"<<"e_Xi[1][2]"<<"\n";
	f<<"e_Xi[2][0]"<<"\n"<<"e_Xi[2][1]"<<"\n"<<"e_Xi[2][2]"<<"\n";
	f.close();

    f.open(path + "/wm_tensor_info", fstream::out | fstream::trunc);
    f<<"xi0"<<"\n"<<"Xi(0)"<<"\n""Xi(1)"<<"\n""Xi(2)"<<"\n";
	f<<"Xi2[0][0]"<<"\n"<<"Xi2[1][1]"<<"\n"<<"Xi2[2][2]"<<"\n";
	f<<"Xi2[0][1]"<<"\n"<<"Xi2[1][2]"<<"\n"<<"Xi2[2][0]"<<"\n";
	f.close();

    f.open(path + "/wm_tensor_error_info", fstream::out | fstream::trunc);
    f<<"e_xi0"<<"\n"<<"e_Xi(0)"<<"\n""e_Xi(1)"<<"\n""e_Xi(2)"<<"\n";
	f<<"e_Xi2[0][0]"<<"\n"<<"e_Xi2[1][1]"<<"\n"<<"e_Xi2[2][2]"<<"\n";
	f<<"e_Xi2[0][1]"<<"\n"<<"e_Xi2[1][2]"<<"\n"<<"e_Xi2[2][0]"<<"\n";
	f.close();

    auto A = sum(regs_t * regs_t.tp())^-1;
    auto B = sum(regs_t * w * regs_t.tp())^-1;
    matrix<3,regs_t.M()> C, e_C;
    C[0] = sum(xi[0] * regs_t.tp()) * A;
    C[1] = sum(xi[1] * regs_t.tp()) * A;
    C[2] = sum((xi * w)[2] * regs_t.tp()) * B;
    cout<<"C = "<<C<<endl;
    xi03 = C(0);
    //cout<<"xi0 = "<<xi0<<endl;
    Xi3 = C.sub<3,3>(0,1);
    //cout<<"Xi = "<<Xi<<endl;

    cout<<"test value @ "<<pull(regs_t,0).sub<3,1>(1,0).tp()<<"\n"<<A * pull(regs_t,0)<<endl;
    cout<<"Linear inverse solution of "<<pull(xi,0).tp()<<":\n"<<(Xi3^-1)*(pull(xi,0)-xi03)<<endl;

    f.open(path + "/wm_matrix", fstream::out | fstream::trunc);
    f<<xi03[0]<<"\t"<<xi03[1]<<"\t"<<xi03[2]<<"\t"
      <<Xi3[0][0]<<"\t"<<Xi3[0][1]<<"\t"<<Xi3[0][2]<<"\t"
      <<Xi3[1][0]<<"\t"<<Xi3[1][1]<<"\t"<<Xi3[1][2]<<"\t"
      <<Xi3[2][0]<<"\t"<<Xi3[2][1]<<"\t"<<Xi3[2][2]<<endl;
    f.close();

    f.open(path + "/wm_tensor", fstream::out | fstream::trunc);
    f<<C<<endl;
    f.close();

    auto sm_xi = C * regs_t;// + B * regs_t2;
    auto res_sm_xi = xi - sm_xi;

    e_C[0] = (A * sum((regs_t * res_sm_xi[0] * res_sm_xi[0] * regs_t.tp())) * A).diag();
    e_C[1] = (A * sum((regs_t * res_sm_xi[1] * res_sm_xi[1] * regs_t.tp())) * A).diag();
    e_C[2] = (B * sum((regs_t * w * res_sm_xi[2] * res_sm_xi[2] * w * regs_t.tp())) * B).diag();
    xi03_e = e_C(0);
    Xi3_e = e_C.sub<3,3>(0,1);

    f.open(path + "/wm_matrix_error", fstream::out | fstream::trunc);
    f<<pow(xi03_e[0],0.5)<<"\t"<<pow(xi03_e[1],0.5)<<"\t"<<pow(xi03_e[2],0.5)<<"\t"
      <<pow(Xi3_e[0][0],0.5)<<"\t"<<pow(Xi3_e[0][1],0.5)<<"\t"<<pow(Xi3_e[0][2],0.5)<<"\t"
      <<pow(Xi3_e[1][0],0.5)<<"\t"<<pow(Xi3_e[1][1],0.5)<<"\t"<<pow(Xi3_e[1][2],0.5)<<"\t"
      <<pow(Xi3_e[2][0],0.5)<<"\t"<<pow(Xi3_e[2][1],0.5)<<"\t"<<pow(Xi3_e[2][2],0.5)<<endl;
    f.close();
    
    f.open(path + "/wm_tensor_error", fstream::out | fstream::trunc);
    f<<e_C<<endl;
    f.close();

    auto cor_sm_xi = sum(res_sm_xi * res_sm_xi.tp()) / sm_xi[0].size();

    cout<<"Coefficient variances of Wind Map fit: "<<endl;
    cout<<e_C<<endl;

    if(out)
    {
        TCanvas* cpts = new TCanvas("cpts","cpts",1000,1000);
        plot2<double> pts(x,y,z,"Chromaticities: Actual vs Linear Map;#xi_{x};#xi_{y};#eta_{1}");
        plot2<double> wmp(sm_xi[0],sm_xi[1],sm_xi[2]);
        pts->Draw("P0");
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
        cpts->SaveAs((path + "/wm_error_plot.pdf").data(),"pdf");

        TCanvas* ceh = new TCanvas("ceh","ceh",1024,768);

        TH2F* hxy = new TH2F("hxy","Distribution of residuals of Linear Model fit;#xi_{x};#xi_{y};Frequency",
                                    101,-5*pow(cor_sm_xi[0][0],0.5),5*pow(cor_sm_xi[0][0],0.5),
                                    101,-5*pow(cor_sm_xi[1][1],0.5),5*pow(cor_sm_xi[1][1],0.5));
        for(int i=0; i<res_sm_xi[0].size(); i++) hxy->Fill(res_sm_xi[0][i],res_sm_xi[1][i]);
        //hxy->Fit("gaus");
        hxy->Draw("colz");
        ceh->SaveAs((path + "/wm_error_dist.pdf(").data(),"pdf");

        TH1F* hx = new TH1F("hx",
            "Distribution of residuals of Linear Model fit;#xi_{x} residual;Frequency",
            51,-5*pow(cor_sm_xi[0][0],0.5),5*pow(cor_sm_xi[0][0],0.5));
        for(int i=0; i<res_sm_xi[0].size(); i++) hx->Fill(res_sm_xi[0][i]);
        hx->Fit("gaus");
        hx->Draw();
        ceh->SaveAs((path + "/wm_error_dist.pdf").data(),"pdf");

        TH1F* hy = new TH1F("hy",
            "Distribution of residuals of Linear Model fit;#xi_{y} residual;Frequency",
            51,-5*pow(cor_sm_xi[1][1],0.5),5*pow(cor_sm_xi[1][1],0.5));
        for(int i=0; i<res_sm_xi[1].size(); i++) hy->Fill(res_sm_xi[1][i]);
        hy->Fit("gaus");
        hy->Draw();
        ceh->SaveAs((path + "/wm_error_dist.pdf").data(),"pdf");

        TH1F* hz = new TH1F("hz",
            "Distribution of residuals of Linear Model fit;#eta_{1} residual;Frequency",
            51,-5*pow(cor_sm_xi[2][2],0.5),5*pow(cor_sm_xi[2][2],0.5));
        for(int i=0; i<res_sm_xi[2].size(); i++) hz->Fill(res_sm_xi[2][i]);
        hz->Fit("gaus");
        hz->Draw();
        ceh->SaveAs((path + "/wm_error_dist.pdf)").data(),"pdf");
    }
    return make_tuple(xi03,Xi3);
}