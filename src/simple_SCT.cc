#include<proton_RF.hh>

array<double,2> proton_RF::simple_SCT(string path)
{
    Data<double> D(path);
    vector<double> sx = D[1], sy = D[2], sz = D[3];
    double N = D.order()[0];

    path = path.substr(0,path.find_last_of("/"));

    valarray<double> mag(0.0,N);
    for(int i=0; i<N; i++) mag[i] = pow(pow(sx[i],2) + pow(sy[i],2) + pow(sz[i],2),0.5);
    
    valarray<double> th(0.0,N), m(0.0,N), n(0.0,N);
    for(int i=0; i<N; i++) m[i] = sx[i] / sz[i];
    m = (m.shift(1) - m) / (1 + m * m.shift(1));
    m = atan(m);

    cout<<"using mode: "<<mode<<endl;

	m = m + 2*pi*nu0[ID][mode][0] + 2*pi*nu0[ID][mode][1] + 2*pi*nu0[ID][mode][2];

    th[0] = 0; n[0] = 0;
    for(int i=0; i<N-1; i++)
    {
        th[i+1] = th[i] + m[i];
        n[i+1] = i+1;
    }

    double E = pow(pow(p[mode],2) + pow(m0[ID],2),0.5);
    double to_t = L[mode] / (p[mode]/E * cel);

    
    matrix<2,1,valarray<double>> vA = {pow(n * 2 * pi,0),pow(n * 2 * pi,1)};
    auto A = (sum(vA * vA.tp())^-1) * sum(vA * th);
    auto B = A.tp() * vA;
 
    double nu_s = -A[1];
    
    matrix<2,1,valarray<double>> vC = {pow(n * to_t,0),pow(n * to_t,2)};
    auto C = (sum(vC * vC.tp())^-1) * sum(vC * mag);
    auto F = C.tp() * vC;
    
    double t = pow((exp(-1) - C[0])/(C[1]),0.5);

    vector<double> rfi, rth, rt, rmag, rB, rF;
    for(int i=1; i<N; i=i+(N>100000 ? N/100000 : 1))
    {
        rfi.push_back(2*pi*n[i]); rth.push_back(th[i]);
        rt.push_back(n[i]*to_t); rmag.push_back(mag[i]);
        rB.push_back(B[i]); rF.push_back(F[i]);
    }
    //valarray<double> rfi = 2*pi*n, &rth = th, rt = n*to_t, &rmag = mag, &rB = B, &rF = F;

    auto mg_sct = new TMultiGraph("SCT",
        "Magnitude of resultant spin vector;time (seconds);#left|P(x)#right| (radians)");

    plot<double> mag_data(rt,rmag);
    mag_data->SetMarkerStyle(1);
    mag_data->SetLineColor(10);
    mg_sct->Add(mag_data.get(),"P");
 
    stringstream ss;
    ss<<"#left|P(x)#right| = "<<(C[1]<0 ? "- " : "")<<abs(C[1])<<" t^2"
        <<(C[0]<0 ? " - " : " + ")<<abs(C[0])<<endl;
    string ts; getline(ss,ts);
 
    plot<double> fi_mag(rt,rF,ts);
    fi_mag->SetMarkerSize(1);
    fi_mag->SetLineColor(2);
    fi_mag->SetLineWidth(2);
    mg_sct->Add(fi_mag.get(),"L");

    auto mg_nus = new TMultiGraph("nu_s",
        "Relative angle of spin vector;#phi (radians);#theta_{x} (radians)");

    plot<double> th_data(rfi,rth,"#theta_{x} at given turn");
    th_data->SetMarkerStyle(1);
    th_data->SetLineColor(10);
    mg_nus->Add(th_data.get(),"P");

    ss.str() = "";
    ss<<"#theta_{x} = "<<(A[0]<0 ? "- " : "")<<abs(A[0])<<" #phi"
                    <<(A[1]<0 ? " - " : " + ")<<abs(A[1])<<endl;
    getline(ss,ts);

    plot<double> fi_th(rfi,rB,ts);
    fi_th->SetMarkerSize(1);
    fi_th->SetLineColor(2);
    fi_th->SetLineWidth(2);
    mg_nus->Add(fi_th.get(),"L");

    TCanvas* csct = new TCanvas("csct","t_mag");
    mg_sct->Draw("AP");
    mg_sct->GetYaxis()->SetTitleOffset(0.9);
    auto leg1 = csct->BuildLegend(0.11,0.25,0.64,0.11);
    leg1->SetFillColorAlpha(10,0.2);
    csct->SetGrid();
    csct->SaveAs((path + "/SCT.pdf").data(),"pdf");

    cout<<"Spin Coherence Time: "<<t<<" seconds"<<endl;
    
    TCanvas* cth = new TCanvas("cth","fi_th");
    mg_nus->Draw("AP");
    mg_nus->GetYaxis()->SetTitleOffset(0.9);
    auto leg2 = cth->BuildLegend(0.11,0.25,0.64,0.11);
    leg2->SetFillColorAlpha(10,0.2);
    cth->SetGrid();
    cth->SaveAs((path + "/theta_fit.pdf").data(),"pdf");

    cout<<"Spin tune: "<<nu_s<<endl;

    return {t, nu_s};
}