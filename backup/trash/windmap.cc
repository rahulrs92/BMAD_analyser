#include<lv1.hh>

using namespace std;

tuple<matrix<2,1>,matrix<2,3>> lv1::windmap(string path, bool pdf)
{
    int nc = bash::to_int("ls -d " + path + "/c* | wc -l");
	Data<double> D = proton::line(path,nc,pdf);
    
    int N = D.order()[0];
    
    valarray<double> tf = D.col("t_{f}");
    valarray<double> td = D.col("t_{d}");
    valarray<double> ts = D.col("t_{s}");
    valarray<double> x = D.col("#xi_{x}");
    valarray<double> y = D.col("#xi_{y}");

    matrix<4,1,valarray<double>> vA = {pow(tf,0),tf,td,ts};
    matrix<2,1,valarray<double>> vB = {x,y};

    auto A = sum(vB * vA.tp()) * (sum(vA * vA.tp())^-1);
    cout<<"A = "<<A<<endl;
    matrix<2,1> xi0 = A(0);
    cout<<"xi0 = "<<xi0<<endl;
    matrix<2,3> Xi = A.sub<2,3>(0,1);
    
    cout<<"Xi = "<<Xi<<endl;

    auto B = A * vA;

    if(pdf)
    {
        TCanvas* cpts = new TCanvas("cpts","cpts",1000,1000);
        plot<double> pts(x,y,"Chromaticities: Actual vs Linear Map;#xi_{x};#xi_{y}");
        plot<double> wmp(B[0],B[1]);
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
        cpts->SaveAs((path + "/wm_error_plot.pdf").data(),"pdf");
    }
    return make_tuple(xi0,Xi);
}