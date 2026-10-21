#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"

#include<lv2.hh>

int main_w_vs_chromaticity(int argc, char **argv)
{
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/escan_C_RF";

	auto ssrc = bash::to_int("find " + path + " -type d -name \"r*\" | wc -l");
	cout<<"ssrc = "<<ssrc<<endl;

	fstream f(path + "/h_info", fstream::out | fstream::trunc);
	f<<"#xi_{x}"<<"\n"<<"#xi_{y}"<<"\n";
	f<<"t_{f}"<<"\n"<<"t_{d}"<<"\n"<<"t_{s}"<<"\n";
	f<<"w[0]"<<"\n"<<"w[1]"<<"\n"<<"w[2]";
	f.close();
	f.open(path + "/h", fstream::out | fstream::trunc);
	
	lv2 az;

	for(int i=1; i<=ssrc; i++)
	{
		string path1 = path + "/r" + to_string(i);
		cout<<"Currently working on "<<path1<<endl;
		
		auto nc = bash::to_int("find " + path1 + " -type d -name \"c*\" | wc -l");
		cout<<"nc = "<<nc<<endl;
		
		auto ln = proton_RF::line(path1,nc,1);
		f<<ln.col("#xi_{x}")[0]<< "\t" << ln.col("#xi_{y}")[0] << "\t";
		f<<ln.col("t_{f}")[0]<< "\t" << ln.col("t_{d}")[0] << "\t" << ln.col("t_{s}")[0] << "\t";
		az.tune_vs_emittance_RF(path1,1);
		f<<az.w[0]<<"\t"<<az.w[1]<<"\t"<<az.w[2]<<endl;
	}
	
	f.close();

	Data<double> D(path + "/h");
	valarray<double> dt = D.col("#delta_{0}");
	matrix<2,1,valarray<double>> xi = {D.col("#xi_{x}"), D.col("#xi_{y}")};
	matrix<3,1,valarray<double>> t = {D.col("t_{f}"), D.col("t_{d}"), D.col("t_{s}")};
	matrix<2,1,valarray<double>> w = {D.col("w[0]"), D.col("w[1]")};
	valarray<double> dn = D.col(0);
	
	int n = xi[0].size();
	cout<<"n = "<<n<<endl;
	//cout<<"c = "<<c<<endl;
	//cout<<"eps = "<<eps<<endl;


	valarray<double> sm_xi_x = valarray<double>(pow(n,2)), sm_xi_y = valarray<double>(pow(n,2));
	/*
	double d_x = (t[0].max() - t[0].min())/(n);
	double d_y = (t[1].max() - t[1].min())/(n);
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++)
		{
			sm_xi_x[i*n + j] = t[0].min() + i * d_x;
			sm_xi_y[i*n + j] = t[1].min() + j * d_y;
		}   
	}
	*/
	double d_x = (xi[0].max() - xi[0].min())/(n);
	double d_y = (xi[1].max() - xi[1].min())/(n);
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++)
		{
			sm_xi_x[i*n + j] = xi[0].min() + i * d_x;
			sm_xi_y[i*n + j] = xi[1].min() + j * d_y;
		}   
	}
	
	matrix<2,1> xi0 = {-4.24664, -3.93523};
	matrix<2,1,valarray<double>> sm_xi = {sm_xi_x, sm_xi_y};
	matrix<2,2> Xi = {2.62575, 2.91311, -1.86668, -13.2338};

	//sm_xi = xi0 + Xi * sm_xi;

	auto sm_t = (Xi^-1) * (sm_xi - xi0);

	
	//cout<<"sm_eps = "<<sm_eps<<endl;

	matrix<3,1,valarray<double>> regs_xi = {pow(xi[0],0),xi[0],xi[1]};
	matrix<3,1,valarray<double>> sm_xi_r = {pow(sm_xi[0],0),sm_xi[0],sm_xi[1]};
	auto A = sum(w * regs_xi.tp()) * (sum(regs_xi * regs_xi.tp())^-1);
	//matrix<3,1,valarray<double>> sm_h = A*sm_xi_r;
	matrix<2,1,valarray<double>> sm_w = A*sm_xi_r;
	cout<<"A = "<<A<<endl;

	TApplication app("app", &argc, argv);

	TCanvas *chcr = new TCanvas("canvas2", "h vs. chromaticity", 1600, 600);
	chcr->Divide(2,1);

	//matrix<1,2,valarray<double>> eps1 = eps.sub<1,2>();

	for(int i=0; i<2; i++)
	{
		chcr->cd(i+1);

		valarray<double> fit_vals = (A*regs_xi)[i];
        //valarray<double> fit_vals = A[i][0] + A[i][1]*xi[0] + A[i][2]*xi[1];
        valarray<double> ydata = w[i];
		cout<<"fit_vals = "<<fit_vals<<endl;

		string plot_str = i<2 ? "w_{" + to_string(i+1) + "}^{*}" : "#sigma_{1}";
        plot2<double> fit2(sm_xi[0],sm_xi[1],sm_w[i],plot_str + " vs chromaticity @ " + 
            proton_RF::mode + " (" + proton_RF::ID + ");#xi_{x};#xi_{y};" + plot_str);
        fit2->SetMarkerColorAlpha(0,0);
        fit2->SetLineColor(2);

        plot2<double> graph2(xi[0],xi[1],ydata,plot_str);
        graph2->SetMarkerStyle(20);
        graph2->SetMarkerColor(4);
        graph2->SetLineColor(4);

		string pad_name = "pad_w" + to_string(i) + "_1";
        TPad* pad1 = new TPad(pad_name.data(), "Main plot", 0, 0.25, 1, 1.0);
        pad1->SetBottomMargin(0.07);
        pad1->Draw();

        pad_name = "pad_w" + to_string(i) + "_2";
        TPad* pad2 = new TPad(pad_name.data(), "Residuals", 0, 0.0, 1, 0.25);
        pad2->SetTopMargin(0.1);
        pad2->SetBottomMargin(0.3);
        pad2->Draw();

        pad1->cd();
        fit2->Draw("surf1");
        //graph2->Draw("AP0");
        graph2->Draw("P0 same");
        fit2->GetXaxis()->SetTitleOffset(2);
        fit2->GetYaxis()->SetTitleOffset(2);

        pad2->cd();

        valarray<double> residuals = ydata - fit_vals;
        valarray<double> zero(0.0,residuals.size());
        plot<double> gres(dn+1, residuals, "");
        plot<double> gzero(dn+1, zero, "Zero Line");
        gres->Sort();
        gres->SetMarkerStyle(7);
        gres->SetLineColor(2);
        gres->GetXaxis()->SetTitle("Data Point Index");
        gres->GetYaxis()->SetTitle("Residuals");
        gres->GetYaxis()->SetTitleSize(0.1);
        gres->GetYaxis()->SetLabelSize(0.10);
        gres->GetYaxis()->SetTitleOffset(0.48);
        gres->GetXaxis()->SetTitleSize(0.12);
        gres->GetXaxis()->SetLabelSize(0.10);
        gres->GetXaxis()->SetLimits(0 - 0.01*dn.size(), dn.size() + 0.01*dn.size());
        //gStyle->SetTitleFontSize(0.1);
        gres->Draw("APL");
        gzero->SetLineColor(1);
        gzero->SetLineWidth(2);
        gzero->Draw("L");
        pad2->SetGrid();

	}
	chcr->SaveAs((path + "/w_vs_chromaticity.pdf").data(), "pdf");

	//az.tune_vs_emittance(path + "/r10",1);

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	return 0;
}
