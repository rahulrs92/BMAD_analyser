#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"

#include<lv2.hh>

int main_c_vs_emittance(int argc, char **argv)
{
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/escan_C";

	auto ssrc = bash::to_int("find " + path + " -type d -name \"r*\" | wc -l");
	cout<<"ssrc = "<<ssrc<<endl;

	fstream f(path + "/h_info", fstream::out | fstream::trunc);
	f<<"#epsilon_{x}"<<"\n"<<"#epsilon_{y}"<<"\n"<<"#delta_{0}"<<"\n";
	f<<"#nu_{s0}"<<"\n"<<"c[0]"<<"\n"<<"c[1]";
	f.close();
	f.open(path + "/h", fstream::out | fstream::trunc);
	
	lv2 az;

	for(int i=1; i<=ssrc; i++)
	{
		string path1 = path + "/r" + to_string(i);
		cout<<"Currently working on "<<path1<<endl;
		
		auto nc = bash::to_int("find " + path1 + " -type d -name \"c*\" | wc -l");
		cout<<"nc = "<<nc<<endl;
		
		auto ln = proton::line(path1,nc,1);
		f<<ln.col("#epsilon_{x0}")[0]<< "\t" << ln.col("#epsilon_{y0}")[0] << "\t"
			 <<ln.col("#delta_{0}")[0] << "\t";
		az.tune_vs_chromaticity(path1,1);
		f<<az.nu_s0<<"\t"<<az.c[0]<<"\t"<<az.c[1]<<endl;
	}
	
	f.close();

	Data<double> D(path + "/h");
	valarray<double> dt = D.col("#delta_{0}");
	matrix<1,3,valarray<double>> eps = {D.col("#epsilon_{x}"), D.col("#epsilon_{y}"), pow(dt,2)};
	matrix<2,1,valarray<double>> xi = {D.col("#xi_{x}"), D.col("#xi_{y}")};
	matrix<3,1,valarray<double>> c = {D.col("#nu_{s0}"), D.col("c[0]"), D.col("c[1]")};
	valarray<double> dn = D.col(0);
	
	int n = eps[0].size();
	cout<<"n = "<<n<<endl;
	//cout<<"c = "<<c<<endl;
	//cout<<"eps = "<<eps<<endl;
	
	/*
	matrix<2,1,valarray<double>> sm_eps = valarray<double>(pow(n,2));
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++)
		{
			sm_eps[0][i*n + j] = eps[0][j];
			sm_eps[1][i*n + j] = eps[1][i];
		}   
	}
	//cout<<"sm_eps = "<<sm_eps<<endl;
	
	//matrix<2,1,valarray<double>> sm_xi = xi;
	*/

	matrix<2,1,valarray<double>> regs_eps = {eps[0],eps[1]};
	matrix<2,1,valarray<double>> sm_eps = regs_eps;
	auto A = sum(c * regs_eps.tp()) * (sum(regs_eps * regs_eps.tp())^-1);
	matrix<1,2> h0 = A[0];
	matrix<2,2> H = A.sub<2,2>(1,0);
	auto sm_c = A*sm_eps;

	TApplication app("app", &argc, argv);

	TCanvas *chcr = new TCanvas("canvas2", "Spin tune vs. Emittance", 2400, 600);
	chcr->Divide(3,1);

	matrix<1,2,valarray<double>> eps1 = eps.sub<1,2>();

	for(int i=0; i<3; i++)
	{
		chcr->cd(i+1);

		valarray<double> fit_vals = (A*regs_eps)[i];
        valarray<double> ydata = c[i];

		string plot_str = (i==0 ? "#nu_{s} @ zero chromaticity" : "c_{" + to_string(i) + "}");
        plot2<double> fit2(sm_eps[0],sm_eps[1],sm_c[i],plot_str + " vs emittance @ " + 
            proton_RF::mode + " (" + proton_RF::ID + ");#epsilon_{x};#epsilon_{y};c_{" + to_string(i+1) + "}");
        fit2->SetMarkerColorAlpha(0,0);
        fit2->SetLineColor(2);

        plot2<double> graph2(eps[0],eps[1],ydata,plot_str);
        graph2->SetMarkerStyle(20);
        graph2->SetMarkerColor(4);
        graph2->SetLineColor(4);

		string pad_name = "pad_c" + to_string(i) + "_1";
        TPad* pad1 = new TPad(pad_name.data(), "Main plot", 0, 0.25, 1, 1.0);
        pad1->SetBottomMargin(0.07);
        pad1->Draw();

        pad_name = "pad_c" + to_string(i) + "_2";
        TPad* pad2 = new TPad(pad_name.data(), "Residuals", 0, 0.0, 1, 0.25);
        pad2->SetTopMargin(0.1);
        pad2->SetBottomMargin(0.3);
        pad2->Draw();

        pad1->cd();
        fit2->Draw("surf1");
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

	}
	chcr->SaveAs((path + "/c_vs_emittance.pdf").data(), "pdf");

	//az.tune_vs_emittance(path + "/r10",1);

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	return 0;
}
