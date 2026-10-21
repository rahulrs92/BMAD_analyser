#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TH1D.h"

#include<lv2.hh>

int main_tracking_chromaticity(int argc, char **argv)
{
	//string path = "/Users/utente/JEDI/PTEDM/1.0/source_prototype2/build/data/sct2";
	string path = "/Users/utente/JEDI/PTEDM/1.0/source_COSY/build/data/Deuteron/sphere/r1";
	TApplication app("app", &argc, argv);

	Data<double> ln(path + "/line");
	ln = ln + Data<double>(path + "/line_error");
	ln.clean_nans();
	valarray<double> Cx = ln.col("#xi_{x}"), Cy = ln.col("#xi_{y}");
	valarray<double> Cxm = ln.col("#xi_{xm}"), Cym = ln.col("#xi_{ym}");
	valarray<double> Cxe = ln.col("#Delta#xi_{xm}"), Cye = ln.col("#Delta#xi_{ym}");

	cout<<"Chromaticity: "<<endl;
	for(size_t i=0; i<10; i++)
	{
		cout<<"Cx = "<<Cx[i]<<", Cxm = "<<Cxm[i]<<" +/- "<<Cxe[i]<<", Cy = "<<Cy[i]<<", Cym = "<<Cym[i]<<" +/- "<<Cye[i]<<endl;
	}
	valarray<double> res_Cx = Cxm - Cx, res_Cy = Cym - Cy;
	valarray<double> w_Cx = pow(Cxe,-2), w_Cy = pow(Cye,-2);
	w_Cx = w_Cx/w_Cx.sum() * w_Cx.size();
	w_Cy = w_Cy/w_Cy.sum() * w_Cy.size();
	// w_Cx = pow(w_Cx,0); w_Cy = pow(w_Cy,0); // to switch between weighted and unweighted residuals
	double std_res_Cx = pow((pow(res_Cx,2) * w_Cx).sum()/res_Cx.size(),0.5);
	double std_res_Cy = pow((pow(res_Cy,2) * w_Cy).sum()/res_Cy.size(),0.5);
	cout<<"Standard deviation of residuals: "<<std_res_Cx<<" (x), "<<std_res_Cy<<" (y)"<<endl;
	cout<<"Max of residuals: "<<res_Cx.max()<<" (x), "<<res_Cy.max()<<" (y)"<<endl;
	TH1D *hCx = new TH1D("hCx","hCx",100,-2*std_res_Cx,2*std_res_Cx);
	TH1D *hCy = new TH1D("hCy","hCy",100,-2*std_res_Cy,2*std_res_Cy);
	for(int i=0; i<Cxm.size(); i++)
	{
		hCx->Fill(res_Cx[i], w_Cx[i]);//
		hCy->Fill(res_Cy[i], w_Cy[i]);//
	}
	TCanvas* c1 = new TCanvas("c1","c1",800,600);
	hCx->SetTitle("Chromaticity: Machine vs Tracking;#xi_{xM} - #xi_{xT};Frequency");
	hCx->SetLineColor(1);
	hCx->Fit("gaus");
	hCx->Draw("HIST");
	// hCx->SetFillColor(2);
	hCx->Draw("FUNC SAME");
	TCanvas* c2 = new TCanvas("c2","c2",800,600);
	hCy->SetTitle("Chromaticity: Machine vs Tracking;#xi_{yM} - #xi_{yT};Frequency");
	hCy->SetLineColor(1);
	hCy->Fit("gaus");
	hCy->Draw("HIST");
	// hCy->SetFillColor(2);
	hCy->Draw("FUNC SAME");

	TRootCanvas *rc = (TRootCanvas *)gPad->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	
	return 0;
}