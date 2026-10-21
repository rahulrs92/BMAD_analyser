#include<lv2.hh>
#include<algorithm>

array<double,3> lv2::optimize(Data<double>& D,Data<double>& M,string path,bool out)
{
	fstream f(path + "/optimization",fstream::out | fstream::trunc);

	TCanvas* copt1 = new TCanvas("copt1","copt1");
	TMultiGraph* mg = new TMultiGraph("mg",
				("Spin Tune vs Third Sextupole Field Strength @ " + 
				proton_RF::mode + " (" + proton_RF::ID + ");t_{s};#nu_{s}").data());
	plot<double> pl(D.col("t_{s}"),D.col("#nu_{s}"),"Spin Tune at Partially Optimized Settings");
	pl->SetMarkerStyle(20);
	pl->SetLineColor(10);
	pl->Sort();
	mg->Add(pl.get(),"P");

	valarray<double> ts = valarray<double>(pl->GetX(),pl->GetN());
	valarray<double> nu = valarray<double>(pl->GetY(),pl->GetN());

	matrix<1,9,valarray<double>> vM;
	for(int i=0; i<9; i++) vM[i] = M.col(i+1);

	matrix<1,2,valarray<double>> regs_ts = {pow(ts,0),pow(ts,1)};
	auto sl = (sum(regs_ts.tp() * regs_ts)^-1) * sum(regs_ts.tp() * nu);
	auto m = (sum(regs_ts.tp() * regs_ts)^-1) * sum(regs_ts.tp() * vM);
	
	stringstream ss;
	ss<<"#nu_{s} = "<<(sl[1]<0 ? "- " : "")<<abs(sl[1])<<" t_{s}"
		<<(sl[0]<0 ? " - " : " + ")<<abs(sl[0])<<endl;
	string tts; getline(ss,tts);
	plot<double> plsl(ts,(regs_ts * sl),tts);
	//plpb->Sort();
	plsl->SetLineColor(2);
	plsl->SetLineWidth(2);
	plsl->SetMarkerColorAlpha(0,0);
	mg->Add(plsl.get(),"L");

	double x0 = - sl[0] / sl[1];
	cout<<"line fit: "<<regs_ts * sl<<endl;
	cout<<"X-intercept at "<<x0<<endl;
	f<<"line fit: "<<regs_ts * sl<<endl;
	f<<"X-intercept at "<<x0<<endl;
	
	auto opt_m = matrix<1,2>{1,x0} * m;
	
	cout<<"opt_m: "<<opt_m<<endl;
	
	matrix<3,1> opt_w;
	matrix<3,2> opt_W;
	for(int i=0; i<3; i++) opt_w[i] = opt_m[i];
	for(int i=0; i<6; i++) opt_W[i/2][i%2] = opt_m[i+3];
	auto opt_xi = ((opt_W.tp() * opt_W)^-1) * opt_W.tp() * (-opt_w);
	
	cout<<"opt_xi = "<<opt_xi<<endl;
	
	mg->Draw("AP");
	auto leg = copt1->BuildLegend(0.13,0.25,0.55,0.13);
	leg->SetFillColorAlpha(10,0.8);
	//leg->AddEntry((TObject*)0,"@ Qx=3.6322, Qy=3.7182","");
	copt1->SetGrid();
	copt1->SaveAs((path + "/optimization.pdf(").data(),"pdf");
	
	TCanvas* copt2 = new TCanvas("copt2","copt2");
	plot<double> pl2(D.col("t_{s}"),D.col("area"),
	"Optimization for Qx=3.6322, Qy=3.7182;t_{s};Area of Triangle");
	pl2->Sort();
	pl2->Draw("AP");
	
	valarray<double> area(pl2->GetY(),pl2->GetN());
	valarray<double> sgnu = nu / abs(nu);
	
	matrix<1,3,valarray<double>> regs_ts2 = {pow(ts,0),pow(ts,1),pow(ts,2)};
	
	/*
	matrix<1,2,valarray<double>> regs_ts = {pow(ts,0),pow(ts,1)};
	valarray<valarray<double>> yy(sgnu * pow(area,0.5),area.size());
	for(int i=0; i<yy.size(); i++) { yy[i] *= sgnu[i] * pow(area[i],0.5); }
	cout<<"yy = "<<endl;
	for(auto a : yy) cout<<a<<endl;
	cout<<endl;
	matrix<2,1,valarray<double>> left = {area,area};
	for(int i=0; i<area.size(); i++)
	{
		auto pb1 = (sum(regs_ts.tp() * regs_ts)^-1) * sum(regs_ts.tp() * yy[i]);
		left[0][i] = pb1[0]; left[1][i] = pb1[1];
	}
	auto pb = sum(left * regs_ts) * (sum(regs_ts.tp() * regs_ts)^-1);
	cout<<"pb = "<<pb<<endl;
	cout<<"det(pb) = "<<det(pb)<<endl;
	cout<<"Vertex at "<< - pb[0][1] / pb[1][1]<<endl;
	*/

	auto pb = (sum(regs_ts2.tp() * regs_ts2)^-1) * sum(regs_ts2.tp() * area);

	plot<double> plpb(ts,(regs_ts2 * pb));
	//plpb->Sort();
	plpb->SetLineColor(2);
	plpb->SetLineWidth(2);
	plpb->Draw("L same");

	cout<<"parabola fit: "<<regs_ts2 * pb<<endl;
	cout<<"Vertex at "<< - pb[1] / (2*pb[2]) <<endl;
	f<<"parabola fit: "<<regs_ts2 * pb<<endl;
	f<<"Vertex at "<< - pb[1] / (2*pb[2]) <<endl;
	f.close();

	copt2->SetGrid();
	copt2->SaveAs((path + "/optimization.pdf)").data(),"pdf");

	cout<<"t_{s}: "<<ts<<endl;
	cout<<"areas: "<<area<<endl;

	return {opt_xi[0],opt_xi[1],x0};
}
