#include<proton.hh>

Data<double> proton_RF::line(string path, int N, bool pdf)
{
	fstream f;
	f.open(path + "/line",fstream::in);
	if(f.is_open())
	{
		f.close();
		Data<double> oD(path + "/line");
		return oD;
	}
	f.close();
	f.open((path + "/line").data(),fstream::out | fstream::trunc);
	for(int i=0; i<N; i++)
	{
		proton_RF P(path + "/c" + to_string(i+1));
		auto CSx = P.Courant_Snyder(0,pdf);
		auto CSy = P.Courant_Snyder(1,pdf);
		auto syn = P.longitudinal(pdf);
		auto eta = P.phase_slip(pdf);
		//cout<<"T = "<<eta[0]<<endl;

		double nu_s = P.spin_tune(pdf);

		f<<P.Cx()<<'\t'<<P.Cy()<<'\t'<<eta[2]<<'\t'
			<<CSx[0]<<'\t'<<CSy[0]<<'\t'<<syn[3]<<'\t'<<syn[4]<<'\t'
			<<eta[0]<<'\t'<<eta[1]<<'\t'<<nu_s<<endl;
	}
	f.close();
	Data<double> oD(path + "/line");
	oD.label(1) = "#xi_{x}";
	oD.label(2) = "#xi_{y}";
	oD.label(3) = "#eta_{1}";
	oD.label(4) = "#epsilon_{x}";
	oD.label(5) = "#epsilon_{y}";
	oD.label(6) = "#delta_{a}";
	oD.label(7) = "#omega_{s}";
	oD.label(8) = "T_{model}";
	oD.label(9) = "#eta_{0}";
	oD.label(10) = "#nu_{s}";

	f.open((path + "/line_info").data(),fstream::out | fstream::trunc);
	for(int i=0; i<oD.order()[1]; i++) f<<oD.label(i+1)<<endl;
	f.close();

	return oD;
}

Data<double> proton::line(string path, int N, bool pdf)
{
	fstream f;
	f.open(path + "/line",fstream::in);
	if(f.is_open())
	{
		f.close();
		Data<double> oD(path + "/line");
		return oD;
	}
	f.close();
	f.open((path + "/line").data(),fstream::out | fstream::trunc);
	for(int i=0; i<N; i++)
	{
		proton P(path + "/c" + to_string(i+1));
		auto CSx = P.Courant_Snyder(0,pdf);
		auto CSy = P.Courant_Snyder(1,pdf);

		auto eta = P.phase_slip(pdf);
		//cout<<"T = "<<eta[0]<<endl;
		//cout<<"Phase slip per turn = "<<eta[1]+eta[0]<<endl;

		double nu_s = P.spin_tune(pdf);

		f<<P.Cx()<<'\t'<<P.Cy()<<'\t'
			<<CSx[0]<<'\t'<<CSy[0]<<'\t'<<P.pz[0]<<'\t'
			<<eta[0]<<'\t'<<eta[1]+eta[0]<<'\t'<<nu_s<<endl;
	}
	f.close();
	Data<double> oD(path + "/line");
	oD.label(1) = "#xi_{x}";
	oD.label(2) = "#xi_{y}";
	oD.label(3) = "#epsilon_{x}";
	oD.label(4) = "#epsilon_{y}";
	oD.label(5) = "#delta_{0}";
	oD.label(6) = "T_{model}";
	oD.label(7) = "T";
	oD.label(8) = "#nu_{s}";

	f.open((path + "/line_info").data(),fstream::out | fstream::trunc);
	for(int i=0; i<oD.order()[1]; i++) f<<oD.label(i+1)<<endl;
	f.close();

	return oD;
}

