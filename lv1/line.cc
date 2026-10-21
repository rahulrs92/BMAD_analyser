#include<lv1.hh>
#include<set>
#include<chrono>

void lv1::line(string path,uint q0,bool out)
{
    fstream f, fe;
	f.open(path + "/line",fstream::in);
	if(f.is_open() && q0==0)
	{
		f.close();
		return;
	}
	f.close();

	f.open(path + "/line_info",fstream::out | fstream::trunc);
	f<<"t_{f}\nt_{d}\nt_{s}\n#xi_{x}\n#xi_{y}\n";
    f<<"#xi_{xm}\n#xi_{ym}\n#eta_{1}\n";
	f<<"#epsilon_{x}\n#epsilon_{y}\n";
	f<<"#delta_{a}\n#omega_{s}\nT\n#nu_{s}"<<endl;
	f.close();

	f.open(path + "/line_error_info",fstream::out | fstream::trunc);
	f<<"#Deltat_{f}\n#Deltat_{d}\n#Deltat_{s}\n#Delta#xi_{x}\n#Delta#xi_{y}\n";
    f<<"#Delta#xi_{xm}\n#Delta#xi_{ym}\n#Delta#eta_{1}\n";
	f<<"#Delta#epsilon_{x}\n#Delta#epsilon_{y}\n";
	f<<"#Delta#delta_{a}\n#Delta#omega_{s}\n#DeltaT\n#Delta#nu_{s}"<<endl;
	f.close();

	int NL=0, Nc=0, Nr=0, Ncr=0;
	string s_NL = "cat $(find " + path + " -type f -name 'spin.dat' | sed -n '1p') | wc -l";
	string s_Nc = "find " + path + " -type d -name 'c*' | wc -l";
	NL = bash::to_int(s_NL);
	Nc = bash::to_int(s_Nc);
	cout<<"NL = "<<NL<<endl;
	cout<<"Nc = "<<Nc<<endl;

	int nth = 1, q = 0;

    vector<thread> vth;
    vector<array<double,14>> M(nth), e_M(nth);
    vector<acquisition> A;
    vector<proton_RF*> P_vec(nth,nullptr);
    //vector<vector<vector<double>>> P(nth,vector<vector<double>>(NL,vector<double>(9,0.0)));

    for(int i=0; i<Nc; i++)   
    {
        array<double,14> &ar = M[q%nth], &e_ar = e_M[q%nth];
        proton_RF* &P = P_vec[q%nth];

        auto func = [&ar,&e_ar,&P,&path,i,out,this]()
        {
            cout<<"Currently working on c"<<i+1<<endl;
            P = new proton_RF(path + "/c" + to_string(i+1));
            auto CSx = P->Courant_Snyder(0,out);
            auto CSy = P->Courant_Snyder(1,out);
            auto syn = P->longitudinal(out);
            auto ps = P->phase_slip(out);
            // auto BTx = P->betatron_tune(0,out);
            // auto BTy = P->betatron_tune(1,out);
            double nus = P->spin_tune(out);

            double kt1 = P->Kt1(), kt2 = P->Kt2(), kt3 = P->Kt3();
            double Cx = P->Cx(), Cy = P->Cy();
            double epx = CSx[0], epy = CSy[0];
            double dta = syn[3], om = syn[4];
            double eta0 = ps[1], eta1 = ps[2];
            double Cxm = P->Cx_m(), Cym = P->Cy_m();

            auto e_CSx = P->Courant_Snyder_errors(0,out);
            auto e_CSy = P->Courant_Snyder_errors(1,out);
            auto e_syn = P->longitudinal_errors(out);
            auto e_PSF = P->phase_slip_errors(out);
            auto e_BTx = P->betatron_tune_error(0,out);
            auto e_BTy = P->betatron_tune_error(1,out);
            double e_nus = P->spin_tune_error(out);
            double epx_e = e_CSx[0], epy_e = e_CSy[0];
            double dta_e = e_syn[3], om_e = e_syn[4];
            double eta0_e = e_PSF[1], eta1_e = e_PSF[2];
            double Cxm_e = e_BTx[1], Cym_e = e_BTy[1];

            ar =    {kt1, kt2, kt3, Cx,  Cy,  Cxm,   Cym,   eta1,   epx,   epy,   dta,   om,   eta0,   nus};
            e_ar =  {0.0, 0.0, 0.0, 0.0, 0.0, Cxm_e, Cym_e, eta1_e, epx_e, epy_e, dta_e, om_e, eta0_e, e_nus};

            delete P;
        };

        if(q>=q0)
        {
            vth.push_back(thread(func));
            //bash::to_int("sleep 0.1");
            this_thread::sleep_for(chrono::milliseconds(10));
        }
        else
        {
            cout<<"Thread "<<q%nth<<" not started: skipping line "<<q<<endl;
        }
        q++;
        if(q%nth == 0 || q == Nc)
        {
            f.open(path + "/line",fstream::out | fstream::app);
            fe.open(path + "/line_error",fstream::out | fstream::app);
            //cout<<"Size of M: "<<M.size()<<endl;
            //cout<<"Size of vth: "<<vth.size()<<endl;
            for(int k=0; k<vth.size(); k++) vth[k].join();
            cout<<vth.size()<<" threads joined!"<<endl;
            for(int k=0; k<vth.size(); k++)
            {
                for(int l=0; l<M[k].size()-1; l++) f<<M[k][l]<<"\t";
                f<<M[k][M[k].size()-1]<<endl;
                // cout<<M[k][5]<<"\t"<<M[k][6]<<endl;
                for(int l=0; l<e_M[k].size()-1; l++) fe<<e_M[k][l]<<"\t";
                fe<<e_M[k][e_M[k].size()-1]<<endl;
            }
            f.close();
            fe.close();
            //M.clear();
            M = vector<array<double,14>>(nth);
            vth.clear();
        }
    }
}