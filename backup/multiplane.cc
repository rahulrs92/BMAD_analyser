#include<lv1.hh>
#include<set>

void lv1::multiplane(acquisition &ia,string sample,string path,uint q0,bool out)
{
    fstream f;
    f.open(path + "/multiplane",fstream::in);
    if(f.is_open())
    {
        f.close();
        return;
    }
    f.close();
    f.open(path + "/multiplane_info",fstream::out | fstream::trunc);
    f<<"t_{f}\nt_{d}\nt_{s}\n#xi_{x}\n#xi_{y}\n";
    f<<"#epsilon_{x}\n#epsilon_{y}\n";
    f<<"#delta_{a}\n#omega_{s}\nT\n#nu_{s}"<<endl;
    f.close();

    int q = 0;
    int nth = 8;
    
    
    Data<double> D(sample);
    map<double,size_t> mr;
    int N = ia.Ncr;
    string line;
    f.open(sample,fstream::in);
    for(auto &a : D.col(3)) mr[a]++;
    cout<<"N = "<<N<<endl;
    cout<<"mr = ";
    for(auto &a : mr) cout<<a.first<<":"<<a.second<<" ";
    cout<<endl;
    f.close();
    int Nr = mr.size();
    cout<<"Nr = "<<Nr<<endl;
    
    vector<thread> vth;
    vector<array<double,11>> M(nth);
    vector<acquisition> A;
    for(int k=0; k<nth; k++) A.push_back(acquisition(ia,k));

    f.open(path + "/multiplane",fstream::out | fstream::trunc);
    int i = 0, j = 0;
    for(auto a : mr)    
    {
        valarray<double> vr = D.col(0)[D.col(3) == a.first];
        vr = vr + 1;
        //cout<<"r"<<i+1<<" contains \n"<<vr<<endl;
        for(auto c : vr)
        {
            array<double,11> &ar = M[q%nth];
            acquisition &aq = A[q%nth];
            if(q>=q0)
            {

                vth.push_back
                (
                    thread
                    (
                        [&ar,&aq,i,c,out,this]()
                        {
                            ar = point_RF(aq,i+1,c,out);
                            //cout<<"Contents of D: ";
                            //for(int l=0; l<4; l++) cout<<ar[l]<<" ";
                            //cout<<ar[4]<<endl;
                        }
                    )
                );
                bash::to_int("sleep 0.1");
            }
            else
            {
                cout<<"Thread "<<q%nth<<" not started: skipping line "<<q<<endl;
            }
            q++;
			if(q%nth == 0 || q == N)
			{
				//cout<<"Size of M: "<<M.size()<<endl;
				//cout<<"Size of vth: "<<vth.size()<<endl;
				for(int k=0; k<vth.size(); k++) vth[k].join();
                cout<<vth.size()<<" threads joined!"<<endl;
				for(int k=0; k<vth.size(); k++)
				{
					for(int l=0; l<10; l++) f<<M[k][l]<<"\t";
					f<<M[k][10]<<endl;
				}
				//M.clear();
                M = vector<array<double,11>>(nth);
				vth.clear();
			}
            j++;
        }
        j=0; i++;
    }
    f.close();
}