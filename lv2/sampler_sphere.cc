#include <lv2.hh>

matrix<3,1,valarray<double>> lv2::sampler_sphere(int N, double r, matrix<3,1> o, string filename, bool out)
{
    matrix<3,1,valarray<double>> points(valarray<double>(0.0, N));
    
    // Generate points with random volume distribution in a cube, select those in the inscribed sphere
    int count = 0;
    while (count < N)
    {
        double x = (2.0 * rand() / RAND_MAX - 1.0) * r;
        double y = (2.0 * rand() / RAND_MAX - 1.0) * r;
        double z = (2.0 * rand() / RAND_MAX - 1.0) * r;
        if (x*x + y*y + z*z <= r*r)
        {
            points[0][count] = x;
            points[1][count] = y;
            points[2][count] = z;
            ++count;
        }
    }

    // Shift points to be centered around o
    points = points + o;

    // Transformation to field coordinates (COSY default working point)
    // matrix<3,1> xi0 = {-4.24664,-3.93523,0.488338};
    // matrix<3,3> Xi = {2.62575,2.91311,9.09487,-1.86668,-13.2338,-5.81579,-0.11602,-1.30336,-5.03159}, E;
    
    // Transform to field coordinates (comment for a sphere in field space)
    // points = (Xi^-1) * (points - xi0);     

    // valarray<double> e(3);
    // Xi.Eigens(e,E);
    // cout<<"Eigenvectors of Xi: \n"<<E<<endl;

    
	if(out)
	{
		fstream f;
		f.open(filename.data(), fstream::out | fstream::trunc);
		for(int i=0; i<N; i++)
			f<<points[0][i]<<"\t"<<points[1][i]<<"\t"<<points[2][i]<<endl;
		f.close();
		f.open((filename + "_stats").data(), fstream::out | fstream::trunc);
		f<<"size of sample: "<<points[0].size()<<endl;
		f<<endl;
		f<<"sample x_max: "<<points[0].max()<<endl;
		f<<"sample y_max: "<<points[1].max()<<endl;
		f<<"sample z_max: "<<points[2].max()<<endl;
		f<<endl;
		f<<"sample x_min: "<<points[0].min()<<endl;
		f<<"sample y_min: "<<points[1].min()<<endl;
		f<<"sample z_min: "<<points[2].min()<<endl;
		f.close();

		new TCanvas("clatplot","clatplot",1900,1000);
		plot2<double> lat_plot(points[0],points[1],points[2],"Sample Lattice",
										"t_{s} (Tesla/m^{2})","t_{l} (Tesla/m^{2})","t_{g} (Tesla/m^{2})");
		//lat_plot.axis_equal();
    	lat_plot->Draw("P0");
	}

    return points;
}