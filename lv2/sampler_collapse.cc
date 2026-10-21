#include <lv2.hh>
#include <random>
#include <numeric> // For std::accumulate

using namespace std;

// Function to generate random sample points uniformly in a sphere, then collapse z to o[2]
matrix<3,1,valarray<double>> lv2::sampler_collapse(int n, array<double,3> o, double r, string filename, bool out)
{
    matrix<3,1,valarray<double>> points = valarray<double>(0.0, n);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist(-1.0, 1.0);
    uniform_real_distribution<> dist01(0.0, 1.0);

    for (int i = 0; i < n; ++i)
    {
        // Marsaglia method for uniform points in a sphere
        double x, y, z;
        while (true) {
            x = dist(gen);
            y = dist(gen);
            z = dist(gen);
            if (x*x + y*y + z*z <= 1.0) break;
        }
        double scale = cbrt(dist01(gen)); // Uniform in volume
        x *= r * scale;
        y *= r * scale;
        z *= r * scale;
        points[0][i] = x + o[0];
        points[1][i] = y + o[1];
        points[2][i] = o[2]; // Collapse z to o[2]
    }

    if(out)
    {
        fstream f;
        f.open(filename.data(), fstream::out | fstream::trunc);
        for(int i=0; i<n; i++)
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

        new TCanvas("clatplot_zcollapse","clatplot_zcollapse",1900,1000);
        plot2<double> lat_plot(points[0],points[1],points[2],"Sample Lattice Collapsed Z",
                                        "t_{s} (Tesla/m^{2})","t_{l} (Tesla/m^{2})","t_{g} (Tesla/m^{2})");
        lat_plot->Draw("P0");
    }
    return points;
}
