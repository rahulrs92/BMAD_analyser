#include <lv2.hh>
#include <random>
#include <numeric> // For std::accumulate

// Function to generate random sample points
matrix<3,1,valarray<double>> lv2::sampler(int n, int l, array<double,3> o, double r, string filename, bool out)
{
    matrix<3,1,valarray<double>> points = valarray<double>(0.0, n);

    // Generate evenly spaced z-values
    std::vector<double> zLv(l);
    for (int i = 0; i < l; ++i)
        zLv[i] = -r + (2.0 * r / (l - 1)) * i;

    //cout << "z-levels: ";
    //for (const auto& a : zLv) cout << a << " ";
    //cout << endl;

    // Compute the area of the circle at each z-level
    std::vector<double> areas(l);
    for (int i = 0; i < l; ++i)
    {
        double z = zLv[i];
        //double rMax = std::sqrt(r * r - z * z); // Radius of the circle at this z-level
        double rMax = std::sqrt(std::max(0.0, r * r - z * z)); // Safe radius calculation
        cout << "z: " << z << ",\trMax: " << rMax << endl;
        areas[i] = M_PI * rMax * rMax; // Area of the circle
    }

    // Normalize the areas to create a probability distribution
    double totalArea = std::accumulate(areas.begin(), areas.end(), 0.0);
    for (double& area : areas) area /= totalArea; // Normalize to sum to 1
    
    // Create a discrete distribution based on the normalized areas
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> zDist(areas.begin(), areas.end());
    std::uniform_real_distribution<> dist(0.0, 1.0); // Uniform distribution for radius within the circle
    std::uniform_real_distribution<> angleDist(0.0, 2 * M_PI);

    for (int i = 0; i < n; i++) {
        // Select a z-level based on the area-weighted probability
        int zIndex = zDist(gen);
        double z = zLv[zIndex];

        // Compute the radius of the circle at this z-level
        //double rMax = std::sqrt(r * r - z * z);
        double rMax = std::sqrt(std::max(0.0, r * r - z * z)); // Safe radius calculation

        // Generate random (x, y) within the circle
        double lr = rMax * std::sqrt(dist(gen)); // sqrt for uniform density
        //cout << "Sample point " << i << ": z = " << z << ", rMax = " << rMax << ", lr = " << lr << endl;
        double theta = angleDist(gen);
        points[0][i] = lr * std::cos(theta);
        points[1][i] = lr * std::sin(theta);
        points[2][i] = z;

        // Apply the offset
        points[0][i] += o[0];
        points[1][i] += o[1];
        points[2][i] += o[2];

        //cout << "Sample point " << i << ": (" << points[0][i] << ", " << points[1][i] << ", " << points[2][i] << ")" << endl;
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

		new TCanvas("clatplot","clatplot",1900,1000);
		plot2<double> lat_plot(points[0],points[1],points[2],"Sample Lattice",
										"t_{s} (Tesla/m^{2})","t_{l} (Tesla/m^{2})","t_{g} (Tesla/m^{2})");
		//lat_plot.axis_equal();
    	lat_plot->Draw("P0");
	}
    return points;
}
