#include <lv2.hh>
#include <vector>
#include <cmath>

// ============================================================================
// Tensor Solver Function
// ============================================================================

/**
 * @brief Solves y^i = y_0^i + A^i_j x^j + B^i_jk x^j x^k for multiple x solutions.
 * * Strategy:
 * 1. Parse parameters y0, A, B from the input data.
 * 2. Generate a grid of initial guesses around the origin.
 * 3. Run Newton-Raphson for each guess.
 * 4. Store unique valid solutions.
 * 5. Return a matrix of valarrays containing all found solutions.
 */
matrix<3, 1, valarray<double>>& lv1::solve_tensor_system(matrix<3, 10> data, matrix<3, 1> target_y)
{
    // --- 1. Parse Input Data ---
    matrix<3, 1> y0;
    matrix<3, 3> A;
    matrix<3, 3> B_slices[3]; 

    for(int i = 0; i < 3; i++)
    {
        y0[i] = data[i][0];
        A[i][0] = data[i][1]; A[i][1] = data[i][2]; A[i][2] = data[i][3];

        // Note: Using your specific symmetric logic (/ 2.0)
        double b11 = data[i][4];
        double b22 = data[i][5];
        double b33 = data[i][6];
        double b12 = data[i][7] / 2.0; 
        double b23 = data[i][8] / 2.0; 
        double b31 = data[i][9] / 2.0; 

        matrix<3, 3> Bi;
        Bi[0][0] = b11; Bi[0][1] = b12; Bi[0][2] = b31;
        Bi[1][0] = b12; Bi[1][1] = b22; Bi[1][2] = b23;
        Bi[2][0] = b31; Bi[2][1] = b23; Bi[2][2] = b33;
        B_slices[i] = Bi;
    }

    // --- 2. Setup Grid Search ---
    vector<double> us_x, us_y, us_z;
    
    // Grid parameters
    double range_min = -5.0;
    double range_max = 5.0;
    double step = 1.0; 
    
    // Helper to check uniqueness
    // Iterate by const-ref to avoid unnecessary copies
    auto is_unique = [&](matrix<3, 1> sol) -> bool
    {
        for(int i = 0; i < us_x.size(); i++)
            if((matrix<3,1>{us_x[i], us_y[i], us_z[i]} - sol).norm() < 1e-6) return false; // Too close to existing solution
        return true;
    };

    // --- 3. Grid Loop ---
    for(double g1 = range_min; g1 <= range_max; g1 += step) {
        for(double g2 = range_min; g2 <= range_max; g2 += step) {
            for(double g3 = range_min; g3 <= range_max; g3 += step) {
                
                // Initial guess
                matrix<3, 1> x; 
                x[0] = g1; x[1] = g2; x[2] = g3;

                // Newton-Raphson variables
                int max_iter = 50; // Iterations per guess
                double tol = 1e-9;
                bool converged = false;

                for(int iter = 0; iter < max_iter; iter++) {
                    // Calculate Residual f(x)
                    matrix<3, 1> Ax = A * x;
                    matrix<3, 1> Bxx;
                    for(int i = 0; i < 3; i++) Bxx[i] = x.tp() * B_slices[i] * x;

                    matrix<3, 1> f = y0 + Ax + Bxx - target_y;

                    // Check Convergence
                    double err = 0;
                    for(int k=0; k<3; k++) err += f[k]*f[k];
                    if(std::sqrt(err) < tol) {
                        converged = true;
                        break;
                    }

                    // Jacobian J(x)
                    matrix<3, 3> J;
                    for(int i = 0; i < 3; i++) {
                        matrix<3, 1> Bx = B_slices[i] * x;
                        for(int j = 0; j < 3; j++) {
                            J[i][j] = A[i][j] + 2.0 * Bx[j];
                        }
                    }

                    // Update Step: dx = - J^-1 * f
                    // Using try-catch or safe inversion check if needed, 
                    // though your header likely handles singular matrices via throw/error or 0-div
                    try {
                        matrix<3, 1> dx = -(J ^ -1) * f; 
                        x = x + dx;
                    } catch(...) {
                        converged = false;
                        cerr << "Singular Jacobian encountered, aborting this guess." << endl;
                        break; // Singular matrix, abort this guess
                    }
                }

                // Store if valid and unique
                if(converged && is_unique(x))
                {
                    // Ensure we push an owning copy (binds its internal storage correctly)
                    cout<<"Found unique solution: \n"<<x<<endl;
                    us_x.push_back(x[0]);
                    us_y.push_back(x[1]);
                    us_z.push_back(x[2]);
                }
            }
        }
    }

    // --- 4. Pack Results into Output Matrix ---
    size_t count = us_x.size();
    
    // Debug print
    std::cout << "Found " << count << " unique solutions." << std::endl;
    for(size_t s = 0; s < count; s++) cout << "Solution " << s+1 << ":\n" << matrix<3,1>{us_x[s], us_y[s], us_z[s]} << endl;
    
    // Prepare valarrays for x1, x2, x3 components
    std::valarray<double> vx1(count), vx2(count), vx3(count);

    matrix<3, 1, std::valarray<double>> &output = tensor_solutions; // Assuming this is a member variable to hold results
    if(count == 0)
    {
        std::cout << "Warning: No valid solutions found for the given tensor system." << std::endl;
        output = matrix<3, 1, valarray<double>>(valarray<double>(0.0,1));
    }
    else
    {
        for(size_t s = 0; s < count; s++)
        {
            vx1[s] = us_x[s];
            vx2[s] = us_y[s];
            vx3[s] = us_z[s];
        }

        output[0] = vx1;
        output[1] = vx2;
        output[2] = vx3;
    }
    return output;
}