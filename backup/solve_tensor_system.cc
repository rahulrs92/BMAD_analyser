#include <lv2.hh> // Your provided header

// ============================================================================
// Tensor Solver Function
// ============================================================================

/**
 * @brief Solves y^i = y_0^i + A^i_j x^j + B^i_jk x^j x^k for x.
 * * @param data A 3x10 matrix containing the system parameters.
 * Rows correspond to i=1, 2, 3.
 * Cols: 0:y0, 1-3:A, 4-9:B components.
 * @param target_y The target vector y (3x1).
 * @return matrix<3, 1> The solution vector x.
 */
matrix<3, 1> lv2::solve_tensor_system(matrix<3, 10> data, matrix<3, 1> target_y)
{
    // --- 1. Parse Input Data ---
    // We decompose the 3x10 input into standard algebraic structures
    // to make the math readable.
    
    matrix<3, 1> y0;
    matrix<3, 3> A;
    // B is a rank-3 tensor. We store it as an array of 3 matrices (slices),
    // where B_slices[i] represents B^i_{jk}.
    matrix<3, 3> B_slices[3]; 

    for(int i = 0; i < 3; i++)
    {
        // Extract y0 (Column 0)
        y0[i] = data[i][0];

        // Extract A (Columns 1, 2, 3)
        A[i][0] = data[i][1];
        A[i][1] = data[i][2];
        A[i][2] = data[i][3];

        // Extract B (Columns 4-9) and enforce symmetry
        // Mapping: 4->11, 5->22, 6->33, 7->12, 8->23, 9->31
        double b11 = data[i][4];
        double b22 = data[i][5];
        double b33 = data[i][6];
        double b12 = data[i][7] / 2.0; // Symmetric
        double b23 = data[i][8] / 2.0; // Symmetric
        double b31 = data[i][9] / 2.0; // Symmetric

        // Construct the 3x3 matrix for this specific i
        matrix<3, 3> Bi;
        Bi[0][0] = b11; Bi[0][1] = b12; Bi[0][2] = b31;
        Bi[1][0] = b12; Bi[1][1] = b22; Bi[1][2] = b23;
        Bi[2][0] = b31; Bi[2][1] = b23; Bi[2][2] = b33;
        
        B_slices[i] = Bi;
    }

    // --- 2. Newton-Raphson Loop ---
    
    matrix<3, 1> x = tensor_Xi(0);// (A^-1) * (target_y - y0); // Initial guess (ignores B)
    int max_iter = 100;
    double tol = 1e-9;

    for(int iter = 0; iter < max_iter; iter++) {
        
        // -- Calculate Residual f(x) --
        // f^i = y_0^i + A^i_j x^j + x^T * B^i * x - target_y^i
        // Note: x^T * B^i * x is a quadratic form (scalar result)
        
        matrix<3, 1> Ax = A * x;
        matrix<3, 1> Bxx; // This will hold the vector of quadratic results
        
        for(int i = 0; i < 3; i++) Bxx[i] = x.tp() * B_slices[i] * x;

        matrix<3, 1> f = y0 + Ax + Bxx - target_y;

        // -- Check Convergence --
        // Simple Euclidean norm check
        double err = 0;
        for(int k=0; k<3; k++) err += f[k]*f[k];
        if(sqrt(err) < tol)
        {
            // cout << "Final y: " << y0 + Ax + Bxx << endl;
            cout << "data = \n" << data << endl;
            return x; // Converged
        }

        // -- Calculate Jacobian J(x) --
        // J^i_j = A^i_j + 2 * B^i_jk * x^k
        // In matrix terms: Row_i(J) = Row_i(A) + 2 * (B^i * x)^T
        
        matrix<3, 3> J;
        for(int i = 0; i < 3; i++) {
            // Calculate derivative of quadratic part: 2 * B^i * x
            matrix<3, 1> Bx = B_slices[i] * x;
            
            // Fill the i-th row of J
            for(int j = 0; j < 3; j++) {
                J[i][j] = A[i][j] + 2.0 * Bx[j];
            }
        }

        // -- Update Step --
        // dx = - J^-1 * f
        // Your library supports inv() for square matrices
        matrix<3, 1> dx = -(J^-1) * f;
        
        x = x + dx;
    }

    cerr << "Warning: Newton solver did not converge." << endl;
    return x;
}

// ============================================================================
// Main Driver
// ============================================================================

/*
int main() {
    // 1. Load Data from file into a 3x10 matrix
    matrix<3, 10> data;
    ifstream file("data.txt");
    
    if(!file.is_open()) {
        cerr << "Error: data.txt not found." << endl;
        return 1;
    }

    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 10; j++) {
            file >> data[i][j];
        }
    }
    file.close();

    // 2. Define target Y (Solving for roots, so y = 0)
    matrix<3, 1> target_y = {0.0, 0.0, 0.0};

    // 3. Solve
    cout << "Solving system..." << endl;
    matrix<3, 1> result = solve_tensor_system(data, target_y);

    // 4. Output
    cout << "Found Solution x:" << endl;
    cout << result << endl;

    return 0;
}
*/