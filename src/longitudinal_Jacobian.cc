#include <proton_RF.hh>

// --- Jacobian Calculation ---

/*
 * Calculates the Jacobian Matrix J (5x5).
 * Returns a matrix<5, 5> object.
 * Rows: z_m, z_a, delta_m, delta_a, omega_s
 * Cols: A, B, C, D, eta_0
 */
matrix<5, 5> proton_RF::longitudinal_Jacobian(matrix<5,1> M)
{
    double A = M[0]; double B = M[1]; double C = M[2]; double D = M[3]; double eta_0 = M[4];
    
    // --- 1. Calculate Intermediates and Values ---
    
    // Safety checks for denominators and square roots
    // if (A <= 0 || B <= 0 || eta_0 == 0) {
    //     cerr << "Error: Invalid inputs (A, B must be > 0; eta_0 != 0)." << endl;
    //     return matrix<5, 5>();
    // }

    double Z = 1.0 + (C * C) / (4.0 * A) + (D * D) / (4.0 * B);
    
    // if (Z < 0) {
    //     cerr << "Error: Z is negative, cannot take square root." << endl;
    //     return matrix<5, 5>();
    // }

    double z_a = sqrt(Z / A);
    double delta_a = sqrt(Z / B);
    
    // Needed for derivatives of z_a and delta_a
    // d(sqrt(u))/dx = (1 / 2*sqrt(u)) * du/dx = (1 / 2*Output) * du/dx
    double prefactor_za = 1.0 / (2.0 * z_a);
    double prefactor_delta_a = 1.0 / (2.0 * delta_a);

    // --- 2. Initialize Jacobian (5 rows, 5 cols) ---
    matrix<5, 5> J; 

    // ------------------------------------------------------------------------
    // Row 0: z_m = -C / (2A)
    // ------------------------------------------------------------------------
    // d(z_m)/dA = -C/2 * (-1/A^2) = C / (2A^2)
    J[0][0] = C / (2.0 * A * A);      
    // d(z_m)/dC = -1 / (2A)
    J[0][2] = -1.0 / (2.0 * A);       
    // d/dB, d/dD, d/deta_0 = 0

    // ------------------------------------------------------------------------
    // Row 1: z_a = sqrt(Z/A)
    // ------------------------------------------------------------------------
    // Let U = Z/A = 1/A + C^2/(4A^2) + D^2/(4AB)
    // dz_a/dX = prefactor_za * dU/dX
    
    double dU_dA = -1.0/(A*A) - (2.0*C*C)/(4.0*A*A*A) - (D*D)/(4.0*A*A*B); // dU/dA
    double dU_dB = -(D*D) / (4.0*A*B*B);                                   // dU/dB
    double dU_dC = (2.0*C) / (4.0*A*A);                                    // dU/dC = C/(2A^2)
    double dU_dD = (2.0*D) / (4.0*A*B);                                    // dU/dD = D/(2AB)
    
    J[1][0] = prefactor_za * dU_dA;
    J[1][1] = prefactor_za * dU_dB;
    J[1][2] = prefactor_za * dU_dC;
    J[1][3] = prefactor_za * dU_dD;
    // d/deta_0 = 0

    // ------------------------------------------------------------------------
    // Row 2: delta_m = -D / (2B)
    // ------------------------------------------------------------------------
    // d(delta_m)/dB = -D/2 * (-1/B^2) = D / (2B^2)
    J[2][1] = D / (2.0 * B * B);
    // d(delta_m)/dD = -1 / (2B)
    J[2][3] = -1.0 / (2.0 * B);
    // d/dA, d/dC, d/deta_0 = 0

    // ------------------------------------------------------------------------
    // Row 3: delta_a = sqrt(Z/B)
    // ------------------------------------------------------------------------
    // Let V = Z/B = 1/B + C^2/(4AB) + D^2/(4B^2)
    // d(delta_a)/dX = prefactor_delta_a * dV/dX

    double dV_dA = -(C*C) / (4.0*A*A*B);                                   // dV/dA
    double dV_dB = -1.0/(B*B) - (C*C)/(4.0*A*B*B) - (2.0*D*D)/(4.0*B*B*B); // dV/dB
    double dV_dC = (2.0*C) / (4.0*A*B);                                    // dV/dC = C/(2AB)
    double dV_dD = (2.0*D) / (4.0*B*B);                                    // dV/dD = D/(2B^2)

    J[3][0] = prefactor_delta_a * dV_dA;
    J[3][1] = prefactor_delta_a * dV_dB;
    J[3][2] = prefactor_delta_a * dV_dC;
    J[3][3] = prefactor_delta_a * dV_dD;
    // d/deta_0 = 0

    // ------------------------------------------------------------------------
    // Row 4: omega_s = (c / eta_0) * sqrt(A/B)
    // ------------------------------------------------------------------------
    // omega_s = c * eta_0^-1 * A^0.5 * B^-0.5

    // d/dA = (c/eta_0) * B^-0.5 * 0.5 * A^-0.5 = (c / 2*eta_0 * sqrt(AB))
    J[4][0] = L[mode] / (2.0 * eta_0 * sqrt(A * B));
    
    // d/dB = (c/eta_0) * A^0.5 * -0.5 * B^-1.5 = - (c * sqrt(A)) / (2*eta_0 * B^1.5)
    J[4][1] = - L[mode] * sqrt(A/B) / (2.0 * eta_0 * B);
    
    // d/deta_0 = c * sqrt(A/B) * -1 * eta_0^-2 = -omega_s / eta_0
    double omega_s = (L[mode] / eta_0) * sqrt(A / B);
    J[4][4] = -omega_s / eta_0;
    
    // d/dC, d/dD = 0

    return J;
}

/*
int main() {
    // 1. Define Inputs and Constants
    // Replace these with your actual operational points
    double A = 10.0;
    double B = 5.0;
    double C = 2.0;
    double D = 3.0;
    double eta_0 = 377.0; // Impedance of free space approx
    const double c = 3.0e8;

    cout << "Inputs: A=" << A << ", B=" << B << ", C=" << C << ", D=" << D << ", eta_0=" << eta_0 << ", c=" << c << endl << endl;

    // 2. Define Input Covariance Matrix (Sigma_u)
    // Order: A, B, C, D, eta_0
    // Example values - Replace with your data
    matrix<5, 5> Sigma_u = {
        1e-4, 0.0,  0.0,  0.0,  0.0, // A
        0.0,  1e-4, 0.0,  0.0,  0.0, // B
        0.0,  0.0,  1e-4, 0.0,  0.0, // C
        0.0,  0.0,  0.0,  1e-4, 0.0, // D
        0.0,  0.0,  0.0,  0.0,  1e-2 // eta_0
    };

    // 3. Compute Jacobian
    matrix<5, 5> J = calculateJacobian(A, B, C, D, eta_0, c);
    
    // Check for validity
    if (J[0][0] == 0 && J[0][2] == 0 && J[1][0] == 0) {
        // Simple heuristic check if J is all zeros due to error return
         // In a real app, handle the empty matrix return explicitly
    }

    // 4. Compute Output Covariance: Sigma_v = J * Sigma_u * J^T
    matrix<5, 5> Sigma_v = J * Sigma_u * J.tp();

    // 5. Output Results 
    cout << "--- Jacobian Matrix J (5x5) ---" << endl;
    cout << J << endl;

    cout << "--- Output Covariance Matrix Sigma_v (5x5) ---" << endl;
    cout << Sigma_v << endl;

    cout << "--- Estimated Uncertainties (Standard Deviations) ---" << endl;
    const char* names[] = {"z_m", "z_a", "delta_m", "delta_a", "omega_s"};
    
    auto variances = Sigma_v.diag();
    
    for(int i=0; i<5; ++i) {
        cout << setw(10) << names[i] << ": " << scientific << setprecision(4) << sqrt(variances[i]) << endl;
    }

    return 0;
}
*/