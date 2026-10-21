#include <proton_RF.hh>

// --- Jacobian Calculation ---

/*
 * Calculates the Jacobian Matrix J (6x5).
 * Returns a matrix<6, 5> object.
 * Rows: alpha, beta, gamma, h, k, epsilon
 * Cols: A, B, C, D, E
 */
matrix<6,5> proton_RF::Courant_Snyder_Jacobian(matrix<5,1> M)
{
    double A = M[0]; double B = M[1]; double C = M[2]; double D = M[3]; double E = M[4];
    // 1. Calculate Intermediates for Z
    double termZ = A * C - (B * B) / 4.0;
    
    if (termZ <= 0) {
        cerr << "Error: Complex number detected (AC - B^2/4 <= 0). Derivatives undefined." << endl;
        return matrix<6, 5>(); // Returns zero matrix
    }

    double Z = pow(termZ, -0.5);

    // --- NEW: Calculate h and k using the corrected matrix form ---
    // [h; k] = - [2A B; B 2C]^-1 * [D; E]
    // Determinant (Det) = 4AC - B^2
    // Inverse = (1/Det) * [2C -B; -B 2A]
    // h = - (1/Det) * (2CD - BE) = (BE - 2CD) / Det
    // k = - (1/Det) * (-BD + 2AE) = (BD - 2AE) / Det
    
    double Det = 4.0 * A * C - B * B;
    double invDet = 1.0 / Det;
    double invDet2 = invDet * invDet; // For quotient rule: 1/Det^2

    double Num_h = B * E - 2.0 * C * D;
    double h = Num_h * invDet;

    double Num_k = B * D - 2.0 * A * E;
    double k = Num_k * invDet;

    // 2. Initialize Jacobian
    matrix<6, 5> J; 

    // 3. Partial Derivatives of Z
    // dZ/dX = -0.5 * Z^3 * d(termZ)/dX
    double preZ = -0.5 * pow(Z, 3);
    double dZ_dA = preZ * C;
    double dZ_dB = preZ * (-0.5 * B);
    double dZ_dC = preZ * A;
    double dZ_dD = 0;
    double dZ_dE = 0;

    // --- Row 0: Alpha (alpha = Z*B / 2) ---
    J[0][0] = (B / 2.0) * dZ_dA;              // dAlpha/dA
    J[0][1] = (B / 2.0) * dZ_dB + (Z / 2.0);  // dAlpha/dB
    J[0][2] = (B / 2.0) * dZ_dC;              // dAlpha/dC
    // J[0][3], J[0][4] are 0

    // --- Row 1: Beta (beta = Z*C) ---
    J[1][0] = C * dZ_dA;            // dBeta/dA
    J[1][1] = C * dZ_dB;            // dBeta/dB
    J[1][2] = C * dZ_dC + Z;        // dBeta/dC

    // --- Row 2: Gamma (gamma = Z*A) ---
    J[2][0] = A * dZ_dA + Z;        // dGamma/dA
    J[2][1] = A * dZ_dB;            // dGamma/dB
    J[2][2] = A * dZ_dC;            // dGamma/dC

    // --- Rows 3 & 4: h and k (Using Quotient Rule) ---
    // Rule: d(Num/Det)/dx = (dNum/dx * Det - Num * dDet/dx) / Det^2
    // Common derivatives for Det = 4AC - B^2
    double dDet_dA = 4.0 * C;
    double dDet_dB = -2.0 * B;
    double dDet_dC = 4.0 * A;
    // dDet_dD and dDet_dE are 0

    auto calc_quotient = [&](double Num, double dNum, double Det, double dDet) {
        return (dNum * Det - Num * dDet) * invDet2;
    };

    // Derivatives for h (Num_h = BE - 2CD)
    J[3][0] = calc_quotient(Num_h, 0.0, Det, dDet_dA);           // dh/dA
    J[3][1] = calc_quotient(Num_h, E, Det, dDet_dB);             // dh/dB
    J[3][2] = calc_quotient(Num_h, -2.0 * D, Det, dDet_dC);      // dh/dC
    J[3][3] = calc_quotient(Num_h, -2.0 * C, Det, 0.0);          // dh/dD
    J[3][4] = calc_quotient(Num_h, B, Det, 0.0);                 // dh/dE

    // Derivatives for k (Num_k = BD - 2AE)
    J[4][0] = calc_quotient(Num_k, -2.0 * E, Det, dDet_dA);      // dk/dA
    J[4][1] = calc_quotient(Num_k, D, Det, dDet_dB);             // dk/dB
    J[4][2] = calc_quotient(Num_k, 0.0, Det, dDet_dC);           // dk/dC
    J[4][3] = calc_quotient(Num_k, B, Det, 0.0);                 // dk/dD
    J[4][4] = calc_quotient(Num_k, -2.0 * A, Det, 0.0);          // dk/dE

    // --- Row 5: Epsilon ---
    // Epsilon depends on h and k, so we use the chain rule with the new dh/dx and dk/dx values
    double Bracket = h * D + k * E;
    auto calc_dEp = [&](double dZ, double dh, double dk, double dD, double dE) {
        double dBracket = dh * D + h * dD + dk * E + k * dE;
        return dZ - 0.5 * (dZ * Bracket + Z * dBracket);
    };

    J[5][0] = calc_dEp(dZ_dA, J[3][0], J[4][0], 0, 0); // dEp/dA
    J[5][1] = calc_dEp(dZ_dB, J[3][1], J[4][1], 0, 0); // dEp/dB
    J[5][2] = calc_dEp(dZ_dC, J[3][2], J[4][2], 0, 0); // dEp/dC
    J[5][3] = calc_dEp(dZ_dD, J[3][3], J[4][3], 1, 0); // dEp/dD
    J[5][4] = calc_dEp(dZ_dE, J[3][4], J[4][4], 0, 1); // dEp/dE

    return J;
}

/*
int main() {
    // 1. Define Inputs
    double A = 10.0;
    double B = 2.0;
    double C = 5.0;
    double D = 0.5;
    double E = 0.8;
    
    cout << "Inputs: A=" << A << ", B=" << B << ", C=" << C << ", D=" << D << ", E=" << E << endl << endl;
    
    // 2. Define Input Covariance Matrix (Sigma_u)
    // Using initializer list supported by your header
    matrix<5, 5> Sigma_u = {
        0.01, 0.00, 0.00, 0.00, 0.00,
        0.00, 0.02, 0.00, 0.00, 0.00,
        0.00, 0.00, 0.01, 0.00, 0.00,
        0.00, 0.00, 0.00, 0.05, 0.00,
        0.00, 0.00, 0.00, 0.00, 0.05
    };
    
    // 3. Compute Jacobian
    matrix<6, 5> J = calculateJacobian(A, B, C, D, E);
    
    // 4. Compute Output Covariance: Sigma_v = J * Sigma_u * J^T
    // Using .tp() from your header for transpose
    matrix<6, 6> Sigma_v = J * Sigma_u * J.tp();
    
    // 5. Output Results using your header's << operators
    cout << "--- Jacobian Matrix J ---" << endl;
    cout << J << endl;
    
    cout << "--- Output Covariance Matrix Sigma_v ---" << endl;
    cout << Sigma_v << endl;
    
    cout << "--- Estimated Uncertainties (Standard Deviations) ---" << endl;
    const char* names[] = {"alpha", "beta", "gamma", "h", "k", "epsilon"};
    
    // Using .diag() from your header
    auto variances = Sigma_v.diag();
    
    for(int i=0; i<6; ++i) {
        cout << setw(10) << names[i] << ": " << sqrt(variances[i]) << endl;
    }
    
    return 0;
}
*/