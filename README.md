#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

// Check if the matrix is symmetric
bool isSymmetric(const vector<vector<double>>& A, double tol = 1e-9) {
    int n = (int)A.size();
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (fabs(A[i][j] - A[j][i]) > tol) return false;
    return true;
}

// Check if the matrix is positive definite using Cholesky
bool isPositiveDefinite(const vector<vector<double>>& A, double tol = 1e-12) {
    int n = (int)A.size();
    vector<vector<double>> L(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0.0;
            for (int k = 0; k < j; ++k) sum += L[i][k] * L[j][k];
            if (i == j) {
                double val = A[i][i] - sum;
                if (val <= tol) return false;
                L[i][j] = sqrt(val);
            } else {
                if (fabs(L[j][j]) <= tol) return false;
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }
    return true;
}

// Cholesky decomposition (fills lower-triangular L)
void choleskyDecomposition(const vector<vector<double>>& A, vector<vector<double>>& L) {
    int n = (int)A.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0.0;
            for (int k = 0; k < j; ++k) sum += L[i][k] * L[j][k];
            if (i == j) L[i][j] = sqrt(A[i][i] - sum);
            else L[i][j] = (A[i][j] - sum) / L[j][j];
        }
        for (int j = i+1; j < n; ++j) L[i][j] = 0.0; // explicitly zero above diag
    }
}

// Forward substitution: L * y = b
vector<double> forwardSubstitution(const vector<vector<double>>& L, const vector<double>& b) {
    int n = (int)L.size();
    vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int k = 0; k < i; ++k) sum += L[i][k] * y[k];
        y[i] = (b[i] - sum) / L[i][i];
    }
    return y;
}

// Back substitution: L^T * x = y  (use L[j][i] for L^T entries)
vector<double> backSubstitutionLt(const vector<vector<double>>& L, const vector<double>& y) {
    int n = (int)L.size();
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += L[j][i] * x[j]; // L^T[i][j] = L[j][i]
        x[i] = (y[i] - sum) / L[i][i];
    }
    return x;
}

int main() {
    cout << fixed << setprecision(6);
    int N;
    cout << "Enter the size of the matrix N: ";
    if (!(cin >> N) || N <= 0) {
        cout << "Invalid size.\n";
        return 1;
    }

    vector<vector<double>> A(N, vector<double>(N));
    cout << "Enter the elements of matrix A (row by row):\n";
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            cin >> A[i][j];

    vector<double> b(N);
    cout << "Enter the right-hand side vector b (N elements):\n";
    for (int i = 0; i < N; ++i) cin >> b[i];

    // checks
    if (!isSymmetric(A)) {
        cout << "\nERROR: The matrix is NOT symmetric. Cholesky decomposition cannot be applied.\n";
        return 0;
    }

    if (!isPositiveDefinite(A)) {
        cout << "\nERROR: The matrix is NOT positive definite. Cholesky decomposition cannot be applied.\n";
        return 0;
    }

    cout << "\nThe matrix is symmetric and positive definite. Proceeding with Cholesky...\n";

    // Cholesky
    vector<vector<double>> L(N, vector<double>(N, 0.0));
    choleskyDecomposition(A, L);

    cout << "\nLower triangular matrix L:\n";
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) cout << L[i][j] << "\t";
        cout << "\n";
    }

    // Solve Ly = b
    vector<double> y = forwardSubstitution(L, b);
    cout << "\nSolution of L * y = b -> y:\n";
    for (double v : y) cout << v << "\t";
    cout << "\n";

    // Solve L^T x = y
    vector<double> x = backSubstitutionLt(L, y);
    cout << "\nSolution of A * x = b -> x:\n";
    for (double v : x) cout << v << "\t";
    cout << "\n";
    // Verification:

compute L * L^T
    cout << "\nVerification (L * L^T):\n";
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double sum = 0.0;
            for (int k = 0; k < N; ++k) sum += L[i][k] * L[j][k];
            cout << sum << "\t";
        }
        cout << "\n";
    }

    return 0;
}
