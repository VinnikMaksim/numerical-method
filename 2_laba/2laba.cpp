#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

const double E1 = 1e-9;
const double E2 = 1e-9;


double f1(double x1, double x2);
double f1dx1(double x1, double x2);
double f1dx2(double x1, double x2);
double f2(double x1, double x2);
double f2dx1(double x1, double x2);
double f2dx2(double x1, double x2);

vector<vector<double>> calculateMatrixJ(double x1, double x2, double k);
void solveNewtonsMethod(double x1, double x2, const double E1, const double E2, const int max_iter, const double M = 0);

vector<double> operator / (const vector<double>& v1, const vector<double>& v2);
vector<double> operator * (const vector<vector<double>>& m1, const vector<double>& v2);
vector<double> operator - (const vector<double>& v);
ostream& operator << (ostream& os, const vector<vector<double>>& A);
ostream& operator << (ostream& os, const vector<double>& A);

bool isVectorWithSameCoordinates(vector<double> v);
bool isMatrixWithProportionalityRows(vector<vector<double>> A, vector<double> b);
vector<double> solveGaussMethod(vector<vector<double>> A, vector<double> b);
void forwardGaussMethod(vector<vector<double>>& A, vector<double>& b);
void backwardGaussMethod(vector<vector<double>> A, vector<double> b, vector<double>& x);

double f1(double x1, double x2) { return sin(x1) - x2 - 1.32; }
double f1dx1(double x1, double x2) { return cos(x1); }
double f1dx2(double x1, double x2) { return -1.0; }
double f2(double x1, double x2) { return cos(x2) - x1 + 0.85; }
double f2dx1(double x1, double x2) { return -1; }
double f2dx2(double x1, double x2) { return -sin(x2); }

vector<double> operator / (const vector<double>& v1, const vector<double>& v2) {
    vector<double> result(v2.size());

    for (int i = 0; i < v2.size(); i++) {
        result[i] = v1[i] / v2[i];
    }
    return result;
}
vector<double> operator * (const vector<vector<double>>& m1, const vector<double>& v2) {
    vector<double> result(v2.size());

    for (int i = 0; i < result.size(); i++) {
        result[i] = 0;
        for (int j = 0; j < result.size(); j++) {
            result[i] += m1[i][j] * v2[j];
        }
    }
    return result;
}
vector<double> operator - (const vector<double>& v) {
    vector<double> result = v;
    for (int i = 0; i < v.size(); i++) {
        result[i] = -result[i];
    }
    return result;
}
ostream& operator << (ostream& os, const vector<vector<double>>& A)
{
    int n = A.size();

    cout << endl;
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[0].size(); j++) {
            os << A[i][j] << " ";
        }
        cout << endl;
    }
    return os;
}
ostream& operator << (ostream& os, const vector<double>& A)
{
    int n = A.size();

    os << "(" << A[0];
    for (int i = 1; i < n; i++) {
        os << "; " << A[i];
    }
    os << ")" << endl << endl;
    return os;
}

bool isVectorWithSameCoordinates(vector<double> v)
{
    double ram = v[0];
    for (int i = 1; i < v.size(); i++) {
        if (ram == v[i])
            continue;
        else
            return false;
    }
    return true;
}
bool isMatrixWithProportionalityRows(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    vector<vector<double>> A_b(n, vector<double>(n + 1, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) {
            if (j == n)
                A_b[i][j] = b[i];
            else
                A_b[i][j] = A[i][j];
        }
    }

    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (isVectorWithSameCoordinates(A_b[i] / A_b[j]))
                return true;
        }
    }

    return false;
}
vector<double> solveGaussMethod(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    if (isMatrixWithProportionalityRows(A, b)) {
        exit(-1);
    }

    forwardGaussMethod(A, b);

    vector<double> x(n);
    backwardGaussMethod(A, b, x);
    return x;
}
void forwardGaussMethod(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        int k = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(A[j][i]) > abs(A[k][i])) {
                k = j;
            }
        }

        swap(A[i], A[k]);
        swap(b[i], b[k]);

        double div = A[i][i];
        for (int j = i; j < n; j++) {
            A[i][j] /= div;
        }
        b[i] /= div;

        for (int j = i + 1; j < n; j++) {
            double mult = A[j][i];
            for (int k = i; k < n; k++) {
                A[j][k] -= mult * A[i][k];
            }
            b[j] -= mult * b[i];
        }
    }
}
void backwardGaussMethod(vector<vector<double>> A, vector<double> b, vector<double>& x) {
    int n = A.size();

    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
    }
}

vector<vector<double>> calculateMatrixJ(double x1, double x2, double k) {
    return vector<vector<double>> { {(f1(x1 + x1 * k, x2) - f1(x1, x2)) / k / x1, (f1(x1, x2 + x2 * k) - f1(x1, x2)) / k / x2},
        { (f2(x1 + x1 * k, x2) - f2(x1, x2)) / k / x1, (f2(x1, x2 + x2 * k) - f2(x1, x2)) / k / x2 } };
}

void solveNewtonsMethod(double x1, double x2, const double E1, const double E2, const int max_iter, const double M) {
    double delta1 = max(abs(f1(x1, x2)), abs(f2(x1, x2)));
    double delta2 = 1;

    cout << "Initial approximation : ( " << x1 << " ; " << x2 << " )\n";
    if (M != 0) { cout << "Relative increment = " << M << "\n"; }

    int iter = 0;
    while ((delta1 > E1 || delta2 > E2) && iter < max_iter) {
        iter++;

        cout << iter << "|\t" << "delta1| " << delta1 << "\t\t" << "delta2| " << delta2 << endl;

        vector<double> F{ f1(x1, x2), f2(x1, x2) };
        vector<vector<double>> J;

        if (M == 0) {
            J = { {f1dx1(x1, x2), f1dx2(x1, x2)}, {f2dx1(x1, x2), f2dx2(x1, x2)} };
        }
        else {
            J = calculateMatrixJ(x1, x2, M);
        }

        vector<double> solution = solveGaussMethod(J, -F);
        x1 += solution[0];
        x2 += solution[1];

        delta1 = abs(F[0]);
        for (int i = 1; i < F.size(); i++) {
            if (delta1 < abs(F[i])) {
                delta1 = abs(F[i]);
            }
        }
        double v1 = abs(x1) < 1 ? abs(solution[0]) : abs(solution[0] / x1);
        double v2 = abs(x2) < 1 ? abs(solution[1]) : abs(solution[1] / x2);
        delta2 = max(v1, v2);
    }
    cout << "\n" << iter << "|\t" << "x1| " << x1 << "\t" << "x2| " << x2 << "\n";
    cout << "____________________________________________________________________________\n\n";

    if (iter == max_iter)
    {
        cout << "limit has been reached !!!!!!" << endl;
    }
}

int main()
{
    cout.precision(15);
    double x1 = 1, x2 = 5;
    int max_iter = 100;


    // analytical jacobian

    // numerical jacobian M = 0.01, M = 0.05 è M = 0.1
    vector<double> vec_M = {0.01, 0.05, 0.1};
    for(auto m : vec_M){
        solveNewtonsMethod(x1, x2, E1, E2, max_iter, m);
    }

    return 0;
}