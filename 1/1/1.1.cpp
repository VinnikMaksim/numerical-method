#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

const int N = 3;

double A[N][N + 1] = { {2.75,1.78,1.11,15.71}, {3.28,0.71,1.15,43.78}, {1.15,2.7,3.58,37.11} };

double A_n[N][N] = { {2.75,1.78,1.11}, {3.28,0.71,1.15}, {1.15,2.7,3.58} };
double B[N] = { 15.71, 43.78, 37.11 };



void multiplyMatrixVector(double matrix[3][3], double vector[3], double result[3], int size = 3) {
    for (int i = 0; i < size; i++) {
        result[i] = 0;
        for (int j = 0; j < size; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

void ABinOneMatrix(double matrix[3][3], double vector[3], double result[3][4], int size = 3) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size+1; j++) {
            if (j == 3)
                result[i][j] = vector[i];
            else 
                result[i][j] = matrix[i][j];
        }
    }
}

void printMatrix(double matrix[N][N + 1])
{
    cout << setw(10);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            cout << matrix[i][j] << setw(10);
        }
        cout << endl;
    }
}

void swapRows(double matrix[N][N + 1], int row1, int row2) 
{
    for (int i = 0; i < N + 1; i++)
    {
        double temp = matrix[row1][i];
        matrix[row1][i] = matrix[row2][i];
        matrix[row2][i] = temp;
    }
}

void gauss(double matrix[N][N + 1])
{
    for (int i = 0; i < N - 1; i++)
    {

        int max_row = i;
        double max_val = abs(matrix[i][i]);
        for (int j = i + 1; j < N; j++)
        {
            if (abs(matrix[j][i]) > max_val)
            {
                max_row = j;
                max_val = abs(matrix[j][i]);
            }
        }


        if (max_row != i)
        {
            swapRows(matrix, i, max_row);
        }


        for (int j = i + 1; j < N; j++)
        {
            double coeff = -matrix[j][i] / matrix[i][i];
            for (int k = i; k < N + 1; k++)
            {
                matrix[j][k] += coeff * matrix[i][k];
            }
        }
    }
}

void backGauss(double matrix[N][N + 1], double solution[N])
{
    for (int i = N - 1; i >= 0; i--)
    {
        solution[i] = matrix[i][N];
        for (int j = i + 1; j < N; j++)
        {
            solution[i] -= matrix[i][j] * solution[j];
        }
        solution[i] /= matrix[i][i];
    }
}

void calcResidual(double matrix[N][N + 1], double solution[N], double residual[N])
{
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < N; j++) {
            sum += matrix[i][j] * solution[j];
        }
        residual[i] = sum - matrix[i][3];
    }
}

double findNormaOfRV(double F[N]) {
    int n = N;
    double norma = 0;

    for (int i = 0; i < n; i++) {
        norma += pow(F[i], 2);
    }
    norma = sqrt(norma);
    return norma;
}

double findCalcError(double x1[N], double x2[N]) {
    int n = N;
    double calc_error = 0;

    double max1 = 0, max2 = 0;
    for (int i = 0; i < n; i++) {
        if (x2[i] - x1[i] > max1)
            max1 = x2[i] - x1[i];
        if (x1[i] > max2)
            max2 = x1[i];
    }

    calc_error = max1 / max2;
    return calc_error;
}

int main()
{
    double c[N][N + 1];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            c[i][j] = A[i][j];
        }
    }
    setlocale(LC_ALL, "ru");
    double solution[N];
    double residual[N];

    cout << "Исходная матрица:" << endl;
    printMatrix(A);

    gauss(A);

    cout << "Треугольная матрица:" << endl;
    printMatrix(A);

    backGauss(A, solution);

    cout << "Решение системы:" << endl;
    for (int i = 0; i < N; i++)
    {
        cout << "x" << i + 1 << " = " << solution[i] << endl;
    }

    cout << "nevyazka: " << endl;
    calcResidual(c, solution, residual);
    for (int i = 0; i < N; i++)
    {
        cout << residual[i] << setprecision(20) << endl;
    }

    cout << "norma = " << findNormaOfRV(residual) << endl;


    double B_new[3];
    multiplyMatrixVector(A_n, solution, B_new);

    double AB[3][4];
    ABinOneMatrix(A_n, B_new, AB, 3);
    gauss(AB);

    double solution_1[N];
    backGauss(A, solution_1);

    cout << "Решение системы:" << endl;
    for (int i = 0; i < N; i++)
    {
        cout << "x" << i + 1 << " = " << solution_1[i] << endl;
    }


    cout << " cal error = " << findCalcError(solution, solution_1);



    return 0;
}




