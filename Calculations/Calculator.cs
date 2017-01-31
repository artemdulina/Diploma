using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Calculations
{
    public static class Calculator
    {
        public static double[,] GetSystemOfEquations(string leftExpression, Func<double, double> rightFunction)
        {
            return new double[5, 5];
        }

        public static double Integrate(double a, double b, double p, Func<double, double> firstFunction, Func<double, double> secondFunction)
        {
            double result = 0;
            double h = (b - a) / p;

            for (int i = 1; i < p; i++)
            {
                result += firstFunction(a + h * i) * secondFunction(a + h * i);
            }
            result = h * (result + (firstFunction(a) * secondFunction(a) + firstFunction(b) * secondFunction(b)) / 2);

            return result;
        }

        public static double[] FivediagonalMatrixAlgorithm(double[][] matrix, double[] b, int n)
        {
            double[] answer = new double[n + 1];
            double[] p = new double[n + 2];
            double[] d = new double[n + 2];
            double[] q = new double[n + 2];
            double[] r = new double[n + 2];

            p[1] = -matrix[0][1] / matrix[0][0];
            d[1] = matrix[1][1] + matrix[1][0] * p[1];
            q[1] = matrix[0][2] / matrix[0][0];
            r[1] = b[0] / matrix[0][0];

            p[2] = (-matrix[1][2] + q[1] * matrix[1][0]) / d[1];
            r[2] = (b[1] - matrix[1][0] * r[1]) / d[1];
            d[2] = matrix[2][2] - matrix[2][0] * q[1] + p[2] * (matrix[2][0] * p[1] + matrix[2][1]);
            q[2] = matrix[1][3] / d[1];

            int t = 3;
            for (int i = 3; i < n; i++)
            {
                p[i] = (-matrix[i - 1][3] + q[i - 1] * (matrix[i - 1][0] * p[i - 2] + matrix[i - 1][1])) / d[i - 1];
                r[i] = (b[i - 1] - matrix[i - 1][0] * r[i - 2] - r[i - 1] * (matrix[i - 1][0] * p[i - 2] + matrix[i - 1][1])) / d[i - 1];
                d[i] = matrix[i][2] - matrix[i][0] * q[i - 1] + p[i] * (matrix[i][0] * p[i - 1] + matrix[i][1]);
                q[i] = matrix[i - 1][4] / d[i - 1];
                t++;
            }

            p[t] = (-matrix[t - 1][3] + q[t - 1] * (matrix[t - 1][0] * p[t - 2] + matrix[t - 1][1])) / d[t - 1];
            r[t] = (b[t - 1] - matrix[t - 1][0] * r[t - 2] - r[t - 1] * (matrix[t - 1][0] * p[t - 2] + matrix[t - 1][1])) / d[t - 1];
            d[t] = matrix[t][2] - matrix[t][0] * q[t - 1] + p[t] * (matrix[t][0] * p[t - 1] + matrix[t][1]);

            t++;
            r[t] = (b[t - 1] - matrix[t - 1][0] * r[t - 2] - r[t - 1] * (matrix[t - 1][0] * p[t - 2] + matrix[t - 1][1])) / d[t - 1];

            answer[n] = r[n + 1];
            answer[n - 1] = p[n] * answer[n] + r[n];

            for (int i = n - 2; i >= 0; i--)
            {
                answer[i] = p[i + 1] * answer[i + 1] - q[i + 1] * answer[i + 2] + r[i + 1];
            }

            return answer;
        }
    }
}



/*double u(double x)
{
    return (pow(x, 5) + 256 * pow(x, 2) - 517 * x - 212) / 60.;
}

double fx(double x)
{
    return x * x;
}

double eZeroone(double x)
{
    return 2 - x * x / (h * h);
}

double eOneone(double x)
{
    return (x * x - 4 * h * x + 4 * h * h) / (h * h);
}

double eZerotwo(double x)
{
    return -x * x / (2 * h * h) + 1.5;
}

double eOnetwo(double x)
{
    return x * x / (4 * h * h) - 3 * x / (2 * h) + 9. / 4;
}

double eNminusOneone(double x)
{
    return (x * x - 2 * h * (n - 3) * x + pow((n - 3) * h, 2)) / (4 * h * h);
}

double eNone(double x)
{
    return (x + h - (n - 1) * h) / h;
}

double eNminusOnetwo(double x)
{
    return (x * x - 2 * h * (n - 2) * x + pow((n - 2) * h, 2)) / (h * h);
}

double eNtwo(double x)
{
    return (2 * x + h - 2 * (n - 1) * h) / h;
}

double integrate(double a, double b, double p, double(*functionA)(double), double(* functionB)(double))
{
	double result = 0;
double h = (b - a) / p;
	for (int i = 1; i<p; i++)
	{
		result += functionA(a + h* i)* functionB(a + h* i);
}
result = h*(result + (functionA(a)* functionB(a) + functionA(b)* functionB(b)) / 2);
	return result;
}

vector<double> gauss(vector<vector<double>> matrix, int n)
{
    vector<double> answer;
    vector<double> wrong(1, -1000000);
    bool flag;
    for (int i = n - 1; i >= 0; i--)
    {
        double temp = 0;
        double tempAns = 0;
        int t = 0;
        flag = false;
        for (int j = n - 1; j >= i; j--)
        {
            if (i == n - 1)
            {
                if (fabs(matrix[i][j] - 0) > 1.0e-12)
                {
                    answer.push_back(matrix[i][n] / matrix[i][j]);
                    flag = true;
                }
                else
                    return wrong;
                break;
            }
            temp += matrix[i][j] * answer[t];
            t++;
            if (j == i + 1)
            {
                if ((fabs(temp - 0) < 1.0e-12 && matrix[i][n] != 0 && fabs(matrix[i][j - 1] - 0) < 1.0e-12) || (fabs(temp - 0) < 1.0e-12 && fabs(matrix[i][n] - 0) < 1.0e-12 && fabs(matrix[i][j - 1] - 0) < 1.0e-12))
                    return wrong;
                tempAns = (matrix[i][n] - temp) / matrix[i][j - 1];
                break;
            }
        }
        if (!flag)
            answer.push_back(tempAns);
    }
    return answer;
}

vector<vector<double>> triangleView(vector<vector<double>> matrix, int n)
{
    for (int i = 0; i < n - 1; i++)
    {
        double max = 0;
        int index = i;
        for (int k = i; k < n; k++)
        {
            if (fabs(matrix[k][i]) > fabs(max))
            {
                max = matrix[k][i];
                index = k;
            }
        }
        if (index != i)
        {
            vector<double> temp = matrix[i];
            matrix[i] = matrix[index];
            matrix[index] = temp;
        }
        for (int t = i + 1; t < n; t++)
        {
            double first = matrix[t][i] / matrix[i][i];
            for (int p = i; p < n + 1; p++)
            {
                matrix[t][p] -= first * matrix[i][p];
            }
        }
    }
    return matrix;
}

double integrateI(double a, double b, double p, double(*functionA)(double), double(* functionB)(double, int), int t)
{
	double result = 0;
double h = (b - a) / p;
	for (int i = 1; i<p; i++)
	{
		result += functionA(a + h* i)* functionB(a + h* i, t);
	}
	result = h*(result + (functionA(a)* functionB(a, t) + functionA(b)* functionB(b, t)) / 2);
	return result;
}

double eione(double x, int i)
{
    return (x * x - 2 * (i - 2) * h * x + pow((i - 2) * h, 2)) / (h * h);
}

double eitwo(double x, int i)
{
    return (-2 * x * x + 4 * (h + (i - 1) * h) * x + (i + 1) * h * h - 2 * (i + 1) * h * (i - 1) * h - (i - 1) * h * h) / (2 * h * h);
}

double eithree(double x, int i)
{
    return (x * x - 2 * (i + 2) * h * x + pow((i + 2) * h, 2)) / (h * h);
}

double integrateEi(int i)
{
    return integrateI((i - 2) * h, (i - 1) * h, p, fx, eione, i) + integrateI((i - 1) * h, (i + 1) * h, p, fx, eitwo, i)
        + integrateI((i + 1) * h, (i + 2) * h, p, fx, eithree, i);
}
*/
