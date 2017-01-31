using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Calculations;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace Diploma
{
    public class Program
    {
        public static int n = 10;
        public static double h = (double)1 / n;
        public static double p = 1000;

        public static void PrintMatrix<T>(T[][] matrix)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix[i].Length; j++)
                {
                    Console.Write(matrix[i][j] + " ");
                }
                Console.WriteLine();
            }
        }

        public static string LeftEquation = "u''";
        public static Func<double, double> RightFunction = value => Math.Pow(value, 2) - 1;

        static void Main(string[] args)
        {
            double[,] equation = Calculator.GetSystemOfEquations(LeftEquation, RightFunction);

            double[] rightValues = Calculator.GetRightValues();

            var left = Matrix<double>.Build.DenseOfArray(equation);
            var right = Vector<double>.Build.Dense(results);
            Console.WriteLine(left.Solve(right));


            /*
            double[][] matrix = new double[n + 1][];
            for (int i = 0; i < n + 1; i++)
            {
                matrix[i] = new double[n + 2];
            }

            matrix[0][0] = 2 * (1 + 1 / (h * h));
            matrix[0][1] = -4 / (h * h);
            matrix[0][2] = 2 / (h * h);
            matrix[0][n + 1] = 10 + integrate(0, h, p, fx, eZeroone) + integrate(h, 2 * h, p, fx, eOneone);
            matrix[1][0] = 1.5 + 1 / (h * h);
            matrix[1][1] = -3 / (2 * h * h);
            matrix[1][3] = 1 / (2 * h * h);
            matrix[1][n + 1] = 7.5 + integrate(0, h, p, fx, eZerotwo) + integrate(h, 3 * h, p, fx, eOnetwo);
            for (int i = 2; i < n - 1; i++)
            {
                matrix[i][i - 2] = -2 / (h * h);
                matrix[i][i - 1] = 4 / (h * h);
                matrix[i][i + 1] = -4 / (h * h);
                matrix[i][i + 2] = 2 / (h * h);
                matrix[i][n + 1] = integrateEi(i);
            }
            matrix[n - 1][n - 3] = -1 / (2 * h * h);
            matrix[n - 1][n - 1] = 1 / (2 * h * h);
            matrix[n - 1][n] = -2;
            matrix[n - 1][n + 1] = -2 + integrate((n - 3) * h, (n - 1) * h, p, fx, eNminusOneone) + integrate((n - 1) * h, n * h, p, fx, eNone);
            matrix[n][n - 2] = -2 / (h * h);
            matrix[n][n - 1] = 2 / (h * h);
            matrix[n][n] = -3;
            matrix[n][n + 1] = -3 + integrate((n - 2) * h, (n - 1) * h, p, fx, eNminusOnetwo) + integrate((n - 1) * h, n * h, p, fx, eNtwo);

            Console.WriteLine("Exact solutions: ");
            for (int i = 0; i <= n; i++)
            {
                Console.WriteLine(u(h * i));
            }
            Console.WriteLine("\nApproximate solutions found by the method of Gauss:");

            vector<vector<double>> triangle = triangleView(matrix, n + 1);
            vector<double> answers = gauss(triangle, n + 1);

            for (int i = n; i >= 0; i--)
            {
                Console.Write("answers[i] ");
            }
            Console.WriteLine("\nApproximate solutions found by a five-point sweep:");

            vector<vector<double>> matrixf(n + 1);
            vector<double> b;
            matrixf[0].push_back(2 * (1 + 1 / (h * h)));
            matrixf[0].push_back(-4 / (h * h));
            matrixf[0].push_back(2 / (h * h));
            b.push_back(10 + integrate(0, h, p, fx, eZeroone) + integrate(h, 2 * h, p, fx, eOneone));
            matrixf[1].push_back(1.5 + 1 / (h * h));
            matrixf[1].push_back(-3 / (2 * h * h));
            matrixf[1].push_back(0);
            matrixf[1].push_back(1 / (2 * h * h));
            b.push_back(7.5 + integrate(0, h, p, fx, eZerotwo) + integrate(h, 3 * h, p, fx, eOnetwo));
            for (int i = 2; i < n - 1; i++)
            {
                matrixf[i].push_back(-2 / (h * h));
                matrixf[i].push_back(4 / (h * h));
                matrixf[i].push_back(0);
                matrixf[i].push_back(-4 / (h * h));
                matrixf[i].push_back(2 / (h * h));
                b.push_back(integrateEi(i));
            }
            matrixf[n - 1].push_back(-1 / (2 * h * h));
            matrixf[n - 1].push_back(0);
            matrixf[n - 1].push_back(1 / (2 * h * h));
            matrixf[n - 1].push_back(-2);
            b.push_back(-2 + integrate((n - 3) * h, (n - 1) * h, p, fx, eNminusOneone) + integrate((n - 1) * h, n * h, p, fx, eNone));
            matrixf[n].push_back(-2 / (h * h));
            matrixf[n].push_back(2 / (h * h));
            matrixf[n].push_back(-3);
            b.push_back(-3 + integrate((n - 2) * h, (n - 1) * h, p, fx, eNminusOnetwo) + integrate((n - 1) * h, n * h, p, fx, eNtwo));

            vector<double> answersf = fivediagonalMatrixAlgorithm(matrixf, b, n);
            for (int i = 0; i < n + 1; i++)
            {
                Console.WriteLine(answersf[i]);
            }
            */

        }
    }
}
