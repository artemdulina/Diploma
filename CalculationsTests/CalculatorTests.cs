using System;
using System.Collections;
using NUnit.Framework;
using Calculations;

namespace CalculationsTests
{
    [TestFixture]
    public class CalculatorTests
    {
        public static int p = 10000;
        public static double delta = 0.01;

        public static IEnumerable IntegrateTestCases
        {
            get
            {
                yield return new TestCaseData();
            }
        }

        [Test]
        public void Integrate_TwoFunctionsMultiplyAllNegativeRange_ShouldReturnCorrect()
        {
            double actual = Calculator.Integrate(-5, -8, p, x => Math.Pow(x, 2), x => Math.Pow(x, 3));

            double expected = 41086.5;

            Assert.AreEqual(expected, actual, delta);
        }

        [Test]
        public void Integrate_TwoFunctionsMultiplyNegativeRange_ShouldReturnCorrect()
        {
            double actual = Calculator.Integrate(-5, 7, p, x => Math.Pow(x, 2), x => Math.Pow(x, 3));

            double expected = 17004;

            Assert.AreEqual(expected, actual, delta);
        }

        [Test]
        public void Integrate_TwoFunctionsMultiply_ShouldReturnCorrect()
        {
            double actual = Calculator.Integrate(1, 2, p, x => Math.Pow(x, 2), x => Math.Pow(x, 3));

            double expected = 10.5;

            Assert.AreEqual(expected, actual, delta);
        }


        public void FivediagonalMatrixAlgorithm()
        {
            Assert.Fail();
        }
    }
}