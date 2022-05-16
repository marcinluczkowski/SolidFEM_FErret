using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace SolidFEM.Classes
{
    public static class FEM
    {
            public static Vector<double> GetShapeFunctions(double r, double s, double t, string elType)
            {
                // Shapefunctions on matrix form
                Vector<double> shapeFunctionsValues = DenseVector.OfArray(new double[] { });
                if (elType == "Hex8")
                {
                    Vector<double> N = DenseVector.OfArray(new double[]
                    {
                    (1-r)*(1-s)*(1-t),
                    (1+r)*(1-s)*(1-t),
                    (1+r)*(1+s)*(1-t),
                    (1-r)*(1+s)*(1-t),
                    (1-r)*(1-s)*(1+t),
                    (1+r)*(1-s)*(1+t),
                    (1+r)*(1+s)*(1+t),
                    (1-r)*(1+s)*(1+t)
                    });
                    N = N.Multiply(0.125);
                    shapeFunctionsValues = N;
                }
                else if (elType == "Tet10")
                {
                    double N1 = 2 * r * r + 4 * r * s + 4 * r * t - 3 * r + 2 * s * s + 4 * s * t - 3 * s + 2 * t * t - 3 * t + 1;
                    double N2 = 2 * r * r - r;
                    double N3 = 2 * s * s - s;
                    double N4 = 2 * t * t - t;
                    double N5 = -4 * r * r - 4 * r * s - 4 * r * t + 4 * r;
                    double N6 = 4 * r * s;
                    double N7 = -4 * r * s - 4 * s * s - 4 * s * t + 4 * s;
                    double N8 = -4 * r * t - 4 * s * t - 4 * t * t + 4 * t;
                    double N9 = 4 * r * t;
                    double N10 = 4 * s * t;

                    Vector<double> N = DenseVector.OfArray(new[] {N1, N2, N3, N4});
                    shapeFunctionsValues = N;
                }
                else if (elType == "Tet4")
                {
                    double N1 = 1 - r - s - t;
                    double N2 = r;
                    double N3 = s;
                    double N4 = t;

                    Vector<double> N = DenseVector.OfArray(new[] { N1, N2, N3, N4 });
                    shapeFunctionsValues = N;
            }

                return shapeFunctionsValues;
            }



            public static Matrix<double> DerivateWithNatrualCoordinates(double r, double s, double t, int nodeDOFS)
            {
                if (nodeDOFS == 2)
                {
                    Matrix<double> shapeFunctionsDerivatedNatural = DenseMatrix.OfArray(new double[,]
                    {
                    {-(1-s), (1-s), (1+s), -(1+s)},
                    {-(1-r), -(1+r), (1+r), (1-r)}
                    });
                    shapeFunctionsDerivatedNatural = shapeFunctionsDerivatedNatural.Multiply(0.25);
                    return shapeFunctionsDerivatedNatural;

                }
                else
                {
                    Matrix<double> shapeFunctionsDerivatedNatural = DenseMatrix.OfArray(new double[,]
                    {
                        {-(1-s)*(1-t), (1-s)*(1-t), (1+s)*(1-t),-(1+s)*(1-t),-(1-s)*(1+t),(1-s)*(1+t),(1+s)*(1+t),-(1+s)*(1+t)},
                        {-(1-r)*(1-t), -(1+r)*(1-t), (1+r)*(1-t),(1-r)*(1-t),-(1-r)*(1+t),-(1+r)*(1+t),(1+r)*(1+t),(1-r)*(1+t)},
                        {-(1-r)*(1-s), -(1+r)*(1-s), -(1+r)*(1+s),-(1-r)*(1+s),(1-r)*(1-s),(1+r)*(1-s),(1+r)*(1+s),(1-r)*(1+s)},
                    });
                    shapeFunctionsDerivatedNatural = shapeFunctionsDerivatedNatural.Multiply(0.125);
                    return shapeFunctionsDerivatedNatural;
                }

            }
            public static Matrix<double> GetNaturalCoordinate(double scale)  // NOT IN USE, 
            {
                Matrix<double> nodes = DenseMatrix.OfArray(new double[,]{});

                
                    double gp = scale;
                    
                    Matrix<double> natNodes = DenseMatrix.OfArray(new double[,]
                    {
                        {-gp, -gp, -gp},
                        {gp, -gp, -gp},
                        {gp, gp, -gp},
                        {-gp, gp, -gp},
                        {-gp, -gp, gp},
                        {gp, -gp, gp},
                        {gp, gp, gp},
                        {-gp, gp, gp},
                    });
                    nodes = natNodes;

                    return nodes;
            }
        }
    
}
