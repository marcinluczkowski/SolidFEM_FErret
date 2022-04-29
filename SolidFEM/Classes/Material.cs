using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace SolidFEM.Classes
{
    public class Material
    {
        public double E { get; set; }
        public double PossionRatio { get; set; }
        public double YieldingStress { get; set; }
        public double Weight { get; set; }
        // Constructor
        public Material()
        {
            // empty constructor
        }

        public Material(double _E, double _possionRatio, double _yieldingStress)
        {
            E = _E;
            PossionRatio = _possionRatio;
            YieldingStress = _yieldingStress;
        }

        // method
        public Matrix<double> GetMaterialConstant()
        {
            double possionRatio = this.PossionRatio;
            double E = this.E;

                double youngModulus = E;
                Matrix<double> C = CreateMatrix.DenseOfArray(new double[,]
                            {
                {1-possionRatio, possionRatio, possionRatio, 0, 0, 0},
                {possionRatio, 1-possionRatio, possionRatio, 0, 0, 0},
                {possionRatio, possionRatio, 1- possionRatio, 0, 0, 0},
                {0, 0, 0, (1-2*possionRatio)/(double)2, 0, 0},
                {0, 0, 0, 0, (1-2*possionRatio)/(double)2, 0},
                {0, 0, 0, 0, 0, (1-2*possionRatio)/(double)2},
                            });
                C = C.Multiply((double)youngModulus / (double)((1 + possionRatio) * (1 - 2 * possionRatio)));
                return C;
            
        }

    }
}
