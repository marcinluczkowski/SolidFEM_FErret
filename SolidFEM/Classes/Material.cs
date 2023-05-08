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
        public string Name { get; set; }
        public double Weight { get; set; }

        public double YieldingStress { get; set; }

        public Material()
        {
            // empty constructor
        }
        // methods
        public virtual Matrix<double> GetMaterialConstant()
        {
            return null;
        }


    }

    public class MaterialIso : Material
    {
        public double YoungModulus { get; set; }
        public double PossionRatio { get; set; }
        


        // Constructor
        public MaterialIso()
        {
            // empty constructor
        }

        public MaterialIso(double _youngModulus, double _possionRatio, double _yieldingStress)
        {
            YoungModulus = _youngModulus;
            PossionRatio = _possionRatio;
            YieldingStress = _yieldingStress;
        }

        public override Matrix<double> GetMaterialConstant()
        {
            double possionRatio = this.PossionRatio;
            double youngModulus = this.YoungModulus;
            Matrix<double> C = CreateMatrix.DenseOfArray(new double[,]
            {
                {1-possionRatio, possionRatio, possionRatio, 0, 0, 0},
                {possionRatio, 1-possionRatio, possionRatio, 0, 0, 0},
                {possionRatio, possionRatio, 1-possionRatio, 0, 0, 0},
                {0, 0, 0, (1-2*possionRatio)/(double)2, 0, 0},
                {0, 0, 0, 0, (1-2*possionRatio)/(double)2, 0},
                {0, 0, 0, 0, 0, (1-2*possionRatio)/(double)2},
            });
            C = C.Multiply((double)youngModulus / (double)((1 + possionRatio) * (1 - 2 * possionRatio)));
            return C;
        }
    }

    public class MaterialOrto : Material
    {
        public double E1 { get; set; }
        public double E2 { get; set; }
        public double E3 { get; set; }
        public double v12 { get; set; }
        public double v13 { get; set; }
        public double v23 { get; set; }
        public double G12 { get; set; }
        public double G13 { get; set; }
        public double G23 { get; set; }

        
        // Constructor
        public MaterialOrto()
        {
            // empty constructor
        }

        public MaterialOrto(double _E1, double _E2, double _E3, double _v12, double _v13, double _v23, double _G12, double _G13, double _G23, double _yieldingStress)
        {
            E1 = _E1;
            E2 = _E2;
            E3 = _E3;
            v12 = _v12;
            v13 = _v13;
            v23 = _v23;
            G12 = _G12;
            G13 = _G13;
            G23 = _G23;


            YieldingStress = _yieldingStress;
        }

        // method
        public override Matrix<double> GetMaterialConstant()
        {
            double E1 = this.E1;
            double E2 = this.E2;
            double E3 = this.E3;

            double v12 = this.v12;
            double v13 = this.v13;
            double v23 = this.v23;

            double G12 = this.G12;
            double G13 = this.G13;
            double G23 = this.G23;


            double v21 = v12 * E2 / E1;
            double v31 = v13 * E3 / E1;
            double v32 = v23 * E3 / E2;

            double D = (1 / (E1 * E2 * E3)) * (1 - 2 * v21 * v13 * v32 - v23 * v32 - v12 * v21 - v13 * v31);

            Matrix<double> C = CreateMatrix.DenseOfArray(new double[,]
            {
                {(1-v23*v32)/(E2*E3*D), (v21 + v31*v23)/(E2*E3*D), (v31+v21*v32)/(E2*E3*D), 0, 0, 0},
                {(v21 + v31*v23)/(E2*E3*D), (1-v13*v31)/(E1*E3*D), (v32+v12*v31)/(E1*E3*D), 0, 0, 0},
                {(v31+v21*v32)/(E2*E3*D), (v32+v12*v31)/(E1*E3*D), (1-v12*v21)/(E1*E2*D), 0, 0, 0},
                {0, 0, 0, G23, 0, 0},
                {0, 0, 0, 0, G13, 0},
                {0, 0, 0, 0, 0, G12},
            });
            return C;
        }

    }


    //TODO Join this and Isotrop material through an abstract parent class. 
}
