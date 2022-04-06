using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics;

namespace SolidFEM.Classes
{
    public class Material
    {
        public double YoungModulus { get; set; }

        public double Exx { get; set; } = 0;
        public double Eyy { get; set; } = 0;
        public double Ezz { get; set; } = 0;
        public double PossionRatio { get; set; }

        public double nu_xy { get; set; } = 0;
        public double nu_yz { get; set; } = 0;
        public double nu_zx { get; set; } = 0;

        public double Gxy { get; set; } = 0;
        public double Gyz { get; set; } = 0;
        public double Gzx { get; set; } = 0;

        public double YieldingStress { get; set; }
        public double Weight { get; set; }
        // Constructor
        public Material()
        {
            // empty constructor
        }

        public Material(double _youngModulus, double _possionRatio, double _yieldingStress)
        {
            YoungModulus = _youngModulus;
            PossionRatio = _possionRatio;
            YieldingStress = _yieldingStress;
        }

        public Material(double _Exx, double _Eyy, double _Ezz, double _nu_xy, double _nu_yz, double _nu_zx, double _Gxy, double _Gzx, double _Gyz, double _yieldingStress)
        {
            Exx = _Exx;
            Eyy = _Eyy;
            Ezz = _Ezz;
            nu_xy = _nu_xy;
            nu_yz = _nu_yz;
            nu_zx = _nu_zx;
            Gxy = _Gxy;
            Gyz = _Gyz;
            Gzx = _Gzx;
            YieldingStress = _yieldingStress;
        }

        // method
        public Matrix<double> GetMaterialConstant()
        {
            if (Exx == 0)
            {
                double possionRatio = this.PossionRatio;
                double youngModulus = this.YoungModulus;
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
            else
            {
                double Exx = this.Exx;
                double Eyy = this.Eyy; 
                double Ezz = this.Ezz;
                double nu_xy = this.nu_xy;
                double nu_yz = this.nu_yz;
                double nu_zx = this.nu_zx;
                double nu_yx = Eyy / Exx * nu_xy;
                double nu_xz = Ezz / Exx * nu_zx;
                double nu_zy = Ezz / Eyy * nu_yz;
                double G_xy = this.Gxy;
                double G_yz = this.Gyz;
                double G_zx = this.Gzx;

                //x = 1, y = 2, z = 3

                double delta = (1 - nu_xy * nu_yx - nu_xz * nu_zx - nu_yz * nu_zy - 2 * nu_yx * nu_xz * nu_zy) / (Exx * Eyy * Ezz);
                double C11 = (1 - nu_yz * nu_zy) / (Eyy * Ezz * delta);
                double C12 = (nu_xy + nu_zx * nu_zy) / (Eyy * Ezz * delta);
                double C13 = (nu_zx + nu_xy * nu_yz) / (Eyy * Ezz * delta);
                double C21 = (nu_xy + nu_zx * nu_zy) / (Exx * Ezz * delta);
                double C22 = (1 - nu_zx * nu_xz) / (Exx * Ezz * delta);
                double C23 = (nu_yz + nu_zx * nu_yx) / (Exx * Ezz * delta);
                double C31 = (nu_zx + nu_xy * nu_yz) / (Exx * Eyy * delta);
                double C32 = (nu_yz + nu_zx * nu_yz) / (Exx * Eyy * delta);
                double C33 = (1 - nu_xy * nu_yx) / (Exx * Eyy * delta);
                double C44 = 2*G_xy;
                double C55 = 2*G_zx;
                double C66 = 2*G_yz;

                Matrix<double> C = CreateMatrix.DenseOfArray(new double[,]
                {
                    {C11, C12, C13,   0,   0,   0},
                    {C21, C22, C23,   0,   0,   0},
                    {C31, C32, C33,   0,   0,   0},
                    {  0,   0,   0, C44,   0,   0},
                    {  0,   0,   0,   0, C55,   0},
                    {  0,   0,   0,   0,   0, C66},

                });
                return C;
            }
        }

    }
}
