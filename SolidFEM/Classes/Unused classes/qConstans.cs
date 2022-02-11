using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SolidFEM.Classes
{
    class qConstants
    {
        public qConstants()
        {
            //Empty constructor
        }

        public double GetThetaTolerance()
        {
            return Math.PI / (double)6;
        }
        public double GetThetaToleranceForClosing()
        {
            return 1.5 * this.GetThetaTolerance();
        }

        public double GetTransitionTolerance()
        {
            return 2.5;
        }

        public double GetEpsilon1()
        {
            return 0.04 * Math.PI;
        }
        public double GetEpsilon2()
        {
            return 0.09 * Math.PI;
        }

    }
}
