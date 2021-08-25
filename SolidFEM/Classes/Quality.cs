using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SolidFEM.Classes
{
    class Quality
    {
        //Properties
        public double AspectRatio { get; set; }
        public double Skewness { get; set; }
        public double JacobianRatio { get; set; }
        public Element element { get; set; }


        //Constructor
        public Quality()
        {
            //Empty constructor
        }

        public Quality(Element _elem, double _aspectRatio, double _skewness, double _jacobianRatio)
        {
            AspectRatio = _aspectRatio;
            Skewness = _skewness;
            JacobianRatio = _jacobianRatio;
            element = _elem;
        }
    }
}
