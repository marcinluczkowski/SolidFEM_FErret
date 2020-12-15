using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace SolidFEM.Classes
{
    class Node
    {
        public int ID;
        public string name;
        public Point3d point;

        public Node()
        {
            //empty constructor
        }
        
    }
}
