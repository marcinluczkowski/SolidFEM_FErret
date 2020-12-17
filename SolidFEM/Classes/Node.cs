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

        public Node(int _id, string _name, Point3d _point)
        {
            ID = _id;
            name = _name;
            point = _point;
        }

    }
}
