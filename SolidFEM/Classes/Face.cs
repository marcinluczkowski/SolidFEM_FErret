using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Grasshopper.Kernel;
using Rhino.Geometry;


namespace SolidFEM.Classes
{
    public class Face
    {
        public string Name;
        public int ID;
        public List<Node> Nodes;
    }
}