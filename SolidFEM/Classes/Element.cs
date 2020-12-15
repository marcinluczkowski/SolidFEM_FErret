using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace SolidFEM.Classes
{
    class Element
    {
        public int ID;
        public string name;
        public Mesh element_mesh;
        public List<Node> nodes;


        public Element()
        {
            //empty constructor
        }

        public static List<Node> sortNodes(Element e)
        {
            //sorting algorithm
            List<Node> sortedNodes = new List<Node>();

            List<Node> unsortedNodes = e.nodes;


            return sortedNodes;
        }
    }
}
