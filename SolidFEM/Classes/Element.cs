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

        public static List<Node> sortOrthogonalNodes(Element e)
        {
            //sorting algorithm
            List<Node> sortedNodes = new List<Node>();

            List<Node> unsortedNodes = e.nodes;

            // Calculate the center point
            double sumX = 0;
            double sumY = 0;
            double sumZ = 0; 
            foreach(Node n in unsortedNodes)
            {
                sumX += n.point.X;
                sumY += n.point.Y;
                sumZ += n.point.Z;
            }
            Point3d centerPt = new Point3d( sumX/ unsortedNodes.Count, sumY / unsortedNodes.Count, sumZ / unsortedNodes.Count) ;

            // If points are below centerPt,
            var bottomNodes = new List<Node>();
            var topNodes = new List<Node>();

            // Assign the nodes in top and bottom list
            foreach(Node n in unsortedNodes)
            {
                if (n.point.Z > centerPt.Z)
                {
                    topNodes.Add(n);
                }
                else
                    bottomNodes.Add(n);
            }

            //Sort bottom nodes

            // Calculate the center point
            double sumX_bottom = 0;
            double sumY_bottom = 0;
            double sumZ_bottom = 0;
            foreach (Node n in bottomNodes)
            {
                sumX_bottom += n.point.X;
                sumY += n.point.Y;
                sumZ += n.point.Z;
            }
            Point3d bottomCenter = new Point3d(sumX_bottom / unsortedNodes.Count, sumY_bottom / unsortedNodes.Count, sumZ_bottom / unsortedNodes.Count);

            //Sort top nodes
            double sumX_top = 0;
            double sumY_top = 0;
            double sumZ_top = 0;
            foreach (Node n in topNodes)
            {
                sumX_bottom += n.point.X;
                sumY += n.point.Y;
                sumZ += n.point.Z;
            }
            Point3d topCenter = new Point3d(sumX_bottom / unsortedNodes.Count, sumY_bottom / unsortedNodes.Count, sumZ_bottom / unsortedNodes.Count);



            Node[] bottom = bottomNodes.ToArray();

            double[] angles = new double[bottom.Length];

            for (int i = 0; i < bottom.Length; i++)

                angles[i] = Math.Atan2(bottom[i].point.Y, bottom[i].point.X);

            Array.Sort(angles, bottom); // Sorting the bottom nodes



            Node[] top = bottomNodes.ToArray();

            double[] anglesTop = new double[bottom.Length];

            for (int i = 0; i < top.Length; i++)

                angles[i] = Math.Atan2(top[i].point.Y, top[i].point.X);

            Array.Sort(angles, top); // Sorting the top nodes. 

            // Join the bottom and top nodes. 
            sortedNodes.AddRange(bottom.ToList());
            sortedNodes.AddRange(top.ToList());
            return sortedNodes;
        }
    }
}
