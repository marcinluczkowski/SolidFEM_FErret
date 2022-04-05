using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace SolidFEM.Classes
{

    public class Element
    {
        public int ID;
        public string name;
        public List<Node> Nodes { get; set; }
        public List<int> Connectivity { get; set; }
        public string Type { get; set; }
        //public Quality MeshQuality { get; set; }
        public int Id { get; set; }
        public Mesh ElementMesh { get; set; }       
        public List<Point3d> TopologyVertices
        {
            get
            {
                List<Point3d> vertices = new List<Point3d>();
                foreach (Node n in Nodes)
                {
                    vertices.Add(n.Coordinate);
                }
                return vertices;
            }
        }      
        public Matrix<double> localK; //local stiffness matrix
        public List<Matrix<double>> LocalB { get; set; } // Local B matrix for the element evaluated in all Gauss points


        // -- constructors ---
        public Element()
        {
            //empty constructor
        }
        public Element(Element element)
        {
            Id = element.Id;
            name = element.name;
            Nodes = element.Nodes;
            Connectivity = element.Connectivity;
            Type = element.Type;
            //MeshQuality = element.MeshQuality;
            ID = element.ID;
            ElementMesh = element.ElementMesh;
        }
        public Element(List<Node> nodeList, int id)
        {
            ID = id;
            Nodes = nodeList;
            name = "Element: " + ID.ToString();
        }
        public Element(int _id, List<Node> _nodes, List<int> _connectivity)
        {
            Id = _id;
            Nodes = _nodes;
            Connectivity = _connectivity;
            GetElementType();
            GetElementMesh();
        }
        //Sort the vertices of the 
        
        // -- methods ---

        public void SortVerticesByGrahamScan()
        {
            List<Point3d> vertices = TopologyVertices; // Get all the topology vertices. 
            //Subtract top and bottom vertices
            #region Split into top and bottom vertices
            List<Point3d> sortedNodes = new List<Point3d>();
            // Calculate the center point
            double sumX = 0;
            double sumY = 0;
            double sumZ = 0;
            foreach (Point3d pt in vertices)
            {
                sumX += pt.X;
                sumY += pt.Y;
                sumZ += pt.Z;
            }
            Point3d centerPt = new Point3d(sumX / vertices.Count, sumY / vertices.Count, sumZ / vertices.Count);
            // If points are below centerPt,
            var bottomNodes = new List<Point3d>();
            var topNodes = new List<Point3d>();
            // Assign the Nodes in top and bottom list
            foreach (Point3d pt in vertices)
            {
                if (pt.Z > centerPt.Z)
                {
                    topNodes.Add(pt);
                }
                else
                    bottomNodes.Add(pt);
            }
            #endregion
            #region Sort top and bottom nodes
            topNodes = topNodes.OrderBy(pt => pt.Y).ThenBy(pt => pt.X).ToList(); // Is it working if several points has the same y-value?
            topNodes = topNodes.OrderBy(pt => Math.Atan2(pt.Y - topNodes[0].Y, pt.X - topNodes[0].X)).ToList();
            List<Point3d> sortedTop = new List<Point3d>();
            while (topNodes.Count > 0)
            {
                GrahamScan(ref topNodes, ref sortedTop);
            }
            bottomNodes = bottomNodes.OrderBy(pt => pt.Y).ThenBy(pt => pt.X).ToList();
            bottomNodes = bottomNodes.OrderBy(pt => Math.Atan2(pt.Y - bottomNodes[0].Y, pt.X - bottomNodes[0].X)).ToList();
            List<Point3d> sortedBottom = new List<Point3d>();
            while (bottomNodes.Count > 0)
            {
                GrahamScan(ref bottomNodes, ref sortedBottom);
            }
            #endregion
            #region Modify element vertices with new list.
            List<Point3d> sortedVertices = new List<Point3d>();
            sortedVertices.AddRange(sortedBottom);
            sortedVertices.AddRange(sortedTop);

            for (int i = 0; i < Nodes.Count; i++)
            {
                Nodes[i].Coordinate = sortedVertices[i];
                Nodes[i].Coordinate = sortedVertices[i];
            }
            #endregion

        }
        public void GetElementType()
        {
            string type = "null";
            switch (Nodes.Count)
            {
                case 3:
                    type = "Triangle";
                    break;
                case 4:
                    type = "Quad";
                    break;
                case 6:
                    type = "Tet";
                    break;
                case 8:
                    type = "Hex";
                    break;
            }
            Type = type;
        }

        public List<List<Node>> GetFaces()
        {
            List<List<Node>> Faces = new List<List<Node>>();
            List<Node> nodes = Nodes;

            if (Type == "Quad" | Nodes.Count == 3) // surface
            {
                Faces.Add(nodes);
            }
            else // solid
            {
                List<List<int>> nodeIndex = new List<List<int>>
                {
                      new List<int> {0, 1, 5, 4}, new List<int> {1, 2, 6, 5}, new List<int> {2, 3, 7, 6},
                      new List<int> { 3, 0, 4, 7 }, new List<int> { 0, 3, 2, 1 },  new List<int> { 4, 5, 6, 7 }
                };

                foreach (List<int> indices in nodeIndex)
                {
                    List<Node> nodesOfFace = new List<Node>();
                    foreach (int i in indices)
                    {
                        nodesOfFace.Add(nodes[i]);
                    }
                    Faces.Add(nodesOfFace);
                }
            }
            return Faces;
        }

        public void GetElementMesh()
        {
            Mesh mesh = new Mesh();
            foreach (Node node in Nodes)
            {
                mesh.Vertices.Add(node.Coordinate);
            }
            if (Type == "Quad")
            {
                mesh.Faces.AddFace(0, 1, 2, 3);
            }
            else if (Type == "Hex")
            {
                mesh.Faces.AddFace(0, 1, 5, 4);
                mesh.Faces.AddFace(1, 2, 6, 5);
                mesh.Faces.AddFace(2, 3, 7, 6);
                mesh.Faces.AddFace(3, 0, 4, 7);
                mesh.Faces.AddFace(0, 1, 2, 3);
                mesh.Faces.AddFace(4, 5, 6, 7);
            }
            else if (Type == "Triangle")
            {
                mesh.Faces.AddFace(0, 1, 2);
            }
            else if (Type == "Tet")
            {
                mesh.Faces.AddFace(0, 1, 2);
                mesh.Faces.AddFace(3, 4, 5);
                mesh.Faces.AddFace(0, 1, 4, 3);
                mesh.Faces.AddFace(1, 2, 5, 4);
                mesh.Faces.AddFace(2, 0, 3, 5);
            }

            mesh.Compact(); //to ensure that it calculate
            mesh.FaceNormals.ComputeFaceNormals();
            mesh.UnifyNormals();
            ElementMesh = mesh;
        }

        private void GrahamScan(ref List<Point3d> pts, ref List<Point3d> selPts)
        {
            if (pts.Count > 0)
            {
                var pt = pts[0];
                if (selPts.Count <= 1)
                {
                    selPts.Add(pt);
                    pts.RemoveAt(0);
                }
                else
                {
                    var pt1 = selPts[selPts.Count - 1];
                    var pt2 = selPts[selPts.Count - 2];
                    Vector3d dir1 = pt1 - pt2;
                    Vector3d dir2 = pt - pt1;
                    var crossProd = Vector3d.CrossProduct(dir1, dir2);
                    if (crossProd.Z < 0) // Check if the turn is clockwise or counter-clockwise
                    {
                        selPts.RemoveAt(selPts.Count - 1); //
                    }
                    else
                    {
                        selPts.Add(pt);
                        pts.RemoveAt(0);
                    }
                }
            }
        }
    }
}
