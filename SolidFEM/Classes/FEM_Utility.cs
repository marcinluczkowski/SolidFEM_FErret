using System;
using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry;
using CSD = CSparse.Double;

namespace SolidFEM.Classes
{
    static class FEM_Utility
    {
        public static void ElementsFromMeshList(List<Mesh> mList, List<Point3d> globalNodePts,out List<Element> femElements)
        {
            List<Element> elements = new List<Element>(); // all the elements of the mesh
            List<Node> globalNodes = new List<Node>(); // all the nodes of the mesh
            
            int IDCounter = 0; // id for each element         
            

            // iterate through all mesh elements
            foreach (var mEl in mList)
            {
                Element el = new Element();
                List<Node> elNodes = new List<Node>();
                el.ID = IDCounter; // assign the ID
                int localNodeID = 0;
                List<int> connectivity = new List<int>();

                // iterate through the vertices
                var vertices = mEl.Vertices.ToPoint3dArray(); // an array of the vertices of each mesh element
                for (int i = 0; i < vertices.Length; i++)
                {
                    Point3d mPt = vertices[i]; // point to create node from

                    // is the point already registered as a global node
                    int globalInd = globalNodePts.IndexOf(mPt); // find the global ID of the element node. -1 if not present
                    Node n;
                    if (globalInd == -1) 
                    {
                        throw new Exception("There is a mismatch between the elements' node and the nodes of the system... ");
                    }
                    // create a new node
                    n = new Node(globalInd, mPt);
                    
                    n.ID = localNodeID; // 
                    ;
                    elNodes.Add(n); // add node to list of the element's nodes
                    connectivity.Add(globalInd); // add the global index to the connectivity list

                    localNodeID++;
                }
                
                el.Nodes = elNodes; // add nodes to the elements
                el.Connectivity = connectivity; // Add the connectivity

                elements.Add(el);
                IDCounter++;
            }

            femElements = elements;

        }


        public static List<Point3d> GetMeshNodes(List<Mesh> meshList)
        {
            List<Point3d> uniquePts = new List<Point3d>();

            foreach (Mesh mesh in meshList)
            {
                Point3d[] vertices = mesh.Vertices.ToPoint3dArray();
                foreach (Point3d point3D in vertices)
                {
                    int index = uniquePts.IndexOf(point3D); // returns -1 if the point is not already added to the list
                    if (index == -1)
                    {
                        uniquePts.Add(point3D); // add the new point to the list
                    }
                }
            }
            return uniquePts;
        }
    }
}
