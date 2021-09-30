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

        public static CSD.DenseMatrix GetBodyForceVector(Material material, List<Element> elements, int numGlobalNodes)
        {
            // Initiate the empty body force vector: 
            CSD.DenseMatrix F_body = new CSD.DenseMatrix(numGlobalNodes * 3, 1);



            for (int i = 0; i < elements.Count; i++)
            {
                // get the current element
                Element el = elements[i];

                // first, get the global coordinates
                CSD.DenseMatrix globalElementCoordinate = new CSD.DenseMatrix(el.Nodes.Count, 3); // one column for x,y, and z coordinate
                for (int j = 0; j < el.Nodes.Count; j++)
                {
                    Node n = el.Nodes[j];
                    globalElementCoordinate[j, 0] = n.Coordinate.X;
                    globalElementCoordinate[j, 1] = n.Coordinate.Y;
                    globalElementCoordinate[j, 2] = n.Coordinate.Z;
                }

                // get the matrix of natural coordinates in gauss points
                var gaussCoordinates = GetGaussPointMatrix(); // by defaul we have a 2x2x2 integration of Hex8 element

                // element force vector
                CSD.DenseMatrix elForceVec = new CSD.DenseMatrix(el.Nodes.Count, 1);

                // iterate through each node
                for (int j = 0; j < gaussCoordinates.RowCount; j++)
                {
                    double r = gaussCoordinates[j, 0];
                    double s = gaussCoordinates[j, 1];
                    double t = gaussCoordinates[j, 2];

                    // get the partial derivatives evaluated at a gauss point
                    var partialDerivatives = PartialDerivateShapeFunctions(r, s, t, "Hex8");

                    // get the jacobian

                    // get the H matrix from the shape functions
                }


                // calculate the elements Jacobian. 
                         

            }

            // I need the shape functions for the element. 8 noded element mean eight nodes. Should ensure that this is applicable for other types of elements as well. 


            return F_body;
        }

        /// <summary>
        /// Returns the natural coordinate of Gauss points in a matrix
        /// </summary>
        /// <param name="numPoints">e.g. 2 for 2x2x2x gauss points</param>
        /// <param name="elType"> Which element to work on use "Hex8" for now. </param>
        /// <returns></returns>
        public static CSD.DenseMatrix GetGaussPointMatrix(int numPoints = 2, string elType = "Hex8")
        {
            // for now we only have Hex8 element implementation
            if (elType == "Hex8")
            {
                if(numPoints == 2)
                {
                    double gaussPoint = 0.57735; // Numerical value of the +- natural coordinate to use in gauss point integration for 2 sampling point. From Bathes book Table 5.6
                    double[] gaussArray = new double[]
                    {
                        -gaussPoint, -gaussPoint, -gaussPoint ,
                        gaussPoint, -gaussPoint, -gaussPoint,
                        gaussPoint, gaussPoint, -gaussPoint,
                        -gaussPoint, gaussPoint, -gaussPoint,
                        -gaussPoint, -gaussPoint, gaussPoint,
                        gaussPoint, -gaussPoint, gaussPoint,
                        gaussPoint, gaussPoint, gaussPoint,
                        -gaussPoint, gaussPoint, gaussPoint

                    };
                    var naturalCoordinatesGauss = new CSD.DenseMatrix(8, 3, gaussArray);
                    return naturalCoordinatesGauss;

                }
                else
                {
                    throw new NotImplementedException("The integration type is not yet implemented.");
                }

            }
            else
            {
                throw new NotImplementedException("This method is not yet implemented");
            }


            return null;
        }

        /// <summary>
        /// Return a vector of a shape function evaluated in the r, s, and t natural coordinates.
        /// </summary>
        /// <param name="r"></param>
        /// <param name="s"></param>
        /// <param name="t"></param>
        /// <param name="elType"></param>
        /// <returns></returns>
        public static CSparse.Storage.DenseColumnMajorStorage<double> GetShapeFunctions(double r, double s, double t ,string elType = "Hex8")
        {
            if (elType == "Hex8")
            {
                double N1 = 0.125 * (1 - r) * (1 - s) * (1- t);
                double N2 = 0.125 * (1 + r) * (1 - s) * (1 - t);
                double N3 = 0.125 * (1 + r) * (1 + s) * (1 - t);
                double N4 = 0.125 * (1 - r) * (1 + s) * (1 - t);
                double N5 = 0.125 * (1 - r) * (1 - s) * (1 + t);
                double N6 = 0.125 * (1 + r) * (1 - s) * (1 + t);
                double N7 = 0.125 * (1 + r) * (1 + s) * (1 + t);
                double N8 = 0.125 * (1 - r) * (1 + s) * (1 + t);

                double[] shapeArray;
            }

        }

        /// <summary>
        /// Returns a matrix of the partial derivatives of the shapefunctions at the input 
        /// natural coordinates r, s, and t.
        /// </summary>
        /// <param name="elType"></param>
        /// <param name="r"></param>
        /// <param name="s"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public static CSparse.Storage.DenseColumnMajorStorage<double> PartialDerivateShapeFunctions(double r, double s, double t, string elType = "Hex8")
        {
            if (elType == "Hex8")
            {
                double c = 0.125;
                double[,] derivateArray = new double[,]
                {
                    {-(1-s)*(1-t)*c, (1-s)*(1-t)*c, (1+s)*(1-t)*c,-(1+s)*(1-t)*c, -(1-s)*(1+t)*c, (1-s)*(1+t)*c,(1+s)*(1+t)*c,-(1+s)*(1+t)*c},
                    {-(1-r)*(1-t)*c, -(1+r)*(1-t)*c, (1+r)*(1-t)*c,(1-r)*(1-t)*c,-(1-r)*(1+t)*c,-(1+r)*(1+t)*c,(1+r)*(1+t)*c,(1-r)*(1+t)*c},
                    {-(1-r)*(1-s)*c, -(1+r)*(1-s)*c, -(1+r)*(1+s)*c,-(1-r)*(1+s)*c,(1-r)*(1-s)*c,(1+r)*(1-s)*c,(1+r)*(1+s)*c,(1-r)*(1+s)*c}
                };
                var derivativeMatrix = CSD.DenseMatrix.OfArray(derivateArray);
                //derivativeMatrix = derivativeMat
                return derivativeMatrix; 

            }
            else
            {
                throw new NotImplementedException("The selected element type is not implemented.");
            }
        }
    }
}
