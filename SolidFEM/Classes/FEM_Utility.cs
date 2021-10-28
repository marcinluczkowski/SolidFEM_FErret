using System;
using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry;
using CSD = CSparse.Double;
using LA = MathNet.Numerics.LinearAlgebra;
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

                // sort the vertices by a Graham 
                //List<Point3d> sortedVertices = SortVerticesByGrahamScan(vertices.ToList());
               

                for (int i = 0; i < vertices.Length; i++)
                {
                    Point3d mPt = vertices[i]; // point to create node from

                    // is the point already registered as a global node
                    int globalInd = -1;

                    
                    for (int j = 0; j < globalNodePts.Count; j++)
                    {
                        double dist = globalNodePts[j].DistanceToSquared(mPt);
                        if (dist < 0.01) // random tolerance
                        {
                            globalInd = j;
                        }
                    }
                    
                    //int globalInd = globalNodePts.IndexOf(mPt); // find the global ID of the element node. -1 if not present
                    if (globalInd == -1) 
                    {
                        throw new Exception("There is a mismatch between the elements' node and the nodes of the system... ");
                    }
                    // create a new node
                    //n = new Node(globalInd, mPt);
                    Node n = new Node(globalInd, globalNodePts[globalInd]);                    
                    n.ID = localNodeID; // 
                    
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

            // add the first points

            var firstPt = meshList[0].Vertices.ToPoint3dArray()[0];
            uniquePts.Add(new Point3d(Math.Round(firstPt.X, 4), Math.Round(firstPt.Y, 4), Math.Round(firstPt.Z, 4))); 
            foreach (Mesh mesh in meshList)
            {
                Point3d[] vertices = mesh.Vertices.ToPoint3dArray();
                foreach (Point3d point3D in vertices)
                {
                    double dist = uniquePts.Min(pt=> pt.DistanceToSquared(point3D));
                    if(dist > 0.00001)
                    {
                        //uniquePts.Add(point3D); // if the minimum distance between the current point and all other unique points are greater than a tolerance, it is not already in the list
                        uniquePts.Add(new Point3d(Math.Round(point3D.X, 4), Math.Round(point3D.Y, 4), Math.Round(point3D.Z, 4)));
                    }
                    /*
                    int index = uniquePts.IndexOf(point3D); // returns -1 if the point is not already added to the list
                    if (index == -1)
                    {
                        uniquePts.Add(point3D); // add the new point to the list
                    }*/
                }
            }
            return uniquePts;
        }
        public static List<Point3d> LocalCartesianCoordinates(Element el)
        {
            List<Point3d> localCoord = new List<Point3d>();
            List<Node> nodes = el.Nodes;
            Point3d newOrigin = nodes[0].Coordinate; // this is the new local origin of the element. 
            localCoord.Add(new Point3d(0, 0, 0)); // set the first point to 0,0,0 in the "local system"
            for (int i = 1; i < nodes.Count; i++)
            {
                double newX = nodes[i].Coordinate.X - newOrigin.X;
                double newY = nodes[i].Coordinate.Y - newOrigin.Y;
                double newZ = nodes[i].Coordinate.Z - newOrigin.Z;

                localCoord.Add(new Point3d(newX, newY, newZ));
            }

            return localCoord;
        }

        public static LA.Vector<double> GetBodyForceVector(Material material, List<Element> elements, int numGlobalNodes)
        {
            // Initiate the empty body force vector: 
            //CSD.DenseMatrix F_body = new CSD.DenseMatrix(numGlobalNodes * 3, 1);
            LA.Vector<double> F_body = LA.Vector<double>.Build.Dense(numGlobalNodes * 3); // create the empty load vector
            LA.Vector<double> bodyLoadVector = LA.Vector<double>.Build.DenseOfArray(new double[] {0, 0, -  (material.Weight * 9.81 * Math.Pow(10, -9) ) });
            double elementJacobianTest = 0;
            for (int i = 0; i < elements.Count; i++)
            {
                // get the current element
                Element el = elements[i];

                // first, get the global coordinates
                LA.Matrix<double> globalElementCoordinate = LA.Matrix<double>.Build.Dense(el.Nodes.Count, 3);
                //CSD.DenseMatrix globalElementCoordinate = new CSD.DenseMatrix(el.Nodes.Count, 3); // one column for x,y, and z coordinate
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
                LA.Vector<double> elForceVec = LA.Vector<double>.Build.Dense(el.Nodes.Count * 3);

                // iterate through each gauss node
                for (int j = 0; j < gaussCoordinates.RowCount; j++)
                {
                    double r = gaussCoordinates[j, 0];
                    double s = gaussCoordinates[j, 1];
                    double t = gaussCoordinates[j, 2];

                    // get the partial derivatives evaluated at a gauss point
                    var partialDerivatives = PartialDerivateShapeFunctions(r, s, t, "Hex8");

                    // get the jacobian
                    var jacobianOperator = partialDerivatives.Multiply(globalElementCoordinate);
                    double jacobianDeterminant = jacobianOperator.Determinant();
                    elementJacobianTest += jacobianDeterminant;
                    //double jacobianDet = jacobianOperator. // how to calculate the  determinant????
                    // get the H matrix from the shape functions
                    var shapeFunctions = GetShapeFunctions(r, s, t, "Hex8");

                    // the H-matrix is the displacement interpolation matrix. 
                    var interpolationMatrix = DisplacementInterpolationMatrix(shapeFunctions, 3);
                    double alpha_ijk = 1.0; // should probably be a vector in a complete solver
                    double t_ijk = 1.0; // not sure what this is yet
                    // with a 2x2x2 integration scheme the integration constant is 1.0 for all gauss points. 
                    var gaussPointLoadVector = interpolationMatrix.TransposeThisAndMultiply(bodyLoadVector);
                    gaussPointLoadVector = gaussPointLoadVector.Multiply( jacobianDeterminant * alpha_ijk * t_ijk );


                    // add the vector from the gauss point to the element load vector
                    elForceVec.Add(gaussPointLoadVector, elForceVec);
                    //elForceVec = elForceVec.Add(gaussPointLoadVector, elForceVec);
                }


                // add the element force vector to the global force vector
                for (int j = 0; j < el.Nodes.Count; j++)
                {

                    int globalInd = el.Connectivity[j]; // index of the global node
                    double deltaR_z = elForceVec[j * 3 + 2]; // since the gravity load is zero in x and y direction these are neglected. 

                    F_body[globalInd * 3 + 2] += deltaR_z; // add the local contribution;
                }
                         

            }

            // I need the shape functions for the element. 8 noded element mean eight nodes. Should ensure that this is applicable for other types of elements as well. 
            // control the sum of the element load. Should be equal to the total weight
            double numericalWeight = F_body.Sum();

            return F_body;
        }

        /// <summary>
        /// Returns the natural coordinate of Gauss points in a matrix
        /// </summary>
        /// <param name="numPoints">e.g. 2 for 2x2x2x gauss points</param>
        /// <param name="elType"> Which element to work on use "Hex8" for now. </param>
        /// <returns></returns>
        public static LA.Matrix<double> GetGaussPointMatrix(int numPoints = 2, string elType = "Hex8")
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
                    //var naturalCoordinatesGauss = new CSD.DenseMatrix(8, 3, gaussArray);
                    var naturalCoordinatesGauss = LA.Matrix<double>.Build.DenseOfRowMajor(8, 3, gaussArray);
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

                double[] shapeArray = new double[] { N1, N2, N3, N4, N5, N6, N7, N8 };
                var shapeFunctions = CSD.DenseMatrix.OfColumnMajor(1, shapeArray.Length, shapeArray);
                return shapeFunctions;
            }
            else
            {
                throw new NotImplementedException("The selected element is not yet implemented.");
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
        public static LA.Matrix<double> PartialDerivateShapeFunctions(double r, double s, double t, string elType = "Hex8")
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
                //var derivativeMatrix = CSD.DenseMatrix.OfArray(derivateArray);
                var derivativeMatrix = LA.Matrix<double>.Build.DenseOfArray(derivateArray);
                //derivativeMatrix = derivativeMat
                return derivativeMatrix; 

            }
            else
            {
                throw new NotImplementedException("The selected element type is not implemented.");
            }
        }

        public static LA.Matrix<double> DisplacementInterpolationMatrix(CSparse.Storage.DenseColumnMajorStorage<double> shapeFunctions, int dofs)
        {
            var interpolationMatrix = LA.Matrix<double>.Build.Dense(dofs, shapeFunctions.ColumnCount * 3);

            for (int i = 0; i < shapeFunctions.ColumnCount; i++)
            {
                double[] shapeFunc = new double[dofs];
                for (int j = 0; j < dofs; j++){ shapeFunc[j] = shapeFunctions[0, i];}

                //var subMat = CSD.DenseMatrix.OfDiagonalArray(shapeFunc);
                var subMat = LA.Matrix<double>.Build.DenseOfDiagonalArray(shapeFunc);
                interpolationMatrix.SetSubMatrix(0, i * dofs, subMat);
            }
            return interpolationMatrix;

        }

        public static List<Point3d> SortedVerticesGraham(Mesh mesh)
        {
            // -- find the lowest and highest face of the mesh --
            int lowface = 0;
            int highFace = 0;


            int count = 0;

            // get the first centroid: 
            var p0 = mesh.TopologyVertices[mesh.Faces[0].A];
            var p1 = mesh.TopologyVertices[mesh.Faces[0].B];
            var p2 = mesh.TopologyVertices[mesh.Faces[0].C];
            var p3 = mesh.TopologyVertices[mesh.Faces[0].D];

            List<Point3d> pts = new List<Point3d>() { p0, p1, p2, p3 };
            double x = 0; double y = 0; double z = 0;
            foreach (Point3d point in pts)
            {
                x += point.X;
                y += point.Y;
                z += point.Z;
            }
            Point3d currentFaceCenter = new Point3d(x / pts.Count, y / pts.Count, z / pts.Count);

            for (int i = 1; i < mesh.Faces.Count; i++)
            {
                // get the first centroid: 
                p0 = mesh.TopologyVertices[mesh.Faces[i].A];
                p1 = mesh.TopologyVertices[mesh.Faces[i].B];
                p2 = mesh.TopologyVertices[mesh.Faces[i].C];
                p3 = mesh.TopologyVertices[mesh.Faces[i].D];

                pts = new List<Point3d>() { p0, p1, p2, p3 };
                x = 0; y = 0; z = 0;
                foreach (Point3d point in pts)
                {
                    x += point.X;
                    y += point.Y;
                    z += point.Z;
                }
                Point3d faceCenter = new Point3d(x / pts.Count, y / pts.Count, z / pts.Count);

                if (faceCenter.Z > currentFaceCenter.Z)
                {
                    highFace = i;
                }
                else if (faceCenter.Z < currentFaceCenter.Z)
                {
                    lowface = i;
                }
            }

            // -- get the vertices from the lowest and highest face --
            List<Point3f> lowPts = new List<Point3f>() { mesh.TopologyVertices[mesh.Faces[lowface].A], mesh.TopologyVertices[mesh.Faces[lowface].B], mesh.TopologyVertices[mesh.Faces[lowface].C], mesh.TopologyVertices[mesh.Faces[lowface].D] };
            List<Point3f> highPts = new List<Point3f>() { mesh.TopologyVertices[mesh.Faces[highFace].A], mesh.TopologyVertices[mesh.Faces[highFace].B], mesh.TopologyVertices[mesh.Faces[highFace].C], mesh.TopologyVertices[mesh.Faces[highFace].D] };



            
            List<Point3d> sortedVertices = new List<Point3d>();

            Stack<Point3d> stack = new Stack<Point3d>();


            return sortedVertices;
        }
        public static List<Point3d> SortVerticesByGrahamScan(List<Point3d> vertices)
        {
            
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
            //topNodes = topNodes.OrderBy(pt => Math.Atan2(pt.Y - topNodes[0].Y, pt.X - topNodes[0].X)).ToList();
            List<Point3d> sortedTop = new List<Point3d>();
            while (topNodes.Count > 0)
            {
                GrahamScan(ref topNodes, ref sortedTop);
            }
            bottomNodes = bottomNodes.OrderBy(pt => pt.Y).ThenBy(pt => pt.X).ToList();
            //bottomNodes = bottomNodes.OrderBy(pt => Math.Atan2(pt.Y - bottomNodes[0].Y, pt.X - bottomNodes[0].X)).ToList();
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
            #endregion
            return sortedVertices;
        }

        private static void GrahamScan(ref List<Point3d> pts, ref List<Point3d> selPts)
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
