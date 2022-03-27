using System;
using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry;
using LA = MathNet.Numerics.LinearAlgebra;
using System.Drawing;
using CSparse;
using CSD = CSparse.Double;
using CSparse.Double.Factorization;
using CSparse.Storage;
using Grasshopper;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry.Collections;
using Matrix = Accord.Math.Matrix;
using Point = Rhino.Geometry.Point;


namespace SolidFEM.Classes
{
    static class FEM_Utility
    {
        public static Mesh AddMidEdgeNodes(Mesh mesh)
        {
            /*for (int i = 0; i < mesh.TopologyEdges.Count; i++)
            {
                Line meshEdge = mesh.TopologyEdges.EdgeLine(i);
                double spX = meshEdge.FromX;
                double spY = meshEdge.FromY;
                double spZ = meshEdge.FromZ;
                double epX = meshEdge.ToX;
                double epY = meshEdge.ToY;
                double epZ = meshEdge.ToZ;

                Point3d midPoint = new Point3d((spX + epX)/2, (spY + epY)/2, (spZ + epZ)/2);
                mesh.Vertices.Add(midPoint);
            }*/
            var meshPts = mesh.TopologyVertices;
            int[,] midNodeIndices = new int[,]  //Array for getting the correct ordering of midside nodes
            {
                {0,1},  // node 5 is between node 1 and 2
                {0,2},  // node 6 is between node 1 and 3
                {0,3},  // node 7 and so on ...
                {1,2},  // node 8...
                {2, 3}, // node 9...
                {1,3}   // node 10...
            };

            for (int i = 0; i < Matrix.Rows(midNodeIndices); i++)
            {
                Point3d pt1 = mesh.TopologyVertices[midNodeIndices[i, 0]];
                Point3d pt2 = mesh.TopologyVertices[midNodeIndices[i, 1]];
                double x1 = pt1.X;
                double x2 = pt2.X;
                double y1 = pt1.Y;
                double y2 = pt2.Y;
                double z1 = pt1.Z;
                double z2 = pt2.Z;

                Point3d midPoint = new Point3d((x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2);
                mesh.Vertices.Add(midPoint);
            }
            return mesh;
        }





        public static void ElementsFromMeshList(List<Mesh> mList, List<Point3d> globalNodePts,out List<Element> femElements)
        {
            List<Element> elements = new List<Element>(); // all the elements of the mesh
            List<Node> globalNodes = new List<Node>(); // all the nodes of the mesh // not in use
            
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
                el.ElementMesh = mList[IDCounter];
                if(elNodes.Count == 8)
                {
                    el.Type = "Hex8";
                }
                else if(elNodes.Count == 4)
                {
                    el.Type = "Tet4";
                }
                else if (elNodes.Count == 10)
                {
                    el.Type = "Tet10";
                }

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

        public static LA.Vector<double> GetBodyForceVector(Material material, List<Element> elements, int numGlobalNodes, FEMLogger logger)
        {
            // Initiate the empty body force vector: 
            var F_body1 = new CSD.DenseMatrix(numGlobalNodes * 3, 1);
            //CSD.DenseMatrix F_body = new CSD.DenseMatrix(numGlobalNodes * 3, 1);
            LA.Vector<double> F_body = LA.Vector<double>.Build.Dense(numGlobalNodes * 3); // create the empty load vector
            LA.Vector<double> bodyLoadVector = LA.Vector<double>.Build.DenseOfArray(new double[] {0, 0, -  (material.Weight * 9.81 * Math.Pow(10, -9) ) });
            double elementJacobianTest = 0;
            for (int i = 0; i < elements.Count; i++)
            {
                // get the current element
                Element el = elements[i];

                // first, get the global coordinates
                LA.Matrix<double> globalElementCoordinates = LA.Matrix<double>.Build.Dense(el.Nodes.Count, 3);
                //LA.Matrix<double> globalElementCoordinate = LA.Matrix<double>.Build.Dense(el.Nodes.Count, 3);
                //CSD.DenseMatrix globalElementCoordinate = new CSD.DenseMatrix(el.Nodes.Count, 3); // one column for x,y, and z coordinate
                for (int j = 0; j < el.Nodes.Count; j++)
                {
                    Node n = el.Nodes[j];
                    globalElementCoordinates[j, 0] = n.Coordinate.X;
                    globalElementCoordinates[j, 1] = n.Coordinate.Y;
                    globalElementCoordinates[j, 2] = n.Coordinate.Z;
                }

                // get the matrix of natural coordinates in gauss points
                int order = 2;  //Order for gauss integration

                // Create from, CSparse
                var gaussCoordinates = GetGaussPointMatrix(order, el.Type); // by defaul we have a 2x2x2 integration of Hex8 element

                // element force vector. Create from CSparse
                LA.Vector<double> elForceVec = LA.Vector<double>.Build.Dense(el.Nodes.Count * 3);

                // iterate through each gauss node
                for (int j = 0; j < gaussCoordinates.RowCount; j++)
                {
                    double r = gaussCoordinates[j, 0];
                    double s = gaussCoordinates[j, 1];
                    double t = gaussCoordinates[j, 2];

                    // get the partial derivatives evaluated at a gauss point
                    var partialDerivatives = PartialDerivateShapeFunctions(r, s, t, el.Type);

                    // get the jacobian
                    var jacobianOperator = partialDerivatives.Multiply(globalElementCoordinates);
                    double jacobianDeterminant = jacobianOperator.Determinant();
                    //double jacobianDet = FEM_Matrices.GetDeterminantJacobi(jacobianOperator, logger);
                    
                    elementJacobianTest += jacobianDeterminant;
                    //double jacobianDet = jacobianOperator. // how to calculate the  determinant????
                    // get the H matrix from the shape functions
                    var shapeFunctions = GetShapeFunctions(r, s, t, el.Type);

                    // the H-matrix is the displacement interpolation matrix. 
                    var interpolationMatrix = DisplacementInterpolationMatrix(shapeFunctions, 3);

                    //Get weights for gauss integration // Can later change GetGaussPointMatrix() to also give the wheights as output

                    double alpha_ijk = 1.0;
                    if (el.Type == "Hex8")
                    {
                        alpha_ijk = 1.0;    // with a 2x2x2 integration scheme the integration constant is 1.0 for all gauss points. 
                    }
                    else if (el.Type == "Tet4" || el.Type == "Tet10")
                    {
                        alpha_ijk = 0.25;   // with a 4 point integration scheme the integration constant is 0.25 for all gauss points
                    }

                    double t_ijk = 1.0; // not sure what this is yet
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
            //double numericalWeight = F_body.Sum();

            return F_body;
        }

        /// <summary>
        /// Returns the natural coordinate of Gauss points in a matrix
        /// </summary>
        /// <param name="order">e.g. 2 for quadratic gauss integration</param>
        /// <param name="elType"> Which element to work on use "Hex8" for now. </param>
        /// <returns></returns>
        public static LA.Matrix<double> GetGaussPointMatrix(int order, string elType)//= "Hex8")
        {
            // for now we only have Hex8 element implementation
            if (elType == "Hex8")
            {
                if(order == 2)
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
            else if (elType == "Tet4" || elType == "Tet10")
            {   if(order == 2)
                {
                    double alpha = 0.58541;
                    double beta = 0.13820;
                    double[] gaussArray = new double[]
                    {
                        alpha, beta, beta ,
                        beta, alpha, beta,
                        beta, beta, alpha,
                        beta, beta, beta
                    };
                    var naturalCoordinatesGauss = LA.Matrix<double>.Build.DenseOfRowMajor(4, 3, gaussArray);
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
        }

        /// <summary>
        /// Return a vector of a shape function evaluated in the r, s, and t natural coordinates.
        /// </summary>
        /// <param name="r"></param>
        /// <param name="s"></param>
        /// <param name="t"></param>
        /// <param name="elType"></param>
        /// <returns></returns>
        public static CSparse.Storage.DenseColumnMajorStorage<double> GetShapeFunctions(double r, double s, double t ,string elType)
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
            else if(elType == "Tet4")
            {
                double N1 = r;
                double N2 = s;
                double N3 = t;
                double N4 = 1 - r - s - t;

                double[] shapeArray = new double[] { N1, N2, N3, N4};
                var shapeFunctions = CSD.DenseMatrix.OfColumnMajor(1, shapeArray.Length, shapeArray);
                return shapeFunctions;
            }
            else if (elType == "Tet10")
            {
                double c = 1 - r - s - t;
                double N1 = r*(2*r - 1);
                double N2 = s*(2*s - 1);
                double N3 = t * (2 * t - 1);
                double N4 = c * (2 * c - 1);
                double N5 = 4 * r * s;
                double N6 = 4 * r * t;
                double N7 = 4 * r * c;
                double N8 = 4 * s * t;
                double N9 = 4 * t * c;
                double N10 = 4 * s * c;

                double[] shapeArray = new double[] { N1, N2, N3, N4, N5, N6, N7, N8, N9, N10 };
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
        public static LA.Matrix<double> PartialDerivateShapeFunctions(double r, double s, double t, string elType)
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

                var derivativeMatrix = LA.Matrix<double>.Build.DenseOfArray(derivateArray);
                //var derivativeMatrix = CSD.DenseMatrix.OfArray(derivateArray);
                //derivativeMatrix = derivativeMat
                return derivativeMatrix; 

            }
            else if (elType == "Tet4")
            {
                double[,] derivateArray = new double[,]
                {
                    {1, 0, 0, -1 },
                    {0, 1, 0, -1 },
                    {0, 0, 1, -1 }
                };
                var derivativeMatrix = LA.Matrix<double>.Build.DenseOfArray(derivateArray);
                //var derivativeMatrix = CSD.DenseMatrix.OfArray(derivateArray);
                return derivativeMatrix;
            }
            else if (elType == "Tet10")
            {
                double[,] derivateArray = new double[,]
                {
                    {4*r - 1, 0, 0, 4*r + 4*s + 4*t - 3, 4*s, 4*t, 4 - 8*r - 4*s - 4*t, 0, -4*t, -4*s},
                    {0, 4*s - 1, 0, 4*r + 4*s + 4*t - 3, 4*r, 0, -4*r, 4*t, -4*t, 4 - 4*r - 8*s - 4*t},
                    {0, 0, 4*t - 1, 4*r + 4*s + 4*t - 3, 0, 4*r, -4*r, 4*s, 4 - 4*r - 4*s -8*t, -4*s}
                };

                var derivativeMatrix = LA.Matrix<double>.Build.DenseOfArray(derivateArray);
                //var derivativeMatrix = CSD.DenseMatrix.OfArray(derivateArray);
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


        /// <summary>
        /// Include boundary conditions, reduce matrices and solve for displacement. 
        /// </summary>
        /// <returns> List of nodal displacement. </returns>
        public static CSD.DenseMatrix CalculateDisplacementCSparse(double[,] K_gl, CSD.DenseMatrix R_gl, List<List<int>> applyBCToDOF, ref FEMLogger logger)
        {
            var timer = new System.Diagnostics.Stopwatch();

            // Make list of boundary condistions
            logger.AddInfo("--- Displacement calculations ---");
            
            
            timer.Start();
            List<int> BCList = new List<int>();
            for (int i = 0; i < applyBCToDOF.Count; i++)
            {
                for (int j = 0; j < applyBCToDOF[0].Count; j++)
                {
                    BCList.Add(applyBCToDOF[i][j]); // list of 0 and 1 values for boundary condition for dof; true = 1, false = 0
                }
            }
            timer.Stop();
            logger.AddInfo("Time to add Boundary conditions: " + timer.ElapsedMilliseconds + " ms"); timer.Reset();


            // Apply boundary conditions to movement
            timer.Start();
            for (int i = 0; i < BCList.Count; i++)
            {
                for (int j = 0; j < BCList.Count; j++)
                {
                    if (BCList[i] == 1)
                    {
                        if (i != j)
                        {
                            //K_gl.Row(i).SetValue(0, j);
                            K_gl[i, j] = 0;
                            K_gl[j, i] = 0;
                        }
                        else
                        {
                            //K_gl.Row(i).SetValue(1, j);
                            K_gl[i, j] = 1;
                            R_gl[i, 0] = 0;
                        }
                    }
                }
            }
            timer.Stop();
            logger.AddInfo("Time to restrain boundary conditions global stiffness and load matrix: " + timer.ElapsedMilliseconds + " ms"); timer.Reset();


            // to do: Time this function. Could be super slow. Better ways to get the array?
            /*
            timer.Start();
            var stiffnessArray = K_gl.AsColumnMajorArray();
            timer.Stop();
            logger.AddInfo("Transforming MathNet Matrix to Column Major Array: " + timer.ElapsedMilliseconds + " ms"); timer.Reset();

            timer.Start();
            //CompressedColumnStorage<double> CCS = CSD.SparseMatrix.OfColumnMajor(K_gl.RowCount, K_gl.ColumnCount, stiffnessArray);
            timer.Stop();
            logger.AddInfo("Create compressed column storage from array: " + timer.ElapsedMilliseconds + " ms"); timer.Reset();
            */
            timer.Start();
            var CCS =  CSD.SparseMatrix.OfArray(K_gl); // Try this instead. Need to convert the K_gl from MathNet to CSparse. 
            timer.Stop(); logger.AddInfo("Create compressed Column Storage of K: " + timer.ElapsedMilliseconds + " ms"); timer.Reset();
            #region Testing problems
            
            /*
            timer.Start();
            var R_LA = LA.Double.DenseVector.Build.DenseOfArray(R_gl.Values);
            var u_LA = K_gl.Solve(R_LA);
            timer.Stop();
            logger.AddInfo("Time spent on testing a MathNet solver: " + timer.ElapsedMilliseconds + " ms"); timer.Reset();
            */
            #endregion


            timer.Start();
            SparseLU CS_K_global = SparseLU.Create(CCS, ColumnOrdering.MinimumDegreeAtPlusA, 0.0);
            //SparseLDL CS_K_global = SparseLDL.Create(CCS, ColumnOrdering.MinimumDegreeAtPlusA); // Try an LDL system instead to test the speed
            timer.Stop();
            logger.AddInfo("Create LU decomposition of stiffness matrix: " + timer.ElapsedMilliseconds + " ms"); timer.Reset();

            //double[] CS_u = CSD.Vector.Create(K_global_red.RowCount * 1, 0.0);

            timer.Start();
            double[] CS_u = CSD.Vector.Create(K_gl.GetLength(0), 0.0);
            //double[] CS_R = R_red.Column(0).ToArray();
            double[] CS_R = R_gl.Column(0).ToArray();
            timer.Stop();
            logger.AddInfo("Prepare u and R arrays before solving for u: " + timer.ElapsedMilliseconds + " ms"); timer.Reset();

            /*
            for (int i = 0; i < CS_R.Length; i++)
            {
                CS_R[i] = Math.Round(CS_R[i], 6);
            }
            */

            timer.Start();
            CS_K_global.Solve(CS_R, CS_u);
            timer.Stop();
            logger.AddInfo("Time spent solving the equations using LU: " + timer.ElapsedMilliseconds + " ms"); timer.Reset();
            //logger.AddInfo("Time spent solving the equations using LDL: " + timer.ElapsedMilliseconds + " ms"); timer.Reset();

            timer.Start();
            var u = new CSD.DenseMatrix(CS_u.Length, 1, CS_u);
            //var u = new CSD.DenseMatrix(CS_u.Length, 1, u_LA.ToArray());
            timer.Stop();
            logger.AddInfo("Making a CSparse matrix from u-array: " + timer.ElapsedMilliseconds + "ms"); timer.Reset();
            //LA.Matrix<double> u = LA.Double.DenseMatrix.OfColumnArrays(CS_u);
            logger.AddInfo("----- Displacement calculation done -----");
            return u;


        }

        public static LA.Matrix<double> CalculateDisplacement(LA.Matrix<double> K_gl, LA.Matrix<double> R_gl, List<List<int>> applyBCToDOF, ref FEMLogger logger)
        {
            // summary: include boundary condistions and calculate global displacement
            var timer = new System.Diagnostics.Stopwatch();

            // Make list of boundary condistions
            logger.AddInfo("--- Displacement calculations ---");
            timer.Start();
            List<int> BCList = new List<int>();
            for (int i = 0; i < applyBCToDOF.Count; i++)
            {
                for (int j = 0; j < applyBCToDOF[0].Count; j++)
                {
                    BCList.Add(applyBCToDOF[i][j]); // list of 0 and 1 values for boundary condition for dof; true = 1, false = 0
                }
            }

            for (int i = 0; i < BCList.Count; i++)
            {
                for (int j = 0; j < BCList.Count; j++)
                {

                    if (BCList[i] == 1)
                    {
                        if (i != j)
                        {
                            K_gl[i, j] = 0;
                        }
                        else
                        {
                            K_gl[i, j] = 1;
                            R_gl[i, 0] = 0;
                        }
                    }
                }
            }

            timer.Stop();
            logger.AddInfo($"Time elapsed during boundary conditions: {timer.ElapsedMilliseconds} ms"); timer.Reset();

            // Reduce K_global and R
            timer.Start();

            // -- EDIT SVERRE ---
            //int numRows = K_gl.RowCount;
            //ReduceMatrices(BCList, numRows, ref K_gl,ref R_gl);

            //var reducedData = ReduceKandR(K_gl, R_gl, BCList);

            //LA.Matrix<double> K_global_red = reducedData.Item1;
            //LA.Matrix<double> R_red = reducedData.Item2;
            timer.Stop();
            logger.AddInfo($"Time elapse for reducing K and R: {timer.ElapsedMilliseconds} ms"); timer.Reset();

            // Time recorder
            var sw0 = new System.Diagnostics.Stopwatch();
            timer.Start();
            // Mathnet.Numerics to CSparse
            //var CMA = K_global_red.Storage.ToColumnMajorArray();
            var CMA = K_gl.Storage.ToColumnMajorArray();
            //CompressedColumnStorage<double> CCS = CSD.SparseMatrix.OfColumnMajor(K_global_red.RowCount, K_global_red.ColumnCount, CMA);
            CompressedColumnStorage<double> CCS = CSD.SparseMatrix.OfColumnMajor(K_gl.RowCount, K_gl.ColumnCount, CMA);
            timer.Stop();
            logger.AddInfo($"Convert MathNet to CSparse: {timer.ElapsedMilliseconds} ms"); timer.Reset();


            SparseLU CS_K_global = SparseLU.Create(CCS, ColumnOrdering.MinimumDegreeAtPlusA, 0.0);
            //double[] CS_u = CSD.Vector.Create(K_global_red.RowCount * 1, 0.0);
            double[] CS_u = CSD.Vector.Create(K_gl.RowCount * 1, 0.0);
            //double[] CS_R = R_red.Column(0).ToArray();
            double[] CS_R = R_gl.Column(0).ToArray();

            timer.Start();
            sw0.Start();
            CS_K_global.Solve(CS_R, CS_u);
            sw0.Stop();
            timer.Stop();
            logger.AddInfo($"Solve the system for displacements: {timer.ElapsedMilliseconds} ms"); timer.Reset();
            //Rhino.RhinoApp.WriteLine($"### {K_global_red.RowCount} x {K_global_red.ColumnCount} Matrix. CSparse Elapsed [msec] = {sw0.Elapsed.TotalMilliseconds}");
            //Rhino.RhinoApp.WriteLine($"### {K_gl.RowCount} x {K_gl.ColumnCount} Matrix. CSparse Elapsed [msec] = {sw0.ElapsedMilliseconds} "); sw0.Reset();

            // CSparse to Mathnet.Numerics
            timer.Start();
            LA.Matrix<double> u = LA.Double.DenseMatrix.OfColumnArrays(CS_u);
            timer.Stop();
            logger.AddInfo($"Convert back to MathNet: {timer.ElapsedMilliseconds} ms"); timer.Reset();
            logger.AddInfo("--- End displacement calculations ---");


            // Get total displacement
            // comment this out when skipping the row and column reduction of stiffness matrix
            /*
            LA.Vector<double> insertVec = DenseVector.Build.Dense(1);

            for (int i = 0; i < BCList.Count; i++)
            {
                if (BCList[i] == 1)
                {
                    u = u.InsertRow(i, insertVec);
                }
            }
            */
            return u;
        }

        /// <summary>
        /// Calculate a list of strain and stress vectors for each node in a element.
        /// </summary>
        /// <returns> Strain and stress vectors for each node in a element. </returns>
        public static Tuple<LA.Matrix<double>, LA.Matrix<double>> CalculateElementStrainStress(Element element, LA.Matrix<double> u, Material material, ref FEMLogger logger)
        {
            LA.Matrix<double> C = material.GetMaterialConstant();


            //List<LA.Matrix<double>> B_local = FEM_Matrices.CalculateElementMatrices(element, material, ref logger).Item2; // this can be changed to save time.. No need to establish the stiffness matrix of an element for this
            List<LA.Matrix<double>> B_local = element.LocalB;
            LA.Matrix<double> elementGaussStrain = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            LA.Matrix<double> elementGaussStress = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            LA.Matrix<double> elementStrain = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            LA.Matrix<double> elementStress = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            LA.Matrix<double> localDeformation = LA.Double.DenseMatrix.Build.Dense(3 * element.Nodes.Count, 1); // Could this be attached to the FE_Element class as well?

            // get deformation of nodes connected to element
            for (int i = 0; i < element.Connectivity.Count; i++)
            {
                localDeformation[3 * i, 0] = u[3 * element.Connectivity[i], 0];
                localDeformation[3 * i + 1, 0] = u[3 * element.Connectivity[i] + 1, 0];
                localDeformation[3 * i + 2, 0] = u[3 * element.Connectivity[i] + 2, 0];
            }
            if (element.Type == "Hex8")
            {
                // get gauss strain and stress
                for (int n = 0; n < B_local.Count; n++)
                {
                    // B-matrix is calculated from gauss points
                    LA.Matrix<double> gaussStrain = B_local[n].Multiply(localDeformation);
                    LA.Matrix<double> gaussStress = C.Multiply(gaussStrain);
                    //LA.Matrix<double> gaussStress = C.Multiply(B_local[n]).Multiply(localDeformation);


                    // Use add column here instead of a loop?
                    for (int i = 0; i < B_local[0].RowCount; i++)
                    {
                        elementGaussStrain[i, n] = gaussStrain[i, 0]; // Should there not be += here? Right now they overwrite the Gauss Strainf for each nodal evaluation
                        elementGaussStress[i, n] = gaussStress[i, 0];
                    }
                }

                // get node strain and stress by extrapolation
                LA.Matrix<double> extrapolationNodes = FEM.GetNaturalCoordinate(Math.Sqrt(3), 3);

                for (int n = 0; n < B_local.Count; n++)
                {
                    // get stress and strain in nodes
                    var r = extrapolationNodes.Row(n)[0];
                    var s = extrapolationNodes.Row(n)[1];
                    double t = extrapolationNodes.Row(n)[2];

                    LA.Vector<double> shapefunctionValuesInNode = FEM.GetShapeFunctions(r, s, t, 3);
                    LA.Vector<double> nodeStrain = elementGaussStrain.Multiply(shapefunctionValuesInNode);
                    LA.Vector<double> nodeStress = elementGaussStress.Multiply(shapefunctionValuesInNode);
                    for (int i = 0; i < B_local[0].RowCount; i++)
                    {
                        elementStrain[i, n] = nodeStrain[i];
                        elementStress[i, n] = nodeStress[i];
                    }
                }
            }
            else if (element.Type == "Tet4")    // no need for gauss strains because B-matrix is constant
            {

                for (int i = 0; i < element.Nodes.Count; i++)
                {

                    //Get deformation at node
                    LA.Matrix<double> nodalDeformation = LA.Double.DenseMatrix.Build.Dense(3, 1);
                    LA.Matrix<double> nodeStrain = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, 1);
                    LA.Matrix<double> nodeStress = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, 1);


                    nodeStrain = (B_local[0]).Multiply(localDeformation);
                    nodeStress.Add(C.Multiply(nodeStrain), nodeStress);

                    elementStrain.SetSubMatrix(0, i, nodeStrain);
                    elementStress.SetSubMatrix(0, i, nodeStress);

                }
            }

            return Tuple.Create(elementStrain, elementStress);
        }







        /// <summary>
        /// Assemble element stress and get global stress and mises stress,
        /// </summary>
        /// <returns> Nodal global stress, node mises stress and element mises stress. </returns>
        /// 
        public static Tuple<LA.Matrix<double>, LA.Vector<double>, LA.Vector<double>> CalculateGlobalStress(List<Element> elements, LA.Matrix<double> u, Material material, ref FEMLogger logger)
        {
            int numNodes = u.RowCount / 3;
            int stressRowDim = 6;
            LA.Matrix<double> globalStress = LA.Double.DenseMatrix.Build.Dense(stressRowDim, numNodes);
            LA.Matrix<double> counter = LA.Double.DenseMatrix.Build.Dense(stressRowDim, numNodes);
            List<LA.Matrix<double>> elementStressList = new List<LA.Matrix<double>>();
            foreach (Element element in elements)
            {
                LA.Matrix<double> elementStress = CalculateElementStrainStress(element, u, material, ref logger).Item2;

                List<int> connectivity = element.Connectivity;

                for (int i = 0; i < elementStress.RowCount; i++) // loop the stress
                {
                    for (int j = 0; j < elementStress.ColumnCount; j++) // loop the element nodes
                    {
                        globalStress[i, connectivity[j]] = globalStress[i, connectivity[j]] + elementStress[i, j];
                        counter[i, connectivity[j]]++;
                    }
                }
                elementStressList.Add(elementStress);
            }

            // get average
            for (int i = 0; i < globalStress.RowCount; i++) // loop the stress
            {
                for (int j = 0; j < globalStress.ColumnCount; j++) // loop the element nodes
                {
                    if (counter[i, j] > 1)
                    {
                        globalStress[i, j] = globalStress[i, j] / (double)counter[i, j];
                        counter[i, j] = 0;
                    }
                }
            }

            // Nodal Mises
            LA.Vector<double> mises = LA.Double.DenseVector.Build.Dense(numNodes);
            for (int i = 0; i < numNodes; i++)
            {
                LA.Vector<double> nodeStress = globalStress.Column(i);
                double Sxx = nodeStress[0];
                double Syy = nodeStress[1];
                double Szz = nodeStress[2];
                double Sxy = nodeStress[3];
                double Sxz = nodeStress[4];
                double Syz = nodeStress[5];
                mises[i] = Math.Sqrt(0.5 * (Math.Pow(Sxx - Syy, 2) + Math.Pow(Syy - Szz, 2) + Math.Pow(Szz - Sxx, 2)) + 3 * (Math.Pow(Sxy, 2) + Math.Pow(Sxz, 2) + Math.Pow(Syz, 2)));
            }

            // Element mises
            LA.Vector<double> elementMises = LA.Double.DenseVector.Build.Dense(elements.Count);
            for (int i = 0; i < elementStressList.Count; i++)
            {
                for (int j = 0; j < elements[i].Nodes.Count; j++)
                {
                    LA.Vector<double> nodeStress = elementStressList[i].Column(j);
                    double Sxx = nodeStress[0];
                    double Syy = nodeStress[1];
                    double Szz = nodeStress[2];
                    double Sxy = nodeStress[3];
                    double Sxz = nodeStress[4];
                    double Syz = nodeStress[5];
                    elementMises[i] += Math.Sqrt(0.5 * (Math.Pow(Sxx - Syy, 2) + Math.Pow(Syy - Szz, 2) + Math.Pow(Szz - Sxx, 2)) + 3 * (Math.Pow(Sxy, 2) + Math.Pow(Sxz, 2) + Math.Pow(Syz, 2))) / elements[i].Nodes.Count;
                }
                // elementMises[i] = elementMises[i] / elements[i].Nodes.Count; // get average of nodal mises
            }


            return Tuple.Create(globalStress, mises, elementMises);
        }

        /// <summary>
        /// Color mesh after nodal stress values.
        /// </summary>
        public static void ColorMeshAfterStress(SmartMesh mesh, LA.Vector<double> mises, Material material)
        {
            double maxValue = material.YieldingStress;
            double minValue = 0;
            Color color = Color.White;

            double range = (maxValue - minValue) / (double)13;
            for (int i = 0; i < mesh.Nodes.Count; i++)
            {
                // to do: Reference Synne, same color mapping
                if (mises[i] < minValue + range) color = Color.Blue;
                else if (mises[i] < minValue + 2 * range) color = Color.RoyalBlue;
                else if (mises[i] < minValue + 3 * range) color = Color.DeepSkyBlue;
                else if (mises[i] < minValue + 4 * range) color = Color.Cyan;
                else if (mises[i] < minValue + 5 * range) color = Color.PaleGreen;
                else if (mises[i] < minValue + 6 * range) color = Color.LimeGreen;
                else if (mises[i] < minValue + 7 * range) color = Color.Lime;
                else if (mises[i] < minValue + 8 * range) color = Color.Lime;
                else if (mises[i] < minValue + 9 * range) color = Color.GreenYellow;
                else if (mises[i] < minValue + 10 * range) color = Color.Yellow;
                else if (mises[i] < minValue + 11 * range) color = Color.Orange;
                else if (mises[i] < minValue + 12 * range) color = Color.OrangeRed;
                else color = Color.Red;

                //mesh.Mesh.VertexColors.Add(color);
            }
        }

        public static List<Point3d> SortedVerticesGraham(Mesh mesh)
        {
            // -- find the lowest and highest face of the mesh --
            int lowface = 0;
            int highFace = 0;
            

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
