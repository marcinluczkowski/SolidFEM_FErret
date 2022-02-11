using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Drawing;
using System.Linq;

// Csparse
using LA = MathNet.Numerics.LinearAlgebra;
using CSparse;
using CSD = CSparse.Double;
using CSparse.Double.Factorization;
using CSparse.Storage;


namespace SolidFEM.FiniteElementMethod
{
    public class NEW_FEMSolver : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMsolver class.
        /// </summary>
        public NEW_FEMSolver()
          : base("NEW FEM Solver", "Solver",
              "Solver for FEM problems with Rhino Mesh" +
                "Uses 3 translation DOFS pr node, linear shape functions, two Gauss Points and full integration.",
              "SmartMesh", "FEM-Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh", "m", "Normal mesh as a list of elements.", GH_ParamAccess.list); // 0
            pManager.AddGenericParameter("Loads", "load", "Load vector from FEM Load.", GH_ParamAccess.list); // 1
            pManager.AddGenericParameter("Boundary Conditions", "BC", "Boundary conditions from FEM Boundary Condtion.", GH_ParamAccess.list); // 2
            pManager.AddGenericParameter("Material", "material", "Material from FEM Material.", GH_ParamAccess.item); // 3
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("u1", "u1", "Displacement of nodes in x-direction.", GH_ParamAccess.list);
            pManager.AddGenericParameter("u2", "u2", "Displacement of nodes in y-direction.", GH_ParamAccess.list);
            pManager.AddGenericParameter("u3", "u3", "Displacement of nodes in z-direction.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Element Stress", "element mises", "List of Von Mises stress in element.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Node Stress", "node mises", "List of von Mises stress in node.", GH_ParamAccess.list);

            pManager.AddTextParameter("Diagonstics", "text", "List of information on the components performance", GH_ParamAccess.list);
            pManager.AddGenericParameter("FE_Mesh", "femesh", "The FE_Mesh containing results from the analysis.", GH_ParamAccess.item) ;
            pManager.AddGenericParameter("Global K", "", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            
            // stopwatch
            var watch = new System.Diagnostics.Stopwatch();
            var infoList = new List<string>(); // list to input information

            List<Mesh> meshList = new List<Mesh>();
            SmartMesh smartMesh = new SmartMesh();
            List<double> loads = new List<double>();
            //List<List<int>> boundaryConditions = new List<List<int>>();
            List<Support> supports = new List<Support>();
            Material material = new Material();

            DA.GetDataList(0, meshList);
            DA.GetDataList(1, loads);
            //DA.GetDataList(2, boundaryConditions);
            DA.GetDataList(2, supports);
            DA.GetData(3, ref material);


            // 0. Initial step

            //List<Node> nodes = smartMesh.Nodes;
            //List<Element> elements = smartMesh.Elements;

            // clean the mesh and sort nodes
            var newMeshList = new List<Mesh>();
            int c = 0; // delete after testing
            foreach (Mesh mesh in meshList)
            {
                if (mesh.Vertices.Count == 8)
                {
                    Mesh nM = GrahamScan.DoGrahamScan(mesh);

                    if (nM.IsValid)
                    {
                        newMeshList.Add(nM);
                    }
                    else newMeshList.Add(mesh);
                    c++;
                }
                else
                {
                    newMeshList.Add(mesh);
                }
            }


            
            List<Element> elements;
            List<Node> femNodes = new List<Node>();

            List<Point3d> nodePos = FEM_Utility.GetMeshNodes(newMeshList);

            for (int i = 0; i < nodePos.Count; i++)
            {
                Node node = new Node(i, nodePos[i]);
                femNodes.Add(node);
            }
            

            int numNodes = nodePos.Count;
            FEM_Utility.ElementsFromMeshList(newMeshList, nodePos ,out elements);

            // 1. Get global stiffness matrix

            watch.Start();
            //LA.Matrix<double> K_global = CalculateGlobalStiffnessMatrix(elements, numNodes, material);
            var K_globalC = GlobalStiffnessCSparse(elements, numNodes, material);
            var sumStiffness = K_globalC.AsColumnMajorArray();
            infoList.Add($"The sum of all elements in the stiffness matrix from mesh elements:  {sumStiffness.Sum()}");

            //var K_sparse = new CSD.SparseMatrix()
            watch.Stop();
            infoList.Add($"Time used calculating global stiffness matrix: {watch.ElapsedMilliseconds} ms"); watch.Reset();

            // 3. Get load vector
            watch.Start();

            
            // self weight
            var selfWeight = FEM_Utility.GetBodyForceVector(material, elements, numNodes);
            CSD.DenseMatrix R_self = new CSD.DenseMatrix(numNodes * 3, 1, selfWeight.ToArray());
            CSD.DenseMatrix R_external = new CSD.DenseMatrix(numNodes * 3, 1);
            
            //LA.Matrix<double> R = LA.Double.DenseMatrix.Build.Dense(numNodes * 3, 1);
            for (int i = 0; i < loads.Count; i++)
            {
                //R[i, 0] = loads[i];
                R_external[i,0] = loads[i];
            }

            CSD.DenseMatrix R = (CSD.DenseMatrix)R_self.Add(R_external);
            
            watch.Stop();
            infoList.Add($"Time used to establish global load vector: {watch.ElapsedMilliseconds} ms"); watch.Reset();

            // 5. Fix BoundaryConditions
            //boundaryConditions = FixBoundaryConditions(boundaryConditions, smartMesh.Nodes.Count);
            watch.Start();
            List<List<int>> boundaryConditions = FixBoundaryConditionsSverre(supports, nodePos);
            infoList.Add($"Time used on boundary conditions: {watch.ElapsedMilliseconds} ms"); watch.Reset();

            // 6. Calculate displacement 
            watch.Start();
            CSD.DenseMatrix u_CSparse = CalculateDisplacementCSparse(K_globalC, R, boundaryConditions, ref infoList);
            infoList.Add($"Time used on displacement calculations with CSparce: {watch.ElapsedMilliseconds} ms"); watch.Reset();

            

            watch.Start();
            //LA.Matrix<double> u = CalculateDisplacement(K_global, R, boundaryConditions, ref infoList);
            infoList.Add($"Time used on regular displacement calculations: {watch.ElapsedMilliseconds} ms"); watch.Reset();

            // 7. Calculate stress

            // convert CSparse to double and MathNet
            double[] u_val = u_CSparse.Values;
            var uCS2MN = new LA.Double.DenseMatrix(u_val.Length, 1, u_val);

            watch.Start();
            var stress = CalculateGlobalStress(elements, uCS2MN, material); // make this compatible with the CSparse matrix as well.
            infoList.Add($"Time used on stress calculations: {watch.ElapsedMilliseconds} ms"); watch.Reset();

            LA.Matrix<double> globalStress = stress.Item1;
            LA.Vector<double> mises = stress.Item2;
            LA.Vector<double> misesElement = stress.Item3;
            watch.Start();
            ColorMeshAfterStress(smartMesh, mises, material);
            infoList.Add($"Time used on mesh colouring: {watch.ElapsedMilliseconds} ms"); watch.Reset();

            // 8. Prepare output
            watch.Start();
            List<double> u1 = new List<double>();
            List<double> u2 = new List<double>();
            List<double> u3 = new List<double>();
            List<double> nodalMises = new List<double>();
            List<double> elementMises = new List<double>();

            for (int i = 0; i < elements.Count; i++)
            {
                elementMises.Add(misesElement[i]);
            }

            for (int i = 0; i < numNodes; i++)
            {
                u1.Add(uCS2MN[i * 3, 0]);
                u2.Add(uCS2MN[i * 3 + 1, 0]);
                u3.Add(uCS2MN[i * 3 + 2, 0]);

                nodalMises.Add(mises[i]);
            }

            // test - delete after
            double maxDisp = u3.Max(x => Math.Abs(x));

            List<double[]> nodalStress = new List<double[]>();
            for (int i = 0; i < globalStress.ColumnCount; i++)
            {
                nodalStress.Add(globalStress.Column(i).ToArray());
            }
            infoList.Add($"Time used on output preparation: {watch.ElapsedMilliseconds} ms"); watch.Reset();

            // TEMPORARY WHILE WAITING FOR MARCIN's MESH CLASS
            TempFE_Mesh outMesh = new TempFE_Mesh(meshList, femNodes, elements, elementMises, nodalMises ,u1, u2, u3, material,
                globalStress.Row(0).ToList(), globalStress.Row(1).ToList(), globalStress.Row(2).ToList(),
                globalStress.Row(3).ToList(), globalStress.Row(4).ToList(), globalStress.Row(5).ToList());
            List<double> test = globalStress.Row(5).ToList();
            // Output
            DA.SetDataList(0, u1);
            DA.SetDataList(1, u2);
            DA.SetDataList(2, u3);
            DA.SetDataList(3, elementMises);
            DA.SetDataList(4, nodalMises);
            // temporary information
            DA.SetDataList(5, infoList);
            DA.SetData(6, outMesh);
            DA.SetDataList(7, K_globalC.AsColumnMajorArray());
        }

        #region Methods
        /// <summary>
        /// Sort multiple boundary conditions to one total boundary condition list. 
        /// </summary>
        /// <returns> A list of boolean values for each dof of each node. </returns>
        private List<List<int>> FixBoundaryConditions(List<List<int>> boundaryConditions, int numNodes)
        {
            List<List<int>> totalBC = new List<List<int>>();
            for (int i = 0; i < numNodes; i++) // loop number nodes
            {
                List<int> dofList = new List<int>(boundaryConditions[i]);  // get dofList of first input list of BC
                for (int j = 0; j < dofList.Count; j++) // loop dofs
                {
                    for (int k = 1; k < boundaryConditions.Count / numNodes; k++) // loop the remaining inout list of BC 
                    {
                        dofList[j] = dofList[j] + boundaryConditions[i + k * numNodes][j];
                    }
                }
                totalBC.Add(dofList);
            }
            return totalBC;
        }
        private List<List<int>> FixBoundaryConditionsSverre(List<Support> sups, List<Point3d> nodesPos)
        {
            List<List<int>> totalBC = new List<List<int>>(); // just using the names from the above function

            // iterate through each support in the mesh
            foreach (Point3d nodePos in nodesPos)
            {
                List<int> bc = new List<int>() { 0, 0, 0 };
                var sup = sups.FirstOrDefault(n => n.Position.DistanceToSquared(nodePos) < 0.001);
                if (!(sup is null))
                {
                    // if there is a support on this node
                    if (sup.Tx == true) { bc[0] = 1; }
                    if (sup.Ty == true){bc[1] = 1;}
                    if (sup.Tz == true) { bc[2] = 1; }

                }
                totalBC.Add(bc);
            }
            /*
            foreach (Support support in sups)
            {
                List<Node> meshNodes = sMesh.Nodes;
                var supNode = meshNodes.FirstOrDefault(n => n.Coordinate.DistanceToSquared(support.Position) < 0.001); // find a node (if present) where
                if(!(supNode is null))
                {
                    // if there is a node close enough: 

                }

            }*/


            return totalBC;
        }

        /// <summary>
        /// Calculate element stifness matrix and element strain matrix.
        /// </summary>
        /// <returns> Element stiffness and strain matrix.</returns>
        private Tuple<LA.Matrix<double>, List<LA.Matrix<double>>> CalculateElementMatrices(Element element, Material material)
        {
            // summary: calculate local K and B matrix
            int roundrecisionBMatrix = 6;
            int rpb = roundrecisionBMatrix;
            // material
            LA.Matrix<double> C = material.GetMaterialConstant();

            // shapefunction
            // create local stiffness matrix
            int numElementNodes = element.Nodes.Count;
            LA.Matrix<double> K_local = LA.Matrix<double>.Build.Dense(3 * numElementNodes, 3 * numElementNodes);

            // create local deformation matrix
            List<LA.Matrix<double>> B_local = new List<LA.Matrix<double>>();

            // Global coordinates of the (corner) nodes of the actual element
            LA.Matrix<double> globalCoordinates = LA.Matrix<double>.Build.Dense(numElementNodes, 3);

            List<Point3d> localCoordinates = FEM_Utility.LocalCartesianCoordinates(element);
            /*
            // this does not seem to have any effect on the global displacements
            for (int i = 0; i < localCoordinates.Count; i++)
            {
                globalCoordinates[i, 0] = localCoordinates[i].X;
                globalCoordinates[i, 1] = localCoordinates[i].Y;
                globalCoordinates[i, 2] = localCoordinates[i].Z;
            }*/
            
            for (int i = 0; i < numElementNodes; i++)
            {
                globalCoordinates[i, 0] = Math.Round(element.Nodes[i].Coordinate.X, rpb); // column of x coordinates
                globalCoordinates[i, 1] = Math.Round(element.Nodes[i].Coordinate.Y, rpb) ; // column of y coordinates
                globalCoordinates[i, 2] = Math.Round(element.Nodes[i].Coordinate.Z, rpb); // colum of z coordinates
            }

            // Different methods for Hex8 and Tet4. Tet4 doesn't need gauss integration because B and Jacobian are constant!
            if (element.Type == "Hex8")
            {
                //Numerical integration
                //LA.Matrix<double> gaussNodes = FEM.GetNaturalCoordinate((double)Math.Sqrt((double)1 / (double)3), 3);
                var gaussCoordinates = FEM_Utility.GetGaussPointMatrix(2, element.Type); // by defaul we have a 2x2x2 integration of Hex8 element
                for (int n = 0; n < gaussCoordinates.RowCount; n++)  // loop gauss nodes
                {
                    // Substitute the natural coordinates into the symbolic expression
                    var r = gaussCoordinates.Row(n)[0];
                    var s = gaussCoordinates.Row(n)[1];
                    var t = gaussCoordinates.Row(n)[2];

                    // Partial derivatives of the shape functions
                    //LA.Matrix<double> shapeFunctionsDerivatedNatural = FEM.DerivateWithNatrualCoordinates(r, s, t, 3);
                    var partialDerivatives = FEM_Utility.PartialDerivateShapeFunctions(r, s, t, element.Type);

                    // Calculate Jacobian matrix
                    LA.Matrix<double> jacobianMatrix = partialDerivatives.Multiply(globalCoordinates);

                    // round the jacobian
                    var jColMajArr = jacobianMatrix.AsColumnMajorArray();
                    for (int k = 0; k < jColMajArr.Length; k++)
                    {
                        jColMajArr[k] = Math.Round(jColMajArr[k], rpb);
                    }
                    jacobianMatrix = LA.Matrix<double>.Build.DenseOfColumnMajor(3, 3, jColMajArr);
                    // Calculate B - LA.Matrix
                    LA.Matrix<double> shapeFuncDerivatedCartesian = (jacobianMatrix.Inverse()).Multiply(partialDerivatives);

                    double jacobianDeterminant = jacobianMatrix.Determinant();
                    if (jacobianDeterminant < 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Negativ jac det"); }
                    int dimRowB = 6;

                    // establish the B-matrix
                    LA.Matrix<double> B_i = LA.Double.DenseMatrix.Build.Dense(dimRowB, 3 * numElementNodes);

                    for (int i = 0; i < numElementNodes; i++)
                    {

                        // with the shape functions derivated with respect to the cartesian coordinates the rotated and unrotated element vectors are not the same... This is the correct one according to the formulas
                        var B_i_sub = DenseMatrix.Build.DenseOfRowMajor(6, 3, new double[] {
                            shapeFuncDerivatedCartesian.Row(0)[i], 0, 0,
                            0, shapeFuncDerivatedCartesian.Row(1)[i], 0,
                            0, 0, shapeFuncDerivatedCartesian.Row(2)[i],
                            shapeFuncDerivatedCartesian.Row(1)[i], shapeFuncDerivatedCartesian.Row(0)[i], 0,
                            shapeFuncDerivatedCartesian.Row(2)[i], 0, shapeFuncDerivatedCartesian.Row(0)[i],
                            0, shapeFuncDerivatedCartesian.Row(2)[i], shapeFuncDerivatedCartesian.Row(1)[i]
                            });
                        /*
                        // with the shape functions derivated with respect to the cartesian coordinates the rotated and unrotated element vectors are not the same... This is the correct one according to the formulas
                        var B_i_sub = DenseMatrix.Build.DenseOfRowMajor(6, 3, new double[] {
                            Math.Round(shapeFuncDerivatedCartesian.Row(0)[i], rpb), 0, 0,
                            0, Math.Round(shapeFuncDerivatedCartesian.Row(1)[i], rpb), 0,
                            0, 0, Math.Round(shapeFuncDerivatedCartesian.Row(2)[i], rpb),
                            Math.Round(shapeFuncDerivatedCartesian.Row(1)[i], rpb), Math.Round(shapeFuncDerivatedCartesian.Row(0)[i],rpb), 0,
                            Math.Round(shapeFuncDerivatedCartesian.Row(2)[i], rpb), 0, Math.Round(shapeFuncDerivatedCartesian.Row(0)[i],rpb),
                            0, Math.Round(shapeFuncDerivatedCartesian.Row(2)[i],rpb), Math.Round(shapeFuncDerivatedCartesian.Row(1)[i],rpb)
                            });
                        */

                        B_i.SetSubMatrix(0, i * 3, B_i_sub);

                    }

                    B_local.Add(B_i);
                    //K_local += ((B_i.Transpose()).Multiply(C).Multiply(B_i)).Multiply(jacobianDeterminant);
                    var k_i = (B_i.Transpose()).Multiply(C.Multiply(B_i)).Multiply(jacobianDeterminant);

                    K_local.Add(k_i, K_local);
                }
            }
            else if (element.Type == "Tet4")
            {
                var partialDerivatives = FEM_Utility.PartialDerivateShapeFunctions(1, 1, 1, "Tet4");    //Random coordinates (1,1,1) because the method requires coordinate inputs
                
                // Calculate Jacobian matrix
                LA.Matrix<double> jacobianMatrix = partialDerivatives.Multiply(globalCoordinates);

                // round the jacobian
                var jColMajArr = jacobianMatrix.AsColumnMajorArray();
                for (int k = 0; k < jColMajArr.Length; k++)
                {
                    jColMajArr[k] = Math.Round(jColMajArr[k], rpb);
                }
                jacobianMatrix = LA.Matrix<double>.Build.DenseOfColumnMajor(3, 3, jColMajArr);

                // Calculate B - LA.Matrix
                LA.Matrix<double> shapeFuncDerivatedCartesian = (jacobianMatrix.Inverse()).Multiply(partialDerivatives);
                
                int dimRowB = 6;
                LA.Matrix<double> B_i = LA.Double.DenseMatrix.Build.Dense(dimRowB, 3 * numElementNodes);

                for (int i = 0; i < numElementNodes; i++)
                {

                    // with the shape functions derivated with respect to the cartesian coordinates the rotated and unrotated element vectors are not the same... This is the correct one according to the formulas
                    var B_i_sub = DenseMatrix.Build.DenseOfRowMajor(6, 3, new double[] {
                            shapeFuncDerivatedCartesian.Row(0)[i], 0, 0,
                            0, shapeFuncDerivatedCartesian.Row(1)[i], 0,
                            0, 0, shapeFuncDerivatedCartesian.Row(2)[i],
                            shapeFuncDerivatedCartesian.Row(1)[i], shapeFuncDerivatedCartesian.Row(0)[i], 0,
                            shapeFuncDerivatedCartesian.Row(2)[i], 0, shapeFuncDerivatedCartesian.Row(0)[i],
                            0, shapeFuncDerivatedCartesian.Row(2)[i], shapeFuncDerivatedCartesian.Row(1)[i]
                            });
                   

                    B_i.SetSubMatrix(0, i * 3, B_i_sub);

                }
                B_local.Add(B_i);
                // Get volume of Tetrahedra
                Brep triangle1 = Brep.CreateFromCornerPoints(element.TopologyVertices[0], element.TopologyVertices[1], element.TopologyVertices[2], 0.0001);
                Brep triangle2 = Brep.CreateFromCornerPoints(element.TopologyVertices[0], element.TopologyVertices[1], element.TopologyVertices[3], 0.0001);
                Brep triangle3 = Brep.CreateFromCornerPoints(element.TopologyVertices[0], element.TopologyVertices[2], element.TopologyVertices[3], 0.0001);
                Brep triangle4 = Brep.CreateFromCornerPoints(element.TopologyVertices[1], element.TopologyVertices[2], element.TopologyVertices[3], 0.0001);

                List<Brep> triangles = new List<Brep> { triangle1, triangle2, triangle3, triangle4 };

                Brep[] tetra = Brep.CreateSolid(triangles, 0.0001);

                VolumeMassProperties vmp = VolumeMassProperties.Compute(tetra[0]);
                double V = vmp.Volume;

                var k_i = V*(B_i.Transpose()).Multiply(C.Multiply(B_i));

                K_local.Add(k_i, K_local);
            }
            return Tuple.Create(K_local, B_local);
        }

        /// <summary>
        /// Construct global stiffness matrix by assembling element stiffness matrices.
        /// </summary>
        /// <returns> Global stiffness matrix. </returns>
        
        
        private LA.Matrix<double> GlobalStiffnessCSparse(List<Element> elements, int numNode, Material material)
        {
            // Initiate empty matrix
            //CSD.DenseMatrix m = new CSD.DenseMatrix(numNode * 3, numNode * 3);
            LA.Matrix<double> m = LA.Matrix<double>.Build.Dense(numNode * 3, numNode * 3);
            List<double> elementSums = new List<double>();
            
            // test the functionality: 
            int[] testElInds = new int[] { 0, 21, 84, 105, 168, 189, 252, 273 };
            List<LA.Matrix<double>> testElsMat = new List<LA.Matrix<double>>();
            List<Element> testEls = new List<Element>();
            int count = 0; 
            // delete the above after finished test
            foreach (Element element in elements)
            {
                List<int> con = element.Connectivity; // get the connectivity of each element

                // iterate over the connectivity indices
                LA.Matrix<double> K_local = CalculateElementMatrices(element, material).Item1;
                elementSums.Add(K_local.AsColumnMajorArray().Sum(x => Math.Abs(x))) ; // the sum of each element

                if (testElInds.Contains(count))
                {
                    testElsMat.Add(K_local);
                    testEls.Add(element);
                }

                // loop nodes of elements
                for (int i = 0; i < con.Count; i++)
                {
                    for (int j = 0; j < con.Count; j++)
                    {
                        // loop relevant local stiffness contribution
                        for (int dofRow = 0; dofRow < 3; dofRow++)
                        {
                            for (int dofCol = 0; dofCol < 3; dofCol++)
                            {
                                m[3 * con[i] + dofRow, 3 * con[j] + dofCol] += K_local[3 * i + dofRow, 3 * j + dofCol];
                            }
                        }
                    }
                }
                count++; // delete after testing
            }

            // small code to round the values of the elements
            /*
            double[] colMayArray = m.AsColumnMajorArray();
            for (int i = 0; i < colMayArray.Length; i++)
            {
                colMayArray[i] = Math.Round(colMayArray[i], 14);
            }
            var globalK = LA.Matrix<double>.Build.DenseOfColumnMajor(numNode * 3, numNode * 3, colMayArray);
            */
            var sum_Element = m.AsColumnMajorArray().Sum();
            return m;
        }

        private LA.Matrix<double> CalculateGlobalStiffnessMatrix(List<Element> elements, int numNode, Material material)
        {            

            // create stiffness matrix
            LA.Matrix<double> K_global = LA.Matrix<double>.Build.Dense(numNode * 3, numNode * 3);
            foreach (Element element in elements)
            {
                
                List<int> con = element.Connectivity;

                
                LA.Matrix<double> K_local = CalculateElementMatrices(element, material).Item1;

                // loop nodes of elements
                for (int i = 0; i < con.Count; i++)
                {
                    for (int j = 0; j < con.Count; j++)
                    {
                        // loop relevant local stiffness contribution
                        for (int dofRow = 0; dofRow < 3; dofRow++)
                        {
                            for (int dofCol = 0; dofCol < 3; dofCol++)
                            {
                                K_global[3 * con[i] + dofRow, 3 * con[j] + dofCol] += + K_local[3 * i + dofRow, 3 * j + dofCol];
                            }
                        }
                    }
                }

            }
            return K_global;
        }

        /// <summary>
        /// Include boundary conditions, reduce matrices and solve for displacement. 
        /// </summary>
        /// <returns> List of nodal displacement. </returns>
        private CSD.DenseMatrix CalculateDisplacementCSparse(LA.Matrix<double> K_gl, CSD.DenseMatrix R_gl, List<List<int>> applyBCToDOF, ref List<string> info)
        {
            var timer = new System.Diagnostics.Stopwatch();

            // Make list of boundary condistions
            info.Add("--- Displacement calculations ---");
            timer.Start();
            List<int> BCList = new List<int>();
            for (int i = 0; i < applyBCToDOF.Count; i++)
            {
                for (int j = 0; j < applyBCToDOF[0].Count; j++)
                {
                    BCList.Add(applyBCToDOF[i][j]); // list of 0 and 1 values for boundary condition for dof; true = 1, false = 0
                }
            }

            // Apply boundary conditions to movement
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

            
            CompressedColumnStorage<double> CCS = CSD.SparseMatrix.OfColumnMajor(K_gl.RowCount, K_gl.ColumnCount, K_gl.AsColumnMajorArray());

            #region Testing problems
            var R_LA = LA.Double.DenseVector.Build.DenseOfArray(R_gl.Values);

            var u_LA = K_gl.Solve(R_LA);
            #endregion


            SparseLU CS_K_global = SparseLU.Create(CCS, ColumnOrdering.MinimumDegreeAtPlusA, 0.0);
            //double[] CS_u = CSD.Vector.Create(K_global_red.RowCount * 1, 0.0);
            double[] CS_u = CSD.Vector.Create(K_gl.RowCount * 1, 0.0);
            //double[] CS_R = R_red.Column(0).ToArray();
            double[] CS_R = R_gl.Column(0).ToArray();
            for (int i = 0; i < CS_R.Length; i++)
            {
                CS_R[i] = Math.Round(CS_R[i], 6);
            }

            CS_K_global.Solve(CS_R, CS_u);
            

            var u = new CSD.DenseMatrix(CS_u.Length, 1, CS_u);
            //var u = new CSD.DenseMatrix(CS_u.Length, 1, u_LA.ToArray());

            //LA.Matrix<double> u = LA.Double.DenseMatrix.OfColumnArrays(CS_u);

            return u;


        }
        private LA.Matrix<double> CalculateDisplacement( LA.Matrix<double> K_gl, LA.Matrix<double> R_gl, List<List<int>> applyBCToDOF, ref List<string> info)
        {
            // summary: include boundary condistions and calculate global displacement
            var timer = new System.Diagnostics.Stopwatch();

            // Make list of boundary condistions
            info.Add("--- Displacement calculations ---");
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
            info.Add($"Time elapsed during boundary conditions: {timer.ElapsedMilliseconds} ms"); timer.Reset();

            // Reduce K_global and R
            timer.Start();

            // -- EDIT SVERRE ---
            //int numRows = K_gl.RowCount;
            //ReduceMatrices(BCList, numRows, ref K_gl,ref R_gl);

            //var reducedData = ReduceKandR(K_gl, R_gl, BCList);
            
            //LA.Matrix<double> K_global_red = reducedData.Item1;
            //LA.Matrix<double> R_red = reducedData.Item2;
            timer.Stop();
            info.Add($"Time elapse for reducing K and R: {timer.ElapsedMilliseconds} ms"); timer.Reset();
            
            // Time recorder
            var sw0 = new System.Diagnostics.Stopwatch();
            timer.Start();
            // Mathnet.Numerics to CSparse
            //var CMA = K_global_red.Storage.ToColumnMajorArray();
            var CMA = K_gl.Storage.ToColumnMajorArray();
            //CompressedColumnStorage<double> CCS = CSD.SparseMatrix.OfColumnMajor(K_global_red.RowCount, K_global_red.ColumnCount, CMA);
            CompressedColumnStorage<double> CCS = CSD.SparseMatrix.OfColumnMajor(K_gl.RowCount, K_gl.ColumnCount, CMA);
            timer.Stop();
            info.Add($"Convert MathNet to CSparse: {timer.ElapsedMilliseconds} ms"); timer.Reset();

            
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
            info.Add($"Solve the system for displacements: {timer.ElapsedMilliseconds} ms"); timer.Reset();
            //Rhino.RhinoApp.WriteLine($"### {K_global_red.RowCount} x {K_global_red.ColumnCount} Matrix. CSparse Elapsed [msec] = {sw0.Elapsed.TotalMilliseconds}");
            Rhino.RhinoApp.WriteLine($"### {K_gl.RowCount} x {K_gl.ColumnCount} Matrix. CSparse Elapsed [msec] = {sw0.ElapsedMilliseconds} "); sw0.Reset();

            // CSparse to Mathnet.Numerics
            timer.Start();
            LA.Matrix<double> u = LA.Double.DenseMatrix.OfColumnArrays(CS_u);
            timer.Stop();
            info.Add($"Convert back to MathNet: {timer.ElapsedMilliseconds} ms"); timer.Reset();
            info.Add("--- End displacement calculations ---");
            
            
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
        /// Reduce stiffness matrix and load vector for fixed boundary conditions.
        /// </summary>
        /// <returns> Reduced stiffness matrix and load vector. </returns>
        private Tuple<LA.Matrix<double>, LA.Matrix<double>> ReduceKandR(LA.Matrix<double> K_global,LA.Matrix<double> R, List<int> BC)
        {
            int removeIndex = 0;
            for (int i = 0; i < K_global.RowCount; i++)
            {
                if (BC[i] == 1)
                {
                    K_global = K_global.RemoveRow(removeIndex);
                    K_global = K_global.RemoveColumn(removeIndex);
                    R = R.RemoveRow(removeIndex);
                    removeIndex--;
                }
                removeIndex++;
            }
            return Tuple.Create(K_global, R);
        }
        private void ReduceMatrices(List<int> BC, int rowCount, ref LA.Matrix<double> K, ref LA.Matrix<double> R)
        {
            int subtract = 0;
            for (int i = 0; i < BC.Count; i++)
            {
                if (BC[i] == 1)
                {
                    K = K.RemoveColumn(i - subtract);
                    K = K.RemoveRow(i - subtract);
                    
                    R = R.RemoveRow(i - subtract);
                    subtract--;
                }
            }


        }

        /// <summary>
        /// Calculate a list of strain and stress vectors for each node in a element.
        /// </summary>
        /// <returns> Strain and stress vectors for each node in a element. </returns>
        private Tuple<LA.Matrix<double>, LA.Matrix<double>> CalculateElementStrainStress(Element element, LA.Matrix<double> u, Material material)
        {
            LA.Matrix<double> C = material.GetMaterialConstant();

            
            List<LA.Matrix<double>> B_local = CalculateElementMatrices(element, material).Item2; // this can be changed to save time.. No need to establish the stiffness matrix of an element for this
            LA.Matrix<double> elementGaussStrain = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            LA.Matrix<double> elementGaussStress = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            LA.Matrix<double> elementStrain = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            LA.Matrix<double> elementStress = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, element.Nodes.Count);
            LA.Matrix<double> localDeformation = LA.Double.DenseMatrix.Build.Dense(3 * element.Nodes.Count, 1);

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
                    LA.Matrix<double> gaussStress = C.Multiply(B_local[n]).Multiply(localDeformation);

                    for (int i = 0; i < B_local[0].RowCount; i++)
                    {
                        elementGaussStrain[i, n] = gaussStrain[i, 0];
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
                
                for(int i = 0; i < element.Nodes.Count; i++)
                {
                    
                    //Get deformation at node
                    LA.Matrix<double> nodalDeformation = LA.Double.DenseMatrix.Build.Dense(3 , 1);
                    LA.Matrix<double> nodeStrain = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, 1);
                    LA.Matrix<double> nodeStress = LA.Double.DenseMatrix.Build.Dense(B_local[0].RowCount, 1);
                    

                    nodeStrain = (B_local[0]).Multiply(localDeformation);
                    nodeStress.Add(C.Multiply(nodeStrain));

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

        private void ElementStrainStressCSparse(Element element, CSD.DenseMatrix u, Material material)
        {

        }
        /*
        private void GlobalStressesCSparse(List<Element> elements, CSD.DenseMatrix u, Material material, out CSD.DenseMatrix GlobalStress, out CSD.DenseMatrix Mises, out CSD.DenseMatrix GlobalMises)
        {
            int numNodes = u.RowCount;
            int stressRowDim = 6;

            CSD.DenseMatrix stress = new CSD.DenseMatrix(stressRowDim, numNodes);
            CSD.DenseMatrix counter = new CSD.DenseMatrix(stressRowDim, numNodes);

            List<CSD.DenseMatrix> elementStressList = new List<CSD.DenseMatrix>();

            // iterate through elements
            foreach (Element el in elements)
            {

            }



        }
        */
        private Tuple<LA.Matrix<double>, LA.Vector<double>, LA.Vector<double>> CalculateGlobalStress(List<Element> elements, LA.Matrix<double> u, Material material)
        {
            int numNodes = u.RowCount / 3;
            int stressRowDim = 6;
            LA.Matrix<double> globalStress = LA.Double.DenseMatrix.Build.Dense(stressRowDim, numNodes);
            LA.Matrix<double> counter = LA.Double.DenseMatrix.Build.Dense(stressRowDim, numNodes);
            List<LA.Matrix<double>> elementStressList = new List<LA.Matrix<double>>();
            foreach (Element element in elements)
            {
                LA.Matrix<double> elementStress = CalculateElementStrainStress(element, u, material).Item2;

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
            LA.Vector<double> mises = DenseVector.Build.Dense(numNodes);
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
            LA.Vector<double> elementMises = DenseVector.Build.Dense(elements.Count);
            for (int i = 0; i < elementStressList.Count; i++)
            {
                for (int j = 0; j < elements[0].Nodes.Count; j++)
                {
                    LA.Vector<double> nodeStress = elementStressList[i].Column(j);
                    double Sxx = nodeStress[0];
                    double Syy = nodeStress[1];
                    double Szz = nodeStress[2];
                    double Sxy = nodeStress[3];
                    double Sxz = nodeStress[4];
                    double Syz = nodeStress[5];
                    elementMises[i] += Math.Sqrt(0.5 * (Math.Pow(Sxx - Syy, 2) + Math.Pow(Syy - Szz, 2) + Math.Pow(Szz - Sxx, 2)) + 3 * (Math.Pow(Sxy, 2) + Math.Pow(Sxz, 2) + Math.Pow(Syz, 2)));
                }
                elementMises[i] = elementMises[i] / (double)8; // get average of nodal mises
            }


            return Tuple.Create(globalStress, mises, elementMises);
        }

        /// <summary>
        /// Color mesh after nodal stress values.
        /// </summary>
        private void ColorMeshAfterStress(SmartMesh mesh, LA.Vector<double> mises, Material material)
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
        #endregion

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                //return Resources.IconForThisComponent;
                //return FEMLoad;
                return Properties.Resources.FESolver;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("b1911217-5eeb-4722-8226-ff72e3b8fbb0"); }
        }
    }
}