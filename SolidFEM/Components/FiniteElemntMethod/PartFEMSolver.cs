using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Linq;
using ClosedXML.Excel;

// Csparse
using LA = MathNet.Numerics.LinearAlgebra;
using CSparse;
using CSD = CSparse.Double;
using CSparse.Double.Factorization;
using CSparse.Storage;
using Grasshopper;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;


namespace SolidFEM.FiniteElementMethod
{
    public class PartFEMSolver : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMsolver class.
        /// </summary>
        public PartFEMSolver()
          : base("Part FEM Solver", "Part Solver",
              "Part Solver for FEM problems with Rhino Mesh" +
                "Uses 3 translation DOFS pr node, linear shape functions, two Gauss Points and full integration.",
              "SolidFEM", "FEM-Mesh")
        {
            //Logger = new FEMLogger();
        }

        // -- logger for component -- 
        private static FEMLogger Logger;


        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh", "m", "Mesh List", GH_ParamAccess.list); // 0
            pManager.AddIntegerParameter("Element order", "elOrder",
                "The order of the finite element to be used. Currently working for 1st and 2nd order element.",
                GH_ParamAccess.item, 1); // 1
            pManager.AddGenericParameter("Loads", "load", "Load vector from FEM Load.", GH_ParamAccess.list); // 2
            pManager.AddGenericParameter("Boundary Conditions", "BC", "Boundary conditions from FEM Boundary Condtion.", GH_ParamAccess.list); // 3
            pManager.AddGenericParameter("Material", "Mat", "Material for each element.", GH_ParamAccess.list); // 4
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
            pManager.AddGenericParameter("FE_Mesh", "femesh", "The FE_Mesh containing results from the analysis.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Points", "pt", "Nodes in mesh", GH_ParamAccess.list);
            pManager.AddTextParameter("Global K", "", "", GH_ParamAccess.list);
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

            Logger = new FEMLogger();
            Logger.AddInfo("Component started.");

            List<Mesh> meshList = new List<Mesh>();
            int degree = 1;
            List<double> loads = new List<double>();
            List<Support> supports = new List<Support>();
            List<Material> material = new List<Material>();

            DA.GetDataList(0, meshList); // 0
            DA.GetData(1, ref degree); // 1
            DA.GetDataList(2, loads); // 2
            DA.GetDataList(3, supports); // 3
            DA.GetDataList(4, material); // 4


            List<string> info = new List<string>();


            // 0. Initial step

            // clean the mesh and sort nodes
            var newMeshList = new List<Mesh>();
            int c = 0; // delete after testing
            string elementType = "";


            foreach (Mesh mesh in meshList)
            {
                if (mesh.Vertices.Count == 8 && degree == 1)
                {
                    Mesh nM = GrahamScan.DoGrahamScan(mesh);

                    if (nM.IsValid)
                    {
                        newMeshList.Add(nM);
                    }
                    else
                    {
                        newMeshList.Add(mesh);
                    }
                    c++;
                }
                else if (mesh.Vertices.Count == 8 && degree == 2)
                {
                    Mesh nM = GrahamScan.DoGrahamScan(mesh);
                    if (nM.IsValid)
                    {
                        Mesh hex20 = FEM_Utility.AddMidEdgeNodes(nM); // Modify this method to make it apply for HEX elements
                        newMeshList.Add(hex20);
                    }
                    else
                    {
                        Mesh hex20 = FEM_Utility.AddMidEdgeNodes(mesh); // Modify this method to make it apply for HEX elements
                        newMeshList.Add(hex20);
                    }
                }
                /*
                else if (mesh.Vertices.Count == 20)
                {
                    newMeshList.Add(mesh);
                }*/
                else if (mesh.Vertices.Count == 4 && degree == 2)
                {
                    Mesh nM = FEM_Utility.AddMidEdgeNodes(mesh);
                    newMeshList.Add(nM);
                }
                else // if tet4
                {
                    newMeshList.Add(mesh);
                }

            }

            List<Element> elements;
            List<Node> femNodes = new List<Node>();

            List<Point3d> nodePos = FEM_Utility.GetMeshNodes(newMeshList); // This needs to be modified? Or can I keep it? Then use the old mesh later? Make sure the element have correct names based on number of vertices.

            for (int i = 0; i < nodePos.Count; i++)
            {
                Node node = new Node(i, nodePos[i]);
                femNodes.Add(node);
            }



            int numNodes = nodePos.Count;
            FEM_Utility.ElementsFromMeshList(newMeshList, nodePos, out elements); // Control if this is working properly. 


            // 1. Get global stiffness matrix
            watch.Start();
            var K_globalC = FEM_Matrices_Part.GlobalStiffnessCSparse(ref elements, numNodes, material, ref Logger);

            //Delete after testing, output K matrix in gh

            LA.Matrix<double> kglobal = LA.Matrix<double>.Build.DenseOfArray(K_globalC);
            info.Add(kglobal.ToString());

            DA.SetDataList(8, info);


            watch.Stop();
            Logger.AddInfo($"Time used calculating global stiffness matrix: {watch.ElapsedMilliseconds} ms"); watch.Reset();

            // 3. Get load vector
            watch.Start();


            // self weight
            var selfWeight = FEM_Utility_Part.GetBodyForceVector(material, elements, numNodes, Logger); //Find method to find body force with different materials

            double weight = 0;

            for (int i = 0; i < selfWeight.Count; i++)
            {
                weight += selfWeight[i];
            }


            CSD.DenseMatrix R_self = new CSD.DenseMatrix(numNodes * 3, 1, selfWeight.ToArray());
            CSD.DenseMatrix R_external = new CSD.DenseMatrix(numNodes * 3, 1);

            for (int i = 0; i < loads.Count; i++)
            {
                R_external[i, 0] = loads[i];
            }

            CSD.DenseMatrix R = (CSD.DenseMatrix)R_self.Add(R_external);
            //CSD.DenseMatrix R = R_external;

            watch.Stop();
            Logger.AddInfo($"Time used to establish global load vector: {watch.ElapsedMilliseconds} ms"); watch.Reset();

            // 5. Fix BoundaryConditions
            watch.Start();
            List<List<int>> boundaryConditions = FixBoundaryConditionsSverre(supports, nodePos);


            Logger.AddInfo($"Time used on boundary conditions: {watch.ElapsedMilliseconds} ms"); watch.Reset();



            // 6. Calculate displacement 
            watch.Start();
            CSD.DenseMatrix u_CSparse = FEM_Utility.CalculateDisplacementCSparse(K_globalC, R, boundaryConditions, ref Logger);
            Logger.AddInfo($"Time used on displacement calculations with CSparce: {watch.ElapsedMilliseconds} ms"); watch.Reset();

            // 7. Calculate stress

            // convert CSparse to double and MathNet
            double[] u_val = u_CSparse.Values;
            var uCS2MN = new LA.Double.DenseMatrix(u_val.Length, 1, u_val);

            watch.Start();
            var stress = FEM_Utility_Part.CalculateGlobalStress(elements, uCS2MN, material, ref Logger); // make this compatible with the CSparse matrix as well.
            Logger.AddInfo($"Time used on stress calculations: {watch.ElapsedMilliseconds} ms"); watch.Reset();

            LA.Matrix<double> globalStress = stress.Item1;
            LA.Vector<double> mises = stress.Item2;
            LA.Vector<double> misesElement = stress.Item3;
            watch.Start();
            Logger.AddInfo($"Time used on mesh colouring: {watch.ElapsedMilliseconds} ms"); watch.Reset();

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
            Logger.AddInfo($"Time used on output preparation: {watch.ElapsedMilliseconds} ms"); watch.Reset();


            // Write logger to .txt file
            Logger.LogToTextFile();




            // TEMPORARY WHILE WAITING FOR MARCIN's MESH CLASS  ////MATERIAL IS SET TO MATERIAL[0]!!
            TempFE_Mesh outMesh = new TempFE_Mesh(meshList, femNodes, elements, elementMises, nodalMises, u1, u2, u3, material[0],
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
            DA.SetDataList(5, Logger.LogList);
            DA.SetData(6, outMesh);
            //DA.SetDataList(7, elements);
            DA.SetDataList(7, nodePos);

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
        private List<List<int>>
            FixBoundaryConditionsSverre(List<Support> sups, List<Point3d> nodesPos)
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
                    if (sup.Ty == true) { bc[1] = 1; }
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
            get { return new Guid("8817A220-EDC9-4275-B328-23E32A98766B"); }
        }
    }
}