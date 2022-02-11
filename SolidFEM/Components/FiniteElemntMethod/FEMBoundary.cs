using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;

namespace SolidFEM.FiniteElementMethod
{
    public class FEMBC : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMBoundaryCondition class.
        /// </summary>
        public FEMBC()
          : base("FEM Boundary Condtion", "BC",
              "Create boundary condition for the FEM Solver.",
              "SmartMesh", "FEM")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "SM", "SmartMesh.", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Faces", "faces", "Index of faces with boundary conditions.", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Edges", "edges", "Index of edges with boundary conditions.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Tx", "Tx", "True if no translation in x-direction.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Ty", "Ty", "True if no translation in y-direction.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Tz", "Tz", "True if no translation in z-direction.", GH_ParamAccess.item);

            pManager[1].Optional = true;
            pManager[2].Optional = true;
            pManager[3].Optional = true;
            pManager[4].Optional = true;
            pManager[5].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Boundary conditions", "BC", "List of DOFS that are fixed.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Points", "pts", "List of points with boundary conditions.", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input

            SmartMesh smartMesh = new SmartMesh();
            List<int> indicesOfFaceWithBC = new List<int>();
            List<int> indicesOfEdgeWithBC = new List<int>();
            bool Tx = false;
            bool Ty = false;
            bool Tz = false;

            DA.GetData(0, ref smartMesh);
            DA.GetDataList(1, indicesOfFaceWithBC);
            DA.GetDataList(2, indicesOfEdgeWithBC);
            DA.GetData(3, ref Tx);
            DA.GetData(4, ref Ty);
            DA.GetData(5, ref Tz);

            // Code

            if (!DA.GetData(0, ref smartMesh)) { return; }
            List<bool> applyBC = new List<bool>() { Tx, Ty, Tz };
            List<List<int>> applyBCToDOF = new List<List<int>>();
            bool nodeIsOnGeometry = false;
            List<Point3d> pointsWithBC = new List<Point3d>();

            // Loop each dof for each node
            for (int i = 0; i < smartMesh.Nodes.Count; i++)
            {
                Node node = smartMesh.Nodes[i];

                if (indicesOfEdgeWithBC.Count > 0) // get nodes on relevant edges
                {
                    for (int n = 0; n < indicesOfEdgeWithBC.Count; n++)
                    {
                        BrepEdge edge = smartMesh.Geometry.Edges[indicesOfEdgeWithBC[n]];
                        if (IsPointOnEdge(node.Coordinate, edge))
                        {
                            nodeIsOnGeometry = true;
                            pointsWithBC.Add(node.Coordinate);
                            break;
                        }
                    }
                }
                else // get nodes on relevant faces
                {
                    for (int n = 0; n < indicesOfFaceWithBC.Count; n++)
                    {
                        BrepFace face = smartMesh.Geometry.Faces[indicesOfFaceWithBC[n]];
                        if (IsPointOnFace(node.Coordinate, face))
                        {
                            nodeIsOnGeometry = true;
                            pointsWithBC.Add(node.Coordinate);
                            break;
                        }
                    }
                }

                // Loop each DOF of node. Add 1 if BC and on geometry, else 0
                List<int> BCNode = new List<int>();
                for (int j = 0; j < 3; j++)
                {
                    if (nodeIsOnGeometry & applyBC[j] == true)
                    {
                        BCNode.Add(1);
                        continue;
                    }
                    BCNode.Add(0); // add information of a single dof of the node
                }
                applyBCToDOF.Add(BCNode); // add information of the node
                nodeIsOnGeometry = false; // reset
            }

            // Output

            DA.SetDataList(0, applyBCToDOF);
            DA.SetDataList(1, pointsWithBC);
        }

        #region Methods
        /// <summary>
        /// Check if point is on a BrepFace.
        /// </summary>
        /// <returns> Boolean value.</returns>
        private bool IsPointOnFace(Point3d point, BrepFace face)
        {
            bool nodeIsOnGeometry = false;
            face.ClosestPoint(point, out double PointOnCurveU, out double PointOnCurveV);
            Point3d testPoint = face.PointAt(PointOnCurveU, PointOnCurveV);  // make test point 
            double distanceToFace = testPoint.DistanceTo(point); // calculate distance between testPoint and node
            if (distanceToFace <= 0.0001 & distanceToFace >= -0.0001) // if distance = 0: node is on edge
            {
                nodeIsOnGeometry = true;
            }
            return nodeIsOnGeometry;
        }

        /// <summary>
        /// Check if point is on a BrepEdge.
        /// </summary>
        /// <returns> Boolean value.</returns>
        private bool IsPointOnEdge(Point3d point, BrepEdge edge)
        {
            bool nodeIsOnGeometry = false;
            edge.ClosestPoint(point, out double PointOnCurve);
            Point3d testPoint = edge.PointAt(PointOnCurve);  // make test point 
            double distanceToEdge = testPoint.DistanceTo(point); // calculate distance between testPoint and node
            if (distanceToEdge <= 0.0001 & distanceToEdge >= -0.0001) // if distance = 0: node is on edge
            {
                nodeIsOnGeometry = true;
            }
            return nodeIsOnGeometry;
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
                // return Resources.IconForThisComponent;
                return null;//  return Properties.Resources.Icon_BC;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("7268f076-dbf3-4dbd-a363-41cd9966ad96"); }
        }
    }
}