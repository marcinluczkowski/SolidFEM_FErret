using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;

namespace SolidFEM.FiniteElementMethod
{
    public class FEM_BoundaryOnPoints_MESH : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEM_BoundaryOnPoints class.
        /// </summary>
        public FEM_BoundaryOnPoints_MESH()
          : base("FEM_BoundaryOnPoints", "Nickname",
              "Description",
              "SmartMesh", "FEM-Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh", "mesh", "The mesh to apply boundary conditions on", GH_ParamAccess.list); // 0
            pManager.AddPointParameter("SupportPoints", "pt", "Position to place boundary points", GH_ParamAccess.list); // 1
            pManager.AddSurfaceParameter("SupportSurface", "Surf", "", GH_ParamAccess.list); //2
            pManager.AddIntegerParameter("Type of support", "T", "Type of support (1 = points, 2 = surface)", GH_ParamAccess.item, 1); //3
            pManager.AddBooleanParameter("Tx", "", "", GH_ParamAccess.item, true); // 4
            pManager.AddBooleanParameter("Ty", "", "", GH_ParamAccess.item, true); // 5
            pManager.AddBooleanParameter("Ty", "", "", GH_ParamAccess.item, true); // 6
            

            pManager[1].Optional = true;
            pManager[2].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Supports", "sup", "Support to apply to model", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // -- variables -- 
            bool tx = true; // 2
            bool ty = true; // 3
            bool tz = true; // 4
            List<Mesh> meshList = new List<Mesh>(); // 0
            //SmartMesh sMesh = new SmartMesh(); // 0
            List<Point3d> positions = new List<Point3d>(); // 1
            List<Brep> supSurf = new List<Brep>();
            int type = 1;


            // - input --

            //if (!DA.GetDataList(0, meshList)) return;
            DA.GetDataList(0, meshList);
            DA.GetDataList(1, positions);
            DA.GetDataList(2, supSurf);
            DA.GetData(3, ref type);
            DA.GetData(4, ref tx);
            DA.GetData(5, ref ty);
            DA.GetData(6, ref tz);
            


            // -- method --
            // clean the mesh and sort nodes, if Hex20 add midside nodes
            var newMeshList = new List<Mesh>();
            int c = 0; // delete after testing
            string elementType = "Tet10";

            foreach (Mesh mesh in meshList)
            {
                if (mesh.Vertices.Count == 8)
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
                else if (elementType == "Tet10")
                {
                    Mesh nM = FEM_Utility.AddMidEdgeNodes(mesh);
                    newMeshList.Add(nM);
                }

                else if(mesh.Vertices.Count == 20)
                {
                        newMeshList.Add(mesh);
                }
                else
                {
                    newMeshList.Add(mesh);
                }


            }

            List<Point3d> nodePts = FEM_Utility.GetMeshNodes(newMeshList);
            List<Support> supportList = new List<Support>();

            if (type == 1)
            {

                foreach (var pt in positions)
                {
                    int nodeIndex = GetClosestNodeIndex(nodePts, pt);
                    if (nodeIndex != -1)
                    {
                        Support sup = new Support(nodePts[nodeIndex], tx, ty, tz);
                        supportList.Add(sup);
                    }

                }
                
            }
            else if (type == 2)
            {
                foreach (Brep surf in supSurf)
                {
                    Surface supSurface = surf.Surfaces[0];
                    for (int i = 0; i < nodePts.Count; i++)
                    {
                        Point3d point = nodePts[i];
                        double u = new double();
                        double v = new double();
                        supSurface.ClosestPoint(point, out u, out v);

                        Point3d surfPt = supSurface.PointAt(u, v);

                        double tol = 0.001;

                        if (surfPt.DistanceTo(point) < tol)
                        {
                            Support sup = new Support(point, tx, ty, tz);
                            supportList.Add(sup);
                        }
                    }
                }
            }
            
            
            

            // -- output --
            DA.SetDataList(0, supportList);

        }
        private int GetClosestNodeIndex(List<Point3d> nodePoints, Point3d loadPosition)
        {
            int nodeIndex = -1;

            for (int i = 0; i < nodePoints.Count; i++)
            {
                double distance = nodePoints[i].DistanceToSquared(loadPosition);
                if (distance < 0.00001) // set tolerance
                {
                    nodeIndex = i;
                    break;
                }
            }
            return nodeIndex;
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                //return Resources.IconForThisComponent;
                return Properties.Resources.FESupport;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("6f64f9c3-3741-4655-ba50-9d0b2f221844"); }
        }
    }
}