using System;
using System.Collections.Generic;
using SolidFEM.Classes;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace SolidFEM.Deconstructors
{
    public class NodalResults : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public NodalResults()
          : base("NodalResults", "nodeRes",
              "Preview nodal stresses and displacements. ",
              "SmartMesh", "Utilities")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("FEMesh", "feMesh", "An instance of a solid mesh", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Points", "pts", "All global node points", GH_ParamAccess.list); // 0
            pManager.AddNumberParameter("U1", "disp1", "Displacment in U1", GH_ParamAccess.list); // 1
            pManager.AddNumberParameter("U2", "disp2", "Displacment in U2", GH_ParamAccess.list); // 2
            pManager.AddNumberParameter("U3", "disp3", "Displacment in U3", GH_ParamAccess.list); // 3
            pManager.AddNumberParameter("Mises", "m", "Nodal mises stresses", GH_ParamAccess.list); // 4
            pManager.AddNumberParameter("xx", "xx", "Nodal stress in xx", GH_ParamAccess.list); // 5
            pManager.AddNumberParameter("yy", "yy", "Nodal stress in yy", GH_ParamAccess.list); // 6
            pManager.AddNumberParameter("zz", "zz", "Nodal stress in zz", GH_ParamAccess.list); // 7
            pManager.AddNumberParameter("xy", "xy", "Nodal stress in xy", GH_ParamAccess.list); // 8
            pManager.AddNumberParameter("xz", "xz", "Nodal stress in xz", GH_ParamAccess.list); // 9
            pManager.AddNumberParameter("yz", "yz", "Nodal stress in yz", GH_ParamAccess.list); // 10



        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            TempFE_Mesh mesh = new TempFE_Mesh();

            if (!DA.GetData(0, ref mesh)) return;

            List<Point3d> points = new List<Point3d>();
            foreach (Node node in mesh.MeshNodes)
            {
                points.Add(node.Coordinate);
            }

           // -- output -- 
           DA.SetDataList(0, points); // 0
           DA.SetDataList(1, mesh.dU); // 1
           DA.SetDataList(2, mesh.dV); // 2
           DA.SetDataList(3, mesh.dW); // 3
           DA.SetDataList(4, mesh.NodelMisesStresses); // 4
           DA.SetDataList(5, mesh.Sigma_xx); // 5
           DA.SetDataList(6, mesh.Sigma_yy); // 6
           DA.SetDataList(7, mesh.Sigma_zz); // 7
           DA.SetDataList(8, mesh.Sigma_xy); // 8
           DA.SetDataList(9, mesh.Sigma_xz); // 9
           DA.SetDataList(10, mesh.Sigma_yz); // 9
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("a7f47096-f22b-4a6c-a277-4550c0a76073"); }
        }
    }
}