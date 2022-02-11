using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;

namespace SolidFEM.FiniteElementMethod
{
    public class FEM_BoundaryOnPoints : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEM_BoundaryOnPoints class.
        /// </summary>
        public FEM_BoundaryOnPoints()
          : base("FEM_BoundaryOnPoints", "Nickname",
              "Description",
              "SmartMesh", "Support")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "mesh", "The mesh to apply boundary conditions to", GH_ParamAccess.item); // 0
            pManager.AddPointParameter("Point", "pt", "Position to place boundary points", GH_ParamAccess.item); // 1
            pManager.AddBooleanParameter("Tx", "", "", GH_ParamAccess.item, true); // 2
            pManager.AddBooleanParameter("Ty", "", "", GH_ParamAccess.item, true); // 3
            pManager.AddBooleanParameter("Ty", "", "", GH_ParamAccess.item, true); // 4
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Supports", "sup", "Support to apply to model", GH_ParamAccess.item);
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
            SmartMesh sMesh = new SmartMesh(); // 0
            Point3d pos = new Point3d(); // 1

            // - input --
            if (!DA.GetData(0, ref sMesh)) return;
            if (!DA.GetData(1, ref pos)) return;
            DA.GetData(2, ref tx);
            DA.GetData(3, ref ty);
            DA.GetData(4, ref tz);


            // -- method --

            Support sup = new Support(pos, tx, ty, tz);


            // -- output --
            DA.SetData(0, sup);

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
            get { return new Guid("46459817-bc55-4b8e-9045-b098663a89cd"); }
        }
    }
}