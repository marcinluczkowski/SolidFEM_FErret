using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints.FiniteElementMethod
{
    public class FEMMaterial : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMMaterial class.
        /// </summary>
        public FEMMaterial()
          : base("FEM Material", "Material",
              "Create material for the FEM Solver.",
              "SmartMesh", "FEM")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Young modulus", "E", "Young modulus [MPa]. Default value: 210000 MPa.", GH_ParamAccess.item, 210000);
            pManager.AddNumberParameter("Possion Ratio", "nu", "Possion Ratio [-]. Default value: 0.3.", GH_ParamAccess.item, 0.3);
            pManager.AddNumberParameter("Yielding stress", "fy", "Yield stress [MPa]. Default value: 355 MPa.", GH_ParamAccess.item, 355);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Material", "mat", "Material.", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            double Emodul = 210000;
            double nu = 0.3;
            double fy = 355;
            DA.GetData(0, ref Emodul);
            DA.GetData(1, ref nu);
            DA.GetData(2, ref fy);


            // Code
            Material material = new Material(Emodul, nu, fy);

            // Output
            DA.SetData(0, material);
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
                return Properties.Resources.Icon_Material;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("6198afae-8cd5-4b28-a841-5b1a48de93ec"); }
        }
    }
}