using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;

namespace SolidFEM.FiniteElementMethod
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
            pManager.AddNumberParameter("Young modulus (isotropic)", "E", "Young modulus [MPa]. Default value: 210000 MPa.", GH_ParamAccess.item, 210000);
            pManager.AddNumberParameter("Young modulus xx-dir (orthotropic)", "Exx", "Young modulus [MPa] in the xx-direction for orthotropic materials", GH_ParamAccess.item, 0);
            pManager.AddNumberParameter("Young modulus yy-dir (orthotropic)", "Eyy", "Young modulus [MPa] in the yy-direction for orthotropic materials", GH_ParamAccess.item, 0);
            pManager.AddNumberParameter("Young modulus zz-dir (orthotropic)", "Ezz", "Young modulus [MPa] in the zz-direction for orthotropic materials", GH_ParamAccess.item, 0);
            pManager.AddNumberParameter("Possion Ratio", "nu", "Possion Ratio [-]. Default value: 0.3.", GH_ParamAccess.item, 0.3);
            pManager.AddNumberParameter("Yielding stress", "fy", "Yield stress [MPa]. Default value: 355 MPa.", GH_ParamAccess.item, 355);
            pManager.AddNumberParameter("Material weight", "rho", "The weight of the material [kg/m^3]", GH_ParamAccess.item, 7850);

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
            double Emodul = 0.0;
            double Exx = 0.0;
            double Eyy = 0.0;
            double Ezz = 0.0;
            double nu = 0.0;
            double fy = 0.0;
            double weight = 0.0;

            DA.GetData(0, ref Emodul);
            DA.GetData(1, ref Exx);
            DA.GetData(2, ref Eyy);
            DA.GetData(3, ref Ezz);
            DA.GetData(4, ref nu);
            DA.GetData(5, ref fy);
            DA.GetData(6, ref weight);


            // Code
            Material material = new Material(Emodul, nu, fy);
            material.Weight = weight;

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
                return Properties.Resources.FEMaterial;// return Properties.Resources.Icon_Material;
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