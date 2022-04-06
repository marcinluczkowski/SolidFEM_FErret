using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;

namespace SolidFEM.FiniteElementMethod
{
    public class FEMMaterial_Orthotropic : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMMaterial class.
        /// </summary>
        public FEMMaterial_Orthotropic()
          : base("FEM Material - Orthotropic", "Orthotropic Material",
              "Create orthotropic material for the FEM Solver.",
              "SmartMesh", "FEM")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Young modulus in xx-direction", "Exx", "Young modulus in xx-direction [MPa].", GH_ParamAccess.item);
            pManager.AddNumberParameter("Young modulus in yy-direction", "Eyy", "Young modulus in yy-direction [MPa].", GH_ParamAccess.item);
            pManager.AddNumberParameter("Young modulus in zz-direction", "Ezz", "Young modulus in zz-direction [MPa].", GH_ParamAccess.item);
            pManager.AddNumberParameter("Possion Ratio in xy-direction", "nu_xy", "Possion Ratio in xy-direction [-]. Default value 0.", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Possion Ratio in yz-direction", "nu_yz", "Possion Ratio in yz-direction [-]. Default value 0.", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Possion Ratio in zx-direction", "nu_zx", "Possion Ratio in zx-direction [-]. Default value 0.", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Shear modulus in xy-direction", "Gxy", "Shear modulus in xy-direction [MPa].", GH_ParamAccess.item);
            pManager.AddNumberParameter("Shear modulus in zx-direction", "Gzx", "Shear modulus in zx-direction [MPa].", GH_ParamAccess.item);
            pManager.AddNumberParameter("Shear modulus in yz-direction", "Gyz", "Shear modulus in yz-direction [MPa].", GH_ParamAccess.item);
            pManager.AddNumberParameter("Yielding stress", "fy", "Yield stress [MPa]. Default value: 355 MPa.", GH_ParamAccess.item);
            pManager.AddNumberParameter("Material weight", "rho", "The weight of the material [kg/m^3]", GH_ParamAccess.item);

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
            double Exx = 0.0;
            double Eyy = 0.0;
            double Ezz = 0.0;
            double nu_xy = 0.0;
            double nu_yz = 0.0;
            double nu_zx = 0.0;
            double Gxy = 0.0;
            double Gzx = 0.0;
            double Gyz = 0.0;
            double fy = 0.0;
            double weight = 0.0;
            DA.GetData(0, ref Exx);
            DA.GetData(1, ref Eyy);
            DA.GetData(2, ref Ezz);
            DA.GetData(3, ref nu_xy);
            DA.GetData(4, ref nu_yz);
            DA.GetData(5, ref nu_zx);
            DA.GetData(6, ref Gxy);
            DA.GetData(7, ref Gzx);
            DA.GetData(8, ref Gyz);
            DA.GetData(9, ref fy);
            DA.GetData(10, ref weight);


            // Code
            Material material = new Material(Exx, Eyy, Ezz, nu_xy, nu_yz, nu_zx, Gxy, Gzx, Gyz, fy);
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
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("05ED8A20-FFD0-4A93-95B7-7A13BE8FAB16"); }
        }
    }
}