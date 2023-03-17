using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using SolidFEM.Classes;

namespace SolidFEM.Components.FiniteElemntMethod
{
    public class FEMMaterialOrto : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMMaterialOrto class.
        /// </summary>
        public FEMMaterialOrto()
          : base("FEM Material Orto", "MaterialOrto",
              "Create orthotropic material for the FEM Solver.",
              "SolidFEM", "FEM")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Young modulus E1", "E1", "Young modulus [MPa]. Default value: 10 000 MPa.", GH_ParamAccess.item, 10000);
            pManager.AddNumberParameter("Young modulus E2", "E2", "Young modulus [MPa]. Default value: 800 MPa.", GH_ParamAccess.item, 800);
            pManager.AddNumberParameter("Young modulus E3", "E3", "Young modulus [MPa]. Default value: 400 MPa.", GH_ParamAccess.item, 400);
            pManager.AddNumberParameter("Shear modulus G1", "G1", "Shear modulus [MPa]. Default value: 600 MPa.", GH_ParamAccess.item, 600);
            pManager.AddNumberParameter("Shear modulus G2", "G2", "Shear modulus [MPa]. Default value: 600 MPa.", GH_ParamAccess.item, 600);
            pManager.AddNumberParameter("Shear modulus G3", "G3", "Shear modulus [MPa]. Default value: 30 MPa.", GH_ParamAccess.item, 30);
            pManager.AddNumberParameter("Possion Ratio v12", "v12", "Possion Ratio [-] plane 12. Default value: 0.5.", GH_ParamAccess.item, 0.5);
            pManager.AddNumberParameter("Possion Ratio v13", "v13", "Possion Ratio [-] plane 13. Default value: 0.6.", GH_ParamAccess.item, 0.6);
            pManager.AddNumberParameter("Possion Ratio v23", "v23", "Possion Ratio [-] plane 23. Default value: 0.6.", GH_ParamAccess.item, 0.6);

            pManager.AddNumberParameter("Yielding stress", "fy", "Yield stress [MPa]. Default value: 30 MPa.", GH_ParamAccess.item, 30);
            pManager.AddNumberParameter("Material weight", "rho", "The weight of the material [kg/m^3]", GH_ParamAccess.item, 300);

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
            double E1 = 0.0;
            double E2 = 0.0;
            double E3 = 0.0;
            double G1 = 0.0;
            double G2 = 0.0;
            double G3 = 0.0;
            double v12 = 0.0;
            double v13 = 0.0;
            double v23 = 0.0;

            double fy = 0.0;
            double weight = 0.0;
            DA.GetData(0, ref E1);
            DA.GetData(1, ref E2);
            DA.GetData(2, ref E3);
            DA.GetData(3, ref G1);
            DA.GetData(4, ref G2);
            DA.GetData(5, ref G3);
            DA.GetData(6, ref v12);
            DA.GetData(7, ref v13);
            DA.GetData(8, ref v23);

            DA.GetData(9, ref fy);
            DA.GetData(10, ref weight);


            // Code
            MaterialOrto material = new MaterialOrto(E1, E2, E3, v12, v13, v23, G1, G2, G3, fy);
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
            get { return new Guid("0734A8EB-DAAD-46BC-88D3-949865FADA2B"); }
        }
    }
}