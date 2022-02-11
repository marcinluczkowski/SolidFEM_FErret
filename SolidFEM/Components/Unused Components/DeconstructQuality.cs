using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;

namespace SolidFEM.DeconstructClasses
{
    public class DeconstructQuality : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DeconstructQuality class.
        /// </summary>
        public DeconstructQuality()
          : base("Deconstruct Quality", "quality",
              "Deconstructing the Quality Class.",
              "SmartMesh", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Quality", "q", "Quality Class.", GH_ParamAccess.item);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Element", "e", "Element.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Aspect Ratio", "AR", "Ratio between shortest and longest element edge.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Skewness", "SK", "Angle ratio of the element.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Jacobian Ratio", "JR", "Jacobian ratio of the element.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //input
            Quality q = new Quality();
            DA.GetData(0, ref q);

            //output
            DA.SetData(0, q.element);
            DA.SetData(1, q.AspectRatio);
            DA.SetData(2, q.Skewness);
            DA.SetData(3, q.JacobianRatio);
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
                return null;// return Properties.Resources.Icon_DecQuality;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("ef2f555e-2d22-4c47-8f10-eb6f7ce8ff52"); }
        }
    }
}