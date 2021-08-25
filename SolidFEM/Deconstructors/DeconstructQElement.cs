using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints.DeconstructClasses
{
    public class DeconstructQElement : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DeconstructQElement class.
        /// </summary>
        public DeconstructQElement()
          : base("Deconstruct qElement", "qElement",
              "Deconstructing the qElement Class.",
              "SmartMesh", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("qElement", "qEl", "qElement Class.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("qEdges", "qE", "Element edges.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Angles", "qE", "Element angles.", GH_ParamAccess.list);
            pManager.AddGenericParameter("IsQuad", "qE", "True if element is a quad.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Distortion Metric", "dM", "Distortion metric of element.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Contour", "c", "Contour of element.", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            qElement element = new qElement();
            DA.GetData(0, ref element);

            DA.SetDataList(0, element.EdgeList);
            DA.SetDataList(1, element.AngleList);
            DA.SetData(2, element.IsQuad);
            DA.SetData(3, element.DistortionMetric);
            DA.SetDataList(4, element.Contour);
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
                return Properties.Resources.Icon_DecQElement;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("6c7a5125-b4e9-45f1-992b-1bcfe1a9274e"); }
        }
    }
}