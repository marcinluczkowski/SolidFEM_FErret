using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;

namespace SolidFEM.DeconstructClasses
{
    public class DeconstructQNode : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DeconstructQNode class.
        /// </summary>
        public DeconstructQNode()
          : base("Deconstruct qNode", "qNode",
              "Deconstructing qNode Class.",
              "SmartMesh", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("qNode", "qN", "qNode Class.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Point", "pt", "Coordinate of node.", GH_ParamAccess.item);
            pManager.AddGenericParameter("IsBoundaryNode", "Boundary", "If true, node is boundary node.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            qNode node = new qNode();
            DA.GetData(0, ref node);
            DA.SetData(0, node.Coordinate);
            DA.SetData(1, node.BoundaryNode);
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
                return null;// return Properties.Resources.Icon_DecQNode;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("ca709769-2507-4b33-880c-faf88fc0d3ce"); }
        }
    }
}