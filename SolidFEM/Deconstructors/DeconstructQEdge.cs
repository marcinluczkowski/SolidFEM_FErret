using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using MeshPoints.Classes;

namespace MeshPoints.DeconstructClasses
{
    public class DeconstructQEdge : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DeconstructQEdge class.
        /// </summary>
        public DeconstructQEdge()
          : base("Deconstruct qEdge", "qEdge",
              "Deconstructing the qEdge Class.",
              "SmartMesh", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("qEdge", "qEg", "qEdge Class.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Start node", "sN", "Start node of edge.", GH_ParamAccess.item);
            pManager.AddGenericParameter("End node", "eN", "End node of edge.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Element 1", "e1", "Adjacent element.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Element 2", "e2", "Adjacent element.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Left front neighbor ", "lF", "Adjacent front egde to left.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Right front neighbor ", "rF", "Adjacent front egde to right.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Line", "l", "Edge line.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Level", "l", "Level.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            qEdge edge = new qEdge();
            DA.GetData(0, ref edge);
            DA.SetData(0, edge.StartNode);
            DA.SetData(1, edge.EndNode);
            DA.SetData(2, edge.Element1);
            DA.SetData(3, edge.Element2);
            DA.SetData(4, edge.LeftFrontNeighbor);
            DA.SetData(5, edge.RightFrontNeighbor);
            DA.SetData(6, edge.EdgeLine);
            DA.SetData(7, edge.Level);
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
                return Properties.Resources.Icon_DecQEdge;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("ece619a7-7e52-461c-ae93-841ce95f6ac3"); }
        }
    }
}