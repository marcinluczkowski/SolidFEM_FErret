using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;

namespace SolidFEM.DeconstructClasses
{
    public class Deconstruct_Node : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Deconstruct_Node class.
        /// </summary>
        public Deconstruct_Node()
          : base("Deconstruct Node", "node",
              "Deconstructing the Node Class",
              "SmartMesh", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Node", "n", "Node Class", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Global Id", "id", "Global Id", GH_ParamAccess.item);
            pManager.AddGenericParameter("Coordinate", "coord", "Coordinate of node.", GH_ParamAccess.item);
            pManager.AddGenericParameter("BC u-dir", "bcu", "Boundary condtion in u-direction", GH_ParamAccess.item);
            pManager.AddGenericParameter("BC v-dir", "bcv", "Boundary condtion in v-direction", GH_ParamAccess.item);
            pManager.AddGenericParameter("BC w-dir", "bcw", "Boundary condtion in w-direction", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //input
            Node n = new Node();
            DA.GetData(0, ref n);

            //output
            DA.SetData(0, n.GlobalId);
            DA.SetData(1, n.Coordinate);
            DA.SetData(2, n.BC_U);
            DA.SetData(3, n.BC_V);
            DA.SetData(4, n.BC_W);
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
                return null;// Properties.Resources.Icon_DecNode;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("3dd1c36b-9f95-4d31-aa83-ebea7b7d0ea9"); }
        }
    }
}