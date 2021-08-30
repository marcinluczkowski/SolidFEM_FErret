using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;

namespace SolidFEM.DeconstructClasses
{
    public class DeconstructSmartMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DeconstructMesh3d class.
        /// </summary>
        public DeconstructSmartMesh()
          : base("Deconstruct SmartMesh", "SmartMesh",
              "Deconstructing the SmartMesh Class.",
              "SmartMesh", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SmartMesh", "SM", "SmartMesh Class.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Elements", "e", "List of elements.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Nodes", "n", "List of nodes.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Geometry", "geo", "Geometry information.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh", "mesh", "Mesh.", GH_ParamAccess.item);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            SmartMesh mesh = new SmartMesh();
            DA.GetData(0, ref mesh);

            // Output
            DA.SetDataList(0, mesh.Elements);
            DA.SetDataList(1, mesh.Nodes);
            DA.SetData(2, mesh.Geometry);
            DA.SetData(3, mesh.Mesh);
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
                return null; // return Properties.Resources.Icon_DecSmartMesh;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("97c30c27-48c9-41ac-b09d-d02f80e806f6"); }
        }
    }
}