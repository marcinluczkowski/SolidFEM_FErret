using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;

namespace SolidFEM.DeconstructClasses
{
    public class DeconstructGeometry : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DeconstructGeometry class.
        /// </summary>
        public DeconstructGeometry()
          : base("Deconstruct Geometry", "geo",
              "Deconstructing the Geometry Class.",
              "SmartMesh", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Geometry", "geo", "Geometry Class", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Faces", "f", "List of faces.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Edges", "e", "List of edges.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Vertices", "v", "List of vertices.", GH_ParamAccess.list);
        }


        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            Geometry geometry = new Geometry();
            DA.GetData(0, ref geometry);

            // Output
            DA.SetDataList(0, geometry.Faces);
            DA.SetDataList(1, geometry.Edges);
            DA.SetDataList(2, geometry.Vertices);
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
                return null;// return Properties.Resources.Icon_DecGeometry;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("bc020810-0e46-4545-a123-a06ed738cb51"); }
        }
    }
}