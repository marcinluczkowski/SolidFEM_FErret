using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;

namespace SolidFEM.DeconstructClasses
{
    public class DeconstructFEMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DeconstructGeometry class.
        /// </summary>
        public DeconstructFEMesh()
          : base("Deconstruct FEMesh", "decMesh",
              "Deconstructing the FEMesh class.",
              "SmartMesh", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("FEMesh", "femesh", "Deconstructs the FEMesh", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Nodes", "f", "List of faces.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Elements", "e", "List of edges.", GH_ParamAccess.list);
        }


        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            TempFE_Mesh feMesh = new TempFE_Mesh();
            if(!DA.GetData(0, ref feMesh)) return;

            // Output
            DA.SetDataList(0, feMesh.MeshNodes);
            DA.SetDataList(1, feMesh.MeshElements);
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
            get { return new Guid("d33184d3-646b-44a6-b78a-02841973f386"); }
        }
    }
}