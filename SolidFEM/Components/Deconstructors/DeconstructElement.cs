using System;
using System.Collections.Generic;
using SolidFEM.Classes;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace SolidFEM.Deconstructors
{
    public class DeconstructElement : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public DeconstructElement()
          : base("DeconstructEl", "DeEl",
              "Deconstruct basic element class",
              "SmartMesh", "Deconstructors")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Element", "el", "solid fem element", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "m", "The elements mesh", GH_ParamAccess.item);
            pManager.AddGenericParameter("Nodes", "ns", "solid fem nodes", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Element el = new Element();
            if(!DA.GetData(0, ref el)) return;
            
            

            
            DA.SetData(0, el.ElementMesh);
            DA.SetDataList(1, el.TopologyVertices);
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
            get { return new Guid("51fe47af-a45c-400b-8511-c947dcffc77b"); }
        }
    }
}