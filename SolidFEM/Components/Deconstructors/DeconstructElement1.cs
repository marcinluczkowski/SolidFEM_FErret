using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;

namespace SolidFEM.DeconstructClasses
{
    public class DeconstructElement : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        /// <summary>
        /// Initializes a new instance of the Deconstruct_Node class.
        /// </summary>
        public DeconstructElement()
        : base("Deconstruct Element", "element",
              "Deconstructing the Element Class.",
              "SmartMesh", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Element", "e", "Element Class.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Nodes", "n", "Nodes of element.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Connectivity", "con", "Connectivity of local to global nodes.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Type", "t", "Element type.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Id", "id", "Element Id.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Mesh", "mesh", "Element mesh.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Determinant", "det", "Jacobian ratio of mesh", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            Element e = new Element();
            DA.GetData(0, ref e);

            // Output
            DA.SetDataList(0, e.Nodes);
            DA.SetDataList(1, e.Connectivity);
            DA.SetData(2, e.Type);
            DA.SetData(3, e.Id);
            DA.SetData(4, e.ElementMesh);

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
                return null; // Properties.Resources.Icon_DecElement;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("df3f3c5f-c32c-44ce-83f1-831a94edd1d8"); }
        }
    }
}