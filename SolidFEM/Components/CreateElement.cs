using System;
using System.Collections.Generic;
using SolidFEM.Classes;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace SolidFEM.Components
{
    public class CreateElement : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public CreateElement()
          : base("CreateEkement", "cel",
              "Create single FEelement",
              "FEM", "Element")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddBrepParameter("singleBrep","brp","box brep to make single FEelement",GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Informations","","",GH_ParamAccess.list);
            pManager.AddGenericParameter("Element", "el", "solid fem element", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Brep b = new Brep();
            DA.GetData(0, ref b);
            var info = new List<string>();

            List<Node> nodes = new List<Node>();
            int i = 0;
            foreach (var v in b.Vertices)
            {
                Point3d pt = v.Location;
                Node n = new Node();
                n.ID = i++;
                n.Name = "3D node";
                n.Point = pt;
                nodes.Add(n);
            }

            Element el = new Element();
            el.ID = 0;
            el.name = "8node3D";
            el.nodes = nodes;


            info.Add("change brep box into fel"); ;
            DA.SetDataList(0, info);
            DA.SetData(1, el);
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
            get { return new Guid("5438539b-fb3c-4dc1-81d0-bb59141146a9"); }
        }
    }
}