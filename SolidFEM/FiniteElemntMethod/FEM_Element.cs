using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;
using System.Linq;

namespace SolidFEM.FiniteElemntMethod
{
    public class FEM_Element : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEM_Element class.
        /// </summary>
        public FEM_Element()
          : base("FEM_Element", "Nickname",
              "Creates solid elements from breps",
              "SmartMesh", "FEM")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Input", "in", "", GH_ParamAccess.list); // 0
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Output", "out", "", GH_ParamAccess.item); // 1
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // --- variables ---
            List<Brep> breps = new List<Brep>(); // 0

            // --- declare input ---
            if (!DA.GetDataList(0, breps)) return;


            // --- method ---

            // create elements from breps
            List<Element> elements = new List<Element>();
            List<Node> all_nodes = new List<Node>();
            
            for(int i = 0; i < breps.Count; i++)
            {
                Point3d[] vertices = breps[i].DuplicateVertices();
                List<Node> el_nodes = new List<Node>();
                
                for(int j = 0; i < vertices.Length; i++)
                {
                    Node n = new Node(vertices[j], j);
                    el_nodes.Add(n);
                    all_nodes.Add(n);
                }
                elements.Add(new Element(el_nodes, i));


            }

            // create SmartMesh from elements
            SmartMesh mesh = new SmartMesh(breps);
            //  SmartMesh mesh = new SmartMesh(all_nodes, elements, "Solid");


            /// <summary 
            /// Want to create a component which creates a solid mesh from a series of input breps?
            /// Test with a straigh line of box breps (HEX-elements)
            /// </summary>


            // --- output
            DA.SetData(0, mesh);

        }

        /// <summary>
        /// Additional functions for the mesh
        /// </summary>
        private Brep CreateGeometry(List<Brep> breps)
        {
            Brep brep = breps[0];
            //breps[0].





            return brep;
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
            get { return new Guid("8530dba7-febd-4c12-bb61-9db15696c1d4"); }
        }
    }
}