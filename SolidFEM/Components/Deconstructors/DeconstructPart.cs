using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
//using Rhino.Render.ChangeQueue;
using SolidFEM.Classes;

namespace SolidFEM.Components.Deconstructors
{
    public class DeconstructPart : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DeconstructPart class.
        /// </summary>
        public DeconstructPart()
          : base("Deconstruct Part", "decPart",
              "Deconstructing the Part class.",
              "SolidFEM", "Deconstruct")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Part", "part", "part to deconstruct", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "mesh", "Mesh from part", GH_ParamAccess.list);
            pManager.AddGenericParameter("Material", "mat", "Material from part", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<PartClass> part = new List<PartClass>();
            DA.GetDataList(0, part);

            List<Mesh> mesh = new List<Mesh>();
            List<Material> material = new List<Material>();

            for(int i=0;i<part.Count;i++) 
            {
                for (int j = 0; j < part[i].mesh.Count; j++)
                { 
                    mesh.Add(part[i].mesh[j]);
                    material.Add(part[i].material[j]);
                }

                //material.Add(part[i].material);
            }

            DA.SetDataList(0, mesh);
            DA.SetDataList(1, material);
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
            get { return new Guid("6E1381A4-85B5-457E-9871-5874658127CD"); }
        }
    }
}