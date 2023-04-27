using System;
using System.Collections.Generic;
using Accord.Collections;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using SolidFEM.Classes;

namespace SolidFEM.Components
{
    public class FEMPart : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Part class.
        /// </summary>
        public FEMPart()
          : base("FEM Part", "Part",
              "Create part for the FEM Solver.",
              "SolidFEM", "FEM-Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh List", "mList", "List of mesh representing solid elements.", GH_ParamAccess.list); // 0
            pManager.AddGenericParameter("Material", "mat", "List of materials", GH_ParamAccess.list); // 1
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Meshes", "Meshes", "List of meshes representing the elements", GH_ParamAccess.list);
            pManager.AddGenericParameter("Materials", "Materials", "List of the material for each element", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Mesh> meshList = new List<Mesh>();
            List<Material> material = new List<Material>();
            DA.GetDataList(0, meshList);
            DA.GetDataList(1, material);

            List<Material> materialList = new List<Material>();

            for (int i = 0; i < material.Count;i++)
            {
                for (int j = 0; j < meshList.Count; j++)
                { 
                    materialList.Add(material[i]); 
                }
            }

            DA.SetDataList(0, meshList);
            DA.SetDataList(1, materialList);
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
            get { return new Guid("C6D18CA2-1031-4409-A6BA-D2EAF0DAB2F9"); }
        }
    }
}