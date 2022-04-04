using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Linq;
using ClosedXML.Excel;

// Csparse
using LA = MathNet.Numerics.LinearAlgebra;
using CSparse;
using CSD = CSparse.Double;
using CSparse.Double.Factorization;
using CSparse.Storage;
using Grasshopper;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Point = Rhino.Geometry.Point;

namespace SolidFEM
{
    public class ConvertMeshTo20Node : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public ConvertMeshTo20Node()
          : base("ConvertMeshTo20Node", "",
              "Converts a 8-node mesh to 20-node mesh",
              "Category", "Subcategory")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "The 8-node mesh", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "The 20-node mesh", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {

            List<Mesh> meshList = new List<Mesh>();
            DA.GetDataList(0, meshList);

            List<Mesh> newMesh = new List<Mesh>();
            
            foreach (Mesh mesh in meshList)
            {
                if (mesh.Vertices.Count == 8)
                {
                    //Mesh nM = GrahamScan.DoGrahamScan(mesh);
                    Mesh nM = mesh;
                    //if (nM.IsValid)
                    //{

                        var vertices_array = nM.Vertices.ToPoint3dArray(); // an array of the vertices of each mesh element

                        List<Point3d> vertices = vertices_array.ToList();


                        nM.Vertices.Add((vertices[0] + vertices[1]) / 2);
                        nM.Vertices.Add((vertices[1] + vertices[2]) / 2);
                        nM.Vertices.Add((vertices[2] + vertices[3]) / 2);
                        nM.Vertices.Add((vertices[3] + vertices[0]) / 2);
                        nM.Vertices.Add((vertices[4] + vertices[5]) / 2);
                        nM.Vertices.Add((vertices[5] + vertices[6]) / 2);
                        nM.Vertices.Add((vertices[6] + vertices[7]) / 2);
                        nM.Vertices.Add((vertices[7] + vertices[4]) / 2);
                        nM.Vertices.Add((vertices[0] + vertices[4]) / 2);
                        nM.Vertices.Add((vertices[1] + vertices[5]) / 2);
                        nM.Vertices.Add((vertices[2] + vertices[6]) / 2);
                        nM.Vertices.Add((vertices[3] + vertices[7]) / 2);


                        newMesh.Add(nM);
                    //}
                    /*
                    else
                    {

                        var vertices_array = mesh.Vertices.ToPoint3dArray(); // an array of the vertices of each mesh element

                        List<Point3d> vertices = vertices_array.ToList();


                        mesh.Vertices.Add((vertices[0] + vertices[1]) / 2);
                        mesh.Vertices.Add((vertices[1] + vertices[2]) / 2);
                        mesh.Vertices.Add((vertices[2] + vertices[3]) / 2);
                        mesh.Vertices.Add((vertices[3] + vertices[0]) / 2);
                        mesh.Vertices.Add((vertices[4] + vertices[5]) / 2);
                        mesh.Vertices.Add((vertices[5] + vertices[6]) / 2);
                        mesh.Vertices.Add((vertices[6] + vertices[7]) / 2);
                        mesh.Vertices.Add((vertices[7] + vertices[4]) / 2);
                        mesh.Vertices.Add((vertices[0] + vertices[4]) / 2);
                        mesh.Vertices.Add((vertices[1] + vertices[5]) / 2);
                        mesh.Vertices.Add((vertices[2] + vertices[6]) / 2);
                        mesh.Vertices.Add((vertices[3] + vertices[7]) / 2);


                        newMesh.Add(mesh);
                    }
                    */
                }
            }


            DA.SetDataList(0, newMesh);

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
            get { return new Guid("F4130C34-5BE4-4A78-9761-075B8E4C99A2"); }
        }
    }
}