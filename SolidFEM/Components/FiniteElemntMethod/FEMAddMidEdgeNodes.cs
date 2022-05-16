using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Matrix = Accord.Math.Matrix;

namespace SolidFEM.Components.FiniteElemntMethod
{
    public class MyComponent1 : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public MyComponent1()
          : base("AddMidEdgeNodes", "TET10",
              "Takes a TET4 mesh and returns a TET10 mesh with nodes on edge midpoints",
              "SmartMesh", "FEM-Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "mesh", "TET4 mesh", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "mesh", "TET10 mesh", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // Input
            Mesh mesh = new Mesh();

            DA.GetData(0, ref mesh);


            // Main work

            var meshPts = mesh.TopologyVertices;
            Mesh newMesh = new Mesh();

            //newMesh.Vertices.AddVertices(meshPts);

            newMesh.Vertices.Add(meshPts[0]);
            newMesh.Vertices.Add(meshPts[1]);
            newMesh.Vertices.Add(meshPts[2]);
            newMesh.Vertices.Add(meshPts[3]);

            newMesh.Faces.AddFace(0, 2, 1);
            newMesh.Faces.AddFace(0, 1, 3);
            newMesh.Faces.AddFace(1, 2, 3);
            newMesh.Faces.AddFace(2, 0, 3);

            newMesh.Compact();
            newMesh.FaceNormals.ComputeFaceNormals();

            var newMeshPts = newMesh.TopologyVertices;


            int[,] midNodeIndices = new int[,]  //Array for getting the correct ordering of midside nodes
            {
                {0,1},  // node 5 is between node 1 and 2
                {1,2},  // node 6 is between node 2 and 3
                {0,2},  // node 7 and so on ...
                {0,3},  // node 8...
                {1,3},  // node 9...
                {2,3}   // node 10... 
            };

            for (int i = 0; i < Matrix.Rows(midNodeIndices); i++)
            {
                Point3d pt1 = newMeshPts[midNodeIndices[i, 0]];
                Point3d pt2 = newMeshPts[midNodeIndices[i, 1]];
                int rp = 6; // rounding precision
                double x1 = Math.Round(pt1.X, rp);
                double x2 = Math.Round(pt2.X, rp);
                double y1 = Math.Round(pt1.Y, rp);
                double y2 = Math.Round(pt2.Y, rp);
                double z1 = Math.Round(pt1.Z, rp);
                double z2 = Math.Round(pt2.Z, rp);

                Point3d midPoint = new Point3d((x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2);
                newMesh.Vertices.Add(midPoint);
            }

            DA.SetData(0, newMesh);

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
            get { return new Guid("6933132B-1D8B-413C-98A4-929081013EB2"); }
        }
    }
}