using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;
using System.Drawing;
    

namespace SolidFEM.Preview
{
    public class MeshPreview : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MeshPreview class.
        /// </summary>
        public MeshPreview()
          : base("MeshPreview", "preview",
              "Preview the results from the analysis.",
              "SmartMesh", "Utilities")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Analysed FE-Mesh","rslt","The resulting mesh from the analysis",GH_ParamAccess.item); // 0
            pManager.AddIntegerParameter("Type", "t", "Type of results to preview. 1 = displacement; 2 = Mises stresses; 3 = utilisation.", GH_ParamAccess.item, 0); // 1
            pManager.AddNumberParameter("Displacement scaling", "scaling", "The scaling factor of the results. 1.0 gives real displacements", GH_ParamAccess.item, 1.0); // 2
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Deformed model", "FE-mesh", "The FE-mesh class with deformed geometry", GH_ParamAccess.item); // 0
            // the final component should output
            pManager.AddMeshParameter("Deformed mesh", "mesh", "The mesh representation of the deformed model", GH_ParamAccess.list);// 1
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // -- variables --
            TempFE_Mesh feMesh = new TempFE_Mesh(); ; // 0
            int type = 0; // 1
            double scaling = 0.0; // 2            

            // -- input --
            if (!DA.GetData(0, ref feMesh)) return;
            DA.GetData(1, ref type);
            DA.GetData(2, ref scaling);


            // -- solve --

            TempFE_Mesh copyMesh = feMesh.DeepCopy();
            List<Mesh> newMeshList = new List<Mesh>();
            List<Node> newNodes = new List<Node>();
            List<Element> newElements = new List<Element>();

            


            // start by creating vectors of scaling
            List<Vector3d> scaleVecs = new List<Vector3d>(); // instantiate empty list
            for (int i = 0; i < copyMesh.MeshNodes.Count; i++)
            {
                double du = copyMesh.dU[i];
                double dv = copyMesh.dV[i];
                double dw = copyMesh.dW[i];

                Vector3d vec =  Vector3d.Multiply( scaling, (new Vector3d(du, dv, dw)) ); // the scaled vector
                scaleVecs.Add(vec);


                Point3d oPt = feMesh.MeshNodes[i].Coordinate;
                Transform translate = Transform.Translation(scaleVecs[i]);
                oPt.Transform(translate);
                newNodes.Add(new Node(i, oPt)); // transform each mesh node
            }
            copyMesh.MeshNodes = newNodes;

            // update the elements positions
            foreach (Element element in feMesh.MeshElements)
            {
                Element newElement = new Element(element);
                List<int> con = element.Connectivity;                
                for (int j = 0; j < con.Count; j++)
                {
                    int nodeInd = con[j]; // the global index of the coordinate
                    newElement.Nodes[j].Coordinate = newNodes[nodeInd].Coordinate; // update the Element nodes
                }
                newElements.Add(newElement);
                newMeshList.Add(MeshFromElement(element));
            }
            copyMesh.MeshElements = newElements;
            copyMesh.MeshList = newMeshList;

            // to do: create a new FE-Mesh with the new lists
            List<Point3d> newPts = new List<Point3d>();
            foreach (Node node in newNodes)
            {
                newPts.Add(node.Coordinate);
            }

            // if von Mises stresses are to be output
            if (type == 2)
            {
                ColorMeshAfterStress(copyMesh);
            }

            // -- output --
            DA.SetData(0, copyMesh);
            DA.SetDataList(1, newMeshList);
        }

        private void ColorMeshAfterStress(TempFE_Mesh mesh)
        {
            Material material = mesh.Material;
            List<double> mises = mesh.MisesStress;
            double maxValue = material.YieldingStress;
            double minValue = 0;
            Color color = Color.White;
            double range = (maxValue - minValue) / (double)13;

            // iterate through each element
            for (int i = 0; i < mesh.MeshElements.Count; i++)
            {
                Element element = mesh.MeshElements[i];
                List<int> con = element.Connectivity; // list of global node indices

                // iterate through each node
                for (int j = 0; j < element.Nodes.Count; j++)
                {
                    int nodeInd = con[j]; // the global index of the node

                    // find the correct colour:
                    if (mises[i] < minValue + range) color = Color.Blue;
                    else if (mises[i] < minValue + 2 * range) color = Color.RoyalBlue;
                    else if (mises[i] < minValue + 3 * range) color = Color.DeepSkyBlue;
                    else if (mises[i] < minValue + 4 * range) color = Color.Cyan;
                    else if (mises[i] < minValue + 5 * range) color = Color.PaleGreen;
                    else if (mises[i] < minValue + 6 * range) color = Color.LimeGreen;
                    else if (mises[i] < minValue + 7 * range) color = Color.Lime;
                    else if (mises[i] < minValue + 8 * range) color = Color.Lime;
                    else if (mises[i] < minValue + 9 * range) color = Color.GreenYellow;
                    else if (mises[i] < minValue + 10 * range) color = Color.Yellow;
                    else if (mises[i] < minValue + 11 * range) color = Color.Orange;
                    else if (mises[i] < minValue + 12 * range) color = Color.OrangeRed;
                    else color = Color.Red;

                    // add the colour to the element
                    bool colourVertex = mesh.MeshList[i].VertexColors.SetColor(j, color);
                    if (!colourVertex)
                    {
                        throw new IndexOutOfRangeException("The vertex index for colouring is out of range");
                    }

                }
            }

            
        }
        private Mesh MeshFromElement(Element el)
        {
            List<Point3d> points = el.TopologyVertices;

            Mesh mesh = new Mesh();

            // create vertices
            mesh.Vertices.AddVertices(points);

            // create faces
            mesh.Faces.AddFace(new MeshFace(0, 1, 2, 3));
            mesh.Faces.AddFace(new MeshFace(4, 5, 6, 7));
            mesh.Faces.AddFace(new MeshFace(0, 1, 5, 4));
            mesh.Faces.AddFace(new MeshFace(1, 2, 6, 5));
            mesh.Faces.AddFace(new MeshFace(2, 3, 7, 6));
            mesh.Faces.AddFace(new MeshFace(3, 0, 4, 7));

            // clean 
            mesh.FaceNormals.ComputeFaceNormals();
            mesh.Compact();

            return mesh;
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
            get { return new Guid("9b5f0406-4f61-486f-a3f3-c614497de36d"); }
        }
    }
}