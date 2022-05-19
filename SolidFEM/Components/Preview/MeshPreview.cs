using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using SolidFEM.Classes;
using System.Drawing;
using System.Linq;

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
            pManager.AddIntegerParameter("Type", "t", "Type of results to preview. 1 = displacement; 2 = Mises stresses; 3 = utilisation; 4 = sigma_xx; 5 = sigma_yy; 6 = sigma_zz.", GH_ParamAccess.item, 0); // 1
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

                //copyMesh.MeshNodes[i].Coordinate = Point3d.Add(copyMesh.MeshNodes[i].Coordinate, vec); // Move the Nodal coordinates
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
                    newElement.Nodes[j] = copyMesh.MeshNodes[nodeInd]; // update the Element nodes
                }
                newElements.Add(newElement);
                newMeshList.Add(MeshFromElement(newElement));
            }
            copyMesh.MeshElements = newElements;
            copyMesh.MeshList = newMeshList;

            // to do: create a new FE-Mesh with the new lists
            List<Point3d> newPts = new List<Point3d>();
            foreach (Node node in newNodes)
            {
                newPts.Add(node.Coordinate);
            }
            if (type == 1)
            {
                // preview displacement as mesh colours
                ColorMeshAfterDisplacements(copyMesh);
            }

            // if von Mises stresses are to be output
            else if (type == 2)
            {
                MeshColorMises(copyMesh);
            }
            else if (type == 3)
            {
               MeshColorUtilisation(copyMesh);
            }
            else if (type == 4)
            {
                MeshColorStresses(copyMesh, copyMesh.Sigma_xx);
            }
            else if (type == 5)
            {
                MeshColorStresses(copyMesh, copyMesh.Sigma_yy);
            }
            else if (type == 6)
            {
                MeshColorStresses(copyMesh, copyMesh.Sigma_zz);
            }
            else if (type == 7)
            {
                MeshColorStresses(copyMesh, copyMesh.Sigma_xy);
            }
            else if (type == 8)
            {
                MeshColorStresses(copyMesh, copyMesh.Sigma_xz);
            }
            else if (type == 9)
            {
                MeshColorStresses(copyMesh, copyMesh.Sigma_yz);
            }
            // -- output --
            DA.SetData(0, copyMesh);
            DA.SetDataList(1, newMeshList);
        }
        private void ColorMeshAfterDisplacements(TempFE_Mesh mesh)
        {
            List<Color> colourPalette = DisplacementColors();

            // get max displacement
            List<double> displacements = new List<double>();
            for (int i = 0; i < mesh.dU.Count; i++)
            {
                double disp = Math.Sqrt( Math.Pow(mesh.dU[i], 2) + Math.Pow(mesh.dV[i], 2) + Math.Pow(mesh.dW[i], 2)); // calculate absolute displacement
                displacements.Add(disp);
            }

            double maxDisp = displacements.Max() + 0.001;
            double min = 0;
            double stepSize = (maxDisp - 0) / colourPalette.Count;

            // iterate through all elements
            for (int i = 0; i < mesh.MeshElements.Count; i++)
            {
                Element el = mesh.MeshElements[i];
                List<int> con = el.Connectivity;

                for (int j = 0; j < con.Count; j++)
                {
                    int globalInd = con[j]; // index of the global node to take the displacement from
                    int colInd = (int)Math.Truncate(displacements[globalInd] / stepSize); // the index to use for vertex colours
                    mesh.MeshList[i].VertexColors.SetColor(j, colourPalette[colInd]); // colour the mesh vertex
                }
            }
        }
        private void ColorMeshAfterMisesStress2(TempFE_Mesh mesh)
        {
            List<Color> colourPalette = MisesColors();
            List<double> misesValues = mesh.MisesStress; // the element mises stresses


            // -- calculate the average mises stresses of each node --

            double[] nodalStresses = new double[mesh.MeshNodes.Count]; // new list double
            double[] nodalCount = new double[mesh.MeshNodes.Count]; // number of element connected to each node
            // iterate through each element
            for (int i = 0; i < mesh.MeshElements.Count; i++)
            {
                Element el = mesh.MeshElements[i];
                List<int> con = el.Connectivity;

                for (int j = 0; j < con.Count; j++)
                {
                    double elStress = misesValues[i];
                    nodalStresses[con[j]] += elStress; // add the element stress to the node
                    nodalCount[con[j]] += 1; // add one node to the nodal count
                }
            }
            // calculate the average stresses
            List<double> avgMises = new List<double>();
            for (int i = 0; i < nodalStresses.Length; i++)
            {
                avgMises.Add( nodalStresses[i] / nodalCount[i] ); // the average mises stresses of each node
            }

            // -- colour the mesh vertices --

            double maxMises = avgMises.Max() + 0.001;
            double stepSize = (maxMises - 0) / colourPalette.Count;
            for (int i = 0; i < mesh.MeshElements.Count; i++)
            {
                Element el = mesh.MeshElements[i];
                List<int> con = el.Connectivity;

                for (int j = 0; j < con.Count; j++)
                {
                    
                    int globalInd = con[j]; // global node index
                    int colourInd = (int)Math.Truncate(avgMises[globalInd] / stepSize); // the index to use for vertex colours
                    mesh.MeshList[i].VertexColors.SetColor(j, colourPalette[colourInd]); // set the mesh vertex colour
                }
            }

        }
        private void MeshColorUtilisation(TempFE_Mesh mesh)
        {
            List<double> nodalMises = mesh.NodelMisesStresses;
            double yieldStress = mesh.Material.YieldingStress;
            List<double> utilisation = new List<double>();
            // divide all stresses with yield stress
            for (int i = 0; i < nodalMises.Count; i++)
            {
                utilisation.Add(nodalMises[i] / yieldStress);
            }
            // colours
            List<Color> positiveCols = PositiveStressColors();
            List<Color> negativeCols = NegativeStressColors();

            double stepSize = 1.001 / positiveCols.Count;

            for (int i = 0; i < mesh.MeshElements.Count; i++)
            {
                Element el = mesh.MeshElements[i];
                List<int> con = el.Connectivity;

                for (int j = 0; j < con.Count; j++)
                {
                    int globalInd = con[j];
                    double util = utilisation[globalInd];

                    if (Math.Abs(util) > 1.0)
                    {
                        // the allowed stress is too high. For now, all nodes with utilisation above 1.0 is black regardless of positive or negative value
                        mesh.MeshList[i].VertexColors.SetColor(j, Color.Black);
                        continue;
                    }

                    if (util > 0.0)
                    {
                        int colInd = (int)Math.Truncate(util / stepSize);
                        mesh.MeshList[i].VertexColors.SetColor(j, positiveCols[colInd]);
                    }
                    else
                    {
                        int colInd = (int)Math.Truncate(Math.Abs(util) / stepSize);
                        mesh.MeshList[i].VertexColors.SetColor(j, negativeCols[colInd]);
                    }


                }
            }
            
        }
        private void MeshColorMises(TempFE_Mesh mesh)
        {
            List<Color> cPalette = MisesColors();
            List<double> nodalMises = mesh.NodelMisesStresses;

            double maxStress = nodalMises.Max() + 0.001;
            double stepSize = (maxStress - 0.0 ) / cPalette.Count;

            for (int i = 0; i < mesh.MeshElements.Count; i++)
            {
                Element el = mesh.MeshElements[i];
                List<int> con = el.Connectivity;

                for (int j = 0; j < con.Count; j++)
                {
                    int globalInd = con[j];
                    int colInd = (int)Math.Truncate( nodalMises[globalInd] / stepSize );
                    mesh.MeshList[i].VertexColors.SetColor(j, cPalette[colInd]);
                }
            }
        }
        
        private void MeshColorStresses(TempFE_Mesh mesh, List<double> stresses)
        {
            List<Color> posColors = PositiveStressColors();
            List<Color> negColors = NegativeStressColors();

            double maxStress = stresses.Max() + 0.001;
            double minStress = stresses.Min() - 0.001;

            double posStep = (maxStress - 0.0) / posColors.Count;
            double negStep = (minStress - 0.0) / negColors.Count;

            for (int i = 0; i < mesh.MeshElements.Count; i++)
            {
                Element el = mesh.MeshElements[i];
                List<int> con = el.Connectivity;

                for (int j = 0; j < con.Count; j++)
                {
                    int globalInd = con[j];
                    double stress = stresses[globalInd];

                    if (stress > 0)
                    {
                        int colInd = (int)Math.Truncate(stress / posStep);
                        mesh.MeshList[i].VertexColors.SetColor(j, posColors[colInd]);
                    }
                    else
                    {
                        int colInd = (int)Math.Truncate(stress / negStep);
                        mesh.MeshList[i].VertexColors.SetColor(j, negColors[colInd]);
                    }

                }
            }
        }
        /// <summary>
        /// Create a Rhino mesh from an element
        /// </summary>
        /// <param name="el">Rhino element to create mesh from. Currently HEX8 and TET4 are implemented</param>
        /// <returns></returns>
        private Mesh MeshFromElement(Element el)
        {
            List<Point3d> points = el.TopologyVertices;
            Mesh oldM = el.ElementMesh;
            
            Mesh mesh = new Mesh();

            // create vertices
            mesh.Vertices.AddVertices(points);

            if (points.Count == 8 || points.Count == 20) // HEX8 element
            {
                // create faces for HEX8
                mesh.Faces.AddFace(new MeshFace(0, 1, 2, 3));
                mesh.Faces.AddFace(new MeshFace(4, 5, 6, 7));
                mesh.Faces.AddFace(new MeshFace(0, 1, 5, 4));
                mesh.Faces.AddFace(new MeshFace(1, 2, 6, 5));
                mesh.Faces.AddFace(new MeshFace(2, 3, 7, 6));
                mesh.Faces.AddFace(new MeshFace(3, 0, 4, 7));
            }
            else if (points.Count == 4)
            {
                // create faces for TET4
                foreach (MeshFace face in oldM.Faces)
                {
                    mesh.Faces.AddFace(new MeshFace(face.A, face.B, face.C));
                }
                /*
                mesh.Faces.AddFace(new MeshFace(0, 1, 3));
                mesh.Faces.AddFace(new MeshFace(1, 2, 3));
                mesh.Faces.AddFace(new MeshFace(0, 2, 3));
                mesh.Faces.AddFace(new MeshFace(0, 1, 2));
                */
            }
            


            // clean 
            mesh.FaceNormals.ComputeFaceNormals();
            //mesh.Compact();

            return mesh;
        }

        private List<Color> DisplacementColors()
        {
            List<Color> colours = new List<Color>(); //

            // add colour gradient from "https://colordesigner.io/gradient-generator"
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#fef1fc")); // 0
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#f7e0f4")); // 1
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#efd0ec")); // 2
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#e7c0e5")); // 3
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#deb0de")); // 4
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#d5a0d7")); // 5
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#cb90d0")); // 6
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#c180ca")); // 7
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#b671c4")); // 8
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#ab61be")); // 9
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#a052b8")); // 10
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#9442b2")); // 11
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#8731ad")); // 12
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#791ea7")); // 13
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#6b00a2")); // 14

            return colours;
        }

        private List<Color> PositiveStressColors()
        {
            List<Color> colours = new List<Color>();

            List<string> colNames = new List<string>() { "#ffffff", "#ecdeff","#d7bdff",  "#c19cff","#a77cff", "#8a5bff", "#6438ff", "#1f02fd" };

            foreach (string colName in colNames)
            {
                colours.Add(ColorTranslator.FromHtml(colName));
            }
            colours.Add(ColorTranslator.FromHtml("#1f02fd")); // 7

            return colours;
        }
        private List<Color> NegativeStressColors()
        {
            List<Color> colours = new List<Color>();

            List<string> colNames = new List<string>() { "#ffffff", "#ffe4da", "#ffc8b6", "#ffac93", "#ff8f71", "#ff704f", "#ff4b2d", "#fd0202" };

            foreach (string colName in colNames)
            {
                colours.Add(ColorTranslator.FromHtml(colName));
            }
            colours.Add(ColorTranslator.FromHtml("#1f02fd")); // 7

            return colours;
        }

        public static List<Color> GetGradientColors(Color start, Color end, int steps, int firstStep, int lastStep)
        {
            var colorList = new List<Color>();
            if (steps <= 0 || firstStep < 0 || lastStep > steps - 1)
                return colorList;

            double aStep = (end.A - start.A) / steps;
            double rStep = (end.R - start.R) / steps;
            double gStep = (end.G - start.G) / steps;
            double bStep = (end.B - start.B) / steps;

            for (int i = firstStep; i < lastStep; i++)
            {
                var a = start.A + (int)(aStep * i);
                var r = start.R + (int)(rStep * i);
                var g = start.G + (int)(gStep * i);
                var b = start.B + (int)(bStep * i);
                colorList.Add(Color.FromArgb(a, r, g, b));
            }

            return colorList;
        }

        private List<Color> MisesColors()
        {
            List<Color> colours = new List<Color>();


            // add colour gradient from "https://colordesigner.io/gradient-generator"
          

            colours.Add(System.Drawing.ColorTranslator.FromHtml("#fef1fc")); // 0
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#f7e0f4")); // 1
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#efd0ec")); // 2
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#e7c0e5")); // 3
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#deb0de")); // 4
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#d5a0d7")); // 5
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#cb90d0")); // 6
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#c180ca")); // 7
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#b671c4")); // 8
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#ab61be")); // 9
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#a052b8")); // 10
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#9442b2")); // 11
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#8731ad")); // 12
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#791ea7")); // 13
            colours.Add(System.Drawing.ColorTranslator.FromHtml("#6b00a2")); // 14



            return colours;


        }
        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.FEPreview;
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