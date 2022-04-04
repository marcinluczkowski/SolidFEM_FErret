using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using LA = MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics;
using SolidFEM.Classes;
using CSD = CSparse.Double;
using System.Linq;


namespace SolidFEM.FiniteElementMethod
{
    public class FEMLoad_MESH : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the FEMLoad class.
        /// </summary>
        public FEMLoad_MESH()
          : base("FEM Load", "Load",
              "Create load for the FEM Solver.",
              "SmartMesh", "FEM-Mesh")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh List", "mList", "List of mesh representing solid elements.", GH_ParamAccess.list); // 0
            pManager.AddIntegerParameter("Type", "type", "Load type: Point load = 1, Surface load = 2.", GH_ParamAccess.item); // 1
            pManager.AddGenericParameter("Position", "pos", "If type = 1: List of coordinates for point loads.", GH_ParamAccess.list); // 2
            pManager.AddSurfaceParameter("Surface", "srf", "If type = 2: surface to apply load", GH_ParamAccess.list); // 3
            pManager.AddGenericParameter("Vector", "vec", "List of vectors of the loads. If surface load, only one vector.", GH_ParamAccess.list); // 4

            pManager[0].Optional = true;
            pManager[1].Optional = true;
            pManager[2].Optional = true;
            pManager[3].Optional = true;


        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Load", "load", "List of residual forces (R).", GH_ParamAccess.list);
            pManager.AddGenericParameter("Points", "pts", "List of points subjected to load.", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Resultant load and area", "R", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {

            // Input

            //SmartMesh smartMesh = new SmartMesh();
            List<Mesh> meshList = new List<Mesh>(); // 0
            int loadType = 0; // 1
            List<Vector3d> loadVectors = new List<Vector3d>();
            List<Point3d> loadPosition = new List<Point3d>();
            List<Surface> surfaces = new List<Surface>();


            DA.GetDataList(0, meshList);
            DA.GetData(1, ref loadType);
            DA.GetDataList(2, loadPosition);
            DA.GetDataList(3, surfaces);
            DA.GetDataList(4, loadVectors);


            // Code
            // clean the mesh and sort nodes
            var newMeshList = new List<Mesh>();
            int c = 0; // delete after testing
            foreach (Mesh mesh in meshList)
            {
                if (mesh.Vertices.Count == 8)
                {
                    Mesh nM = GrahamScan.DoGrahamScan(mesh);

                    if (nM.IsValid)
                    {
                        
                        newMeshList.Add(nM);
                    }
                    else newMeshList.Add(mesh);
                    c++;
                }
                else if (elementType == "Tet10")
                {
                    Mesh nM = FEM_Utility.AddMidEdgeNodes(mesh);
                    newMeshList.Add(nM);
                }
                else
                {
                    newMeshList.Add(mesh);
                }

            }
                    else
                    {
                        newMeshList.Add(mesh);
                    }
                    c++;

                }
                else if (mesh.Vertices.Count == 20)
                {
                    newMeshList.Add(mesh);
                }
                else
                {
                    newMeshList.Add(mesh);
                }
            }


            double totalArea = 0;

            List<Point3d> pointsWithLoad = new List<Point3d>();

            List<Point3d> nodePts = FEM_Utility.GetMeshNodes(newMeshList);

            int numNodes = nodePts.Count;

            double[] residualForces = new double[numNodes * 3]; // number of nodes times the three global directions

            // Assign external load
            if (loadType == 1) // point load
            {
                // control the length of load vector- and position list
                int loadLength = loadVectors.Count;
                int posLength = loadPosition.Count;

                if (loadLength != posLength) 
                {
                    // ensure that longest list principle applies
                    if (posLength > loadLength)
                    {
                        int diff = posLength - loadLength;
                        var last = loadVectors[loadLength - 1];
                        for (int i = 0; i < diff; i++)
                        {
                            loadVectors.Add(last);
                        }
                    }

                }
                
                for (int i = 0; i < loadVectors.Count; i++)
                {
                    //int nodeIndex = GetClosestNodeIndex(smartMesh.Nodes, loadPosition[i]);
                    //int nodeIndex = GetClosestNodeIndex(mesh)
                    int nodeIndex = GetClosestNodeIndex(nodePts, loadPosition[i]);
                    if (nodeIndex == -1)
                    {
                        continue; // the point is not on the mesh
                    }

                    // Deconstruct load vector
                    double xLoad = loadVectors[i].X;
                    double yLoad = loadVectors[i].Y;
                    double zLoad = loadVectors[i].Z;



                    // Construct residual force list
                    residualForces[nodeIndex * 3] += xLoad;
                    residualForces[nodeIndex * 3 + 1] += yLoad;
                    residualForces[nodeIndex * 3 + 2] += zLoad;
                    pointsWithLoad.Add(nodePts[nodeIndex]);
                }
            }
            else if(loadType == 2)
            {
                for (int i = 0; i < surfaces.Count; i++)
                {
                    foreach (Mesh mesh in newMeshList)
                    {
                        if (mesh.Vertices.Count == 8)
                        {
                            List<Point3d> pts = new List<Point3d>();
                            foreach (Point3d pt in mesh.Vertices)
                            {
                                // Find closest point on load surface
                                double u = 0;
                                double v = 0;
                                surfaces[i].ClosestPoint(pt, out u, out v);
                                if (pt.DistanceTo(surfaces[i].PointAt(u, v)) < 0.001)   // if closer than a tolerance add it to meshpoints
                                {
                                    pts.Add(pt);
                                }
                            }

                            if (pts.Count > 2)
                            {
                                NurbsSurface meshSrf = NurbsSurface.CreateFromCorners(pts[0], pts[1], pts[3], pts[2]);  // may need to do a graham scan on meshSrf
                                AreaMassProperties amp = AreaMassProperties.Compute(meshSrf);
                                double A = amp.Area;
                                Vector3d nodalForce = A * loadVectors[i] / pts.Count;

                                // Deconstruct load vector
                                double xLoad = nodalForce.X;
                                double yLoad = nodalForce.Y;
                                double zLoad = nodalForce.Z;

                                foreach (var pt in pts)
                                {
                                    int nodeIndex = GetClosestNodeIndex(nodePts, pt);
                                    if (nodeIndex == -1)
                                    {
                                        continue; // the point is not on the mesh
                                    }

                                    // Construct residual force list
                                    residualForces[nodeIndex * 3] += xLoad;
                                    residualForces[nodeIndex * 3 + 1] += yLoad;
                                    residualForces[nodeIndex * 3 + 2] += zLoad;
                                    pointsWithLoad.Add(nodePts[nodeIndex]);
                                }
                            }
                        }
                        else if (mesh.Vertices.Count == 4 || mesh.Vertices.Count == 10)
                        {
                            List<Point3d> pts = new List<Point3d>();
                            foreach (Point3d pt in mesh.Vertices)
                            {
                                double u = 0;
                                double v = 0;
                                surfaces[i].ClosestPoint(pt, out u, out v);
                                if (pt.DistanceTo(surfaces[i].PointAt(u, v)) < 0.001)
                                {
                                    pts.Add(pt);
                                }
                            }

                            if (pts.Count > 2)
                            {
                                NurbsSurface meshSrf = NurbsSurface.CreateFromCorners(pts[0], pts[1], pts[2]);
                                AreaMassProperties amp = AreaMassProperties.Compute(meshSrf);
                                double A = amp.Area;
                                Vector3d nodalForce = A * loadVectors[i] / pts.Count;

                                // Deconstruct load vector
                                double xLoad = nodalForce.X;
                                double yLoad = nodalForce.Y;
                                double zLoad = nodalForce.Z;

                                foreach (var pt in pts)
                                {
                                    int nodeIndex = GetClosestNodeIndex(nodePts, pt);
                                    if (nodeIndex == -1)
                                    {
                                        continue; // the point is not on the mesh
                                    }

                                    // Construct residual force list
                                    residualForces[nodeIndex * 3] += xLoad;
                                    residualForces[nodeIndex * 3 + 1] += yLoad;
                                    residualForces[nodeIndex * 3 + 2] += zLoad;
                                    pointsWithLoad.Add(nodePts[nodeIndex]);
                                }
                            }
                        }
                        else if (mesh.Vertices.Count == 20)
                        {
                            int method = 2;
                            List<Point3d> pts = new List<Point3d>();
                            foreach (Point3d pt in mesh.Vertices)
                            {
                                // Find closest point on load surface
                                double u = 0;
                                double v = 0;
                                surfaces[i].ClosestPoint(pt, out u, out v);
                                if (pt.DistanceTo(surfaces[i].PointAt(u, v)) < 0.001)   // if closer than a tolerance add it to meshpoints
                                {
                                    pts.Add(pt);

                                }
                            }

                            if (pts.Count > 4 && method == 0)
                            {
                                int order = 2;
                                string eltype = "Q8";

                                // first, get the global coordinates
                                LA.Matrix<double> globalCoordinates = LA.Matrix<double>.Build.Dense(8, 3);


                                List<int> nodal_ID = new List<int>();
                                for (int j = 0; j < mesh.Vertices.Count; j++)
                                {
                                    Point3d mesh_pt = mesh.Vertices[j];
                                    for (int k = 0; k < pts.Count; k++)
                                    {
                                        if (mesh_pt.DistanceTo(pts[k]) < 0.001)
                                        {
                                            nodal_ID.Add(j);
                                        }
                                        globalCoordinates[k, 0] = pts[k].X;
                                        globalCoordinates[k, 1] = pts[k].Y;
                                        globalCoordinates[k, 2] = pts[k].Z;
                                    }
                                }

                                var gaussCoordinates = FEM_Utility.GetGaussPointMatrix(order, eltype);

                                for (int j = 0; j < gaussCoordinates.RowCount; j++)
                                {
                                    double r = gaussCoordinates[j, 0];
                                    double s = gaussCoordinates[j, 1];

                                    double a = 0.25; double b = 0.5;
                                    var shapeVec = new CSD.DenseMatrix(1, 8);

                                    var derivativeMatrix = LA.Matrix<double>.Build.Dense(3, 8);

                                    double[] natCoordsArray = new double[]
                                        {
                                            -1, -1,
                                            1, -1,
                                            1, 1,
                                            -1, 1,
                                            0, -1,
                                            1, 0,
                                            0, 1,
                                            -1, 0
                                        };

                                    var naturalCoordinates = LA.Matrix<double>.Build.DenseOfRowMajor(8, 2, natCoordsArray);

                                    for (int k = 0; k < naturalCoordinates.RowCount; k++)
                                    {
                                        var x = naturalCoordinates[k, 0];
                                        var y = naturalCoordinates[k, 1];

                                        if (k < 4)
                                        {
                                            shapeVec[0, k] = a * (1 + r * x) * (1 + s * y) * (r * x + s * y - 1);

                                            derivativeMatrix[0, k] = b * x * (1 + s * y) * (r * x + b * s * y);
                                            derivativeMatrix[1, k] = b * y * (1 + r * x) * (s * y + b * r * x);
                                            derivativeMatrix[2, k] = 0.00001;
                                        }
                                        else if (k == 4 || k == 6)
                                        {
                                            shapeVec[0, k] = b * (1 - Math.Pow(r, 2)) * (1 - s * y);

                                            derivativeMatrix[0, k] = -r * (1 + s * y);
                                            derivativeMatrix[1, k] = 0.5 * x * (1 - Math.Pow(s, 2));
                                            derivativeMatrix[2, k] = 0.00001;
                                        }
                                        else if (k == 5 || k == 7)
                                        {
                                            shapeVec[0, k] = b * (1 + r * x) * (1 - Math.Pow(s, 2));

                                            derivativeMatrix[0, k] = b * y * (1 - Math.Pow(r, 2));
                                            derivativeMatrix[1, k] = -s * (1 + x * r);
                                            derivativeMatrix[2, k] = 0.00001;
                                        }
                                    }

                                    var jacobianOperator = derivativeMatrix.Multiply(globalCoordinates);
                                    double jacobianDeterminant = jacobianOperator.Determinant();

                                    totalArea += jacobianDeterminant;

                                    var shapeMat = FEM_Utility.DisplacementInterpolationMatrix(shapeVec, 3);
                                    var loadvec = LA.Vector<double>.Build.DenseOfArray(new double[] { loadVectors[i].X, loadVectors[i].Y, loadVectors[i].Z });

                                    var gaussPointLoadVector = shapeMat.TransposeThisAndMultiply(loadvec);

                                    double w_r = 1; double w_s = 1; double w_t = 1;

                                    var globalLoadVector = gaussPointLoadVector.Multiply(jacobianDeterminant) * w_r * w_s * w_t;

                                    int cnt = 0;
                                    foreach (Point3d pt in pts)
                                    {
                                        int globalID = GetClosestNodeIndex(nodePts, pt);
                                        int localID = cnt;

                                        // Construct residual force list
                                        residualForces[globalID * 3] += globalLoadVector[localID * 3];
                                        residualForces[globalID * 3 + 1] += globalLoadVector[localID * 3 + 1];
                                        residualForces[globalID * 3 + 2] += globalLoadVector[localID * 3 + 2];
                                        pointsWithLoad.Add(nodePts[globalID]);
                                        cnt++;

                                    }
                                }


                            }

                            else if (pts.Count > 4 && method == 1)
                            {
                                int order = 3;
                                string eltype = "Hex20";

                                // first, get the global coordinates
                                LA.Matrix<double> globalCoordinates = LA.Matrix<double>.Build.Dense(20, 3);


                                List<int> nodal_ID = new List<int>();
                                for (int j = 0; j < mesh.Vertices.Count; j++)
                                {
                                    Point3d mesh_pt = mesh.Vertices[j];
                                    for (int k = 0; k < pts.Count; k++)
                                    {
                                        if (mesh_pt.DistanceTo(pts[k]) < 0.001)
                                        {
                                            nodal_ID.Add(j);
                                        }
                                    }
                                }

                                for (int j = 0; j < pts.Count; j++)
                                {
                                    int ID = nodal_ID[j];
                                    globalCoordinates[ID, 0] = pts[j].X;
                                    globalCoordinates[ID, 1] = pts[j].Y;
                                    globalCoordinates[ID, 2] = pts[j].Z;
                                }

                                var gaussCoordinates = FEM_Utility.GetGaussPointMatrix(order, eltype);

                                NurbsSurface meshSrf = NurbsSurface.CreateFromCorners(pts[0], pts[1], pts[3], pts[2]);  // may need to do a graham scan on meshSrf
                                AreaMassProperties amp = AreaMassProperties.Compute(meshSrf);
                                double A = (amp.Area);



                                foreach (int j in nodal_ID)
                                {
                                    double r = gaussCoordinates[j, 0];
                                    double s = gaussCoordinates[j, 1];
                                    double t = gaussCoordinates[j, 2];

                                    double a = 0.125; double b = 0.25;
                                    var shapeVec = new CSD.DenseMatrix(1, 20);


                                    double c1 = 0.125; double c2 = 0.25; double c3 = 0.5;
                                    var derivativeMatrix = LA.Matrix<double>.Build.Dense(3, 20);

                                    var ind1 = new List<int>() { 0, 1, 2, 3, 4, 5, 6, 7 };
                                    var ind2 = new List<int>() { 8, 10, 12, 14 };
                                    var ind3 = new List<int>() { 9, 11, 13, 15 };
                                    var ind4 = new List<int>() { 16, 17, 18, 19 };

                                    int gp = 1;

                                    double[] gaussArr = new double[]
                                        {
                                            -gp, -gp, -gp,
                                            gp, -gp, -gp,
                                            gp, gp, -gp,
                                            -gp, gp, -gp,
                                            -gp, -gp, gp,
                                            gp, -gp, gp,
                                            gp, gp, gp,
                                            -gp, gp, gp,
                                            0, -gp, -gp,
                                            gp, 0, -gp,
                                            0, gp, -gp,
                                            -gp, 0, -gp,
                                            0, -gp, gp,
                                           gp, 0, gp,
                                            0, gp, gp,
                                            -gp, 0, gp,
                                            -gp, -gp, 0,
                                            gp, -gp, 0,
                                            gp, gp, 0,
                                            -gp, gp, 0
                                        };

                                    var naturalCoordinates = LA.Matrix<double>.Build.DenseOfRowMajor(20, 3, gaussArr);

                                    foreach (int ID in nodal_ID)
                                    {
                                        var natcor_r = naturalCoordinates[ID, 0];
                                        var natcor_s = naturalCoordinates[ID, 1];
                                        var natcor_t = naturalCoordinates[ID, 2];

                                        if (ind1.Contains(ID))
                                        {
                                            shapeVec[0, ID] = a * (1 + natcor_r * r) * (1 + natcor_s * s) * (1 + natcor_t * t) * (natcor_r * r + natcor_s * s + natcor_t * t - 2);

                                            derivativeMatrix[0, ID] = c1 * natcor_r * (1 + natcor_s * s) * (1 + natcor_t * t) * (2 * natcor_r * r + natcor_s * s + natcor_t * t - 1);
                                            derivativeMatrix[1, ID] = c1 * natcor_s * (1 + natcor_r * r) * (1 + natcor_t * t) * (natcor_r * r + 2 * natcor_s * s + natcor_t * t - 1);
                                            derivativeMatrix[2, ID] = c1 * natcor_t * (1 + natcor_r * r) * (1 + natcor_s * s) * (natcor_r * r + natcor_s * s + 2 * natcor_t * t - 1);
                                        }
                                        else if (ind2.Contains(ID))
                                        {
                                            shapeVec[0, ID] = b * (1 - Math.Pow(r, 2)) * (1 + natcor_s * s) * (1 + natcor_t * t);

                                            derivativeMatrix[0, ID] = -c3 * r * (1 + natcor_s * s) * (1 + natcor_t * t);
                                            derivativeMatrix[1, ID] = c2 * natcor_s * (1 - Math.Pow(r, 2)) * (1 + natcor_t * t);
                                            derivativeMatrix[2, ID] = c2 * natcor_t * (1 - Math.Pow(r, 2)) * (1 + natcor_s * s);
                                        }
                                        else if (ind3.Contains(ID))
                                        {
                                            shapeVec[0, ID] = b * (1 - Math.Pow(s, 2)) * (1 + natcor_r * r) * (1 + natcor_t * t);

                                            derivativeMatrix[0, ID] = c2 * natcor_r * (1 - Math.Pow(s, 2)) * (1 + natcor_t * t);
                                            derivativeMatrix[1, ID] = -c3 * s * (1 + natcor_r * r) * (1 + natcor_t * t);
                                            derivativeMatrix[2, ID] = c2 * natcor_t * (1 - Math.Pow(s, 2)) * (1 + natcor_r * r);
                                        }
                                        else if (ind4.Contains(ID))
                                        {
                                            shapeVec[0, ID] = b * (1 - Math.Pow(t, 2)) * (1 + natcor_r * r) * (1 + natcor_s * s);

                                            derivativeMatrix[0, ID] = c2 * natcor_r * (1 - Math.Pow(t, 2)) * (1 + natcor_s * s);
                                            derivativeMatrix[1, ID] = c2 * natcor_s * (1 - Math.Pow(t, 2)) * (1 + natcor_r * r);
                                            derivativeMatrix[2, ID] = -c3 * t * (1 + natcor_r * r) * (1 + natcor_s * s);
                                        }
                                    }

                                    var jacobianOperator = derivativeMatrix.Multiply(globalCoordinates);
                                    double jacobianDeterminant = jacobianOperator.Determinant();

                                    totalArea += jacobianDeterminant;

                                    var shapeMat = FEM_Utility.DisplacementInterpolationMatrix(shapeVec, 3);
                                    var loadvec = LA.Vector<double>.Build.DenseOfArray(new double[] { loadVectors[i].X, loadVectors[i].Y, loadVectors[i].Z });


                                    var gaussPointLoadVector = shapeMat.TransposeThisAndMultiply(loadvec);


                                    //double w_r = 1; double w_s = 1; double w_t = 1;

                                    double w_r = 0.555556; double w_s = 0.555556; double w_t = 0.555556;

                                    if (r == 0)
                                    {
                                        w_r = 0.888889;
                                    }
                                    if (s == 0)
                                    {
                                        w_s = 0.888889;
                                    }
                                    if (t == 0)
                                    {
                                        w_t = 0.888889;
                                    }


                                    var globalLoadVector = gaussPointLoadVector.Multiply(A) * w_r * w_s * w_t;

                                    int cnt = 0;
                                    foreach (Point3d pt in pts)
                                    {
                                        int globalID = GetClosestNodeIndex(nodePts, pt);
                                        int localID = nodal_ID[cnt];

                                        // Construct residual force list
                                        residualForces[globalID * 3] += globalLoadVector[localID * 3];
                                        residualForces[globalID * 3 + 1] += globalLoadVector[localID * 3 + 1];
                                        residualForces[globalID * 3 + 2] += globalLoadVector[localID * 3 + 2];
                                        pointsWithLoad.Add(nodePts[globalID]);
                                        cnt++;

                                    }

                                }

                            }

                            else if (pts.Count > 4 && method == 2)
                            {
                                NurbsSurface meshSrf = NurbsSurface.CreateFromCorners(pts[0], pts[1], pts[2], pts[3]);  // may need to do a graham scan on meshSrf
                                AreaMassProperties amp = AreaMassProperties.Compute(meshSrf);
                                double A = amp.Area;
                                totalArea += A;

                                int cnt = 0;

                                foreach (var pt in pts)
                                {
                                    double xLoad = 0;
                                    double yLoad = 0;
                                    double zLoad = 0;
                                    if (cnt < 4)
                                    {
                                        Vector3d nodalForce = -A * loadVectors[i] / 12;

                                        // Deconstruct load vector
                                        xLoad = nodalForce.X;
                                        yLoad = nodalForce.Y;
                                        zLoad = nodalForce.Z;
                                    }
                                    else
                                    {
                                        Vector3d nodalForce = A * loadVectors[i] / 3;

                                        // Deconstruct load vector
                                        xLoad = nodalForce.X;
                                        yLoad = nodalForce.Y;
                                        zLoad = nodalForce.Z;
                                    }
                                    

                                    int nodeIndex = GetClosestNodeIndex(nodePts, pt);
                                    if (nodeIndex == -1)
                                    {
                                        continue; // the point is not on the mesh
                                    }

                                    // Construct residual force list
                                    residualForces[nodeIndex * 3] += xLoad;
                                    residualForces[nodeIndex * 3 + 1] += yLoad;
                                    residualForces[nodeIndex * 3 + 2] += zLoad;
                                    pointsWithLoad.Add(nodePts[nodeIndex]);

                                    cnt++;
                                }
                            }
                            else if (pts.Count > 0 && method == 3)
                            {
                                NurbsSurface meshSrf = NurbsSurface.CreateFromCorners(pts[0], pts[1], pts[2], pts[3]);  // may need to do a graham scan on meshSrf
                                AreaMassProperties amp = AreaMassProperties.Compute(meshSrf);
                                double A = amp.Area;
                                totalArea += A;

                                int cnt = 0;

                                foreach (var pt in pts)
                                {
                                    double xLoad = 0;
                                    double yLoad = 0;
                                    double zLoad = 0;

                                    Vector3d nodalForce = -A * loadVectors[i] / 8;

                                    // Deconstruct load vector
                                    xLoad = nodalForce.X;
                                    yLoad = nodalForce.Y;
                                    zLoad = nodalForce.Z;
                                    
                                    int nodeIndex = GetClosestNodeIndex(nodePts, pt);
                                    if (nodeIndex == -1)
                                    {
                                        continue; // the point is not on the mesh
                                    }

                                    // Construct residual force list
                                    residualForces[nodeIndex * 3] += xLoad;
                                    residualForces[nodeIndex * 3 + 1] += yLoad;
                                    residualForces[nodeIndex * 3 + 2] += zLoad;
                                    pointsWithLoad.Add(nodePts[nodeIndex]);

                                    cnt++;
                                }
                            }
                        }
                    }
                }

            }

            List<double> res = new List<double>();
            res.Add( residualForces.Sum());
            res.Add(totalArea);
            
            // Output

            DA.SetDataList(0, residualForces);
            DA.SetDataList(1, pointsWithLoad);
            DA.SetDataList(2, res);
        }

        #region Methods

        /// <summary>
        /// Get the index of the node closest to the load position.
        /// </summary>
        /// <returns> Index of node closest to load position.</returns>
        private int GetClosestNodeIndex(List<Point3d> nodePoints, Point3d loadPosition)
        {
            int nodeIndex = -1;

            for (int i = 0; i < nodePoints.Count; i++)
            {
                double distance = nodePoints[i].DistanceToSquared(loadPosition);
                if (distance < 0.00001) // set tolerance
                {
                    nodeIndex = i;
                    break;
                }
            }
            return nodeIndex;
        }

        /// <summary>
        /// Get the index of the nodes on the surface.
        /// </summary>
        /// <returns> List of index of nodes on the surface.</returns>
        private List<int> GetNodeIndexOnSurface(List<Node> nodes, BrepFace surface)
        {
            List<int> nodeIndexOnSurface = new List<int>();
            for (int i = 0; i < nodes.Count; i++)
            {
                if (nodes[i].IsOnFace(surface))
                {
                    nodeIndexOnSurface.Add(i);
                }
            }
            return nodeIndexOnSurface;
        }


        #endregion

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return Properties.Resources.FELoad;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("e933626f-742e-42d6-bb97-30c4e2cca685"); }
        }
    }
}