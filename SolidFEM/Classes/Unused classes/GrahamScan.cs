using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace SolidFEM
{
    public static class GrahamScan
    {
        public static Point3f p0;
        private static Mesh SortedMesh(List<Point3f> pts)
        {
            Mesh mesh = new Mesh();

            // add vertices
            foreach (Point3f p in pts)
            {
                mesh.Vertices.Add(p);
            }

            // add faces
            mesh.Faces.AddFace(0, 1, 2, 3);
            mesh.Faces.AddFace(4, 5, 6, 7);
            mesh.Faces.AddFace(0, 1, 5, 4);
            mesh.Faces.AddFace(1, 2, 6, 5);
            mesh.Faces.AddFace(2, 3, 7, 6);
            mesh.Faces.AddFace(3, 0, 4, 7);

            return mesh;
        }

        private static List<List<Point3f>> GetTopAndBottomNodes(Mesh mesh)
        {

            var returnList = new List<List<Point3f>>();
            // -- find the lowest and highest face of the mesh --
            int lowface = 0;
            int highFace = 0;

            // get the first centroid:
            var p0 = mesh.TopologyVertices[mesh.Faces[0].A];
            var p1 = mesh.TopologyVertices[mesh.Faces[0].B];
            var p2 = mesh.TopologyVertices[mesh.Faces[0].C];
            var p3 = mesh.TopologyVertices[mesh.Faces[0].D];

            List<Point3d> pts = new List<Point3d>() { p0, p1, p2, p3 };
            double x = 0; double y = 0; double z = 0;
            foreach (Point3d point in pts)
            {
                x += point.X;
                y += point.Y;
                z += point.Z;
            }
            Point3d currentLowFaceCenter = new Point3d(x / pts.Count, y / pts.Count, z / pts.Count);
            Point3d currentHigFaceCenter = new Point3d(x / pts.Count, y / pts.Count, z / pts.Count);

            for (int i = 1; i < mesh.Faces.Count; i++)
            {
                // get the first centroid:
                p0 = mesh.TopologyVertices[mesh.Faces[i].A];
                p1 = mesh.TopologyVertices[mesh.Faces[i].B];
                p2 = mesh.TopologyVertices[mesh.Faces[i].C];
                p3 = mesh.TopologyVertices[mesh.Faces[i].D];

                pts = new List<Point3d>() { p0, p1, p2, p3 };
                x = 0; y = 0; z = 0;
                foreach (Point3d point in pts)
                {
                    x += point.X;
                    y += point.Y;
                    z += point.Z;
                }
                Point3d faceCenter = new Point3d(x / pts.Count, y / pts.Count, z / pts.Count);

                if (faceCenter.Z > currentHigFaceCenter.Z)
                {
                    highFace = i;
                    currentHigFaceCenter = new Point3d(faceCenter);

                }
                else if (faceCenter.Z < currentLowFaceCenter.Z)
                {
                    lowface = i;
                    currentLowFaceCenter = new Point3d(faceCenter);
                }
            }

            // -- get the vertices from the lowest and highest face --
            List<Point3f> lowPts = new List<Point3f>() { mesh.TopologyVertices[mesh.Faces[lowface].A], mesh.TopologyVertices[mesh.Faces[lowface].B], mesh.TopologyVertices[mesh.Faces[lowface].C], mesh.TopologyVertices[mesh.Faces[lowface].D] };
            List<Point3f> highPts = new List<Point3f>() { mesh.TopologyVertices[mesh.Faces[highFace].A], mesh.TopologyVertices[mesh.Faces[highFace].B], mesh.TopologyVertices[mesh.Faces[highFace].C], mesh.TopologyVertices[mesh.Faces[highFace].D] };

            returnList.Add(lowPts);
            returnList.Add(highPts);

            return returnList;

        }

        public static Mesh DoGrahamScan(Mesh oldMesh)
        {
            // find top and bottom mesh
            List<List<Point3f>> topAndBottom = GetTopAndBottomNodes(oldMesh);

            List<Point3f> bottomPts = topAndBottom[0];
            List<Point3f> topPts = topAndBottom[1];

            List<Point3f> sortedBottom = GrahamScanFace(bottomPts);
            List<Point3f> sortedTop = GrahamScanFace(topPts);
            if (sortedBottom == null || sortedTop == null)
            {
                return oldMesh;
            }

            List<Point3f> sortedPts = new List<Point3f>( sortedBottom);
            sortedPts.AddRange(sortedTop);
            // with the sorted points, create an updated mesh
            Mesh mesh = SortedMesh(sortedPts);

            return mesh;

        }
        public static Mesh DoGrahamScanOriginal(Mesh oldMesh)
        {

            // find top and bottom mesh
            List<List<Point3f>> topAndBottom = GetTopAndBottomNodes(oldMesh);

            List<Point3f> bottomPts = topAndBottom[0];
            List<Point3f> topPts = topAndBottom[1];


            List<Point3f> pts = oldMesh.TopologyVertices.ToList();


            int n = pts.Count();
            Stack<Point3f> pointStack = new Stack<Point3f>();
            List<Point3f> sortedPts = new List<Point3f>();

            // find point with lowest y- and x- coordinate:
            var sortedXY = pts.OrderBy(pt => pt.Y).ThenBy(pt => pt.X).ToList();
            p0 = sortedXY[0]; // the first point
                                     // remove the first point from the sorted point.
                                     //sortedXY.RemoveAt(0);

            // sort point by polar angle with P0
            sortedXY.Sort(Compare);

            // in the case of two or more points bein colinear, we only want to keep the one at the furthest distance
            int m = 1;
            for (int i = 0; i < n; i++)
            {

                while ((i < n - 1) && (Orientation(p0, sortedXY[i], sortedXY[i + 1])) == 0)
                {
                    i += 1;
                }
                sortedXY[m] = sortedXY[i];
                m += 1;
            }

            if (m < 3) return null; // only three points exists, return

            // create empty stack and push first three points to it
            Stack<Point3f> stack = new Stack<Point3f>();
            stack.Push(sortedXY[0]);
            stack.Push(sortedXY[1]);
            stack.Push(sortedXY[2]);

            

            // process for remaing n-3 points
            for (int i = 3; i < m; i++)
            {
                // keep removing top while the angle formed by points next-to-top and sortedXY[i] makes a non-left turn
                while (stack.Count > 1 && Orientation(NextToPop(stack), stack.Peek(), sortedXY[i]) != 2)
                {
                    stack.Pop();
                }
                stack.Push(sortedXY[i]);
            }

            sortedPts = stack.ToList();
            sortedPts.Reverse();


            // with the sorted points, create an updated mesh
            Mesh mesh = SortedMesh(sortedPts);

            return mesh;
        }

        private static List<Point3f> GrahamScanFace(List<Point3f> pts)
        {

            int n = pts.Count();
            Stack<Point3f> pointStack = new Stack<Point3f>();
            List<Point3f> sortedPts = new List<Point3f>();

            // find point with lowest y- and x- coordinate:
            var sortedXY = pts.OrderBy(pt => pt.Y).ThenBy(pt => pt.X).ToList();
            p0 = sortedXY[0]; // the first point
                              // remove the first point from the sorted point.
                              //sortedXY.RemoveAt(0);

            // sort point by polar angle with P0
            sortedXY.Sort(Compare);

            // in the case of two or more points bein colinear, we only want to keep the one at the furthest distance
            int m = 1;
            for (int i = 0; i < n; i++)
            {

                while ((i < n - 1) && (Orientation(p0, sortedXY[i], sortedXY[i + 1])) == 0)
                {
                    i += 1;
                }
                sortedXY[m] = sortedXY[i];
                m += 1;
            }

            if (m < 3) return null; // only three points exists, return

            // create empty stack and push first three points to it
            Stack<Point3f> stack = new Stack<Point3f>();
            stack.Push(sortedXY[0]);
            stack.Push(sortedXY[1]);
            stack.Push(sortedXY[2]);



            // process for remaing n-3 points
            for (int i = 3; i < m; i++)
            {
                // keep removing top while the angle formed by points next-to-top and sortedXY[i] makes a non-left turn
                while (stack.Count > 1 && Orientation(NextToPop(stack), stack.Peek(), sortedXY[i]) != 2)
                {
                    stack.Pop();
                }
                stack.Push(sortedXY[i]);
            }

            sortedPts = stack.ToList();
            sortedPts.Reverse();

            return sortedPts;
        }
        private static Point3f NextToPop(Stack<Point3f> stack)
        {
            Point3f p = stack.Pop();
            Point3f next = stack.Pop();

            // replace the elements in stack
            stack.Push(next);
            stack.Push(p);
            return next;
        }

        private static void Swap(ref Point3f p1, ref Point3f p2)
        {
            Point3f temp = new Point3f(p1.X, p1.Y, p1.Z);
            p1 = p2;
            p2 = temp;
        }
        private static int Compare(Point3f p1, Point3f p2)
        {
            // find orientation
            var o = Orientation(p0, p1, p2);

            if (o == 0)
            {
                if (p0.DistanceTo(p2) >= p0.DistanceTo(p1)) return -1;
                else return 1;
            }
            else
            {
                if (o == 2) return -1;
                else return 1;
            }

        }

        private static int Orientation(Point3f p1, Point3f p2, Point3f p3)
        {
            double val = (p2.Y - p1.Y) * (p3.X - p2.X) - (p2.X - p1.X) * (p3.Y - p2.Y);
            if (val == 0) return 0;

            if (val > 0) return 1;
            else return 2;

        }


    }
}
