using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace SolidFEM.Classes
{
    class Geometry
    {
        public Brep Brep { get; set; }
        public List<BrepFace> Faces { get; set; }
        public List<BrepEdge> Edges { get; set; }
        public List<BrepVertex> Vertices { get; set; }

        // Constructors
        public Geometry()
        {
            // empty constructor
        }
        public Geometry(Brep _brep, List<BrepFace> _faces, List<BrepEdge> _edges, List<BrepVertex> _vertices) // surface and sweep
        {
            Brep = _brep;
            Faces = _faces;
            Edges = _edges;
            Vertices = _vertices;
        }
        public Geometry(Brep _brep, int bottomFace) // solid
        {
            Brep = _brep;
            Faces = SortBrepFaces(_brep, bottomFace);
            Edges = SortBrepEdges(_brep, bottomFace);
            Vertices = SortBrepVertex(_brep, bottomFace);
        }
        private List<BrepFace> SortBrepFaces(Brep brep, int bottomFace)
        {
            List<BrepFace> faceSorted = new List<BrepFace>();
            // Find top and bottom edge
            List<BrepFace> brepFace = brep.Faces.ToList();
            List<int> indexAdjecentFaces = (brepFace[bottomFace].AdjacentFaces()).ToList();
            indexAdjecentFaces.Add(bottomFace);
            for (int i = 0; i < brepFace.Count; i++)
            {
                if (!indexAdjecentFaces.Contains(brepFace.IndexOf(brepFace[i])))
                {
                    faceSorted.Add(brepFace[bottomFace]); // bottom face
                    faceSorted.Add(brepFace[i]); // top face
                    continue;
                }
            }
            indexAdjecentFaces.Remove(bottomFace);
            foreach (int index in indexAdjecentFaces) { faceSorted.Add(brepFace[index]); }
            return faceSorted;
        }
        private List<BrepEdge> SortBrepEdges(Brep brep, int bottomFace)
        {
            List<BrepFace> faceSorted = SortBrepFaces(brep, bottomFace);
            List<BrepEdge> edgeSorted = new List<BrepEdge>();
            List<int> indexAdjecentEdges = new List<int>();

            // Find edges connected to top and bottom face.
            indexAdjecentEdges.AddRange(faceSorted[0].AdjacentEdges().ToList());
            indexAdjecentEdges.AddRange(faceSorted[1].AdjacentEdges().ToList());

            // Add edges to list.
            foreach (int index in indexAdjecentEdges) { edgeSorted.Add(brep.Edges[index]); }

            // Find rest of edges
            List<BrepEdge> brepEdgesCopy = brep.Edges.ToList();
            List<Curve> rails = new List<Curve>(brepEdgesCopy);
            foreach (int index in indexAdjecentEdges) { rails.Remove(brepEdgesCopy[index]); }

            // Add rest of edges to list.
            foreach (BrepEdge edge in rails) { edgeSorted.Add(edge); }
            return edgeSorted;
        }
        public List<BrepVertex> SortBrepVertex(Brep brep, int bottomFace)
        {
            List<BrepVertex> vertexSorted = new List<BrepVertex>();
            List<BrepFace> brepFaces = SortBrepFaces(brep, bottomFace);
            List<BrepVertex> brepVertex = brep.Vertices.ToList();

            foreach (BrepVertex vertex in brepVertex)
            {
                bool isOnBottomFace = IsVertexOnFace(vertex.Location, brepFaces[0]);
                if (isOnBottomFace) { vertexSorted.Add(vertex); }
            }

            foreach (BrepVertex vertex in brepVertex)
            {
                bool isOnTopFace = IsVertexOnFace(vertex.Location, brepFaces[1]);
                if (isOnTopFace) { vertexSorted.Add(vertex); }
            }
            return vertexSorted;
        }
        public bool IsVertexOnFace(Point3d point, BrepFace face)
        {
            bool isOnFace = false;

            face.ClosestPoint(point, out double PointOnCurveU, out double PointOnCurveV);
            Point3d testPoint = face.PointAt(PointOnCurveU, PointOnCurveV);  // make test point 
            double distanceToFace = (testPoint - point).Length; // calculate distance between testPoint and node
            if (distanceToFace <= 0.0001 & distanceToFace >= -0.0001) // if distance = 0: node is on edge
            {
                isOnFace = true;
            }
            return isOnFace;
        }
    }
}
