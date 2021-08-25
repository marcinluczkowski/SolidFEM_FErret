using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using Rhino.Geometry;

namespace SolidFEM.Classes
{
    class qNode
    {
        public Point3d Coordinate { get; set; }
        public bool BoundaryNode { get; set; }
        public bool OBS { get; set; }

        // Constructors:
        public qNode()
        {
            // empty constructor
        }

        public qNode(Point3d _coordinate, bool _boundaryNode)
        {
            Coordinate = _coordinate;
            BoundaryNode = _boundaryNode;
        }

        // Methods:
        public List<qEdge> GetConnectedEdges(List<qEdge> globalEdgeList)
        {
            // summary: get connected edges to a node
            List<qEdge> connectedEdges = new List<qEdge>();
            foreach (qEdge edge in globalEdgeList)
            {
                if ((edge.StartNode == this) | (edge.EndNode == this))
                {
                    connectedEdges.Add(edge);
                }
            }
            return connectedEdges;
        }
        public List<qElement> GetQuadsConnectedToNode(List<qEdge> globalEdgeList)
        {
            // summary: get all quad elements connected to a node
            List<qElement> quadElements = new List<qElement>();
            var connectedEdges = this.GetConnectedEdges(globalEdgeList);

            foreach (qEdge edge in connectedEdges)
            {
                List<qElement> connectedElements = edge.GetConnectedElements();
                foreach (qElement element in connectedElements)
                {
                    if (element.IsQuad) { quadElements.Add(element); }
                }
            }

            // delete dublicates
            List<qElement> quadElementsNoDublicates = new List<qElement>();
            if (quadElements.Count > 0)
            {
                foreach (qElement element in quadElements)
                {
                    if (!quadElementsNoDublicates.Contains(element)) { quadElementsNoDublicates.Add(element); }
                }
            }

            return quadElementsNoDublicates;
        }
        public bool IsFrontNode(List<qEdge> frontEdges)
        {
            // summary: check if node is a front node
            bool isFrontNode = false;
            foreach (qEdge front in frontEdges)
            {
                if (this == front.StartNode | this == front.EndNode)
                {
                    isFrontNode = true;
                    break;
                }
            }
            return isFrontNode;
        }
        public List<qElement> GetTrianglesConnectedToNode(List<qEdge> globalEdgeList)
        {
            // summary: get all quad elements connected to a node
            List<qElement> triangleElements = new List<qElement>();
            var connectedEdges = this.GetConnectedEdges(globalEdgeList);

            foreach (qEdge edge in connectedEdges)
            {
                List<qElement> connectedElements = edge.GetConnectedElements();
                foreach (qElement element in connectedElements)
                {
                    if (!element.IsQuad) { triangleElements.Add(element); }
                }
            }

            // delete dublicates
            List<qElement> triangleElementsNoDublicates = new List<qElement>();
            if (triangleElements.Count > 0)
            {
                foreach (qElement element in triangleElements)
                {
                    if (!triangleElementsNoDublicates.Contains(element)) { triangleElementsNoDublicates.Add(element); }
                }
            }

            return triangleElementsNoDublicates;
        }
        public List<qElement> GetConnectedElements(List<qEdge> globalEdgeList)
        {
            qNode node = this;
            List<qElement> connectedElements = new List<qElement>();
            List<qEdge> connectedEdges = node.GetConnectedEdges(globalEdgeList);
            foreach (qEdge edge in connectedEdges)
            {
                if (!connectedElements.Contains(edge.Element1)) { connectedElements.Add(edge.Element1); }
                if (!connectedElements.Contains(edge.Element2) & edge.Element2 != null) { connectedElements.Add(edge.Element2); }
            }
            return connectedElements;
        }
        public bool IsOnEdge(BrepEdge edge)
        {
            Point3d point = this.Coordinate;
            bool isOnEdge = false;

            edge.ClosestPoint(point, out double PointOnCurve);
            Point3d testPoint = edge.PointAt(PointOnCurve);  // make test point 
            double distanceToEdge = (testPoint - point).Length; // calculate distance between testPoint and node
            if (distanceToEdge <= 0.0001 & distanceToEdge >= -0.0001) // if distance = 0: node is on edge
            {
                isOnEdge = true;
            }
            return isOnEdge;
        }
    }
}
