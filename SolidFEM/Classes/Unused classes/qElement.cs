using System;
using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;

namespace SolidFEM.Classes
{
    class qElement
    {
        public List<qEdge> EdgeList { get; set; }
        public List<double> AngleList { get; set; } // todo: when angle is larger than pi it does not work..
        public List<Line> Contour { get; set; }
        public bool IsQuad { get; set; }
        public double DistortionMetric { get; set; }


        // Constructors:
        public qElement()
        {
            // empty constructor
        }
        public qElement(List<qEdge> _edgeList)
        {
            EdgeList = _edgeList;
            CalculateAngles();
            GetContourOfElement();

            if (_edgeList.Count == 4) { IsQuad = true; }
            else { IsQuad = false; }
            if (! IsQuad)
            {
                FixEdgeOrder();
            }
        }



        // Methods
        public List<double> CalculateAnglesOld(List<qEdge> _edgeList)// to do: slett
        {
            Vector3d vec1 = Vector3d.Zero;
            Vector3d vec2 = Vector3d.Zero;
            double ang = 0;
            List<double> angList = new List<double>();
            List<qEdge> edgeListCopy = new List<qEdge>(_edgeList);
            edgeListCopy.Add(edgeListCopy[0]);

            // todo: change calculation of angles AnalyzeTriangle (do not work for chevron)
            for (int i = 0; i < edgeListCopy.Count - 1; i++)
            {
                Point3d start1 = edgeListCopy[i].StartNode.Coordinate;
                Point3d end1 = edgeListCopy[i].EndNode.Coordinate;
                Point3d start2 = edgeListCopy[i + 1].StartNode.Coordinate;
                Point3d end2 = edgeListCopy[i + 1].EndNode.Coordinate;

                vec1 = start1 - end1;
                vec2 = end2 - start2;

                ang = Vector3d.VectorAngle(vec1, vec2); // radian
                angList.Add(ang);
            }
            return angList;
        }

        public void CalculateAngles()
        {
            List<qEdge> _edgeList =  EdgeList;
            List<double> angList = new List<double>();

            List<qEdge> edgeListCopy = new List<qEdge>(_edgeList);
            edgeListCopy.Add(edgeListCopy[0]);

            var vectors1 = _edgeList[0].CalculateVectorsFromSharedNode(_edgeList[1]);
            angList.Add(Vector3d.VectorAngle(vectors1.Item1, vectors1.Item2, Vector3d.CrossProduct(vectors1.Item1, vectors1.Item2)));

            var vectors2 = _edgeList[0].CalculateVectorsFromSharedNode(_edgeList[2]);
            angList.Add(Vector3d.VectorAngle(vectors2.Item1, vectors2.Item2, Vector3d.CrossProduct(vectors2.Item1, vectors2.Item2)));

            if (_edgeList.Count == 3)
            {
                var vectors3Tri = _edgeList[1].CalculateVectorsFromSharedNode(_edgeList[2]);
                angList.Add(Vector3d.VectorAngle(vectors3Tri.Item1, vectors3Tri.Item2, Vector3d.CrossProduct(vectors3Tri.Item1, vectors3Tri.Item2)));
            }
            else
            {
                var vectors3 = _edgeList[1].CalculateVectorsFromSharedNode(_edgeList[3]);
                angList.Add(Vector3d.VectorAngle(vectors3.Item1, vectors3.Item2, Vector3d.CrossProduct(vectors3.Item1, vectors3.Item2)));

                var vectors4 = _edgeList[2].CalculateVectorsFromSharedNode(_edgeList[3]);
                angList.Add(Vector3d.VectorAngle(vectors4.Item1, vectors4.Item2, Vector3d.CrossProduct(vectors4.Item1, vectors4.Item2)));
            }

             AngleList = angList;
        }

        public void GetContourOfElement()
        {
            List<Line> contour = new List<Line>();
            foreach (qEdge edge in  EdgeList)
            {
                Line contourLine = new Line(edge.StartNode.Coordinate, edge.EndNode.Coordinate);
                contour.Add(contourLine);
            }
             Contour = contour;
        }
        public Point3d GetElementCenter()
        {
            // summary: get center of an element
            double sx = 0;
            double sy = 0;
            double sz = 0;

            List<qEdge> edgeList =  EdgeList;
            foreach (qEdge edge in edgeList)
            {
                Point3d startPoint = edge.StartNode.Coordinate;
                Point3d endPoint = edge.EndNode.Coordinate;
                List<Point3d> pts = new List<Point3d>() { startPoint, endPoint };
                foreach (Point3d pt in pts)
                {
                    sx = sx + pt.X;
                    sy = sy + pt.Y;
                    sz = sz + pt.Z;
                }
            }
            int n = edgeList.Count * 2;
            Point3d centerPt = new Point3d(sx / n, sy / n, sz / n);
            return centerPt;
        }
        public List<qNode> GetNodesOfElement()
        {
            // summary: get nodes of an element
            List<qNode> nodeList = new List<qNode>();
            if (! IsQuad)
            {
                List<qEdge> triEdges =  EdgeList;
                qEdge baseEdge = triEdges[0];
                qEdge rightEdge = triEdges[1];
                qEdge leftEdge = triEdges[2];

                qNode node1 = new qNode();
                qNode node2 = new qNode();
                qNode node3 = new qNode();

                if (baseEdge.StartNode == leftEdge.StartNode | baseEdge.StartNode == leftEdge.EndNode)
                {
                    node1 = baseEdge.StartNode;
                    node2 = baseEdge.EndNode;
                    if (baseEdge.StartNode == leftEdge.StartNode) { node3 = leftEdge.EndNode; }
                    if (baseEdge.StartNode == leftEdge.EndNode) { node3 = leftEdge.StartNode; }
                }
                else if (baseEdge.EndNode == leftEdge.StartNode | baseEdge.EndNode == leftEdge.EndNode)
                {
                    node1 = baseEdge.EndNode;
                    node2 = baseEdge.StartNode;
                    if (baseEdge.EndNode == leftEdge.StartNode) { node3 = leftEdge.EndNode; }
                    if (baseEdge.EndNode == leftEdge.EndNode) { node3 = leftEdge.StartNode; }
                }
                else if (baseEdge.StartNode == rightEdge.StartNode | baseEdge.StartNode == rightEdge.EndNode)
                {
                    node2 = baseEdge.StartNode;
                    node1 = baseEdge.EndNode;
                    if (baseEdge.StartNode == rightEdge.StartNode) { node3 = rightEdge.EndNode; }
                    if (baseEdge.StartNode == rightEdge.EndNode) { node3 = rightEdge.StartNode; }
                }
                else if (baseEdge.EndNode == rightEdge.StartNode | baseEdge.EndNode == rightEdge.EndNode)
                {
                    node2 = baseEdge.EndNode;
                    node1 = baseEdge.StartNode;
                    if (baseEdge.EndNode == rightEdge.StartNode) { node3 = rightEdge.EndNode; }
                    if (baseEdge.EndNode == rightEdge.EndNode) { node3 = rightEdge.StartNode; }
                }

                nodeList = new List<qNode> { node1, node2, node3 }; // n1: bottom left, n2: bottom right, n3: top 

                /* todo: slett
                foreach (qEdge edge in  EdgeList)
                {
                    if (!nodeList.Contains(edge.StartNode))
                    {
                        nodeList.Add(edge.StartNode);
                    }
                    
                    if (!nodeList.Contains(edge.EndNode))
                    {
                        nodeList.Add(edge.EndNode);
                    }
                }*/
            }
            else if ( IsQuad)
            {
                List<qEdge> quadEdges =  EdgeList;

                qEdge baseEdge = quadEdges[0];
                qEdge rightEdge = quadEdges[1];
                qEdge leftEdge = quadEdges[2];
                qEdge topEdge = quadEdges[3];

                qNode node1 = new qNode();
                qNode node2 = new qNode();
                qNode node3 = new qNode();
                qNode node4 = new qNode();

                if (baseEdge.StartNode == leftEdge.StartNode | baseEdge.StartNode == leftEdge.EndNode)
                {
                    node1 = baseEdge.StartNode;
                    node2 = baseEdge.EndNode;
                }
                else if (baseEdge.EndNode == leftEdge.StartNode | baseEdge.EndNode == leftEdge.EndNode)
                {
                    node1 = baseEdge.EndNode;
                    node2 = baseEdge.StartNode;
                }
                else if (baseEdge.StartNode == rightEdge.StartNode | baseEdge.StartNode == rightEdge.EndNode)
                {
                    node2 = baseEdge.StartNode;
                    node1 = baseEdge.EndNode;
                }
                else if (baseEdge.EndNode == rightEdge.StartNode | baseEdge.EndNode == rightEdge.EndNode)
                {
                    node2 = baseEdge.EndNode;
                    node1 = baseEdge.StartNode;
                }

                if (topEdge.StartNode == leftEdge.StartNode | topEdge.StartNode == leftEdge.EndNode)
                {
                    node3 = topEdge.EndNode;
                    node4 = topEdge.StartNode;
                }
                else if (topEdge.EndNode == leftEdge.StartNode | topEdge.EndNode == leftEdge.EndNode)
                {
                    node3 = topEdge.StartNode;
                    node4 = topEdge.EndNode;
                }
                else if (topEdge.StartNode == rightEdge.StartNode | topEdge.StartNode == rightEdge.EndNode)
                {
                    node4 = topEdge.EndNode;
                    node3 = topEdge.StartNode;
                }
                else if (topEdge.EndNode == rightEdge.StartNode | topEdge.EndNode == rightEdge.EndNode)
                {
                    node4 = topEdge.StartNode;
                    node3 = topEdge.EndNode;
                }

                nodeList = new List<qNode> { node1, node2, node3, node4 }; // n1: bottom left, n2: bottom right, n3: top right, n4: top left
            }
            return nodeList;
        }
        public void FixEdgeOrder()
        {
            // summary: fix edge order of triangle elements

            qEdge edge =  EdgeList[0]; // assume this to be the E_front
            qEdge edgeConnectedToStartNode = new qEdge();
            qEdge edgeConnectedToEndNode = new qEdge();
            qEdge edgeNotConnected = new qEdge();

            for (int i = 1; i < EdgeList.Count; i++)
            {
                if (edge.StartNode == EdgeList[i].StartNode | edge.StartNode == EdgeList[i].EndNode)
                {
                    edgeConnectedToStartNode = EdgeList[i];
                }
                else if (edge.EndNode == EdgeList[i].StartNode | edge.EndNode == EdgeList[i].EndNode)
                {
                    edgeConnectedToEndNode = EdgeList[i];
                }
                else
                {
                    edgeNotConnected = EdgeList[i];
                }
            }
            Point3d midPointEdg = 0.5 * (edge.StartNode.Coordinate + edge.EndNode.Coordinate); // mid point of edge
            Point3d centerPoint = GetElementCenter();
            Vector3d centerToMidVector = midPointEdg - centerPoint;
            Vector3d centerToEndNodeVector = edge.EndNode.Coordinate - centerPoint;
            Vector3d centerToStartNodeVector = edge.StartNode.Coordinate - centerPoint;
            double startAngle = Vector3d.VectorAngle(centerToMidVector, centerToStartNodeVector, Vector3d.ZAxis); // todo: make normal more general
            double endAngle = Vector3d.VectorAngle(centerToMidVector, centerToEndNodeVector, Vector3d.ZAxis); // todo: make normal more general

            if (endAngle < startAngle)
            {
                if (! IsQuad)
                {
                     EdgeList = new List<qEdge>() { edge, edgeConnectedToEndNode, edgeConnectedToStartNode };
                }
                else
                {
                     EdgeList = new List<qEdge>() { edge, edgeConnectedToEndNode, edgeConnectedToStartNode, edgeNotConnected };
                }
            }
            else
            {
                if (! IsQuad)
                {
                     EdgeList = new List<qEdge>() { edge, edgeConnectedToStartNode, edgeConnectedToEndNode };
                }
                else
                {
                     EdgeList = new List<qEdge>() { edge, edgeConnectedToStartNode, edgeConnectedToEndNode, edgeNotConnected };
                }
            }



            /* to do: fixslett
            var sideEdges = edge.OrientateNeigborEdges(edgeConnectedToStartNode, edgeConnectedToEndNode);
            qEdge leftSide = sideEdges.Item1;
            qEdge rightSide = sideEdges.Item2;

            if (! IsQuad)
            {
                 EdgeList = new List<qEdge>() { edge, rightSide, leftSide };
            }
            else
            {
                 EdgeList = new List<qEdge>() { edge, rightSide, leftSide, edgeNotConnected };
            }*/

        }

        public bool IsChevron()
        {
            if ( AngleList.Max() > (double)200 / (double)180 * Math.PI) { return true; }
            else { return false; }
        }
        public bool IsInverted()
        {
            // summary: check if a triangle or quad element is inverted
            bool isInverted = false;

            if (! IsQuad)
            {
                Point3d A =  EdgeList[0].GetSharedNode( EdgeList[1]).Coordinate;
                Point3d B =  EdgeList[1].GetSharedNode( EdgeList[2]).Coordinate;
                Point3d C =  EdgeList[2].GetSharedNode( EdgeList[0]).Coordinate;

                // check area
                double area = 0.5 * (A.X * (B.Y - C.Y) + B.X * (C.Y - A.Y) + C.X * (A.Y - B.Y));
                if (area <= 0) { isInverted = true; }
            }
            else // todo: fix inverted elements for quads.
            {
                List<qNode> elementNodes =  GetNodesOfElement();
                NurbsCurve line1 = new Line(elementNodes[0].Coordinate, elementNodes[3].Coordinate).ToNurbsCurve();
                NurbsCurve line2 = new Line(elementNodes[1].Coordinate, elementNodes[2].Coordinate).ToNurbsCurve();
                var placesWithIntersection = Intersection.CurveCurve(line1, line2, 0.0001, 0.0001);
                if (placesWithIntersection.Count == 0)
                {
                    isInverted = false;
                }
                else
                {
                    isInverted = true;
                }
                /*
                var nodes =  GetNodesOfElement();
                Point3d A1 = nodes[0].Coordinate;
                Point3d B1 = nodes[1].Coordinate;
                Point3d C1 = nodes[3].Coordinate;

                Point3d A2 = nodes[1].Coordinate;
                Point3d B2 = nodes[2].Coordinate;
                Point3d C2 = nodes[3].Coordinate;

                // check area
                double area1 = 0.5 * (A1.X * (B1.Y - C1.Y) + B1.X * (C1.Y - A1.Y) + C1.X * (A1.Y - B1.Y));
                double area2 = 0.5 * (A2.X * (B2.Y - C2.Y) + B2.X * (C2.Y - A2.Y) + C2.X * (A2.Y - B2.Y));
                if (area1 <= 0 | area2 <= 0) { isInverted = true; }*/
            }
            return isInverted;
        }
    }
}
