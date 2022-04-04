using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace SolidFEM.Classes
{
    class qEdge
    {
        public qNode StartNode { get; set; }
        public qNode EndNode { get; set; }
        public double Length { get; set; }
        public Line EdgeLine { get; set; } // for visualization
        public qElement Element1 { get; set; }
        public qElement Element2 { get; set; }
        public qEdge LeftFrontNeighbor { get; set; }
        public qEdge RightFrontNeighbor { get; set; }
        public int Level { get; set; } // level if front edge
        public bool Unselectable { get; set; }



        // Constructors:
        public qEdge()
        {
            // empty constructor
        }

        public qEdge(qNode _startNode, qNode _endNode)
        {
            StartNode = _startNode;
            EndNode = _endNode;
            Length = CalculateLength(_startNode, _endNode);
            EdgeLine = VisualizeLine(_startNode, _endNode);
            Level = 0; // to do: check if needed
        }

        // Methods:
        public double CalculateLength(qNode _startNode, qNode _endNode)
        {
            return _startNode.Coordinate.DistanceTo(_endNode.Coordinate);
        }
        public Line VisualizeLine(qNode _startNode, qNode _endNode)
        {
            return new Line(_startNode.Coordinate, _endNode.Coordinate);
        }
        public bool IsFrontEdge()
        {
            qEdge edge = this;
            // summary: check if an edge is a front edge
            bool check = false;
            List<qElement> connectedElements = this.GetConnectedElements();
            if (connectedElements.Count == 1)
            {
                if (!edge.Element1.IsQuad)
                {
                    check = true;
                }
            }
            else if (connectedElements.Count != 1)
            {
                if (!edge.Element1.IsQuad & edge.Element2.IsQuad) { check = true; }
                else if (edge.Element1.IsQuad & !edge.Element2.IsQuad) { check = true; }
            }

            return check;
        }

        public List<qElement> GetConnectedElements()
        {
            qEdge edge = this;
            // summary: get conneccted elements to an edge. Assume edge has updated elements element 1 and/or element 2.
            List<qElement> connectedElements = new List<qElement>();
            connectedElements.Add(edge.Element1);
            if (edge.Element2 != null)
            {
                connectedElements.Add(edge.Element2);
            }
            return connectedElements;
        }

        public double CalculateAngleOfNeighborFrontEdges(string nodeToCalculate)
        {
            // summary: calculate the angle between front edges with a shared point, i.e. neighbor front edges. NodeToCalculate 0 and 1 refer to left and right respectively.

            qEdge edgeNeighbor = new qEdge();
            double angle = 0;

            // Get neighbor nodes
            if (nodeToCalculate == "left") { edgeNeighbor = this.LeftFrontNeighbor; }
            else { edgeNeighbor = this.RightFrontNeighbor; }

            // Create vectors from shared node
            var vectors = this.CalculateVectorsFromSharedNode(edgeNeighbor);
            Vector3d vec1 = vectors.Item1;
            Vector3d vec2 = vectors.Item2;

            if (nodeToCalculate == "left")
            {
                angle = Vector3d.VectorAngle(vec1, vec2, Vector3d.ZAxis); // to do: make normal more general
            }
            else
            {
                angle = Vector3d.VectorAngle(vec2, vec1, Vector3d.ZAxis); // to do: make normal more general
            }

            /* to do: CHECK IF NEEDED
            qNode sharedNode = vectors.Item3;
            Vector3d V_k = vec1.Length * vec2 + vec2.Length * vec1; // angle bisector
            if (V_k == Vector3d.Zero) // create V_k if zero vector
            {
                if (nodeToCalculate == 0) { V_k = vec1; }
                else { V_k = vec2; }

                Vector3d rotationAxis = Vector3d.CrossProduct(vec1, vec2);
                V_k.Rotate(0.5 * Math.PI, rotationAxis);
            }
            V_k = V_k / V_k.Length; // normalize

            // check with domain
            Point3d endPointV_k = Point3d.Add(sharedNode.Coordinate, V_k); // endpoint of V_k from sharedNode
            Point3d endPointV_kNegative = Point3d.Add(sharedNode.Coordinate, -V_k); // endpoint of -V_k from sharedNode

            // distance to domain
            Point3d edgeFaceCenter = GetElementCenter(GetFrontElement(edge));
            Point3d edgeNeighborFaceCenter = GetElementCenter(GetFrontElement(edgeNeighbor));
            Point3d midFaceCenter = (edgeFaceCenter + edgeNeighborFaceCenter) / 2; // mid face center

            double distanceEndPoint = midFaceCenter.DistanceTo(endPointV_k); // todo: ok for all cases?
            double distanceEndPointNegative = midFaceCenter.DistanceTo(endPointV_kNegative); // todo: ok for all cases?
            double alpha = Vector3d.VectorAngle(vec1, vec2);

            if (distanceEndPoint < distanceEndPointNegative)
            {
                // V_k is inside domain
                angle = alpha;
            }
            else
            {
                // V_k is outside domain
                angle = 2 * Math.PI - alpha;
            }*/
            return angle;
        }

        public Tuple<Vector3d, Vector3d> CalculateVectorsFromSharedNode(qEdge neighborEdge)
        {
            // summary: calculate vectors for two edges with a shared node with direction away from the shared node.
            qNode sharedNode = this.GetSharedNode(neighborEdge);
            qNode node1 = this.GetOppositeNode(sharedNode);
            qNode node2 = neighborEdge.GetOppositeNode(sharedNode);

            Vector3d vec1 = node1.Coordinate - sharedNode.Coordinate;
            Vector3d vec2 = node2.Coordinate - sharedNode.Coordinate;

            return Tuple.Create(vec1, vec2);
        }

        public qNode GetSharedNode(qEdge neighborEdge)
        {
            // summary: get shared node of edges
            qNode sharedNode = new qNode();

            if (this.StartNode == neighborEdge.StartNode)
            {
                sharedNode = this.StartNode;
            }
            else if (this.StartNode == neighborEdge.EndNode)
            {
                sharedNode = this.StartNode;
            }
            else if (this.EndNode == neighborEdge.StartNode)
            {
                sharedNode = this.EndNode;
            }
            else if (this.EndNode == neighborEdge.EndNode)
            {
                sharedNode = this.EndNode;
            }
            return sharedNode;
        }

        public qNode GetOppositeNode(qNode node)
        {
            qNode oppositeNode = new qNode();
            if (node == this.StartNode) { oppositeNode = this.EndNode; }
            else if (node == this.EndNode) { oppositeNode = this.StartNode; }
            return oppositeNode;
        }

        public qElement GetFrontElement()
        {
            // summary: get the connected triangle element of a front edge
            qElement triangleElement = new qElement();
            if (this.GetConnectedElements().Count == 1) { triangleElement = this.Element1; }
            else if (!this.Element1.IsQuad) { triangleElement = this.Element1; }
            else { triangleElement = this.Element2; }

            return triangleElement;
        }

        public Tuple<qEdge, qEdge> OrientateNeigborEdges(qEdge neigborEdgeToStartNode, qEdge neigborEdgeToEndNode)
        {
            // summary: oridentate the neigbors of an edge. Return the left edge and the right edge, respectively.
            qEdge edge = this;
            Point3d midPointEdg = 0.5 * (edge.StartNode.Coordinate + edge.EndNode.Coordinate); // mid point of edge

            Point3d centerPoint = this.GetFrontElement().GetElementCenter();

            Vector3d centerToMidVector = midPointEdg - centerPoint;

            Vector3d centerToEndNodeVector = edge.EndNode.Coordinate - centerPoint;

            Vector3d centerToStartNodeVector = edge.StartNode.Coordinate - centerPoint;

            double startAngle = Vector3d.VectorAngle(centerToMidVector, centerToStartNodeVector, Vector3d.ZAxis); // todo: make normal more general

            double endAngle = Vector3d.VectorAngle(centerToMidVector, centerToEndNodeVector, Vector3d.ZAxis); // todo: make normal more general

            qEdge leftEdge = new qEdge();
            qEdge rightEdge = new qEdge();

            if (endAngle < startAngle)
            {
                leftEdge = neigborEdgeToStartNode;
                rightEdge = neigborEdgeToEndNode;
            }
            else
            {
                leftEdge = neigborEdgeToEndNode;
                rightEdge = neigborEdgeToStartNode;
            }
            return Tuple.Create(leftEdge, rightEdge);
        }

    }
}
