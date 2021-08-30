using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

    namespace SolidFEM.Classes
    {
        class SmartMesh
        {
            public List<Element> Elements { get; set; } //list of elements
            public List<Node> Nodes { get; set; } //list of nodes
            public Mesh Mesh { get; set; } //mesh
            public int nu { get; set; } //number of nodes in x-dir
            public int nv { get; set; } //number of nodes in y-dir
            public int nw { get; set; } //number of nodes in z-dir
            public string Type { get; set; } // to do: inplementer
            public Geometry Geometry { get; set; }
            public List<List<List<Point3d>>> GridInformation { get; set; }

            // Constructors
            public SmartMesh()
            {
                //Empty constructor
            }

            public SmartMesh(int _nu, int _nv, List<Node> _nodes, List<Element> _elements, Mesh _mesh) // for surface mesh
            {
                nu = _nu;
                nv = _nv;
                nw = 1;
                Nodes = _nodes;
                Elements = _elements;
                Mesh = _mesh;
                Type = "Surface";
            }
            public SmartMesh(int _nu, int _nv, int _nw, List<Node> _nodes, List<Element> _elements, Mesh _mesh) // for solid mesh
            {
                nu = _nu;
                nv = _nv;
                nw = _nw;
                Nodes = _nodes;
                Elements = _elements;
                Mesh = _mesh;
                Type = "Solid";
            }
            public SmartMesh(List<Node> _nodes, List<Element> _elements, String _type) // for unstructured surface mesh
            {
                Nodes = _nodes;
                Elements = _elements;
                Type = _type;
                CreateMesh();
            }

            // Methods
            public void CreateNodes(List<Point3d> meshPoints, int u, int v, int w)
            {
                List<Node> nodes = new List<Node>();
                int id = 0;
                int nu = u + 1; // number nodes in u dir
                int nv = v + 1; // number nodes in v dir 
                int nw = w + 1; // number nodes in w dir 

                for (int i = 0; i < nw; i++)
                {
                    int row = 0;
                    int column = 0;
                    for (int j = 0; j < (nu * nv); j++)
                    {
                        bool BC_U = false;
                        bool BC_V = false;
                        bool BC_W = false;

                        if (column == 0 | column == nu - 1) { BC_U = true; } // assign BCU
                        if (row == 0 | row == nv - 1) { BC_V = true; } // assign BCV
                        if (i == 0 | i == nw - 1) { BC_W = true; } // assign BCW

                        Node node = new Node(id, meshPoints[j + i * (nu * nv)], BC_U, BC_V, BC_W);
                        id++;

                        column++;
                        if (column == nu)
                        {
                            row++;
                            column = 0;
                        }
                        nodes.Add(node);
                    }
                }
                this.Nodes = nodes;
            }
            public void CreateQuadElements()
            {
                // Create quad elements of a structured SmartMesh given nodes are assigned
                List<Node> nodes = this.Nodes;
                int nu = this.nu;
                int nv = this.nv;

                List<Element> elements = new List<Element>();
                int uSequence = 0;
                int counter = 0;

                for (int i = 0; i < (nu - 1) * (nv - 1); i++) // loop elements
                {
                    List<Node> elementNodes = new List<Node>();
                    List<int> connectivity = new List<int>();
                    connectivity.Add(counter);
                    connectivity.Add(counter + 1);
                    connectivity.Add(counter + nu + 1);
                    connectivity.Add(counter + nu);

                    foreach (int id in connectivity)
                    {
                        elementNodes.Add(nodes[id]);
                    };

                    Element element = new Element(i, elementNodes, connectivity);
                    elements.Add(element); // add element to list of elements

                    counter++;
                    uSequence++;
                    if (uSequence == (nu - 1)) // check if done with a v sequence
                    {
                        counter++;
                        uSequence = 0; // new v sequence
                    }
                }
                this.Elements = elements;
            }
            public void CreateHexElements()
            {
                // Create hex elements of a structured SmartMesh given nodes are assigned

                int nu = this.nu;
                int nv = this.nv;
                int nw = this.nw;
                List<Node> nodes = this.Nodes;
                List<Element> elements = new List<Element>();

                int elemId = 0;

                for (int i = 0; i < nw - 1; i++)  // loop levels
                {
                    int sequence = 0;
                    int counter = (nu * nv) * i;

                    for (int j = 0; j < (nu * nv) - nu - 1; j++) // loop elements in a level
                    {
                        List<Node> elementNodes = new List<Node>();
                        List<int> connectivity = new List<int>();

                        if (sequence < nu - 1)
                        {
                            connectivity.Add(counter);
                            connectivity.Add(counter + 1);
                            connectivity.Add(counter + nu + 1);
                            connectivity.Add(counter + nu);
                            connectivity.Add(counter + nu * nv);
                            connectivity.Add(counter + 1 + nu * nv);
                            connectivity.Add(counter + nu + 1 + nu * nv);
                            connectivity.Add(counter + nu + nu * nv);

                            foreach (int id in connectivity)
                            {
                                elementNodes.Add(nodes[id]);
                            }

                            Element element = new Element(elemId, elementNodes, connectivity);
                            elements.Add(element);

                            sequence++;
                            elemId++;
                            counter++;
                        }
                        else { sequence = 0; counter++; }
                    }
                }
                this.Elements = elements;
            }
            public void CreateMesh()
            {
                Mesh mesh = new Mesh();

                // Create mesh vertices from node coordinates        
                foreach (Node node in this.Nodes)
                {
                    mesh.Vertices.Add(node.Coordinate);
                }

                foreach (Element element in this.Elements)
                {
                    // Create mesh faces from element connectivity
                    if (element.Type == "Quad")
                    {

                        mesh.Faces.AddFace(element.Connectivity[0], element.Connectivity[1], element.Connectivity[2], element.Connectivity[3]);
                    }
                    else if (element.Type == "Hex")
                    {
                        int a = element.Connectivity[0];
                        int b = element.Connectivity[1];
                        int c = element.Connectivity[2];
                        int d = element.Connectivity[3];
                        int e = element.Connectivity[4];
                        int f = element.Connectivity[5];
                        int g = element.Connectivity[6];
                        int h = element.Connectivity[7];

                        mesh.Faces.AddFace(a, b, f, e);
                        mesh.Faces.AddFace(b, c, g, f);
                        mesh.Faces.AddFace(c, d, h, g);
                        mesh.Faces.AddFace(d, a, e, h);
                        mesh.Faces.AddFace(a, b, c, d);
                        mesh.Faces.AddFace(e, f, g, h);
                    }
                    else if (element.Type == "Triangle")
                    {
                        mesh.Faces.AddFace(element.Connectivity[0], element.Connectivity[1], element.Connectivity[2]);
                    }
                    else if (element.Type == "Tet")
                    {
                        int a = element.Connectivity[0];
                        int b = element.Connectivity[1];
                        int c = element.Connectivity[2];
                        int d = element.Connectivity[3];
                        int e = element.Connectivity[4];
                        int f = element.Connectivity[5];

                        mesh.Faces.AddFace(a, b, c);
                        mesh.Faces.AddFace(d, e, f);
                        mesh.Faces.AddFace(a, b, e, d);
                        mesh.Faces.AddFace(b, c, f, e);
                        mesh.Faces.AddFace(c, a, d, f);
                    }
                }


                mesh.Normals.ComputeNormals();
                mesh.Compact();
                mesh.FaceNormals.ComputeFaceNormals();
                this.Mesh = mesh;
            }
        }
    }


