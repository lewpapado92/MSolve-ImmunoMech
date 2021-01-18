using System;
using System.Linq;
using System.Text;
using System.Collections.Generic;
using System.Globalization;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.FEM.Elements.BoundaryConditionElements;

namespace ISAAR.MSolve.FEM.Readers
{
    public class ComsolMeshReader
    {
        public Model Model { get; private set; }
        public IList<IList<Node>> nodeBoundaries;
        public IList<IList<Element>> elementBoundaries;
        public IList<IList<Element>> elementDomains;
        public double diffusionCeoff;
        enum Attributes
        {
            sdim = 1001,
            Meshpointcoordinates = 1002,
            numberofmeshpoints = 1003,
            numberofelementtypes = 1004,
            vtx = 1006,
            edg = 1007,
            tri = 1008,
            tet = 1009,
            hex = 1010,
            quad = 1011,
        }

        public string Filename { get; private set; }

        public ComsolMeshReader(string filename)
        {
            Filename = filename;
        }

        int NumberOfNodes;
        int NumberOfTriElements;
        int NumberOfTetElements;
        int NumberOfHexElements;
        int NumberOfQuadElements;

        public Model CreateModelFromFile()
        {
            double[] U = {1.0,1,1};
            double k = 1.0;
            double L = .0;
            diffusionCeoff = k;
            var elementFactory3D = new ConvectionDiffusionElement3DFactory(new ConvectionDiffusionMaterial(k, U, L));
            var boundaryFactory3D = new SurfaceBoundaryFactory3D(0, new ConvectionDiffusionMaterial(k, U, L));
            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(0);
            // Material

            char[] delimeters = { ' ', '=', '\t' };
            Attributes? name = null;

            String[] text = System.IO.File.ReadAllLines(Filename);
            elementBoundaries = new List<IList<Element>>();
            elementDomains = new List<IList<Element>>();
            nodeBoundaries = new List<IList<Node>>();
            for (int i = 0; i < 10; i++)
            {
                elementBoundaries.Add(new List<Element>());
                nodeBoundaries.Add(new List<Node>());
            }
            for (int i = 0; i < 2; i++)
            {
                elementDomains.Add(new List<Element>());
            }

            for (int i = 0; i < text.Length; i++)
            {
                String[] line = text[i].Split(delimeters, StringSplitOptions.RemoveEmptyEntries);
                StringBuilder comparisonString = new StringBuilder(null, 50);
                if (line.Length == 0)
                {
                    continue;
                }
                if (line[0] == "3" & line.Length > 1)
                {
                    try
                    {
                        name = (Attributes)Enum.Parse(typeof(Attributes), line[1]);
                    }
                    catch (Exception exception)
                    {
                        name = null;
                    }
                }
                else
                {
                    for (int linePosition = 0; linePosition < line.Length; linePosition++)
                    {
                        if (line[linePosition] == "#")
                        {
                            for (int ij = linePosition + 1; ij < line.Length; ij++)
                            {
                                comparisonString.Append(line[ij]);
                            }
                            try
                            {
                                name = (Attributes)Enum.Parse(typeof(Attributes), comparisonString.ToString());
                            }
                            catch (Exception exception)
                            {
                                name = null;
                            }
                        }
                    }
                }
                switch (name)
                {
                    case Attributes.sdim:
                        int NumberOfDimensions = Int32.Parse(line[0]);
                        break;
                    case Attributes.numberofmeshpoints:
                        NumberOfNodes = Int32.Parse(line[0]);
                        break;
                    case Attributes.Meshpointcoordinates:
                        IList<Node> nodelist = new List<Node> { null };
                        for (int j = 0; j < NumberOfNodes; j++)
                        {
                            i++;
                            line = text[i].Split(delimeters);
                            int nodeGlobalID = j;
                            double x = Double.Parse(line[0], CultureInfo.InvariantCulture);
                            double y = Double.Parse(line[1], CultureInfo.InvariantCulture);
                            double z = Double.Parse(line[2], CultureInfo.InvariantCulture);
                            Node node = new Node (nodeGlobalID, x, y, z);
                            model.NodesDictionary.Add(nodeGlobalID, node);
                            nodelist.Add(node);
                        }

                        break;
                    case Attributes.tri:
                        do
                        {
                            i++;
                            line = text[i].Split(delimeters);
                        } while (line[0] == "");
                        i++;
                        line = text[i].Split(delimeters);
                        NumberOfTriElements = Int32.Parse(line[0]);
                        i++;
                        IList<Node> nodesCollection = new List<Node>();
                        for (int TriID = 0; TriID < NumberOfTriElements; TriID++)
                        {
                            i++;
                            line = text[i].Split(delimeters);

                            for (int j = 0; j < (line.Length - 1); j++)
                            {
                                nodesCollection.Add(model.NodesDictionary[Int32.Parse(line[j])]);
                            }

                            IReadOnlyList<Node> nodes = new List<Node>
                            {
                                model.NodesDictionary[Int32.Parse(line[2])],
                                model.NodesDictionary[Int32.Parse(line[1])],
                                model.NodesDictionary[Int32.Parse(line[0])]
                            };

                            var Tri3 = boundaryFactory3D.CreateElement(CellType.Tri3, nodes);
                            var element = new Element();
                            element.ID = TriID;
                            element.ElementType = Tri3;
                            model.SubdomainsDictionary[0].Elements.Add(element);
                            model.ElementsDictionary.Add(TriID, element);
                            double r3 = 0;
                            foreach (Node node in nodes)
                            {
                                element.AddNode(node);
                                r3 += Math.Sqrt(Math.Pow(node.X, 2) + Math.Pow(node.Y, 2) + Math.Pow(node.Z, 2));
                            }
                            if (element.Nodes[0].X == 0 && r3 < 3 * 5e-4)
                            {
                                elementBoundaries[0].Add(model.ElementsDictionary[TriID]);
                                //model.NodesDictionary[element.Nodes[0].ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 10 });
                            }
                            else if (element.Nodes[0].Y == 0 && r3 < 3 * 5e-4)
                            {
                                elementBoundaries[1].Add(model.ElementsDictionary[TriID]);
                            }
                            else if (element.Nodes[0].Z == 0 && r3 < 3 * 5e-4)
                            {
                                elementBoundaries[2].Add(model.ElementsDictionary[TriID]);
                            }
                            else if (element.Nodes[0].X == 0 && r3 > 3 * 5e-4)
                            {
                                elementBoundaries[3].Add(model.ElementsDictionary[TriID]);
                            }
                            else if (element.Nodes[0].Y == 0 && r3 > 3 * 5e-4)
                            {
                                elementBoundaries[4].Add(model.ElementsDictionary[TriID]);
                            }
                            else if (r3 == 3 * 5e-4)
                            {
                                elementBoundaries[5].Add(model.ElementsDictionary[TriID]);
                            }
                            else if (element.Nodes[0].Z == 0.1 && element.Nodes[1].Z == 0.1 && element.Nodes[2].Z == 0.1)
                            {
                                elementBoundaries[6].Add(model.ElementsDictionary[TriID]);
                            }
                            else if (element.Nodes[0].Z == 0 && r3 > 3 * 5e-4)
                            {
                                elementBoundaries[7].Add(model.ElementsDictionary[TriID]);
                            }
                            else if (element.Nodes[0].Y == 0.1 && element.Nodes[1].Y == 0.1 && element.Nodes[2].Y == 0.1)
                            {
                                elementBoundaries[8].Add(model.ElementsDictionary[TriID]);
                            }
                            else if (element.Nodes[0].X == 0.1 && element.Nodes[1].X == 0.1 && element.Nodes[2].X == 0.1)
                            {
                                elementBoundaries[9].Add(model.ElementsDictionary[TriID]);
                            }
                        }
                        nodesCollection = nodesCollection.Distinct().ToList();
                        foreach (Node node in nodesCollection)
                        {
                            double r = Math.Sqrt(Math.Pow(node.X, 2) + Math.Pow(node.Y, 2) + Math.Pow(node.Z, 2));
                            if (node.X == 0 && node.Y <= 5e-4 && node.Z <= 5e-4)
                            {
                                nodeBoundaries[0].Add(model.NodesDictionary[node.ID]);
                                //model.NodesDictionary[node.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 10 });
                            }
                            else if (node.X <= 5e-4 && node.Y == 0 && node.Z <= 5e-4)
                            {
                                nodeBoundaries[1].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.X <= 5e-4 && node.Y <= 5e-4 && node.Z == 0)
                            {
                                nodeBoundaries[2].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.X == 0 && r > 5e-4)
                            {
                                nodeBoundaries[3].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.Y == 0 && r > 5e-4)
                            {
                                nodeBoundaries[4].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (r == 5e-4)
                            {
                                nodeBoundaries[5].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.Z == 0.1)
                            {
                                nodeBoundaries[6].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.Z == 0 && r > 5e-4)
                            {
                                nodeBoundaries[7].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.Y == 0.1)
                            {
                                nodeBoundaries[8].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.X == 0.1)
                            {
                                nodeBoundaries[9].Add(model.NodesDictionary[node.ID]);
                            }
                            //else
                            //{
                            //    model.Loads.Add(new Load() { Node = model.NodesDictionary[node.ID], DOF = ThermalDof.Temperature, Amount = 10 });
                            //}
                        }
                        break;
                    case Attributes.quad:
                        do
                        {
                            i++;
                            line = text[i].Split(delimeters);
                        } while (line[0] == "");
                        i++;
                        line = text[i].Split(delimeters);
                        NumberOfTriElements = Int32.Parse(line[0]);
                        i++;
                        IList<Node> quadNodesCollection = new List<Node>();
                        for (int QuadID = 0; QuadID < NumberOfTriElements; QuadID++)
                        {
                            i++;
                            line = text[i].Split(delimeters);

                            for (int j = 0; j < (line.Length - 1); j++)
                            {
                                quadNodesCollection.Add(model.NodesDictionary[Int32.Parse(line[j])]);
                            }

                            IReadOnlyList<Node> nodes = new List<Node>
                            {
                                model.NodesDictionary[Int32.Parse(line[3])],
                                model.NodesDictionary[Int32.Parse(line[2])],
                                model.NodesDictionary[Int32.Parse(line[1])],
                                model.NodesDictionary[Int32.Parse(line[0])]
                            };

                            var Quad = boundaryFactory3D.CreateElement(CellType.Quad4, nodes);
                            var element = new Element();
                            element.ID = QuadID;
                            element.ElementType = Quad;
                            model.SubdomainsDictionary[0].Elements.Add(element);
                            model.ElementsDictionary.Add(QuadID, element);
                            //double r3 = 0;
                            //foreach (Node node in nodes)
                            //{
                            //    element.AddNode(node);
                            //    r3 += Math.Sqrt(Math.Pow(node.X, 2) + Math.Pow(node.Y, 2) + Math.Pow(node.Z, 2));
                            //}
                            //if (element.Nodes[0].X == 0 && r3 < 3 * 5e-4)
                            //{
                            //    elementBoundaries[0].Add(model.ElementsDictionary[QuadID]);
                            //    //model.NodesDictionary[element.Nodes[0].ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 10 });
                            //}
                            //else if (element.Nodes[0].Y == 0 && r3 < 3 * 5e-4)
                            //{
                            //    elementBoundaries[1].Add(model.ElementsDictionary[QuadID]);
                            //}
                            //else if (element.Nodes[0].Z == 0 && r3 < 3 * 5e-4)
                            //{
                            //    elementBoundaries[2].Add(model.ElementsDictionary[QuadID]);
                            //}
                            //else if (element.Nodes[0].X == 0 && r3 > 3 * 5e-4)
                            //{
                            //    elementBoundaries[3].Add(model.ElementsDictionary[QuadID]);
                            //}
                            //else if (element.Nodes[0].Y == 0 && r3 > 3 * 5e-4)
                            //{
                            //    elementBoundaries[4].Add(model.ElementsDictionary[QuadID]);
                            //}
                            //else if (r3 == 3 * 5e-4)
                            //{
                            //    elementBoundaries[5].Add(model.ElementsDictionary[QuadID]);
                            //}
                            //else if (element.Nodes[0].Z == 0.1 && element.Nodes[1].Z == 0.1 && element.Nodes[2].Z == 0.1)
                            //{
                            //    elementBoundaries[6].Add(model.ElementsDictionary[QuadID]);
                            //}
                            //else if (element.Nodes[0].Z == 0 && r3 > 3 * 5e-4)
                            //{
                            //    elementBoundaries[7].Add(model.ElementsDictionary[QuadID]);
                            //}
                            //else if (element.Nodes[0].Y == 0.1 && element.Nodes[1].Y == 0.1 && element.Nodes[2].Y == 0.1)
                            //{
                            //    elementBoundaries[8].Add(model.ElementsDictionary[QuadID]);
                            //}
                            //else if (element.Nodes[0].X == 0.1 && element.Nodes[1].X == 0.1 && element.Nodes[2].X == 0.1)
                            //{
                            //    elementBoundaries[9].Add(model.ElementsDictionary[QuadID]);
                            //}
                        }
                        quadNodesCollection = quadNodesCollection.Distinct().ToList();
                        foreach (Node node in quadNodesCollection)
                        {
                            double r = Math.Sqrt(Math.Pow(node.X, 2) + Math.Pow(node.Y, 2) + Math.Pow(node.Z, 2));
                            if (node.X == 0 && node.Y <= 5e-4 && node.Z <= 5e-4)
                            {
                                nodeBoundaries[0].Add(model.NodesDictionary[node.ID]);
                                //model.NodesDictionary[node.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 10 });
                            }
                            else if (node.X <= 5e-4 && node.Y == 0 && node.Z <= 5e-4)
                            {
                                nodeBoundaries[1].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.X <= 5e-4 && node.Y <= 5e-4 && node.Z == 0)
                            {
                                nodeBoundaries[2].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.X == 0 && r > 5e-4)
                            {
                                nodeBoundaries[3].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.Y == 0 && r > 5e-4)
                            {
                                nodeBoundaries[4].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (r == 5e-4)
                            {
                                nodeBoundaries[5].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.Z == 0.1)
                            {
                                nodeBoundaries[6].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.Z == 0 && r > 5e-4)
                            {
                                nodeBoundaries[7].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.Y == 0.1)
                            {
                                nodeBoundaries[8].Add(model.NodesDictionary[node.ID]);
                            }
                            else if (node.X == 0.1)
                            {
                                nodeBoundaries[9].Add(model.NodesDictionary[node.ID]);
                            }
                            //else
                            //{
                            //    model.Loads.Add(new Load() { Node = model.NodesDictionary[node.ID], DOF = ThermalDof.Temperature, Amount = 10 });
                            //}
                        }
                        break;
                    case Attributes.tet:
                        do
                        {
                            i++;
                            line = text[i].Split(delimeters);
                        } while (line[0] == "");
                        i++;
                        line = text[i].Split(delimeters);
                        NumberOfTetElements = Int32.Parse(line[0]);
                        i++;
                        for (int TetID = NumberOfTriElements + 1; TetID < NumberOfTriElements + NumberOfTetElements; TetID++)
                        {
                            i++;
                            line = text[i].Split(delimeters);

                            IReadOnlyList<Node> nodes = new List<Node>
                            {
                                model.NodesDictionary[Int32.Parse(line[3])],
                                model.NodesDictionary[Int32.Parse(line[2])],
                                model.NodesDictionary[Int32.Parse(line[1])],
                                model.NodesDictionary[Int32.Parse(line[0])]
                            };
                            var Tet4 = elementFactory3D.CreateElement(CellType.Tet4, nodes);
                            var element = new Element();
                            element.ID = TetID;
                            element.ElementType = Tet4;
                            double r4 = 0;
                            foreach (Node node in nodes)
                            {
                                element.AddNode(node);
                                r4 += Math.Sqrt(Math.Pow(node.X, 2) + Math.Pow(node.Y, 2) + Math.Pow(node.Z, 2));
                            }
                            model.SubdomainsDictionary[0].Elements.Add(element);
                            model.ElementsDictionary.Add(TetID, element);
                            if (r4 < 4*5e-4)
                                elementDomains[0].Add(element);
                            else
                                elementDomains[1].Add(element);
                        }
                        name = null;
                        break;
                }
            }
            return model;
        }

    }
}
