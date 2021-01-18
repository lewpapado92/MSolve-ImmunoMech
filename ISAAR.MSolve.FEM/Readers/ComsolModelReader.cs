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
using ISSAR.MSolve.Discretization.Loads;

namespace ISAAR.MSolve.FEM.Readers
{
    public class ComsolModelReader
    {
        public Model Model { get; private set; }
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
        }

        public string Filename { get; private set; }

        public ComsolModelReader(string filename)
        {
            Filename = filename;
        }

        int NumberOfNodes;
        int NumberOfTriElements;
        int NumberOfElements;

        public Model CreateModelFromFile()
        {
            double density = 1.0;
            double k = 1.0;
            double c = 1.0;
            var elementFactory = new ThermalElement3DFactory(new ThermalMaterial(1,density, c, k));
            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(0);
            // Material

            char[] delimeters = { ' ', '=', '\t' };
            Attributes? name = null;

            String[] text = System.IO.File.ReadAllLines(Filename);

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
                            double x = Double.Parse(line[2], CultureInfo.InvariantCulture);
                            double y = Double.Parse(line[1], CultureInfo.InvariantCulture);
                            double z = Double.Parse(line[0], CultureInfo.InvariantCulture);
                            Node node = new Node(nodeGlobalID, x, y, z);
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
                        }
                        nodesCollection = nodesCollection.Distinct().ToList();
                        foreach (Node node in nodesCollection)
                        {
                            if (node.X <= 5e-4 && node.Y <= 5e-4 && node.Z <= 5e-4)
                            {
                                model.NodesDictionary[node.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 10 });
                            }
                            else
                            {
                                model.Loads.Add(new Load() { Node = model.NodesDictionary[node.ID], DOF = ThermalDof.Temperature, Amount = 10 });
                            }
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
                        NumberOfElements = Int32.Parse(line[0]);
                        i++;
                        for (int TetID = 0; TetID < NumberOfElements; TetID++)
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
                            var Tet4 = elementFactory.CreateElement(CellType.Tet4, nodes);
                            var element = new Element();
                            element.ID = TetID;
                            element.ElementType = Tet4;
                            foreach (Node node in nodes)
                            {
                                element.AddNode(node);
                            }
                            model.SubdomainsDictionary[0].Elements.Add(element);
                            model.ElementsDictionary.Add(TetID, element);
                        }
                        break;
                }
            }
            return model;
        }

    }
}
