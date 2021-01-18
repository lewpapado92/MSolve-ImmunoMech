using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Readers;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;
using System.IO;
using ISAAR.MSolve.FEM.Elements.BoundaryConditionElements;
using ISSAR.MSolve.Discretization.Loads;

namespace ISAAR.MSolve.Tests.FEM
{
    public class TumorCellsConcentrationHexa27
    {
        private const int subdomainID = 0;

        [Fact]
        private static void RunTest()
        {
            Model model = CreateModel();
            IVectorView solution = SolveModel(model);
            Assert.True(CompareResults(solution));
        }

        private static bool CompareResults(IVectorView solution)
        {
            var comparer = new ValueComparer(1E-5);

            //                                                   dofs:   1,   2,   4,   5,   7,   8
            var expectedSolution = Vector.CreateFromArray(new double[] { 150, 200, 150, 200, 150, 200 });
            int numFreeDofs = 6;
            if (solution.Length != 6) return false;
            for (int i = 0; i < numFreeDofs; ++i)
            {
                if (!comparer.AreEqual(expectedSolution[i], solution[i])) return false;
            }
            return true;
        }

        private static Model CreateModel()
        {
            IList<IList<Node>> nodeBoundaries;
            IList<IList<Element>> elementBoundaries;
            double density = 1.0;
            double c = 1.0;
            double k = 1.0;
            double[] U = {1.0,1,1};
            double h = .0;
            var elementFactory3D = new ConvectionDiffusionElement3DFactory(new ConvectionDiffusionMaterial(k, U, 0));
            var boundaryFactory3D = new SurfaceBoundaryFactory3D(0, new ConvectionDiffusionMaterial(k, U, 0));
            var model = new Model();

            model.SubdomainsDictionary[0] = new Subdomain(0);
            elementBoundaries = new List<IList<Element>>();
            nodeBoundaries = new List<IList<Node>>();
            for (int i = 0; i < 10; i++)
            {
                elementBoundaries.Add(new List<Element>());
                nodeBoundaries.Add(new List<Node>());
            }

            double[,] nodeData = new double[,]
            {
                {0,0,0}, {0,.1,0}, {0,.1,.1}, {0,0,.1}, {.1,0,0}, {.1,.1,0}, {.1,.1,.1}, {.1,0,.1}, 
                {0,0.05,0}, {0,0,0.05}, {.05,0,0}, {0,.1,0.05}, {.05,.1,0}, {0,0.05,.1}, {.05,.1,.1},
                {.05,0,.1}, {.1,0.05,0}, {.1,0,0.05}, {.1,.1,0.05}, {.1,0.05,.1}, {0,0.05,0.05},
                {.05,0.05,0}, {.05,0,0.05}, {.05,.1,0.05}, {.05,0.05,.1}, {.1,0.05,0.05},
                {.05,0.05,0.05},    
                      
            };

            for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nNode, new Node(id: nNode, x: nodeData[nNode, 0], y: nodeData[nNode, 1], z: nodeData[nNode, 2]));

            }
            IReadOnlyList<Node> nodes = new List<Node>
            {
                model.NodesDictionary[0],model.NodesDictionary[1],model.NodesDictionary[2],
                model.NodesDictionary[3],model.NodesDictionary[4],model.NodesDictionary[5],
                model.NodesDictionary[6],model.NodesDictionary[7],model.NodesDictionary[8],
                model.NodesDictionary[9],model.NodesDictionary[10],model.NodesDictionary[11],
                model.NodesDictionary[12],model.NodesDictionary[13],model.NodesDictionary[14],
                model.NodesDictionary[15],model.NodesDictionary[16],model.NodesDictionary[17],
                model.NodesDictionary[18],model.NodesDictionary[19],model.NodesDictionary[20],
                model.NodesDictionary[21],model.NodesDictionary[22],model.NodesDictionary[23],
                model.NodesDictionary[24],model.NodesDictionary[25],model.NodesDictionary[26],
            };
            IList<IReadOnlyList<Node>> face = new List<IReadOnlyList<Node>>();
            //face = new List<IReadOnlyList<Node>>();
            for (int i = 0; i < 6; i++)
            {
                face.Add(new List<Node>());
            }

            face[0] = new List<Node>
            {
                model.NodesDictionary[0],
                model.NodesDictionary[1],
                model.NodesDictionary[2],
                model.NodesDictionary[3],
                model.NodesDictionary[4],
                model.NodesDictionary[5],
                model.NodesDictionary[6],
                model.NodesDictionary[7],
                model.NodesDictionary[8],
            };
            face[1] = new List<Node>
            {
                model.NodesDictionary[0],
                model.NodesDictionary[1],
                model.NodesDictionary[2],
                model.NodesDictionary[9],
                model.NodesDictionary[10],
                model.NodesDictionary[11],
                model.NodesDictionary[18],
                model.NodesDictionary[19],
                model.NodesDictionary[20],
            };
            face[2] = new List<Node>
            {
                model.NodesDictionary[18],
                model.NodesDictionary[19],
                model.NodesDictionary[20],
                model.NodesDictionary[21],
                model.NodesDictionary[22],
                model.NodesDictionary[23],
                model.NodesDictionary[24],
                model.NodesDictionary[25],
                model.NodesDictionary[26],
            };
            face[3] = new List<Node>
            {
                model.NodesDictionary[6],
                model.NodesDictionary[7],
                model.NodesDictionary[8],
                model.NodesDictionary[15],
                model.NodesDictionary[16],
                model.NodesDictionary[17],
                model.NodesDictionary[24],
                model.NodesDictionary[25],
                model.NodesDictionary[26],
            };
            face[4] = new List<Node>
            {
                model.NodesDictionary[0],
                model.NodesDictionary[3],
                model.NodesDictionary[6],
                model.NodesDictionary[9],
                model.NodesDictionary[12],
                model.NodesDictionary[15],
                model.NodesDictionary[18],
                model.NodesDictionary[21],
                model.NodesDictionary[24],
            };
            face[5] = new List<Node>
            {
                model.NodesDictionary[2],
                model.NodesDictionary[5],
                model.NodesDictionary[8],
                model.NodesDictionary[11],
                model.NodesDictionary[14],
                model.NodesDictionary[17],
                model.NodesDictionary[20],
                model.NodesDictionary[23],
                model.NodesDictionary[26],
            };

            var Hexa27 = elementFactory3D.CreateElement(CellType.Hexa27, nodes);
            var element = new Element();
            element.ID = 0;
            element.ElementType = Hexa27;
            foreach (Node node in nodes)
            {
                element.AddNode(node);
            }
            model.SubdomainsDictionary[0].Elements.Add(element);
            model.ElementsDictionary.Add(0, element);
            for (int i = 1; i < 7; i++)
            {
                var Quad = boundaryFactory3D.CreateElement(CellType.Quad9, face[i-1]);
                var faceElement = new Element();
                faceElement.ID = i;
                faceElement.ElementType = Quad;
                foreach (Node node in face[i-1])
                {
                    faceElement.AddNode(node);
                }
                model.SubdomainsDictionary[0].Elements.Add(faceElement);
                model.ElementsDictionary.Add(i, faceElement);
            }
            //constraints
            foreach (int i in new int[] { 2, 6, 7 })
            {
                model.NodesDictionary[i].Constraints.Add(new Constraint()
                {
                    Amount = 0,
                    DOF = ThermalDof.Temperature
                });
            }

            // fortish korufhs
            Load load1;
            foreach (int i in new int[] {0, 1, 3, 4, 5})
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[i],
                    DOF = ThermalDof.Temperature,
                    Amount = 100
                };
                model.Loads.Add(load1);
            }


            return model;
        }

        private static IVectorView SolveModel(Model model)
        {
            double[] temp0 = new double[] { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            Vector initialTemp = Vector.CreateFromArray(temp0);
            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);
            var provider = new ProblemConvectionDiffusion(model, solver);

            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new ConvectionDiffusionExplicitDynamicAnalyzer(model, solver, provider, childAnalyzer, 0.05, 5, initialTemp);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.LinearSystems[subdomainID].Solution;
        }
    }
}

