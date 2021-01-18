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
using ISAAR.MSolve.FEM.Loading.SurfaceLoads;
using static ISAAR.MSolve.FEM.Loading.SurfaceLoads.WeakDirichlet;
using ISAAR.MSolve.FEM.Loading;

namespace ISAAR.MSolve.Tests.FEM
{
    public class TumorCellsConcentrationHexa8
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
            double[] U = {2,2,2};
            double k = 1.0;
            double L = .0;
            var elementFactory3D = new ConvectionDiffusionElement3DFactory(new ConvectionDiffusionMaterial(k, U, L));
            var boundaryFactory3D = new SurfaceBoundaryFactory3D(0, new ConvectionDiffusionMaterial(k, U, L));
            var model = new Model();

            model.SubdomainsDictionary[0] = new Subdomain(0);
            elementBoundaries = new List<IList<Element>>();
            nodeBoundaries = new List<IList<Node>>();
            for (int i = 0; i < 10; i++)
            {
                elementBoundaries.Add(new List<Element>());
                nodeBoundaries.Add(new List<Node>());
            }

            double[,] nodeData = new double[,] { {0,0,0},
            {0,1,0},
            {0,1,1},
            {0,0,1},
            {1,0,0},
            {1,1,0},
            {1,1,1},
            {1,0,1} };

            for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nNode, new Node(id: nNode, x: nodeData[nNode, 0], y: nodeData[nNode, 1], z: nodeData[nNode, 2]));

            }
            IReadOnlyList<Node> nodes = new List<Node>
            {
                model.NodesDictionary[0],
                model.NodesDictionary[1],
                model.NodesDictionary[2],
                model.NodesDictionary[3],
                model.NodesDictionary[4],
                model.NodesDictionary[5],
                model.NodesDictionary[6],
                model.NodesDictionary[7],
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
            };
            face[1] = new List<Node>
            {
                model.NodesDictionary[4],
                model.NodesDictionary[5],
                model.NodesDictionary[6],
                model.NodesDictionary[7],
            };
            face[2] = new List<Node>
            {
                model.NodesDictionary[0],
                model.NodesDictionary[1],
                model.NodesDictionary[5],
                model.NodesDictionary[4],
            };
            face[3] = new List<Node>
            {
                model.NodesDictionary[1],
                model.NodesDictionary[2],
                model.NodesDictionary[6],
                model.NodesDictionary[5],
            };
            face[4] = new List<Node>
            {
                model.NodesDictionary[2],
                model.NodesDictionary[6],
                model.NodesDictionary[7],
                model.NodesDictionary[3],
            };
            face[5] = new List<Node>
            {
                model.NodesDictionary[0],
                model.NodesDictionary[4],
                model.NodesDictionary[7],
                model.NodesDictionary[3],
            };
            int[] elementData = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };// the last line will not be used. We assign only one element

            var Hexa8 = elementFactory3D.CreateElement(CellType.Hexa8, nodes);
            var element = new Element();
            element.ID = 0;
            element.ElementType = Hexa8;
            for (int j = 0; j < 8; j++)
            {
                element.NodesDictionary.Add(elementData[j], model.NodesDictionary[elementData[j]]);
            }
            model.SubdomainsDictionary[0].Elements.Add(element);
            model.ElementsDictionary.Add(0, element);

            var flux1 = new FluxLoad(1);
            var flux2 = new FluxLoad(10);
            var dir1 = new DirichletDistribution(list => {
                return Vector.CreateWithValue(list.Count, 0);
            });
            var dir2 = new DirichletDistribution(list => {
                return Vector.CreateWithValue(list.Count, 10);
            });
            var weakDirichlet1 = new WeakDirichlet(dir1,k);
            var weakDirichlet2 = new WeakDirichlet(dir2,k);

            var dirichletFactory1 = new SurfaceLoadElementFactory(weakDirichlet1);
            var dirichletFactory2 = new SurfaceLoadElementFactory(weakDirichlet2);
            var fluxFactory1 = new SurfaceLoadElementFactory(flux1);
            var fluxFactory2 = new SurfaceLoadElementFactory(flux2);
            //foreach (int i in new int[] { 1 })
            //{
            //    var Quad = boundaryFactory3D.CreateElement(CellType.Quad4, face[i - 1]);
            //    var faceElement = new Element();
            //    faceElement.ID = i;
            //    faceElement.ElementType = Quad;
            //    foreach (Node node in face[i - 1])
            //    {
            //        faceElement.AddNode(node);
            //    }
            //    model.SubdomainsDictionary[0].Elements.Add(faceElement);
            //    model.ElementsDictionary.Add(i, faceElement);

            //    var dirichletElement = dirichletFactory1.CreateElement(CellType.Quad4, face[i - 1]);
            //    //var fluxElement = fluxFactory.CreateElement(CellType.Quad4, face[i - 1]);

            //    model.SurfaceLoads.Add(dirichletElement);
            //    //model.SurfaceLoads.Add(fluxElement);
            //}

            foreach (int i in new int[] { 0 })
            {
                //var Quad = boundaryFactory3D.CreateElement(CellType.Quad4, face[i]);
                //var faceElement = new Element();
                //faceElement.ID = 1;
                //faceElement.ElementType = Quad;
                //foreach (Node node in face[i])
                //{
                //    faceElement.AddNode(node);
                //}
                //model.SubdomainsDictionary[0].Elements.Add(faceElement);
                //model.ElementsDictionary.Add(1, faceElement);

                //var dirichletElement = dirichletFactory2.CreateElement(CellType.Quad4, face[i]);
                //model.SurfaceLoads.Add(dirichletElement);

                var fluxElement = fluxFactory1.CreateElement(CellType.Quad4, face[i]);
                model.SurfaceLoads.Add(fluxElement);
            }

            foreach (int i in new int[] { 1 })
            {
                //var Quad = boundaryFactory3D.CreateElement(CellType.Quad4, face[i]);
                //var faceElement = new Element();
                //faceElement.ID = 2;
                //faceElement.ElementType = Quad;
                //foreach (Node node in face[i])
                //{
                //    faceElement.AddNode(node);
                //}
                //model.SubdomainsDictionary[0].Elements.Add(faceElement);
                //model.ElementsDictionary.Add(2, faceElement);

                //var dirichletElement = dirichletFactory2.CreateElement(CellType.Quad4, face[i]);
                //model.SurfaceLoads.Add(dirichletElement);

                var fluxElement = fluxFactory2.CreateElement(CellType.Quad4, face[i]);
                model.SurfaceLoads.Add(fluxElement);
            }

            //constraints
            //foreach (int i in new int[] { 4, 5, 6, 7 })
            //{
            //    model.NodesDictionary[i].Constraints.Add(new Constraint()
            //    {
            //        Amount = 0,
            //        DOF = ThermalDof.Temperature
            //    });
            //}

            //// loads
            //Load load1;
            //foreach (int i in new int[] {0, 1, 2, 3})
            //{
            //    load1 = new Load()
            //    {
            //        Node = model.NodesDictionary[i],
            //        DOF = ThermalDof.Temperature,
            //        Amount = 12.5
            //    };
            //    model.Loads.Add(load1);
            //}


            return model;
        }



        private static IVectorView SolveModel(Model model)
        {
            double[] temp0 = new double[] { 0, 0, 0, 0, 0, 0, 0, 0 };
            Vector initialTemp = Vector.CreateFromArray(temp0);
            DenseMatrixSolver solver = (new DenseMatrixSolver.Builder()).BuildSolver(model);
            var provider = new ProblemConvectionDiffusion(model, solver);

            var childAnalyzer = new LinearAnalyzer(model, solver, provider); // NonlinearAnalyzer
            var parentAnalyzer = new ConvectionDiffusionExplicitDynamicAnalyzer(model, solver, provider, childAnalyzer, 2.5e-4, 2.5, initialTemp); // ConvectionDiffusionStaticAnalyzer

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.LinearSystems[subdomainID].Solution;
        }
    }
}

