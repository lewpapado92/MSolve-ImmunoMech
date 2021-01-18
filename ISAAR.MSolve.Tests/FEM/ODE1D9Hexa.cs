using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Readers;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;
using System.IO;
using ISAAR.MSolve.FEM.Loading.SurfaceLoads;
using static ISAAR.MSolve.FEM.Loading.SurfaceLoads.WeakDirichlet;
using ISAAR.MSolve.FEM.Loading;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System;
using ISAAR.MSolve.FEM.Elements.BoundaryConditionElements;
using ISAAR.MSolve.FEM.Loading.BodyLoads;

namespace ISAAR.MSolve.Tests.FEM
{
    public class ODE1D9Hexa
    {
        private const int subdomainID = 0;
        [Fact]
        private static void RunTest()
        {
            Model model = CreateModel(1, 1).Item1;
            ComsolMeshReader3 modelReader = CreateModel(1, 1).Item2;
            IVectorView solution = SolveModel(model, modelReader);
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

        private static Tuple<Model, ComsolMeshReader3> CreateModel(double C1, double C0)
        {
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "9hexa.mphtxt");
            ComsolMeshReader3 modelReader = new ComsolMeshReader3(filename, C1, C0);
            Model model = modelReader.CreateModelFromFile();
            //Boundary Conditions
            var bodyLoad = new ODEDomainLoad(new ODEMaterial(C1, C0), 1, ThermalDof.Temperature);
            var flux1 = new FluxLoad(1);
            var flux2 = new FluxLoad(10);
            var bodyLoadElementFactory = new BodyLoadElementFactory(bodyLoad, model);
            var fluxFactory1 = new SurfaceLoadElementFactory(flux1);
            var fluxFactory2 = new SurfaceLoadElementFactory(flux2);

            int[] domainIDs = new int[] { 0, };
            foreach (int domainID in domainIDs)
            {
                foreach (Element element in modelReader.elementDomains[domainID])
                {
                    //IReadOnlyList<Node> nodes = (IReadOnlyList<Node>)element.Nodes;
                    //var fluxElement1 = fluxFactory1.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(fluxElement1);
                    //var bodyLoadElementCellType = element.ElementType.CellType;
                    var nodes = (IReadOnlyList<Node>)element.Nodes;
                    var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Hexa8, nodes);
                    model.BodyLoads.Add(bodyLoadElement);
                    //var surfaceElement = new SurfaceLoadElement();
                    //element.ID = TriID;
                    //surfaceElement.ElementType = DirichletElement1;
                    //model.SubdomainsDictionary[0].Elements.Add(dirichletElement1);
                    //model.ElementsDictionary.Add(TriID, surfaceElement);

                    //model.NodesDictionary[surfaceElement.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100 });
                }
            }
            int[] boundaryIDs = new int[] { 0, };
            int QuadID = model.ElementsDictionary.Count + 1;
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    //IReadOnlyList<Node> nodes = (IReadOnlyList<Node>)element.Nodes;
                    //var fluxElement1 = fluxFactory1.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(fluxElement1);
                    //var dirichletElement1 = dirichletFactory1.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(dirichletElement1);
                    //var surfaceBoundaryElement = boundaryFactory3D.CreateElement(CellType.Quad4, nodes);
                    //var element = new Element();
                    //element.ID = QuadID;
                    //element.ElementType = surfaceBoundaryElement;
                    //model.SubdomainsDictionary[0].Elements.Add(element);
                    //model.ElementsDictionary.Add(QuadID, element);
                    //foreach (Node node in nodes)
                    //{
                    //    element.AddNode(node);
                    //}
                    //QuadID += 1;
                    // var surfaceElement = new SurfaceLoadElement();
                    //element.ID = TriID;
                    //surfaceElement.ElementType = DirichletElement1;
                    //model.SubdomainsDictionary[0].Elements.Add(dirichletElement1);
                    //model.ElementsDictionary.Add(TriID, surfaceElement);

                }
                //foreach (Node node in modelReader.nodeBoundaries[boundaryID])
                //{
                //    model.NodesDictionary[node.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 0 });
                //}
            }
            boundaryIDs = new int[] { 5 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    //IReadOnlyList<Node> nodes = (IReadOnlyList<Node>)element.Nodes;
                    //var fluxElement2 = fluxFactory2.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(fluxElement2);
                    //var dirichletElement2 = dirichletFactory2.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(dirichletElement2);
                    //var surfaceBoundaryElement = boundaryFactory3D.CreateElement(CellType.Quad4, nodes);
                    //var element = new Element();
                    //element.ID = QuadID;
                    //element.ElementType = surfaceBoundaryElement;
                    //model.SubdomainsDictionary[0].Elements.Add(element);
                    //model.ElementsDictionary.Add(QuadID, element);
                    //foreach (Node node in nodes)
                    //{
                    //    element.AddNode(node);
                    //}
                    //QuadID += 1;
                    // var surfaceElement = new SurfaceLoadElement();
                    //element.ID = TriID;
                    //surfaceElement.ElementType = DirichletElement1;
                    //model.SubdomainsDictionary[0].Elements.Add(dirichletElement1);
                    //model.ElementsDictionary.Add(TriID, surfaceElement);

                }
                //foreach (Node node in modelReader.nodeBoundaries[boundaryID])
                //{
                //    model.NodesDictionary[node.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 0 });
                //}
            }

            //boundaryIDs = new int[] { 1, 2, 3, 4 };
            //foreach (int boundaryID in boundaryIDs)
            //{
            //    foreach (Element element in modelReader.elementBoundaries[boundaryID])
            //    {
            //        IReadOnlyList<Node> nodes = (IReadOnlyList<Node>)element.Nodes;
            //        var fluxElement = fluxFactory.CreateElement(CellType.Quad4, nodes);
            //        model.SurfaceLoads.Add(fluxElement);
            //    }
            //}

            return new Tuple<Model, ComsolMeshReader3>(model, modelReader);
        }

        private static IVectorView SolveModel(Model model, ComsolMeshReader3 modelReader)
        {
            //var initialTemp = Vector.CreateZero(model.Nodes.Count);
            double[] value0 = new double[model.Nodes.Count];
            int[] boundaryIDs = new int[] { 0 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
                    {
                        value0[node.ID] = 0;
                    }
                }
            }
            boundaryIDs = new int[] { 5 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
                    {
                        value0[node.ID] = 0;
                    }
                }
            }
            Vector initialValue = Vector.CreateFromArray(value0);
            var builder = new SkylineSolver.Builder();
            //builder.IsMatrixPositiveDefinite = false;
            var solver = builder.BuildSolver(model);
            var provider = new ProblemODE(model, solver);

            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new ODEDynamicAnalyzer(model, solver, provider, childAnalyzer, .5,  .1, 3, initialValue);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.LinearSystems[subdomainID].Solution;
        }
    }
}

