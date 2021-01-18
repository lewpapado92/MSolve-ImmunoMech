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

namespace ISAAR.MSolve.Tests.FEM
{
    public class TumorCellsConcentrationDynamic1
    {
        private const int subdomainID = 0;
        [Fact]
        private static void RunTest()
        {
            Model model = CreateModel().Item1;
            ComsolMeshReader modelReader = CreateModel().Item2;
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

        private static Tuple<Model,ComsolMeshReader> CreateModel()
        {
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "mesh1.mphtxt");
            ComsolMeshReader modelReader = new ComsolMeshReader(filename);
            Model model = modelReader.CreateModelFromFile();
            //Boundary Conditions
            var flux = new FluxLoad(10);
            var dir1 = new DirichletDistribution(list => {
                return Vector.CreateWithValue(list.Count, 0);
            });
            var dir2 = new DirichletDistribution(list => {
                return Vector.CreateWithValue(list.Count, 1);
            });
            var weakDirichlet1 = new WeakDirichlet(dir1,modelReader.diffusionCeoff);
            var weakDirichlet2 = new WeakDirichlet(dir2, modelReader.diffusionCeoff);

            var dirichletFactory1 = new SurfaceLoadElementFactory(weakDirichlet1);
            var dirichletFactory2 = new SurfaceLoadElementFactory(weakDirichlet2);
            var fluxFactory = new SurfaceLoadElementFactory(flux);

            int[] boundaryIDs = new int[] { 2, 7, 8, 9 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Element element in modelReader.elementBoundaries[boundaryID])
                {
                    IReadOnlyList<Node> nodes = (IReadOnlyList<Node>)element.Nodes; 
                    var dirichletElement1 = dirichletFactory1.CreateElement(CellType.Tri3, nodes);
                    model.SurfaceLoads.Add(dirichletElement1);
                    //var surfaceElement = new SurfaceLoadElement();
                    //element.ID = TriID;
                    //surfaceElement.ElementType = DirichletElement1;
                    //model.SubdomainsDictionary[0].Elements.Add(dirichletElement1);
                    //model.ElementsDictionary.Add(TriID, surfaceElement);

                    //model.NodesDictionary[surfaceElement.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100 });
                }
            }
            boundaryIDs = new int[]{ 0, 1, 3, 4, 6 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Element element in modelReader.elementBoundaries[boundaryID])
                {
                    IReadOnlyList<Node> nodes = (IReadOnlyList<Node>)element.Nodes;
                    var dirichletElement2 = dirichletFactory2.CreateElement(CellType.Tri3, nodes);
                    model.SurfaceLoads.Add(dirichletElement2);
                }
            }

            return new Tuple<Model,ComsolMeshReader>(model,modelReader);
        }

        private static IVectorView SolveModel(Model model, ComsolMeshReader modelReader)
        {
            double[] temp0 = new double[model.Nodes.Count];
            int[] boundaryIDs = new int[] {2, 7, 8, 9 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Node node in modelReader.nodeBoundaries[boundaryID])
                {
                    temp0[node.ID] = 0;
                }
            }
            boundaryIDs = new int[] { 0, 1, 3, 4, 6 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Node node in modelReader.nodeBoundaries[boundaryID])
                {
                    temp0[node.ID] = 1;
                }
            }
            Vector initialTemp = Vector.CreateFromArray(temp0);
            DenseMatrixSolver solver = (new DenseMatrixSolver.Builder()).BuildSolver(model);
            var provider = new ProblemConvectionDiffusion(model, solver);

            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new ConvectionDiffusionExplicitDynamicAnalyzer(model, solver, provider, childAnalyzer, 5e-8, 1e-3, initialTemp);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return parentAnalyzer.temperature[subdomainID];
 //           return solver.LinearSystems[subdomainID].Solution;
        }
    }
}

