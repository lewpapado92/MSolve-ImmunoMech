using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;
using Accord;
using ISAAR.MSolve.Discretization.Interfaces;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Interfaces;
using System.Linq;
using ISAAR.MSolve.Solvers;
using ISSAR.MSolve.Discretization.Loads;

namespace ISAAR.MSolve.Tests.FEM
{
    public class ConvectionDiffusionMulti1DBenchmark
    {
        private const int subdomainID = 0;

        [Fact]
        private static void RunTest()
        {
            var models = new[] { CreateModel(1, new double[] { 2, 2, 2 }, 0, 1), CreateModel(0, new double[] { 2, 2, 2 }, 0, 1) };
            IVectorView[] solutions = SolveModels(models);
            Assert.True(CompareResults(solutions[0]));
        }

        private static void UpdateModels(Dictionary<int, IVector>[] solutions, IStructuralModel[] modelsToReplace, ISolver[] solversToReplace,
            IConvectionDiffusionIntegrationProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
        {
            modelsToReplace[0] = CreateModel(1, new double[] { 2, 2, 2 }, 0, 1);
            modelsToReplace[1] = CreateModel(0, new double[] { 2, 2, 2 }, 0, 1);

            for (int i = 0; i < modelsToReplace.Length; i++)
            {
                solversToReplace[i] = (new DenseMatrixSolver.Builder()).BuildSolver(modelsToReplace[i]);
                providersToReplace[i] = new ProblemConvectionDiffusion((Model)modelsToReplace[i], solversToReplace[i]);
                childAnalyzersToReplace[i] = new LinearAnalyzer(modelsToReplace[i], solversToReplace[i], providersToReplace[i]);
            }
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

        private static Model CreateModel(double k, double[] U, double L, double h)
        {
            var model = new Model();

            // Subdomains
            model.SubdomainsDictionary.Add(0, new Subdomain(subdomainID));

            double crossSectionArea = 1;
            var material = new ConvectionDiffusionMaterial(k, U, L);

            // Nodes
            int numNodes = 10;
            for (int i = 0; i < numNodes; i++)
            {
                model.NodesDictionary.Add(i, new Node(id: i, x: i * h));
            }

            //Elements
            var numElements = numNodes - 1;
            for (int i = 0; i < numElements; i++)
            {
                Node[] startEndNodes = { model.NodesDictionary[i], model.NodesDictionary[i+1] };
                var elementType = new ConvectionDiffusionRod(startEndNodes, crossSectionArea, material);
                var elementWrapper = new Element() { ID = i, ElementType = elementType };
                foreach (var node in startEndNodes) elementWrapper.AddNode(node);
                model.ElementsDictionary[elementWrapper.ID] = elementWrapper;
                model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
            }

            // Dirichlet BC
            //double kappa = 1;
            //double startNodeTemperature = 1;
            //double endNodeTemperature = 0;
            //model.NodesDictionary[0].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = kappa * k / h *startNodeTemperature });
            //model.NodesDictionary[9].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = kappa * k / h * endNodeTemperature });

            // Dirichlet BC
            double kappa = 1;
            double startNodeTemperature = 1;
            double endNodeTemperature = 0;
            model.Loads.Add(new Load() { Amount = kappa * k / h * startNodeTemperature, Node = model.NodesDictionary[0], DOF = ThermalDof.Temperature });
            model.Loads.Add(new Load() { Amount = kappa * k / h * endNodeTemperature, Node = model.NodesDictionary[9], DOF = ThermalDof.Temperature });
            //model.Loads.Add(new Load() { Amount = q / 2.0, Node = model.NodesDictionary[8], DOF = ThermalDof.Temperature });


            ////////model.TimeDependentNodalLoads.Add(new SteadyNodalLoad(q / 2.0) { Node = model.NodesDictionary[2], DOF = DOFType.Temperature });

            return model;
        }

        private static IVectorView[] SolveModels(Model[] models)
        {
            Vector[] initialTemps = new Vector[models.Length];
            DenseMatrixSolver[] solvers = new DenseMatrixSolver[models.Length];
            IConvectionDiffusionIntegrationProvider[] providers = new IConvectionDiffusionIntegrationProvider[models.Length];
            IChildAnalyzer[] childAnalyzers = new IChildAnalyzer[models.Length];
            for (int i = 0; i < models.Length; i++)
            {
                initialTemps[i] = Vector.CreateFromArray(new double[] { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
                solvers[i] = (new DenseMatrixSolver.Builder()).BuildSolver(models[i]);
                providers[i] = new ProblemConvectionDiffusion(models[i], solvers[i]);
                childAnalyzers[i] = new LinearAnalyzer(models[i], solvers[i], providers[i]);
            }

            var parentAnalyzer = new ConvectionDiffusionDynamicAnalyzerMultiModel_Beta(UpdateModels, models, solvers, providers, childAnalyzers, 0.05, 5, initialTemperature: initialTemps);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solvers.Select(x => x.LinearSystems[subdomainID].Solution).ToArray();
        }
    }
}