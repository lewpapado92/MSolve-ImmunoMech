using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISSAR.MSolve.Discretization.Loads;
using Xunit;

namespace ISAAR.MSolve.Tests
{
    public class Beam2DNewmarkDynamicAanalysisMultiTest
    {
        private const int subdomainID = 0;
        private static SkylineSolver.Builder builder = new SkylineSolver.Builder();

        [Fact]
        private static void RunTest()
        {
            //var models = new[] { CreateModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item1, CreateModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item1 };
            //var modelReaders = new[] { CreateModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item2, CreateModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item2 };

            var models = new[] { CreateModel(21000, 91.04, 1000.0), CreateModel(21000, 91.04, 1000.0) };
            IVectorView[] solutions = SolveModels(models);
            Assert.True(CompareResults(solutions[0]));
        }

        private static void UpdateModels(Dictionary<int, IVector>[] solutions, IStructuralModel[] modelsToReplace, ISolver[] solversToReplace,
            IImplicitIntegrationProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
        {
            double[] sol0 = solutions[0][0].CopyToArray();
            double[] sol1 = solutions[1][0].CopyToArray();
            modelsToReplace[0] = CreateModel(21000, 91.04, 1000.0);
            modelsToReplace[1] = CreateModel(21000, 91.04, 1000.0);

            for (int i = 0; i < modelsToReplace.Length; i++)
            {
                //providersToReplace[i].
                solversToReplace[i] = builder.BuildSolver(modelsToReplace[i]);
                providersToReplace[i] = new ProblemStructural((Model)modelsToReplace[i], solversToReplace[i]);
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

        private static Model CreateModel(double youngModulus, double area, double nodalLoad)
        {
            //double youngModulus = 21000;
            double poissonRatio = 0.3;
            //double area = 91.04;
            double inertiaY = 2843.0;
            double inertiaZ = 8091.0;
            double density = 7.85;
            //double nodalLoad = 1000.0;
            int totalNodes = 2;

            var material = new ElasticMaterial()
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Node creation
            IList<Node> nodes = new List<Node>();
            Node node0 = new Node(id: 0, x: 0.0, y: 0.0, z: 0.0);
            Node node1 = new Node(id: 1, x: 300.0, y: 0.0, z: 0.0);
            nodes.Add(node0);
            nodes.Add(node1);

            // Model creation
            Model model = new Model();

            // Add a single subdomain to the model
            model.SubdomainsDictionary.Add(0, new Subdomain(0));

            // Add nodes to the nodes dictonary of the model
            for (int i = 0; i < nodes.Count; ++i)
            {
                model.NodesDictionary.Add(i, nodes[i]);
            }

            // Constrain bottom nodes of the model
            model.NodesDictionary[0].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
            model.NodesDictionary[0].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
            model.NodesDictionary[0].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });

            // Create a new Beam2D element
            var beam = new EulerBeam2D(youngModulus)
            {
                Density = density,
                SectionArea = area,
                MomentOfInertia = inertiaZ
            };

            var element = new Element()
            {
                ID = 0,
                ElementType = beam
            };

            // Add nodes to the created element
            element.AddNode(model.NodesDictionary[0]);
            element.AddNode(model.NodesDictionary[1]);

            // Element Stiffness Matrix
            var a = beam.StiffnessMatrix(element);
            var b = beam.MassMatrix(element);

            // Add Hexa element to the element and subdomains dictionary of the model
            model.ElementsDictionary.Add(element.ID, element);
            model.SubdomainsDictionary[0].Elements.Add(element);

            // define loads
            model.Loads.Add(new Load { Amount = nodalLoad, Node = model.NodesDictionary[totalNodes-1], DOF = StructuralDof.TranslationY });

            return model;
        }


        private static IVectorView[] SolveModels(Model[] models)
        {

            SkylineSolver[] solvers = new SkylineSolver[models.Length];
            IImplicitIntegrationProvider[] providers = new IImplicitIntegrationProvider[models.Length];
            IChildAnalyzer[] childAnalyzers = new IChildAnalyzer[models.Length];
            for (int i = 0; i < models.Length; i++)
            {
                //var builder = new DenseMatrixSolver.Builder();
                solvers[i] = builder.BuildSolver(models[i]);
                providers[i] = new ProblemStructural(models[i], solvers[i]);
                childAnalyzers[i] = new LinearAnalyzer(models[i], solvers[i], providers[i]);
            }

            var parentAnalyzer = new NewmarkDynamicAnalyzerMultiModel(UpdateModels, models, solvers, providers, childAnalyzers, 0.28, 3.36, .25, .5);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solvers.Select(x => x.LinearSystems[subdomainID].Solution).ToArray();
        }
    }
}