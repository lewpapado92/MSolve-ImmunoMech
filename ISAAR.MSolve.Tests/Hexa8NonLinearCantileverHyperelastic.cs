using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISSAR.MSolve.Discretization.Loads;
using Xunit;

namespace ISAAR.MSolve.Tests.FEM
{
    public static class Hexa8NonLinearCantileverDefGradHyperelastic
    {
        private const int subdomainID = 1;

        private static bool AreDisplacementsSame(IReadOnlyList<Dictionary<int, double>> expectedDisplacements, TotalDisplacementsPerIterationLog computedDisplacements)
        {
            var comparer = new ValueComparer(1E-13);
            for (int iter = 0; iter < expectedDisplacements.Count; ++iter)
            {
                foreach (int dof in expectedDisplacements[iter].Keys)
                {
                    if (!comparer.AreEqual(expectedDisplacements[iter][dof], computedDisplacements.GetTotalDisplacement(iter, subdomainID, dof)))
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        [Fact]
        private static void RunTest()
        {
            IReadOnlyList<Dictionary<int, double>> expectedDisplacements = GetExpectedDisplacements();
            TotalDisplacementsPerIterationLog computedDisplacements = SolveModel();
            Assert.True(AreDisplacementsSame(expectedDisplacements, computedDisplacements));
        }

        private static IReadOnlyList<Dictionary<int, double>> GetExpectedDisplacements()
        {
            var expectedDisplacements = new Dictionary<int, double>[6]; //TODO: this should be 11 EINAI ARRAY APO DICTIONARIES

            expectedDisplacements[0] = new Dictionary<int, double> {
      {0, -0.012630163044092212}, {5,-0.052941878275004031} };
            expectedDisplacements[1] = new Dictionary<int, double> {
    { 0,-0.012356629242693141 }, {5,-0.05341157042273393 }};
            expectedDisplacements[2] = new Dictionary<int, double> {
    { 0,-0.012357034444187773 }, {5,-0.053421028639746654}};
            expectedDisplacements[3] = new Dictionary<int, double> {
    { 0,-0.024423182208720692 }, {5,-0.10728250571358597}};
            expectedDisplacements[4] = new Dictionary<int, double> {
    { 0,-0.02409681387875337 }, {5,-0.10763333288739485 }};
            expectedDisplacements[5] = new Dictionary<int, double> {
    { 0,-0.024097071796207306 }, {5,-0.10764279522486279 }};
            


            return expectedDisplacements;
        }

        private static TotalDisplacementsPerIterationLog SolveModel()
        {
            //VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain( subdomainID));

            BuildCantileverModel(model, -(0.15*2)/20);

            //model.ConnectDataStructures();

            //var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
            //linearSystems[subdomainID] = new SkylineLinearSystem(subdomainID, model.Subdomains[0].Forces);

            //ProblemStructural provider = new ProblemStructural(model, linearSystems);

            //var solver = new SolverSkyline(linearSystems[subdomainID]);
            //var linearSystemsArray = new[] { linearSystems[subdomainID] };
            //var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            //var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };

            //var increments = 2;
            //var childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Problem type
            var provider = new ProblemStructural(model, solver);

            //var solver = new SolverSkyline(linearSystems[subdomainID]);
            //var linearSystemsArray = new[] { linearSystems[subdomainID] };

            //var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            //var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };

            var increments = 2;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            childAnalyzerBuilder.ResidualTolerance = 1E-8;
            childAnalyzerBuilder.MaxIterationsPerIncrement = 100;
            childAnalyzerBuilder.NumIterationsForMatrixRebuild = 1;
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            var watchDofs = new Dictionary<int, int[]>();
            watchDofs.Add(subdomainID, new int[2] { 0, 5 }); //[12] { 0, 1,2,3,4,5,6,7,8,9,10,11});
            var log1 = new TotalDisplacementsPerIterationLog(watchDofs);
            childAnalyzer.TotalDisplacementsPerIterationLog = log1;


            //childAnalyzer.SetMaxIterations = 100;
            //childAnalyzer.SetIterationsForMatrixRebuild = 1;

            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            //parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();


            return log1;
        }

        private static void BuildCantileverModel(Model model, double load_value)
        {
            //xrhsimopoiithike to  Hexa8NonLinearCantileverDefGrad
            // allagh tou material


            IContinuumMaterial3DDefGrad material1 = new HyperElasticMaterial3DDefGrad() { C1 = 0.035, C2 = 0.057, k_cons = 1 };

            double[,] nodeData = new double[,] { {1,-1,-1},
            {1,1,-1},
            {1,-1,1},
            {1,1,1},
            {-1,-1,-1},
            {-1,1,-1},
            {-1,-1,1},
            {-1,1,1} };

            int[,] elementData = new int[,] {{1,4,8,7,3,2,6,5,1},
            {2,12,11,9,10,8,7,5,6} };// the last line will not be used. We assign only one element

            // orismos shmeiwn
            for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nNode + 1, new Node(id: nNode + 1, x: nodeData[nNode, 0], y:  nodeData[nNode, 1], z: nodeData[nNode, 2] ));

            }

            // orismos elements 
            Element e1;
            int subdomainID = 1;
            for (int nElement = 0; nElement < elementData.GetLength(0)-1; nElement++)
            {
                e1 = new Element()
                {
                    ID = nElement + 1,
                    ElementType = new Hexa8NonLinearDefGrad(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8                    
                };
                for (int j = 0; j < 8; j++)
                {
                    e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
                model.SubdomainsDictionary[subdomainID].Elements.Add( e1);
            }

            // constraint to to deftero miso apo th list twn nodes
            foreach(int k in new int[] { 5,6,7,8})
            {
                model.NodesDictionary[k].Constraints.Add(new Constraint()
                {
                    Amount = 0,
                    DOF = StructuralDof.TranslationX
                });
                model.NodesDictionary[k].Constraints.Add(new Constraint()
                {
                    Amount = 0,
                    DOF = StructuralDof.TranslationY
                });
                model.NodesDictionary[k].Constraints.Add(new Constraint()
                {
                    Amount = 0,
                    DOF = StructuralDof.TranslationZ
                });
            }

            // fortish korufhs
            Load load1;
            for (int k = 4; k < 5; k++)
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[k],
                    DOF = StructuralDof.TranslationZ,
                    Amount = 1 * load_value
                };
                model.Loads.Add(load1);
            }
        }
    }

}
