namespace ISAAR.MSolve.Tests
{
	using System.Collections.Generic;
	using System.IO;
	using ISAAR.MSolve.Discretization;
	using ISAAR.MSolve.Discretization.Integration.Quadratures;
	using ISAAR.MSolve.Discretization.Mesh;
	using ISAAR.MSolve.FEM.Entities;
	using ISAAR.MSolve.FEM.Loading;
	using ISAAR.MSolve.FEM.Interpolation;
	using ISAAR.MSolve.Solvers.Direct;
	using Xunit;
    using ISAAR.MSolve.Materials;
    using ISAAR.MSolve.Analyzers.Dynamic;
    using ISAAR.MSolve.Analyzers;
    using ISAAR.MSolve.Problems;
    using ISAAR.MSolve.Logging;
    using ISAAR.MSolve.Analyzers.NonLinear;
    using ISAAR.MSolve.Discretization.FreedomDegrees;
    using ISAAR.MSolve.Materials.Interfaces;
    using ISAAR.MSolve.FEM.Elements;
    using ISSAR.MSolve.Discretization.Loads;

    public class ContinuumElements3DDynamicTests
	{
		private const int subdomainID = 0;
		private static int numberOfNodes;
		private const int increments = 1;
		private const double youngModulus = 4.0;
		private const double poissonRatio = 0.4;
		private const double yieldStress = 0.12;
		private const double plasticModulus = 0.1;
		private const double density = 0.001;
		private const double hardeningRation = plasticModulus / youngModulus;
		private const int monitorDof = 11;
		private const double lambdag = 1;

        //Linear Geometry Examples
       [Fact]
        public void ContinuumElement3DElasticMaterialDynamicConsistent()
        {
            IList<Node> nodes = new List<Node>();

            // Create Model
            var model = new Model();

            // Create Subdomain
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

            // Create Elastic Material
            var solidMaterial = new ElasticMaterial3D()
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Create Dynamic Material
            var dynamicMaterial = new DynamicMaterial(density, 0, 0, true);

            DefineContinuumElement3DLinear(model, solidMaterial, dynamicMaterial);

            BoundaryConditions(model);

            // Choose linear equation system solver
            var solverBuilder = new SkylineSolver.Builder();
            SkylineSolver solver = solverBuilder.BuildSolver(model);

            // Choose the provider of the problem -> here a structural problem
            var provider = new ProblemStructural(model, solver);

            // Choose child analyzer -> Linear
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);

            // Choose parent analyzer -> Parent: Static
            var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 1.0, 100.0);
            parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration();
            NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

            // Request output
            childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Check output
            DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0];
            double computedValue = log.DOFValues[monitorDof];
            double expectedValue = 0.531178; // consistent: 0.531178 // lumped: 0.894201
            Assert.Equal(expectedValue, computedValue, 3);
        }

        //[Fact]
        //public void ContinuumElement3DElasticMaterialDynamicLumped()
        //{
        //	IList<Node> nodes = new List<Node>();

        //	// Create Model
        //	var model = new Model();

        //	// Create Subdomain
        //	model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

        //	// Create Elastic Material
        //	var solidMaterial = new ElasticMaterial3D()
        //	{
        //		YoungModulus = youngModulus,
        //		PoissonRatio = poissonRatio,
        //	};

        //	// Create Dynamic Material
        //	var dynamicMaterial = new DynamicMaterial(density, 0, 0, false);

        //	DefineContinuumElement3DLinear(model, solidMaterial, dynamicMaterial);

        //	BoundaryConditions(model);

        //	// Choose linear equation system solver
        //	var solverBuilder = new SkylineSolver.Builder();
        //	SkylineSolver solver = solverBuilder.BuildSolver(model);

        //	// Choose the provider of the problem -> here a structural problem
        //	var provider = new ProblemStructural(model, solver);

        //	// Choose child analyzer -> Linear
        //	var childAnalyzer = new LinearAnalyzer(model, solver, provider);

        //	// Choose parent analyzer -> Parent: Static
        //	var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 1.0, 100.0);
        //	parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration();
        //	NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

        //	// Request output
        //	childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

        //	// Run the analysis
        //	parentAnalyzer.Initialize();
        //	parentAnalyzer.Solve();

        //	// Check output
        //	DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0];
        //	double computedValue = log.DOFValues[monitorDof];
        //	double expectedValue = 0.894201;
        //	Assert.Equal(expectedValue, computedValue, 3);
        //}

        [Fact]
        public void ContinuumElement3DVonMisesMaterialDynamicConsistent()
        {
            IList<Node> nodes = new List<Node>();

            // Create Model
            var model = new Model();

            // Create Subdomain
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

            // Create von Mises Plastic Material
            var solidMaterial = new VonMisesMaterial3D(youngModulus, poissonRatio, yieldStress, plasticModulus);

            // Create Dynamic Material
            var dynamicMaterial = new DynamicMaterial(density, 0, 0, true);

            DefineContinuumElement3DLinear(model, solidMaterial, dynamicMaterial);

            BoundaryConditionsNLM(model);

            // Choose linear equation system solver
            var solverBuilder = new SkylineSolver.Builder();
            SkylineSolver solver = solverBuilder.BuildSolver(model);

            // Choose the provider of the problem -> here a structural problem
            var provider = new ProblemStructural(model, solver);

            // Choose child analyzer -> Child: NewtonRaphson - LoadControl
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) };
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments)
            {
                MaxIterationsPerIncrement = 50,
                NumIterationsForMatrixRebuild = 1,
                ResidualTolerance = 1E-06,
            };
            var childAnalyzer = childAnalyzerBuilder.Build();

            // Choose parent analyzer -> Parent: Static
            var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 1.0, 100.0);
            parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration();
            NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

            // Request output
            childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Check output
            DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0];
            double computedValue = log.DOFValues[monitorDof];
            double expectedValue = 3.47644;
            Assert.Equal(expectedValue, computedValue, 2);
        }

        //[Fact]
        //public void ContinuumElement3DVonMisesMaterialDynamicLumped()
        //{
        //	IList<Node> nodes = new List<Node>();

        //	// Create Model
        //	var model = new Model();

        //	// Create Subdomain
        //	model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

        //	// Create von Mises Plastic Material
        //	var solidMaterial = new VonMisesMaterial3D(youngModulus, poissonRatio, yieldStress, plasticModulus);

        //	// Create Dynamic Material
        //	var dynamicMaterial = new DynamicMaterial(density, 0, 0, false);

        //	DefineContinuumElement3DLinear(model, solidMaterial, dynamicMaterial);

        //	BoundaryConditionsNLM(model);

        //	// Choose linear equation system solver
        //	var solverBuilder = new SkylineSolver.Builder();
        //	SkylineSolver solver = solverBuilder.BuildSolver(model);

        //	// Choose the provider of the problem -> here a structural problem
        //	var provider = new ProblemStructural(model, solver);

        //	// Choose child analyzer -> Child: NewtonRaphson - LoadControl
        //	var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) };
        //	var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments)
        //	{
        //		MaxIterationsPerIncrement = 50,
        //		NumIterationsForMatrixRebuild = 1,
        //		ResidualTolerance = 1E-06,
        //	};
        //	var childAnalyzer = childAnalyzerBuilder.Build();

        //	// Choose parent analyzer -> Parent: Dynamic
        //	var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 1.0, 100.0);
        //	parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration();
        //	NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

        //	// Request output
        //	childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

        //	// Run the analysis
        //	parentAnalyzer.Initialize();
        //	parentAnalyzer.Solve();

        //	// Check output
        //	DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0];
        //	double computedValue = log.DOFValues[monitorDof];
        //	double expectedValue = 3.78687;
        //	Assert.Equal(expectedValue, computedValue, 3);
        //}

        // NonLinear Geometry Examples
        [Fact]
		public void ContinuumElement3DNonLinearElasticMaterialDynamicConsistent()
		{
			IList<Node> nodes = new List<Node>();

			// Create Model
			var model = new Model();

			// Create Subdomain
			model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

			// Create Elastic Material
			var solidMaterial = new ElasticMaterial3D()
			{
				YoungModulus = youngModulus,
				PoissonRatio = poissonRatio,
			};

			// Create Dynamic Material
			var dynamicMaterial = new DynamicMaterial(density, 0, 0, true);

			DefineContinuumElement3DNonLinear(model, solidMaterial, dynamicMaterial);

			BoundaryConditionsNLM(model);

			// Choose linear equation system solver
			var solverBuilder = new SkylineSolver.Builder();
			SkylineSolver solver = solverBuilder.BuildSolver(model);

			// Choose the provider of the problem -> here a structural problem
			var provider = new ProblemStructural(model, solver);

			// Choose child analyzer -> Child: NewtonRaphson - LoadControl
			var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) };
			var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments)
			{
				MaxIterationsPerIncrement = 50,
				NumIterationsForMatrixRebuild = 1,
				ResidualTolerance = 1E-06,
			};
			var childAnalyzer = childAnalyzerBuilder.Build();

			// Choose parent analyzer -> Parent: Static
			var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 1.0, 100.0);
			parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration();
			NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

			// Request output
			childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			// Check output
			DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0];
			double computedValue = log.DOFValues[monitorDof];
			double expectedValue = 0.500629;
			Assert.Equal(expectedValue, computedValue, 3);
		}

		[Fact]
		public void ContinuumElement3DNonLinearElasticMaterialDynamicLumped()
		{
			IList<Node> nodes = new List<Node>();

			// Create Model
			var model = new Model();

			// Create Subdomain
			model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

			// Create Elastic Material
			var solidMaterial = new ElasticMaterial3D()
			{
				YoungModulus = youngModulus,
				PoissonRatio = poissonRatio,
			};

			// Create Dynamic Material
			var dynamicMaterial = new DynamicMaterial(density, 0, 0, false);

			DefineContinuumElement3DNonLinear(model, solidMaterial, dynamicMaterial);

			BoundaryConditionsNLM(model);

			// Choose linear equation system solver
			var solverBuilder = new SkylineSolver.Builder();
			SkylineSolver solver = solverBuilder.BuildSolver(model);

			// Choose the provider of the problem -> here a structural problem
			var provider = new ProblemStructural(model, solver);

			// Choose child analyzer -> Child: NewtonRaphson - LoadControl
			var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) };
			var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments)
			{
				MaxIterationsPerIncrement = 50,
				NumIterationsForMatrixRebuild = 1,
				ResidualTolerance = 1E-06,
			};
			var childAnalyzer = childAnalyzerBuilder.Build();

			// Choose parent analyzer -> Parent: Static
			var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 1.0, 100.0);
			parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration();
			NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

			// Request output
			childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			// Check output
			DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0];
			double computedValue = log.DOFValues[monitorDof];
			double expectedValue = 2.97084;
			Assert.Equal(expectedValue, computedValue, 3);
		}
		[Fact]
		public void ContinuumElement3DNonLinearVonMisesMaterialDynamicConsistent()
		{
			IList<Node> nodes = new List<Node>();

			// Create Model
			var model = new Model();

			// Create Subdomain
			model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

			// Create von Mises Plastic Material
			var solidMaterial = new VonMisesMaterial3D(youngModulus, poissonRatio, yieldStress, plasticModulus);

			// Create Dynamic Material
			var dynamicMaterial = new DynamicMaterial(density, 0, 0, true);

			DefineContinuumElement3DNonLinear(model, solidMaterial, dynamicMaterial);

			BoundaryConditionsNLM(model);

			// Choose linear equation system solver
			var solverBuilder = new SkylineSolver.Builder();
			SkylineSolver solver = solverBuilder.BuildSolver(model);

			// Choose the provider of the problem -> here a structural problem
			var provider = new ProblemStructural(model, solver);

			// Choose child analyzer -> Child: NewtonRaphson - LoadControl
			var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) };
			var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments)
			{
				MaxIterationsPerIncrement = 50,
				NumIterationsForMatrixRebuild = 1,
				ResidualTolerance = 1E-06,
			};
			var childAnalyzer = childAnalyzerBuilder.Build();

			// Choose parent analyzer -> Parent: Static
			var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 1.0, 100.0);
			parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration();
			NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

			// Request output
			childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			// Check output
			DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0];
			double computedValue = log.DOFValues[monitorDof];
			double expectedValue = 1.93737;
			Assert.Equal(expectedValue, computedValue, 2);
		}

		[Fact]
		public void ContinuumElement3DNonLinearDefGradDynamicConsistent()
		{
			IList<Node> nodes = new List<Node>();

			// Create Model
			var model = new Model();

			// Create Subdomain
			model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

			// Create von Mises Plastic Material
			//var solidMaterial = new VonMisesMaterial3D(youngModulus, poissonRatio, yieldStress, plasticModulus);
			var solidMaterial = new HyperElasticMaterial3DDefGrad() { C1 = 1, C2 = 1, k_cons = 1 };

			// Create Dynamic Material
			var dynamicMaterial = new DynamicMaterial(density, 0, 0, true);

			//DefineContinuumElement3DNonLinear(model, solidMaterial, dynamicMaterial);
			DefineContinuumElement3DNonLinearDefGrad(model, solidMaterial, dynamicMaterial);

			BoundaryConditionsNLM(model);

			// Choose linear equation system solver
			var solverBuilder = new SkylineSolver.Builder();
			SkylineSolver solver = solverBuilder.BuildSolver(model);

			// Choose the provider of the problem -> here a structural problem
			var provider = new ProblemStructural(model, solver);

			// Choose child analyzer -> Child: NewtonRaphson - LoadControl
			var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) };
			var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, 2)
			{
				MaxIterationsPerIncrement = 50,
				NumIterationsForMatrixRebuild = 1,
				ResidualTolerance = 1E-06,
			};
			var childAnalyzer = childAnalyzerBuilder.Build();

			// Choose parent analyzer -> Parent: Static
			var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 1.0, 100.0);
			parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration();
			NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

			// Request output
			childAnalyzer.LogFactories[subdomainID] = new LinearAnalyzerLogFactory(new int[] { monitorDof });

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			// Check output
			DOFSLog log = (DOFSLog)childAnalyzer.Logs[subdomainID][0];
			double computedValue = log.DOFValues[monitorDof];
			double expectedValue = 2.32803;
			Assert.Equal(expectedValue, computedValue, 2);
		}
		// DEFINE ELEMENTS AND BOUNDARY CONDITIONS
		internal static int DefineContinuumElement3DLinear(Model model, IIsotropicContinuumMaterial3D material3D, IDynamicMaterial dynamicMaterial)
		{
			// Node creation
			var node1 = new Node(id: 1, x: 0.0, y: 0.0, z: 0.0);
			var node2 = new Node(id: 2, x: 100.0, y: 0.0, z: 0.0);
			var node3 = new Node(id: 3, x: 0.0, y: 100.0, z: 0.0);
			var node4 = new Node(id: 4, x: 100.0, y: 100.0, z: 0.0);
			var node5 = new Node(id: 5, x: 0.0, y: 0.0, z: 100.0);
			var node6 = new Node(id: 6, x: 100.0, y: 0.0, z: 100.0);
			var node7 = new Node(id: 7, x: 0.0, y: 100.0, z: 100.0);
			var node8 = new Node(id: 8, x: 100.0, y: 100.0, z: 100.0);

			// Create List of nodes
			IList<Node> nodes = new List<Node>();
			nodes.Add(node1);
			nodes.Add(node2);
			nodes.Add(node3);
			nodes.Add(node4);
			nodes.Add(node5);
			nodes.Add(node6);
			nodes.Add(node7);
			nodes.Add(node8);

			numberOfNodes = nodes.Count;

			// Add nodes to the nodes dictonary of the model
			for (int i = 0; i < numberOfNodes; ++i)
			{
				model.NodesDictionary.Add(i, nodes[i]);
			}

			int[][] connectivityMatrix = new int[8][];
			connectivityMatrix[0] = new int[] { 1, 1, 2, 4, 3, 5, 6, 8, 7 };

			var factory = new ContinuumElement3DFactory(material3D, dynamicMaterial);

			for (int i = 0; i < 1; i++)
			{
				List<Node> elementNodeSet = new List<Node>(8);
				for (int j = 1; j < 9; j++)
				{
					elementNodeSet.Add((Node)model.NodesDictionary[connectivityMatrix[i][j] - 1]);
				}

				var hexa8element = new Element()
				{
					ID = connectivityMatrix[i][0],
					ElementType = factory.CreateElement(CellType.Hexa8, elementNodeSet),
				};

				for (int j = 1; j < 9; j++)
				{
					hexa8element.AddNode(model.NodesDictionary[connectivityMatrix[i][j] - 1]);
				}

				model.ElementsDictionary.Add(hexa8element.ID, hexa8element);
				model.SubdomainsDictionary[0].Elements.Add(hexa8element);
			}

			return numberOfNodes;
		}

		internal static int DefineContinuumElement3DNonLinear(Model model, IIsotropicContinuumMaterial3D material3D, IDynamicMaterial dynamicMaterial)
		{
			// Node creation
			var node1 = new Node(id: 1, x: 0.0, y: 0.0, z: 0.0);
			var node2 = new Node(id: 2, x: 100.0, y: 0.0, z: 0.0);
			var node3 = new Node(id: 3, x: 0.0, y: 100.0, z: 0.0);
			var node4 = new Node(id: 4, x: 100.0, y: 100.0, z: 0.0);
			var node5 = new Node(id: 5, x: 0.0, y: 0.0, z: 100.0);
			var node6 = new Node(id: 6, x: 100.0, y: 0.0, z: 100.0);
			var node7 = new Node(id: 7, x: 0.0, y: 100.0, z: 100.0);
			var node8 = new Node(id: 8, x: 100.0, y: 100.0, z: 100.0);

			// Create List of nodes
			IList<Node> nodes = new List<Node>();
			nodes.Add(node1);
			nodes.Add(node2);
			nodes.Add(node3);
			nodes.Add(node4);
			nodes.Add(node5);
			nodes.Add(node6);
			nodes.Add(node7);
			nodes.Add(node8);

			numberOfNodes = nodes.Count;

			// Add nodes to the nodes dictonary of the model
			for (int i = 0; i < numberOfNodes; ++i)
			{
				model.NodesDictionary.Add(i, nodes[i]);
			}

			int[][] connectivityMatrix = new int[8][];
			connectivityMatrix[0] = new int[] { 1, 1, 2, 4, 3, 5, 6, 8, 7 };

			var factory = new ContinuumElement3DNonLinearFactory(material3D, dynamicMaterial);

			for (int i = 0; i < 1; i++)
			{
				List<Node> elementNodeSet = new List<Node>(8);
				for (int j = 1; j < 9; j++)
				{
					elementNodeSet.Add((Node)model.NodesDictionary[connectivityMatrix[i][j] - 1]);
				}

				var hexa8element = new Element()
				{
					ID = connectivityMatrix[i][0],
					ElementType = factory.CreateElement(CellType.Hexa8, elementNodeSet),
				};

				for (int j = 1; j < 9; j++)
				{
					hexa8element.AddNode(model.NodesDictionary[connectivityMatrix[i][j] - 1]);
				}

				model.ElementsDictionary.Add(hexa8element.ID, hexa8element);
				model.SubdomainsDictionary[0].Elements.Add(hexa8element);
			}

			return numberOfNodes;
		}

		internal static int DefineContinuumElement3DNonLinearDefGrad(Model model, IContinuumMaterial3DDefGrad material3D, IDynamicMaterial dynamicMaterial)
		{
			// Node creation
			var node1 = new Node(id: 1, x: 0.0, y: 0.0, z: 0.0);
			var node2 = new Node(id: 2, x: 100.0, y: 0.0, z: 0.0);
			var node3 = new Node(id: 3, x: 0.0, y: 100.0, z: 0.0);
			var node4 = new Node(id: 4, x: 100.0, y: 100.0, z: 0.0);
			var node5 = new Node(id: 5, x: 0.0, y: 0.0, z: 100.0);
			var node6 = new Node(id: 6, x: 100.0, y: 0.0, z: 100.0);
			var node7 = new Node(id: 7, x: 0.0, y: 100.0, z: 100.0);
			var node8 = new Node(id: 8, x: 100.0, y: 100.0, z: 100.0);

			// Create List of nodes
			IList<Node> nodes = new List<Node>();
			nodes.Add(node1);
			nodes.Add(node2);
			nodes.Add(node3);
			nodes.Add(node4);
			nodes.Add(node5);
			nodes.Add(node6);
			nodes.Add(node7);
			nodes.Add(node8);

			numberOfNodes = nodes.Count;

			// Add nodes to the nodes dictonary of the model
			for (int i = 0; i < numberOfNodes; ++i)
			{
				model.NodesDictionary.Add(i, nodes[i]);
			}

			int[][] connectivityMatrix = new int[8][];
			connectivityMatrix[0] = new int[] { 1, 1, 2, 4, 3, 5, 6, 8, 7 };

			var factory = new ContinuumElement3DNonLinearDefGradFactory(material3D, dynamicMaterial, lambdag);

			for (int i = 0; i < 1; i++)
			{
				List<Node> elementNodeSet = new List<Node>(8);
				for (int j = 1; j < 9; j++)
				{
					elementNodeSet.Add((Node)model.NodesDictionary[connectivityMatrix[i][j] - 1]);
				}

				var hexa8element = new Element()
				{
					ID = connectivityMatrix[i][0],
					ElementType = factory.CreateElement(CellType.Hexa8, elementNodeSet),
				};

				for (int j = 1; j < 9; j++)
				{
					hexa8element.AddNode(model.NodesDictionary[connectivityMatrix[i][j] - 1]);
				}

				model.ElementsDictionary.Add(hexa8element.ID, hexa8element);
				model.SubdomainsDictionary[0].Elements.Add(hexa8element);
			}

			return numberOfNodes;
		}
		internal static void BoundaryConditions(Model model)
		{
			// Boundary Condtitions - Rolled Left-Face
			for (int iNode = 1; iNode <= 4; iNode++)
			{
				model.NodesDictionary[iNode - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			}

			// Boundary Condtitions - Rolled Bottom-Face
			model.NodesDictionary[1 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			model.NodesDictionary[2 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			model.NodesDictionary[5 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			model.NodesDictionary[6 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });

			// Boundary Condtitions - Rolled Back-Face
			model.NodesDictionary[1 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			model.NodesDictionary[3 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			model.NodesDictionary[5 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			model.NodesDictionary[7 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });

			// Boundary Condtitions - Prescribed Right-Face
			for (int iNode = 5; iNode <= 8; iNode++)
			{
				model.Loads.Add(new Load() { Node = model.NodesDictionary[iNode - 1], DOF = StructuralDof.TranslationZ, Amount = +50.0 });
			}
		}

		internal static void BoundaryConditionsNLM(Model model)
		{
			// Boundary Condtitions - Rolled Left-Face
			for (int iNode = 1; iNode <= 4; iNode++)
			{
				model.NodesDictionary[iNode - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
			}

			// Boundary Condtitions - Rolled Bottom-Face
			model.NodesDictionary[1 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			model.NodesDictionary[2 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			model.NodesDictionary[5 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
			model.NodesDictionary[6 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });

			// Boundary Condtitions - Rolled Back-Face
			model.NodesDictionary[1 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			model.NodesDictionary[3 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			model.NodesDictionary[5 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
			model.NodesDictionary[7 - 1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });

			// Boundary Condtitions - Prescribed Right-Face
			for (int iNode = 5; iNode <= 8; iNode++)
			{
				model.Loads.Add(new Load() { Node=model.NodesDictionary[iNode - 1], DOF=StructuralDof.TranslationZ, Amount= +200.0 });
			}
		}

	}
}
