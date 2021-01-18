using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;
using ISAAR.MSolve.Discretization.Interfaces;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Interfaces;
using System.Linq;
using ISAAR.MSolve.Solvers;
using System.IO;
using ISAAR.MSolve.FEM.Readers;
using ISAAR.MSolve.FEM.Readers.Interfaces;
using ISAAR.MSolve.FEM.Loading.SurfaceLoads;
using static ISAAR.MSolve.FEM.Loading.SurfaceLoads.WeakDirichlet;
using ISAAR.MSolve.FEM.Loading;
using ISAAR.MSolve.FEM.Elements.BoundaryConditionElements;
using System;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISSAR.MSolve.Discretization.Loads;
using ISAAR.MSolve.FEM.Loading.BodyLoads;
using System.Reflection;
using ISAAR.MSolve.FEM.Elements;

namespace ISAAR.MSolve.Tests.FEM
{
    public class HyperElasticGrowthCube
    {
        private const int subdomainID = 0;
        private static double lambdag = 1;
        private const double a = 1;
        private const double b = 1;
        private static DenseMatrixSolver.Builder builder = new DenseMatrixSolver.Builder();
        private static SkylineSolver.Builder structuralBuilder = new SkylineSolver.Builder();
        private static Dictionary<int, IVector>[] Solutions;
        private static Dictionary<int, IVector> Accelerations;
        private static Dictionary<int, IVector> Velocities;
        private static Dictionary<int, IVector> Displacements;

        [Fact]
        private static void RunTest()
        {
            //var models = new[] { CreateModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item1, CreateModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item1 };
            //var modelReaders = new[] { CreateModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item2, CreateModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0).Item2 };

            //var modelTuples = new[] { CreateConvectionDiffusionModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0), CreateConvectionDiffusionModel(1, new double[] { 2, 2, 2 }, 0, 1, 0, 0, 0) };
            var modelTuple1 = CreateGrowthModel(0, new double[] { 0, 0, 0 }, 1, 0, 0, 1.5);
            //var modelTuple3 = CreateStructuralModel(10e4, 0, new DynamicMaterial(.001, 0, 0, true), 0, new double[] { 0, 0, 0 });
            var models = new[] { modelTuple1.Item1 };
            var modelReaders = new[] { modelTuple1.Item2 };
            //IVectorView[] solutions = SolveModels(models, modelReaders);
            IVectorView[] solutions = SolveModelsWithNewmark(models, modelReaders);

            #region paraview commands
            //string path3 = Path.Combine(Directory.GetCurrentDirectory(), "CD_Structural_9HexaOutput.vtu");
            //var numberOfPoints = models[0].Nodes.Count;
            //var numberOfCells = models[0].Elements.Count;
            //using (StreamWriter outputFile = new StreamWriter(path3))
            //{
            //    outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
            //    outputFile.WriteLine("  <UnstructuredGrid>");
            //    outputFile.WriteLine($"     <Piece NumberOfPoints=\"{numberOfPoints}\" NumberOfCells=\"{numberOfCells}\">");
            //    outputFile.WriteLine("          <Points>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"position\" NumberOfComponents=\"3\" format =\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{models[0].Nodes[i].X} {models[0].Nodes[i].Y} {models[0].Nodes[i].Z} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("          </Points>");
            //    outputFile.WriteLine("          <PointData>");

            //    outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"node_ID\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{i + 1}");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"solution1\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{Solutions[0][0][i]} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"solution2\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{Solutions[1][0][i]} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"displacement\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfPoints; i++)
            //        outputFile.WriteLine($"{Displacements[0][3 * i]} ");
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("          </PointData>");
            //    outputFile.WriteLine("          <CellData>");
            //    outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"element_ID\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfCells; i++)
            //    {
            //        outputFile.WriteLine($"{i + 1}");
            //    }
            //    outputFile.WriteLine("              </DataArray>");
            //    outputFile.WriteLine("          </CellData>");
            //    outputFile.WriteLine("          <Cells>");

            //    outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"connectivity\">");
            //    for (int i = 0; i < numberOfCells; i++)
            //    {
            //        for (int j = 0; j < models[0].Elements[i].Nodes.Count; j++)
            //            outputFile.Write($"{models[0].Elements[i].Nodes[j].ID} ");
            //        outputFile.WriteLine("");
            //    }
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">");
            //    var offset = 0;
            //    for (int i = 0; i < numberOfCells; i++)
            //    {
            //        offset += models[0].Elements[i].Nodes.Count;
            //        outputFile.WriteLine($"{offset} ");
            //    }
            //    outputFile.WriteLine("              </DataArray>");

            //    outputFile.WriteLine("              <DataArray type=\"Int32\" Name =\"types\" NumberOfComponents =\"1\" format=\"ascii\">");
            //    for (int i = 0; i < numberOfCells; i++)
            //    {
            //        if (models[0].Elements[i].Nodes.Count == 8)
            //            outputFile.WriteLine($"{12} ");
            //        else outputFile.WriteLine($"{9} ");
            //    }
            //    outputFile.WriteLine("              </DataArray>");
            //    outputFile.WriteLine("          </Cells>");
            //    outputFile.WriteLine("      </Piece>");
            //    outputFile.WriteLine("  </UnstructuredGrid>");
            //    outputFile.WriteLine("</VTKFile>");
            //}end*/
            #endregion
            Assert.True(CompareResults(solutions[0]));
        }

        private static void UpdateModels(Dictionary<int, IVector>[] solutions, IStructuralModel[] modelsToReplace, ISolver[] solversToReplace,
            IConvectionDiffusionIntegrationProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
        {
            Solutions = solutions;
            double[] lg = solutions[0][0].CopyToArray();
            //double[] sol1 = solutions[1][0].CopyToArray();
            lambdag = lg[0];// (sol0[39] + sol0[38] + sol0[37] + sol1[36])/8 ;
            modelsToReplace[0] = CreateGrowthModel(0, new double[] { 0, 0, 0 }, 1, 0, 0, 1.5).Item1;
            //modelsToReplace[1] = CreateConvectionDiffusionModel2(1, new double[] { 0, 0, 0 }, 0, 0, fluxLoad, new double[] { }).Item1;

            for (int i = 0; i < modelsToReplace.Length; i++)
            {
                solversToReplace[i] = builder.BuildSolver(modelsToReplace[i]);
                providersToReplace[i] = new ProblemConvectionDiffusion2((Model)modelsToReplace[i], solversToReplace[i]);
                childAnalyzersToReplace[i] = new LinearAnalyzer(modelsToReplace[i], solversToReplace[i], providersToReplace[i]);
            }
        }
        private static void ReplaceLambdaGInModel(IStructuralModel model, double lg)
        {
            foreach (var e in model.Elements)
            {
                var et = (ContinuumElement3DNonLinearDefGrad)e.ElementType;
                var bindFlags = BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Static;
                FieldInfo field = typeof(ContinuumElement3DNonLinearDefGrad).GetField("lambdag", bindFlags);
                field.SetValue(et, lg);
            }
        }

        private static void UpdateNewmarkModel(Dictionary<int, IVector> accelerations, Dictionary<int, IVector> velocities, Dictionary<int, IVector> displacements, IStructuralModel[] modelsToReplace,
            ISolver[] solversToReplace, IImplicitIntegrationProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
        {
            double[] sol0 = Solutions[0][0].CopyToArray();
            //double[] sol1 = Solutions[1][0].CopyToArray();
            //var load = new double[] { sol1[39], sol1[39], sol1[39] };
            //IDynamicMaterial commonDynamicMaterialProperties = new DynamicMaterial(.001, 0, 0, true);
            //Accelerations = accelerations;
            //Velocities = velocities;
            //Displacements = displacements;
            ////modelsToReplace[0] = CreateStructuralModel(3e4, 0, commonDynamicMaterialProperties, 0, new double[] { 1000, 0, 0 }, lambdag).Item1;
            ReplaceLambdaGInModel(modelsToReplace[0], lambdag);
            solversToReplace[0] = structuralBuilder.BuildSolver(modelsToReplace[0]);
            providersToReplace[0] = new ProblemStructural(modelsToReplace[0], solversToReplace[0]);
            solversToReplace[0].HandleMatrixWillBeSet();
            //childAnalyzersToReplace[0] = new LinearAnalyzer(modelsToReplace[0], solversToReplace[0], providersToReplace[0]);
            var increments = 2;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(modelsToReplace[0], solversToReplace[0], (INonLinearProvider)providersToReplace[0], increments);
            childAnalyzerBuilder.ResidualTolerance = 1E-6;
            childAnalyzerBuilder.MaxIterationsPerIncrement = 50;
            childAnalyzerBuilder.NumIterationsForMatrixRebuild = 1;
            childAnalyzersToReplace[0] = childAnalyzerBuilder.Build();
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
        private static Tuple<Model, IModelReader> CreateGrowthModel(double k, double[] U, double L, double b, double f, double bl)
        {
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "125HexaHyperelasticCube.mphtxt");
            var modelReader = new ComsolMeshReader2(filename, new double[] { k }, new double[][] { U }, new double[] { L });
            Model model = modelReader.CreateModelFromFile();
            var material = new ConvectionDiffusionMaterial(k, U, L);
            //Boundary Conditions
            //var flux1 = new FluxLoad(f);
            //var dir1 = new DirichletDistribution(list => {
            //    return Vector.CreateWithValue(list.Count, b);
            //});
            //var dir2 = new DirichletDistribution(list =>
            //{
            //    return Vector.CreateWithValue(list.Count, b2);
            //});
            //var weakDirichlet1 = new WeakDirichlet(dir1, k);
            //var weakDirichlet2 = new WeakDirichlet(dir2, k);

            //var dirichletFactory1 = new SurfaceLoadElementFactory(weakDirichlet1);
            //var dirichletFactory2 = new SurfaceLoadElementFactory(weakDirichlet2);
            //var fluxFactory1 = new SurfaceLoadElementFactory(flux1);
            //var boundaryFactory3D = new SurfaceBoundaryFactory3D(0,
            //    new ConvectionDiffusionMaterial(k, new double[] { 0, 0, 0 }, 0));


            int[] domainIDs = new int[] { 0, };
            foreach (int domainID in domainIDs)
            {
                foreach (Element element in modelReader.elementDomains[domainID])
                {
                    ////IReadOnlyList<Node> nodes = (IReadOnlyList<Node>)element.Nodes;
                    //var fluxElement1 = fluxFactory1.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(fluxElement1);
                    var bodyLoadElementCellType = element.ElementType.CellType;
                    var nodes = (IReadOnlyList<Node>)element.Nodes;
                    var domainLoad = new ConvectionDiffusionDomainLoad(material, bl, ThermalDof.Temperature);
                    var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
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
            //foreach (Node node in model.Nodes)
            //{
            //    model.Loads.Add(new Load() { Amount = bl[node.ID], Node = node, DOF = ThermalDof.Temperature });
            //}

            int[] boundaryIDs = new int[] { 0, };
            int QuadID = model.ElementsDictionary.Count + 1;

            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    //var fluxElement1 = fluxFactory1.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(fluxElement1);
                    //var dirichletElement1 = dirichletFactory1.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(dirichletElement1);
                    //var SurfaceBoundaryElement = boundaryFactory3D.CreateElement(CellType.Quad4, nodes);
                    //var element = new Element();
                    //element.ID = QuadID;
                    //element.ElementType = SurfaceBoundaryElement;
                    //model.SubdomainsDictionary[0].Elements.Add(element);
                    //model.ElementsDictionary.Add(QuadID, element);
                    //foreach (Node node in nodes)
                    //{
                    //    element.AddNode(node);
                    //}
                    //QuadID += 1;
                }
            }
            boundaryIDs = new int[] { 5 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    //var fluxElement1 = fluxFactory1.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(fluxElement1);
                    //var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(bodyLoadElement);
                    //var SurfaceBoundaryElement = boundaryFactory3D.CreateElement(CellType.Quad4, nodes);
                    //var element = new Element();
                    //element.ID = QuadID;
                    //element.ElementType = SurfaceBoundaryElement;
                    //model.SubdomainsDictionary[0].Elements.Add(element);
                    //model.ElementsDictionary.Add(QuadID, element);
                    //foreach (Node node in nodes)
                    //{
                    //    element.AddNode(node);
                    //}
                    //QuadID += 1;
                }
            }
            return new Tuple<Model, IModelReader>(model, modelReader);
        }
        private static Tuple<Model, IModelReader> CreateStructuralModel(double[] C1, double[] C2, IDynamicMaterial[] commonDynamicMaterialProperties, double b, double[] l, double lambdag)
        {
            double[] poissonV = new double[C1.Length];
            double[] muLame = new double[C1.Length];
            double[] bulkModulusnew = new double[C1.Length];
            for (int i = 0; i < C1.Length; i++)
            {
                poissonV[i] = 0.2;
                muLame[i] = 2 * C1[i];
                bulkModulusnew[i] = 2 * muLame[i] * (1 + poissonV[i]) / (3 * (1 - 2 * poissonV[i]));
            }
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "125HexaHyperelasticCube.mphtxt");
            var modelReader = new ComsolMeshReader1(filename, C1, C2, bulkModulusnew, commonDynamicMaterialProperties, new double[] { lambdag, 1 });
            Model model = modelReader.CreateModelFromFile();
            //Boundary Conditions
            var lx = l[0];
            var ly = l[1];
            var lz = l[2];
            var distributedLoad = new DistributedLoad(lx, ly, lz);
            //var flux2 = new FluxLoad(f2);
            //var dir1 = new DirichletDistribution(list => {
            //    return Vector.CreateWithValue(list.Count, b1);
            //});
            //var dir2 = new DirichletDistribution(list => {
            //    return Vector.CreateWithValue(list.Count, b2);
            //});
            //var weakDirichlet1 = new WeakDirichlet(dir1, k);
            //var weakDirichlet2 = new WeakDirichlet(dir2, k);

            //var dirichletFactory1 = new SurfaceLoadElementFactory(weakDirichlet1);
            //var dirichletFactory2 = new SurfaceLoadElementFactory(weakDirichlet2);
            var distributedLoadFactory = new SurfaceLoadElementFactory(distributedLoad);
            //var fluxFactory2 = new SurfaceLoadElementFactory(flux2);
            //var boundaryFactory3D = new SurfaceBoundaryFactory3D(0,
            //    new ConvectionDiffusionMaterial(k, new double[] { 0, 0, 0 }, 0));


            int[] boundaryIDs = new int[] { 0 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Node node in modelReader.nodeBoundaries[boundaryID])
                {
                        node.Constraints.Add(new Constraint()
                        {
                            Amount = b,
                            DOF = StructuralDof.TranslationX
                        });
                        //node.Constraints.Add(new Constraint()
                        //{
                        //    Amount = b,
                        //    DOF = StructuralDof.TranslationY
                        //});
                        //node.Constraints.Add(new Constraint()
                        //{
                        //    Amount = b,
                        //    DOF = StructuralDof.TranslationZ
                        //});
                }
            }
            boundaryIDs = new int[] { 1 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Node node in modelReader.nodeBoundaries[boundaryID])
                {
                    //node.Constraints.Add(new Constraint()
                    //{
                    //    Amount = b,
                    //    DOF = StructuralDof.TranslationX
                    //});
                    node.Constraints.Add(new Constraint()
                        {
                            Amount = b,
                            DOF = StructuralDof.TranslationY
                        });
                        //node.Constraints.Add(new Constraint()
                        //{
                        //    Amount = b,
                        //    DOF = StructuralDof.TranslationZ
                        //});
                }
            }
            boundaryIDs = new int[] { 2 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Node node in modelReader.nodeBoundaries[boundaryID])
                {
                    //node.Constraints.Add(new Constraint()
                    //{
                    //    Amount = b,
                    //    DOF = StructuralDof.TranslationX
                    //});
                    //node.Constraints.Add(new Constraint()
                    //{
                    //    Amount = b,
                    //    DOF = StructuralDof.TranslationY
                    //});
                    node.Constraints.Add(new Constraint()
                        {
                            Amount = b,
                            DOF = StructuralDof.TranslationZ
                        });
                }
            }
            boundaryIDs = new int[] { 3 };
            int QuadID = model.ElementsDictionary.Count + 1;
            //foreach (int boundaryID in boundaryIDs)
            //{
            //    foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
            //    {
            //        //var distributedLoadElement = distributedLoadFactory.CreateElement(CellType.Quad4, nodes);
            //        //model.SurfaceLoads.Add(distributedLoadElement);
            //        //var dirichletElement2 = dirichletFactory2.CreateElement(CellType.Quad4, nodes);
            //        //model.SurfaceLoads.Add(dirichletElement2);
            //        //var SurfaceBoundaryElement = boundaryFactory3D.CreateElement(CellType.Quad4, nodes);
            //        //var element = new Element();
            //        //element.ID = QuadID;
            //        //element.ElementType = SurfaceBoundaryElement;
            //        //model.SubdomainsDictionary[0].Elements.Add(element);
            //        //model.ElementsDictionary.Add(QuadID, element);
            //        //foreach (Node node in nodes)
            //        //{
            //        //    element.AddNode(node);
            //        //}
            //        //QuadID += 1;
            //        foreach (Node node in nodes)
            //        {
            //            model.Loads.Add(new Load() { Node = node, DOF = StructuralDof.TranslationZ, Amount = +100.0 });
            //        }
            //    }
            //}
            //int[] nodeIDs = new int[] { 1, 2, 3, 4 };
            int[] nodeIDs = new int[] { 0, 5, 10, 40 };
            foreach (int nodeID in nodeIDs)
            {
                model.Loads.Add(new Load() { Node = model.NodesDictionary[nodeID], DOF = StructuralDof.TranslationZ, Amount = +lz });
            }
            return new Tuple<Model, IModelReader>(model, modelReader);
        }
        private static IVectorView[] SolveModelsWithNewmark(Model[] models, IModelReader[] modelReaders)
        {
            Vector[] initialValues = new Vector[models.Length];
            var value0 = new Dictionary<int, double[]>();
            for (int i = 0; i < models.Length; i++)
            {
                double[] v0 = new double[models[i].Nodes.Count];
                value0.Add(i, v0);
            }
            foreach (Node node in models[0].Nodes)
            {
                value0[0][node.ID] = 1;
            }

            DenseMatrixSolver[] solvers = new DenseMatrixSolver[models.Length];
            IConvectionDiffusionIntegrationProvider[] providers = new IConvectionDiffusionIntegrationProvider[models.Length];
            IChildAnalyzer[] childAnalyzers = new IChildAnalyzer[models.Length];
            for (int i = 0; i < models.Length; i++)
            {
                initialValues[i] = Vector.CreateFromArray(value0[i]);
                //var builder = new DenseMatrixSolver.Builder();
                builder.IsMatrixPositiveDefinite = false;
                solvers[i] = builder.BuildSolver(models[i]);
                providers[i] = new ProblemConvectionDiffusion2(models[i], solvers[i]);
                childAnalyzers[i] = new LinearAnalyzer(models[i], solvers[i], providers[i]);
            }

            const double timestep = .1;
            const double time = 3;
            var parentAnalyzer = new ConvectionDiffusionImplicitDynamicAnalyzerMultiModel(UpdateModels, models, solvers,
                providers, childAnalyzers, timestep, time, initialTemperature: initialValues);
            parentAnalyzer.Initialize();

            //var structuralModel = CreateStructuralModel(10.5e3, 0, new DynamicMaterial(.001, 0, 0, true), 0, new double[] { 0, 0, 25000 }, lambdag).Item1; // new Model();
            DynamicMaterial[] dynamicMaterials = new DynamicMaterial[] { new DynamicMaterial(1000, 0, 0, true)};
            var structuralModel = CreateStructuralModel(new double[] { 1 }, new double[] { 1 }, dynamicMaterials, 0, new double[] { 0, 0, 2500 }, lambdag).Item1; // new Model();
            var structuralSolver = structuralBuilder.BuildSolver(structuralModel);
            var structuralProvider = new ProblemStructural(structuralModel, structuralSolver);
            //var structuralChildAnalyzer = new LinearAnalyzer(structuralModel, structuralSolver, structuralProvider);
            var increments = 2;
            var structuralChildAnalyzerBuilder = new LoadControlAnalyzer.Builder(structuralModel, structuralSolver, structuralProvider, increments);
            structuralChildAnalyzerBuilder.ResidualTolerance = 1E-6;
            structuralChildAnalyzerBuilder.MaxIterationsPerIncrement = 50;
            structuralChildAnalyzerBuilder.NumIterationsForMatrixRebuild = 1;
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer structuralChildAnalyzer = structuralChildAnalyzerBuilder.Build();
            var structuralParentAnalyzer = new NewmarkDynamicAnalyzer(UpdateNewmarkModel, structuralModel, structuralSolver,
                structuralProvider, structuralChildAnalyzer, timestep, time, 0.25, 0.5);
            structuralParentAnalyzer.Initialize();

            for (int i = 0; i < time / timestep; i++)
            {
                parentAnalyzer.SolveTimestep(i);
                structuralParentAnalyzer.SolveTimestep(i);
            }

            return solvers.Select(x => x.LinearSystems[subdomainID].Solution).ToArray();
        }
    }
}