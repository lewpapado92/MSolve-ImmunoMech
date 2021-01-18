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
using ISSAR.MSolve.Discretization.Loads;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using MathNet.Numerics.LinearAlgebra;
using System.Text;
using ISAAR.MSolve.FEM.Readers.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using MathNet.Numerics;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using System.Reflection;

namespace ISAAR.MSolve.Tests.FEM
{
    public class HyperElasticStrechCube
    {
        private const int subdomainID = 0;
        private static double lambdag = 1;
        private static IVector Displacements;
        private static SkylineSolver.Builder builder = new SkylineSolver.Builder();

        [Fact]
        private static void RunTest()
        {

            Model model = CreateModel1(1, 1, new DynamicMaterial[] { new DynamicMaterial(.001, 0, 0, true) }, 0, new double[] { 0, 0, 250 }, lambdag).Item1; ;
            IModelReader modelReader = CreateModel1(1, 1, new DynamicMaterial[] { new DynamicMaterial(.001, 0, 0, true) }, 0, new double[] { 0, 0, 250 }, lambdag).Item2;
            string path0 = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput";
            string path3 = Path.Combine(Directory.GetCurrentDirectory(), "HyperElastiGrowthCantilever.vtu");
            //var path2 = Path.Combine(path0, $"nodes.txt");

            IVectorView solution = SolveModel(model, modelReader);

            var numberOfPoints = model.Nodes.Count;
            var numberOfCells = model.Elements.Count;
            using (StreamWriter outputFile = new StreamWriter(path3))
            {
                outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
                outputFile.WriteLine("  <UnstructuredGrid>");
                outputFile.WriteLine($"     <Piece NumberOfPoints=\"{numberOfPoints}\" NumberOfCells=\"{numberOfCells}\">");
                outputFile.WriteLine("          <Points>");

                outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"position\" NumberOfComponents=\"3\" format =\"ascii\">");
                for (int i = 0; i < numberOfPoints; i++)
                    outputFile.WriteLine($"{model.Nodes[i].X} {model.Nodes[i].Y} {model.Nodes[i].Z} ");
                outputFile.WriteLine("              </DataArray>");

                outputFile.WriteLine("          </Points>");
                outputFile.WriteLine("          <PointData>");

                outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"node_ID\" NumberOfComponents=\"1\" format=\"ascii\">");
                for (int i = 0; i < numberOfPoints; i++)
                    outputFile.WriteLine($"{i + 1}");
                outputFile.WriteLine("              </DataArray>");

                outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"displacement\" NumberOfComponents=\"1\" format=\"ascii\">");
                for (int i = 0; i < 4; i++)
                    outputFile.WriteLine($"{0} ");
                for (int i = 0; i < numberOfPoints - 4; i++)
                    outputFile.WriteLine($"{Displacements[3 * (i + 1) - 1]} ");
                outputFile.WriteLine("              </DataArray>");
                outputFile.WriteLine("          </PointData>");
                outputFile.WriteLine("          <CellData>");
                outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"element_ID\" NumberOfComponents=\"1\" format=\"ascii\">");
                for (int i = 0; i < numberOfCells; i++)
                {
                    outputFile.WriteLine($"{i + 1}");
                }
                outputFile.WriteLine("              </DataArray>");
                outputFile.WriteLine("          </CellData>");
                outputFile.WriteLine("          <Cells>");

                outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"connectivity\">");
                for (int i = 0; i < numberOfCells; i++)
                {
                    for (int j = 0; j < model.Elements[i].Nodes.Count; j++)
                        outputFile.Write($"{model.Elements[i].Nodes[j].ID} ");
                    outputFile.WriteLine("");
                }
                outputFile.WriteLine("              </DataArray>");

                outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">");
                var offset = 0;
                for (int i = 0; i < numberOfCells; i++)
                {
                    offset += model.Elements[i].Nodes.Count;
                    outputFile.WriteLine($"{offset} ");
                }
                outputFile.WriteLine("              </DataArray>");

                outputFile.WriteLine("              <DataArray type=\"Int32\" Name =\"types\" NumberOfComponents =\"1\" format=\"ascii\">");
                for (int i = 0; i < numberOfCells; i++)
                {
                    if (model.Elements[i].Nodes.Count == 8)
                        outputFile.WriteLine($"{12} ");
                    else outputFile.WriteLine($"{9} ");
                }
                outputFile.WriteLine("              </DataArray>");
                outputFile.WriteLine("          </Cells>");
                outputFile.WriteLine("      </Piece>");
                outputFile.WriteLine("  </UnstructuredGrid>");
                outputFile.WriteLine("</VTKFile>");
            }

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
            //double[] disp = displacements[0].CopyToArray();
            //Displacements = Vector.CreateFromArray(disp);
            IDynamicMaterial commonDynamicMaterialProperties = new DynamicMaterial(.001, 0, 0, true);
            //modelsToReplace[0] = CreateModel1(1, 1, new DynamicMaterial(.001, 0, 0, true), 0, new double[] { 0, 0, 200 }, lambdag).Item1;
            ReplaceLambdaGInModel(modelsToReplace[0], lambdag);
            solversToReplace[0] = builder.BuildSolver(modelsToReplace[0]);
            providersToReplace[0] = new ProblemStructural(modelsToReplace[0], solversToReplace[0]);
            var increments = 2;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(modelsToReplace[0], solversToReplace[0], (INonLinearProvider)providersToReplace[0], increments);
            childAnalyzersToReplace[0] = childAnalyzerBuilder.Build();
        }
        private static Tuple<Model, IModelReader> CreateModel1(double C1, double C2, IDynamicMaterial[] commonDynamicMaterialProperties, double b, double[] l, double lambdag)
        {
            double poissonV = 0.2;
            double muLame = 2 * C1;
            double bulkModulus = 2 * muLame * (1 + poissonV) / (3 * (1 - 2 * poissonV));
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "125HexaHyperelasticCube.mphtxt");
            var modelReader = new ComsolMeshReader1(filename, new double[] { C1 }, new double[] { C2 }, new double[] { 1 }, commonDynamicMaterialProperties, new double[] { lambdag });
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
                foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
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
            }
            boundaryIDs = new int[] { 1 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
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
            }
            boundaryIDs = new int[] { 2 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
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
            int[] nodeIDs = new int[] { 0, 5, 10, 40 };
            foreach (int nodeID in nodeIDs)
            {
                model.Loads.Add(new Load() { Node = model.NodesDictionary[nodeID], DOF = StructuralDof.TranslationZ, Amount = +lz });
            }
            return new Tuple<Model, IModelReader>(model, modelReader);
        }
        private static Tuple<Model, IModelReader> CreateModel2(double C1, double C2, IDynamicMaterial[] commonDynamicMaterialProperties, double b, double[] l, double lambdag)
        {
            double poissonV = 0.2;
            double muLame = 2 * C1;
            double bulkModulus = 2 * muLame * (1 + poissonV) / (3 * (1 - 2 * poissonV));
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "9hexa.mphtxt");
            var modelReader = new ComsolMeshReader1(filename, new double[] { C1 }, new double[] { C2 }, new double[] { 1 }, commonDynamicMaterialProperties, new double[] { lambdag });
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


            int[] boundaryIDs = new int[] { 0, 5 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
                    {
                        node.Constraints.Add(new Constraint()
                        {
                            Amount = b,
                            DOF = StructuralDof.TranslationX
                        });
                        node.Constraints.Add(new Constraint()
                        {
                            Amount = b,
                            DOF = StructuralDof.TranslationY
                        });
                        node.Constraints.Add(new Constraint()
                        {
                            Amount = b,
                            DOF = StructuralDof.TranslationZ
                        });
                    }
                }
            }
            return new Tuple<Model, IModelReader>(model, modelReader);
        }
        private static Tuple<Model, IModelReader> CreateModel3(double C1, double C2, IDynamicMaterial[] commonDynamicMaterialProperties, double b, double[] l, double lambdag)
        {
            double poissonV = 0.2;
            double muLame = 2 * C1;
            double bulkModulus = 2 * muLame * (1 + poissonV) / (3 * (1 - 2 * poissonV));
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "9hexa.mphtxt");
            var modelReader = new ComsolMeshReader1(filename, new double[] { C1 }, new double[] { C2 }, new double[] { 1 }, commonDynamicMaterialProperties, new double[] { lambdag });
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


            int[] boundaryIDs = new int[] { 0, 5 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.quadBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
                    {
                        node.Constraints.Add(new Constraint()
                        {
                            Amount = b,
                            DOF = StructuralDof.TranslationX
                        });
                        node.Constraints.Add(new Constraint()
                        {
                            Amount = b,
                            DOF = StructuralDof.TranslationY
                        });
                        node.Constraints.Add(new Constraint()
                        {
                            Amount = b,
                            DOF = StructuralDof.TranslationZ
                        });
                    }
                }
            }
            return new Tuple<Model, IModelReader>(model, modelReader);
        }

        private static IVectorView SolveModel(Model model, IModelReader modelReader)
        {
            var builder = new SkylineSolver.Builder();
            //builder.IsMatrixPositiveDefinite = false;
            var solver = builder.BuildSolver(model);
            const double timestep = 1;
            const double time = 100;

            var provider = new ProblemStructural(model, solver);
            var increments = 2;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            childAnalyzerBuilder.ResidualTolerance = 1E-6;
            childAnalyzerBuilder.MaxIterationsPerIncrement = 50;
            childAnalyzerBuilder.NumIterationsForMatrixRebuild = 1;
            LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();

            //var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, timestep, time);
            //parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration();
            //NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();
            var parentAnalyzer = new NewmarkDynamicAnalyzer(UpdateNewmarkModel, model, solver,
                provider, childAnalyzer, timestep, time, .25, .5);


            parentAnalyzer.Initialize();
            for (int i = 0; i < time / timestep; i++)
            {
                //lambdag = .01 * i + 1;
                parentAnalyzer.SolveTimestep(i);
            }

            return solver.LinearSystems[subdomainID].Solution;
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
                model.NodesDictionary.Add(nNode + 1, new Node(id: nNode + 1, x: nodeData[nNode, 0], y: nodeData[nNode, 1], z: nodeData[nNode, 2]));

            }

            // orismos elements 
            Element e1;
            int subdomainID = 1;
            for (int nElement = 0; nElement < elementData.GetLength(0) - 1; nElement++)
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
                model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
            }

            // constraint to to deftero miso apo th list twn nodes
            foreach (int k in new int[] { 5, 6, 7, 8 })
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

