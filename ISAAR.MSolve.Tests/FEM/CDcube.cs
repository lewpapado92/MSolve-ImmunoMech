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

namespace ISAAR.MSolve.Tests.FEM
{
    public class CDcube
    {
        private const int subdomainID = 0;
        [Fact]
        private static void RunTest()
        {

            Model model = CreateModel(1, new double[] { 2, 2, 2 }, 0).Item1;
            ComsolMeshReader2 modelReader = CreateModel(1, new double[] { 2, 2, 2 }, 0).Item2;
            string path0 = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput";
            string path3 = Path.Combine(Directory.GetCurrentDirectory(), "CDcubeOutput.vtu"); 
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
                    outputFile.WriteLine($"{i+1}");
                outputFile.WriteLine("              </DataArray>");

                outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"solution\" NumberOfComponents=\"1\" format=\"ascii\">");
                for (int i = 0; i < numberOfPoints; i++)
                    outputFile.WriteLine($"{solution[i]} ");
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
                {   if (model.Elements[i].Nodes.Count == 4)
                        outputFile.WriteLine($"{10} ");
                    else outputFile.WriteLine($"{5} ");
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

        private static Tuple<Model, ComsolMeshReader2> CreateModel(double k, double[] U, double L)
        {
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "fullModel.mphtxt");
            ComsolMeshReader2 modelReader = new ComsolMeshReader2(filename, new double[] { k }, new double[][] { U }, new double[] { L });
            Model model = modelReader.CreateModelFromFile();
            //Boundary Conditions
            var bodyLoad = new ConvectionDiffusionDomainLoad(new ConvectionDiffusionMaterial(k, U, L), 1, ThermalDof.Temperature);
            var flux1 = new FluxLoad(1);
            var flux2 = new FluxLoad(10);
            var dir1 = new DirichletDistribution(list => {
                return Vector.CreateWithValue(list.Count, 1);
            });
            var dir2 = new DirichletDistribution(list => {
                return Vector.CreateWithValue(list.Count, 0);
            });
            var weakDirichlet1 = new WeakDirichlet(dir1, k);
            var weakDirichlet2 = new WeakDirichlet(dir2, k);

            var bodyLoadElementFactory = new BodyLoadElementFactory(bodyLoad, model);
            var dirichletFactory1 = new SurfaceLoadElementFactory(weakDirichlet1);
            var dirichletFactory2 = new SurfaceLoadElementFactory(weakDirichlet2);
            var fluxFactory1 = new SurfaceLoadElementFactory(flux1);
            var fluxFactory2 = new SurfaceLoadElementFactory(flux2);
            var boundaryFactory3D = new SurfaceBoundaryFactory3D(0,
                new ConvectionDiffusionMaterial(k, new double[] { 0, 0, 0 }, 0));

            int[] domainIDs = new int[] { 0, };
            foreach (int domainID in domainIDs)
            {
                foreach (Element element in modelReader.elementDomains[domainID])
                {
                    //IReadOnlyList<Node> nodes = (IReadOnlyList<Node>)element.Nodes;
                    //var fluxElement1 = fluxFactory1.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(fluxElement1);
                    //var bodyLoadElementCellType = element.ElementType.CellType;
                    //var nodes = (IReadOnlyList<Node>)element.Nodes;
                    //var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Hexa8, nodes);
                    //model.BodyLoads.Add(bodyLoadElement);
                    // var surfaceElement = new SurfaceLoadElement();
                    //element.ID = TriID;
                    //surfaceElement.ElementType = DirichletElement1;
                    //model.SubdomainsDictionary[0].Elements.Add(dirichletElement1);
                    //model.ElementsDictionary.Add(TriID, surfaceElement);

                    //model.NodesDictionary[surfaceElement.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100 });
                }
            }

            //foreach (Node node in model.Nodes)
            //{
            //    model.Loads.Add(new Load() { Amount = .25, Node = node, DOF = ThermalDof.Temperature });
            //}

            int[] boundaryIDs = new int[] { 0,1,2,3,4,7 };
            int TriID = model.ElementsDictionary.Count + 1;
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.triBoundaries[boundaryID])
                {
                    //IReadOnlyList<Node> nodes = (IReadOnlyList<Node>)element.Nodes;
                    //var fluxElement1 = fluxFactory1.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(fluxElement1);
                    var dirichletElement1 = dirichletFactory1.CreateElement(CellType.Tri3, nodes);
                    model.SurfaceLoads.Add(dirichletElement1);
                    var surfaceBoundaryElement = boundaryFactory3D.CreateElement(CellType.Tri3, nodes);
                    var element = new Element();
                    element.ID = TriID;
                    element.ElementType = surfaceBoundaryElement;
                    model.SubdomainsDictionary[0].Elements.Add(element);
                    model.ElementsDictionary.Add(TriID, element);
                    foreach (Node node in nodes)
                    {
                        element.AddNode(node);
                    }
                    TriID += 1;
                    // var surfaceElement = new SurfaceLoadElement();
                    //element.ID = TriID;
                    //surfaceElement.ElementType = DirichletElement1;
                    //model.SubdomainsDictionary[0].Elements.Add(dirichletElement1);
                    //model.ElementsDictionary.Add(TriID, surfaceElement);

                    //model.NodesDictionary[surfaceElement.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100 });
                }
            }
            boundaryIDs = new int[] { 6, 8, 9 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IReadOnlyList<Node> nodes in modelReader.triBoundaries[boundaryID])
                {
                    //IReadOnlyList<Node> nodes = (IReadOnlyList<Node>)element.Nodes;
                    //var fluxElement2 = fluxFactory2.CreateElement(CellType.Quad4, nodes);
                    //model.SurfaceLoads.Add(fluxElement2);
                    var dirichletElement2 = dirichletFactory2.CreateElement(CellType.Tri3, nodes);
                    model.SurfaceLoads.Add(dirichletElement2);
                    var surfaceBoundaryElement = boundaryFactory3D.CreateElement(CellType.Tri3, nodes);
                    var element = new Element();
                    element.ID = TriID;
                    element.ElementType = surfaceBoundaryElement;
                    model.SubdomainsDictionary[0].Elements.Add(element);
                    model.ElementsDictionary.Add(TriID, element);
                    foreach (Node node in nodes)
                    {
                        element.AddNode(node);
                    }
                    TriID += 1;
                    // var surfaceElement = new SurfaceLoadElement();
                    //element.ID = TriID;
                    //surfaceElement.ElementType = DirichletElement1;
                    //model.SubdomainsDictionary[0].Elements.Add(dirichletElement1);
                    //model.ElementsDictionary.Add(TriID, surfaceElement);

                    //model.NodesDictionary[surfaceElement.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100 });
                }
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

            return new Tuple<Model, ComsolMeshReader2>(model, modelReader);
        }

        private static IVectorView SolveModel(Model model, ComsolMeshReader2 modelReader)
        {
            //var initialTemp = Vector.CreateZero(model.Nodes.Count);
            double[] temp0 = new double[model.Nodes.Count];
            int[] boundaryIDs = new int[] { 0, 1, 2, 3, 4, 7 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IList<Node> nodes in modelReader.triBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
                    {
                        temp0[node.ID] = 1;
                    }
                }
            }
            boundaryIDs = new int[] { 6, 8, 9 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (IList<Node> nodes in modelReader.triBoundaries[boundaryID])
                {
                    foreach (Node node in nodes)
                    {
                        temp0[node.ID] = 0;
                    }
                }
            }
            Vector initialTemp = Vector.CreateFromArray(temp0);
            var builder = new DenseMatrixSolver.Builder();
            builder.IsMatrixPositiveDefinite = false;
            var solver = builder.BuildSolver(model);
            var provider = new ProblemConvectionDiffusion2(model, solver);

            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var BDForder = 1;
            var parentAnalyzer = new BDF(model, solver, provider, childAnalyzer, .01, 2, BDForder, initialTemp);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return parentAnalyzer.temperature[BDForder][subdomainID];
            //           return solver.LinearSystems[subdomainID].Solution;
        }
    }
}

