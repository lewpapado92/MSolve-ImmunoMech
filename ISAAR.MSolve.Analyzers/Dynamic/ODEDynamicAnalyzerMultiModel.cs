using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: Optimization: I could avoid initialization and GC of some vectors by reusing existing ones.
//TODO: Use a base class for implicit time integration methods (perhaps to together with explicit)
namespace ISAAR.MSolve.Analyzers.Dynamic
{
    public class ODEDynamicAnalyzerMultiModel : INonLinearParentAnalyzer //TODO: why is this non linear
    {
        private readonly int maxStaggeredSteps;
        private readonly double beta, timeStep, totalTime, tolerance;
        private IStructuralModel[] models;
        private readonly IReadOnlyDictionary<int, ILinearSystem>[] linearSystems;
        private readonly ISolver[] solvers;
        private readonly IChildAnalyzer[] childAnalyzers;
        private readonly IImplicitIntegrationProvider[] providers;
        private readonly IVector[] initialValue;
        private readonly Dictionary<int, IVector>[] rhs;
        private readonly Dictionary<int, IVector>[] rhsPrevious;//TODO: at the moment domain loads are not implemented in this
        private readonly Dictionary<int, IVector>[] unknown;
        private readonly Dictionary<int, IVector>[] unknownFromPreviousStaggeredStep;
        private readonly Dictionary<int, IVector>[] firstOrderMatrixTimesUnknown;
        private readonly Dictionary<int, IVector>[] secondOrderMatrixTimesUnknown;
        private readonly Action<Dictionary<int, IVector>[], IStructuralModel[], ISolver[], IImplicitIntegrationProvider[], IChildAnalyzer[]> CreateNewModel;

        public ODEDynamicAnalyzerMultiModel(Action<Dictionary<int, IVector>[], IStructuralModel[], ISolver[], IImplicitIntegrationProvider[], IChildAnalyzer[]> modelCreator,
            IStructuralModel[] models, ISolver[] solvers, IImplicitIntegrationProvider[] providers, IChildAnalyzer[] childAnalyzers, double timeStep, double totalTime,
            int maxStaggeredSteps = 100, double tolerance = 1e-3, IVector[] initialValue = null)
        {
            this.CreateNewModel = modelCreator;
            this.maxStaggeredSteps = maxStaggeredSteps;
            this.tolerance = tolerance;
            this.models = models;
            this.linearSystems = new IReadOnlyDictionary<int, ILinearSystem>[solvers.Length];
            rhs = new Dictionary<int, IVector>[solvers.Length];
            rhsPrevious = new Dictionary<int, IVector>[solvers.Length];
            unknown = new Dictionary<int, IVector>[solvers.Length];
            unknownFromPreviousStaggeredStep = new Dictionary<int, IVector>[solvers.Length];
            firstOrderMatrixTimesUnknown = new Dictionary<int, IVector>[solvers.Length];
            secondOrderMatrixTimesUnknown = new Dictionary<int, IVector>[solvers.Length];
            for (int i = 0; i < solvers.Length; i++)
            {
                this.linearSystems[i] = solvers[i].LinearSystems;
                rhs[i] = new Dictionary<int, IVector>();
                rhsPrevious[i] = new Dictionary<int, IVector>();//TODO: at the moment domain loads are not implemented in this
                unknown[i] = new Dictionary<int, IVector>();
                unknownFromPreviousStaggeredStep[i] = new Dictionary<int, IVector>();
                firstOrderMatrixTimesUnknown[i] = new Dictionary<int, IVector>();
                secondOrderMatrixTimesUnknown[i] = new Dictionary<int, IVector>();
            }
            this.solvers = solvers;
            //solver.PreventFromOverwrittingSystemMatrices(); //TODO: If the scheme is purely implicit we can overwrite the matrix.
            this.providers = providers;
            this.ChildAnalyzer = childAnalyzers[0];
            this.childAnalyzers = childAnalyzers;
            this.timeStep = timeStep;
            this.totalTime = totalTime;
            this.ChildAnalyzer.ParentAnalyzer = this;
            this.initialValue = initialValue;
        }

        public Dictionary<int, IAnalyzerLog[]> Logs => null; //TODO: this can't be right
        public Dictionary<int, ImplicitIntegrationAnalyzerLog> ResultStorages { get; }
            = new Dictionary<int, ImplicitIntegrationAnalyzerLog>();

        public IChildAnalyzer ChildAnalyzer { get; }

        public void BuildMatrices()
        {
            var coeffs = new ImplicitIntegrationCoefficients
            {
                Mass = 1 / timeStep,
                Stiffness = beta
            };
            for (int i = 0; i < linearSystems.Length; i++)
            {
                foreach (ILinearSystem linearSystem in linearSystems[i].Values)
                {
                    linearSystem.Matrix = providers[i].LinearCombinationOfMatricesIntoStiffness(coeffs, linearSystem.Subdomain);
                }
            }
        }

        public IVector GetOtherRhsComponents(ILinearSystem linearSystem, IVector currentSolution)
        {
            #region old code
            //// u[id]: old solution
            //// v[id]: current solution
            //// vv: old acceleration
            //// v2: current acceleration
            //// v1: current velocity
            ////double vv = v2[id].Data[j];
            ////v2[id].Data[j] = a0 * (v[id].Data[j] - u[id].Data[j]) - a2 * v1[id].Data[j] - a3 * vv;
            ////v1[id].Data[j] += a6 * vv + a7 * v2[id].Data[j];

            //int id = subdomain.ID;
            //Vector<double> currentAcceleration = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> currentVelocity = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> uu = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> uc = new Vector<double>(subdomain.Solution.Length);
            //for (int j = 0; j < subdomain.Rhs.Length; j++)
            //{
            //    currentAcceleration.Data[j] = a0 * (currentSolution[j] - v[id].Data[j]) - a2 * v1[id].Data[j] - a3 * v2[id].Data[j];
            //    currentVelocity.Data[j] = v1[id].Data[j] + a6 * v2[id].Data[j] + a7 * currentAcceleration.Data[j];
            //    uu.Data[j] = a0 * currentSolution[j] + a2 * currentVelocity.Data[j] + a3 * currentAcceleration.Data[j];
            //    uc.Data[j] = a1 * currentSolution[j] + a4 * currentVelocity.Data[j] + a5 * currentAcceleration.Data[j];
            //}

            //Vector<double> tempResult = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> result = new Vector<double>(subdomain.Solution.Length);
            //provider.MassMatrixVectorProduct(subdomain, uu, tempResult.Data);
            //result.Add(tempResult);

            //provider.DampingMatrixVectorProduct(subdomain, uc, tempResult.Data);
            //result.Add(tempResult);

            //return result.Data;

            //Vector<double> uu = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> uc = new Vector<double>(subdomain.Solution.Length);
            //int id = subdomain.ID;
            //for (int j = 0; j < subdomain.Rhs.Length; j++)
            //{
            //    uu.Data[j] = -a0 * (v[id].Data[j] - currentSolution[j]) - a2 * v1[id].Data[j] - a3 * v2[id].Data[j];
            //    uc.Data[j] = -a1 * (v[id].Data[j] - currentSolution[j]) - a4 * v1[id].Data[j] - a5 * v2[id].Data[j];
            //}
            //provider.MassMatrixVectorProduct(subdomain, uu, tempResult.Data);
            //result.Add(tempResult);
            //provider.DampingMatrixVectorProduct(subdomain, uc, tempResult.Data);
            //result.Add(tempResult);

            ////CalculateRhsImplicit(subdomain, result.Data, false);
            ////result.Scale(-1d);
            #endregion

            // result = M * u
            return providers[0].MassMatrixVectorProduct(linearSystem.Subdomain, currentSolution);
        }

        public void Initialize(bool isFirstAnalysis = true)
        {
            //if (isFirstAnalysis)
            //{
            //    // The order in which the next initializations happen is very important.
            //    model.ConnectDataStructures();
            //    solver.OrderDofs(false);
            //    foreach (ILinearSystem linearSystem in linearSystems.Values)
            //    {
            //        linearSystem.Reset(); // Necessary to define the linear system's size 
            //        linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
            //    }
            //}
            //else
            //{
            //    foreach (ILinearSystem linearSystem in linearSystems.Values)
            //    {
            //        //TODO: Perhaps these shouldn't be done if an analysis has already been executed. The model will not be 
            //        //      modified. Why should the linear system be?
            //        linearSystem.Reset();
            //    }
            //}

            ////TODO: Perhaps this should be called by the child analyzer
            //BuildMatrices();

            //// Loads must be created after building the matrices.
            ////TODO: Some loads may not have to be recalculated each time the stiffness changes.
            //model.AssignLoads(solver.DistributeNodalLoads);
            //foreach (ILinearSystem linearSystem in linearSystems.Values)
            //{
            //    linearSystem.RhsVector = linearSystem.Subdomain.Forces;
            //}

            ////InitializeCoefficients();
            //InitializeInternalVectors();
            ////InitializeMatrices();
            //InitializeRhs();

            //if (ChildAnalyzer == null) throw new InvalidOperationException("Newmark analyzer must contain an embedded analyzer.");
            //ChildAnalyzer.Initialize(isFirstAnalysis);
        }

        private void InitializeInternal()
        {
            if (ChildAnalyzer == null) throw new InvalidOperationException("Newmark analyzer must contain an embedded analyzer.");
            for (int i = 0; i < linearSystems.Length; i++)
            {
                models[i].ConnectDataStructures();
                solvers[i].OrderDofs(false);
                foreach (ILinearSystem linearSystem in linearSystems[i].Values)
                {
                    linearSystem.Reset(); // Necessary to define the linear system's size 
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }
            }

            //TODO: Perhaps this should be called by the child analyzer
            BuildMatrices();

            for (int i = 0; i < linearSystems.Length; i++)
            {
                models[i].AssignLoads(solvers[i].DistributeNodalLoads);
                foreach (ILinearSystem linearSystem in linearSystems[i].Values)
                {
                    linearSystem.RhsVector = linearSystem.Subdomain.Forces;
                }

            }

            InitializeInternalVectors();
            for (int i = 0; i < linearSystems.Length; i++)
            {
                InitializeRhs(i);
                childAnalyzers[i].Initialize(true);
            }
        }

        public void SolveTimestep(int t)
        {
            DateTime start = DateTime.Now;
            Debug.WriteLine("Newmark step: {0}", t);

            int staggeredStep = 0;
            var temperatureNorm = 0d;
            var previousTemperatureNorm = 0d;
            var error = 1d;
            do
            {
                previousTemperatureNorm = temperatureNorm;

                for (int i = 0; i < linearSystems.Length; i++)
                {
                    unknownFromPreviousStaggeredStep[i].Clear();
                    foreach (var linearSystem in linearSystems[i].Values)
                    {
                        if (linearSystem.Solution != null)
                        {
                            unknownFromPreviousStaggeredStep[i].Add(linearSystem.Subdomain.ID, linearSystem.Solution.Copy());
                        }
                        this.linearSystems[i] = solvers[i].LinearSystems;
                    }
                }

                InitializeInternal();
                for (int i = 0; i < linearSystems.Length; i++)
                {
                    IDictionary<int, IVector> rhsVectors = providers[i].GetRhsFromHistoryLoad(t);
                    foreach (var l in linearSystems[i].Values) l.RhsVector = rhsVectors[l.Subdomain.ID];
                    InitializeRhs(i);
                    CalculateRhsImplicit(i);

                    childAnalyzers[i].Solve();
                }

                temperatureNorm = 0;
                for (int i = 0; i < linearSystems.Length; i++)
                {
                    foreach (var linearSystem in linearSystems[i].Values)
                    {
                        temperatureNorm += linearSystem.Solution.Norm2();
                    }
                }
                error = temperatureNorm != 0 ? Math.Abs(temperatureNorm - previousTemperatureNorm) / temperatureNorm : 0;
                Debug.WriteLine("Staggered step: {0} - error {1}", staggeredStep, error);
                staggeredStep++;

                CreateNewModel(unknown, models, solvers, providers, childAnalyzers);
            }
            while (staggeredStep < maxStaggeredSteps && error > tolerance);

            DateTime end = DateTime.Now;
            UpdateUnknown(t);
            UpdateResultStorages(start, end);
            Debug.WriteLine("-------------");
        }

        public void Solve()
        {
            int numTimeSteps = (int)(totalTime / timeStep);
            for (int t = 0; t < numTimeSteps; ++t)
            {
                SolveTimestep(t);
            }
        }

        private void CalculateRhsImplicit(int modelNo)
        {
            foreach (ILinearSystem linearSystem in linearSystems[modelNo].Values)
            {
                linearSystem.RhsVector = CalculateRhsImplicit(linearSystem, modelNo, true);
            }
        }

        private IVector CalculateRhsImplicit(ILinearSystem linearSystem, int modelNo, bool addRhs)
        {
            //TODO: what is the meaning of addRhs? Do we need this when solving dynamic thermal equations?
            //TODO: stabilizingRhs has not been implemented

            // result = -dt(conductuvity*temperature + rhs -dt(stabilizingConductivity*temperature + StabilizingRhs)) 
            int id = linearSystem.Subdomain.ID;
            firstOrderMatrixTimesUnknown[modelNo][id] = providers[modelNo].MassMatrixVectorProduct(linearSystem.Subdomain, unknown[modelNo][id]);
            secondOrderMatrixTimesUnknown[modelNo][id] = providers[modelNo].DampingMatrixVectorProduct(linearSystem.Subdomain, unknown[modelNo][id]);

            IVector rhsResult = rhsPrevious[modelNo][id].LinearCombination(1 - beta, rhs[modelNo][id], beta);
            rhsResult.AxpyIntoThis(firstOrderMatrixTimesUnknown[modelNo][id], 1 / timeStep);
            rhsResult.AxpyIntoThis(secondOrderMatrixTimesUnknown[modelNo][id], -(1 - beta));

            rhsPrevious[modelNo][id] = rhs[modelNo][id];
            return rhsResult;
        }

        private void InitializeInternalVectors()
        {
            for (int i = 0; i < linearSystems.Length; i++)
            {
                //temperature[i].Clear();
                rhs[i].Clear();
                rhsPrevious[i].Clear();
                secondOrderMatrixTimesUnknown[i].Clear();
                firstOrderMatrixTimesUnknown[i].Clear();

                foreach (ILinearSystem linearSystem in linearSystems[i].Values)
                {
                    int id = linearSystem.Subdomain.ID;
                    firstOrderMatrixTimesUnknown[i].Add(id, linearSystem.CreateZeroVector());
                    secondOrderMatrixTimesUnknown[i].Add(id, linearSystem.CreateZeroVector());
                    //temperature.Add(id, linearSystem.CreateZeroVector());
                    rhs[i].Add(id, linearSystem.CreateZeroVector());
                    rhsPrevious[i].Add(id, linearSystem.CreateZeroVector());

                    if (unknown[i].ContainsKey(id) == false)
                        unknown[i].Add(id, initialValue[i].Copy());

                    //// Account for initial conditions coming from a previous solution. 
                    ////TODO: This doesn't work as intended. The solver (previously the LinearSystem) initializes the solution to zero.
                    //if (linearSystem.Solution != null) temperature[i].Add(id, linearSystem.Solution.Copy());
                    //else temperature[i].Add(id, initialTemperature[i].Copy());
                }
            }
        }

        private void InitializeRhs(int modelNo)
        {
            ImplicitIntegrationCoefficients coeffs = new ImplicitIntegrationCoefficients
            {
                Mass = 0,  //never used
                Stiffness = 0
            };
            foreach (ILinearSystem linearSystem in linearSystems[modelNo].Values)
            {
                providers[modelNo].ProcessRhs(coeffs, linearSystem.Subdomain, linearSystem.RhsVector);
                rhs[modelNo][linearSystem.Subdomain.ID] = linearSystem.RhsVector.Copy(); //TODO: copying the vectors is wasteful.
            }
        }

        private void UpdateResultStorages(DateTime start, DateTime end)
        {
            for (int i = 0; i < linearSystems.Length; i++)
            {
                foreach (ILinearSystem linearSystem in linearSystems[i].Values)
                {
                    int id = linearSystem.Subdomain.ID;
                    if (ResultStorages.ContainsKey(id))
                        if (ResultStorages[id] != null)
                            foreach (var l in ChildAnalyzer.Logs[id])
                                ResultStorages[id].StoreResults(start, end, l);
                }
            }
        }

        private void UpdateUnknown(int timeStep)
        {
            for (int i = 0; i < linearSystems.Length; i++)
            {
                foreach (ILinearSystem linearSystem in linearSystems[i].Values)
                {
                    int id = linearSystem.Subdomain.ID;
                    unknown[i][id].CopyFrom(linearSystem.Solution);
                    unknown[i][id].AddIntoThis(linearSystem.Solution);
                    if ((timeStep + 1) % 1 == 0)
                    {
                        string path0 = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput";
                        //string path1 = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput\temperature0.txt";
                        //string path = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput";
                        var path2 = Path.Combine(path0, $"temperature{i}-{timeStep}.txt");
                        var writer = new LinearAlgebra.Output.FullVectorWriter() { ArrayFormat = Array1DFormat.PlainVertical };
                        writer.WriteToFile(unknown[i][id], path2);
                        //writer.WriteToFile(temperature[id][0], path1);

                        //File.AppendAllLines(path1, new string[] { temperature[id][0].ToString() }, Encoding.UTF8);
                        //File.AppendAllLines(path2, new string[] { temperature[id][340].ToString() }, Encoding.UTF8);
                    }
                }
            }
        }
    }
}
