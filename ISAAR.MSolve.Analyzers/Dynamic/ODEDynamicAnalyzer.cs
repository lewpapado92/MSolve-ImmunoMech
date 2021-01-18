using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization.Interfaces;
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
    public class ODEDynamicAnalyzer : INonLinearParentAnalyzer //TODO: why is this non linear
    {
        private readonly double beta, timeStep, totalTime;
        private readonly IStructuralModel model;
        private readonly IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        private readonly ISolver solver;
        private readonly IImplicitIntegrationProvider provider;
        private Dictionary<int, IVector> rhs = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> rhsPrevious = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> unknown = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> firstOrderMatrixTimesUnknown = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> zerothOrderMatrixTimesUnknown = new Dictionary<int, IVector>();

        public ODEDynamicAnalyzer(IStructuralModel model, ISolver solver, IImplicitIntegrationProvider provider,
            IChildAnalyzer childAnalyzer, double beta, double timeStep, double totalTime, IVector initialValue = null)
        {
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.solver = solver;
            //solver.PreventFromOverwrittingSystemMatrices(); //TODO: If the scheme is purely implicit we can overwrite the matrix.
            this.provider = provider;
            this.ChildAnalyzer = childAnalyzer;
            this.beta = beta;
            this.timeStep = timeStep;
            this.totalTime = totalTime;
            this.ChildAnalyzer.ParentAnalyzer = this;
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
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.Matrix = provider.LinearCombinationOfMatricesIntoStiffness(coeffs, linearSystem.Subdomain);
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
            return provider.MassMatrixVectorProduct(linearSystem.Subdomain, currentSolution);
        }

        public void Initialize(bool isFirstAnalysis = true)
        {
            if (isFirstAnalysis)
            {
                // The order in which the next initializations happen is very important.
                model.ConnectDataStructures();
                solver.OrderDofs(false);
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    linearSystem.Reset(); // Necessary to define the linear system's size 
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }
            }
            else
            {
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    //TODO: Perhaps these shouldn't be done if an analysis has already been executed. The model will not be 
                    //      modified. Why should the linear system be?
                    linearSystem.Reset();
                }
            }

            //TODO: Perhaps this should be called by the child analyzer
            BuildMatrices();

            // Loads must be created after building the matrices.
            //TODO: Some loads may not have to be recalculated each time the stiffness changes.
            model.AssignLoads(solver.DistributeNodalLoads);
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.RhsVector = linearSystem.Subdomain.Forces;
            }

            //InitializeCoefficients();
            InitializeInternalVectors();
            //InitializeMatrices();
            InitializeRhs();

            if (ChildAnalyzer == null) throw new InvalidOperationException("Dynamic analyzer must contain an embedded analyzer.");
            ChildAnalyzer.Initialize(isFirstAnalysis);
        }

        public void Solve()
        {
            int numTimeSteps = (int)(totalTime / timeStep);
            for (int t = 0; t < numTimeSteps; ++t)
            {
                Debug.WriteLine("Dynamic Analyzer step: {0}", t);

                IDictionary<int, IVector> rhsVectors = provider.GetRhsFromHistoryLoad(t);
                foreach (var l in linearSystems.Values) l.RhsVector = rhsVectors[l.Subdomain.ID];
                InitializeRhs();
                CalculateRhsImplicit();

                DateTime start = DateTime.Now;
                ChildAnalyzer.Solve();
                DateTime end = DateTime.Now;

                UpdateUnknown(t);
                UpdateResultStorages(start, end);
            }
        }

        private void CalculateRhsImplicit()
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.RhsVector = CalculateRhsImplicit(linearSystem, true);
            }
        }

        private IVector CalculateRhsImplicit(ILinearSystem linearSystem, bool addRhs)
        {
            //TODO: what is the meaning of addRhs? Do we need this when solving dynamic thermal equations?

            // result = (1-b)* rhsPrevious + beta * rhs + 1/dt * 1stOrderMatrix * unknown - (1-b) * 0thOrderMatrix * unknown
            int id = linearSystem.Subdomain.ID;
            firstOrderMatrixTimesUnknown[id] = provider.MassMatrixVectorProduct(linearSystem.Subdomain, unknown[id]);
            zerothOrderMatrixTimesUnknown[id] = provider.DampingMatrixVectorProduct(linearSystem.Subdomain, unknown[id]);

            IVector rhsResult = rhsPrevious[id].LinearCombination(1 - beta, rhs[id], beta);
            rhsResult.AxpyIntoThis(firstOrderMatrixTimesUnknown[id], 1 / timeStep);
            rhsResult.AxpyIntoThis(zerothOrderMatrixTimesUnknown[id], -(1 - beta));

            rhsPrevious[id] = rhs[id];
            return rhsResult;
        }

        private void InitializeInternalVectors()
        {
            unknown.Clear();
            rhs.Clear();
            rhsPrevious.Clear();
            firstOrderMatrixTimesUnknown.Clear();
            zerothOrderMatrixTimesUnknown.Clear();

            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                firstOrderMatrixTimesUnknown.Add(id, linearSystem.CreateZeroVector());
                zerothOrderMatrixTimesUnknown.Add(id, linearSystem.CreateZeroVector());
                //temperature.Add(id, linearSystem.CreateZeroVector());
                rhs.Add(id, linearSystem.CreateZeroVector());
                rhsPrevious.Add(id, linearSystem.CreateZeroVector());

                // Account for initial conditions coming from a previous solution. 
                //TODO: This doesn't work as intended. The solver (previously the LinearSystem) initializes the solution to zero.
                if (linearSystem.Solution != null) unknown.Add(id, linearSystem.Solution.Copy());
                else unknown.Add(id, linearSystem.CreateZeroVector());
            }
        }

        private void InitializeRhs()
        {
            ImplicitIntegrationCoefficients coeffs = new ImplicitIntegrationCoefficients
            {
                Mass = 0,  //never used
                Stiffness = 0
            };
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                provider.ProcessRhs(coeffs, linearSystem.Subdomain, linearSystem.RhsVector);
                rhs[linearSystem.Subdomain.ID] = linearSystem.RhsVector.Copy(); //TODO: copying the vectors is wasteful.
            }
        }

        private void UpdateResultStorages(DateTime start, DateTime end)
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                if (ResultStorages.ContainsKey(id))
                    if (ResultStorages[id] != null)
                        foreach (var l in ChildAnalyzer.Logs[id])
                            ResultStorages[id].StoreResults(start, end, l);
            }
        }

        private void UpdateUnknown(int timeStep)
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                unknown[id].CopyFrom(linearSystem.Solution);
                unknown[id].AddIntoThis(linearSystem.Solution);
                if ((timeStep + 1) % 1 == 0)
                {
                    string path0 = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput";
                    //string path1 = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput\temperature0.txt";
                    //string path = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput";
                    var path2 = Path.Combine(path0, $"unknown-{timeStep}.txt");
                    var writer = new LinearAlgebra.Output.FullVectorWriter() { ArrayFormat = Array1DFormat.PlainVertical };
                    writer.WriteToFile(unknown[id], path2);
                    //writer.WriteToFile(temperature[id][0], path1);

                    //File.AppendAllLines(path1, new string[] { temperature[id][0].ToString() }, Encoding.UTF8);
                    //File.AppendAllLines(path2, new string[] { temperature[id][340].ToString() }, Encoding.UTF8);
                }
            }
        }
    }
}
