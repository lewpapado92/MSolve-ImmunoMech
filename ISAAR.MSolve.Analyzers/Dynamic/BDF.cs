using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
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
    /// <summary>
    /// 
    /// Authors: Odysseas Kokkinos
    /// </summary>
    public class BDF : INonLinearParentAnalyzer //TODO: why is this non linear
    {
        private readonly double timeStep, totalTime;
        private readonly IStructuralModel model;
        private readonly IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        private readonly ISolver solver;
        private readonly IConvectionDiffusionIntegrationProvider provider;
        private readonly IVector initialTemperature;
        private readonly int BDForder;
        private Dictionary<int, IVector> rhs = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> stabilizingRhs = new Dictionary<int, IVector>();//TODO: has to be implemented, pertains to domain loads
        private Dictionary<int, IVector> rhsPrevious = new Dictionary<int, IVector>();//TODO: at the moment domain loads are not implemented in this
        public IList<Dictionary<int, IVector>> temperature = new List<Dictionary<int, IVector>>();
        //private Dictionary<int, IVector> conductivityTimesTemperature = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> capacityTimesTemperature = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> diffusionConductivityTimesTemperature = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> massTransportConductivityTimesTemperature = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> stabilizingConductivityTimesTemperature = new Dictionary<int, IVector>();
        private static double[,] BDFcoeff = new double[,]{{-1,1,1,0,0,0,0,0},
            { 1, -4, 3, 2, 0, 0, 0, 0},
            { -2, 9, -18, 11, 6, 0, 0, 0},
            { 3, -16, 36, -48, 25, 12, 0, 0},
            { -12, 75, -200, 300, -300, 137, 60, 0},
            { 10, -72, 225, -400, 450, -360, 147, 60},};
        //private Matrix BDFmatrix = Matrix.CreateFromArray(BDFcoeff);
        private Dictionary<int, IVector> effectiveTemperature = new Dictionary<int, IVector>();

        public BDF(IStructuralModel model, ISolver solver, IConvectionDiffusionIntegrationProvider provider,
            IChildAnalyzer childAnalyzer, double timeStep, double totalTime,int BDForder, IVector initialTemperature = null)
        {
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.solver = solver;
            //solver.PreventFromOverwrittingSystemMatrices(); //TODO: If the scheme is purely implicit we can overwrite the matrix.
            this.provider = provider;
            this.ChildAnalyzer = childAnalyzer;
            this.timeStep = timeStep;
            this.totalTime = totalTime;
            this.ChildAnalyzer.ParentAnalyzer = this;
            this.BDForder = BDForder;
            this.initialTemperature = initialTemperature;

        }

        public Dictionary<int, IAnalyzerLog[]> Logs => null; //TODO: this can't be right
        public Dictionary<int, ImplicitIntegrationAnalyzerLog> ResultStorages { get; }
            = new Dictionary<int, ImplicitIntegrationAnalyzerLog>();

        public IChildAnalyzer ChildAnalyzer { get; }

        public void BuildMatrices()
        {
            if (timeIncrement < BDForder-1)
            {
                var coeffs = new ImplicitIntegrationCoefficients
                {
                    Mass = 1 / timeStep / timeStep,
                    Damping = -1,
                    Stiffness = 1 / timeStep
                };
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    linearSystem.Matrix = provider.LinearCombinationOfMatricesIntoStiffness(coeffs, linearSystem.Subdomain);
                }
            }
            else
            {
                var coeffs = new ImplicitIntegrationCoefficients
                {
                    Mass = BDFcoeff[BDForder - 1, BDForder] / timeStep / timeStep,
                    Damping = -BDFcoeff[BDForder - 1, BDForder + 1],
                    Stiffness = BDFcoeff[BDForder - 1, BDForder + 1] / timeStep
                };
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    linearSystem.Matrix = provider.LinearCombinationOfMatricesIntoStiffness(coeffs, linearSystem.Subdomain);
                }
            }
        }

        private int timeIncrement = 0;
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

            // result = Κ * u
            return provider.ConductivityMatrixVectorProduct(linearSystem.Subdomain, currentSolution);
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

            if (ChildAnalyzer == null) throw new InvalidOperationException("Newmark analyzer must contain an embedded analyzer.");
            ChildAnalyzer.Initialize(isFirstAnalysis);
        }

        public void Solve()
        {
            int numTimeSteps = (int)(totalTime / timeStep);
            for (timeIncrement = 0; timeIncrement < numTimeSteps; ++timeIncrement)
            {
                Debug.WriteLine("Newmark step: {0}", timeIncrement);

                IDictionary<int, IVector> rhsVectors = provider.GetRhsFromHistoryLoad(timeIncrement);
                foreach (var l in linearSystems.Values) l.RhsVector = rhsVectors[l.Subdomain.ID];
                InitializeRhs();
                CalculateRhsImplicit(timeIncrement);
                if (timeIncrement == BDForder-1) BuildMatrices();

                DateTime start = DateTime.Now;
                ChildAnalyzer.Solve();
                DateTime end = DateTime.Now;

                UpdateTemperature(timeIncrement);
                UpdateResultStorages(start, end);
            }
        }

        private void CalculateRhsImplicit(int currentTimeStep)
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.RhsVector = CalculateRhsImplicit(linearSystem, currentTimeStep);
            }
        }
        private IVector CalculateRhsImplicit(ILinearSystem linearSystem, int ti)
        {
            //TODO: what is the meaning of addRhs? Do we need this when solving dynamic thermal equations?
            //TODO: stabilizingRhs has not been implemented

            // result = (capacity/dt^2-diffusionConductivity/dt-massTransportConductivity/dt+stabilizingConductivity)*temperature - (StabilizingRhs - rhs/dt)) 
            // result = -dt(conductuvity*temperature - rhs -dt(stabilizingConductivity*temperature + StabilizingRhs)) 
            int id = linearSystem.Subdomain.ID;
            double a0 = 1 / Math.Pow(timeStep, 2);
            double a1 = 1 / (2 * timeStep);
            double a2 = 1 / timeStep;
            if (timeIncrement < BDForder-1)
            {
                capacityTimesTemperature[id] = provider.CapacityMatrixVectorProduct(linearSystem.Subdomain, temperature[BDForder][id]);
                capacityTimesTemperature[id].ScaleIntoThis(a0);
                //stabilizingConductivityTimesTemperature[id] = provider.StabilizingConductivityMatrixVectorProduct(linearSystem.Subdomain, temperature[id]);
                rhs[id].ScaleIntoThis(a2);
                var rhsResult = capacityTimesTemperature[id].Subtract(stabilizingRhs[id]);
                rhsResult.AddIntoThis(rhs[id]);
                return rhsResult;
            }
            else
            {
                effectiveTemperature[id].Clear();
                for (int i = 0; i < BDForder; i++)
                {
                    var a = temperature[i][id].Scale(BDFcoeff[BDForder - 1, i]);
                    effectiveTemperature[id].AddIntoThis(a); 
                }
                capacityTimesTemperature[id] = provider.CapacityMatrixVectorProduct(linearSystem.Subdomain, effectiveTemperature[id]);
                capacityTimesTemperature[id].ScaleIntoThis(-a0);
                //stabilizingConductivityTimesTemperature[id] = provider.StabilizingConductivityMatrixVectorProduct(linearSystem.Subdomain, temperature[id]);
                var rhsLoads = rhs[id].LinearCombination(a2, stabilizingRhs[id], 1);
                rhsLoads.ScaleIntoThis(BDFcoeff[BDForder - 1, BDForder + 1]);
                var rhsResult = capacityTimesTemperature[id].Add(rhsLoads);
                return rhsResult;
            }
        }


        private void InitializeInternalVectors()
        {
            for (int i = 0; i <= BDForder; i++)
            {
                temperature.Add(new Dictionary<int, IVector>());
            }
            //temperature.Clear();
            rhs.Clear();
            stabilizingRhs.Clear();
            rhsPrevious.Clear();
            diffusionConductivityTimesTemperature.Clear();
            massTransportConductivityTimesTemperature.Clear();
            stabilizingConductivityTimesTemperature.Clear();
            capacityTimesTemperature.Clear();
            effectiveTemperature.Clear();

            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                diffusionConductivityTimesTemperature.Add(id, linearSystem.CreateZeroVector());
                massTransportConductivityTimesTemperature.Add(id, linearSystem.CreateZeroVector());
                stabilizingConductivityTimesTemperature.Add(id, linearSystem.CreateZeroVector());
                capacityTimesTemperature.Add(id, linearSystem.CreateZeroVector());
                //temperature.Add(id, linearSystem.CreateZeroVector());
                rhs.Add(id, linearSystem.CreateZeroVector());
                stabilizingRhs.Add(id, linearSystem.CreateZeroVector());
                rhsPrevious.Add(id, linearSystem.CreateZeroVector());
                effectiveTemperature.Add(id, linearSystem.CreateZeroVector());

                // Account for initial conditions coming from a previous solution. 
                //TODO: This doesn't work as intended. The solver (previously the LinearSystem) initializes the solution to zero.
                for (int i = 0; i <= BDForder; i++)
                {
                    if (linearSystem.Solution != null) temperature[i].Add(id, linearSystem.Solution.Copy());
                    else temperature[i].Add(id, initialTemperature.Copy());
                }
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
                provider.ProcessRhs(linearSystem.Subdomain, linearSystem.RhsVector);
                rhs[linearSystem.Subdomain.ID] = linearSystem.RhsVector.Copy(); //TODO: copying the vectors is wasteful.
                stabilizingRhs[linearSystem.Subdomain.ID] = provider.StabilizingRhs(linearSystem.Subdomain);
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

        private void UpdateTemperature(int timeStep)
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                temperature[BDForder][id].CopyFrom(linearSystem.Solution);
                for (int i = 0; i < BDForder; i++)
                {
                    temperature[i][id].CopyFrom(temperature[i+1][id]);
                }
                if ((timeStep + 1) % 1 == 0)
                {
                    //string path0 = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput";
                    ////string path1 = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput\temperature0.txt";
                    ////string path = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput";
                    //var path2 = Path.Combine(path0, $"temperature{timeStep}.txt");
                    //var writer = new LinearAlgebra.Output.FullVectorWriter() { ArrayFormat = Array1DFormat.PlainVertical };
                    //writer.WriteToFile(temperature[BDForder][id], path2);
                    ////writer.WriteToFile(temperature[id][0], path1);

                    ////File.AppendAllLines(path1, new string[] { temperature[id][0].ToString() }, Encoding.UTF8);
                    ////File.AppendAllLines(path2, new string[] { temperature[id][340].ToString() }, Encoding.UTF8);
                    string path0 = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput";
                    //string path1 = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput\temperature0.txt";
                    //string path = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput";
                    var path2 = Path.Combine(path0, $"temperature{timeStep}.txt");
                    var writer = new LinearAlgebra.Output.FullVectorWriter() { ArrayFormat = Array1DFormat.PlainHorizontal };
                    writer.WriteToFile(temperature[BDForder][id], path2);
                    //writer.WriteToFile(temperature[id][0], path1);

                    //File.AppendAllLines(path1, new string[] { temperature[id][0].ToString() }, Encoding.UTF8);
                    //File.AppendAllLines(path2, new string[] { temperature[id][340].ToString() }, Encoding.UTF8);
                }
            }
        }
    }
}
