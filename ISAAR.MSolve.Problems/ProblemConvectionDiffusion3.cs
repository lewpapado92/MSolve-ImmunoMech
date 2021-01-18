using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.FEM.Loading.Providers;
using ISAAR.MSolve.FEM.Providers;
using ISSAR.MSolve.Discretization.Loads;

//TODO: Usually the LinearSystem is passed in, but for GetRHSFromHistoryLoad() it is stored as a field. Decide on one method.
//TODO: I am not too fond of the provider storing global sized matrices.
namespace ISAAR.MSolve.Problems
{
    public class ProblemConvectionDiffusion3 : IConvectionDiffusionIntegrationProvider, IStaticProvider, INonLinearProvider
    {
        private Dictionary<int, IMatrix> capacity, diffusionConductivityFreeFree, massTransportConductivityFreeFree, stabilizingConductivity;
        private Dictionary<int, IVector> stabilizingDomainLoad;
        private Dictionary<int, IMatrixView> diffusionConductivityFreeConstr, diffusionConductivityConstrFree, diffusionConductivityConstrConstr,
            massTransportConductivityFreeConstr, massTransportConductivityConstrFree, massTransportConductivityConstrConstr;
        private readonly Model model;
        private readonly ISolver solver;
        private IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        private ElementMassTransportConductivityProvider massTransportConductivityProvider = new ElementMassTransportConductivityProvider();
        private ElementDiffusionConductivityProvider diffusionConductivityProvider = new ElementDiffusionConductivityProvider();
        private ElementStructuralStiffnessProvider conductivityProvider = new ElementStructuralStiffnessProvider();
        private ElementStructuralMassProvider capacityProvider = new ElementStructuralMassProvider();
        private ElementStructuralDampingProvider stabilizingConductivityProvider = new ElementStructuralDampingProvider();

        public ProblemConvectionDiffusion3(Model model, ISolver solver)
        {
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.solver = solver;
            this.DirichletLoadsAssembler = new DirichletEquivalentLoadsStructural(conductivityProvider);
        }

        public IDirichletEquivalentLoadsAssembler DirichletLoadsAssembler { get; }

        private IDictionary<int, IMatrix> Capacity
        {
            get
            {
                if (capacity == null) BuildCapacity();
                return capacity;
            }
        }

        private IDictionary<int, IMatrix> DiffusionConductivity
        {
            get
            {
                if (diffusionConductivityFreeFree == null) BuildDiffusionConductivityFreeFree();
                //else RebuildConductivityMatrices();
                return diffusionConductivityFreeFree;
            }
        }
        private IDictionary<int, IMatrix> MassTransportConductivity
        {
            get
            {
                if (massTransportConductivityFreeFree == null) BuildMassTransportConductivityFreeFree();
                //else RebuildConductivityMatrices();
                return massTransportConductivityFreeFree;
            }
        }

        private IDictionary<int, IMatrix> StabilizingConductivity
        {
            get
            {
                if (stabilizingConductivity == null) BuildStabilizingConductivity();
                return stabilizingConductivity;
            }
        }
        private IDictionary<int, IVector> StabilizingDomainLoad
        {
            get
            {
                if (stabilizingDomainLoad == null) stabilizingDomainLoad = model.BuildGlobalStabilizingBodyLoads(solver.DistributeNodalLoads);
                return stabilizingDomainLoad;
            }
        }

        private void BuildDiffusionConductivityFreeFree()
        {
            double[,] kappa = new double[model.Nodes.Count, model.Nodes.Count];
            kappa[0, 0] = 1;
            kappa[model.Nodes.Count - 1, model.Nodes.Count - 1] = 1;
            IMatrix penalty = Matrix.CreateFromArray(kappa);
            diffusionConductivityFreeFree = solver.BuildGlobalMatrices(diffusionConductivityProvider);
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                diffusionConductivityFreeFree[id].AddIntoThis(penalty);
            }
        }
        private void BuildDiffusionConductivitySubmatrices()
        {
            Dictionary<int, (IMatrix Cff, IMatrixView Cfc, IMatrixView Ccf, IMatrixView Ccc)> matrices =
                solver.BuildGlobalSubmatrices(diffusionConductivityProvider);

            diffusionConductivityFreeFree = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            diffusionConductivityFreeConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            diffusionConductivityConstrFree = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            diffusionConductivityConstrConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                diffusionConductivityFreeFree.Add(id, matrices[id].Cff);
                diffusionConductivityFreeConstr.Add(id, matrices[id].Cfc);
                diffusionConductivityConstrFree.Add(id, matrices[id].Ccf);
                diffusionConductivityConstrConstr.Add(id, matrices[id].Ccc);
            }
        }
        private void BuildMassTransportConductivityFreeFree() => massTransportConductivityFreeFree = solver.BuildGlobalMatrices(massTransportConductivityProvider);

        private void BuildMassTransportConductivitySubmatrices()
        {
            Dictionary<int, (IMatrix Cff, IMatrixView Cfc, IMatrixView Ccf, IMatrixView Ccc)> matrices =
                solver.BuildGlobalSubmatrices(massTransportConductivityProvider);

            massTransportConductivityFreeFree = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            massTransportConductivityFreeConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            massTransportConductivityConstrFree = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            massTransportConductivityConstrConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                massTransportConductivityFreeFree.Add(id, matrices[id].Cff);
                massTransportConductivityFreeConstr.Add(id, matrices[id].Cfc);
                massTransportConductivityConstrFree.Add(id, matrices[id].Ccf);
                massTransportConductivityConstrConstr.Add(id, matrices[id].Ccc);
            }
        }

        private void RebuildDiffusionConductivityFreeFree()
        {
            //TODO: This will rebuild all the stiffnesses of all subdomains, if even one subdomain has MaterialsModified = true.
            //      Optimize this, by passing a flag foreach subdomain to solver.BuildGlobalSubmatrices().

            bool mustRebuild = false;
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                if (subdomain.StiffnessModified)
                {
                    mustRebuild = true;
                    break;
                }
            }
            if (mustRebuild) diffusionConductivityFreeFree = solver.BuildGlobalMatrices(diffusionConductivityProvider);
            foreach (ISubdomain subdomain in model.Subdomains) subdomain.ResetMaterialsModifiedProperty();
        }
        private void RebuildMassTransportConductivityFreeFree()
        {
            //TODO: This will rebuild all the stiffnesses of all subdomains, if even one subdomain has MaterialsModified = true.
            //      Optimize this, by passing a flag foreach subdomain to solver.BuildGlobalSubmatrices().

            bool mustRebuild = false;
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                if (subdomain.StiffnessModified)
                {
                    mustRebuild = true;
                    break;
                }
            }
            if (mustRebuild) massTransportConductivityFreeFree = solver.BuildGlobalMatrices(massTransportConductivityProvider);
            foreach (ISubdomain subdomain in model.Subdomains) subdomain.ResetMaterialsModifiedProperty();
        }

        private void BuildCapacity() => capacity = solver.BuildGlobalMatrices(capacityProvider);
        private void BuildStabilizingConductivity() => stabilizingConductivity = solver.BuildGlobalMatrices(stabilizingConductivityProvider);
        //private void BuildStabilizingDomainLoad() => stabilizingDomainLoad = solver.DistributeNodalLoads(stabilizingDomainLoadProvider.StabilizingLoad);

        #region IAnalyzerProvider Members
        public void ClearMatrices()
        {
            capacity = null;
            diffusionConductivityFreeFree = null;
            diffusionConductivityFreeConstr = null;
            diffusionConductivityConstrFree = null;
            diffusionConductivityConstrConstr = null;
            massTransportConductivityFreeFree = null;
            massTransportConductivityFreeConstr = null;
            massTransportConductivityConstrFree = null;
            massTransportConductivityConstrConstr = null;
        }

        public void Reset()
        {
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                foreach (IElement element in subdomain.Elements)
                {
                    ((IFiniteElement)element.ElementType).ClearMaterialState();
                }
            }

            diffusionConductivityFreeFree = null;
            diffusionConductivityConstrFree = null;
            diffusionConductivityConstrConstr = null;
            massTransportConductivityFreeFree = null;
            massTransportConductivityConstrFree = null;
            massTransportConductivityConstrConstr = null;
            capacity = null;
        }
        #endregion 

        #region IImplicitIntegrationProvider Members

        public IMatrixView LinearCombinationOfMatricesIntoStiffness(ImplicitIntegrationCoefficients coefficients,
            ISubdomain subdomain)
        {
            int id = subdomain.ID;
            var conductivity = DiffusionConductivity[id].Add(MassTransportConductivity[id]);
            conductivity.ScaleIntoThis(coefficients.Stiffness);
            var capacityTimesCoeff = Capacity[id].Scale(coefficients.Mass);
            capacityTimesCoeff.AddIntoThis(conductivity);
            return capacityTimesCoeff.LinearCombination(1,StabilizingConductivity[id],coefficients.Damping);
        }

        public void ProcessRhs(ISubdomain subdomain, IVector rhs)
        {
            // Method intentionally left empty.
        }

        public IDictionary<int, IVector> GetAccelerationsOfTimeStep(int timeStep)
        {
            throw new InvalidOperationException("This is does not make sense in explicit methods for first order equations");
        }

        public IDictionary<int, IVector> GetVelocitiesOfTimeStep(int timeStep)
        {
            throw new InvalidOperationException("This is not needed in explicit methods for first order equations");
        }

        public IDictionary<int, IVector> GetRhsFromHistoryLoad(int timeStep)
        {
            foreach (ISubdomain subdomain in model.Subdomains) subdomain.Forces.Clear(); //TODO: this is also done by model.AssignLoads()

            model.AssignLoads(solver.DistributeNodalLoads); // Time-independent nodal loads
            model.AssignTimeDependentNodalLoads(timeStep, solver.DistributeNodalLoads); // Time-dependent nodal loads

            var rhsVectors = new Dictionary<int, IVector>();
            foreach (ISubdomain subdomain in model.Subdomains) rhsVectors.Add(subdomain.ID, subdomain.Forces.Copy());
            return rhsVectors;
        }
        public IVector ConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
        {
            throw new NotImplementedException();
        }
        public IVector CapacityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.Capacity[subdomain.ID].Multiply(vector);
        public IVector DiffusionConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.DiffusionConductivity[subdomain.ID].Multiply(vector);
        public IVector MassTransportConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.MassTransportConductivity[subdomain.ID].Multiply(vector);

        //TODO: Ok this is weird. These methods should be named Second/First/ZeroOrderCoefficientTimesVector()
        public IVector StabilizingConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.StabilizingConductivity[subdomain.ID].Multiply(vector);
        public IVector StabilizingRhs(ISubdomain subdomain)
            => this.StabilizingDomainLoad[subdomain.ID];

        #endregion

        #region IStaticProvider Members

        public IMatrixView CalculateMatrix(ISubdomain subdomain)
        {
            if (diffusionConductivityFreeFree == null) BuildDiffusionConductivityFreeFree();
            if (massTransportConductivityFreeFree == null) BuildMassTransportConductivityFreeFree();
            return diffusionConductivityFreeFree[subdomain.ID].LinearCombination(1,massTransportConductivityFreeFree[subdomain.ID],1);
        }

        public (IMatrixView matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr) CalculateSubMatrices(ISubdomain subdomain)
        {
            int id = subdomain.ID;
            if ((diffusionConductivityFreeFree == null) || (diffusionConductivityFreeConstr == null)
                || (diffusionConductivityConstrFree == null) || (diffusionConductivityConstrConstr == null))
            {
                BuildDiffusionConductivitySubmatrices();
            }
            if ((massTransportConductivityFreeFree == null) || (massTransportConductivityFreeConstr == null)
                || (massTransportConductivityConstrFree == null) || (massTransportConductivityConstrConstr == null))
            {
                BuildMassTransportConductivitySubmatrices();
            }
            var cff = diffusionConductivityFreeFree[id].LinearCombination(1, massTransportConductivityFreeFree[id], 1);
            var cfc = diffusionConductivityFreeConstr[id].LinearCombination(1, massTransportConductivityFreeConstr[id], 1);
            var ccf = diffusionConductivityConstrFree[id].LinearCombination(1, massTransportConductivityConstrFree[id], 1);
            var ccc = diffusionConductivityConstrConstr[id].LinearCombination(1, massTransportConductivityConstrConstr[id], 1);
            return (cff, cfc, ccf, ccc);
        }
        #endregion

        #region INonLinearProvider Members

        public double CalculateRhsNorm(IVectorView rhs) => rhs.Norm2();

        public void ProcessInternalRhs(ISubdomain subdomain, IVectorView solution, IVector rhs) { }

        #endregion
    }
}
