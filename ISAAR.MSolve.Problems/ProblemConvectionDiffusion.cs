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
using ISSAR.MSolve.Discretization.Loads;

//TODO: Usually the LinearSystem is passed in, but for GetRHSFromHistoryLoad() it is stored as a field. Decide on one method.
//TODO: I am not too fond of the provider storing global sized matrices.
namespace ISAAR.MSolve.Problems
{
    public class ProblemConvectionDiffusion : IConvectionDiffusionIntegrationProvider, IStaticProvider, INonLinearProvider
    {
        private Dictionary<int, IMatrix> capacity, conductivityFreeFree, stabilizingConductivity;
        private Dictionary<int, IVector> stabilizingDomainLoad;
        private Dictionary<int, IMatrixView> conductivityFreeConstr, conductivityConstrFree, conductivityConstrConstr;
        private readonly Model model;
        private readonly ISolver solver;
        private IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        private ElementStructuralStiffnessProvider conductivityProvider = new ElementStructuralStiffnessProvider();
        private ElementStructuralMassProvider capacityProvider = new ElementStructuralMassProvider();
        private ElementStructuralDampingProvider stabilizingConductivityProvider = new ElementStructuralDampingProvider();

        public ProblemConvectionDiffusion(Model model, ISolver solver)
        {
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.solver = solver;
            //this.DirichletLoadsAssembler = new DirichletEquivalentLoadsStructural(conductivityProvider);
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

        private IDictionary<int, IMatrix> Conductivity
        {
            get
            {
                if (conductivityFreeFree == null) BuildConductivityFreeFree();
                //else RebuildConductivityMatrices();
                return conductivityFreeFree;
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

        private void BuildConductivityFreeFree() => conductivityFreeFree = solver.BuildGlobalMatrices(conductivityProvider);

        private void BuildConductivitySubmatrices()
        {
            Dictionary<int, (IMatrix Cff, IMatrixView Cfc, IMatrixView Ccf, IMatrixView Ccc)> matrices =
                solver.BuildGlobalSubmatrices(conductivityProvider);

            conductivityFreeFree = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            conductivityFreeConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            conductivityConstrFree = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            conductivityConstrConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                conductivityFreeFree.Add(id, matrices[id].Cff);
                conductivityFreeConstr.Add(id, matrices[id].Cfc);
                conductivityConstrFree.Add(id, matrices[id].Ccf);
                conductivityConstrConstr.Add(id, matrices[id].Ccc);
            }
        }

        private void RebuildConductivityFreeFree()
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
            if (mustRebuild) conductivityFreeFree = solver.BuildGlobalMatrices(conductivityProvider);
            foreach (ISubdomain subdomain in model.Subdomains) subdomain.ResetMaterialsModifiedProperty();
        }

        private void BuildCapacity() => capacity = solver.BuildGlobalMatrices(capacityProvider);
        private void BuildStabilizingConductivity() => stabilizingConductivity = solver.BuildGlobalMatrices(stabilizingConductivityProvider);
        //private void BuildStabilizingDomainLoad() => stabilizingDomainLoad = solver.DistributeNodalLoads(stabilizingDomainLoadProvider.StabilizingLoad);

        #region IAnalyzerProvider Members
        public void ClearMatrices()
        {
            capacity = null;
            conductivityFreeFree = null;
            conductivityFreeConstr = null;
            conductivityConstrFree = null;
            conductivityConstrConstr = null;
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

            conductivityFreeFree = null;
            conductivityConstrFree = null;
            conductivityConstrConstr = null;
            capacity = null;
        }
        #endregion 

        #region IImplicitIntegrationProvider Members

        public IMatrixView LinearCombinationOfMatricesIntoStiffness(ImplicitIntegrationCoefficients coefficients, ISubdomain subdomain)
        {
            int id = subdomain.ID;
            return Capacity[id];
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
        public IVector CapacityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.Capacity[subdomain.ID].Multiply(vector);

        public IVector ConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.Conductivity[subdomain.ID].Multiply(vector);

        //TODO: Ok this is weird. These methods should be named Second/First/ZeroOrderCoefficientTimesVector()
        public IVector StabilizingConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.StabilizingConductivity[subdomain.ID].Multiply(vector);
        public IVector StabilizingRhs(ISubdomain subdomain)
            => this.StabilizingDomainLoad[subdomain.ID];

        #endregion

        #region IStaticProvider Members

        public IMatrixView CalculateMatrix(ISubdomain subdomain)
        {
            if (conductivityFreeFree == null) BuildConductivityFreeFree();
            return conductivityFreeFree[subdomain.ID];
        }

        public (IMatrixView matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr) CalculateSubMatrices(ISubdomain subdomain)
        {
            int id = subdomain.ID;
            if ((conductivityFreeFree == null) || (conductivityFreeConstr == null)
                || (conductivityConstrFree == null) || (conductivityConstrConstr == null))
            {
                BuildConductivitySubmatrices();
            }
            return (conductivityFreeFree[id], conductivityFreeConstr[id],
                conductivityConstrFree[id], conductivityConstrConstr[id]);
        }
        #endregion

        #region INonLinearProvider Members

        public double CalculateRhsNorm(IVectorView rhs) => rhs.Norm2();

        public void ProcessInternalRhs(ISubdomain subdomain, IVectorView solution, IVector rhs) { }
        public IVector MassTransportConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
        {
            throw new NotImplementedException();
        }

        public IVector DiffusionConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
        {
            throw new NotImplementedException();
        }


        #endregion
    }
}
