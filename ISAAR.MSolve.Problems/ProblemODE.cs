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
using ISSAR.MSolve.Discretization.Loads;

//TODO: Usually the LinearSystem is passed in, but for GetRHSFromHistoryLoad() it is stored as a field. Decide on one method.
//TODO: I am not too fond of the provider storing global sized matrices.
namespace ISAAR.MSolve.Problems
{
    public class ProblemODE : IImplicitIntegrationProvider, IStaticProvider, INonLinearProvider
    {
        private Dictionary<int, IMatrix> firstOrderMatrix, zerothOrderMatrixFreeFree;
        private Dictionary<int, IMatrixView> zerothOrderMatrixFreeConstr, zerothOrderMatrixConstrFree, zerothOrderMatrixConstrConstr;
        private readonly Model model;
        private readonly ISolver solver;
        private IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        private ElementStructuralStiffnessProvider zerothOrderMatrixProvider = new ElementStructuralStiffnessProvider();
        private ElementStructuralMassProvider firstOrderMatrixProvider = new ElementStructuralMassProvider();

        public ProblemODE(Model model, ISolver solver)
        {
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.solver = solver;
            this.DirichletLoadsAssembler = new DirichletEquivalentLoadsStructural(zerothOrderMatrixProvider);
        }

        public IDirichletEquivalentLoadsAssembler DirichletLoadsAssembler { get; }

        private IDictionary<int, IMatrix> FirstOrderMatrix
        {
            get
            {
                if (firstOrderMatrix == null) BuildFirstOrderMatrix();
                return firstOrderMatrix;
            }
        }

        private IDictionary<int, IMatrix> ZerothOrderMatrix
        {
            get
            {
                if (zerothOrderMatrixFreeFree == null) BuildZerothOrderMatrixFreeFree();
                //else RebuildConductivityMatrices();
                return zerothOrderMatrixFreeFree;
            }
        }

        private void BuildZerothOrderMatrixFreeFree() => zerothOrderMatrixFreeFree = solver.BuildGlobalMatrices(zerothOrderMatrixProvider);

        private void BuildConductivitySubmatrices()
        {
            Dictionary<int, (IMatrix Cff, IMatrixView Cfc, IMatrixView Ccf, IMatrixView Ccc)> matrices =
                solver.BuildGlobalSubmatrices(zerothOrderMatrixProvider);

            zerothOrderMatrixFreeFree = new Dictionary<int, IMatrix>(model.Subdomains.Count);
            zerothOrderMatrixFreeConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            zerothOrderMatrixConstrFree = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            zerothOrderMatrixConstrConstr = new Dictionary<int, IMatrixView>(model.Subdomains.Count);
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int id = subdomain.ID;
                zerothOrderMatrixFreeFree.Add(id, matrices[id].Cff);
                zerothOrderMatrixFreeConstr.Add(id, matrices[id].Cfc);
                zerothOrderMatrixConstrFree.Add(id, matrices[id].Ccf);
                zerothOrderMatrixConstrConstr.Add(id, matrices[id].Ccc);
            }
        }

        private void RebuildZerothOrderMatrixFreeFree()
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
            if (mustRebuild) zerothOrderMatrixFreeFree = solver.BuildGlobalMatrices(zerothOrderMatrixProvider);
            foreach (ISubdomain subdomain in model.Subdomains) subdomain.ResetMaterialsModifiedProperty();
        }

        private void BuildFirstOrderMatrix() => firstOrderMatrix = solver.BuildGlobalMatrices(firstOrderMatrixProvider);

        #region IAnalyzerProvider Members
        public void ClearMatrices()
        {
            firstOrderMatrix = null;
            zerothOrderMatrixFreeFree = null;
            zerothOrderMatrixFreeConstr = null;
            zerothOrderMatrixConstrFree = null;
            zerothOrderMatrixConstrConstr = null;
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

            zerothOrderMatrixFreeFree = null;
            zerothOrderMatrixConstrFree = null;
            zerothOrderMatrixConstrConstr = null;
            firstOrderMatrix = null;
        }
        #endregion 

        #region IImplicitIntegrationProvider Members

        public IMatrixView LinearCombinationOfMatricesIntoStiffness(ImplicitIntegrationCoefficients coefficients,
            ISubdomain subdomain)
        {
            // The effective matrix should not overwrite the conductivity matrix. 
            // In a dynamic analysis that is not purely implicit we need the conductivity matrix.
            int id = subdomain.ID;
            return ZerothOrderMatrix[id].LinearCombination(coefficients.Stiffness, FirstOrderMatrix[id], coefficients.Mass);
        }

        public void ProcessRhs(ImplicitIntegrationCoefficients coefficients, ISubdomain subdomain, IVector rhs)
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

        public IVector MassMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.FirstOrderMatrix[subdomain.ID].Multiply(vector);

        //TODO: Ok this is weird. These methods should be named Second/First/ZeroOrderCoefficientTimesVector()
        public IVector DampingMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.ZerothOrderMatrix[subdomain.ID].Multiply(vector);

        #endregion

        #region IStaticProvider Members

        public IMatrixView CalculateMatrix(ISubdomain subdomain)
        {
            if (zerothOrderMatrixFreeFree == null) BuildZerothOrderMatrixFreeFree();
            return zerothOrderMatrixFreeFree[subdomain.ID];
        }

        public (IMatrixView matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr) CalculateSubMatrices(ISubdomain subdomain)
        {
            int id = subdomain.ID;
            if ((zerothOrderMatrixFreeFree == null) || (zerothOrderMatrixFreeConstr == null)
                || (zerothOrderMatrixConstrFree == null) || (zerothOrderMatrixConstrConstr == null))
            {
                BuildConductivitySubmatrices();
            }
            return (zerothOrderMatrixFreeFree[id], zerothOrderMatrixFreeConstr[id],
                zerothOrderMatrixConstrFree[id], zerothOrderMatrixConstrConstr[id]);
        }
        #endregion

        #region INonLinearProvider Members

        public double CalculateRhsNorm(IVectorView rhs) => rhs.Norm2();

        public void ProcessInternalRhs(ISubdomain subdomain, IVectorView solution, IVector rhs) { }

        #endregion
    }
}
