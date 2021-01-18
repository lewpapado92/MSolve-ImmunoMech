using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IConvectionDiffusionIntegrationProvider : IAnalyzerProvider
    {
        //TODO: This should not exist at all. The provider should return the 0th order (stiffness), 1st order (damping) and 2nd
        //      order matrices (or some matrix representations that can be combined between them and multiplied with vectors).
        IMatrixView LinearCombinationOfMatricesIntoStiffness(ImplicitIntegrationCoefficients coefficients, ISubdomain subdomain);

        //TODO: Way too generic name. Probably needs refactoring as well.
        void ProcessRhs(ISubdomain subdomain, IVector rhs);

        //TODO: I think the analyzer is responsible for these. E.g. Newmark has the formulas with beta and gamma, Euler has the
        //      central differences formulas, etc.
        IDictionary<int, IVector> GetAccelerationsOfTimeStep(int timeStep);
        IDictionary<int, IVector> GetVelocitiesOfTimeStep(int timeStep);

        IDictionary<int, IVector> GetRhsFromHistoryLoad(int timeStep);

        //TODO: what about thermal? There is no mass matrix there. Define these as 1st order matrix coeff, 2nd order ...
        IVector ConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector);
        IVector CapacityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector);
        IVector DiffusionConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector);
        IVector MassTransportConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector);
        IVector StabilizingConductivityMatrixVectorProduct(ISubdomain subdomain, IVectorView vector);
        IVector StabilizingRhs(ISubdomain subdomain);
    }
}
