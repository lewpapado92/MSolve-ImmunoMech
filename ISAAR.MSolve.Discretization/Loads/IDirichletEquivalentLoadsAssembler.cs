using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISSAR.MSolve.Discretization.Loads
{
	public interface IDirichletEquivalentLoadsAssembler
	{
		IVector GetEquivalentNodalLoads(ISubdomain subdomain, IVectorView solution, double constraintScalingFactor);
	}
}
