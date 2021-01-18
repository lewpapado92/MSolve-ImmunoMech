using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISSAR.MSolve.Discretization.Loads
{
	public interface ITimeDependentNodalLoad
	{
		INode Node { get; set; }
		IDofType DOF { get; set; }

		double GetLoadAmount(int timeStep);
	}
}
