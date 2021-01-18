using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISSAR.MSolve.Discretization.Loads
{
	public class MassAccelerationLoad
	{
		public IDofType DOF { get; set; }
		public double Amount { get; set; }
	}
}
