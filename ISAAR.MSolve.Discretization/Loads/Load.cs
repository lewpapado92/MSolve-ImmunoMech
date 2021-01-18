using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISSAR.MSolve.Discretization.Loads
{
	public class Load
	{
		public INode Node { get; set; }
		public IDofType DOF { get; set; }
		public double Amount { get; set; }
	}
}
