namespace ISAAR.MSolve.Materials.Interfaces
{
	public interface IDynamicMaterial
	{
		double Density { get; }
		double RayleighCoeffMass { get; }
		double RayleighCoeffStiffness { get; }
	}
}
