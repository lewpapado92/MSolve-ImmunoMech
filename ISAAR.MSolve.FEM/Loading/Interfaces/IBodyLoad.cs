namespace ISAAR.MSolve.FEM.Loading.Interfaces
{
	using System.Collections.Generic;
	using Entities;
	using Interpolation;
    using ISAAR.MSolve.Discretization.Commons;
    using ISAAR.MSolve.Discretization.FreedomDegrees;
    using ISAAR.MSolve.Discretization.Integration.Quadratures;
    using ISAAR.MSolve.Discretization.Interfaces;

    public interface IBodyLoad
	{
        Table<INode, IDofType, double> CalculateBodyLoad(IIsoparametricInterpolation3D interpolation,
            IQuadrature3D integration, IReadOnlyList<Node> nodes);
        Table<INode, IDofType, double> CalculateStabilizingBodyLoad(IIsoparametricInterpolation3D interpolation,
            IQuadrature3D integration, IReadOnlyList<Node> nodes);
    }
}
