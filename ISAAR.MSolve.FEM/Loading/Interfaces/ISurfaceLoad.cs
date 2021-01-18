using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;

namespace ISAAR.MSolve.FEM.Loading.Interfaces
{
    public interface ISurfaceLoad
    {
        Table<INode, IDofType, double> CalculateSurfaceLoad(IIsoparametricInterpolation2D interpolation,
            IQuadrature2D integration, IReadOnlyList<Node> nodes);
    }
}