using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IConvectionDiffusionElement: IFiniteElement
    {
        IMatrix MassTransportConductivityMatrix(IElement element);
        IMatrix DiffusionConductivityMatrix(IElement element);
    }
}
