using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IConvectionDiffusionBoundaryElement : IFiniteElement
    {
        IMatrix RHSPrescribedMatrix(IElement element);
        IMatrix RHSFLuxMatrix(IElement element);
    }
}
