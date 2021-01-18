using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.FEM.Providers
{
    public class ElementDiffusionConductivityProvider : IElementMatrixProvider
    {
        public IMatrix Matrix(IElement element)
        {
            IConvectionDiffusionElement elementType = (IConvectionDiffusionElement)element.ElementType;
            return elementType.DiffusionConductivityMatrix(element);
        }
    }
}
