using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.FEM.Providers
{
    public class ElementMassTransportConductivityProvider : IElementMatrixProvider
    {
        public IMatrix Matrix(IElement element)
        {
            IConvectionDiffusionElement elementType = (IConvectionDiffusionElement)element.ElementType;
            return elementType.MassTransportConductivityMatrix(element);
        }
    }
}
