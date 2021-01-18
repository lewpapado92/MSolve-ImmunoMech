using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.FEM.Loading.Interfaces
{
    public interface ISurfaceLoadElement
    {
        Table<INode, IDofType, double> CalculateSurfaceLoad();
    }
}