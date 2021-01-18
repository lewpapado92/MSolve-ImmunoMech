using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Loading.Interfaces
{
    public interface IBodyLoadElement
    {
        Table<INode, IDofType, double> CalculateBodyLoad();
        Table<INode, IDofType, double> CalculateStabilizingBodyLoad();
    }
}
