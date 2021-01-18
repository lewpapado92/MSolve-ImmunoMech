using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IFiniteElement : IElementType
    {
        int ID { get; }
        ElementDimensions ElementDimensions { get; }
    }
}
