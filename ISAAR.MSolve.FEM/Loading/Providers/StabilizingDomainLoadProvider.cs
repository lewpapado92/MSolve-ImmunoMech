using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Loading.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.FEM.Loading.Providers
{
    public class StabilizingDomainLoadProvider
    {
        private IBodyLoadElement element;
        //public StabilizingDomainLoadProvider(IBodyLoadElement element)
        //{
        //    this.element = element;
        //}
        public Table<INode, IDofType, double> StabilizingLoad => element.CalculateStabilizingBodyLoad();

    }
}
