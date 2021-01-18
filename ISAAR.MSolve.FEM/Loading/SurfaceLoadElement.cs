using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Loading.Interfaces;

namespace ISAAR.MSolve.FEM.Loading
{
    public class SurfaceLoadElement:ISurfaceLoadElement
    {
        private readonly ISurfaceLoad surfaceLoad;
        private readonly IIsoparametricInterpolation2D isoparametricInterpolation2D;
        private readonly IQuadrature2D quadrature2D;
        private readonly IReadOnlyList<Node> nodes;

        public SurfaceLoadElement(ISurfaceLoad surfaceLoad, IIsoparametricInterpolation2D isoparametricInterpolation2D,
            IQuadrature2D quadrature2D, IReadOnlyList<Node> nodes)
        {
            this.surfaceLoad = surfaceLoad;
            this.isoparametricInterpolation2D = isoparametricInterpolation2D;
            this.quadrature2D = quadrature2D;
            this.nodes = nodes;
        }

        public Table<INode, IDofType, double> CalculateSurfaceLoad() =>
            surfaceLoad.CalculateSurfaceLoad(isoparametricInterpolation2D, quadrature2D, nodes);

    }
}
