using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Loading.Interfaces;

namespace ISAAR.MSolve.FEM.Loading
{
    public class SurfaceLoadElementFactory
    {
        private static readonly IReadOnlyDictionary<CellType, IQuadrature2D> integrationForLoad;
        private static readonly IReadOnlyDictionary<CellType, IIsoparametricInterpolation2D> interpolations;
        private readonly ISurfaceLoad _surfaceLoad;

        static SurfaceLoadElementFactory()
        {
            var interpolations = new Dictionary<CellType, IIsoparametricInterpolation2D>();
            var integrationForLoad=new Dictionary<CellType, IQuadrature2D>();

            interpolations.Add(CellType.Quad4, InterpolationQuad4.UniqueInstance);
            integrationForLoad.Add(CellType.Quad4, GaussLegendre2D.GetQuadratureWithOrder(2,2));

            interpolations.Add(CellType.Quad8, InterpolationQuad8.UniqueInstance);
            integrationForLoad.Add(CellType.Quad8, GaussLegendre2D.GetQuadratureWithOrder(3,3));

            interpolations.Add(CellType.Tri3, InterpolationTri3.UniqueInstance);
            integrationForLoad.Add(CellType.Tri3, TriangleQuadratureSymmetricGaussian.Order1Point1);

            interpolations.Add(CellType.Tri6, InterpolationTri6.UniqueInstance);
            integrationForLoad.Add(CellType.Tri6, TriangleQuadratureSymmetricGaussian.Order2Points3);

            SurfaceLoadElementFactory.integrationForLoad = integrationForLoad;
            SurfaceLoadElementFactory.interpolations = interpolations;
        }

        public SurfaceLoadElementFactory(ISurfaceLoad surfaceLoad)
        {
            _surfaceLoad = surfaceLoad;
        }

        public SurfaceLoadElement CreateElement(CellType cellType, IReadOnlyList<Node> nodes)
        {
            return new SurfaceLoadElement(_surfaceLoad,interpolations[cellType], 
                integrationForLoad[cellType],nodes);
        }
    }
}
