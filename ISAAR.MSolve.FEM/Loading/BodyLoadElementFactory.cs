using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Loading
{
	using Entities;
	using Interfaces;
	using Interpolation;
	using Interpolation.GaussPointExtrapolation;
    using ISAAR.MSolve.Discretization.Integration.Quadratures;
    using ISAAR.MSolve.Discretization.Mesh;

    public class BodyLoadElementFactory
	{
		private static readonly Dictionary<CellType, IIsoparametricInterpolation3D> interpolations;
		private static readonly Dictionary<CellType, IQuadrature3D> integrationsForLoad;
		private readonly IBodyLoad _bodyLoad;

		static BodyLoadElementFactory()
		{
			var interpolations = new Dictionary<CellType, IIsoparametricInterpolation3D>();
			var integrationsForLoad = new Dictionary<CellType, IQuadrature3D>();

			interpolations.Add(CellType.Tet4, InterpolationTet4.UniqueInstance);
			integrationsForLoad.Add(CellType.Tet4, TetrahedronQuadrature.Order1Point1);

			interpolations.Add(CellType.Tet10, InterpolationTet10.UniqueInstance);
			integrationsForLoad.Add(CellType.Tet10, TetrahedronQuadrature.Order2Points4);

			interpolations.Add(CellType.Hexa8, InterpolationHexa8.UniqueInstance);
			integrationsForLoad.Add(CellType.Hexa8, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2));

			interpolations.Add(CellType.Hexa20, InterpolationHexa20.UniqueInstance);
			integrationsForLoad.Add(CellType.Hexa20, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3));

			interpolations.Add(CellType.Hexa27, InterpolationHexa27.UniqueInstance);
			integrationsForLoad.Add(CellType.Hexa27, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3));

			interpolations.Add(CellType.Wedge6, InterpolationWedge6.UniqueInstance);
			integrationsForLoad.Add(CellType.Wedge6, WedgeQuadrature.Points6);

			interpolations.Add(CellType.Wedge15, InterpolationWedge15.UniqueInstance);
			integrationsForLoad.Add(CellType.Wedge15, WedgeQuadrature.Points8);

			interpolations.Add(CellType.Wedge18, InterpolationWedge18.UniqueInstance);
			integrationsForLoad.Add(CellType.Wedge18, WedgeQuadrature.Points8);

			interpolations.Add(CellType.Pyra5, InterpolationPyra5.UniqueInstance);
			integrationsForLoad.Add(CellType.Pyra5, PyramidQuadrature.Points5);

			interpolations.Add(CellType.Pyra13, InterpolationPyra13.UniqueInstance);
			integrationsForLoad.Add(CellType.Pyra13, PyramidQuadrature.Points6);

			interpolations.Add(CellType.Pyra14, InterpolationPyra14.UniqueInstance);
			integrationsForLoad.Add(CellType.Pyra14, PyramidQuadrature.Points6);

			BodyLoadElementFactory.interpolations = interpolations;
			BodyLoadElementFactory.integrationsForLoad = integrationsForLoad;
		}

		public BodyLoadElementFactory(IBodyLoad load, Model model)
		{
			_bodyLoad = load;

		}

		public BodyLoadElement CreateElement(CellType cellType, IReadOnlyList<Node> nodes) =>
			new BodyLoadElement(_bodyLoad, interpolations[cellType],
				integrationsForLoad[cellType], nodes);
	}
}
