using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Loading
{
	using Entities;
	using Interfaces;
	using Interpolation;
    using ISAAR.MSolve.Discretization.Commons;
    using ISAAR.MSolve.Discretization.FreedomDegrees;
    using ISAAR.MSolve.Discretization.Integration.Quadratures;
    using ISAAR.MSolve.Discretization.Interfaces;
	public class BodyLoadElement : IBodyLoadElement
	{
		private readonly IBodyLoad _bodyLoad;
		private readonly IIsoparametricInterpolation3D _interpolation;
		private readonly IQuadrature3D _quadrature;
		private readonly IReadOnlyList<Node> _nodes;

		public BodyLoadElement(IBodyLoad bodyLoad, IIsoparametricInterpolation3D interpolation3D,
			IQuadrature3D quadrature3D, IReadOnlyList<Node> nodes)
		{
			_bodyLoad = bodyLoad;
			_interpolation = interpolation3D;
			_quadrature = quadrature3D;
			_nodes = nodes;
		}

		public Table<INode, IDofType, double> CalculateBodyLoad() =>
			_bodyLoad.CalculateBodyLoad(_interpolation, _quadrature, _nodes);

		public Table<INode, IDofType, double> CalculateStabilizingBodyLoad() =>
			_bodyLoad.CalculateStabilizingBodyLoad(_interpolation, _quadrature, _nodes);
	}
}
