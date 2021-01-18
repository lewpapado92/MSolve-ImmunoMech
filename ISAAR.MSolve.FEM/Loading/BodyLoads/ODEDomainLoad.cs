using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Loading.BodyLoads
{
	using Entities;
	using Interfaces;
	using Interpolation;
	using ISAAR.MSolve.Discretization.Commons;
	using ISAAR.MSolve.Discretization.FreedomDegrees;
	using ISAAR.MSolve.Discretization.Integration.Quadratures;
	using ISAAR.MSolve.Discretization.Interfaces;
	using ISAAR.MSolve.FEM.Interpolation.Jacobians;
	using ISAAR.MSolve.LinearAlgebra.Matrices;
	using ISAAR.MSolve.LinearAlgebra.Vectors;
	using ISAAR.MSolve.Materials;

	public class ODEDomainLoad : IBodyLoad
	{
		private readonly ODEMaterial _material;
		protected double _load;
		private readonly IDofType _dofType;

		public ODEDomainLoad(ODEMaterial material, double load, IDofType dofType)
		{
			_material = material;
			_load = load;
			_dofType = dofType;
		}

		public Table<INode, IDofType, double> CalculateBodyLoad(IIsoparametricInterpolation3D interpolation, IQuadrature3D integration, IReadOnlyList<Node> nodes)
		{
			var loadTable = new Table<INode, IDofType, double>();
			IReadOnlyList<Matrix> shapeGradientsNatural =
				interpolation.EvaluateNaturalGradientsAtGaussPoints(integration);
			IReadOnlyList<double[]> shapeFunctionNatural =
				interpolation.EvaluateFunctionsAtGaussPoints(integration);

			for (int gp = 0; gp < integration.IntegrationPoints.Count; gp++)
			{
				var jacobian = new IsoparametricJacobian3D(nodes, shapeGradientsNatural[gp]);
				var jacdet = jacobian.DirectDeterminant;

				var weightFactor = integration.IntegrationPoints[gp].Weight;
				for (int indexNode = 0; indexNode < nodes.Count; indexNode++)
				{
					var node = nodes[indexNode];
					var valueX = _load * shapeFunctionNatural[gp][indexNode] * jacobian.DirectDeterminant *
								 weightFactor;
					if (loadTable.Contains(node, _dofType))
					{
						loadTable[node, _dofType] += valueX;
					}
					else
					{
						loadTable.TryAdd(node, _dofType, valueX);
					}
				}
			}

			return loadTable;
		}

        public Table<INode, IDofType, double> CalculateStabilizingBodyLoad(IIsoparametricInterpolation3D interpolation, IQuadrature3D integration, IReadOnlyList<Node> nodes)
        {
            throw new NotImplementedException();
        }
    }
}
