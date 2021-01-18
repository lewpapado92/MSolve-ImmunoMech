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

	public class ConvectionDiffusionDomainLoad : IBodyLoad
	{
		private readonly ConvectionDiffusionMaterial _material;
		protected double _load;
		private readonly IDofType _dofType;

		public ConvectionDiffusionDomainLoad(ConvectionDiffusionMaterial material, double load, IDofType dofType)
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
			var loadTable = new Table<INode, IDofType, double>();
			IReadOnlyList<Matrix> shapeGradientsNatural =
				interpolation.EvaluateNaturalGradientsAtGaussPoints(integration);
			IReadOnlyList<double[]> shapeFunctionNatural =
				interpolation.EvaluateFunctionsAtGaussPoints(integration);

			for (int gp = 0; gp < integration.IntegrationPoints.Count; gp++)
			{
				var jacobian = new IsoparametricJacobian3D(nodes, shapeGradientsNatural[gp]);
				Matrix shapeGradientsCartesian =
				   jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
					//TODO: isn't this just the transpose of [dNi/dxj]?
				var deformation = Matrix.CreateZero(3, nodes.Count);
				for (int nodeIdx = 0; nodeIdx < nodes.Count; ++nodeIdx)
				{
					deformation[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
					deformation[1, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
					deformation[2, nodeIdx] = shapeGradientsCartesian[nodeIdx, 2];
				}
				Vector deformationX = deformation.GetRow(0);
				Vector deformationY = deformation.GetRow(1);
				Vector deformationZ = deformation.GetRow(2);
				Vector partialK = deformationX.Scale(_material.ConvectionCoeff[0]) +
					deformationY.Scale(_material.ConvectionCoeff[1]) +
					deformationZ.Scale(_material.ConvectionCoeff[2]);
				//loadTable.AxpyIntoThis(partialK, dA);

				var weightFactor = integration.IntegrationPoints[gp].Weight;
				for (int indexNode = 0; indexNode < nodes.Count; indexNode++)
				{
					var node = nodes[indexNode];
					var valueX = -0.5 * _load * partialK[indexNode] * jacobian.DirectDeterminant *
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
	}
}
