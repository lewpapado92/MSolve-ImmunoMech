using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.Jacobians;
using ISAAR.MSolve.FEM.Loading.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.FEM.Loading.SurfaceLoads
{
    public class DistributedLoad:ISurfaceLoad
    {
        private readonly double _loadX;
        private readonly double _loadY;
        private readonly double _loadZ;

        public DistributedLoad(double loadX, double loadY, double loadZ)
        {
            _loadX = loadX;
            _loadY = loadY;
            _loadZ = loadZ;
        }

        public Table<INode, IDofType, double> CalculateSurfaceLoad(IIsoparametricInterpolation2D interpolation,
            IQuadrature2D integration, IReadOnlyList<Node> nodes)
        {
            var loadTable = new Table<INode, IDofType, double>();
            IReadOnlyList<Matrix> shapeGradientsNatural =
                interpolation.EvaluateNaturalGradientsAtGaussPoints(integration);
            IReadOnlyList<double[]> shapeFunctionNatural =
                interpolation.EvaluateFunctionsAtGaussPoints(integration);
			for (int gp = 0; gp < integration.IntegrationPoints.Count; gp++)
			{
				var jacobianMatrix = Matrix.CreateZero(2, 3);
				for (int indexNode = 0; indexNode < nodes.Count; indexNode++)
				{
					jacobianMatrix[0, 0] += shapeGradientsNatural[gp][indexNode, 0] * nodes[indexNode].X;
					jacobianMatrix[0, 1] += shapeGradientsNatural[gp][indexNode, 0] * nodes[indexNode].Y;
					jacobianMatrix[0, 2] += shapeGradientsNatural[gp][indexNode, 0] * nodes[indexNode].Z;

					jacobianMatrix[1, 0] += shapeGradientsNatural[gp][indexNode, 1] * nodes[indexNode].X;
					jacobianMatrix[1, 1] += shapeGradientsNatural[gp][indexNode, 1] * nodes[indexNode].Y;
					jacobianMatrix[1, 2] += shapeGradientsNatural[gp][indexNode, 1] * nodes[indexNode].Z;
				}

				var tangentVector1 = jacobianMatrix.GetRow(0);
				var tangentVector2 = jacobianMatrix.GetRow(1);
				var normalVector = tangentVector1.CrossProduct(tangentVector2);

				Vector surfaceBasisVector1 = Vector.CreateZero(3);
				surfaceBasisVector1[0] = jacobianMatrix[0, 0];
				surfaceBasisVector1[1] = jacobianMatrix[0, 1];
				surfaceBasisVector1[2] = jacobianMatrix[0, 2];

				Vector surfaceBasisVector2 = Vector.CreateZero(3);
				surfaceBasisVector2[0] = jacobianMatrix[1, 0];
				surfaceBasisVector2[1] = jacobianMatrix[1, 1];
				surfaceBasisVector2[2] = jacobianMatrix[1, 2];

				Vector surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);
				var jacdet = surfaceBasisVector3.Norm2();

				var weightFactor = integration.IntegrationPoints[gp].Weight;
				for (int indexNode = 0; indexNode < nodes.Count; indexNode++)
				{
					var node = nodes[indexNode];
					var valueX = _loadX * shapeFunctionNatural[gp][indexNode] * jacdet *
								 weightFactor;
					var valueY = _loadY * shapeFunctionNatural[gp][indexNode] * jacdet *
								 weightFactor;
					var valueZ = _loadZ * shapeFunctionNatural[gp][indexNode] * jacdet *
								 weightFactor;
					if (loadTable.Contains(node, StructuralDof.TranslationX))
					{
						loadTable[node, StructuralDof.TranslationX] += valueX;
					}
					else
					{
						loadTable.TryAdd(node, StructuralDof.TranslationX, valueX);
					}

					if (loadTable.Contains(node, StructuralDof.TranslationY))
					{
						loadTable[node, StructuralDof.TranslationY] += valueY;
					}
					else
					{
						loadTable.TryAdd(node, StructuralDof.TranslationY, valueY);
					}

					if (loadTable.Contains(node, StructuralDof.TranslationZ))
					{
						loadTable[node, StructuralDof.TranslationZ] += valueZ;
					}
					else
					{
						loadTable.TryAdd(node, StructuralDof.TranslationZ, valueZ);
					}
				}
			}

			return loadTable;
        }
    }
}
