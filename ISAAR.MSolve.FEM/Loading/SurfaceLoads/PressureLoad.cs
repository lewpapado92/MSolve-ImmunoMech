using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Loading.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.FEM.Loading.SurfaceLoads
{
    public class PressureLoad:ISurfaceLoad
    {
        private readonly double _pressure;

        public PressureLoad(double pressure)
        {
            _pressure = pressure;
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

                var jacdet = (jacobianMatrix[0, 0] * jacobianMatrix[1, 1])
                                - (jacobianMatrix[1, 0] * jacobianMatrix[0, 1]);

                var weightFactor = integration.IntegrationPoints[gp].Weight;
                for (int indexNode = 0; indexNode < nodes.Count; indexNode++)
                {
                    var node = nodes[indexNode];
                    var valueX = _pressure * normalVector[0] * shapeFunctionNatural[gp][indexNode] * jacdet *
                                 weightFactor;
                    var valueY = _pressure * normalVector[1] * shapeFunctionNatural[gp][indexNode] * jacdet *
                                 weightFactor;
                    var valueZ = _pressure * normalVector[2] * shapeFunctionNatural[gp][indexNode] * jacdet *
                                 weightFactor;
                    if (loadTable.Contains(node,StructuralDof.TranslationX))
                    {
                        loadTable[node, StructuralDof.TranslationX] += valueX;
                    }
                    else
                    {
                        loadTable.TryAdd(node, StructuralDof.TranslationX, valueX);
                    }

                    if (loadTable.Contains(node,StructuralDof.TranslationY))
                    {
                        loadTable[node, StructuralDof.TranslationY] += valueY;
                    }
                    else
                    {
                        loadTable.TryAdd(node, StructuralDof.TranslationY, valueY);
                    }

                    if (loadTable.Contains(node,StructuralDof.TranslationZ))
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
