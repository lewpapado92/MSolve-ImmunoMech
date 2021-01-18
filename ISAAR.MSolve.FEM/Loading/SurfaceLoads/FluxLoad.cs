using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Loading.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;
using System.Collections.Generic;
using System.Reflection.Emit;
using System.Text;

namespace ISAAR.MSolve.FEM.Loading.SurfaceLoads
{
    public class FluxLoad : ISurfaceLoad
    {
        private readonly double _flux;

        public FluxLoad(double flux)
        {
            _flux = flux;
        }

        public Table<INode, IDofType, double> CalculateSurfaceLoad(IIsoparametricInterpolation2D interpolation, IQuadrature2D integration, IReadOnlyList<Node> nodes)
        {
            int numDofs = nodes.Count;
            var stiffness = Vector.CreateZero(numDofs);
            IReadOnlyList<double[]> shapeFunctions =
                interpolation.EvaluateFunctionsAtGaussPoints(integration);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                interpolation.EvaluateNaturalGradientsAtGaussPoints(integration);

            for (int gp = 0; gp < integration.IntegrationPoints.Count; ++gp)
            {
                Vector shapeFunctionVector =Vector.CreateFromArray(shapeFunctions[gp]);
                Matrix jacobianMatrix = Matrix.CreateZero(2, 3);
                for (int k = 0; k < nodes.Count; k++)
                {
                    jacobianMatrix[0, 0] += shapeGradientsNatural[gp][k, 0] * nodes[k].X;
                    jacobianMatrix[0, 1] += shapeGradientsNatural[gp][k, 0] * nodes[k].Y;
                    jacobianMatrix[0, 2] += shapeGradientsNatural[gp][k, 0] * nodes[k].Z;
                    jacobianMatrix[1, 0] += shapeGradientsNatural[gp][k, 1] * nodes[k].X;
                    jacobianMatrix[1, 1] += shapeGradientsNatural[gp][k, 1] * nodes[k].Y;
                    jacobianMatrix[1, 2] += shapeGradientsNatural[gp][k, 1] * nodes[k].Z;
                }
                Vector tangentVector1 = jacobianMatrix.GetRow(0);
                Vector tangentVector2 = jacobianMatrix.GetRow(1);
                Vector normalVector = tangentVector1.CrossProduct(tangentVector2);

                var jacdet = normalVector.Norm2();

                double dA = jacdet * integration.IntegrationPoints[gp].Weight;
                stiffness.AxpyIntoThis(shapeFunctionVector, dA);
            }
            var appliedDisplacements = Vector.CreateWithValue(nodes.Count, _flux);
            var fluxLoad = stiffness.DoEntrywise(appliedDisplacements, (stiffi, appi) => stiffi*appi);
            var table = new Table<INode, IDofType, double>();
            for (int i = 0; i < nodes.Count; i++)
            {
                table.TryAdd(nodes[i], ThermalDof.Temperature, fluxLoad[i]);
            }
            return table;
        }
    }
}
