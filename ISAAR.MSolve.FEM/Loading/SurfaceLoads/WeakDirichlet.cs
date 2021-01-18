using DotNumerics.ODE.DVode;
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
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Loading.SurfaceLoads
{
    public class WeakDirichlet : ISurfaceLoad
    {
        private readonly double _magnitude;
        private readonly double _diffusionCoeff;
        private readonly DirichletDistribution _distribution;

        public delegate Vector DirichletDistribution(IReadOnlyList<Node> list);
        public WeakDirichlet(DirichletDistribution distribution, double diffusionCoeff)
        {
            _distribution = distribution;
            _diffusionCoeff = diffusionCoeff;
        }

        public Table<INode, IDofType, double> CalculateSurfaceLoad(IIsoparametricInterpolation2D interpolation, IQuadrature2D integration, IReadOnlyList<Node> nodes)
        {
            //nodes[0].ElementsDictionary.Values.ToList<>
            int numDofs = nodes.Count;
            double kappaCoef;
            if (numDofs == 4) kappaCoef = 100;
            else kappaCoef = 1;
            var stiffness = Matrix.CreateZero(numDofs, numDofs);
            IReadOnlyList<double[]> shapeFunctions =
                interpolation.EvaluateFunctionsAtGaussPoints(integration);
            IReadOnlyList<Matrix> shapeGradientsNatural =
               interpolation.EvaluateNaturalGradientsAtGaussPoints(integration);

            double[] dist = new double[nodes.Count];
            List<INode> neighborNodes = new List<INode>();
            for (int i = 0; i < nodes.Count; i++)
            {
                neighborNodes.Clear();
                foreach (var element in nodes[i].ElementsDictionary.Values)
                {
                    neighborNodes.AddRange(element.Nodes);
                }
                neighborNodes = neighborNodes.Distinct().ToList();
                neighborNodes.Remove(nodes[i]);
                var minDist = neighborNodes.Select(x => Math.Sqrt(
                      Math.Pow(nodes[i].X - x.X, 2) +
                              Math.Pow(nodes[i].Y - x.Y, 2) + Math.Pow(nodes[i].Z - x.Z, 2))).Min();
                dist[i] = minDist;
            }
            //double kappa = kappaCoef * _diffusionCoeff / dist.Min();
            //double[] kappa = new double[nodes.Count];
            for (int gp = 0; gp < integration.IntegrationPoints.Count; ++gp)
            {
                //var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                
                var shapeFunctionMatrix = Vector.CreateFromArray(shapeFunctions[gp]);
                Matrix jacobianMatrix = Matrix.CreateZero(2, 3);
                for (int k = 0; k < nodes.Count; k++)
                {
                    //xGaussPoint += shapeFunctionMatrix[gp, k] * this.Nodes[k].X;
                    //yGaussPoint += shapeFunctionMatrix[gp, k] * this.Nodes[k].Y;
                    //zGaussPoint += shapeFunctionMatrix[gp, k] * this.Nodes[k].Z;
                    jacobianMatrix[0, 0] += shapeGradientsNatural[gp][k, 0] * nodes[k].X;
                    jacobianMatrix[0, 1] += shapeGradientsNatural[gp][k, 0] * nodes[k].Y;
                    jacobianMatrix[0, 2] += shapeGradientsNatural[gp][k, 0] * nodes[k].Z;
                    jacobianMatrix[1, 0] += shapeGradientsNatural[gp][k, 1] * nodes[k].X;
                    jacobianMatrix[1, 1] += shapeGradientsNatural[gp][k, 1] * nodes[k].Y;
                    jacobianMatrix[1, 2] += shapeGradientsNatural[gp][k, 1] * nodes[k].Z;
                    //kappa[k] = _diffusionCoeff / .05 / dist[k];
                }
                Vector tangentVector1 = jacobianMatrix.GetRow(0);
                Vector tangentVector2 = jacobianMatrix.GetRow(1);
                Vector normalVector = tangentVector1.CrossProduct(tangentVector2);
                var jacdet = normalVector.Norm2(); 
                normalVector.ScaleIntoThis(1/jacdet);
                Matrix jacobianMatrixLeftInverse = jacobianMatrix.Transpose() *
                    (jacobianMatrix * jacobianMatrix.Transpose()).Invert();
                Matrix shapeGradientsCartesian = (shapeGradientsNatural[gp]) * jacobianMatrixLeftInverse.Transpose();
                Matrix deformation = shapeGradientsCartesian.Transpose();
                Vector deformationNormal = normalVector * deformation;
                //Vector deformationX = deformation.GetRow(0);
                //Vector deformationY = deformation.GetRow(1);
                //Vector deformationZ = deformation.GetRow(2);
                //Matrix partial = kappa * shapeFunctionMatrix.TensorProduct(shapeFunctionMatrix)
                //    -_diffusionCoeff * deformationNormal.TensorProduct(shapeFunctionMatrix);
                double kappa = _diffusionCoeff / .05;
                Matrix partial = kappa * shapeFunctionMatrix.TensorProduct(shapeFunctionMatrix);

                //Vector surfaceBasisVector1 = Vector.CreateZero(3);
                //surfaceBasisVector1[0] = jacobianMatrix[0, 0];
                //surfaceBasisVector1[1] = jacobianMatrix[0, 1];
                //surfaceBasisVector1[2] = jacobianMatrix[0, 2];

                //Vector surfaceBasisVector2 = Vector.CreateZero(3);
                //surfaceBasisVector2[0] = jacobianMatrix[1, 0];
                //surfaceBasisVector2[1] = jacobianMatrix[1, 1];
                //surfaceBasisVector2[2] = jacobianMatrix[1, 2];

                //Vector surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);
                //var jacdet = surfaceBasisVector3.Norm2();

                double dA = jacdet * integration.IntegrationPoints[gp].Weight;
                stiffness.AxpyIntoThis(partial, dA);
            }
            //var appliedDisplacements=Vector.CreateWithValue(nodes.Count, _magnitude);
            var appliedDisplacements = _distribution(nodes);
            var weakDirichletForces = stiffness * appliedDisplacements;

            var table = new Table<INode, IDofType, double>();
            for (int i = 0; i < nodes.Count; i++)
            {
                table.TryAdd(nodes[i], ThermalDof.Temperature, weakDirichletForces[i]);
            }
            return table;
        }
    }

}
