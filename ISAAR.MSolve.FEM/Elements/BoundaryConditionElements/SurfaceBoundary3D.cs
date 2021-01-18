using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.FEM.Interpolation.Jacobians;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISSAR.MSolve.Discretization.Loads;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Elements.BoundaryConditionElements
{
    public class SurfaceBoundary3D : IConvectionDiffusionElement,  ICell<Node>
    {
        private readonly IDofType[][] dofTypes; //TODO: this should not be stored for each element. Instead store it once for each Quad4, Tri3, etc. Otherwise create it on the fly.
        private readonly ConvectionDiffusionMaterial material;
        //private readonly Dictionary<GaussPoint, ThermalMaterial> materialsAtGaussPoints;


        public SurfaceBoundary3D(double thickness, IReadOnlyList<Node> nodes, IIsoparametricInterpolation2D interpolation,
            IQuadrature2D quadratureForStiffness, IQuadrature2D quadratureForConsistentMass,
            IGaussPointExtrapolation2D gaussPointExtrapolation,
            ConvectionDiffusionMaterial material)
        {
            this.material = material;
            this.GaussPointExtrapolation = gaussPointExtrapolation;
            this.Nodes = nodes;
            this.Interpolation = interpolation;
            this.QuadratureForConsistentMass = quadratureForConsistentMass;
            this.QuadratureForStiffness = quadratureForStiffness;
            this.Thickness = thickness;

            dofTypes = new IDofType[nodes.Count][];
            for (int i = 0; i < interpolation.NumFunctions; ++i) dofTypes[i] = new IDofType[] { ThermalDof.Temperature };
        }

        public CellType CellType => Interpolation.CellType;
        public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

        public int ID => throw new NotImplementedException(
            "Element type codes should be in a settings class. Even then it's a bad design choice");

        public IGaussPointExtrapolation2D GaussPointExtrapolation { get; }
        public IIsoparametricInterpolation2D Interpolation { get; }
        public IReadOnlyList<Node> Nodes { get; }
        public IQuadrature2D QuadratureForConsistentMass { get; }
        public IQuadrature2D QuadratureForStiffness { get; }
        public double Thickness { get; }

        public bool MaterialModified => throw new NotImplementedException();

        public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();
        //int IElement.ID { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        //public IElementType ElementType => throw new NotImplementedException();

        //IReadOnlyList<INode> IElement.Nodes => throw new NotImplementedException();

        //public ISubdomain Subdomain => throw new NotImplementedException();

        public IMatrix MassMatrix(IElement element)
        {
            return BuildCapacityMatrix();
        }

        public Matrix BuildCapacityMatrix()
        {
            int numDofs = Nodes.Count;
            var capacity = Matrix.CreateZero(numDofs, numDofs);
            //IReadOnlyList<double[]> shapeFunctions =
            //    Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
            //IReadOnlyList<Matrix> shapeGradientsNatural =
            //    Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            //for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
            //{
            //    Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
            //    Matrix partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
            //    var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
            //    double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
            //    capacity.AxpyIntoThis(partial, dA);
            //}

            ////WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            //capacity.ScaleIntoThis(Thickness * material.Density * material.SpecialHeatCoeff);
            return capacity;
        }
        public Matrix BuildDiffusionMatrix()
        {
            int numDofs = Nodes.Count;
            double kappaCoef;
            if (numDofs == 4) kappaCoef = 100;
            else kappaCoef = 1;
            var stiffness = Matrix.CreateZero(numDofs, numDofs);
            IReadOnlyList<double[]> shapeFunctions =
                Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            double[] dist = new double[Nodes.Count];
            List<INode> neighborNodes = new List<INode>();
            for (int i = 0; i < Nodes.Count; i++)
            {
                neighborNodes.Clear();
                foreach (var element in Nodes[i].ElementsDictionary.Values)
                {
                    neighborNodes.AddRange(element.Nodes);
                }
                neighborNodes = neighborNodes.Distinct().ToList();
                neighborNodes.Remove(Nodes[i]);
                var minDist = neighborNodes.Select(x => Math.Sqrt(
                      Math.Pow(Nodes[i].X - x.X, 2) +
                              Math.Pow(Nodes[i].Y - x.Y, 2) + Math.Pow(Nodes[i].Z - x.Z, 2))).Min();
                dist[i] = minDist;
            }
            //double kappa = kappaCoef * material.DiffusionCoeff/ dist.Min();
            double kappa = material.DiffusionCoeff / .05;

            for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
            {
                //var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                Vector shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
                Matrix jacobianMatrix = Matrix.CreateZero(2, 3);
                //Matrix jacobian = Matrix.CreateZero(2, 2);
                for (int k = 0; k < this.Nodes.Count; k++)
                {
                    jacobianMatrix[0, 0] += shapeGradientsNatural[gp][k, 0] * this.Nodes[k].X;
                    jacobianMatrix[0, 1] += shapeGradientsNatural[gp][k, 0] * this.Nodes[k].Y;
                    //jacobian[0, 0] += shapeGradientsNatural[gp][k, 0] * this.Nodes[k].X;
                    //jacobian[0, 1] += shapeGradientsNatural[gp][k, 0] * this.Nodes[k].Y;
                    jacobianMatrix[0, 2] += shapeGradientsNatural[gp][k, 0] * this.Nodes[k].Z;
                    jacobianMatrix[1, 0] += shapeGradientsNatural[gp][k, 1] * this.Nodes[k].X;
                    jacobianMatrix[1, 1] += shapeGradientsNatural[gp][k, 1] * this.Nodes[k].Y;
                    //jacobian[1, 0] += shapeGradientsNatural[gp][k, 1] * this.Nodes[k].X;
                    //jacobian[1, 1] += shapeGradientsNatural[gp][k, 1] * this.Nodes[k].Y;
                    jacobianMatrix[1, 2] += shapeGradientsNatural[gp][k, 1] * this.Nodes[k].Z;
                }
                Matrix jacobianMatrixLeftInverse = jacobianMatrix.Transpose() * (jacobianMatrix * jacobianMatrix.Transpose()).Invert();
                //jacobian = jacobianMatrix * jacobianMatrix.Transpose();
                //jacobian.ScaleIntoThis(20);
                //var jacinv = jacobian.Invert();
                Vector tangentVector1 = jacobianMatrix.GetRow(0);
                Vector tangentVector2 = jacobianMatrix.GetRow(1);
                Vector normalVector = tangentVector1.CrossProduct(tangentVector2);
                var jacdet = normalVector.Norm2();
                normalVector.ScaleIntoThis(1 / jacdet);
                //Matrix shapeGradientsCartesian1 = (shapeGradientsNatural[gp]) * jacobian.Invert();
                Matrix shapeGradientsCartesian = (shapeGradientsNatural[gp]) * jacobianMatrixLeftInverse.Transpose();
                Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);
                Vector deformationNormal = normalVector * deformation;
                //Matrix partial = -2 * (shapeFunctionMatrix.TensorProduct(deformationNormal))
                //    + kappa * shapeFunctionMatrix.TensorProduct(shapeFunctionMatrix);
                Matrix partial = kappa * shapeFunctionMatrix.TensorProduct(shapeFunctionMatrix)
                    - 2* material.DiffusionCoeff*(shapeFunctionMatrix.TensorProduct(deformationNormal));

                //Vector surfaceBasisVector1 = Vector.CreateZero(3);
                //surfaceBasisVector1[0] = jacobianMatrix[0, 0];
                //surfaceBasisVector1[1] = jacobianMatrix[0, 1];
                //surfaceBasisVector1[2] = jacobianMatrix[0, 2];

                //Vector surfaceBasisVector2 = Vector.CreateZero(3);
                //surfaceBasisVector2[0] = jacobianMatrix[1, 0];
                //surfaceBasisVector2[1] = jacobianMatrix[1, 1];
                //surfaceBasisVector2[2] = jacobianMatrix[1, 2];

                //Vector surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);

                double dA = jacdet * QuadratureForStiffness.IntegrationPoints[gp].Weight;
                stiffness.AxpyIntoThis(partial, dA);
            }

            //WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            //stiffness.ScaleIntoThis(1);
            return stiffness;
        }
        public Matrix BuildZeroMatrix()
        {
            int numDofs = Nodes.Count;
            var value = Matrix.CreateZero(numDofs, numDofs);
            return value;
        }

        private Matrix BuildDeformationMatrix(Matrix shapeGradientsCartesian)
        {
            //TODO: isn't this just the transpose of [dNi/dxj]?
            var deformation = Matrix.CreateZero(3, Nodes.Count);
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                deformation[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
                deformation[1, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
                deformation[2, nodeIdx] = shapeGradientsCartesian[nodeIdx, 2];
            }
            return deformation;
        }

        /// <summary>
        /// The shape function matrix is 1-by-n, where n = is the number of shape functions.
        /// </summary>
        public Vector BuildShapeFunctionMatrix(double[] shapeFunctions) //TODO: reconsider this. As it is, it just returns the shape functions in a Matrix
        {
            //var array2D = new double[1, shapeFunctions.Length];
            //for (int i = 0; i < shapeFunctions.Length; ++i)
            //{
            //    array2D[0, i] = shapeFunctions[i];
            //}
            //return Matrix.CreateFromArray(array2D);
            return Vector.CreateFromArray(shapeFunctions);
        }

        public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofTypes;

        public void ResetMaterialModified()
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForces(IElement element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
        {
            throw new NotImplementedException();
        }

        public void SaveMaterialState()
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialState()
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialStresses()
        {
            throw new NotImplementedException();
        }

        public IMatrix StiffnessMatrix(IElement element)
        {
            return BuildDiffusionMatrix();//D*u
        }
        public IMatrix DiffusionConductivityMatrix(IElement element)
        {
            return BuildDiffusionMatrix();//D*u
        }
        public IMatrix MassTransportConductivityMatrix(IElement element)
        {
            return BuildZeroMatrix();
        }

        //public IMatrix RHSFLuxMatrix(IElement element)
        //{
        //    return BuildRHSFluxMatrix();//F*flux
        //}
        //public IMatrix RHSPrescribedMatrix(IElement element)
        //{
        //    return BuildRHSPrescribedMatrix();//P*u0
        //}

        public IMatrix DampingMatrix(IElement element)
        {
            return BuildZeroMatrix();
        }
    }
}
