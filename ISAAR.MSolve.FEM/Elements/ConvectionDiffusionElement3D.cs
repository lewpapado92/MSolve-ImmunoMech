using System;
using System.Collections.Generic;
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

namespace ISAAR.MSolve.FEM.Elements
{
    public class ConvectionDiffusionElement3D : IConvectionDiffusionElement, ICell<Node>
    {
        private readonly static IDofType[] nodalDOFTypes = new IDofType[] { ThermalDof.Temperature };
        private readonly IDofType[][] dofTypes; //TODO: this should not be stored for each element. Instead store it once for each Quad4, Tri3, etc. Otherwise create it on the fly.
        private readonly ConvectionDiffusionMaterial material;

        public ConvectionDiffusionElement3D(IReadOnlyList<Node> nodes, IIsoparametricInterpolation3D interpolation,
        IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
        IGaussPointExtrapolation3D gaussPointExtrapolation, ConvectionDiffusionMaterial material)
        {
            this.material = material;
            this.GaussPointExtrapolation = gaussPointExtrapolation;
            this.Nodes = nodes;
            this.Interpolation = interpolation;
            this.QuadratureForConsistentMass = quadratureForMass;
            this.QuadratureForStiffness = quadratureForStiffness;

            dofTypes = new IDofType[nodes.Count][];
            for (int i = 0; i < interpolation.NumFunctions; ++i) dofTypes[i] = new IDofType[] { ThermalDof.Temperature };
        }
        public CellType CellType => Interpolation.CellType;
        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

        public int ID => throw new NotImplementedException(
            "Element type codes should be in a settings class. Even then it's a bad design choice");

        public IGaussPointExtrapolation3D GaussPointExtrapolation { get; }
        public IIsoparametricInterpolation3D Interpolation { get; }
        public IReadOnlyList<Node> Nodes { get; }
        public IQuadrature3D QuadratureForConsistentMass { get; }
        public IQuadrature3D QuadratureForStiffness { get; }

        public bool MaterialModified => throw new NotImplementedException();

        public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

        public IMatrix MassMatrix(IElement element)
        {
            return BuildCapacityMatrix();
        }

        public Matrix BuildCapacityMatrix()
        {
            int numDofs = Nodes.Count;
            var capacity = Matrix.CreateZero(numDofs, numDofs);
            IReadOnlyList<double[]> shapeFunctions =
                Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
            {
                Vector shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
                Matrix partial = shapeFunctionMatrix.TensorProduct(shapeFunctionMatrix);
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
                capacity.AxpyIntoThis(partial, dA);
            }

            //WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            return capacity;
        }
        public Matrix BuildMassTransportConductivityMatrix()
        {
            int numDofs = Nodes.Count;
            var convection = Matrix.CreateZero(numDofs, numDofs);
            IReadOnlyList<double[]> shapeFunctions =
                Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
            {
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                Vector shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
                Matrix shapeGradientsCartesian =
                   jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
                Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);
                Vector deformationX = deformation.GetRow(0);
                Vector deformationY = deformation.GetRow(1);
                Vector deformationZ = deformation.GetRow(2);
                Matrix partialK = shapeFunctionMatrix.TensorProduct(deformationX) * material.ConvectionCoeff[0] +
                    shapeFunctionMatrix.TensorProduct(deformationY) * material.ConvectionCoeff[1] +
                    shapeFunctionMatrix.TensorProduct(deformationZ) * material.ConvectionCoeff[2];
                double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight;
                convection.AxpyIntoThis(partialK, dA);
            }

            //WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            //convection.Scale(1);
            return convection;
        }
        public Matrix BuildLoadFromUnknownConductivityMatrix()
        {
            int numDofs = Nodes.Count;
            var conductivity = Matrix.CreateZero(numDofs, numDofs);
            IReadOnlyList<double[]> shapeFunctions =
                Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
            {
                Vector shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
                Matrix partial = shapeFunctionMatrix.TensorProduct(shapeFunctionMatrix);
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
                conductivity.AxpyIntoThis(partial, dA * material.LoadFromUnknownCoeff);
            }

            //WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            //conductivity.Scale(1);
            return conductivity;
        }
        public Matrix BuildStabilizingConductivityMatrix()
        {
            int numDofs = Nodes.Count;
            var convection = Matrix.CreateZero(numDofs, numDofs);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

            for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
            {
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                Matrix shapeGradientsCartesian =
                   jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
                Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);
                Vector deformationX = deformation.GetRow(0);
                Vector deformationY = deformation.GetRow(1);
                Vector deformationZ = deformation.GetRow(2);
                Matrix partialK = deformationX.TensorProduct(deformationY) * material.ConvectionCoeff[0] * material.ConvectionCoeff[1] +
                    deformationX.TensorProduct(deformationZ) * material.ConvectionCoeff[0] * material.ConvectionCoeff[2] +
                    deformationY.TensorProduct(deformationZ) * material.ConvectionCoeff[1] * material.ConvectionCoeff[2];
                double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight;
                convection.AxpyIntoThis(partialK, -dA);
            }

            //WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            //convection.Scale(1);
            return convection;
        }

        public Matrix BuildStabilizingLoadFromUnknownConductivityMatrix()
        {
            int numDofs = Nodes.Count;
            var convection = Matrix.CreateZero(numDofs, numDofs);
            IReadOnlyList<double[]> shapeFunctions =
                Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
            {
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                Vector shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
                Matrix shapeGradientsCartesian =
                   jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
                Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);
                Vector deformationX = deformation.GetRow(0);
                Vector deformationY = deformation.GetRow(1);
                Vector deformationZ = deformation.GetRow(2);
                Matrix partialK = shapeFunctionMatrix.TensorProduct(deformationX) * material.ConvectionCoeff[0] * material.LoadFromUnknownCoeff +
                    shapeFunctionMatrix.TensorProduct(deformationY) * material.ConvectionCoeff[1] * material.LoadFromUnknownCoeff +
                    shapeFunctionMatrix.TensorProduct(deformationZ) * material.ConvectionCoeff[2] * material.LoadFromUnknownCoeff;
                double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight;
                convection.AxpyIntoThis(partialK, -0.5 * dA );
            }

            //WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            //convection.Scale(1);
            return convection;
        }

        public IMatrix StiffnessMatrix(IElement element)
        {
            return BuildDiffusionConductivityMatrix() + BuildMassTransportConductivityMatrix() +
                BuildLoadFromUnknownConductivityMatrix();
        }
        public IMatrix DiffusionConductivityMatrix(IElement element)
        {
            return BuildDiffusionConductivityMatrix() + BuildLoadFromUnknownConductivityMatrix();
        }
        public IMatrix MassTransportConductivityMatrix(IElement element)
        {
            return BuildMassTransportConductivityMatrix();
        }

        public Matrix BuildDiffusionConductivityMatrix()
        {
            int numDofs = Nodes.Count;
            var conductivity = Matrix.CreateZero(numDofs, numDofs);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

            for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
            {
                // Calculate the necessary quantities for the integration
                //Matrix constitutive = (Matrix)(materialsAtGaussPoints[gp].ConstitutiveMatrix); // ugly cast will be removed along with the legacy Matrix classes
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                Matrix shapeGradientsCartesian =
                    jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
                Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);

                // Contribution of this gauss point to the element stiffness matrix
                Matrix partialK = deformation.Transpose() * deformation;
                //Matrix partialΚ = deformation.Transpose() * (constitutive * deformation);
                //partialK.Scale(materialsAtGaussPoints[gaussPoint].ThermalConductivity);

                double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight; //TODO: this is used by all methods that integrate. I should cache it.
                conductivity.AxpyIntoThis(partialK, dA * material.DiffusionCoeff);
                //conductivity.AxpyIntoThis(partialK, dA * 1);
            }
            
            //conductivity.Scale(1);
            return conductivity;
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

        public IMatrix DampingMatrix(IElement element)
        {
            return BuildStabilizingConductivityMatrix() + BuildStabilizingLoadFromUnknownConductivityMatrix();
        }
    }
}
