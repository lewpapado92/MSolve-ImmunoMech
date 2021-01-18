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
using ISAAR.MSolve.Materials;
using ISSAR.MSolve.Discretization.Loads;

namespace ISAAR.MSolve.FEM.Elements
{
    public class ODEElement3D : IFiniteElement, ICell<Node>
    {
        private readonly static IDofType[] nodalDOFTypes = new IDofType[] { ThermalDof.Temperature };
        private readonly IDofType[][] dofTypes; //TODO: this should not be stored for each element. Instead store it once for each Quad4, Tri3, etc. Otherwise create it on the fly.
        private readonly ODEMaterial material;

        public ODEElement3D(IReadOnlyList<Node> nodes, IIsoparametricInterpolation3D interpolation,
        IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
        IGaussPointExtrapolation3D gaussPointExtrapolation, ODEMaterial material)
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
            return BuildFirstOrderMatrix();
        }
        public Matrix BuildFirstOrderMatrix()
        {
            int numDofs = Nodes.Count;
            var firstOrderMatrix = Matrix.CreateZero(numDofs, numDofs);
            IReadOnlyList<double[]> shapeFunctions =
                Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
            {
                Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
                Matrix partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
                firstOrderMatrix.AxpyIntoThis(partial, dA * material.FirstOrderCoeff);
            }

            //WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            //firstOrderMatrix.Scale(material.FirstOrderCoeff);
            return firstOrderMatrix;
        }

        public IMatrix StiffnessMatrix(IElement element)
        {
            return BuildZerothOrderMatrix();
        }
        
        public Matrix BuildZerothOrderMatrix()
        {
            int numDofs = Nodes.Count;
            var zerothOrderMatrix = Matrix.CreateZero(numDofs, numDofs);
            IReadOnlyList<double[]> shapeFunctions =
                Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
            {
                Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
                Matrix partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
                zerothOrderMatrix.AxpyIntoThis(partial, dA * material.ZerothOrderCoeff);
            }

            //WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            //zerothOrderMatrix.Scale(material.ZerothOrderCoeff);
            return zerothOrderMatrix;
        }
        
        public Matrix BuildShapeFunctionMatrix(double[] shapeFunctions) //TODO: reconsider this. As it is, it just returns the shape functions in a Matrix
        {
            //var array2D = new double[1, shapeFunctions.Length];
            //for (int i = 0; i < shapeFunctions.Length; ++i)
            //{
            //    array2D[0, i] = shapeFunctions[i];
            //}
            //return Matrix.CreateFromArray(array2D);
            return Matrix.CreateFromArray(shapeFunctions, 1, shapeFunctions.Length);
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
            throw new NotImplementedException();
        }


    }
}