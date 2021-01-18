using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISSAR.MSolve.Discretization.Loads;

namespace ISAAR.MSolve.FEM.Elements.BoundaryConditionElements
{
    /// <summary>
    /// Finite element for implementation of boundary conditions along an edge on a boundary.
    /// </summary>
    public class EdgeBoundary2D : IFiniteElement, IEmbeddedElement
    {
        private const int numNodes = 2;
        private const int numDofs = 2;
        private static readonly IDofType[][] dofTypes = {
            new IDofType[] { ThermalDof.Temperature }, new IDofType[] { ThermalDof.Temperature } }; 

        private readonly ThermalMaterial material;

        public EdgeBoundary2D(IReadOnlyList<Node> nodes, double crossSectionArea, ThermalMaterial material)
        {
            Debug.Assert(nodes.Count == 2, "Thermal rod element must have exactly 2 nodes.");
            this.material = material;
            this.Nodes = nodes;
            this.CrossSectionArea = crossSectionArea;
            this.Length = nodes[0].CalculateEuclidianDistanceFrom(nodes[1]);
        }

        public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

        public int ID => throw new NotImplementedException(
            "Element type codes should be in a settings class. Even then it's a bad design choice");
        public CellType CellType { get; } = CellType.Line;

        public double CrossSectionArea { get; }
        public double Length { get; }
        public IReadOnlyList<Node> Nodes { get; }

        public bool MaterialModified => throw new NotImplementedException();

        public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

        public IList<EmbeddedNode> EmbeddedNodes { get; } = new List<EmbeddedNode>();

        public IMatrix MassMatrix(IElement element)
        {
            return BuildCapacityMatrix();
        }

        public Matrix BuildCapacityMatrix()
        {
            double kdAL = material.SpecialHeatCoeff * material.Density * CrossSectionArea * Length;
            double[,] capacity = { { kdAL / 3.0, kdAL/ 6.0 }, { kdAL / 6.0, kdAL / 3.0 } };
            return Matrix.CreateFromArray(capacity);
        }

        public Matrix BuildDiffusionMatrix()
        {

            double cAoverL = material.ThermalConductivity * CrossSectionArea / Length;
            double[,] conductivity = { { cAoverL, -cAoverL }, { -cAoverL, cAoverL } };
            return Matrix.CreateFromArray(conductivity);
        }
        public Matrix BuildBoundaryConductivityMatrix()
        {
            //Vector normalVector = 
            double cAoverL = material.ThermalConductivity * CrossSectionArea / 2;
            double[,] convection = { { -cAoverL, -cAoverL }, { cAoverL, cAoverL } };
            return Matrix.CreateFromArray(convection);
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
            return DofEnumerator.GetTransformedMatrix(BuildDiffusionMatrix());
        }

        public IMatrix DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        public Dictionary<IDofType, int> GetInternalNodalDOFs(Element element, Node node)
        {
            if (node.ID == this.Nodes[0].ID) return new Dictionary<IDofType, int> { { ThermalDof.Temperature, 0 } };
            else if (node.ID == this.Nodes[1].ID) return new Dictionary<IDofType, int> { { ThermalDof.Temperature, 1 } };
            else throw new ArgumentException($"GetInternalNodalDOFs: Node {node.ID} not found in element {element.ID}.");
        }

        public double[] GetLocalDOFValues(Element hostElement, double[] hostDOFValues)
        {
            return DofEnumerator.GetTransformedDisplacementsVector(hostDOFValues);
        }
    }
}
