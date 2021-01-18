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
using ISAAR.MSolve.Materials;
using ISSAR.MSolve.Discretization.Loads;

namespace ISAAR.MSolve.FEM.Elements
{
    /// <summary>
    /// Finite element for convection-diffusion transfer along a single direction. Can be used for 1D, 2D and 3D problems. Does not take into 
    /// account geometric non-linearities.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ConvectionDiffusionRod : IConvectionDiffusionElement, IEmbeddedElement
    {
        private const int numNodes = 2;
        private const int numDofs = 2;
        private static readonly IDofType[][] dofTypes = {
            new IDofType[] { ThermalDof.Temperature }, new IDofType[] { ThermalDof.Temperature } };

        private readonly ConvectionDiffusionMaterial material;

        public ConvectionDiffusionRod(IReadOnlyList<Node> nodes, double crossSectionArea, ConvectionDiffusionMaterial material)
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
            double kdAL =  CrossSectionArea * Length;
            double[,] capacity = { { kdAL / 3.0, kdAL / 6.0 }, { kdAL / 6.0, kdAL / 3.0 } };
            return Matrix.CreateFromArray(capacity);
        }

        public Matrix BuildDiffusionConductivityMatrix()
        {

            double cAoverL = material.DiffusionCoeff * CrossSectionArea / Length;
            double[,] conductivity = { { cAoverL, -cAoverL }, { -cAoverL, cAoverL } };
            return Matrix.CreateFromArray(conductivity);
        }
        public Matrix BuildMassTransportConductivityMatrix()
        {

            double conA = material.ConvectionCoeff[0] * CrossSectionArea / 2;
            double[,] conductivity = { { -conA, conA }, { -conA, conA } };
            return Matrix.CreateFromArray(conductivity);
        }
        public Matrix BuildLoadFromUnknownConductivityMatrix()
        {

            double conA = material.LoadFromUnknownCoeff * CrossSectionArea * Length;
            double[,] conductivity = { { conA / 3.0, conA / 6.0 }, { conA / 6.0, conA / 3.0 } };
            return Matrix.CreateFromArray(conductivity);
        }
        public Matrix BuildStabilizingConductivityMatrix()
        {

            double cAoverL = -.5 * Math.Pow(material.ConvectionCoeff[0],2) * CrossSectionArea / Length;
            double[,] conductivity = { { cAoverL, -cAoverL }, { -cAoverL, cAoverL } };
            return Matrix.CreateFromArray(conductivity);
        }
        public Matrix BuildStabilizingLoadFromUnknownConductivityMatrix()
        {
            double conA = -.5 * material.LoadFromUnknownCoeff * material.ConvectionCoeff[0] * CrossSectionArea / 2;
            double[,] conductivity = { { -conA, conA }, { -conA, conA } };
            return Matrix.CreateFromArray(conductivity);
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
            return DofEnumerator.GetTransformedMatrix(BuildDiffusionConductivityMatrix() + BuildMassTransportConductivityMatrix() + 
                BuildLoadFromUnknownConductivityMatrix());
            //return BuildDiffusionConductivityMatrix() + BuildMassTransportConductivityMatrix();
        }

        public IMatrix DampingMatrix(IElement element)
        {
            return DofEnumerator.GetTransformedMatrix(BuildStabilizingConductivityMatrix() + 
                BuildStabilizingLoadFromUnknownConductivityMatrix());
            //return BuildStabilizingConductivityMatrix();
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

        public IMatrix MassTransportConductivityMatrix(IElement element)
        {
            return DofEnumerator.GetTransformedMatrix(BuildMassTransportConductivityMatrix());
        }

        public IMatrix DiffusionConductivityMatrix(IElement element)
        {
            return DofEnumerator.GetTransformedMatrix(BuildDiffusionConductivityMatrix());
        }
    }
}
