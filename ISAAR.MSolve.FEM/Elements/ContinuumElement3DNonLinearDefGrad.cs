using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Mesh;
using System.Linq;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.Jacobians;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;
using ISSAR.MSolve.Discretization.Loads;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.Materials;

namespace ISAAR.MSolve.FEM.Elements
{
    /// <summary>
    /// Continuum finite Element for 3d problems with material and geometric nonlinearities
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class ContinuumElement3DNonLinearDefGrad : IStructuralFiniteElement//, IEmbeddedHostElement
    {
        protected readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
        protected readonly IDofType[][] dofTypes;
        private readonly IDynamicMaterial dynamicProperties;
        protected readonly IContinuumMaterial3DDefGrad[] materialsAtGaussPoints;
        protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

        private readonly int nGaussPoints;
        private bool isInitialized = false;

        private double[][] initialCoordinates; //not defined by user. 8 arrays of 3 elements
        private double[][] totalDisplacements;
        private double[] integrationCoeffs;

        private double[][] strainsVec;
        //private double[][] strainsVec_last_converged;
        private double[][] DefGradVec;
        private double lambdag = 1;

        protected ContinuumElement3DNonLinearDefGrad()
        {
        }

        public ContinuumElement3DNonLinearDefGrad(IReadOnlyList<Node> nodes, IContinuumMaterial3DDefGrad material, IQuadrature3D quadratureForStiffness,
             IIsoparametricInterpolation3D interpolation)
        {
            this.nGaussPoints = quadratureForStiffness.IntegrationPoints.Count;
            this.QuadratureForStiffness = quadratureForStiffness;
            this.Interpolation = interpolation;


            materialsAtGaussPoints = new IContinuumMaterial3DDefGrad[nGaussPoints];
            for (int i = 0; i < nGaussPoints; i++)
                materialsAtGaussPoints[i] = (IContinuumMaterial3DDefGrad)material.Clone();

            dofTypes = new IDofType[nodes.Count][];
            for (int i = 0; i < nodes.Count; i++)
            {
                dofTypes[i] = new IDofType[]
                {
                    StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
                };
            }
        }

        public ContinuumElement3DNonLinearDefGrad(IReadOnlyList<Node> nodes, IIsoparametricInterpolation3D interpolation,
            IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
            IGaussPointExtrapolation3D gaussPointExtrapolation,
            IContinuumMaterial3DDefGrad material, IDynamicMaterial dynamicProperties)
        {
            this.dynamicProperties = dynamicProperties;
            //this.materialsAtGaussPoints[] = materialsAtGaussPoints;
            //this.GaussPointExtrapolation = gaussPointExtrapolation;
            this.Nodes = nodes;
            this.Interpolation = interpolation;
            this.nGaussPoints = quadratureForStiffness.IntegrationPoints.Count;
            this.QuadratureForConsistentMass = quadratureForMass;
            this.QuadratureForStiffness = quadratureForStiffness;

            materialsAtGaussPoints = new IContinuumMaterial3DDefGrad[nGaussPoints];
            for (int i = 0; i < nGaussPoints; i++)
                materialsAtGaussPoints[i] = (IContinuumMaterial3DDefGrad)material.Clone();

            dofTypes = new IDofType[nodes.Count][];
            for (int i = 0; i < nodes.Count; i++)
            {
                dofTypes[i] = new IDofType[]
                {
                    StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
                };
            }

            //strainsVec = new double[nGaussPoints][];
            //strainsVecLastConverged = new double[nGaussPoints][];
            //for (int gpoint = 0; gpoint < materialsAtGaussPoints.Count; gpoint++)
            //{
            //    strainsVec[gpoint] = new double[6];
            //    strainsVecLastConverged[gpoint] = new double[6];
            //}
        }
        public ContinuumElement3DNonLinearDefGrad(IReadOnlyList<Node> nodes, IIsoparametricInterpolation3D interpolation,
            IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
            IGaussPointExtrapolation3D gaussPointExtrapolation,
            IContinuumMaterial3DDefGrad material, IDynamicMaterial dynamicProperties, double lambdag)
        {
            this.dynamicProperties = dynamicProperties;
            //this.materialsAtGaussPoints[] = materialsAtGaussPoints;
            //this.GaussPointExtrapolation = gaussPointExtrapolation;
            this.Nodes = nodes;
            this.Interpolation = interpolation;
            this.nGaussPoints = quadratureForStiffness.IntegrationPoints.Count;
            this.QuadratureForConsistentMass = quadratureForMass;
            this.QuadratureForStiffness = quadratureForStiffness;
            this.lambdag = lambdag;

            materialsAtGaussPoints = new IContinuumMaterial3DDefGrad[nGaussPoints];
            for (int i = 0; i < nGaussPoints; i++)
                materialsAtGaussPoints[i] = (IContinuumMaterial3DDefGrad)material.Clone();

            dofTypes = new IDofType[nodes.Count][];
            for (int i = 0; i < nodes.Count; i++)
            {
                dofTypes[i] = new IDofType[]
                {
                    StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
                };
            }

            //strainsVec = new double[nGaussPoints][];
            //strainsVecLastConverged = new double[nGaussPoints][];
            //for (int gpoint = 0; gpoint < materialsAtGaussPoints.Count; gpoint++)
            //{
            //    strainsVec[gpoint] = new double[6];
            //    strainsVecLastConverged[gpoint] = new double[6];
            //}
        }

        public IIsoparametricInterpolation3D Interpolation { get; }
        public IQuadrature3D QuadratureForStiffness { get; }
        public IQuadrature3D QuadratureForConsistentMass { get; }
        public IReadOnlyList<Node> Nodes { get; }

        public int ID => 13;
        public CellType CellType => Interpolation.CellType;

        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

        public IElementDofEnumerator DofEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public IReadOnlyList<IFiniteElementMaterial> Materials => materialsAtGaussPoints;

        public bool MaterialModified
        {
            get
            {
                foreach (IContinuumMaterial3DDefGrad material in materialsAtGaussPoints)
                    if (material.Modified) return true;
                return false;
            }
        }

        private Matrix[] Getbl13DeformationMatrices(IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives)
        {
            Matrix[] bl13Matrices;
            bl13Matrices = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                bl13Matrices[npoint] = Matrix.CreateZero(9, 3 * shapeFunctionNaturalDerivatives[npoint].NumRows);
                for (int m = 0; m < shapeFunctionNaturalDerivatives[npoint].NumRows; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        bl13Matrices[npoint][n, 3 * m + 0] = shapeFunctionNaturalDerivatives[npoint][m, n];
                        bl13Matrices[npoint][n + 3, 3 * m + 1] = shapeFunctionNaturalDerivatives[npoint][m, n];
                        bl13Matrices[npoint][n + 6, 3 * m + 2] = shapeFunctionNaturalDerivatives[npoint][m, n];
                    }
                }
            }
            return bl13Matrices;
        }

        private Matrix[] Getbl11aDeformationMatrices(Matrix[] jacobianInverse)
        {
            Matrix[] bl11aMatrices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                bl11aMatrices[gpoint] = Matrix.CreateZero(6, 9);
                for (int m = 0; m < 3; m++) // calculate first three data lines of the matrix
                {
                    for (int n = 0; n < 3; n++)
                    {
                        bl11aMatrices[gpoint][m, 3 * m + n] = jacobianInverse[gpoint][m, n];
                    }
                }
                for (int n = 0; n < 3; n++)
                {
                    bl11aMatrices[gpoint][3, n] = jacobianInverse[gpoint][1, n]; // calculate 4th data line
                    bl11aMatrices[gpoint][3, 3 + n] = jacobianInverse[gpoint][0, n];
                    bl11aMatrices[gpoint][4, 3 + n] = jacobianInverse[gpoint][2, n]; // calculate 5th data line
                    bl11aMatrices[gpoint][4, 6 + n] = jacobianInverse[gpoint][1, n];
                    bl11aMatrices[gpoint][5, 0 + n] = jacobianInverse[gpoint][2, n]; // calculate 6th data line
                    bl11aMatrices[gpoint][5, 6 + n] = jacobianInverse[gpoint][0, n];
                }
            }

            return bl11aMatrices;
        }

        private Matrix[] GetBL12DeformationMatrices(Matrix[] jacobianInverse)
        {
            Matrix[] bl12Marices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                bl12Marices[gpoint] = Matrix.CreateZero(9, 9);
                for (int m = 0; m < 3; m++) // calculate first three data lines of the matrix
                {
                    for (int n = 0; n < 3; n++)
                    {
                        bl12Marices[gpoint][m, 3 * m + n] = jacobianInverse[gpoint][0, n];
                    }
                }
                for (int m = 0; m < 3; m++) // calculate  data lines 4:6
                {
                    for (int n = 0; n < 3; n++)
                    {
                        bl12Marices[gpoint][3 + m, 3 * m + n] = jacobianInverse[gpoint][1, n];
                    }
                }
                for (int m = 0; m < 3; m++) // calculate  data lines 7:8
                {
                    for (int n = 0; n < 3; n++)
                    {
                        bl12Marices[gpoint][6 + m, 3 * m + n] = jacobianInverse[gpoint][2, n];
                    }
                }

            }

            return bl12Marices;
        }

        private Matrix[] Getbl01MDeformationMatrices(Matrix[] jacobianInverse)
        {
            Matrix[] bl01Matrices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                bl01Matrices[gpoint] = Matrix.CreateZero(6, 9);
                for (int m = 0; m < 3; m++) // calculate first three data lines of the matrix
                {
                    for (int n = 0; n < 3; n++)
                    {
                        bl01Matrices[gpoint][m, 3 * m + n] = jacobianInverse[gpoint][m, n];
                    }
                }
                for (int n = 0; n < 3; n++)
                {
                    bl01Matrices[gpoint][3, n] = jacobianInverse[gpoint][1, n]; // calculate 4th data line
                    bl01Matrices[gpoint][3, 3 + n] = jacobianInverse[gpoint][0, n];
                    bl01Matrices[gpoint][4, 3 + n] = jacobianInverse[gpoint][2, n]; // calculate 5th data line
                    bl01Matrices[gpoint][4, 6 + n] = jacobianInverse[gpoint][1, n];
                    bl01Matrices[gpoint][5, 0 + n] = jacobianInverse[gpoint][2, n]; // calculate 6th data line
                    bl01Matrices[gpoint][5, 6 + n] = jacobianInverse[gpoint][0, n];
                }
            }
            return bl01Matrices;
        }

        private Matrix[] GetAuxilliaryDeformationbnl1Matrices(Matrix[] jacobianInverse)
        {
            Matrix[] bnl1Matrices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                bnl1Matrices[gpoint] = Matrix.CreateZero(9, 9);
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            bnl1Matrices[gpoint][3 * m + n, 3 * m + p] = jacobianInverse[gpoint][n, p];
                        }
                    }
                }
            }
            return bnl1Matrices;
        }

        private void CalculateInitialConfigurationData(IElement element)
        {
            int numNodes = element.Nodes.Count;
            IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            Matrix[] bl13Matrices;
            bl13Matrices = Getbl13DeformationMatrices(shapeFunctionNaturalDerivatives);


            Matrix[] bnl1Matrices;

            initialCoordinates = new double[numNodes][];
            totalDisplacements = new double[numNodes][];

            var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(element.Nodes, x));

            Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            double[] jacobianDeterminants = jacobians.Select(x => x.DirectDeterminant).ToArray();





            integrationCoeffs = new double[nGaussPoints];

            bnl1Matrices = GetAuxilliaryDeformationbnl1Matrices(jacobianInverse);

            for (int j = 0; j < numNodes; j++)
            {
                initialCoordinates[j] = new double[] { element.Nodes[j].X, element.Nodes[j].Y, element.Nodes[j].Z, };
            }

            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                integrationCoeffs[gpoint] = jacobianDeterminants[gpoint] * QuadratureForStiffness.IntegrationPoints[gpoint].Weight;

            }

            totalDisplacements = new double[numNodes][];
            //strainsVec = new double[nGaussPoints][]; //MS
            //strainsVec_last_converged = new double[nGaussPoints][];
            DefGradVec = new double[nGaussPoints][];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                //strainsVec[gpoint] = new double[6]; //MS
                //strainsVec_last_converged[gpoint] = new double[6];
                DefGradVec[gpoint] = new double[9];
            }
            for (int k = 0; k < numNodes; k++)
            {
                totalDisplacements[k] = new double[3];
            }
            isInitialized = true;

        }

        private void UpdateCoordinateData(double[] localdisplacements, out double[][] deformedCoordinates)
        {
            int numNodes = localdisplacements.Length / 3;
            deformedCoordinates = new double[numNodes][];
            for (int j = 0; j < numNodes; j++)
            {
                deformedCoordinates[j] = new double[3];
                for (int k = 0; k < 3; k++)
                {
                    totalDisplacements[j][k] = localdisplacements[3 * j + k];
                    deformedCoordinates[j][k] = initialCoordinates[j][k] + totalDisplacements[j][k];
                }
            }
        }

        private void CalculateStrains(double[] localdisplacements, IElement element, double[][] deformedCoordinates)
        {
            IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(element.Nodes, x));
            Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            double[] jacobianDeterminants = jacobians.Select(x => x.DirectDeterminant).ToArray();
            //TODO: possibility of caching shapeFunctionNaturalDerivatives or J_0inv

            Matrix[] deformationGradientsTransposed = new Matrix[nGaussPoints];
            //Matrix[] GL = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                deformationGradientsTransposed[npoint] = Matrix.CreateZero(3, 3);
                //GL[npoint] = Matrix.CreateZero(3, 3);
            }

            var jacobiansDeformed = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(deformedCoordinates, x, false)).ToArray();
            Matrix[] jacobiansDeformedMatrices = jacobiansDeformed.Select(x => x.DirectMatrix).ToArray();


            //
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                //
                deformationGradientsTransposed[npoint] = jacobianInverse[npoint] * jacobiansDeformedMatrices[npoint];
                DefGradVec[npoint] = new double[9] { deformationGradientsTransposed[npoint][0, 0], deformationGradientsTransposed[npoint][1, 1],
                    deformationGradientsTransposed[npoint][2, 2], deformationGradientsTransposed[npoint][1, 0], deformationGradientsTransposed[npoint][2, 1],
                    deformationGradientsTransposed[npoint][0, 2], deformationGradientsTransposed[npoint][2, 0], deformationGradientsTransposed[npoint][0, 1],
                    deformationGradientsTransposed[npoint][1, 2], };//MS
                ////
                //GL[npoint] = deformationGradientsTransposed[npoint] * deformationGradientsTransposed[npoint].Transpose();
                //for (int m = 0; m < 3; m++)
            }

        }

        private double[] UpdateForces(IElement element)
        {
            //TODO: the gauss point loop should be the outer one


            // Matrices that are not currently cached are calculated here.
            int numNodes = element.Nodes.Count();
            Matrix ll2 = Matrix.CreateZero(numNodes, 3);
            for (int m = 0; m < numNodes; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    ll2[m, n] = totalDisplacements[m][n];
                }
            }
            IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(element.Nodes, x));
            Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            double[] jacobianDeterminants = jacobians.Select(x => x.DirectDeterminant).ToArray();

            Matrix[] bl13Matrices;
            bl13Matrices = Getbl13DeformationMatrices(shapeFunctionNaturalDerivatives);
            Matrix[] bl11aMatrices; // dimension number of gpoints
            Matrix[] bl12Marices;
            Matrix[] bl01Matrices;
            bl11aMatrices = Getbl11aDeformationMatrices(jacobianInverse);
            bl12Marices = GetBL12DeformationMatrices(jacobianInverse);
            bl01Matrices = Getbl01MDeformationMatrices(jacobianInverse);

            //INITIALIZATION of MAtrixes that are currently not cached
            double[][] integrCoeffsTimesStresses = new double[nGaussPoints][];
            Matrix[] blMatrices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                integrCoeffsTimesStresses[gpoint] = new double[6];
                blMatrices[gpoint] = Matrix.CreateZero(6, 3 * numNodes);
            }

            double[][] forces = new double[nGaussPoints + 1][];
            for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
            {
                forces[npoint] = new double[3 * numNodes];
            }

            Matrix[] bl11Matrices = new Matrix[nGaussPoints];
            Matrix[] bL1112Plus01Matrices = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                bl11Matrices[npoint] = Matrix.CreateZero(6, 9);
                bL1112Plus01Matrices[npoint] = Matrix.CreateZero(6, 9);
            }

            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                var stressesElastic = materialsAtGaussPoints[npoint].Stresses;
                var secondPiolaElastic = new double[3, 3] { {stressesElastic[0],stressesElastic[3],stressesElastic[5] },
                    { stressesElastic[3],stressesElastic[1],stressesElastic[4] },
                    {stressesElastic[5],stressesElastic[4],stressesElastic[2] } };
                //DefGradVec[npoint] = new double[9] { deformationGradientsTransposed[npoint][0, 0], 
                //deformationGradientsTransposed[npoint][1, 1], deformationGradientsTransposed[npoint][2, 2], 
                //    deformationGradientsTransposed[npoint][1, 0], deformationGradientsTransposed[npoint][2, 1], 
                //    deformationGradientsTransposed[npoint][0, 2], deformationGradientsTransposed[npoint][2, 0], 
                //    deformationGradientsTransposed[npoint][0, 1], deformationGradientsTransposed[npoint][1, 2], };//MS
                var defGradTransposed = new double[3, 3] {{DefGradVec[npoint][0],
                DefGradVec[npoint][7], DefGradVec[npoint][5] },
                    { DefGradVec[npoint][3], DefGradVec[npoint][1],
                    DefGradVec[npoint][8] }, {DefGradVec[npoint][6],
                    DefGradVec[npoint][4], DefGradVec[npoint][2] }, };
                var defGradElasticTransposed = new double[3, 3] {{DefGradVec[npoint][0] / lambdag,
                DefGradVec[npoint][7] / lambdag, DefGradVec[npoint][5] / lambdag },
                    { DefGradVec[npoint][3] / lambdag, DefGradVec[npoint][1] / lambdag,
                    DefGradVec[npoint][8] / lambdag}, {DefGradVec[npoint][6] / lambdag,
                    DefGradVec[npoint][4] / lambdag, DefGradVec[npoint][2] / lambdag }, };
                var firstPiolaElastic = new double[3, 3];
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            firstPiolaElastic[i, j] += secondPiolaElastic[i, k] * defGradElasticTransposed[k, j];
                        }
                    }
                }
                double defGradDeterminant = 0;
                for (int i = 0; i < 3; i++)
                    defGradDeterminant = defGradDeterminant + (defGradTransposed[0, i] * (defGradTransposed[1, (i + 1) % 3]
                        * defGradTransposed[2, (i + 2) % 3] - defGradTransposed[1, (i + 2) % 3] * defGradTransposed[2, (i + 1) % 3]));

                double defGradElasticDeterminant = 0;
                for (int i = 0; i < 3; i++)
                    defGradElasticDeterminant = defGradElasticDeterminant + (defGradElasticTransposed[0, i] * (defGradElasticTransposed[1, (i + 1) % 3]
                        * defGradElasticTransposed[2, (i + 2) % 3] - defGradElasticTransposed[1, (i + 2) % 3] * defGradElasticTransposed[2, (i + 1) % 3]));

                var firstPiola = new double[3, 3];
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        firstPiola[i, j] = defGradDeterminant / defGradElasticDeterminant * firstPiolaElastic[i, j] / lambdag / defGradElasticDeterminant;
                    }
                }

                double[,] defGradInverseTransposed = new double[3, 3] { { (defGradTransposed[2, 2]*defGradTransposed[1,1]- defGradTransposed[2, 1]
                    * defGradTransposed[1, 2])/defGradDeterminant,
                        (-(defGradTransposed[2, 2] * defGradTransposed[0, 1] - defGradTransposed[2, 1]
                    * defGradTransposed[0, 2]))/defGradDeterminant,
                        (defGradTransposed[1,2] * defGradTransposed[0, 1] - defGradTransposed[1, 1] * defGradTransposed[0, 2])/defGradDeterminant },
                    { (-(defGradTransposed[2,2]*defGradTransposed[1,0]-defGradTransposed[2,0]*defGradTransposed[1,2]))/defGradDeterminant,
                        (defGradTransposed[2,2]*defGradTransposed[0,0]-defGradTransposed[2,0]*defGradTransposed[0,2])/ defGradDeterminant,
                        (-(defGradTransposed[1,2]*defGradTransposed[0,0]-defGradTransposed[1,0]*defGradTransposed[0,2]))/defGradDeterminant },
                    {(defGradTransposed[2,1]*defGradTransposed[1,0]-defGradTransposed[2,0]*defGradTransposed[1,1])/defGradDeterminant,
                        (-(defGradTransposed[2,1]*defGradTransposed[0,0]-defGradTransposed[2,0]*defGradTransposed[0,1]))/defGradDeterminant,
                        (defGradTransposed[1,1]*defGradTransposed[0,0]-defGradTransposed[1,0]*defGradTransposed[0,1])/ defGradDeterminant } };

                var secondPiola = new double[3, 3];
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            secondPiola[i, j] += firstPiola[i, k] * defGradInverseTransposed[k, j];
                        }
                    }
                }

                var stresses = new double[6] { secondPiola[0, 0], secondPiola[1, 1], secondPiola[2, 2], secondPiola[0, 1], secondPiola[1, 2], secondPiola[2, 0] };
                integrCoeffsTimesStresses[npoint] = stresses.Scale(integrationCoeffs[npoint]);

                //
                Matrix lcyrcumflex;//= Matrix.CreateZero(3, 3);
                lcyrcumflex = shapeFunctionNaturalDerivatives[npoint].Transpose() * ll2;

                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            bl11Matrices[npoint][m, n] += bl11aMatrices[npoint][m, p] * lcyrcumflex[p, n];
                            bl11Matrices[npoint][m, 3 + n] += bl11aMatrices[npoint][m, 3 + p] * lcyrcumflex[p, n];
                            bl11Matrices[npoint][m, 6 + n] += bl11aMatrices[npoint][m, 6 + p] * lcyrcumflex[p, n];
                        }
                    }
                }

                //
                bL1112Plus01Matrices[npoint] = bl11Matrices[npoint] * bl12Marices[npoint];
                bL1112Plus01Matrices[npoint].AddIntoThis(bl01Matrices[npoint]);

                // 
                blMatrices[npoint] = bL1112Plus01Matrices[npoint] * bl13Matrices[npoint];

                //              
                forces[npoint] = blMatrices[npoint].Multiply(integrCoeffsTimesStresses[npoint], true);
            }

            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                forces[nGaussPoints].AddIntoThis(forces[npoint]);
            }

            return forces[nGaussPoints];
        }

        private Matrix UpdateKmatrices(IElement element)
        {
            int numNodes = element.Nodes.Count();
            Matrix elementStiffnessMatrix = Matrix.CreateZero(3 * numNodes, 3 * numNodes);


            // initialization of matrices that are not cached currently
            double[][] integrCoeffsTimesSpkvec = new double[nGaussPoints][];
            Matrix[] blMatrices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                integrCoeffsTimesSpkvec[gpoint] = new double[6];
                blMatrices[gpoint] = Matrix.CreateZero(6, 3 * numNodes);

            }
            Matrix totalDisplacementsMatrixReordered = Matrix.CreateZero(numNodes, 3);
            for (int m = 0; m < numNodes; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    totalDisplacementsMatrixReordered[m, n] = totalDisplacements[m][n];
                }
            }
            IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
            shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(element.Nodes, x));
            Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
            double[] jacobianDeterminants = jacobians.Select(x => x.DirectDeterminant).ToArray();

            Matrix[] bl13Matrices;
            bl13Matrices = Getbl13DeformationMatrices(shapeFunctionNaturalDerivatives);
            Matrix[] bl11aMatrices; // dimension: gpoints
            Matrix[] bl12Marices;
            Matrix[] bl01Matrices;
            bl11aMatrices = Getbl11aDeformationMatrices(jacobianInverse);
            bl12Marices = GetBL12DeformationMatrices(jacobianInverse);
            bl01Matrices = Getbl01MDeformationMatrices(jacobianInverse);

            Matrix[] bl11Matrices = new Matrix[nGaussPoints];
            Matrix[] bL1112Plus01Matrices = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                bl11Matrices[npoint] = Matrix.CreateZero(6, 9);
                bL1112Plus01Matrices[npoint] = Matrix.CreateZero(6, 9); //TODO this may be unnescessary
            }



            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                var stressesElastic = materialsAtGaussPoints[npoint].Stresses;
                var secondPiolaElastic = new double[3, 3] { {stressesElastic[0],stressesElastic[3],stressesElastic[5] },
                    { stressesElastic[3],stressesElastic[1],stressesElastic[4] },
                    {stressesElastic[5],stressesElastic[4],stressesElastic[2] } };
                //DefGradVec[npoint] = new double[9] { deformationGradientsTransposed[npoint][0, 0], 
                //deformationGradientsTransposed[npoint][1, 1], deformationGradientsTransposed[npoint][2, 2], 
                //    deformationGradientsTransposed[npoint][1, 0], deformationGradientsTransposed[npoint][2, 1], 
                //    deformationGradientsTransposed[npoint][0, 2], deformationGradientsTransposed[npoint][2, 0], 
                //    deformationGradientsTransposed[npoint][0, 1], deformationGradientsTransposed[npoint][1, 2], };//MS
                var defGradTransposed = new double[3, 3] {{DefGradVec[npoint][0],
                DefGradVec[npoint][7], DefGradVec[npoint][5] },
                    { DefGradVec[npoint][3], DefGradVec[npoint][1],
                    DefGradVec[npoint][8] }, {DefGradVec[npoint][6],
                    DefGradVec[npoint][4], DefGradVec[npoint][2] }, };
                var defGradElasticTransposed = new double[3, 3] {{DefGradVec[npoint][0] / lambdag,
                DefGradVec[npoint][7] / lambdag, DefGradVec[npoint][5] / lambdag },
                    { DefGradVec[npoint][3] / lambdag, DefGradVec[npoint][1] / lambdag,
                    DefGradVec[npoint][8] / lambdag}, {DefGradVec[npoint][6] / lambdag,
                    DefGradVec[npoint][4] / lambdag, DefGradVec[npoint][2] / lambdag }, };
                var firstPiolaElastic = new double[3, 3];
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            firstPiolaElastic[i, j] += secondPiolaElastic[i, k] * defGradElasticTransposed[k, j];
                        }
                    }
                }
                double defGradDeterminant = 0;
                for (int i = 0; i < 3; i++)
                    defGradDeterminant = defGradDeterminant + (defGradTransposed[0, i] * (defGradTransposed[1, (i + 1) % 3]
                        * defGradTransposed[2, (i + 2) % 3] - defGradTransposed[1, (i + 2) % 3] * defGradTransposed[2, (i + 1) % 3]));

                double defGradElasticDeterminant = 0;
                for (int i = 0; i < 3; i++)
                    defGradElasticDeterminant = defGradElasticDeterminant + (defGradElasticTransposed[0, i] * (defGradElasticTransposed[1, (i + 1) % 3]
                        * defGradElasticTransposed[2, (i + 2) % 3] - defGradElasticTransposed[1, (i + 2) % 3] * defGradElasticTransposed[2, (i + 1) % 3]));

                var firstPiola = new double[3, 3];
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        firstPiola[i, j] = defGradDeterminant / defGradElasticDeterminant * firstPiolaElastic[i, j] / lambdag;
                    }
                }

                double[,] defGradInverseTransposed = new double[3, 3] { { (defGradTransposed[2, 2]*defGradTransposed[1,1]- defGradTransposed[2, 1]
                    * defGradTransposed[1, 2])/defGradDeterminant,
                        (-(defGradTransposed[2, 2] * defGradTransposed[0, 1] - defGradTransposed[2, 1]
                    * defGradTransposed[0, 2]))/defGradDeterminant,
                        (defGradTransposed[1,2] * defGradTransposed[0, 1] - defGradTransposed[1, 1] * defGradTransposed[0, 2])/defGradDeterminant },
                    { (-(defGradTransposed[2,2]*defGradTransposed[1,0]-defGradTransposed[2,0]*defGradTransposed[1,2]))/defGradDeterminant,
                        (defGradTransposed[2,2]*defGradTransposed[0,0]-defGradTransposed[2,0]*defGradTransposed[0,2])/ defGradDeterminant,
                        (-(defGradTransposed[1,2]*defGradTransposed[0,0]-defGradTransposed[1,0]*defGradTransposed[0,2]))/defGradDeterminant },
                    {(defGradTransposed[2,1]*defGradTransposed[1,0]-defGradTransposed[2,0]*defGradTransposed[1,1])/defGradDeterminant,
                        (-(defGradTransposed[2,1]*defGradTransposed[0,0]-defGradTransposed[2,0]*defGradTransposed[0,1]))/defGradDeterminant,
                        (defGradTransposed[1,1]*defGradTransposed[0,0]-defGradTransposed[1,0]*defGradTransposed[0,1])/ defGradDeterminant } };

                var secondPiola = new double[3, 3];
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            secondPiola[i, j] += firstPiola[i, k] * defGradInverseTransposed[k, j];
                        }
                    }
                }

                var stresses = new double[6] { secondPiola[0, 0], secondPiola[1, 1], secondPiola[2, 2], secondPiola[0, 1], secondPiola[1, 2], secondPiola[2, 0] };

                // 
                integrCoeffsTimesSpkvec[npoint] = stresses.Scale(integrationCoeffs[npoint]);

                //
                Matrix lcyrcumflex = Matrix.CreateZero(3, 3);
                lcyrcumflex = shapeFunctionNaturalDerivatives[npoint].Transpose() * totalDisplacementsMatrixReordered;

                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            bl11Matrices[npoint][m, n] += bl11aMatrices[npoint][m, p] * lcyrcumflex[p, n];
                            bl11Matrices[npoint][m, 3 + n] += bl11aMatrices[npoint][m, 3 + p] * lcyrcumflex[p, n];
                            bl11Matrices[npoint][m, 6 + n] += bl11aMatrices[npoint][m, 6 + p] * lcyrcumflex[p, n];
                        }
                    }
                }

                // 
                bL1112Plus01Matrices[npoint] = bl11Matrices[npoint] * bl12Marices[npoint];
                bL1112Plus01Matrices[npoint].AddIntoThis(bl01Matrices[npoint]);

                //
                blMatrices[npoint] = bL1112Plus01Matrices[npoint] * bl13Matrices[npoint];

            }
            // TODO: BL and above calculations can cached from calculate forces method

            Matrix[] bnl1Matrices;
            Matrix[] bnlMatrices;
            bnl1Matrices = GetAuxilliaryDeformationbnl1Matrices(jacobianInverse);
            bnlMatrices = new Matrix[nGaussPoints];
            for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
            {
                //bnlMatrices[gpoint] = Matrix.CreateZero(9, 3*numNodes); //todo this may be unnescessary

                bnlMatrices[gpoint] = bnl1Matrices[gpoint] * bl13Matrices[gpoint];

            }


            Matrix[] integrCoeffsTimesStresses = new Matrix[nGaussPoints];
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                integrCoeffsTimesStresses[npoint] = Matrix.CreateZero(3, 3);
            }

            Matrix[] klStiffnessMatrixContributions = new Matrix[nGaussPoints + 1];
            Matrix[] knlStiffnessMatrixContributions = new Matrix[nGaussPoints + 1];
            for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
            {
                klStiffnessMatrixContributions[npoint] = Matrix.CreateZero(3 * numNodes, 3 * numNodes);
                knlStiffnessMatrixContributions[npoint] = Matrix.CreateZero(3 * numNodes, 3 * numNodes);
            }



            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                Matrix integrCoeffsTimesStressesTimesbnlMatrices = Matrix.CreateZero(9, 3 * numNodes); //TODO
                Matrix integrCoeffsTimesConsMatrix = Matrix.CreateZero(6, 6); //TODO
                Matrix integrCoeffTimesConsMatrixTimesBLMatrices = Matrix.CreateZero(6, 3 * numNodes);//TODO

                //
                integrCoeffsTimesStresses[npoint][0, 0] = integrCoeffsTimesSpkvec[npoint][0];
                integrCoeffsTimesStresses[npoint][0, 1] = integrCoeffsTimesSpkvec[npoint][3];
                integrCoeffsTimesStresses[npoint][0, 2] = integrCoeffsTimesSpkvec[npoint][5];
                integrCoeffsTimesStresses[npoint][1, 0] = integrCoeffsTimesSpkvec[npoint][3];
                integrCoeffsTimesStresses[npoint][1, 1] = integrCoeffsTimesSpkvec[npoint][1];
                integrCoeffsTimesStresses[npoint][1, 2] = integrCoeffsTimesSpkvec[npoint][4];
                integrCoeffsTimesStresses[npoint][2, 0] = integrCoeffsTimesSpkvec[npoint][5];
                integrCoeffsTimesStresses[npoint][2, 1] = integrCoeffsTimesSpkvec[npoint][4];
                integrCoeffsTimesStresses[npoint][2, 2] = integrCoeffsTimesSpkvec[npoint][2];

                //
                IMatrixView consDisp = materialsAtGaussPoints[npoint].ConstitutiveMatrix;

                for (int m = 0; m < 6; m++)
                {
                    for (int n = 0; n < 6; n++)
                    {
                        integrCoeffsTimesConsMatrix[m, n] = integrationCoeffs[npoint] * consDisp[m, n];
                    }
                }

                //
                integrCoeffTimesConsMatrixTimesBLMatrices = integrCoeffsTimesConsMatrix * blMatrices[npoint];

                //
                klStiffnessMatrixContributions[npoint] = blMatrices[npoint].Transpose() * integrCoeffTimesConsMatrixTimesBLMatrices;

                //
                for (int m = 0; m < 3; m++) // 3x24 dimensions
                {
                    for (int n = 0; n < 3 * numNodes; n++)
                    {
                        for (int p = 0; p < 3; p++)
                        {
                            integrCoeffsTimesStressesTimesbnlMatrices[m, n] += integrCoeffsTimesStresses[npoint][m, p] * bnlMatrices[npoint][p, n];
                            integrCoeffsTimesStressesTimesbnlMatrices[3 + m, n] += integrCoeffsTimesStresses[npoint][m, p] * bnlMatrices[npoint][3 + p, n];
                            integrCoeffsTimesStressesTimesbnlMatrices[6 + m, n] += integrCoeffsTimesStresses[npoint][m, p] * bnlMatrices[npoint][6 + p, n];
                        }
                    }
                }

                //
                knlStiffnessMatrixContributions[npoint] = bnlMatrices[npoint].Transpose() * integrCoeffsTimesStressesTimesbnlMatrices;
            }

            // Add contributions of each gp on the total element stiffness matrix elementStiffnessMatrix            
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                for (int m = 0; m < 3 * numNodes; m++)
                {
                    for (int n = 0; n < 3 * numNodes; n++)
                    {
                        klStiffnessMatrixContributions[nGaussPoints][m, n] += klStiffnessMatrixContributions[npoint][m, n];
                        knlStiffnessMatrixContributions[nGaussPoints][m, n] += knlStiffnessMatrixContributions[npoint][m, n];
                    }
                }
            }
            for (int m = 0; m < 3 * numNodes; m++)
            {
                for (int n = 0; n < 3 * numNodes; n++)
                {
                    elementStiffnessMatrix[m, n] = klStiffnessMatrixContributions[nGaussPoints][m, n] + knlStiffnessMatrixContributions[nGaussPoints][m, n];
                }
            }

            return elementStiffnessMatrix;
        }

        public Tuple<double[], double[]> CalculateStresses(IElement element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            this.UpdateCoordinateData(localTotalDisplacements, out double[][] deformedCoordinates);
            this.CalculateStrains(localTotalDisplacements, element, deformedCoordinates);
            //double[] strainsVec_strain_minus_last_converged_value = new double[6];
            for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
            {
                //strainsVec_strain_minus_last_converged_value = new double[6] 
                //{
                //    strainsVec[npoint][0]- strainsVec_last_converged[npoint][0],
                //    strainsVec[npoint][1] - strainsVec_last_converged[npoint][1],
                //    strainsVec[npoint][2] - strainsVec_last_converged[npoint][2],
                //    strainsVec[npoint][3]- strainsVec_last_converged[npoint][3],
                //    strainsVec[npoint][4]- strainsVec_last_converged[npoint][4],
                //    strainsVec[npoint][5]- strainsVec_last_converged[npoint][5]
                //};
                //materialsAtGaussPoints[npoint].UpdateMaterial(strainsVec_strain_minus_last_converged_value); 
                // //To update with total strain simplY = materialsAtGaussPoints[npoint].UpdateMaterial(strainsVec[npoint]);
                double[] DefGradVecEl = new double[9];
                for (int i = 0; i < 9; i++)
                {
                    DefGradVecEl[i] = DefGradVec[npoint][i] / lambdag;
                }
                materialsAtGaussPoints[npoint].UpdateMaterial(DefGradVecEl); //MS
            };
            //return new Tuple<double[], double[]>(strainsVec_strain_minus_last_converged_value, 
            //    materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses);

            return new Tuple<double[], double[]>(DefGradVec[materialsAtGaussPoints.Length - 1], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses);

            //TODO return data with total strains data would be:
            //return new Tuple<double[], double[]>(strainsVec[materialsAtGaussPoints.Length - 1], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses);
            //TODO: why return only the strain- stress of the gausspoint that is last on the array, Where is it needed?
        }

        public double[] CalculateForces(IElement element, double[] localTotalDisplacements, double[] localdDisplacements)
            => this.UpdateForces(element);

        public double[] CalculateForcesForLogging(IElement element, double[] localDisplacements)
            => CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);

        public virtual IMatrix StiffnessMatrix(IElement element)
        {
            if (!isInitialized)
            {
                int numNodes = element.Nodes.Count();
                this.CalculateInitialConfigurationData(element);
                var localTotalDisplacements = new double[3 * numNodes];
                this.UpdateCoordinateData(localTotalDisplacements, out double[][] deformedCoordinates);
                this.CalculateStrains(localTotalDisplacements, element, deformedCoordinates);
            }
            Matrix elementStiffness = this.UpdateKmatrices(element);
            //It doesn't implement Iembedded to return dof.Enumerator.GetTransformedMatrix
            return elementStiffness;
        }

        public void ResetMaterialModified()
        {
            foreach (IContinuumMaterial3DDefGrad material in materialsAtGaussPoints) material.ResetModified();
        }

        public void ClearMaterialState()
        {
            //TODO: the next throws an exception. Investigate. Possible changes in Analyzers may be the cause.
            //foreach (IContinuumMaterial3DDefGrad m in materialsAtGaussPoints) m.ClearState();
        }

        public void SaveMaterialState()
        {
            //for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
            //{
            //    for (int i1 = 0; i1 < 6; i1++)
            //    { strainsVec_last_converged[npoint][i1] = strainsVec[npoint][i1]; }
            //}

            foreach (IContinuumMaterial3DDefGrad m in materialsAtGaussPoints) m.SaveState();
        }

        public void ClearMaterialStresses()
        {
            foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearStresses();
        }

        public virtual IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofTypes;
        public virtual IMatrix MassMatrix(IElement element)
        {
            return ((DynamicMaterial)dynamicProperties).UseConsistentMass ? BuildConsistentMassMatrix() : BuildLumpedMassMatrix();
        }

        public Matrix BuildConsistentMassMatrix()
        {
            int numberOfDofs = 3 * Nodes.Count;
            var mass = Matrix.CreateZero(numberOfDofs, numberOfDofs);
            IReadOnlyList<double[]> shapeFunctions =
                Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
            {
                Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
                Matrix partial = shapeFunctionMatrix.MultiplyRight(shapeFunctionMatrix, true, false);
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
                mass.AxpyIntoThis(partial, dA);
            }

            mass.ScaleIntoThis(dynamicProperties.Density);
            return mass;
        }

        public Matrix BuildLumpedMassMatrix()
        {
            int numberOfDofs = 3 * Nodes.Count;
            var lumpedMass = Matrix.CreateZero(numberOfDofs, numberOfDofs);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            double area = 0;
            for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
            {
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                area += jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
            }

            double nodalMass = area * dynamicProperties.Density / Nodes.Count;
            for (int i = 0; i < numberOfDofs; i++)
            {
                lumpedMass[i, i] = nodalMass;
            }

            return lumpedMass;
        }
        private Matrix BuildShapeFunctionMatrix(double[] shapeFunctions)
        {
            var shapeFunctionMatrix = Matrix.CreateZero(3, 3 * shapeFunctions.Length);
            for (int i = 0; i < shapeFunctions.Length; i++)
            {
                shapeFunctionMatrix[0, 3 * i] = shapeFunctions[i];
                shapeFunctionMatrix[1, (3 * i) + 1] = shapeFunctions[i];
                shapeFunctionMatrix[2, (3 * i) + 2] = shapeFunctions[i];
            }

            return shapeFunctionMatrix;
        }

        #region not implemented
        public double[] CalculateAccelerationForces(IElement element, IList<MassAccelerationLoad> loads)
        {
            throw new NotImplementedException();
        }


        public virtual IMatrix DampingMatrix(IElement element)
        {
            IMatrix damping = StiffnessMatrix(element);
            damping.ScaleIntoThis(dynamicProperties.RayleighCoeffStiffness);
            damping.AxpyIntoThis(MassMatrix(element), dynamicProperties.RayleighCoeffMass);
            //IMatrix damping = null;
            return damping;
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }
        #endregion




    }


}