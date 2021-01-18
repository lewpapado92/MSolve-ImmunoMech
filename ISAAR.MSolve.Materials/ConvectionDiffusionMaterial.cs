using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials
{
    public class ConvectionDiffusionMaterial
    {
        public ConvectionDiffusionMaterial(double diffusionCoeff, double[] convectionCoeff, double loadFromUnknownCoeff)
        {
            this.DiffusionCoeff = diffusionCoeff;
            this.ConvectionCoeff = convectionCoeff;
            this.LoadFromUnknownCoeff = loadFromUnknownCoeff;
        }

        public double DiffusionCoeff { get; }
        public double[] ConvectionCoeff { get; }
        public double LoadFromUnknownCoeff { get; }

        public ConvectionDiffusionMaterial Clone() => new ConvectionDiffusionMaterial(DiffusionCoeff, ConvectionCoeff, LoadFromUnknownCoeff);
    }
}
