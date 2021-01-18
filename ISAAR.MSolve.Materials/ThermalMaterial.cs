using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials
{
    public class ThermalMaterial
    {
        public ThermalMaterial(double density, double specialHeatCoeff, double thermalConductivity, double thermalConvection)
        {
            this.Density = density;
            this.SpecialHeatCoeff = specialHeatCoeff;
            this.ThermalConductivity = thermalConductivity;
            this.ThermalConvection = thermalConvection;
        }

        public double Density { get; }
        public double SpecialHeatCoeff { get; }
        public double ThermalConductivity { get; }
        public double ThermalConvection { get; }

        public ThermalMaterial Clone() => new ThermalMaterial(Density, SpecialHeatCoeff, ThermalConductivity, ThermalConvection);
    }
}
