using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials
{
    public class ODEMaterial
    {
        public ODEMaterial(double firstOrderCoeff, double zerothOrderCoeff)
        {
            this.FirstOrderCoeff = firstOrderCoeff;
            this.ZerothOrderCoeff = zerothOrderCoeff;
        }

        public double FirstOrderCoeff { get; }
        public double ZerothOrderCoeff { get; }

        public ODEMaterial Clone() => new ODEMaterial(FirstOrderCoeff, ZerothOrderCoeff);
    }
}
