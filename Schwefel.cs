using System;

namespace MPOA
{
    class Schwefel : ObjFunc
    {
        public double Call(int n, double[] x)
        { // минимизировать
            double x1 = x[0];
            double x2 = x[1];
            return -(x1*Math.Sin(Math.Sqrt(Math.Abs(x1))) + x2*Math.Sin(Math.Sqrt(Math.Abs(x2))));
        }
    }
}
