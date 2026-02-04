using System;

namespace MPOA
{
    class RA
    {
        public int n;
        public double[] center;
        public double radius;

        public double[] fungi;
        public double f;

        public int zapret;

        public bool includePoint(double[] p)
        {
            double sum = 0;
            for (int i = 0; i < n; i++)
            {
                sum += Math.Pow(p[i] - center[i], 2);
            }
            return (Math.Sqrt(sum) <= radius);
        }

        public bool Intersect(double[] p, double p_radius)
        {
            double sum = 0;
            for (int i = 0; i < n; i++)
            {
                sum += Math.Pow(p[i] - center[i], 2);
            }
            return (Math.Sqrt(sum) <= radius + p_radius);
        }

        public void decreaseZapret()
        {
            --zapret;
        }
    }
}
