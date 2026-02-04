using System;

namespace MPOA
{
    public abstract class AlgOpt
    {
        protected int n;
        public ObjFunc function;
        public double[] optimum;
        public double f_opt;
        protected double[][] FS; // Feasible Set
        protected AlgParameter[] parameters;
        public int count_of_parameters;

        public AlgParameter GetAlgParameter(int no)
        {
            if (no < 0 || no >= parameters.Length) throw new IndexOutOfRangeException();
            return parameters[no];
        }
        public void SetAlgParameter(int no, String name, Object value)
        {
            if (no < 0 || no >= parameters.Length) throw new IndexOutOfRangeException();
            parameters[no].Name = name;
            parameters[no].Value = value;
        }
        public abstract void Run();
        public abstract void Run(int time);
        public abstract void RunParallel();
        public abstract void RunParallel(int time);
    }
}
