using System;

namespace MPOA
{
    class Program
    {
        static void Main(string[] args)
        {
            ObjFunc function = new Schwefel();
            double solution = 837.9658;

            double[][] FS = new double[2][];
            FS[0] = new double[2];
            FS[1] = new double[2];
            FS[0][0] = -500; FS[0][1] = 500;
            FS[1][0] = -500; FS[1][1] = 500;

            // гиперпараметры по умолчанию
            AlgOpt alg = new MPOA(2, function, FS);
            alg.Run();

            for (int i = 0; i < 2; i++)
                Console.WriteLine("x" + (i + 1).ToString() + ": " + alg.optimum[i].ToString());
            Console.WriteLine("f(x): " + (-alg.f_opt).ToString());
            Console.WriteLine("Отклонение по модулю: " + (Math.Abs(solution - -(alg.f_opt))).ToString());
        }
    }
}
