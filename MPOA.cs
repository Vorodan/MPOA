using System;
using System.Collections.Generic;
using System.IO;

namespace MPOA
{
    class MPOA : AlgOpt
    {
        int T = 10; // количество проходов локального поиска
        int Np = 5; // количество грибников в группе
        int PC = 10; //количество точек, посещаемых грибником в процессе исследования грибного места
        int RAt = 6; // количество итераций, в течение которых грибникам будет запрещено посещать исследованные грибные места
        double r = 0.1; // радиус грибного места
        double alpha = 1; // величина шага перелёта Леви
        double lambda = 1.5; // параметр распределения Леви
        double ep = 1; // extrasearch possibility
        double D_eps = 1e-6; // нижняя грань показателя рассредоточения
        double rho_eps = 1e-3; // минимально возможный радиус самой большой гиперсферы локального поиска

        double s = 0.2; // масштаб спирали
        double eps1 = 1e-6; // для окончания локального поиска
        double eps2 = 1e-6;

        double[][] pickers;

        double f_max;

        (double[] fungi, double f)[] Basket;

        List<RA> RAs;
        List<RA> new_RAs;

        Random randomd, randomi;

        public MPOA(int n, ObjFunc function, double[][] FS, int T, int Np, int PC, int RAt, double r, double alpha, double lambda, double ep, double D_eps, double rho_eps, double s, double eps1, double eps2)
        {
            this.n = n;
            this.function = function;
            this.T = T;
            this.Np = Np;
            this.PC = PC;
            this.RAt = RAt;
            this.r = r;
            this.alpha = alpha;
            this.lambda = lambda;
            this.ep = ep;
            this.D_eps = D_eps;
            this.rho_eps = rho_eps;
            this.s = s;
            this.eps1 = eps1;
            this.eps2 = eps2;

            this.FS = new double[n][];
            for (int i = 0; i < n; i++)
            {
                this.FS[i] = new double[2];
                this.FS[i][0] = FS[i][0];
                this.FS[i][1] = FS[i][1];
            }

            count_of_parameters = 13;
            parameters = new AlgParameter[count_of_parameters];
            parameters[0] = new AlgParameter("T", T);
            parameters[1] = new AlgParameter("Np", Np);
            parameters[2] = new AlgParameter("PC", PC);
            parameters[3] = new AlgParameter("RAt", RAt);
            parameters[4] = new AlgParameter("r", r);
            parameters[5] = new AlgParameter("alpha", alpha);
            parameters[6] = new AlgParameter("lambda", lambda);
            parameters[7] = new AlgParameter("ep", ep);
            parameters[8] = new AlgParameter("D_eps", D_eps);
            parameters[9] = new AlgParameter("rho_eps", rho_eps);
            parameters[10] = new AlgParameter("s", s);
            parameters[11] = new AlgParameter("eps1", eps1);
            parameters[12] = new AlgParameter("eps2", eps2);

            pickers = new double[Np][];
            for (int j = 0; j < Np; j++)
                pickers[j] = new double[n];

            RAs = new List<RA>();
            new_RAs = new List<RA>();

            Basket = null;
        }

        public MPOA(int n, ObjFunc function, double[][] FS)
        {
            this.n = n;
            this.function = function;
            this.FS = new double[n][];
            for (int i = 0; i < n; i++)
            {
                this.FS[i] = new double[2];
                this.FS[i][0] = FS[i][0];
                this.FS[i][1] = FS[i][1];
            }

            count_of_parameters = 13;
            parameters = new AlgParameter[count_of_parameters];
            parameters[0] = new AlgParameter("T", T);
            parameters[1] = new AlgParameter("Np", Np);
            parameters[2] = new AlgParameter("PC", PC);
            parameters[3] = new AlgParameter("RAt", RAt);
            parameters[4] = new AlgParameter("r", r);
            parameters[5] = new AlgParameter("alpha", alpha);
            parameters[6] = new AlgParameter("lambda", lambda);
            parameters[7] = new AlgParameter("ep", ep);
            parameters[8] = new AlgParameter("D_eps", D_eps);
            parameters[9] = new AlgParameter("rho_eps", rho_eps);
            parameters[10] = new AlgParameter("s", s);
            parameters[11] = new AlgParameter("eps1", eps1);
            parameters[12] = new AlgParameter("eps2", eps2);

            pickers = new double[Np][];
            for (int j = 0; j < Np; j++)
                pickers[j] = new double[n];

            RAs = new List<RA>();
            new_RAs = new List<RA>();

            Basket = null;
        }

        public override void Run()
        {
            randomd = new Random();
            randomi = new Random();

            T = Int32.Parse(parameters[0].Value.ToString());
            Np = Int32.Parse(parameters[1].Value.ToString());
            PC = Int32.Parse(parameters[2].Value.ToString());
            RAt = Int32.Parse(parameters[3].Value.ToString());
            r = Double.Parse(parameters[4].Value.ToString());
            alpha = Double.Parse(parameters[5].Value.ToString());
            lambda = Double.Parse(parameters[6].Value.ToString());
            ep = Double.Parse(parameters[7].Value.ToString());
            D_eps = Double.Parse(parameters[8].Value.ToString());
            rho_eps = Double.Parse(parameters[9].Value.ToString());
            s = Double.Parse(parameters[10].Value.ToString());
            eps1 = Double.Parse(parameters[11].Value.ToString());
            eps2 = Double.Parse(parameters[12].Value.ToString());

            pickers = new double[Np][];
            for (int j = 0; j < Np; j++) pickers[j] = new double[n];

            GenPop();
            RAs.Clear();
            Basket = null;
            int k = 1;

            do
            { // глобальный поиск
                LevyFlight(k);

                SetNewRA(k);

                ClearOldRA();

                Search(k);

                SortBasket();

                // обновление RAs, добавление из newRAs
                UpdateRA();

                k++; // здесь, чтобы не было коллизии в traj3 при k=1
                if (randomd.NextDouble() <= ep)
                    ExtraSearch(k);

                if (GlobalEnd()) break;

            } while (true);

            RAs.Clear();
            SortBasket();

            // локальный поиск
            {
                double rho;
                double rho_k;
                double[] m_tau = new double[n];
                Array.Copy(Basket[0].fungi, m_tau, n);
                double f_tau = Basket[0].f;
                double[] m_tau_prev = new double[n]; double f_tau_prev = Double.MaxValue;

                for (int tau = 0; tau < T; tau++)
                {
                    if (tau > 0)
                    {
                        if (Math.Abs(f_tau - f_tau_prev) < eps1) break;
                        double cur;
                        double max = 0;
                        for (int i = 0; i < n; i++)
                        {
                            if (m_tau[i] == 0)
                            {
                                max = eps2 + 1;
                                break;
                            }
                            cur = Math.Abs((m_tau[i] - m_tau_prev[i]) / m_tau[i]);
                            if (max < cur) max = cur;
                        }
                        if (max < eps2) break;
                    }
                    rho = CalcRho(); // когда равно нулю, значит все грибы в одном месте, скорее всего
                                     // тогда локальный поиск холостой
                    rho_k = rho;
                    k = 1;
                    while (rho_k >= rho_eps)
                    {
                        SetSpheres(rho_k);
                        for (int j = 0; j < Np; j++) DoSphere(j);
                        SortBasket();
                        k++;
                        rho_k /= 2;
                        new_RAs.Clear();
                    }
                    Array.Copy(m_tau, m_tau_prev, n);
                    f_tau_prev = f_tau;
                    Array.Copy(Basket[0].fungi, m_tau, n);
                    f_tau = Basket[0].f;
                }
            }

            optimum = new double[n];
            Array.Copy(ToAbsolutSize(Basket[0].fungi), optimum, n);
            f_opt = Basket[0].f;
        }

        public override void Run(int time)
        {
            throw new NotImplementedException();
        }

        public override void RunParallel()
        {
            throw new NotImplementedException();
        }

        public override void RunParallel(int time)
        {
            throw new NotImplementedException();
        }
        
        void GenPop()
        {
            for (int j = 0; j < Np; j++)
            {
                int eingang_index = randomi.Next(n);
                if (randomi.Next(n) % 2 == 0) pickers[j][eingang_index] = 0;
                else pickers[j][eingang_index] = 1;
                for (int i = 0; i < n; i++)
                {
                    if (i == eingang_index) continue;
                    pickers[j][i] = randomd.NextDouble();
                }
            }
        }

        void LevyFlight(int k)
        {
            if (k == 1) return;
            for (int j = 0; j < Np; j++)
            {
                while (true)
                {
                    double[] new_pos = new double[n];
                    for (int i = 0; i < n; i++)
                    {
                        double Ri = randomd.NextDouble() * (1.0 - 1e-7) + 1e-7;
                        double Thetai = Ri*2*Math.PI;
                        double Li = Math.Pow(Ri, -1 / lambda);
                        if (i < n / 2) new_pos[i] = pickers[j][i] + alpha*Li*Math.Sin(Thetai);
                        else new_pos[i] = pickers[j][i] + alpha*Li*Math.Cos(Thetai);
                        if (new_pos[i] < 0 || new_pos[i] > 1)
                            new_pos[i] = randomd.NextDouble();
                    }
                    if (InRA(new_pos)) continue;
                    Array.Copy(new_pos, pickers[j], n);
                    break;
                }
            }
        }

        void SetNewRA(int k)
        {
            for (int j = 0; j < Np; j++)
            {
                RA new_ra = new RA();
                new_ra.n = n;
                new_ra.center = new double[n];
                Array.Copy(pickers[j], new_ra.center, n);
                new_ra.radius = r;
                new_ra.zapret = RAt;
                new_ra.fungi = new double[n];

                for (int i = 0; i < RAs.Count; i++)
                {
                    if (new_ra.Intersect(RAs[i].center, RAs[i].radius))
                    {
                        new_ra.radius /= 2;
                        i--;
                    }
                }

                if (k > 1)
                    for (int i = 0; i < n; i++)
                    {
                        if (new_ra.center[i] + new_ra.radius > 1 ||
                            new_ra.center[i] - new_ra.radius < 0)
                        {
                            new_ra.radius /= 2;
                            i--;
                        }
                    }

                new_RAs.Add(new_ra);
            }
        }

        void ClearOldRA()
        {
            for (int i = 0; i < RAs.Count; i++)
            {
                RAs[i].zapret -= 1;
                if (RAs[i].zapret == 0)
                {
                    RAs.RemoveAt(i);
                    i--;
                }
            }
        }

        void Search(int k)
        {
            for (int j = 0; j < Np; j++)
            {
                if (k == 1) { DoRandomWalk(j, k); continue; }
                switch (randomi.Next(5))
                {
                    case 0:
                        DoLine(j);
                        break;
                    case 1:
                        DoSphere(j);
                        break;
                    case 2:
                        DoRandomWalk(j, k);
                        break;
                    case 3:
                        DoSpiral(j);
                        break;
                    case 4:
                        DoLoop(j);
                        break;
                }
            }
        }

        bool InRA(double[] pos)
        { // заход в запретную зону
            foreach (RA ra in RAs)
                if (ra.includePoint(pos)) return true;
            return false;
        }

        void GenVec(double[] vec)
        {
            double sum = 0;
            for (int j = 0; j < n; j++)
            {
                vec[j] = randomd.NextDouble();
                if (randomd.NextDouble() > 0.5) vec[j] *= -1;
                sum += Math.Pow(vec[j], 2);
            }
            double length = Math.Sqrt(sum);
            if (length > 0)
            {
                for (int j = 0; j < n; j++)
                    vec[j] /= length;
            }
        }

        void DoLine(int no)
        {
            double radius = new_RAs[no].radius;
            double step = radius / (PC - 1);
            double[] vec = new double[n];
            double length;
            GenVec(vec);
            for (int i = 0; i < n; i++)
                pickers[no][i] = new_RAs[no].center[i] + vec[i] * radius; // первая точка на сфере
            double[] best_point = new double[n];
            Array.Copy(pickers[no], best_point, n);
            double best_f = function.Call(n, ToAbsolutSize(best_point));
            double worst_f = f_max;

            for (int q = 1; q < PC; q++)
            {
                for (int i = 0; i < n; i++)
                {
                    pickers[no][i] -= vec[i] * step;
                    if (pickers[no][i] < 0 || pickers[no][i] > 1)
                        break;
                    if (pickers[no][i] is double.NaN)
                        break;
                }


                double f = function.Call(n, ToAbsolutSize(pickers[no]));
                if (f < best_f)
                {
                    best_f = f;
                    Array.Copy(pickers[no], best_point, n);
                }
                if (f > f_max) f_max = f;
            }

            Array.Copy(best_point, new_RAs[no].fungi, n);
            new_RAs[no].f = best_f;

            for (int i = 0; i < n; i++)
            {
                if (best_point[i] < 0 || best_point[i] > 1)
                    break;

            }
        }

        void DoLineExtra(int alpha)
        { // считаем только точки между alpha и beta, новый гриб не будет в beta
            double radius = new_RAs[0].radius;
            double step = radius/(PC - 1);
            double[] vec = new double[n];
            double length = 0;

            for (int i = 0; i < n; i++)
            { // грибник в alpha
                vec[i] = pickers[0][i] - new_RAs[0].center[i];
                if (vec[i] == 0)
                    vec[i] = 0;
                length += Math.Pow(vec[i], 2);
            }
            if (length == 0)
                return; // это для ситуаций, когда грибы не отличимы
            length = Math.Sqrt(length);
            for (int i = 0; i < n; i++)
                vec[i] /= length;

            double[] best_point = new double[n];
            Array.Copy(pickers[0], best_point, n);
            double best_f = Basket[alpha].f;

            double worst_f = f_max;

            for (int q = 1; q < PC - 1; q++)
            {
                for (int i = 0; i < n; i++)
                {
                    pickers[0][i] -= vec[i] * step;
                    if (pickers[0][i] < 0 || pickers[0][i] > 1)
                        break;
                    if (pickers[0][i] is double.NaN)
                        break;
                }


                double f = function.Call(n, ToAbsolutSize(pickers[0]));
                if (f < best_f)
                {
                    best_f = f;
                    Array.Copy(pickers[0], best_point, n);
                }
                if (f > f_max) f_max = f;
            }

            Array.Copy(best_point, new_RAs[0].fungi, n);
            new_RAs[0].f = best_f;

            for (int i = 0; i < n; i++)
            {
                if (best_point[i] < 0 || best_point[i] > 1)
                    break;

            }
        }

        void DoSphere(int no)
        {
            double[] vec = new double[n];
            double radius = new_RAs[no].radius;

            for (int i = 0; i < n; i++)
                pickers[no][i] = new_RAs[no].center[i] + vec[i] * radius; // первая точка на сфере
            double[] best_point = new double[n];
            Array.Copy(pickers[no], best_point, n);
            double best_f = function.Call(n, ToAbsolutSize(best_point));

            double[] new_pos = new double[n];
            for (int q = 1; q < PC; q++)
            {
                bool _out = false;
                GenVec(vec);
                for (int i = 0; i < n; i++)
                {
                    new_pos[i] = new_RAs[no].center[i] + vec[i] * radius;
                    if (new_pos[i] < 0 || new_pos[i] > 1)
                    { // выход за границу области
                        q--;
                        _out = true;
                        break;
                    }
                }
                if (_out) continue;

                double f = function.Call(n, ToAbsolutSize(new_pos));
                if (f < best_f)
                {
                    best_f = f;
                    Array.Copy(new_pos, best_point, n);
                }
                if (f > f_max) f_max = f;
            }

            Array.Copy(best_point, new_RAs[no].fungi, n);
            Array.Copy(new_pos, pickers[no], n);
            new_RAs[no].f = best_f;

            for (int i = 0; i < n; i++)
            {
                if (best_point[i] < 0 || best_point[i] > 1)
                    break;
            }
        }

        void DoRandomWalk(int no, int k)
        {
            double[] vec = new double[n];
            double radius = new_RAs[no].radius;
            double step = radius / (PC - 1);

            double[] best_point = new double[n];
            Array.Copy(pickers[no], best_point, n);
            double best_f = function.Call(n, ToAbsolutSize(best_point)); // первая точка в центре
            if (k == 1 && no == 0) f_max = best_f; // инициализация f_max (только здесь)

            double[] new_pos = new double[n];
            for (int q = 1; q < PC; q++)
            {
                bool _out = false;
                GenVec(vec);
                for (int i = 0; i < n; i++)
                {
                    new_pos[i] = pickers[no][i] + vec[i] * step;
                    if (new_pos[i] < 0 || new_pos[i] > 1)
                    { // выход за границу области
                        q--;
                        _out = true;
                        break;
                    }
                }
                if (_out) continue;

                Array.Copy(new_pos, pickers[no], n);

                double f = function.Call(n, ToAbsolutSize(new_pos));
                if (f < best_f)
                {
                    best_f = f;
                    Array.Copy(new_pos, best_point, n);
                }
                if (f > f_max) f_max = f;
            }

            Array.Copy(best_point, new_RAs[no].fungi, n);
            new_RAs[no].f = best_f;

            for (int i = 0; i < n; i++)
            {
                if (best_point[i] < 0 || best_point[i] > 1)
                    break;
            }
        }

        void DoLoop(int no)
        {
            double[] vec1 = new double[n];
            GenVec(vec1);
            int z = -1;
            for (int i = 0; i < n; i++)
            {
                if (vec1[i] != 0)
                {
                    z = i;
                    break;
                }
            }

            double[][] vs = new double[n][];
            for (int t = 0; t < n; t++)
                vs[t] = new double[n];
            Array.Copy(vec1, vs[0], n);

            /// построение ортогонального дополнения к vec1 путём решения Однор. СЛАУ
            for (int t = 1; t < n; t++)
            {
                for (int i = 0; i < n; i++)
                {
                    if (i == z && t <= z) vs[t][i] = -vs[0][t - 1] / vs[0][z];
                    else
                    if (i == z && t > z) vs[t][i] = -vs[0][t] / vs[0][z];
                    else
                    if (i < z && t == i + 1 || i > z && t == i) vs[t][i] = 1;
                    else
                        vs[t][i] = 0;
                }
            }
            ///

            /// ортогонализация всей совокупности + нормализация
            for (int t = 0; t < n; t++)
            { // не оптимально
                double sum = 0;
                double[] vt_copy = new double[n];
                Array.Copy(vs[t], vt_copy, n);
                for (int i = 0; i < n; i++)
                {
                    double sum_i = vt_copy[i];
                    for (int m = 0; m < t; m++)
                        sum_i -= DotProduct(vt_copy, vs[m]) * vs[m][i];
                    vs[t][i] = sum_i;
                    sum += Math.Pow(sum_i, 2);
                }
                for (int i = 0; i < n; i++) vs[t][i] /= Math.Sqrt(sum); // нормализация
            }
            ///

            double[] best_point = new double[n];
            double best_f = Double.MaxValue;
            double[] pq = new double[n];
            double[] Spq = new double[n];
            for (int q = 0; q < PC; q++)
            {
                double tq = q * Math.PI / PC;
                pq[0] = new_RAs[no].radius * randomd.NextDouble() * Math.Sin(tq) + DotProduct(new_RAs[no].center, vs[0]);
                pq[1] = new_RAs[no].radius * randomd.NextDouble() * Math.Sin(tq) * Math.Cos(tq) + DotProduct(new_RAs[no].center, vs[1]);
                for (int i = 2; i < n; i++)
                    pq[i] = DotProduct(new_RAs[no].center, vs[i]);

                for (int i = 0; i < n; i++)
                {
                    Spq[i] = 0;
                    for (int m = 0; m < n; m++)
                        Spq[i] += vs[m][i] * pq[m];
                }

                double f = function.Call(n, ToAbsolutSize(Spq));
                if (f < best_f)
                {
                    best_f = f;
                    Array.Copy(Spq, best_point, n);
                }
                if (f > f_max) f_max = f;
            }

            Array.Copy(best_point, new_RAs[no].fungi, n);
            new_RAs[no].f = best_f;

            for (int i = 0; i < n; i++)
            {
                if (best_point[i] < 0 || best_point[i] > 1)
                    break;
            }
        }

        double DotProduct(double[] a, double[] b)
        {
            double sum = 0;
            for (int i = 0; i < a.Length; i++) sum += a[i] * b[i];
            return sum;
        }

        void DoSpiral(int no)
        {
            double radius = new_RAs[no].radius;
            double[] vec = new double[n];
            GenVec(vec);

            double[] p1 = new double[n];
            double[] pPC = new double[n];
            double[] pq = new double[n];

            double[] best_point = new double[n];
            double best_f;

            for (int i = 0; i < n; i++)
            {
                pPC[i] = new_RAs[no].center[i];
                p1[i] = new_RAs[no].center[i] + radius * vec[i];
            }

            best_f = function.Call(n, ToAbsolutSize(p1));
            double f = function.Call(n, ToAbsolutSize(pPC));

            if (f > f_max) f_max = f;
            if (best_f > f_max) f_max = f;

            if (f < best_f)
            {
                best_f = f;
                Array.Copy(pPC, best_point, n);
            }
            else Array.Copy(p1, best_point, n);

            for (int q = 1; q < PC - 1; q++)
            {
                double tq = randomd.NextDouble() * 2 - 2;
                for (int i = 0; i < n; i++)
                    pq[i] = radius * Math.Exp(s * tq) * Math.Cos(2 * Math.PI * tq) + pPC[i];

                f = function.Call(n, ToAbsolutSize(pq));
                if (f < best_f)
                {
                    best_f = f;
                    Array.Copy(pq, best_point, n);
                }
                if (f > f_max) f_max = f;
            }

            Array.Copy(best_point, new_RAs[no].fungi, n);
            new_RAs[no].f = best_f;

            for (int i = 0; i < n; i++)
            {
                if (best_point[i] < 0 || best_point[i] > 1)
                    break;
            }
        }

        void SortBasket()
        {
            if (Basket is null)
            {
                Basket = new (double[] fungi, double f)[2 * Np];
                for (int j = 0; j < Np; j++)
                {
                    Basket[j] = (new_RAs[j].fungi, new_RAs[j].f);
                    Basket[j + Np] = (new_RAs[j].fungi, new_RAs[j].f);
                }
            }
            else
            {
                for (int j = 0; j < new_RAs.Count; j++) // в extrasearch мб всего 1 newRA, остальные null
                {
                    Basket[j + Np] = (new_RAs[j].fungi, new_RAs[j].f);
                }
                Array.Sort(Basket, (fungi1, fungi2) => fungi1.f.CompareTo(fungi2.f));
            }
        }

        void UpdateRA()
        {
            for (int t = 0; t < new_RAs.Count; t++)
            {
                RAs.Add(new_RAs[t]);
                new_RAs.RemoveAt(t);
                t--;
            }
        }

        void ExtraSearch(int k)
        {
            int alpha = randomi.Next(0, Np);
            int beta = randomi.Next(0, Np);
            if (alpha == beta) return;

            RA ra = new RA();
            ra.center = new double[n];
            ra.fungi = new double[n];
            ra.f = Double.MaxValue;
            double radius = 0;
            for (int i = 0; i < n; i++)
            {
                ra.center[i] = (Basket[alpha].fungi[i] + Basket[beta].fungi[i]) / 2;
                pickers[0][i] = Basket[alpha].fungi[i];
                radius += Math.Pow(Basket[alpha].fungi[i] - Basket[beta].fungi[i], 2);
                if (pickers[0][i] is double.NaN)
                    break;
            }
            radius = Math.Sqrt(radius) / 2;
            if (radius == 0)
                return; // ничего не делать для одинаковых грибов

            ra.radius = radius;
            new_RAs.Add(ra);

            DoLineExtra(alpha);

            SortBasket();

            new_RAs.Clear();

        }

        bool GlobalEnd()
        {
            double diversity = 0;
            double volume = 1;
            double[] B_ = new double[n];

            for (int i = 0; i < n; i++)
            {
                volume *= 1;
                B_[i] = 0;
                for (int j = 0; j < Np; j++)
                    B_[i] += Basket[j].fungi[i];
                B_[i] /= Np;
            }

            for (int j = 0; j < Np; j++)
            {
                double sum = 0;
                for (int i = 0; i < n; i++)
                    sum += Math.Pow(Basket[j].fungi[i] - B_[i], 2);
                diversity += Math.Sqrt(sum);
            }

            diversity /= Np*volume;

            return diversity < D_eps;
        }

        double CalcRho()
        {
            double rho = 0;
            for (int j = 1; j < Np; j++)
            {
                double sum = 0;
                for (int i = 0; i < n; i++)
                    sum += Math.Pow((Basket[j].fungi[i] - Basket[0].fungi[i]), 2);
                if (rho < Math.Sqrt(sum)) rho = Math.Sqrt(sum); // выбор максимального ро
            }
            return rho;
        }

        void SetSpheres(double rho_k)
        {
            for (int j = 0; j < Np; j++)
            {
                RA ra = new RA();
                ra.n = n;
                ra.center = new double[n];
                Array.Copy(Basket[0].fungi, ra.center, n);
                ra.radius = (1 - Math.Cos((double)(j + 1) / Np * Math.PI / 2)) * rho_k;
                ra.fungi = new double[n];
                ra.zapret = 0;
                new_RAs.Add(ra);
            }
        }

        void WriteToLog(string str, bool append = true)
        {
            using (StreamWriter outputFile = new StreamWriter("Log.txt", append))
            {
                outputFile.WriteLine(str);
            }
        }

        double[] ToAbsolutSize(double[] point)
        {
            double[] new_point = new double[n];
            for (int i = 0; i < n; i++)
                new_point[i] = FS[i][0] + point[i] * (FS[i][1] - FS[i][0]);
            return new_point;
        }
    }
}
