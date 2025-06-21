using System;
using System.Numerics;
using System.Collections.Generic;

namespace ClothoidX
{
    public class ClothoidSolutionBertolazziFrego
    {
        /// <summary>
        /// Threshold for A to solve the fresnel equation with different solutions
        /// </summary>
        public static readonly double A_THRESHOLD = 1E-2;
        /// <summary>
        /// When A is small the momenta integrals are evaluated with an infinite series.
        /// Turns out this series converges very fast so we only need the first few terms.
        /// This is that number of terms.
        /// </summary>
        public static readonly int A_SMALL_SERIES_SIZE = 3;
        /// <summary>
        /// Root finding tolerance
        /// </summary>
        private static readonly double ROOT_TOLERANCE = 1E-2;
        /// <summary>
        /// These are coefficients used in the initial guess of A in the root finding algorithm
        /// </summary>
        private static readonly double[] CF = new double[6] { 2.989696028701907, 0.716228953608281, -0.458969738821509, -0.502821153340377, 0.261062141752652, -0.045854475238709 };

        public static double[] SolveG1Parameters(HermiteData point1, HermiteData point2)
        {
            ClothoidCurve c = G1(point1.x, point1.z, point1.tangentAngle, point2.x, point2.z, point2.tangentAngle);
            return new double[4] { c[0].Sharpness, c[0].StartCurvature, c.AngleOffset, c.TotalArcLength };
        }

        public static double[] SolveG1Parameters(ClothoidCurve curve)
        {
            return new double[4] { curve[0].Sharpness, curve[0].StartCurvature, curve.AngleOffset, curve.TotalArcLength };
        }

        public static ClothoidCurve G1Spline(Posture[] data)
        {
            ClothoidCurve c = new ClothoidCurve();

            for (int i = 0; i + 1 < data.Length; i++)
            {
                c += G1(data[i].X, data[i].Z, data[i].Angle, data[i + 1].X, data[i + 1].Z, data[i + 1].Angle);
            }

            c.Offset = data[0].Position;
            c.AngleOffset = data[0].Angle;

            return c;
        }

        /// <summary>
        /// Get a G1 clothoid curve using two points and associated tangents.
        /// Angles should be in radians. Like all G1 and G2 solutions, this 
        /// requires returning a ClothoidCurve object in order to leverage
        /// the Offset and AngleOffset properties, which are required for
        /// the endpoints to interpolate.
        /// </summary>
        /// <param name="x0">start point x</param>
        /// <param name="z0">start point z</param>
        /// <param name="t0">start tangent angle in radians</param>
        /// <param name="x1">end point x</param>
        /// <param name="z1">end point z</param>
        /// <param name="t1">end point tangent angle in radians</param>
        /// <param name="addOffsets">If true, the returned curve will be offet and rotated by the first point and its tangent. We might want to leave it false if we are building a G1 spline, in which case we would only offset and rotate the entire G1 curve by the first point after building the whole thing. See <see cref="G1Spline"/></param>
        /// <returns></returns>
        public static ClothoidCurve G1(double x0, double z0, double t0, double x1, double z1, double t1)
        {
            double dx = x1 - x0;
            double dz = z1 - z0;

            double phi = Math.Atan2(dz, dx);

            double r = Math.Sqrt((dx * dx) + (dz * dz));

            double phi0 = NormalizeAngle(t0 - phi);
            double phi1 = NormalizeAngle(t1 - phi);

            double e = 1E-4;
            //UnityEngine.Debug.Log($"phi0: {phi0} | phi1: {phi1}");
            if ((Math.Abs(phi0) < e && Math.Abs(phi1) == 0) || (phi0 + Math.PI < e && phi1 - Math.PI < e) || (phi0 - Math.PI < e && phi1 + Math.PI < e))
            {
                return ClothoidCurve.FromSegments(new ClothoidSegment(0, 0, (float)r));
            }

            double d = phi1 - phi0;

            //Calculate the bounds of the solution A_max and T_max
            //double Tmax = Math.Max(0, PI_2 + (Math.Sign(phi1) * phi0));
            //double Amax = Tmax == 0 ? absd : absd + (2 * Tmax * (1 + Math.Sqrt(1 + (absd / Tmax))));

            double g;
            double dg;
            List<double[]> IntCS;
            int u = 0;

            double A = InitialGuessA(phi0, phi1);
            //UnityEngine.Debug.LogWarning($"A guess: {A}");
            do
            {
                IntCS = GeneralizedFresnelCS(3, 2 * A, d - A, phi0);
                g = IntCS[1][0];
                dg = IntCS[0][2] - IntCS[0][1];
                A -= g / dg;//
                //UnityEngine.Debug.LogWarning($"u: {u} | g: {g} | dg: {dg} | new A -= {g / dg}: {A}");
            } while (++u < 30 && Math.Abs(g) > ROOT_TOLERANCE);

            if (Math.Abs(g) > ROOT_TOLERANCE)
            {
                //Assert.IsTrue(Math.Abs(g) <= ROOT_TOLERANCE);
                //could not find a root
                //return new ClothoidCurve();
            }
            else
            {
                //UnityEngine.Debug.Log($"Convergence iterations: {u}\nAmax: {Amax}\nTmax: {Tmax}\nA: {A}");
            }

            //UnityEngine.Debug.Log($"Number of attempts: {u}");
            double[] intCS;

            intCS = GeneralizedFresnelCS(2 * A, d - A, phi0);
            double s = r / intCS[0];


            //Assert.IsTrue(s > 0);
            if (s <= 0)
            {
                //UnityEngine.Debug.LogWarning($"ArcLength Negative or Zero! s: {s} | r: {r} | intCS[0]: {intCS[0]} | intCS[1]: {intCS[1]} | phi0: {phi0 * 180 / Math.PI} | phi1: {phi1 * 180 / Math.PI} | A: {A} | d: {d * 180 / Math.PI}");
                return G1CurveM(x0, z0, t0, x1, z1, t1);
                //return new();
            }

            double startCurvature = (d - A) / s;
            double sharpness = 2 * A / (s * s);

            ClothoidCurve c = ClothoidCurve.FromSegments(new ClothoidSegment(startCurvature, sharpness, s));
            c.Offset = new Vector3((float)x0, 0, (float)z0);
            c.AngleOffset = (float)t0;
            return c;
        }

        /// <summary>
        /// Solve the G1 curve by reflecting it across the origin, do this by simply swapping the tangent points and solving the G1 curve again,
        /// then setting the offset to the final point, and the angle offset to t1 + 180.
        /// </summary>
        /// <param name="x0"></param>
        /// <param name="z0"></param>
        /// <param name="t0"></param>
        /// <param name="x1"></param>
        /// <param name="z1"></param>
        /// <param name="t1"></param>
        /// <param name="addOffsets"></param>
        /// <returns></returns>
        public static ClothoidCurve G1CurveM(double x0, double z0, double t0, double x1, double z1, double t1)
        /*{
            ClothoidCurve c = G1Curve(x0, z0, t1, x1, z1, t0, false);
            if (addOffsets)
            {
                c.Offset = new UnityEngine.Vector3((float)x1, 0, (float)z1);
                c.AngleOffset = (float)t1 + 180;
            }
            return c;
        }*/
        {
            double dx = x1 - x0;
            double dz = z1 - z0;

            double phi = Math.Atan2(dz, dx);
            double r = Math.Sqrt((dx * dx) + (dz * dz));

            //NOTE: At the moment I have t1 and t0 swapped to test
            //solving the root for the curve rotated by 180 deg. 
            //Then rotating the solution curve by another 180 deg.
            double phi0 = NormalizeAngle(t1 - phi);
            double phi1 = NormalizeAngle(t0 - phi);

            double e = 1E-4;
            if ((Math.Abs(phi0) < e && Math.Abs(phi1) == 0) || (phi0 + Math.PI < e && phi1 - Math.PI < e) || (phi0 - Math.PI < e && phi1 + Math.PI < e))
            {
                //UnityEngine.Debug.Log("G1 Curve is a line!");
                ClothoidCurve cc = new ClothoidCurve().AddLine((float)r);
                cc.Offset = new Vector3((float)x1, 0, (float)z1);
                cc.AngleOffset = (float)t0;
                return cc;
            }

            double d = phi1 - phi0;

            double g;
            double dg;
            List<double[]> IntCS;
            int u = 0;

            double A = InitialGuessA(phi0, phi1);

            do
            {
                IntCS = GeneralizedFresnelCS(3, 2 * A, d - A, phi0);
                g = IntCS[1][0];
                dg = IntCS[0][2] - IntCS[0][1];
                A -= g / dg;//
                //UnityEngine.Debug.LogWarning($"u: {u} | g: {g} | dg: {dg} | A: {A}");
            } while (++u < 20 && Math.Abs(g) > ROOT_TOLERANCE);

            if (Math.Abs(g) > ROOT_TOLERANCE)
            {
                //UnityEngine.Debug.LogWarning($"No root found! (g, tol, tol2): ({g}, {ROOT_TOLERANCE})");
                //Assert.IsTrue(Math.Abs(g) <= ROOT_TOLERANCE);
                //could not find a root
                return new ClothoidCurve();
            }

            //UnityEngine.Debug.Log($"Number of attempts: {u}");
            double[] intCS;

            intCS = GeneralizedFresnelCS(2 * A, d - A, phi0);
            double s = r / intCS[0];


            //Assert.IsTrue(s > 0);
            if (s <= 0)
            {
                //UnityEngine.Debug.LogWarning($"ArcLength Negative or Zero! s: {s} | r: {r} | intCS[0]: {intCS[0]} | intCS[1]: {intCS[1]} | phi0: {phi0 * 180 / Math.PI} | phi1: {phi1 * 180 / Math.PI} | A: {A} | d: {d * 180 / Math.PI}");

                //intCS = GeneralizedFresnelCS(2 * -A, d + A, phi0);
                //s = r / intCS[0];

                //UnityEngine.Debug.LogWarning($"ArcLength Negative or Zero! s: {s} | r: {r} | intCS[0]: {intCS[0]} | intCS[1]: {intCS[1]} | phi0: {phi0 * 180 / Math.PI} | phi1: {phi1 * 180 / Math.PI} | A: {A} | d: {d * 180 / Math.PI}");
                return new ClothoidCurve();
            }

            double startCurvature = (d - A) / s;
            double sharpness = 2 * A / (s * s);

            ClothoidCurve c = ClothoidCurve.FromSegments(new ClothoidSegment((float)startCurvature, (float)sharpness, (float)s));
            c.Offset = new Vector3((float)x1, 0, (float)z1);
            c.AngleOffset = (float)t1 + 180;
            return c;
        }


        /// <summary>
        /// Normalize an angle in radians to be between -pi and pi.
        /// </summary>
        /// <param name="angle"></param>
        /// <returns></returns>
        private static double NormalizeAngle(double angle)
        {
            while (angle > Math.PI) angle -= 2 * Math.PI;
            while (angle < -Math.PI) angle += 2 * Math.PI;
            return angle;
        }

        private static double InitialGuessA(double phi0, double phi1)
        {
            double X = phi0 / Math.PI;
            double Y = phi1 / Math.PI;
            double xy = X * Y;
            X *= X;
            Y *= Y;
            //double c1 = 3.070645;
            //double c2 = 0.947923;
            //double c3 = -0.673029;
            //return (phi0 + phi1) * (3.070645 + (0.947923 * X * Y) + (-0.673029 * ((X * X) + (Y * Y))));
            return (phi0 + phi1) * (CF[0] + xy * (CF[1] + xy * CF[2]) + (CF[3] + xy * CF[4]) * (X + Y) + CF[5] * (X * X + Y * Y));
        }

        public static double RLommel(double mu, double v, double b)
        {

            double t = 1 / ((mu + v + 1) * (mu - v + 1));
            double r = t;
            double n = 1;
            double e = .00000000001;

            while (Math.Abs(t) > e * Math.Abs(r))
            {
                t *= ((-b) / ((2 * n) + mu - v + 1)) * (b / ((2 * n) + mu + v + 1));
                r += t;
                n++;
            }

            return r;
        }

        /// <summary>
        /// Evaluate XY when a is 0
        /// </summary>
        /// <param name="b"></param>
        /// <param name="k"></param>
        /// <param name="e">The allowable error</param>
        /// <returns></returns>
        private static List<double[]> EvalXYaZero(double b, int k)
        {
            double b2 = b * b;
            double sb = Math.Sin(b);
            double cb = Math.Cos(b);

            double[] X = new double[45];
            double[] Y = new double[45];

            if (Math.Abs(b) < 1E-3)
            {
                X[0] = 1 - (b2 / 6) * (1 - (b2 / 20) * (1 - (b2 / 42)));
                Y[0] = (b2 / 2) * (1 - (b2 / 12) * (1 - (b2 / 30)));
            }
            else
            {
                X[0] = sb / b;
                Y[0] = (1 - cb) / b;
            }

            //int m = (int)Math.Min(Math.Max(1, Math.Floor(2 * b)), k);
            int m = (int)Math.Floor(2 * b);

            if (m >= k) m = k - 1;
            if (m < 1) m = 1;

            for (int i = 1; i < m; ++i)
            {
                X[i] = (sb - (i * Y[i - 1])) / b;
                Y[i] = ((i * X[i - 1]) - cb) / b;
            }

            if (m < k)
            {
                double A = b * sb;
                double D = sb - (b * cb);
                double B = b * D;
                double C = -b2 * sb;
                double rLa = RLommel(m + 0.5, 1.5, b);
                double rLd = RLommel(m + 0.5, 0.5, b);
                double rLb;
                double rLc;

                for (int i = m; i < k; ++i)
                {
                    rLb = RLommel(i + 1.5, 0.5, b);
                    rLc = RLommel(i + 1.5, 1.5, b);
                    //UnityEngine.Debug.LogWarning($"i, m, k: {i}, {m}, {k}");
                    X[i] = ((i * A * rLa) + (B * rLb) + cb) / (1 + i);
                    Y[i] = (((C * rLc) + sb) / (2 + i)) + (D * rLd);
                    rLa = rLc;
                    rLd = rLb;
                    //
                }

            }

            return new List<double[]>() { X, Y };
        }

        private static List<double[]> EvalXYaSmall2(int k, double a, double b, int p)
        {
            if (a == 0) return EvalXYaZero(b, k);

            List<double[]> XY = EvalXYaZero(b, k + (4 * p) + 2);
            double[] X0 = XY[0];
            double[] Y0 = XY[1];
            double[] X = new double[3];
            double[] Y = new double[3];

            double term1;
            double term2;
            double term3;
            double z, zz, zzz;
            int n;
            int nk = 0;
            do
            {
                X[nk] = 0;
                Y[nk] = 0;
                for (n = 0; n <= p; n++)
                {
                    z = Math.Pow(a / 2, 2 * n);
                    zz = Math.Pow(-1, n);
                    zzz = Mathc.Fact(2 * n);
                    term1 = z * zz / zzz;
                    //term1 = Math.Pow(a / 2, 2 * n) * Math.Pow(-1, n) / Mathc.Fact(2 * n);
                    term2 = X0[(4 * n) + k] - (a * Y0[(4 * n) + k + 2] / (2 * ((2 * n) + 1)));
                    term3 = Y0[(4 * n) + k] - (a * X0[(4 * n) + k + 2] / (2 * ((2 * n) + 1)));
                    //UnityEngine.Debug.LogWarning($"{z} | {zz} | {zzz}");

                    if (!double.IsNormal(term1) || !double.IsNormal(term2) || !double.IsNormal(term3))
                    {
                        //UnityEngine.Debug.Log($"n: {n} | term1: {term1} | term2: {term2} | term3: {term3}");
                        //UnityEngine.Debug.Log($"k: {k} | a: {a} | b: {b} | p: {p}");
                    }

                    X[nk] += term1 * term2;
                    Y[nk] += term1 * term3;
                }
            } while (++nk < k);

            return new List<double[]>() { X, Y };
        }

        private static List<double[]> EvalXYaSmall(int k, double a, double b, int p)
        {

            int nkk = k + (4 * p) + 2;
            List<double[]> points = EvalXYaZero(b, nkk);
            double[] X0 = points[0];
            double[] Y0 = points[1];
            double[] X = new double[3];
            double[] Y = new double[3];

            for (int j = 0; j < k; j++)
            {
                X[j] = X0[j] - (a / 2) * Y0[j + 2];
                Y[j] = Y0[j] + (a / 2) * X0[j + 2];
            }

            double t = 1;
            double aa = -a * a / 4;
            double bf;
            int jj;
            for (int n = 1; n <= p; n++)
            {
                t *= aa / (2 * n * ((2 * n) - 1));
                bf = a / ((4 * n) + 2);

                for (int j = 0; j < k; j++)
                {
                    jj = (4 * n) + j;
                    X[j] += t * (X0[jj] - (bf * Y0[jj + 2]));
                    Y[j] += t * (Y0[jj] + (bf * X0[jj + 2]));
                }
            }

            return new List<double[]>() { X, Y };
        }

        private static double[] EvalXYaLarge(double a, double b)
        {
            List<double[]> XY = EvalXYaLarge(0, a, b);
            return new double[2] { XY[0][0], XY[1][0] };
        }

        public static List<double[]> EvalXYaLarge(int n, double a, double b)
        {
            if (n > 3 || n < 0) throw new ArgumentOutOfRangeException();

            double s = a > 0 ? 1 : -1;
            double absa = Math.Abs(a);
            double z = Math.Sqrt(absa) / Math.Sqrt(Math.PI);
            double ell = s * b / Math.Sqrt(Math.PI * absa);
            double g = -s * b * b / (2 * absa);
            double cg = Math.Cos(g) / z;
            double sg = Math.Sin(g) / z;
            List<double[]> CS1 = FresnelCS(n, ell);
            List<double[]> CSz = FresnelCS(n, ell + z);
            double[] C1 = CS1[0];
            double[] S1 = CS1[1];
            double[] Cz = CSz[0];
            double[] Sz = CSz[1];

            double[] dC = new double[n];
            double[] dS = new double[n];

            for (int i = 0; i < n; i++)
            {
                dC[i] = Cz[i] - C1[i];
                dS[i] = Sz[i] - S1[i];
            }
            double[] X = new double[3];
            double[] Y = new double[3];

            //UnityEngine.Debug.LogWarning($"Number of elements in CS1: ({C1.Length}, {S1.Length}) | in CSZ: ({Cz.Length}, {Sz.Length})");

            X[0] = (cg * dC[0]) - (s * sg * dS[0]);
            Y[0] = (sg * dC[0]) + (s * cg * dS[0]);
            if (n > 1)
            {
                cg /= z;
                sg /= z;
                double DC = dC[1] - (ell * dC[0]);
                double DS = dS[1] - (ell * dS[0]);

                X[1] = (cg * DC) - (s * sg * DS);
                Y[1] = (sg * DC) + (s * cg * DS);

                if (n > 2)
                {
                    cg /= z;
                    sg /= z;
                    DC = dC[2] + (ell * ((ell * dC[0]) - (2 * dC[1])));
                    DS = dS[2] + (ell * ((ell * dS[0]) - (2 * dS[1])));
                    X[2] = (cg * DC) - (s * sg * DS);
                    Y[2] = (sg * DC) + (s * cg * DS);
                }
            }

            return new List<double[]>() { X, Y };
        }

        /// <summary>
        /// More closely resembles the paper's mathematical description.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        public static List<double[]> EvalXYaLarge2(int n, double a, double b)
        {
            if (n > 3 || n < 0) throw new ArgumentOutOfRangeException();

            if (double.IsNaN(a)) Console.WriteLine($"a is NaN: {a}");//UnityEngine.Debug.LogWarning($"a is NaN: {a}");

            double[] X = new double[3];
            double[] Y = new double[3];

            double s = a > 0 ? 1 : -1;
            double absa = Math.Abs(a);
            double z = s * Math.Sqrt(absa / Math.PI);
            double wm = b / Math.Sqrt(Math.PI * absa);
            double wp = wm + z;
            double t = -b * b / (2 * a);
            double cz = Math.Cos(t) / z;
            double sz = Math.Sin(t) / z;

            List<double[]> CSp = FresnelCS(3, wp);
            List<double[]> CSm = FresnelCS(3, wm);

            double dC0 = CSp[0][0] - CSm[0][0];
            double dS0 = CSp[1][0] - CSm[1][0];

            X[0] = cz * dC0 - (s * sz * dS0);
            Y[0] = sz * dC0 + (s * cz * dS0);

            if (n > 1)
            {
                cz /= z;
                sz /= z;
                dC0 *= -wm;
                dS0 *= -wm;

                double dC1 = CSp[0][1] - CSm[0][1];
                double dS1 = CSp[1][1] - CSm[1][1];

                X[1] = cz * (dC0 + dC1) - (s * sz * (dS0 + dS1));
                Y[1] = sz * (dC0 + dC1) + (s * cz * (dS0 + dS1));

                if (n > 2)
                {
                    cz /= z;
                    sz /= z;
                    dC0 *= -wm;
                    dS0 *= -wm;
                    dC1 *= -2 * wm;
                    dS1 *= -2 * wm;

                    double dC2 = CSp[0][2] - CSm[0][2];
                    double dS2 = CSp[1][1] - CSm[1][1];

                    X[2] = cz * (dC0 + dC1 + dC2) - (s * sz * (dS0 + dS1 + dS2));
                    Y[2] = sz * (dC0 + dC1 + dC2) + (s * cz * (dS0 + dS1 + dS2));
                }
            }

            return new List<double[]> { X, Y };

        }

        /// <summary>
        /// Evaluate the fresnel integral and its momentae at arc length t
        /// </summary>
        /// <param name="n"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        private static List<double[]> FresnelCS(int n, double t)
        {
            double[] C = new double[n];
            double[] S = new double[n];
            Mathc.FresnelCS(t, out double c, out double s);
            C[0] = c;
            S[0] = s;

            if (n > 1)
            {
                C[1] = C1(t);
                S[1] = S1(t);
                /*
                double tt = Math.PI * t * t / 2;
                double stt = Math.Sin(tt);
                double ctt = Math.Cos(tt);
                C[1] = stt / Math.PI;
                S[1] = (1 - ctt) / Math.PI;*/
                if (n > 2)
                {
                    C[2] = C2(t); //((t * stt) - S[0]) / Math.PI;
                    S[2] = S2(t); //(C[0] - (t * ctt)) / Math.PI;
                }
            }

            return new List<double[]>() { C, S };
        }

        /// <summary>
        /// Compute the Generalized Fresnel integral as described by Bertolazzi and Frego.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        private static List<double[]> GeneralizedFresnelCS(int n, double a, double b, double c)
        {
            //UnityEngine.Debug.Log($"nk: {nk}, a: {a}, b: {b}, c: {c}");
            if (n < 1 || n > 3) throw new ArgumentOutOfRangeException($"Value of index: {n} | Expected value between 1 and 3 inclusive.");
            double cc = Math.Cos(c);
            double sc = Math.Sin(c);
            List<double[]> CS;

            if (Math.Abs(a) < A_THRESHOLD)
            {
                //UnityEngine.Debug.Log("A < THRESHOLD");
                CS = EvalXYaSmall(n, a, b, A_SMALL_SERIES_SIZE);
            }
            else
            {
                //UnityEngine.Debug.Log("A >= THRESHOLD");
                CS = EvalXYaLarge(n, a, b);
            }

            double ci;
            double si;

            for (int i = 0; i < n; i++)
            {
                ci = CS[0][i];
                si = CS[1][i];
                CS[0][i] = (ci * cc) - (si * sc);
                CS[1][i] = (ci * sc) + (si * cc);
            }

            return CS;
        }

        private static double[] GeneralizedFresnelCS(double a, double b, double c)
        {
            List<double[]> CS = GeneralizedFresnelCS(1, a, b, c);
            return new double[] { CS[0][0], CS[1][0] };
        }

        private static double C0(double t)
        {
            return Mathc.C((float)t);
        }

        private static double S0(double t)
        {
            return Mathc.S((float)t);
        }

        private static double C1(double t)
        {
            return Math.Sin(Math.PI * t * t / 2) / Math.PI;
        }

        private static double S1(double t)
        {
            return (1 - Math.Cos(Math.PI * t * t / 2)) / Math.PI;
        }

        private static double C2(double t)
        {
            return ((t * Math.Sin(Math.PI * t * t / 2)) - S0(t)) / Math.PI;
        }

        private static double S2(double t)
        {
            return (C0(t) - (t * Math.Cos(Math.PI * t * t / 2))) / Math.PI;
        }

        /// <summary>
        /// Evaluate the clothoid curve along its arc length using the internal GeneralizedFresnelCS function instead of the built in sampling methods in the ClothoidCurve class. 
        /// </summary>
        /// <param name="arcLength"></param>
        /// <param name="curvature"></param>
        /// <param name="sharpness"></param>
        /// <param name="startAngle"></param>
        /// <param name="start"></param>
        /// <param name="numSamples"></param>
        /// <returns></returns>
        public static List<Vector3> Eval(float arcLength, float curvature, float sharpness, float startAngle, Vector3 start, int numSamples = 100)
        {
            double x0 = start.X;
            double y0 = start.Z;
            List<Vector3> points = new List<Vector3>();
            double[] XY;
            double increment = arcLength / numSamples;
            for (double s = 0; s < arcLength; s += increment)
            {
                XY = GeneralizedFresnelCS(sharpness * s * s, curvature * s, startAngle);
                points.Add(new Vector3((float)(x0 + (XY[0] * s)), 0, (float)(y0 + (XY[1] * s))));
            }
            //add final point
            XY = GeneralizedFresnelCS(sharpness * arcLength * arcLength, curvature * arcLength, startAngle);
            points.Add(new Vector3((float)(x0 + (XY[0] * arcLength)), 0, (float)(y0 + (XY[1] * arcLength))));

            return points;
        }
    }
}