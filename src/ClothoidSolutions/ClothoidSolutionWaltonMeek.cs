using System;
using System.Collections.Generic;
using System.Numerics;

namespace ClothoidX
{
    /// <summary>
    /// Solution sets given by D.J. Walton and D.S. Meek
    /// </summary>
    public static class ClothoidSolutionWaltonMeek
    {
        private const int MAXITER = 10;
        private const double TOL = 1E-7;
        private static void Reverse(ref double phi1, ref double phi2, ref Mathc.VectorDouble P1, ref Mathc.VectorDouble D)
        {
            D = -D;
            P1 += D;
            (phi1, phi2) = (-phi2, -phi1);
        }

        private static void Reflect(ref double phi1, ref double phi2)
        {
            phi1 = -phi1;
            phi2 = -phi2;
        }

        private static Mathc.VectorDouble Rotate(Mathc.VectorDouble v, double radians)
        {
            //return v.RotateAboutAxis(Mathc.VectorDouble.UnitY, radians);
            return new Mathc.VectorDouble((v.X * Math.Cos(radians)) - (v.Z * Math.Sin(radians)), 0, (v.X * Math.Sin(radians)) + (v.Z * Math.Cos(radians)));
        }

        internal class FinalParameters
        {
            public Mathc.VectorDouble? P0;
            public Mathc.VectorDouble? T0;
            public Mathc.VectorDouble? N0;
            public double a;
            public double t1;
            public double t2;
            public bool failed;
            public bool reflect;
        }

        public static ClothoidCurve G1Spline(Posture[] data)
        {
            ClothoidCurve c = new ClothoidCurve();
            for (int i = 0; i + 1 < data.Length; i++)
            {
                c += G1Segment(data[i].PositionD, data[i].Tangent, data[i + 1].PositionD, data[i + 1].Tangent, TOL, MAXITER);
                //Console.WriteLine($"G1(X0: {data[i].X}, Z0: {data[i].Z}, T0: {data[i].Angle}, X1: {data[i + 1].X}, Z1: {data[i + 1].Z}, T1: {data[i + 1].Angle})");
            }
            /*
            c.Offset = data[0].PositionD;
            c.AngleOffset = data[0].Angle;*/
            return c;
        }

        public static ClothoidCurve G1Spline(List<Vector3> points, bool loop = false)
        {
            return G1Spline(Posture.CalculatePostures(points, loop).ToArray());
        }


        /// <summary>
        /// G1 Spline with ClothoidSegment2 as the output.
        /// </summary>
        /// <param name="P1"></param>
        /// <param name="T1"></param>
        /// <param name="P2"></param>
        /// <param name="T2"></param>
        /// <param name="tol"></param>
        /// <param name="maxIter"></param>
        /// <returns></returns>
        public static ClothoidSegment G1Segment(Mathc.VectorDouble P1, Mathc.VectorDouble T1, Mathc.VectorDouble P2, Mathc.VectorDouble T2, double tol = TOL, int maxIter = MAXITER)
        {
            Mathc.VectorDouble cacheP1 = new Mathc.VectorDouble(P1.X, 0, P1.Z);
            bool reverseFlag = false;
            Mathc.VectorDouble D = P2 - P1;
            double d = D.Length;

            if (d <= tol)
            {
                //degenerate case
                return new ClothoidSegment(Vector3.Zero, Math.Atan2(T1.Z, T1.X), 0, 0, 0);
            }
            else
            {
                double phi1 = Math.Atan2(Mathc.Cross2(T1, D), Mathc.VectorDouble.Dot(T1, D));
                double phi2 = Math.Atan2(Mathc.Cross2(D, T2), Mathc.VectorDouble.Dot(D, T2));
                //Console.WriteLine($"Phi1: {phi1}, Phi2: {phi2}, D: {D.Length}");

                if (Math.Abs(phi1) > Math.Abs(phi2))
                {
                    reverseFlag = true;
                    Reverse(ref phi1, ref phi2, ref P1, ref D);
                }

                bool reflectFlag = false;
                if (phi2 < 0)
                {
                    reflectFlag = true;
                    Reflect(ref phi1, ref phi2);
                }

                double k;
                double dk;
                double L;

                if ((phi1 == 0 && phi2 == Math.PI) || (phi1 == Math.PI && phi2 == 0) || (phi1 == Math.PI && phi2 == Math.PI))
                {
                    //ambiguous case
                    //perturb t1 or t2 and retry
                    Random p = new Random();
                    //perturbation amount
                    double amt = (p.NextDouble() * .01) - .005;
                    Vector3 v;
                    if (p.NextDouble() > .5)
                    {
                        v = Mathc.VectorDouble.RotateAboutAxis(T1, Mathc.VectorDouble.UnitY, amt);
                        T1 = new Mathc.VectorDouble(v.X, 0, v.Z);
                    }
                    else
                    {
                        v = Mathc.VectorDouble.RotateAboutAxis(T2, Mathc.VectorDouble.UnitY, amt);
                        T2 = new Mathc.VectorDouble(v.X, 0, v.Z);
                    }
                    return G1Segment(P1, T1, P2, T2, tol, maxIter);
                }
                else if (Math.Abs(phi1) <= tol && Math.Abs(phi2) <= tol)
                {
                    //straight line segment
                    //TODO: Rotate the line segment to actually interpolate the points between p1 and p2
                    ClothoidSegment s = new ClothoidSegment(cacheP1, Math.Atan2(T1.Z, T1.X), 0, 0, D.Length);
                    return s;
                }
                else if (Math.Abs(phi2 - phi1) <= tol)
                {
                    d = D.Length;
                    //circular segment
                    k = 2 * Math.Sin(phi1) / d;
                    dk = 0;
                    L = d * phi1 / Math.Sin(phi1);

                    if (reverseFlag) k *= -1;
                    if (reflectFlag) k *= -1;

                    ClothoidSegment s = new ClothoidSegment(cacheP1, Math.Atan2(T1.Z, T1.X), k, dk, L);
                    return s;
                }
                else
                {
                    FinalParameters f = FitEuler(P1, D / d, d, phi1, phi2, reflectFlag, tol, maxIter);
                    L = f.a * Math.Abs(f.t2 - f.t1);
                    dk = Math.PI / (f.a * f.a);
                    if (reverseFlag)
                    {
                        k = -Math.PI * f.t2 / f.a;
                    }
                    else
                    {
                        k = Math.PI * f.t1 / f.a;
                    }
                    if (reflectFlag)
                    {
                        k *= -1;
                        dk *= -1;
                    }

                    ClothoidSegment s = new ClothoidSegment(cacheP1, Math.Atan2(T1.Z, T1.X), k, dk, L);
                    return s;
                }

            }
        }

        /// <summary>
        /// The solution to the G1 clothoid interpolation 
        /// problem as given by Walton & Meek in:
        /// 
        /// An Improved Euler Spiral Algorithm for Shape Completion
        /// DOI 10.1109/CRV.2008.11
        /// 
        /// It improves upon work by Kimea, Frankel and Popescu:
        /// Euler Spiral for Shape Completion
        /// DOI 10.1023/A:1023713602895
        /// 
        /// Author: Greg Pritchard 17/06/2025
        /// </summary>
        /// <param name="x0">start position x</param>
        /// <param name="z0">start position z</param>
        /// <param name="t0">start tangent angle radians</param>
        /// <param name="x1">end position x</param>
        /// <param name="z1">end position z</param>
        /// <param name="t1">end tangent angle</param>
        /// <returns></returns>
        public static ClothoidCurve G1Segment(double x0, double z0, double t0, double x1, double z1, double t1)
        {
            Mathc.VectorDouble P1 = new Mathc.VectorDouble(x0, 0, z0);
            Mathc.VectorDouble T = Mathc.VectorDouble.RotateAboutAxis(Mathc.VectorDouble.UnitX, Mathc.VectorDouble.UnitY, -t0);
            Mathc.VectorDouble T1 = new Mathc.VectorDouble(T.X, 0, T.Z);
            Mathc.VectorDouble P2 = new Mathc.VectorDouble(x1, 0, z1);
            T = Mathc.VectorDouble.RotateAboutAxis(Mathc.VectorDouble.UnitX, Mathc.VectorDouble.UnitY, -t1);
            Mathc.VectorDouble T2 = new Mathc.VectorDouble(T.X, 0, T.Z);

            //Console.WriteLine($"G1(X0: {P1.X}, Z0: {P1.Z}, T0: {T1}, X1: {P2.X}, Z1: {P2.Z}, T1: {T2})");

            return G1(P1, T1, P2, T2, TOL, MAXITER);
        }

        /// <summary>
        /// The solution to the G1 clothoid interpolation 
        /// problem as given by Walton & Meek in:
        /// 
        /// An Improved Euler Spiral Algorithm for Shape Completion
        /// DOI 10.1109/CRV.2008.11
        /// 
        /// It improves upon work by Kimea, Frankel and Popescu:
        /// Euler Spiral for Shape Completion
        /// DOI 10.1023/A:1023713602895
        /// 
        /// Author: Greg Pritchard 17/06/2025
        /// </summary>
        /// <param name="P1">the start position</param>
        /// <param name="T1">the start tangent vector</param>
        /// <param name="P2">the end position</param>
        /// <param name="T2">the end tangent vector</param>
        /// <param name="tol">solution tolerance</param>
        /// <param name="maxIter">maximum newton iterations</param>
        /// <returns></returns>
        public static ClothoidCurve G1(Mathc.VectorDouble P1, Mathc.VectorDouble T1, Mathc.VectorDouble P2, Mathc.VectorDouble T2, double tol, int maxIter)
        {
            Mathc.VectorDouble cacheP1 = new Mathc.VectorDouble(P1.X, 0, P1.Z);
            bool reverseFlag = false;
            Mathc.VectorDouble D = P2 - P1;
            double d = D.Length;

            if (d <= tol)
            {
                //degenerate case
                return new ClothoidCurve();
            }
            else
            {
                double phi1 = Math.Atan2(Mathc.Cross2(T1, D), Mathc.VectorDouble.Dot(T1, D));
                double phi2 = Math.Atan2(Mathc.Cross2(D, T2), Mathc.VectorDouble.Dot(D, T2));
                //Console.WriteLine($"Phi1: {phi1}, Phi2: {phi2}, D: {D.Length}");

                if (Math.Abs(phi1) > Math.Abs(phi2))
                {
                    reverseFlag = true;
                    Reverse(ref phi1, ref phi2, ref P1, ref D);
                }

                bool reflectFlag = false;
                if (phi2 < 0)
                {
                    reflectFlag = true;
                    Reflect(ref phi1, ref phi2);
                }

                double k;
                double dk;
                double L;

                if ((phi1 == 0 && phi2 == Math.PI) || (phi1 == Math.PI && phi2 == 0) || (phi1 == Math.PI && phi2 == Math.PI))
                {
                    //ambiguous case
                    //perturb t1 or t2 and retry
                    Random p = new Random();
                    //perturbation amount
                    double amt = (p.NextDouble() * .01) - .005;
                    Vector3 v;
                    if (p.NextDouble() > .5)
                    {
                        v = Mathc.VectorDouble.RotateAboutAxis(T1, Mathc.VectorDouble.UnitY, amt);
                        T1 = new Mathc.VectorDouble(v.X, 0, v.Z);
                    }
                    else
                    {
                        v = Mathc.VectorDouble.RotateAboutAxis(T2, Mathc.VectorDouble.UnitY, amt);
                        T2 = new Mathc.VectorDouble(v.X, 0, v.Z);
                    }
                    return G1(P1, T1, P2, T2, tol, maxIter);
                }
                else if (Math.Abs(phi1) <= tol && Math.Abs(phi2) <= tol)
                {
                    //straight line segment
                    //TODO: Rotate the line segment to actually interpolate the points between p1 and p2
                    ClothoidCurve c = new ClothoidCurve().AddSegment(new ClothoidSegment(cacheP1, Math.Atan2(T1.Z, T1.X), 0, 0, D.Length));
                    //c.Offset = cacheP1;
                    //c.AngleOffset = Math.Atan2(T1.Z, T1.X);
                    return c;
                }
                else if (Math.Abs(phi2 - phi1) <= tol)
                {
                    d = D.Length;
                    //circular segment
                    k = 2 * Math.Sin(phi1) / d;
                    dk = 0;
                    L = d * phi1 / Math.Sin(phi1);

                    if (reverseFlag) k *= -1;
                    if (reflectFlag) k *= -1;

                    ClothoidCurve c = new ClothoidCurve().AddSegment(new ClothoidSegment(cacheP1, Math.Atan2(T1.Z, T1.X), k, dk, D.Length));
                    //c.Offset = cacheP1;
                    //c.AngleOffset = Math.Atan2(T1.Z, T1.X);
                    return c;
                }
                else
                {
                    FinalParameters f = FitEuler(P1, D / d, d, phi1, phi2, reflectFlag, tol, maxIter);
                    L = f.a * Math.Abs(f.t2 - f.t1);
                    dk = Math.PI / (f.a * f.a);
                    if (reverseFlag)
                    {
                        k = -Math.PI * f.t2 / f.a;
                    }
                    else
                    {
                        k = Math.PI * f.t1 / f.a;
                    }
                    if (reflectFlag)
                    {
                        k *= -1;
                        dk *= -1;
                    }

                    ClothoidCurve c = new ClothoidCurve().AddSegment(new ClothoidSegment(cacheP1, Math.Atan2(T1.Z, T1.X), k, dk, D.Length));
                    //c.Offset = cacheP1;
                    //c.AngleOffset = Math.Atan2(T1.Z, T1.X);
                    return c;
                }

            }
        }

        private static FinalParameters FitEuler(Mathc.VectorDouble P1, Mathc.VectorDouble T, double d, double phi1, double phi2, bool reflect, double tol, int maxIter)
        {
            bool failed = true;
            double theta = 0;
            int segno = 1;
            double t1 = 0;
            double t2 = Math.Sqrt(2 * (phi1 + phi2) / Math.PI);
            Mathc.FresnelCS(t2, out double C, out double S);
            double h = (S * Math.Cos(phi1)) - (C * Math.Sin(phi1));

            double theta0;
            double theta1;
            if (phi1 > 0 && h <= 0)
            {
                //C shaped
                if (h > tol)
                {
                    //solution theta = 0
                }
                else
                {
                    double l = (1 - Math.Cos(phi1)) / (1 - Math.Cos(phi2));
                    double l2 = l * l;
                    theta0 = (l2 / (1 - l2)) * (phi1 + phi2);
                    failed = Solve(0, theta0, phi1, phi2, segno, tol, maxIter, out theta);
                }
            }
            else
            {
                segno = -1;
                theta0 = Math.Max(0, -phi1);
                theta1 = (Math.PI / 2) - phi1;
                failed = Solve(theta0, theta1, phi1, phi2, segno, tol, maxIter, out theta);
            }

            t1 = segno * Math.Sqrt(2 * theta / Math.PI);
            t2 = Math.Sqrt(2 * (theta + phi1 + phi2) / Math.PI);
            Mathc.FresnelCS(t1, out double C1, out double S1);
            Mathc.FresnelCS(t2, out double C2, out double S2);
            double phi = phi1 + theta;
            double a = d / (((S2 - S1) * Math.Sin(phi)) + ((C2 - C1) * Math.Cos(phi)));

            Mathc.VectorDouble T0;
            Mathc.VectorDouble N0;
            Mathc.VectorDouble P0;
            if (reflect)
            {
                T0 = Rotate(T, phi);
                N0 = Rotate(T0, -Math.PI / 2);
            }
            else
            {
                T0 = Rotate(T, -phi);
                N0 = Rotate(T0, Math.PI / 2);
            }

            P0 = P1 - (a * (C1 * T0) + (S1 * N0));

            return new FinalParameters()
            {
                P0 = P0,
                T0 = T0,
                N0 = N0,
                a = a,
                t1 = t1,
                t2 = t2,
                failed = failed,
                reflect = reflect,
            };
        }

        private static bool Solve(double a, double b, double phi1, double phi2, int segno, double tol, int maxIter, out double theta)
        {
            Evaluate(a, phi1, phi2, segno, out double fa, out double dfa);
            Evaluate(b, phi1, phi2, segno, out double fb, out double dfb);
            theta = (a + b) / 2;
            Evaluate(theta, phi1, phi2, segno, out double f, out double df);
            double err = b - a;
            int u = 0;
            bool newtonFail;
            double solution;
            double d;
            while (err > tol && ++u < maxIter)
            {
                newtonFail = true;
                if (Math.Abs(df) > tol)
                {
                    solution = theta - (f / df);
                    d = Math.Abs(theta - solution);
                    if (solution > a && solution < b && (d < .5 * err))
                    {
                        newtonFail = false;
                        theta = solution;
                        err = d;
                    }
                }

                if (newtonFail)
                {
                    if (fa * f < 0)
                    {
                        b = theta;
                        fb = f;
                    }
                    else
                    {
                        a = theta;
                        fa = f;
                    }
                    theta = (a + b) / 2;
                    err = b - a;
                }

                Evaluate(theta, phi1, phi2, segno, out f, out df);
            }

            if (u >= maxIter) return false;
            return true;
        }

        private static void Evaluate(double theta, double phi1, double phi2, double segno, out double f, out double df)
        {
            Fresnel2(theta, out double C1, out double S1);
            Fresnel2(theta + phi1 + phi2, out double C2, out double S2);
            double c = Math.Cos(theta + phi1);
            double s = Math.Sin(theta + phi1);
            f = Math.Sqrt(2 * Math.PI) * ((c * (S2 - (segno * S1))) - (s * (C2 - (segno * C1))));
            df = (Math.Sin(phi2) / Math.Sqrt(theta + phi1 + phi2)) + (segno * Math.Sin(phi1) / Math.Sqrt(theta)) - (Math.Sqrt(2 * Math.PI) * ((s * (S2 - (segno * S1))) + (c * (C2 - (segno * C1)))));
        }

        private static void Fresnel2(double t, out double C, out double S)
        {
            Mathc.FresnelCS(Math.Sqrt(2 * Math.Abs(t) / Math.PI), out C, out S);
            C *= Math.Sign(t);
            S *= Math.Sign(t);
        }
    }
}