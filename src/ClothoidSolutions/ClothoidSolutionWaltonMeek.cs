using System;
using System.Linq;
using System.Numerics;
using System.Collections.Generic;

namespace ClothoidX
{
    /// <summary>
    /// Solution sets given by D.J. Walton and D.S. Meek
    /// </summary>
    public static class ClothoidSolutionWaltonMeek3
    {
        private static readonly int MAXITER = 10;
        private static readonly double TOL = 1E-7;
        private static void Reverse(ref double phi1, ref double phi2, ref Vector2 P1, ref Vector2 D)
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

        private static Vector2 Rotate(Vector2 v, double radians)
        {
            return new Vector2((v.X * (float)Math.Cos(radians)) - (v.Y * (float)Math.Sin(radians)), (v.X * (float)Math.Sin(radians)) + (v.Y * (float)Math.Cos(radians)));
        }

        internal class FinalParameters
        {
            public Vector2 P0;
            public Vector2 T0;
            public Vector2 N0;
            public double a;
            public double t1;
            public double t2;
            public bool failed;
            public bool reflect;
        }

        /// <summary>
        /// Extent a clothoid curve with a new endpoint and tangent angle (in radians)
        /// </summary>
        /// <param name="curve"></param>
        /// <param name="x">New curve endpoint x</param>
        /// <param name="z">new curve endpoint z</param>
        /// <param name="t">new curve endpoint tangent in radians</param>
        /// <returns></returns>
        public static ClothoidCurve G1Extend(ClothoidCurve curve, double x, double z, double t)
        {
            return new ClothoidCurve();
        }

        /// <summary>
        /// Helper method to get G1 curve using Posture data.
        /// </summary>
        /// <param name="data"></param>
        /// <returns></returns>
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

        public static ClothoidCurve G1Spline(List<Vector3> inputPolyline)
        {
            return G1Spline(Posture.CalculatePostures(inputPolyline).ToArray());
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
        public static ClothoidCurve G1(double x0, double z0, double t0, double x1, double z1, double t1)
        {
            Vector2 P1 = new Vector2((float)x0, (float)z0);
            Vector3 T = ClothoidSegment.RotateAboutAxisRad(Vector3.UnitX, Vector3.UnitY, -t0);
            Vector2 T1 = new Vector2(T.X, T.Z);
            Vector2 P2 = new Vector2((float)x1, (float)z1);
            T = ClothoidSegment.RotateAboutAxisRad(Vector3.UnitX, Vector3.UnitY, -t1);
            Vector2 T2 = new Vector2(T.X, T.Z);

            return G1(P1, T1, P2, T2, .01, 20);
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
        /// <param name="iterLim">maximum newton iterations</param>
        /// <returns></returns>
        public static ClothoidCurve G1(Vector2 P1, Vector2 T1, Vector2 P2, Vector2 T2, double tol, double iterLim)
        {
            Vector3 cacheP1 = new Vector3(P1.X, 0, P1.Y);
            bool reverseFlag = false;
            Vector2 D = P2 - P1;
            float d = D.Length();

            if (d <= tol)
            {
                //degenerate case
                return new ClothoidCurve();
            }
            else
            {
                double phi1 = Math.Atan2(Mathc.Cross(T1, D), Vector2.Dot(T1, D));
                double phi2 = Math.Atan2(Mathc.Cross(D, T2), Vector2.Dot(D, T2));
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
                        v = ClothoidSegment.RotateAboutAxisRad(new Vector3(T1.X, 0, T1.Y), Vector3.UnitY, amt);
                        T1 = new Vector2(v.X, v.Z);
                    }
                    else
                    {
                        v = ClothoidSegment.RotateAboutAxisRad(new Vector3(T2.X, 0, T2.Y), Vector3.UnitY, amt);
                        T2 = new Vector2(v.X, v.Z);
                    }
                    return G1(P1, T1, P2, T2, tol, iterLim);
                }
                else if (Math.Abs(phi1) <= tol && Math.Abs(phi2) <= tol)
                {
                    //straight line segment
                    ClothoidCurve c = new ClothoidCurve().AddLine(D.Length());
                    c.Offset = cacheP1;
                    c.AngleOffset = Math.Atan2(T1.Y, T1.X);
                    return c;
                }
                else if (Math.Abs(phi2 - phi1) <= tol)
                {
                    d = D.Length();
                    //circular segment
                    k = 2 * Math.Sin(phi1) / d;
                    dk = 0;
                    L = d * phi1 / Math.Sin(phi1);

                    if (reverseFlag) k *= -1;
                    if (reflectFlag) k *= -1;

                    ClothoidCurve c = new ClothoidCurve() + new ClothoidSegment(k, dk, L);
                    c.Offset = cacheP1;
                    c.AngleOffset = Math.Atan2(T1.Y, T1.X);
                    return c;
                }
                else
                {
                    FinalParameters f = FitEuler(P1, Vector2.Normalize(D), d, phi1, phi2, reflectFlag);
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

                    ClothoidCurve c = new ClothoidCurve() + new ClothoidSegment(k, dk, L);
                    c.Offset = cacheP1;
                    c.AngleOffset = Math.Atan2(T1.Y, T1.X);
                    return c;
                }

            }
        }

        private static FinalParameters FitEuler(Vector2 P1, Vector2 T, double d, double phi1, double phi2, bool reflect)
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
                if (h > TOL)
                {
                    //solution theta = 0
                }
                else
                {
                    double l = (1 - Math.Cos(phi1)) / (1 - Math.Cos(phi2));
                    double l2 = l * l;
                    theta0 = (l2 / (1 - l2)) * (phi1 + phi2);
                    failed = Solve(0, theta0, phi1, phi2, segno, out theta);
                }
            }
            else
            {
                segno = -1;
                theta0 = Math.Max(0, -phi1);
                theta1 = (Math.PI / 2) - phi1;
                failed = Solve(theta0, theta1, phi1, phi2, segno, out theta);
            }

            t1 = segno * Math.Sqrt(2 * theta / Math.PI);
            t2 = Math.Sqrt(2 * (theta + phi1 + phi2) / Math.PI);
            Mathc.FresnelCS(t1, out double C1, out double S1);
            Mathc.FresnelCS(t2, out double C2, out double S2);
            double phi = phi1 + theta;
            double a = d / (((S2 - S1) * Math.Sin(phi)) + ((C2 - C1) * Math.Cos(phi)));

            Vector2 T0;
            Vector2 N0;
            Vector2 P0;
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

            P0 = P1 - ((float)a * ((float)C1 * T0) + ((float)S1 * N0));

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

        private static bool Solve(double a, double b, double phi1, double phi2, int segno, out double theta)
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
            while (err > TOL && ++u < MAXITER)
            {
                newtonFail = true;
                if (Math.Abs(df) > TOL)
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

            if (u >= MAXITER) return false;
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
    public static class ClothoidSolutionWaltonMeek
    {
        public static ClothoidCurve G1Spline(List<Vector3> inputPolyline)
        {
            ClothoidCurve c = new ClothoidCurve();
            AsymmetricParameters p;
            Vector2 p1;
            Vector2 p2;
            Vector2 p3;
            for (int i = 0; i + 2 < inputPolyline.Count; i++)
            {
                p1 = new Vector2(inputPolyline[i].X, inputPolyline[i].Z);
                p2 = new Vector2(inputPolyline[i + 1].X, inputPolyline[i + 1].Z);
                p3 = new Vector2(inputPolyline[i + 2].X, inputPolyline[i + 2].Z);
                p = AsymmetricClothoid(new Vector2[3] { p1, p2, p3 }, 0.5f);
                //TODO: add t1 and t2 parameters to AsymmetricParameters so we can 
                //calculate L, k and dk to add segments the normal way.

                //c += new ClothoidSegment()
            }

            c.Offset = inputPolyline[0];
            return c;
        }

        private static double Angle(Vector2 v1, Vector2 v2)
        {
            return Math.Acos(Vector2.Dot(v1, v2)) / (v1.Length() * v2.Length());
        }

        public static double ThetaToT(double theta) => Math.Sqrt(2 * theta / Math.PI);
        public static double TToTheta(double t) => t * t * Math.PI / 2;

        public static Vector2[] MakePointsFeasible(Vector2[] points, double tau)
        {
            double alpha = Angle(points[1] - points[0], points[2] - points[1]);
            double t = ThetaToT(alpha);
            double limit = Mathc.C((float)t) / Mathc.S((float)t);
            double a = (points[1] - points[0]).Length();
            double b = (points[2] - points[1]).Length();
            double g;
            double h;

            if (a > b)
            {
                g = a;
                h = b;
            }
            else
            {
                (points[0], points[2]) = (points[2], points[0]);
                g = b;
                h = a;
            }

            Vector2 T0 = (points[1] - points[0]) / (float)g;
            Vector2 newPoint = points[0];

            if (((g / h) + Math.Cos(alpha)) / Math.Sin(alpha) > limit)
            {
                double gLim = h * ((limit * Math.Sin(alpha)) - Math.Cos(alpha));
                double tLerp = ((1 - tau) * h) + (tau * gLim);
                newPoint = points[1] - ((float)tLerp * T0);
            }
            List<Vector2> ps = points.ToList();
            ps.Insert(1, newPoint);
            return ps.ToArray();
        }

        private static double F(double theta, double alpha, double k)
        {
            double t0 = ThetaToT(theta);
            double t1 = ThetaToT(alpha - theta);
            return (Math.Sqrt(theta) * ((Mathc.C(t0) * Math.Sin(alpha)) - (Mathc.S(t0) * (k + Math.Cos(alpha))))) + (Math.Sqrt(alpha - theta) * ((Mathc.S(t1) * (1 + (k * Math.Cos(alpha)))) - (k * Mathc.C(t1) * Math.Sin(alpha))));
        }

        private static double FP(double theta, double alpha, double k)
        {
            double t0 = ThetaToT(theta);
            double t1 = ThetaToT(alpha - theta);
            double st0 = Mathc.S(t0);
            double st1 = Mathc.S(t1);
            double sa = Math.Sin(alpha);
            double ca = Math.Cos(alpha);

            return (st0 * sa / (2 * Math.Sqrt(theta))) *
            ((Mathc.C(t0) / st0) - ((k + ca) / sa)) +
            (k * st1 * sa / (2 * Math.Sqrt(alpha - theta))) *
            ((Mathc.C(t1) / st1) - ((1 + (k * ca)) / (k * sa)));
        }

        private static AsymmetricParameters AsymmetricClothoid(Vector2[] points, double tau)
        {
            double omega = Angle(points[1] - points[0], points[1] - points[2]);
            if (Math.Abs(omega - Math.PI) < 1E-6) throw new InvalidClothoidSegmentException("angle too small");

            Vector2[] p = MakePointsFeasible(points, tau);
            double alpha = Angle(p[2] - p[1], p[3] - p[2]);
            double a = (p[2] - p[1]).Length();
            double b = (p[3] - p[2]).Length();
            double g = Math.Max(a, b);
            double h = Math.Min(a, b);
            double theta0 = alpha / 2;
            double f;
            int u = 0;
            do
            {
                f = F(theta0, alpha, g / h);
                theta0 -= f / FP(theta0, alpha, g / h);
            } while (Math.Abs(f) > 1E-3 && ++u <= 20);

            double theta1 = alpha - theta0;
            double t0 = ThetaToT(theta0);
            double t1 = ThetaToT(theta1);

            double num = g + (h * Math.Cos(alpha));
            double scale = Math.Sqrt(theta1 / theta0);
            double denom = Mathc.C(t0) + (scale * Mathc.C(t1) * Math.Cos(alpha)) + (scale * Mathc.S(t1) * Math.Sin(alpha));
            double a0 = num / denom;
            double a1 = a0 * scale;

            Vector2 P0 = p[1];
            Vector2 T0 = (p[2] - p[1]) / (float)g;
            Vector2 P1 = p[3];
            Vector2 T1 = (p[2] - p[3]) / (float)h;

            Vector2 N0 = new Vector2(-T0.Y, T0.X);
            Vector2 N1 = new Vector2(-T1.Y, T1.X);

            Func<Vector2, Vector2, double> cross = (Vector2 v1, Vector2 v2) => (v1.X * v2.Y) - (v1.Y * v2.X);
            double crossT1T0 = cross(T1, T0);
            if (Math.Sign(cross(T0, N0)) != Math.Sign(crossT1T0)) N0 *= -1;
            if (Math.Sign(cross(N1, T1)) != Math.Sign(crossT1T0)) N1 *= -1;

            return new AsymmetricParameters()
            {
                point0 = points[0],
                point1 = points[1],
                P0 = P0,
                T0 = T0,
                N0 = N0,
                P1 = P1,
                T1 = T1,
                N1 = N1,
                a0 = a0,
                t0 = t0,
                a1 = a1,
                t1 = t1
            };
        }

        internal struct AsymmetricParameters
        {
            public Vector2 point0, point1;
            public Vector2 P0, T0, N0, P1, T1, N1;
            public double a0, t0, a1, t1;
        }

        private static double GetCurvature(double t, double a)
        {
            return Math.Sqrt(2 * Math.PI * t) / a;
        }
    }
    /*
    public class ClothoidSolutionWaltonMeek2 : ClothoidSolution
    {
        private float W { get; set; }

        private List<Posture> postures;

        public override ClothoidCurve CalculateClothoidCurve(List<Vector3> inputPolyline, float allowableError = 0.1f, float endpointWeight = 1)
        {
            this.postures = Posture.CalculatePostures(inputPolyline);
            Vector3 R;
            Vector3 startPosition;
            Vector3 endPosition;
            Vector3 endTangent;
            Posture sp;
            Posture ep;
            for (int i = 0; i + 1 < postures.Count; i++)
            {
                //postures only record the total offsets, we need each pair to be shifted to where the position of the first one is at the origin,
                //and then we need to rotate it so that the initial tangent is aligned with the positive x axis.
                //this solution is only valid when the final tangent angle is less than 180 degrees. Therefore we also need to possibly mirror the
                //solution along the x axis if after shifting and rotating the final point has a negative z value. Once the calculation is finished
                //we can mirror -> rotate -> translate to realign it with the rest of the curve. But the clothoid segment object handles that for us, 
                //so we just need to add clothoid segments in standard position.
                sp = postures[i];
                ep = postures[i + 1];
                //shift the end position by start position
                startPosition = sp.Position - sp.Position;
                endPosition = ep.Position - sp.Position;
                //rotate the end tangent by the start tangent angle
                endTangent = ClothoidSegment.RotateAboutAxisRad(ep.Tangent, Vector3.up, -sp.Angle);
                float endAngle;

                //check if the points can be defined with circular segments or line segments. The clothoid segment object can handle these
                //with just the start and end curvature values so we just check here.
                switch (ClothoidSegment.GetLineTypeFromCurvatureDiff(sp.Curvature, ep.Curvature))
                {
                    case LineType.LINE:
                        break;
                    case LineType.CIRCLE:
                        break;
                    case LineType.CLOTHOID:
                        //check if the postures in standard form are going into the the negative z subplane. if so flag it and mirror them.
                        if (startPosition.z > endPosition.z)
                        {
                            endPosition = new Vector3(endPosition.x, endPosition.y, -endPosition.z);
                            endTangent = new Vector3(endTangent.x, endTangent.y, -endTangent.z);
                        }
                        else
                        {
                        }

                        endAngle = Mathf.Atan2(endTangent.z, endTangent.x);
                        R = CalculateR(endAngle, ep.Curvature);

                        if (IsInGamma(endPosition, R, sp.Curvature, ep.Curvature))
                        {

                        }
                        break;
                }
            }
            throw new NotImplementedException();
        }

        public override List<Vector3> GetFitSamples(int numSaples)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// This returns parameters for the gamma region. 
        /// If startCurvature = 0: 
        /// Item1: slope from origin to point R
        /// Item2: y value of point R
        /// Item3: Vector3.zero
        /// 
        /// If startCurvature > 0:
        /// Item1: slope from origin to point R
        /// Item2: radius of bounding circle
        /// Item3: center of bounding circle
        /// </summary>
        /// <param name="startCurvature"></param>
        /// <param name="endCurvature"></param>
        /// <param name="endTangent">the tangent in radians</param>
        /// <returns></returns>
        public static (float, float, Vector3) GetGammaRegionParameters(float startCurvature, float endCurvature, float endTangent)
        {
            Vector3 R = CalculateR(endTangent, endCurvature);
            float slope = R.z / R.x;
            float radius;
            Vector3 center;
            if (startCurvature > 0)
            {
                radius = (1 / startCurvature) - (1 / endCurvature);
                center = new Vector3(R.x, 0, R.z + (1 / startCurvature) - (1 / endCurvature));
            }
            else
            {
                radius = R.z;
                center = Vector3.zero;
            }

            return (slope, radius, center);
        }

        public static float CalculateW(float endTangent)
        {
            return endTangent * 2 / Mathf.PI;
        }

        public static Vector3 CalculateR(float endTangent, float endCurvature)
        {
            float x = (float)Math.Sin(endTangent);
            float z = 1 - (float)Math.Cos(endTangent);
            //Debug.Log($"X: {x} | Z: {z}");
            return new Vector3(x, 0, z) / endCurvature;
        }


        /// <summary>
        /// Given A and C on a monotonically increasing curve, calculate which side of the curve Q is on.
        /// </summary>
        /// <param name="A">arc length at point A</param>
        /// <param name="C">arc length at point C</param>
        /// <param name="Q">the point we are testing</param>
        /// <returns>0: Q is on the convex side of the curve. 1: Q is on the concave side of the curve. 2: Q is on the curve itself.</returns>
        private int WhichSideD(double A, double C, Vector3 Q)
        {
            //Calculate B
            double B = (A + C) / 2;
            return 0;
        }

        public static bool IsInGamma(Vector3 point, float startCurvature, float endCurvature, float endTangent)
        {
            Vector3 R = CalculateR(endTangent, endCurvature);
            return IsInGamma(point, R, startCurvature, endCurvature);
        }

        /// <summary>
        /// Checks if the point is in an acceptable region where we can find a solution.
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        public static bool IsInGamma(Vector3 point, Vector3 R, float startCurvature, float endCurvature)
        {
            // check if point.z > R.z, the x value is a bit more complicated.
            // Draw a ray from the origin that passes through R (R is in quadrant 1)
            // if our point.x > ray.x at point.z then we are in bounds.
            if (point.z < R.z)
            {
                return false;
            }
            // y = mx + b where b = 0, m = R.y/R.x => x = y / m
            //check if slope of Q is greater than slope of R
            if (point.z / point.x >= R.z / R.x)
            {
                return false;
            }

            if (startCurvature > 0)
            {
                // Additionally check if we are in the defined circle
                float radius = (1 / startCurvature) - (1 / endCurvature);
                Vector3 center = new Vector3(R.x, 0, (1 / startCurvature) - R.z);
                float a = point.x - center.x;
                float b = point.z - center.z;
                if ((a * a) + (b * b) >= radius * radius)
                {
                    return false;
                }
            }

            return true;
        }

        public static List<Vector3> CalculateDCurve(float W, float Kp, float Kq)
        {

            List<Vector3> points = new List<Vector3>();
            int numSamples = 250;
            float Kq2 = Kq * Kq;
            float Kp2 = Kp * Kp;
            float B;
            float KqB;
            float KpB;
            float omega;
            if (Kp > 0)
            {
                for (float t = 0; t < W; t += W / numSamples)
                {
                    B = CalculateB(W, t, Kq2, Kp2);
                    KqB = Kq * B;
                    KpB = Kp * B;
                    omega = -(Kp * Kp * B * B);
                    points.Add(SampleD(t, B, W, omega, Kq, KqB, KpB));
                }
                //add last point manually
                B = CalculateB(W, W, Kq2, Kp2);
                KqB = Kq * B;
                KpB = Kp * B;
                omega = -(KpB * KpB);
                points.Add(SampleD(W, B, W, omega, Kq, KqB, KpB));
            }
            else
            {
                Vector3 R = CalculateR(W * Mathf.PI / 2, Kq);
                for (float t = 0; t < W; t += W / numSamples)
                {
                    points.Add(SampleD2(t, W, Kq, R));
                }
                //add last point manually
                points.Add(SampleD2(W, W, Kq, R));
            }

            return points;
        }

        public static List<Vector3> CalculateECurve(float W, float Kp, float Kq)
        {
            float omega;
            int numPoints = 250;
            List<Vector3> points = new List<Vector3>();
            float B;
            Vector3 a;
            Vector3 b;
            Vector3 S;
            if (Kp > 0)
            {
                for (float t = 0; t < W; t += W / numPoints)
                {
                    B = CalculateB(W, t, Kq * Kq, Kp * Kp);
                    omega = t - (Kp * B * Kp * B);
                    a = GammaProduct(omega, Kp * B, Kq * B, B);
                    b = new Vector3(Mathf.Sin(Mathf.PI * t / 2), 0, 1 - Mathf.Cos(Mathf.PI * t / 2)) / Kp;
                    points.Add(a + b);
                }
                B = CalculateB(W, W, Kq * Kq, Kp * Kp);
                omega = W - (Kp * B * Kp * B);
                a = GammaProduct(omega, Kp * B, Kq * B, B);
                b = new Vector3(Mathf.Sin(Mathf.PI * W / 2), 0, 1 - Mathf.Cos(Mathf.PI * W / 2)) / Kp;
                points.Add(a + b);

            }
            else
            {
                S = SampleD2(0, W, Kq, CalculateR(W * Mathf.PI / 2, Kq));
                points.Add(S);
                points.Add(S + (Vector3.right * 100));
            }
            return points;
        }

        //TODO: Return D(t) for Kp > 0
        /// <summary>
        /// Sample D(t) when Kp > 0
        /// </summary>
        /// <param name="t"></param>
        /// <param name="B"></param>
        /// <param name="W"></param>
        /// <param name="omega"></param>
        /// <param name="Kq"></param>
        /// <param name="KqB"></param>
        /// <param name="KpB"></param>
        /// <returns></returns>
        private static Vector3 SampleD(float t, float B, float W, float omega, float Kq, float KqB, float KpB)
        {
            float x = Mathf.PI * W / 2;
            float y = Mathf.PI * (W - t) / 2;
            float vx = -Mathf.Sin(x) + Mathf.Sin(y);
            float vz = Mathf.Cos(x) - Mathf.Cos(y);
            Vector3 a = GammaProduct(omega, KpB, KqB, B);
            return a - (new Vector3(vx, 0, vz) / Kq);
        }

        /// <summary>
        /// Sample D when Kp = 0
        /// </summary>
        /// <param name="t"></param>
        /// <param name="W"></param>
        /// <param name="Kq"></param>
        /// <returns></returns>
        private static Vector3 SampleD2(float t, float W, float Kq, Vector3 R)
        {
            return R + new Vector3((float)Mathc.C(Mathf.Sqrt(W - t)), 0, (float)Mathc.S(Mathf.Sqrt(W - t))) * Mathf.PI / Kq;
        }

        private bool IsInGammaA(Vector3 point, Vector3 R, float W, float startCurvature, float endCurvature)
        {
            return false;
        }

        private bool IsInGammaB(Vector3 point, Vector3 R, float W, float startCurvature, float endCurvature)
        {
            return false;
        }

        private static double[][] RotationMatrix(float omega)
        {
            float w = Mathf.PI * omega / 2;
            return new double[2][] {
                new double[] {Mathf.Cos(w), -Mathf.Sin(w)},
                new double[] {Mathf.Sin(w), Mathf.Cos(w)}
            };
        }

        private static double[] RotationMatrix(double omega, double x, double y)
        {
            double w = System.Math.PI * omega / 2;
            return new double[2] { (x * Math.Cos(w)) - (y * Math.Sin(w)), (x * Math.Sin(w)) + (y * Math.Cos(w)) };
        }

        private static double[] Vector(float KpB, float KqB)
        {
            return new double[2] { Mathc.C(KqB) - Mathc.C(KpB), Mathc.S(KqB) - Mathc.S(KpB) };
        }

        public static Vector3 GammaProduct(float omega, float KpB, float KqB, float B)
        {
            //double[][] a = Clothoid.Mathc.SVDJacobiProgram.MatProduct(RotationMatrix(omega), Vector(KpB, KqB));
            double[] v = Vector(KpB, KqB);
            double[] r = RotationMatrix(omega, v[0], v[1]);
            return B * Mathf.PI * new Vector3((float)r[0], 0, (float)r[1]);
        }

        private static float CalculateB(float W, float t, float Kq2, float Kp2)
        {
            return Mathf.Sqrt((W - t) / (Kq2 - Kp2));
        }
    }*/
}