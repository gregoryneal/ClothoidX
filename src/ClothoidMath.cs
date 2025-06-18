using System;
using System.Numerics;
using System.Collections.Generic;

namespace ClothoidX
{

    public static class Mathc
    {
        static readonly double[] fn = new double[11] {
            0.49999988085884732562,
            1.3511177791210715095,
            1.3175407836168659241,
            1.1861149300293854992,
            0.7709627298888346769,
            0.4173874338787963957,
            0.19044202705272903923,
            0.06655998896627697537,
            0.022789258616785717418,
            0.0040116689358507943804,
            0.0012192036851249883877
        };
        static readonly double[] fd = new double[12] {
            1.0,
            2.7022305772400260215,
            4.2059268151438492767,
            4.5221882840107715516,
            3.7240352281630359588,
            2.4589286254678152943,
            1.3125491629443702962,
            0.5997685720120932908,
            0.20907680750378849485,
            0.07159621634657901433,
            0.012602969513793714191,
            0.0038302423512931250065
        };
        static readonly double[] gn = new double[11] {
            0.50000014392706344801,
            0.032346434925349128728,
            0.17619325157863254363,
            0.038606273170706486252,
            0.023693692309257725361,
            0.007092018516845033662,
            0.0012492123212412087428,
            0.00044023040894778468486,
            -8.80266827476172521e-6,
            -1.4033554916580018648e-8,
            2.3509221782155474353e-10
        };
        static readonly double[] gd = new double[12]{
            1.0,
            2.0646987497019598937,
            2.9109311766948031235,
            2.6561936751333032911,
            2.0195563983177268073,
            1.1167891129189363902,
            0.57267874755973172715,
            0.19408481169593070798,
            0.07634808341431248904,
            0.011573247407207865977,
            0.0044099273693067311209,
            -0.00009070958410429993314
        };

        /// <summary>
        /// Compute the Fresnel integral using a method defined by Venkata Sivakanth Telasula ~2005.
        /// </summary>
        /// <param name="arcLength"></param>
        /// <param name="C"></param>
        /// <param name="S"></param>
        /// <exception cref="Exception"></exception>
        public static void FresnelCS(double arcLength, out double C, out double S)
        {
            double x = Math.Abs(arcLength);
            double EPS = 1E-15;
            double PI_2 = Math.PI / 2;

            if (x < 1.0)
            {
                double term, sum;
                double s = PI_2 * (x * x);
                double t = -s * s;

                // Cosine integral series
                double twofn = 0.0;
                double fact = 1.0;
                double denterm = 1.0;
                double numterm = 1.0;
                sum = 1.0;
                do
                {
                    twofn += 2.0;
                    fact *= twofn * (twofn - 1.0);
                    denterm += 4.0;
                    numterm *= t;
                    term = numterm / (fact * denterm);
                    sum += term;
                } while (Math.Abs(term) > EPS * Math.Abs(sum));

                C = x * sum;

                // Sine integral series
                twofn = 1.0;
                fact = 1.0;
                denterm = 3.0;
                numterm = 1.0;
                sum = 1.0 / 3.0;
                do
                {
                    twofn += 2.0;
                    fact *= twofn * (twofn - 1.0);
                    denterm += 4.0;
                    numterm *= t;
                    term = numterm / (fact * denterm);
                    sum += term;
                } while (Math.Abs(term) > EPS * Math.Abs(sum));

                S = PI_2 * sum * (x * x * x);
            }
            else if (x < 6.0)
            {
                // Rational approximation for f
                double sumn = 0.0;
                double sumd = fd[11];
                for (int k = 10; k >= 0; --k)
                {
                    sumn = fn[k] + x * sumn;
                    sumd = fd[k] + x * sumd;
                }
                double f = sumn / sumd;

                // Rational approximation for g
                sumn = 0.0;
                sumd = gd[11];
                for (int k = 10; k >= 0; --k)
                {
                    sumn = gn[k] + x * sumn;
                    sumd = gd[k] + x * sumd;
                }
                double g = sumn / sumd;

                double U = PI_2 * (x * x);
                double SinU = Math.Sin(U);
                double CosU = Math.Cos(U);

                C = 0.5 + f * SinU - g * CosU;
                S = 0.5 - f * CosU - g * SinU;
            }
            else
            {
                // x >= 6.0; asymptotic expansions for f and g
                double s = Math.PI * x * x;
                double t = -1 / (s * s);

                // Expansion for f
                double numterm = -1.0;
                double term = 1.0;
                double sum = 1.0;
                double oldterm = 1.0;
                double eps10 = 0.1 * EPS;
                double absterm;

                do
                {
                    numterm += 4.0;
                    term *= numterm * (numterm - 2.0) * t;
                    sum += term;
                    absterm = Math.Abs(term);
                    if (oldterm < absterm)
                    {
                        throw new Exception($"In FresnelCS f not converged to eps, x = {x}, oldterm = {oldterm}, absterm = {absterm}");
                    }
                    oldterm = absterm;
                } while (absterm > eps10 * Math.Abs(sum));

                double f = sum / (Math.PI * x);

                // Expansion for g
                numterm = -1.0;
                term = 1.0;
                sum = 1.0;
                oldterm = 1.0;

                do
                {
                    numterm += 4.0;
                    term *= numterm * (numterm + 2.0) * t;
                    sum += term;
                    absterm = Math.Abs(term);
                    if (oldterm < absterm)
                    {
                        throw new Exception($"In FresnelCS g not converged to eps, x = {x}, oldterm = {oldterm}, absterm = {absterm}");
                    }
                    oldterm = absterm;
                } while (absterm > eps10 * Math.Abs(sum));

                double g = sum / ((Math.PI * x) * (Math.PI * x * x));

                double U = PI_2 * (x * x);
                double SinU = Math.Sin(U);
                double CosU = Math.Cos(U);

                C = 0.5 + f * SinU - g * CosU;
                S = 0.5 - f * CosU - g * SinU;
            }

            if (arcLength < 0)
            {
                C = -C;
                S = -S;
            }
        }

        /// <summary>
        /// This is the Fresnel Cosine Integral approximation, as given by:
        /// 
        /// HEALD M. A.: Rational approximations for the
        /// fresnel integrals. Mathematics of Computation 44, 170 (1985), 459–461.
        ///
        /// </summary>
        /// <param name="arcLength"></param>
        /// <returns></returns>
        /*public static float C(float arcLength)
        {
            arcLength /= ClothoidSegment.CURVE_FACTOR;
            int scalar = 1;
            if (arcLength < 0)
            {
                scalar = -1;
                arcLength *= -1;
            }
            float val = .5f - (R(arcLength) * (float)System.Math.Sin(System.Math.PI * 0.5f * (A(arcLength) - (arcLength * arcLength))));
            return scalar * val * ClothoidSegment.CURVE_FACTOR;
        }*/

        public static double C(double arcLength)
        {
            arcLength /= ClothoidSegment.CURVE_FACTOR;
            int scalar = 1;
            if (arcLength < 0)
            {
                scalar = -1;
                arcLength *= -1;
            }
            double val = .5f - (R(arcLength) * System.Math.Sin(System.Math.PI * 0.5f * (A(arcLength) - (arcLength * arcLength))));
            return scalar * val * ClothoidSegment.CURVE_FACTOR;
        }

        /// <summary>
        /// The Fresnel Sine Integral approximation, given by HEALD M. A.
        /// </summary>
        /// <param name="arcLength"></param>
        /// <returns></returns>
        /*public static float S(float arcLength)
        {
            arcLength /= ClothoidSegment.CURVE_FACTOR;
            int scalar = 1;
            if (arcLength < 0)
            {
                scalar = -1;
                arcLength *= -1;
            }
            float val = 0.5f - (R(arcLength) * (float)System.Math.Cos(System.Math.PI * 0.5f * (A(arcLength) - (arcLength * arcLength))));
            return scalar * val * ClothoidSegment.CURVE_FACTOR;
        }*/

        public static double S(double arcLength)
        {
            arcLength /= ClothoidSegment.CURVE_FACTOR;
            int scalar = 1;
            if (arcLength < 0)
            {
                scalar = -1;
                arcLength *= -1;
            }
            double val = 0.5f - (R(arcLength) * (float)System.Math.Cos(System.Math.PI * 0.5f * (A(arcLength) - (arcLength * arcLength))));
            return scalar * val * ClothoidSegment.CURVE_FACTOR;
        }

        /// <summary>
        /// Helper function used by the Fresnel Integral approximations
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        private static double R(double t)
        {
            double num = (.506 * t) + 1;
            double denom = (1.79 * t * t) + (2.054 * t) + (float)System.Math.Sqrt(2);
            return num / denom;
        }

        /// <summary>
        /// Helper function used by the Fresnel Integral approximations.
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        private static double A(double t)
        {
            return 1 / ((.803 * t * t * t) + (1.886 * t * t) + (2.524 * t) + 2);
        }

        /// <summary>
        /// Get the tangent vector from an angle in radians.
        /// </summary>
        /// <param name="radians"></param>
        /// <returns></returns>
        public static Vector3 GetTangent(double radians)
        {
            //it is negative because a "positive" angle should rotate in a counterclockwise direction
            return ClothoidSegment.RotateAboutAxisRad(new Vector3(1, 0, 0), Vector3.UnitY, -radians);
        }

        public static Vector3 GetTangentDeg(double degrees)
        {
            return GetTangent(degrees * Math.PI / 180);
        }

        /// <summary>
        /// Approximate the curvature of the circum-circle that passes through all three points.
        /// If they are collinear, return 0, if any of the points are coincident, also return 0.
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <param name="point3"></param>
        /// <returns>A positive only value curvature.</returns>
        public static float MengerCurvature(Vector3 point1, Vector3 point2, Vector3 point3)
        {
            if (point1 == point2 || point2 == point3 || point1 == point3) return 0;
            float a = Vector3.Distance(point1, point2);
            float b = Vector3.Distance(point2, point3);
            float c = Vector3.Distance(point1, point3);
            float s = (a + b + c) / 2; //half perimeter
            float area = (float)Math.Sqrt(s * (s - a) * (s - b) * (s - c)); //herons forumula
            if (area == 0) return 0; //collinear points
            float numerator = 4 * area;
            float denominator = a * b * c;
            return numerator / denominator;
        }

        /// <summary>
        /// Discrete curvature estimation given by Moreton and Sequin in 1992.
        /// This one is the most useful since is distinguishes between positive and 
        /// negative curvature. Positive curvature corresponds to a left turn, and
        /// negative curvature corresponds to a right turn.
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <param name="point3"></param>
        /// <returns></returns>
        public static double MoretonSequinCurvature(Vector3 point1, Vector3 point2, Vector3 point3)
        {
            Vector3 v1 = point2 - point1;
            Vector3 v2 = point3 - point2;
            double det = (v1.X * v2.Z) - (v2.X * v1.Z);
            return 2 * det / (v1.Length() * v2.Length() * (v1 + v2).Length());
        }

        /// <summary>
        /// Basically the same as the absolute value of the MoretonSequinCurvature()
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <param name="point3"></param>
        /// <returns></returns>
        public static double FrenetSerretCurvature(Vector3 point1, Vector3 point2, Vector3 point3)
        {
            Vector3 v1 = point2 - point1;
            Vector3 v2 = point3 - point2;
            float mag1 = v1.Length();
            float mag2 = v2.Length();
            if (mag1 == 0 || mag2 == 0) return 0;
            double theta = Math.Acos(Vector3.Dot(Vector3.Normalize(v1), Vector3.Normalize(v2)));
            return 2 * Math.Sin(theta / 2) / Math.Sqrt(mag1 * mag2);
        }

        /// <summary>
        /// Checks if the three provided points are colinear in the XZ plane.
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <param name="point3"></param>
        /// <returns></returns>
        public static bool AreColinearPoints(Vector3 point1, Vector3 point2, Vector3 point3, float minError = 0.01f)
        {
            float dx1 = point2.X - point1.X;
            float dx2 = point3.X - point2.X;

            //Infinite slope1 and maybe infinite slope 2
            if (dx1 == 0) return dx2 == 0;
            //Infinite slope2 but finite slope1
            else if (dx2 == 0) return false;
            //Slope1 and slope2 are finite
            else return Math.Abs(((point2.Z - point1.Z) / dx1) - ((point3.Z - point2.Z) / dx2)) < minError;
        }

        /// <summary>
        /// Find a circle that goes through three non-colinear points by finding the intersection of the perpendicular bisectors of two cords.
        /// Be sure to handle colinearity seperately.
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <param name="point3"></param>
        /// <returns></returns>
        public static bool CenterOfCircleOfThreePoints(out Vector3 center, Vector3 point1, Vector3 point2, Vector3 point3)
        {
            center = point1;
            //chord of point1 and 2
            (float, float) line1 = EquationOfLineFromTwoPoints(point1, point2);
            //chord of point2 and 3
            (float, float) line2 = EquationOfLineFromTwoPoints(point2, point3);
            //get perpendicular slope
            float slopePerpLine1 = float.IsFinite(line1.Item1) ? -1f / line1.Item1 : 0;
            float slopePerpLine2 = float.IsFinite(line2.Item1) ? -1f / line2.Item1 : 0;
            //get the midpoint of the chord, the line formed by this point and the perpendicular slope passes through the center of the circle
            Vector3 midway1And2 = Vector3.Lerp(point1, point2, 0.5f);
            Vector3 midway2And3 = Vector3.Lerp(point2, point3, 0.5f);
            //find the point where the chord perpendicular bisectors intersect, that is our circle center.
            if (LineLineIntersection(out Vector3 intersection, midway1And2, slopePerpLine1, midway2And3, slopePerpLine2))
            {
                center = intersection;
                return true;
            }
            //Debug.LogError($"Error with chord bisection algorithm => Point1: {point1} | Point2: {point2} | Point3: {point3}");
            return false;
        }

        /// <summary>
        /// Calculate a line intersection based on two points and their respective slopes.
        /// </summary>
        /// <param name="intersection"></param>
        /// <param name="point1"></param>
        /// <param name="slope1"></param>
        /// <param name="point2"></param>
        /// <param name="slope2"></param>
        /// <returns></returns>
        public static bool LineLineIntersection(out Vector3 intersection, Vector3 point1, float slope1, Vector3 point2, float slope2)
        {
            intersection = Vector3.Zero;
            if (Math.Abs(slope1) == Math.Abs(slope2)) return false;

            float int1 = InterceptFromPointAndSlope(point1, slope1);
            float int2 = InterceptFromPointAndSlope(point2, slope2); //y = mx + b => m1x + b1 = m2x + b2 at intersection so x = (b2 - b1) / (m1 - m2)

            float xValue;
            float yValue;

            // the slopes can be infinite and if thats the case the x intercept will be returned instead of the y, we need to check for that.
            if (!float.IsFinite(slope1))
            {
                //slope2 is not infinity, we need the intersection of line two with a vertical line instead
                xValue = point1.X;
                yValue = (slope2 * xValue) + int2;
            }
            else if (!float.IsFinite(slope2))
            {
                xValue = point2.X;
                yValue = (slope1 * xValue) + int1;
            }
            else if (slope1 == 0 && slope2 != 0)
            {
                //if a slope is 0 we need to check for intersection between horizontal line and line with slope
                yValue = point1.Z;
                xValue = (yValue - int2) / slope2;
            }
            else if (slope2 == 0 && slope1 != 0)
            {
                yValue = point2.Z;
                xValue = (yValue - int1) / slope1;
            }
            else if (slope1 == 0 && slope2 == 0)
            {
                if (int1 == int2)
                {
                    xValue = point1.X;
                    yValue = point1.Z;
                }
                else
                {
                    //Debug.LogError($"Slope1: {slope1} | Slope2: {slope2} | Intercept1: {int1} | Intercept2: {int2}");
                    return false;
                } //different horizontal lines
            }
            else
            {
                xValue = (int2 - int1) / (slope1 - slope2);
                yValue = (slope1 * xValue) + int1; //y = mx + b
            }

            /*if (float.IsNaN(xValue) || float.IsNaN(yValue))
            {
                Debug.Log($"Error nan value: intercept 1: {int1} | intercept 2: {int2} | slope 1: {slope1} | slope 2: {slope2} | point 1: {point1} | point 2: {point2}");
            }*/

            intersection = new Vector3(xValue, 0, yValue);

            return true;
        }

        public static bool LineCircleIntersection(out List<Vector3> intersections, Vector3 point, float slope, float circleRadius, Vector3 circleCenter)
        {
            intersections = new List<Vector3>();
            double zIntercept = InterceptFromPointAndSlope(point, slope);
            double discriminant;
            double a;
            double b;
            double c;
            double x;
            double y;
            if (!double.IsFinite(slope))
            {
                //handle case where line is vertical line x = k
                a = 1;
                b = -2 * circleCenter.Z;
                c = (circleCenter.X * circleCenter.X) + (circleCenter.Z * circleCenter.Z) - (circleRadius * circleRadius) - (2 * zIntercept * circleCenter.X) + (zIntercept * zIntercept);

                discriminant = (b * b) - (4 * a * c);

                if (discriminant < 0)
                {
                    return false;
                }
                else if (discriminant > 0)
                {
                    y = (-b + Math.Sqrt(discriminant)) / (2 * a);
                    intersections.Add(new Vector3((float)zIntercept, 0, (float)y));
                    y = (-b - Math.Sqrt(discriminant)) / (2 * a);
                    intersections.Add(new Vector3((float)zIntercept, 0, (float)y));
                    return true;
                }
                else
                {
                    y = -b / (2 * a);
                    intersections.Add(new Vector3((float)zIntercept, 0, (float)y));
                    return true;
                }

            }
            else
            {
                a = (slope * slope) + 1;
                b = 2 * ((slope * zIntercept) - (slope * circleCenter.Z) - circleCenter.X);
                c = (circleCenter.Z * circleCenter.Z) - (circleRadius * circleRadius) + (circleCenter.X * circleCenter.X) - (2 * zIntercept * circleCenter.Z) + (zIntercept * zIntercept);

                discriminant = (b * b) - (4 * a * c);

                if (discriminant < 0)
                {
                    //imaginary solutions 
                    return false;
                }
                else if (discriminant > 0)
                {
                    //2 solutions
                    x = (-b + Math.Sqrt(discriminant)) / (2 * a);
                    intersections.Add(new Vector3((float)x, 0, (float)((slope * x) + zIntercept)));
                    x = (-b - Math.Sqrt(discriminant)) / (2 * a);
                    intersections.Add(new Vector3((float)x, 0, (float)((slope * x) + zIntercept)));
                    return true;
                }
                else
                {
                    //1 solution
                    x = -b / (2 * a);
                    intersections.Add(new Vector3((float)x, 0, (float)((slope * x) + zIntercept)));
                    return true;
                }
            }
        }

        /// <summary>
        /// Returns a pair of floats (x, z) where the x value is the slope and the z value is the z intercept of the line passing through the two points given.
        /// If the slope is infinite, the second value will be the x intercept.
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <returns></returns>
        public static (float, float) EquationOfLineFromTwoPoints(Vector3 point1, Vector3 point2)
        {
            float diffZ = point2.Z - point1.Z;
            float diffX = point2.X - point1.X;
            float slope;

            if (diffX == 0)
            {
                slope = float.PositiveInfinity;
            }
            else
            {
                slope = diffZ / diffX;
            }

            return (slope, InterceptFromPointAndSlope(point1, slope));
        }

        /// <summary>
        /// Returns the intercept from a point and a slope. If the slope is infinity, 
        /// the x intercept will be returned, otherwise the z intercept will be returned.
        /// </summary>
        /// <param name="point"></param>
        /// <param name="slope"></param>
        /// <returns></returns>
        public static float InterceptFromPointAndSlope(Vector3 point, float slope)
        {
            if (!float.IsFinite(slope)) return point.X;
            return point.Z - (slope * point.X); // b = y - mx
        }

        internal class SVDJacobiProgram
        {
            public static void Test()
            {
                //Console.WriteLine("\nBegin SVD decomp via Jacobi algorithm using C#");

                double[][] A = new double[4][];
                A[0] = new double[] { 1, 2, 3 };
                A[1] = new double[] { 5, 0, 2 };
                A[2] = new double[] { 8, 5, 4 };
                A[3] = new double[] { 6, 9, 7 };

                //Console.WriteLine("\nSource matrix: ");
                MatShow(A, 1, 5);

                double[][] U;
                double[][] Vh;
                double[] s;

                //Console.WriteLine("\nPerforming SVD decomposition ");
                SVD_Jacobi(A, out U, out Vh, out s);

                //Console.WriteLine("\nResult U = ");
                MatShow(U, 4, 9);
                //Console.WriteLine("\nResult s = ");
                VecShow(s, 4, 9);
                //Console.WriteLine("\nResult Vh = ");
                MatShow(Vh, 4, 9);

                double[][] S = MatDiag(s);
                double[][] US = MatProduct(U, S);
                double[][] USVh = MatProduct(US, Vh);
                //Console.WriteLine("\nU * S * Vh = ");
                MatShow(USVh, 4, 9);

                //Console.WriteLine("\nEnd demo ");
                Console.ReadLine();
            } // Main

            public static bool SVD_Jacobi(double[][] M, out double[][] U, out double[][] Vh, out double[] s)
            {
                double DBL_EPSILON = 1.0e-15;

                double[][] A = MatCopy(M); // working U
                int m = A.Length;
                int n = A[0].Length;
                double[][] Q = MatIdentity(n); // working V
                double[] t = new double[n];    // working s

                // initialize counters
                int count = 1;
                int sweep = 0;
                //int sweepmax = 5 * n;

                double tolerance = 10 * m * DBL_EPSILON; // heuristic

                // Always do at least 12 sweeps
                int sweepmax = System.Math.Max(5 * n, 12); // heuristic

                // store the column error estimates in St for use
                // during orthogonalization

                for (int j = 0; j < n; ++j)
                {
                    double[] cj = MatGetColumn(A, j);
                    double sj = VecNorm(cj);
                    t[j] = DBL_EPSILON * sj;
                }

                // orthogonalize A by plane rotations
                while (count > 0 && sweep <= sweepmax)
                {
                    // initialize rotation counter
                    count = n * (n - 1) / 2;

                    for (int j = 0; j < n - 1; ++j)
                    {
                        for (int k = j + 1; k < n; ++k)
                        {
                            double cosine, sine;

                            double[] cj = MatGetColumn(A, j);
                            double[] ck = MatGetColumn(A, k);

                            double p = 2.0 * VecDot(cj, ck);
                            double a = VecNorm(cj);
                            double b = VecNorm(ck);

                            double q = a * a - b * b;
                            double v = Hypot(p, q);

                            // test for columns j,k orthogonal,
                            // or dominant errors 
                            double abserr_a = t[j];
                            double abserr_b = t[k];

                            bool sorted = (a >= b);
                            bool orthog = (System.Math.Abs(p) <=
                                tolerance * (a * b));
                            bool noisya = (a < abserr_a);
                            bool noisyb = (b < abserr_b);

                            if (sorted && (orthog ||
                                noisya || noisyb))
                            {
                                --count;
                                continue;
                            }

                            // calculate rotation angles
                            if (v == 0 || !sorted)
                            {
                                cosine = 0.0;
                                sine = 1.0;
                            }
                            else
                            {
                                cosine = System.Math.Sqrt((v + q) / (2.0 * v));
                                sine = p / (2.0 * v * cosine);
                            }

                            // apply rotation to A (U)
                            for (int i = 0; i < m; ++i)
                            {
                                double Aik = A[i][k];
                                double Aij = A[i][j];
                                A[i][j] = Aij * cosine + Aik * sine;
                                A[i][k] = -Aij * sine + Aik * cosine;
                            }

                            // update singular values
                            t[j] = System.Math.Abs(cosine) * abserr_a +
                                System.Math.Abs(sine) * abserr_b;
                            t[k] = System.Math.Abs(sine) * abserr_a +
                                System.Math.Abs(cosine) * abserr_b;

                            // apply rotation to Q (V)
                            for (int i = 0; i < n; ++i)
                            {
                                double Qij = Q[i][j];
                                double Qik = Q[i][k];
                                Q[i][j] = Qij * cosine + Qik * sine;
                                Q[i][k] = -Qij * sine + Qik * cosine;
                            } // i
                        } // k
                    } // j

                    ++sweep;
                } // while

                //    compute singular values
                double prevNorm = -1.0;

                for (int j = 0; j < n; ++j)
                {
                    double[] column = MatGetColumn(A, j);
                    double norm = VecNorm(column);

                    // determine if singular value is zero
                    if (norm == 0.0 || prevNorm == 0.0
                        || (j > 0 &&
                            norm <= tolerance * prevNorm))
                    {
                        t[j] = 0.0;
                        for (int i = 0; i < m; ++i)
                            A[i][j] = 0.0;
                        prevNorm = 0.0;
                    }
                    else
                    {
                        t[j] = norm;
                        for (int i = 0; i < m; ++i)
                            A[i][j] = A[i][j] * 1.0 / norm;
                        prevNorm = norm;
                    }
                }

                U = A;
                Vh = MatTranspose(Q);
                s = t;

                // to sync with default np.linalg.svd() shapes:
                // if m "lt" n, extract 1st m columns of U
                //     extract 1st m values of s
                //     extract 1st m rows of Vh

                if (m < n)
                {
                    U = MatExtractFirstColumns(U, m);
                    s = VecExtractFirst(s, m);
                    Vh = MatExtractFirstRows(Vh, m);
                }

                if (count > 0)
                {
                    //Console.WriteLine("Jacobi iterations did not converge");
                    return false;
                }

                return true;

            } // SVD_Jacobi()

            // === helper functions =================================
            //
            // MatMake, MatCopy, MatIdentity, MatGetColumn,
            // MatExtractFirstColumns, MatExtractFirstRows,
            // MatTranspose, MatDiag, MatProduct, VecNorm, VecDot,
            // Hypot, VecExtractFirst, MatShow, VecShow
            //
            // ======================================================

            public static double[][] MatMake(int r, int c)
            {
                double[][] result = new double[r][];
                for (int i = 0; i < r; ++i)
                    result[i] = new double[c];
                return result;
            }

            static double[][] MatCopy(double[][] m)
            {
                int r = m.Length; int c = m[0].Length;
                double[][] result = MatMake(r, c);
                for (int i = 0; i < r; ++i)
                    for (int j = 0; j < c; ++j)
                        result[i][j] = m[i][j];
                return result;
            }

            static double[][] MatIdentity(int n)
            {
                double[][] result = MatMake(n, n);
                for (int i = 0; i < n; ++i)
                    result[i][i] = 1.0;
                return result;
            }

            static double[] MatGetColumn(double[][] m, int j)
            {
                int rows = m.Length;
                double[] result = new double[rows];
                for (int i = 0; i < rows; ++i)
                    result[i] = m[i][j];
                return result;
            }

            static double[][] MatExtractFirstColumns(double[][] src,
                int n)
            {
                int nRows = src.Length;
                // int nCols = src[0].Length;
                double[][] result = MatMake(nRows, n);
                for (int i = 0; i < nRows; ++i)
                    for (int j = 0; j < n; ++j)
                        result[i][j] = src[i][j];
                return result;
            }

            static double[][] MatExtractFirstRows(double[][] src,
                int n)
            {
                // int nRows = src.Length;
                int nCols = src[0].Length;
                double[][] result = MatMake(n, nCols);
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < nCols; ++j)
                        result[i][j] = src[i][j];
                return result;
            }

            public static double[][] MatTranspose(double[][] m)
            {
                int r = m.Length;
                int c = m[0].Length;
                double[][] result = MatMake(c, r);
                for (int i = 0; i < r; ++i)
                    for (int j = 0; j < c; ++j)
                        result[j][i] = m[i][j];
                return result;
            }

            static double[][] MatDiag(double[] vec)
            {
                int n = vec.Length;
                double[][] result = MatMake(n, n);
                for (int i = 0; i < n; ++i)
                    result[i][i] = vec[i];
                return result;
            }

            public static double[][] MatProduct(double[][] matA, double[][] matB)
            {
                int aRows = matA.Length;
                int aCols = matA[0].Length;
                int bRows = matB.Length;
                int bCols = matB[0].Length;
                if (aCols != bRows)
                {
                    MatShow(matA);
                    MatShow(matB);
                    throw new Exception("Non-conformable matrices");
                }

                double[][] result = MatMake(aRows, bCols);

                for (int i = 0; i < aRows; ++i)
                    for (int j = 0; j < bCols; ++j)
                        for (int k = 0; k < aCols; ++k)
                            result[i][j] += matA[i][k] * matB[k][j];

                return result;
            }

            static double VecNorm(double[] vec)
            {
                double sum = 0.0;
                int n = vec.Length;
                for (int i = 0; i < n; ++i)
                    sum += vec[i] * vec[i];
                return System.Math.Sqrt(sum);
            }

            static double VecDot(double[] v1, double[] v2)
            {
                int n = v1.Length;
                double sum = 0.0;
                for (int i = 0; i < n; ++i)
                    sum += v1[i] * v2[i];
                return sum;
            }

            static double Hypot(double x, double y)
            {
                // fancy sqrt(x^2 + y^2)
                double xabs = System.Math.Abs(x);
                double yabs = System.Math.Abs(y);
                double min, max;

                if (xabs < yabs)
                {
                    min = xabs; max = yabs;
                }
                else
                {
                    min = yabs; max = xabs;
                }

                if (min == 0)
                    return max;

                double u = min / max;
                return max * System.Math.Sqrt(1 + u * u);
            }

            static double[] VecExtractFirst(double[] vec, int n)
            {
                double[] result = new double[n];
                for (int i = 0; i < n; ++i)
                    result[i] = vec[i];
                return result;
            }

            // ------------------------------------------------------

            public static void MatShow(double[][] m, int decimalPlaces = 4, int padWidth = 9)
            {
                for (int i = 0; i < m.Length; ++i)
                {
                    string rowString = "";
                    for (int j = 0; j < m[0].Length; ++j)
                    {
                        double v = m[i][j];
                        if (System.Math.Abs(v) < 1.0e-8) v = 0.0;    // hack
                        rowString += v.ToString("F" + decimalPlaces).PadLeft(padWidth);
                    }
                    //Console.WriteLine(rowString);
                }
            }

            // ------------------------------------------------------

            public static void VecShow(double[] vec, int dec, int wid)
            {
                string rowString = "";
                for (int i = 0; i < vec.Length; ++i)
                {
                    double x = vec[i];
                    if (System.Math.Abs(x) < 1.0e-8) x = 0.0;
                    rowString += x.ToString("F" + dec).PadLeft(wid);
                }
                //Console.WriteLine(rowString);
            }
        }

        /// <summary>
        /// Approximate any single valued integral using the box method and some number of boxes.
        /// This is here because I tried to approximate the actual Fresnel Integrals at first since they are
        /// single valued and simple, but the above approximation is way more accurate and computationally less
        /// expensive. This is here to remind myself how simple some things really are. 
        /// </summary>
        /// <param name="from"></param>
        /// <param name="to"></param>
        /// <param name="function"></param>
        /// <param name="numIntervals">The number of boxes. The higher this number is the more accurate the approximation is.</param>
        /// <returns></returns>
        public static float IntegralApproximation(float from, float to, Func<float, float> function, int numIntervals = 20)
        {
            float integral = 0;
            float dx = (from - to) / numIntervals;
            for (int i = 0; i < numIntervals; i++)
            {
                float x = i * dx;
                integral += function(x) * dx;
            }
            return integral;
        }

        /// <summary>
        /// Approximate a function using simpsons method
        /// </summary>
        /// <param name="from"></param>
        /// <param name="to"></param>
        /// <param name="func"></param>
        /// <returns></returns>
        public static double SimpsonApproximation(double from, double to, Func<double, double> func, int numIntervals)
        {
            if (numIntervals % 2 != 0) numIntervals++;
            double intervalLength = (to - from) / numIntervals;
            double coeff = intervalLength / 6;
            double sum = 0;
            for (double i = from; i + intervalLength <= to; i += intervalLength)
            {
                //new from is i, new to is i+intervalLength
                sum += coeff * (func(i) + (4 * func((i + i + intervalLength) / 2)) + func(i + intervalLength));
            }
            return sum;
        }

        public static BigInteger GetBinCoeff(long N, long K)
        {
            // This function gets the total number of unique combinations based upon N and K.
            // N is the total number of items.
            // K is the size of the group.
            // Total number of unique combinations = N! / ( K! (N - K)! ).
            // This function is less efficient, but is more likely to not overflow when N and K are large.
            // Taken from:  http://blog.plover.com/math/choose.html
            //
            BigInteger r = 1;
            long d;
            if (K > N) return 0;
            for (d = 1; d <= K; d++)
            {
                r *= N--;
                r /= d;
            }
            return r;
        }

        /// <summary>
        /// Factorial an integer
        /// </summary>
        /// <param name="f"></param>
        /// <returns></returns>
        public static int Fact(int f)
        {
            if (f == 0)
                return 1;
            else
                return f * Fact(f - 1);
        }

        /// <summary>
        /// Calculate the cross product of two Vector2s.
        /// </summary>
        /// <param name="rhs"></param>
        /// <param name="lhs"></param>
        /// <returns></returns>
        public static double Cross(Vector2 lhs, Vector2 rhs)
        {
            return (lhs.X * rhs.Y) - (lhs.Y * rhs.X);
        }
    }
}