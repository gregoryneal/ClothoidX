using System;
using System.Collections.Generic;

namespace ClothoidX
{
    /// <summary>
    /// A solution to solve a general clothoid curve with the constraint that every point in the input polyline must be passed through by the resulting curve.
    /// 
    /// It works by creating a Posture object for every sequence of three input nodes, and then drawing a curve between each Posture. We will take left turns to correspond
    /// to negative radii/curvature, and right turns to correspond to positive radii/curvature.
    /// 
    /// This method was described by Dong Hun Shin and Sajiv Singh in their December 1990 paper titled "Path Generation for Robot Vehicles Using Composite Clothoid Segments"
    /// Copyright Carnegie-Mellon University.
    /// </summary>
    public class ClothoidSolutionShinSingh
    {
        /// <summary>
        /// Create 3 clothoid segments for each successive pair of postures depending on the orientation of the two postures.
        /// Let Pi, Ci denote the first posture and circle, and let Pf, Cf denote the second. The two circles will overlap in two places,
        /// marking the start and end nodes of the three clothoid segments. 
        /// 
        /// First, figure out which sign the curvature takes at the start and end node. Visually this involves drawing the Ci and Cf, with 
        /// the tangents visible. If the tangent at Pf points outside of Ci, set the sign of the sharpness of the first clothoid segment so
        /// that it stays within Ci. If the tangent at Pf points inside of Ci, set the sign of the sharpness of the first clothoid segment 
        /// so that it goes outside of Ci. Otherwise set the sharpness to 0 so that the curve follows the circle Ci (Ci = Cf).
        /// 
        /// Now to find the sign of the final (third) clothoid segment we do something very similar. If the tangent of Pi lies outside of Cf,
        /// set the sign of the sharpness of the last clothoid segment so that the curve lies outside Cf. If the tangent at Pi lies inside Cf,
        /// set the sign of the sharpness of the last clothoid segment so that the curve lies inside Cf. Otherwise set the sharpness to 0.
        /// 
        /// Here is a helpful table, Xi and Xf represent the initial and final sharpness sign
        /// and in this solution positive values denote left turns and negative values denote right turns:
        /// 
        /// ki > kf > 0 -> Xi, Xf > 0 
        /// ki > 0 > kf -> Xi, Xf > 0 
        /// 0 > ki > kf -> Xi, Xf < 0 
        /// 0 > kf > ki -> Xi, Xf > 0
        /// kf > 0 > ki -> Xi, Xf < 0
        /// kf > ki > 0 -> Xi, Xf < 0
        /// ki > kf = 0 -> Xi, Xf > 0
        /// 0 = ki > kf -> Xi, Xf > 0
        /// 0 = kf > ki -> Xi, Xf < 0
        /// kf > ki = 0 -> Xi, Xf < 0
        /// kf = ki = 0 -> Xi, Xf = 0
        /// 
        /// Once we figure this out we then need to build three clothoid segments between each posture pair.
        /// They will have the (sharpness, arcLength) values of (x, l1), (-x, l2), (x, l3) where the sign of
        /// x depends on the postures being compared using the table above.
        /// 
        /// We would need to solve fresnel integrals with arc length l1, l2, l3 to find values for x, l1, l2,
        /// and l3. So instead we use numerical methods to compute them.
        /// 
        /// Let l1 and l2 be 1/3 the average length of the (small) circular arcs that connect the initial and final posture.
        /// Those values are:
        /// 
        /// si = ri * theta1 => ri is radius of posture Pi.
        /// theta1 = angle between line CiPi and CiPf (this time Pi and Pf represent the points (Pi.X, Pi.Z) and (Pf.X, Pf.Z))
        /// 
        /// you can derive theta2 similarly as
        /// theta2 = angle between line CfPi and CfPf
        /// 
        /// Therefore set l1 = l2 = (si + sf) / (2 * 3) = (si + sf) / 6
        /// 
        /// Then follow a numerical method iteratively guess the real values of x, l1, l2 and l3:
        /// 
        /// Using l1 and l2, solve for l3, x with the equations:
        /// 
        /// kf = ki + xl1 - xl2 + xl3
        /// -> x = (kf - ki) / (l1 - l2 + l3)
        /// -> l3 = (kf - ki - xl1 + xl2) / x                                                              (1)
        /// a denotes the angle of the posture.
        /// 
        /// af = ai + ki(l1 + l2 + l3) + x(l1l2 - l2l3 + l1l3) + (x/2)(l1^2 - l2^2 + l3^2) 
        /// -> 
        /// x = (af - ai - ki(l1 + l2 + l3)) / ((l1l2 - l2l3 + l1l3) + ((l1^2 - l2^2 + l3^2) / 2))         (2)
        /// 
        /// substitution eq (1) into eq (2) and simplifying yields the following quadratic on x:
        /// 
        /// x^2(l2^2 - l1l2) + x(af - ai - 2l2ki) - (kf^2 - ki^2) / 2 = 0
        /// 
        /// choose the solution with the sign according to the above convention, then substitute back into
        /// eq (1) to solve for l3. Now we have estimated the quadruple (k, l1, l2, l3) which defines the 
        /// 3 clothoid segments connecting the two postures. To check if this solution is valid calculate
        /// the endpoint of the three clothoid segments and see if it nearby the second posture position.
        /// If it is then we have our parameters, if not then we must adjust l1 and l2. To make this numerical
        /// procedure a little faster, we can follow this algorithm:
        /// 
        /// guess l1, l2
        /// -> fix l2 and randomly alter l1 5 times
        ///     -> calculate final position for each pair
        ///     -> choose new l1 that minimizes distance to final posture
        /// -> fix l1 and randomly alter l2 5 times
        ///     -> calculate final position for each pair
        ///     -> choose new l2 that minimizes distance to final posture
        /// -> if distance is within bound exit, otherwise repeat loop.
        /// 
        /// </summary>
        

        /// <summary>
        /// Each posture overlaps its neighbors at two points since they both share two points. Given two neighboring 
        /// Postures, sample points along the small circular arcs between them. 
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static List<Mathc.VectorDouble>[] GetSmallArcsThatConnectPostures(Posture p1, Posture p2)
        {
            List<Mathc.VectorDouble> arc1;
            List<Mathc.VectorDouble> arc2;
            int numSamples = 100;
            double thetao = Math.Atan2(p1.Position.Z - p1.CircleCenter.Z, p1.Position.X - p1.CircleCenter.X);
            double thetaf = Math.Atan2(p2.Position.Z - p1.CircleCenter.Z, p2.Position.X - p1.CircleCenter.X);
            arc1 = p1.GetSamplesD(numSamples, thetaf, thetao);
            thetao = Math.Atan2(p1.Position.Z - p2.CircleCenter.Z, p1.Position.X - p2.CircleCenter.X);
            thetaf = Math.Atan2(p2.Position.Z - p2.CircleCenter.Z, p2.Position.X - p2.CircleCenter.X);
            arc2 = p2.GetSamplesD(numSamples, thetaf, thetao);

            return new List<Mathc.VectorDouble>[2] { arc1, arc2 };
        }

        public Func<double, double> SetupCosTheta1(double thetaI, double curvature, double sharpness)
        {
            return x => Math.Cos(SetupThetaE1(thetaI, curvature, sharpness)(x));
        }

        public Func<double, double> SetupCosTheta2(double thetaI, double curvature, double sharpness, double arcLength1)
        {
            return x => Math.Cos(SetupThetaE2(thetaI, curvature, sharpness, arcLength1)(x));
        }

        public Func<double, double> SetupCosTheta3(double thetaI, double curvature, double sharpness, double arcLength1, double arcLength2)
        {
            return x => Math.Cos(SetupThetaE3(thetaI, curvature, sharpness, arcLength1, arcLength2)(x));
        }

        public Func<double, double> SetupSinTheta1(double thetaI, double curvature, double sharpness)
        {
            return x => Math.Sin(SetupThetaE1(thetaI, curvature, sharpness)(x));
        }

        public Func<double, double> SetupSinTheta2(double thetaI, double curvature, double sharpness, double arcLength1)
        {
            return x => Math.Sin(SetupThetaE2(thetaI, curvature, sharpness, arcLength1)(x));
        }

        public Func<double, double> SetupSinTheta3(double thetaI, double curvature, double sharpness, double arcLength1, double arcLength2)
        {
            return x => Math.Sin(SetupThetaE3(thetaI, curvature, sharpness, arcLength1, arcLength2)(x));
        }

        private Func<double, double> SetupThetaE1(double thetaI, double curvature, double sharpness)
        {
            return x =>
            {
                return thetaI +
                (curvature * x) +
                (sharpness * x * x / 2);
            };
        }

        private Func<double, double> SetupThetaE2(double thetaI, double curvature, double sharpness, double arcLength1)
        {
            return x =>
            {
                return thetaI +
                (curvature * arcLength1) +
                (sharpness * arcLength1 * arcLength1 / 2) +
                ((curvature + (sharpness * arcLength1)) * x) -
                (sharpness * x * x / 2);
            };
        }

        private Func<double, double> SetupThetaE3(double thetaI, double curvature, double sharpness, double arcLength1, double arcLength2)
        {
            return x =>
            {
                return thetaI +
                (curvature * (arcLength1 + arcLength2)) +
                (sharpness * arcLength1 * arcLength2) +
                (sharpness * ((arcLength1 * arcLength1) - (arcLength2 * arcLength2)) / 2) +
                ((curvature + (sharpness * (arcLength1 - arcLength2))) * x) +
                (sharpness * x * x / 2);
            };
        }



        public class Guess : IComparable
        {
            public double s1;
            public double s2;
            public double s3;
            public double x;
            public double ki;
            public double kf;
            public double dist;
            public Guess(double s1, double s2, double s3, double x, double ki, double kf, double dist)
            {
                this.s1 = s1;
                this.s2 = s2;
                this.s3 = s3;
                this.x = x;
                this.ki = ki;
                this.kf = kf;
                this.dist = dist;
            }

            public int CompareTo(object obj)
            {
                if (obj.GetType() != typeof(Guess)) throw new InvalidOperationException();
                return dist.CompareTo(((Guess)obj).dist);
            }
        }


    }
}
