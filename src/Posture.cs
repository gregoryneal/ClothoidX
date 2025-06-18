using System;
using System.Numerics;
using System.Collections.Generic;

namespace ClothoidX
{

    /// <summary>
    /// A Posture is a method to generate target tangents and curvatures for a set of input nodes.
    /// For every subsequent three points a, b, c the posture position P = point b (except for the first
    /// and last posture, which use P = a and P = c respectively). The curvature is the radius of the 
    /// circle C that passes through all three points (it handles collinear points as well), the sign of which
    /// denotes the direction of curvature. The angle is defined by the tangent of C at P. 
    /// </summary>
    public class Posture
    {
        public bool isMirroredX = false;
        public Vector3 Position => new Vector3((float)X, 0, (float)Z);
        public Vector3 Tangent { get; private set; }

        /// <summary>
        /// The point on the circumference of the circle represented by this posture.
        /// </summary>
        public double X { get; private set; }
        public double Z { get; private set; }

        private double _angle = double.NaN;

        /// <summary>
        /// The the tangent angle in radians of the circle at point (x, z)
        /// </summary>
        public double Angle
        {
            get
            {
                if (double.IsNaN(_angle))
                {
                    _angle = Math.Atan2(Tangent.Z, Tangent.X);
                    //while (_angle < 0) _angle += 360;
                }
                return _angle;
            }
        }
        public double Curvature { get; private set; }
        public double Radius { get { return 1 / Curvature; } }
        public Vector3 CircleCenter { get; private set; }
        /*
                public Posture(float x, float z, float angle, float curvature, Vector3 circleCenter) {
                    this.X = x;
                    this.Z = z;
                    this.Angle = angle;
                    this.Curvature = curvature;
                    this.CircleCenter = circleCenter;
                }
        */
        public Posture(double x, double z, double curvature, Vector3 tangent, Vector3 circleCenter)
        {
            this.X = x;
            this.Z = z;
            this.Curvature = curvature;
            this.Tangent = tangent;
            this.CircleCenter = circleCenter;
        }

        public override string ToString()
        {
            return $"position: ({X}, {Z}) | angle: {Angle} | curvature: {Curvature}";
        }

        public List<Vector3> GetSamples(int numSamples, double startAngle = 0, double endAngle = 360)
        {
            //Debug.Log(posture.ToString());
            if (startAngle > endAngle)
            {
                //swap values without introducing temp variable
                (endAngle, startAngle) = (startAngle, endAngle);
            }
            if (Curvature != 0)
            {
                Vector3 sampleCircle = new Vector3((float)Math.Abs(Radius), 0, 0);
                List<Vector3> samples = new List<Vector3>();
                double sampleAngleDeg = (endAngle - startAngle) / numSamples;
                for (double angleDeg = -startAngle; angleDeg > -endAngle; angleDeg -= sampleAngleDeg)
                {
                    samples.Add(ClothoidSegment.RotateAboutAxisDeg(sampleCircle, Vector3.UnitY, (float)angleDeg) + CircleCenter);
                }
                //add last point manually
                samples.Add(ClothoidSegment.RotateAboutAxisDeg(sampleCircle, Vector3.UnitY, -(float)endAngle) + CircleCenter);
                return samples;
            }
            else
            {
                //CircleCenter is the endpoint if the curvature is 0
                return new List<Vector3>() { Position, CircleCenter };
            }
        }

        /// <summary>
        /// Get the Arc Length given an angle in degrees
        /// </summary>
        /// <param name="angleDiff"></param>
        /// <returns></returns>
        public double GetArcLength(double angleDiff)
        {
            //S = R * Theta
            return Math.Abs(Radius * angleDiff) * Math.PI / 180;
        }

        public static List<Posture> CalculatePostures(List<Vector3> polyline)
        {
            List<Posture> postures = new List<Posture>();
            if (polyline.Count < 3) return postures;

            for (int i = 0; i < polyline.Count; i++)
            {
                if (i == 0)
                {
                    postures.Add(CalculatePosture(polyline[0], polyline[1], polyline[2], 1));
                    continue;
                }
                if (i == polyline.Count - 1)
                {
                    postures.Add(CalculatePosture(polyline[^3], polyline[^2], polyline[^1], 3));
                    continue;
                }
                postures.Add(CalculatePosture(polyline[i - 1], polyline[i], polyline[i + 1]));
            }

            return postures;
        }

        /// <summary>
        /// Calulates a posture on point{I}, where point1 is the first point and point3 is the last point of the subsequence in the input polyline.
        /// The Posture is simply a list of useful properties of each point on the polyline, position, heading, curvature, etc. 
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <param name="point3"></param>
        /// <param name="I">The index that should be considered the position of the posture, for tangent calculations.</param>
        /// <returns></returns>
        private static Posture CalculatePosture(Vector3 point1, Vector3 point2, Vector3 point3, int I = 2)
        {
            float x;
            float z;
            if (I == 1)
            {
                x = point1.X;
                z = point1.Z;
            }
            else if (I == 2)
            {
                x = point2.X;
                z = point2.Z;
            }
            else
            {
                x = point3.X;
                z = point3.Z;
            }

            double curvature;
            Vector3 tangent;
            Vector3 circleCenter;
            bool areCollinear = Mathc.AreColinearPoints(point1, point2, point3);
            bool foundCenter = false;
            if (areCollinear)
            {
                curvature = 0;
                //tangent = point3 - point1;
                //use the adjacent point to calculate tangent and final position
                if (I == 1)
                {
                    circleCenter = point2;
                    tangent = point2 - point1;
                }
                else if (I == 2)
                {
                    circleCenter = point3;
                    tangent = point3 - point2;
                }
                else
                {
                    circleCenter = point2;
                    tangent = point3 - point2;
                }
            }
            else
            {
                curvature = Mathc.MoretonSequinCurvature(point1, point2, point3);
                if (Mathc.CenterOfCircleOfThreePoints(out Vector3 center, point1, point2, point3))
                {
                    circleCenter = center;
                    //tangent slope is negative reciprocal of normal slope since they are orthogonal.
                    tangent = new Vector3(z - circleCenter.Z, 0, -(x - circleCenter.X));
                    //Debug.Log($"New Tangent: {tangent}");
                    if (curvature > 0) tangent *= -1;

                }
                else
                {
                    //we shouldn't get here so we can use this as debug.
                    Console.WriteLine($"Warning, something is wrong with the center of the circle calculations! Points: ({point1.X}, {point1.Z}) | ({point2.X}, {point2.Z}) | ({point3.X}, {point3.Z})  || Colinear: {areCollinear}");
                    circleCenter = point1;
                    tangent = Vector3.Zero;
                }
            }

            //angle = Mathf.Abs(angle);

            if (foundCenter)
            {
                //Debug.Log($"We found the center! point1: {point1} | point2: {point2} | point3: {point3} | areCollinear: {areCollinear} | curvature: {curvature} | circleCenter: {circleCenter} | angle: {angle}");
            }

            return new Posture(x, z, curvature, tangent, circleCenter);
        }

        /// <summary>
        /// Given two postures, return two of the same postures but in standard form. 
        /// The first Posture will be centered on the origin with its tangent aligned with the +X-Axis.
        /// The second posture will be relative to that. If the second posture ends up in the fourth 
        /// quadrant, it will be mirrored along the X-Axis and it will be flagged isMirroredX.
        /// </summary>
        /// <param name="postures"></param>
        /// <returns></returns>
        public static (Posture, Posture) InStandardForm(Posture posture1, Posture posture2)
        {
            Posture newPosture1;
            Posture newPosture2;

            //Debug.Log($"Rotate by {posture1.Angle} degrees");
            //Debug.Log($"Rotate 2 by {posture2.Angle} degrees");
            //rotate relative positions by the start tangent angle
            Vector3 endPosition = ClothoidSegment.RotateAboutAxisRad(posture2.Position - posture1.Position, Vector3.UnitY, posture1.Angle);
            Vector3 circle1 = ClothoidSegment.RotateAboutAxisRad(posture1.CircleCenter - posture1.Position, Vector3.UnitY, posture1.Angle);
            Vector3 circle2 = ClothoidSegment.RotateAboutAxisRad(posture2.CircleCenter - posture1.Position, Vector3.UnitY, posture1.Angle);
            Vector3 endTangent = ClothoidSegment.RotateAboutAxisRad(posture2.Tangent, Vector3.UnitY, posture1.Angle);

            newPosture1 = new Posture(0, 0, Math.Abs(posture1.Curvature), Vector3.UnitX * posture1.Tangent.Length(), circle1);
            newPosture2 = new Posture(endPosition.X, Math.Abs(endPosition.Z), Math.Abs(posture2.Curvature), endTangent, circle2);

            if (endPosition.Z < 0)
            {
                //Debug.Log($"IS MIRRORED: {endPosition}");
                newPosture2.isMirroredX = true;
                newPosture2.Tangent = new Vector3(newPosture2.Tangent.X, newPosture2.Tangent.Y, -newPosture2.Tangent.Z);

                newPosture1.CircleCenter = new Vector3(newPosture1.CircleCenter.X, newPosture1.CircleCenter.Y, -newPosture1.CircleCenter.Z);
                newPosture2.CircleCenter = new Vector3(newPosture2.CircleCenter.X, newPosture2.CircleCenter.Y, -newPosture2.CircleCenter.Z);
            }

            return (newPosture1, newPosture2);
        }

        public static IEnumerator<(Posture, Posture)> GetStandardPostures(List<Posture> postures)
        {
            for (int i = 0; i + 1 < postures.Count; i++)
            {
                yield return InStandardForm(postures[i], postures[i + 1]);
            }
        }
    }
}