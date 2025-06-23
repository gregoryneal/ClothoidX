using System;
using System.Collections.Generic;
using System.Numerics;

namespace ClothoidX
{
    /// <summary>
    /// A version of ClothoidSegment with an updated sampling method that ensures G0 continuity of all input polyline nodes.
    /// </summary>
    [Serializable]
    public class ClothoidSegment2
    {
        public static double CURVE_FACTOR = Math.Sqrt(Math.PI / 2);
        //if the curvature is less than this we will approximate it as a line segment.
        public static double MIN_CURVATURE_DIFF = 1E-15;
        /// <summary>
        /// Represents the arc length on the parent ClothoidCurve where this segment starts.
        /// </summary>
        public double ArcLengthStart { get; protected set; }
        /// <summary>
        /// Arc length on the parent ClothoidCurve where this segment ends.
        /// </summary>
        public double ArcLengthEnd { get; protected set; }
        /// <summary>
        /// The total arc length of this segment.
        /// </summary>
        public double TotalArcLength => ArcLengthEnd - ArcLengthStart;
        /// <summary>
        /// Curvature in radians at the start of the segment. Positive -> right hand turn, negative -> left hand turn.
        /// </summary>
        public double StartCurvature { get; protected set; }
        protected double _StartCurvature { get; set; }
        /// <summary>
        /// Curvature in radians at the end of the segment. Positive -> right hand turn, negative -> left hand turn.
        /// </summary>
        public double EndCurvature { get; protected set; }
        protected double _EndCurvature { get; set; }
        /// <summary>
        /// Scale factor for the entire curve, can be negative.
        /// </summary>
        public double B { get; protected set; }
        /// <summary>
        /// The scale factor as it was originally calculated, so we don't forget.
        /// </summary>
        protected double _B { get; set; }

        public double Sharpness
        {
            get
            {
                if (EndCurvature == StartCurvature) return 0;
                else if (TotalArcLength > 0) return (EndCurvature - StartCurvature) / TotalArcLength;
                else return 0;
            }
        }

        /// <summary>
        /// The xz offset of the endpoint of this segment, relative to the origin.
        /// </summary>
        public Vector3 RelativeEndpoint { get; private set; }
        /// <summary>
        /// The start position of this segment.
        /// </summary>
        public Vector3 Position { get; private set; }

        public Vector3 Endpoint { get => Position + RelativeEndpoint; }

        /// <summary>
        /// The angle of the endpoint tangent, in radians.
        /// </summary>
        public double Rotation { get; private set; }

        /// <summary>
        /// What type of segment is this? A line, a circular arc, or a clothoid segment.
        /// </summary>
        public LineType LineType { get; private set; }

        /// <summary>
        /// This determines the clothoid segment parameters. All clothoids are calculated first in "standard" form. 
        /// That is, from the origin with the initial tangent positioned along the positive x axis. 
        /// To connect them we must apply a rotation, a mirror, and then a translation.
        /// </summary>
        /// <param name="arcLengthStart"></param>
        /// <param name="arcLengthEnd"></param>
        /// <param name="startCurvature"></param>
        /// <param name="endCurvature"></param>
        /// <param name="B"></param>
        public ClothoidSegment2(Vector3 position, double rotation, double arcLengthStart, double arcLengthEnd, double startCurvature, double endCurvature, double B)
        {
            this.ArcLengthStart = arcLengthStart;
            this.ArcLengthEnd = arcLengthEnd;
            this.StartCurvature = startCurvature;
            this._StartCurvature = startCurvature;
            this.EndCurvature = endCurvature;
            this._EndCurvature = endCurvature;
            this.B = B;
            this._B = B;
            this.LineType = GetLineTypeFromCurvatureDiff(StartCurvature, EndCurvature);
            this.Position = position;
            this.Rotation = rotation;
            //Console.WriteLine($"I am a {this.LineType} - segment!");
            CalculateOffsetAndRotation();
        }

        public ClothoidSegment2(Vector3 position, double rotation, double arcLengthStart, double arcLengthEnd, double startCurvature, double endCurvature) : this(position, rotation, arcLengthStart, arcLengthEnd, startCurvature, endCurvature, CalculateB(arcLengthStart, arcLengthEnd, startCurvature, endCurvature)) { }

        public ClothoidSegment2(Vector3 position, double rotation, double startCurvature, double sharpness, double totalArcLength) : this(position, rotation, 0, totalArcLength, startCurvature, (sharpness * totalArcLength) + startCurvature) { }

        public ClothoidSegment2(ClothoidSegment2 segment) : this(segment.Position, segment.Rotation, segment.StartCurvature, segment.Sharpness, segment.TotalArcLength) { }

        /// <summary>
        /// Calcluate the scaling factor for the constrained SinghMcCrae clothoid segment.
        /// </summary>
        /// <returns></returns>
        public static double CalculateB(double arcLengthStart, double arcLengthEnd, double startCurvature, double endCurvature)
        {
            double arcLength = arcLengthEnd - arcLengthStart;
            double B = Math.Sqrt(arcLength / (Math.PI * Math.Abs(endCurvature - startCurvature)));
            return B;
        }

        /// <summary>
        /// Sample this segment instance in local space (without offsets and rotations added).
        /// </summary>
        /// <param name="arcLength"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        public Vector3 SampleSegmentByTotalArcLength(double arcLength)
        {
            Vector3 v;
            if (arcLength > ArcLengthEnd || arcLength < ArcLengthStart) throw new ArgumentOutOfRangeException();
            double interp = (arcLength - ArcLengthStart) / TotalArcLength;
            switch (this.LineType)
            {
                case LineType.LINE:
                    v = new Vector3((float)(arcLength - ArcLengthStart), 0, 0);
                    break;
                case LineType.CIRCLE:
                    //we build the circle constructively to find the point and rotation angle.
                    //start with a positive or negative radius, and build a vector centered on the origin with the z value as the radius.
                    //rotate the point by the desired theta, positive theta rotates in the clockwise direction, negative values are ccw.
                    //now subtract the final z value by the radius (point.z -= radius) to get the final offset. 
                    double radius = -2f / (StartCurvature + EndCurvature); //this might be positive or negative
                    double circumference = 2 * Math.PI * radius;
                    double fullSweepAngle_deg = 360f * TotalArcLength / circumference;
                    double rotationAngle = interp * fullSweepAngle_deg; //same here, value in degrees
                    //Vector3 vector = RotateAboutAxisDeg(new Vector3(0, 0, (float)radius), Vector3.UnitY, rotationAngle);
                    Vector3 vector = RotateAboutAxisRad(new Vector3(0, 0, (float)radius), Vector3.UnitY, interp * TotalArcLength / radius);
                    v = new Vector3(vector.X, vector.Y, (float)(vector.Z - radius));
                    //return new Vector3(vector.x, vector.y, isMirroredX ? radius - vector.z : vector.z - radius);
                    break;
                default:
                    //clothoid
                    v = this.SampleClothoidLocal(interp);
                    break;
            }

            //if (Sharpness == 0) v = new Vector3(v.x, v.y, -v.z);
            return v;
        }

        public Mathc.VectorDouble SampleSegmentByTotalArcLengthDouble(double arcLength)
        {
            if (arcLength > ArcLengthEnd || arcLength < ArcLengthStart) throw new ArgumentOutOfRangeException();
            double interp = (arcLength - ArcLengthStart) / TotalArcLength;
            switch (this.LineType)
            {
                case LineType.LINE:
                    return new Mathc.VectorDouble(arcLength - ArcLengthStart, 0, 0);
                case LineType.CIRCLE:
                    double radius = -2f / (StartCurvature + EndCurvature); //this might be positive or negative
                    double circumference = 2 * Math.PI * radius;
                    //double fullSweepAngle_deg = 360.0 * TotalArcLength / circumference;
                    //double rotationAngle = interp * fullSweepAngle_deg; //same here, value in radians
                    double rot = interp * TotalArcLength / radius;
                    Mathc.VectorDouble vector = Mathc.VectorDouble.RotateAboutAxis(new Mathc.VectorDouble(0, 0, radius), Mathc.VectorDouble.UnitY, rot);
                    return new Mathc.VectorDouble(vector.X, vector.Y, vector.Z - radius);
                default:
                    return SampleClothoidD(interp);
            }
        }

        /// <summary>
        /// Shifts the start arc length to newStartArcLength, shifts the end arc length value to newStartArcLength + TotalArcLength of the segment at creation.
        /// </summary>
        /// <param name="newStartArcLength"></param>
        public void ShiftStartArcLength(double newStartArcLength)
        {
            double arcLength = TotalArcLength;
            this.ArcLengthStart = newStartArcLength;
            this.ArcLengthEnd = newStartArcLength + arcLength;
        }

        /// <summary>
        /// Sample this segment by a value constrained between 0 and 1.
        /// </summary>  
        /// <param name="value"></param>
        /// <returns></returns>
        public Vector3 SampleSegmentByRelativeLength(double value)
        {
            if (value < 0 || value > 1) throw new ArgumentOutOfRangeException();
            double totalArcLength = ArcLengthStart + (value * (ArcLengthEnd - ArcLengthStart));
            return SampleSegmentByTotalArcLength(totalArcLength);
        }

        public List<Vector3> GetSamples(int numSamples)
        {
            List<Vector3> points = new List<Vector3>();
            for (int i = 0; i < numSamples - 1; i++)
            {
                points.Add(SampleSegmentByRelativeLength(i / (numSamples - 1)));
            }
            //do the endpoint manually
            points.Add(SampleSegmentByRelativeLength(1));
            return points;
        }

        /// <summary>
        /// Get the tangent at a specific arc length, in radians.
        /// </summary>
        /// <param name="arcLength"></param>
        /// <returns></returns>
        public double Tangent(double arcLength)
        {
            switch (this.LineType)
            {
                case LineType.LINE:
                    Console.WriteLine("Line segment");
                    return 0;
                case LineType.CIRCLE:
                    Console.WriteLine("Circle segment");
                    double radius = 2f / (StartCurvature + EndCurvature); //this might be positive or negative
                    double circumference = 2 * Math.PI * radius;
                    return 2 * Math.PI * arcLength / circumference;
                default:
                    Console.WriteLine("Clothoid segment");
                    return (Sharpness * arcLength * arcLength / 2) + (StartCurvature * arcLength);
            }
        }

        /// <summary>
        /// Reflect the segment along the common axis between the start and end points.
        /// Note, this will affect the ClothoidCurve that contains this segment,
        /// it will affect the endpoint tangent where this segment ends.
        /// </summary>
        public void Reflect()
        {
            StartCurvature *= -1;
            EndCurvature *= -1;
            CalculateOffsetAndRotation();
        }

        /// <summary>
        /// Reverse the segment. Note this affect the ClothoidCurve at
        /// the point where this segment begins.
        /// </summary>
        public void Reverse()
        {
            (StartCurvature, EndCurvature) = (-EndCurvature, -StartCurvature);
            CalculateOffsetAndRotation();
        }

        /// <summary>
        /// Here we calculate the offset and rotation for the segment endpoint. All of the shapes can be thought of as being generated from the origin.
        /// Lines are generated in the positive X direction, circles with negative radius (Left turning) are generate in the positive z planar subsection,
        /// centered on the y axis. Clothoids are generated in the positive X direction, but can turn into the positive or negative z direction.
        /// 
        /// Rotations are calculated as follows: for lines there is no rotation, so it is 0. For circles it is calculated by the ratio 360deg * (arcLength / circumference).
        /// Clothoids are calculated with (EndCurvature - StartCurvature) * ArcLengthEnd * ArcLengthEnd / (2 * (ArcLengthEnd - ArcLengthStart)) 
        /// </summary>
        protected void CalculateOffsetAndRotation()
        {
            switch (this.LineType)
            {
                case LineType.LINE:
                    this.RelativeEndpoint = new Vector3((float)TotalArcLength, 0, 0);
                    this.Rotation = 0;
                    break;
                case LineType.CIRCLE:
                    //we build the circle constructively to find the point and rotation angle.
                    //start with a positive or negative radius, and build a vector centered on the origin with the z value as the radius.
                    //rotate the point by the desired theta, positive theta rotates in the clockwise direction, negative values are ccw.
                    //now subtract the final z value by the radius (point.z -= radius) to get the final offset. 
                    double radius = -2f / (StartCurvature + EndCurvature); //this might be positive or negative

                    bool negativeCurvature = radius < 0;
                    double circumference = 2f * Math.PI * Math.Abs(radius);
                    double rotationAngle = TotalArcLength * 2 * Math.PI / circumference; //same here, value in degrees
                    Vector3 vector;

                    if (!negativeCurvature)
                    {
                        vector = RotateAboutAxisRad(new Vector3(0, 0, (float)radius), Vector3.UnitY, rotationAngle);
                        //Console.WriteLine($"Circle offset before radius subtraction: {vector}");
                        vector = new Vector3(vector.X, vector.Y, (float)(vector.Z - radius));
                        this.Rotation = -rotationAngle;
                        this.RelativeEndpoint = vector;
                    }
                    else
                    {
                        radius = Math.Abs(radius);
                        vector = RotateAboutAxisRad(new Vector3(0, 0, (float)-radius), Vector3.UnitY, -rotationAngle);
                        //Console.WriteLine($"Circle offset before radius addition: {vector}");
                        vector = new Vector3(vector.X, vector.Y, (float)(vector.Z + radius));
                        this.Rotation = rotationAngle;
                        this.RelativeEndpoint = vector;
                        //Console.WriteLine($"Circle Offset: {Offset} | Rotation: {Rotation}");
                    }
                    //Console.WriteLine($"Circle stats: radius: {radius}, negativeCurvature: {negativeCurvature}, circumference: {circumference}...");
                    //Console.WriteLine($"Circle stats cont: rotationAngle: {rotationAngle}, offset: {vector}");
                    break;
                default:
                    //clothoid
                    //the angle is represented by the difference of squares of the scaled curvature parameters.
                    this.RelativeEndpoint = SampleClothoidLocal(1);
                    //float c = (float)System.Math.Sqrt(System.Math.PI/2);
                    double t1 = StartCurvature * B * CURVE_FACTOR;
                    double t2 = EndCurvature * B * CURVE_FACTOR;
                    //this.Rotation = - TotalArcLength * EndCurvature * 180 / (2 * Mathf.PI);

                    if (t2 > t1) this.Rotation = (t2 * t2) - (t1 * t1);// * 180f / Mathf.PI;
                    else this.Rotation = (t1 * t1) - (t2 * t2);// * 180f / Mathf.PI;
                    break;
            }

            //Console.WriteLine($"New Segment End tangent angle: {this.Rotation}");

            /*if (isMirroredX) {
                this.Rotation += 180;
            }*/
        }

        /// <summary>
        /// Sample a clothoid segment given its arc length, curvature difference and an interpolations value
        /// </summary>
        /// <param name="arcLengthStart"></param>
        /// <param name="arcLengthEnd"></param>
        /// <param name="curvatureStart"></param>
        /// <param name="curvatureEnd"></param>
        /// <param name="interpolation"></param>
        /// <returns></returns>
        public static Vector3 SampleClothoidSegment2(double arcLengthStart, double arcLengthEnd, double curvatureStart, double curvatureEnd, double interpolation)
        {
            double arcLength = arcLengthEnd - arcLengthStart;
            double curveDiff = curvatureEnd - curvatureStart;
            double B = Math.Sqrt(arcLength / (Math.PI * Math.Abs(curveDiff)));
            double t1 = curvatureStart * B * CURVE_FACTOR;
            double t2 = curvatureEnd * B * CURVE_FACTOR;
            double interpolatedT = t1 + (interpolation * (t2 - t1));
            Vector3 output = (float)(B * Math.PI) * SampleClothoidSegment(t1, t2, interpolatedT);
            if (float.IsNaN(output.X) || float.IsNaN(output.Y) || float.IsNaN(output.Z))
            {
                Console.WriteLine("NAN OUTPUT");
                Console.WriteLine($"arclengthstart: {arcLengthStart}");
                Console.WriteLine($"arcLegnthEnd: {arcLengthEnd}");
                Console.WriteLine($"curvatureStart: {curvatureStart}");
                Console.WriteLine($"curvatureEnd: {curvatureEnd}");
                Console.WriteLine($"interpolation: {interpolation}");
                Console.WriteLine($"curveDiff: {curveDiff}");
                Console.WriteLine($"B: {B}");
                Console.WriteLine($"t1: {t1}");
                Console.WriteLine($"t2: {t2}");
                Console.WriteLine($"interpolatedT: {interpolatedT}");
                Console.WriteLine($"output: {output}");
            }
            return output;
        }

        /// <summary>
        /// Sample a clothoid between two curvature values
        /// </summary>
        /// <param name="t1"></param>
        /// <param name="t2"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public static Vector3 SampleClothoidSegment(double t1, double t2, double t)
        {
            Vector3 point = new Vector3((float)Mathc.C(t), 0, (float)Mathc.S(t)) - new Vector3((float)Mathc.C(t1), 0, (float)Mathc.S(t1));
            if (t2 > t1)
            {
                point = RotateAboutAxisRad(point, Vector3.UnitY, t1 * t1);
            }
            else
            {
                point = RotateAboutAxisRad(point, Vector3.UnitY, t1 * t1 + Math.PI);
                point = new Vector3(point.X, point.Y, -point.Z);
            }
            return point / (float)CURVE_FACTOR;
        }

        public static Mathc.VectorDouble SampleClothoidSegmentD(double t1, double t2, double t)
        {
            Mathc.VectorDouble point = new Mathc.VectorDouble(Mathc.C(t), 0, Mathc.S(t)) - new Mathc.VectorDouble(Mathc.C(t1), 0, Mathc.S(t1));
            if (t2 > t1)
            {
                point = Mathc.VectorDouble.RotateAboutAxis(point, Mathc.VectorDouble.UnitY, t1 * t1);
            }
            else
            {
                point = Mathc.VectorDouble.RotateAboutAxis(point, Mathc.VectorDouble.UnitY, t1 * t1 + Math.PI);
                point = new Mathc.VectorDouble(point.X, point.Y, -point.Z);
            }
            return point / CURVE_FACTOR;
        }

        /// <summary>
        /// Sample the clothoid segment in local space, with an interpolation value between 0 and 1. 
        /// </summary>
        /// <param name="interpolation"></param>
        /// <returns></returns>
        public Vector3 SampleClothoidLocal(double interpolation)
        {
            double t1 = StartCurvature * B * CURVE_FACTOR;
            double t2 = EndCurvature * B * CURVE_FACTOR;
            double t = t1 + (interpolation * (t2 - t1));
            Vector3 s = SampleClothoidSegment(t1, t2, t);

            //if (isMirroredX) return B * Mathf.PI * new Vector3(s.x, s.y, -s.z);
            /*else */
            return (float)(B * Math.PI) * s;
        }

        public Mathc.VectorDouble SampleClothoidD(double interpolation)
        {
            double t1 = StartCurvature * B * CURVE_FACTOR;
            double t2 = EndCurvature * B * CURVE_FACTOR;
            double t = t1 + (interpolation * (t2 - t1));
            Mathc.VectorDouble s = SampleClothoidSegmentD(t1, t2, t);
            return B * Math.PI * s;
        }

        /// <summary>
        /// Rotates a vector about an axis by an angle in degrees. If the angle is positive, this results in a clockwise rotation.
        /// If the angle is negative, this results in counterclockwise rotation.
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="axis"></param>
        /// <param name="degrees"></param>
        /// <returns></returns>
        public static Vector3 RotateAboutAxisDeg(Vector3 vector, Vector3 axis, double degrees)
        {
            return RotateAboutAxisRad(vector, axis, degrees * Math.PI / 180);
        }

        public static Vector3 RotateAboutAxisRad(Vector3 vector, Vector3 axis, double radians)
        {
            Quaternion rotation = Quaternion.CreateFromAxisAngle(axis, (float)radians);
            return Vector3.Transform(vector, rotation);
        }

        /// <summary>
        /// Helper function to convert a list of N Vector3(arcLength, 0, curvature) into a list of N-1 ClothoidSegment2s
        /// </summary>
        /// <param name="lkNodes"></param>
        /// <returns></returns>
        public static List<ClothoidSegment2> GenerateSegmentsFromLKGraph(List<Vector3> lkNodes)
        {
            List<ClothoidSegment2> segments = new List<ClothoidSegment2>();
            for (int segmentIndex = 0; segmentIndex + 1 < lkNodes.Count; segmentIndex++)
            {
                float arcLengthStart = lkNodes[segmentIndex].X;
                float arcLengthEnd = lkNodes[segmentIndex + 1].X;
                float startCurvature = lkNodes[segmentIndex].Z;
                float endCurvature = lkNodes[segmentIndex + 1].Z;
                ClothoidSegment2 segment = new ClothoidSegment2(Vector3.Zero, arcLengthStart, arcLengthEnd, startCurvature, endCurvature);
                segments.Add(segment);
                //Console.WriteLine(segment.ToString());
            }
            return segments;
        }

        /// <summary>
        /// The start and end curvature values uniquely determine the type of shape to be drawn, these are enumerated in the LineType object as the straight line segment,
        /// the circular arc segment, and the clothoid segment.
        /// </summary>
        /// <param name="startCurvature"></param>
        /// <param name="endCurvature"></param>
        /// <returns></returns>
        public static LineType GetLineTypeFromCurvatureDiff(double startCurvature, double endCurvature)
        {
            double curveDiff = endCurvature - startCurvature;
            //Console.WriteLine($"endCurvature - startCurvature :=> {endCurvature} - {startCurvature} = {curveDiff}");
            //Console.WriteLine($"MIN_CURVATURE_DIFF: {MIN_CURVATURE_DIFF}");
            if (Math.Abs(curveDiff) >= MIN_CURVATURE_DIFF)
            {
                return LineType.CLOTHOID;
            }
            else if (Math.Abs(curveDiff) < MIN_CURVATURE_DIFF && (Math.Abs(endCurvature) < MIN_CURVATURE_DIFF || Math.Abs(startCurvature) < MIN_CURVATURE_DIFF))
            {
                return LineType.LINE;
            }
            else
            {
                return LineType.CIRCLE;
            }
        }

        public override string ToString()
        {
            return $"LineType {LineType} | ArcLengthStart {ArcLengthStart} | ArcLengthEnd {ArcLengthEnd} | StartCurvature {StartCurvature} | EndCurvature {EndCurvature} | B {B}";
        }

        public string Description()
        {
            string s = "";

            switch (this.LineType)
            {
                case LineType.LINE:
                    s += $"Line with length {TotalArcLength}";
                    break;
                case LineType.CIRCLE:
                    double radius = 2f / (StartCurvature + EndCurvature);
                    string b = radius < 0 ? "negative" : "positive";
                    s += $"Circle with a {b} radius of {Math.Abs(radius)}";
                    break;
                case LineType.CLOTHOID:
                    string c = StartCurvature >= 0 ? "positive" : "negative";
                    string d = EndCurvature >= 0 ? "positive" : "negative";
                    s += $"Clothoid with {c} start curvature and {d} end curvature";
                    break;
            }

            s += $"\nMy total arc length is {TotalArcLength}.";
            return s;
        }
    }
}