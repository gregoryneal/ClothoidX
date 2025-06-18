using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;

namespace ClothoidX
{

    /// <summary>
    /// A class that holds a sequence of ClothoidSegments and can calculate any position given an arc length (optionally parameterized to the range [0f,1f])
    /// </summary>
    public class ClothoidCurve
    {

        /// <summary>
        /// Overall Offset vector for the curve. If the curve passes through all the polyline nodes, this should be the position of the first node.
        /// </summary>
        public Vector3 Offset = Vector3.Zero;
        /// <summary>
        /// Overall angle offset of the curve in radians. This should be close to if not equal to the first Posture's tangent angle. 
        /// </summary>
        public double AngleOffset = 0;

        /// <summary>
        /// Number of segments in this curve. 
        /// </summary>
        public int Count { get { return segments.Count; } }

        public ClothoidSegment this[int index]
        {
            get => this.segments[index];
        }

        public Vector3 Endpoint => SampleCurveFromArcLength(TotalArcLength);

        public int PolylineCount { get { return this.inputPolyline.Count; } }

        public double StartCurvature
        {
            get
            {
                if (segments.Count == 0) return 0;
                return segments[0].StartCurvature;
            }
        }

        public double EndCurvature
        {
            get
            {
                if (segments.Count == 0) return 0;
                return segments[^1].EndCurvature;
            }
        }

        public double TotalArcLength
        {
            get
            {
                if (segments.Count == 0) return 0;
                return segments[^1].ArcLengthEnd - segments[0].ArcLengthStart;
            }
        }
        protected List<ClothoidSegment> segments = new List<ClothoidSegment>();
        protected List<Vector3> inputPolyline = new List<Vector3>();

        /// <summary>
        /// If the curve position and rotation are approximated, this will be the center of mass of the generated curve in local space.
        /// </summary>
        protected Vector3 curveCM = Vector3.Zero;

        /// <summary>
        /// If the curve position and rotation are approximated, this will be the center of mass of the input polyline nodes.
        /// </summary>
        protected Vector3 polylineCM = Vector3.Zero;

        /// <summary>
        /// If the curve position and rotation are approximated, this will be a 3 column row vector that the curve sample point (3 row column vector)
        /// gets multiplied by  after offsetting the point by the overall curve center of mass.
        /// </summary>
        public double[][] rotationMatrix = Mathc.SVDJacobiProgram.MatMake(1, 3);// new double[3][] { new double[] { 0, 0, 0 } };

        /// <summary>
        /// Create a clothoid curve with just a curve factor.
        /// This curve factor is used to scale down the curve before sampling the fresnel integrals. 
        /// Then the point is scaled back up by the same factor.
        /// As far as I know it is only used in the SinghMcCrae solution.
        /// </summary>
        /// <param name="CURVE_FACTOR"></param>
        public ClothoidCurve()
        {
            this.segments = new List<ClothoidSegment>();
            this.inputPolyline = new List<Vector3>();
        }

        public ClothoidCurve(List<Vector3> inputPolyline)
        {
            this.segments = new List<ClothoidSegment>();
            this.inputPolyline = inputPolyline;
        }

        public ClothoidCurve(List<ClothoidSegment> orderedSegments, List<Vector3> inputPolyline)
        {
            this.segments = orderedSegments;
            this.inputPolyline = inputPolyline;
        }

        /// <summary>
        /// Create a copy of a ClothoidCurve object. 
        /// </summary>
        /// <param name="curve"></param>
        public ClothoidCurve(ClothoidCurve curve) : this(CopySegments(curve), curve.inputPolyline)
        {
            Offset = curve.Offset;
            AngleOffset = curve.AngleOffset;
        }

        /// <summary>
        /// Get the tangent angle from the positive X axis at a specific arc length (in radians).
        /// </summary>
        /// <param name="arcLength"></param>
        /// <returns></returns>
        public double Tangent(double arcLength)
        {
            if (arcLength < 0 || arcLength > TotalArcLength) throw new ArgumentOutOfRangeException("arc length needs to be positive and less than the total arc length of the curve");
            //TODO: this doesn't provide a good value after the first segment, fix ur shit bro
            double angle = AngleOffset;
            for (int i = 0; i < Count; i++)
            {
                if (arcLength >= this[i].ArcLengthStart && arcLength <= this[i].ArcLengthEnd)
                {
                    double relativeAL = arcLength - this[i].ArcLengthStart; //arc length relative to this segment
                    //this works because segments are only represented in the standard frame, with their initial tangent oriented along the +x axis
                    angle += this[i].Tangent(relativeAL);
                    break;
                }
                else
                {
                    angle += this[i].Rotation;
                }
            }
            return angle;
        }

        /// <summary>
        /// Best fit translation and rotation are used when the curve is an approximation of your input polyline. 
        /// This as of now is only used in the SinghMcCrae solution. But if you can't calculate an offset for your
        /// clothoid curve (typically with the offset and tangent of the fist polyline node), you can solve for a best fit
        /// using the SinghMcCrae solution and then running this method. 
        /// 
        /// TODO: See if this can be replaced with the Offset and AngleOffset properties. 
        /// </summary>
        /// <param name="curveCM"></param>
        /// <param name="polylineCM"></param>
        /// <param name="bestRotate"></param>
        public void AddBestFitTranslationRotation(Vector3 curveCM, Vector3 polylineCM, double[][] bestRotate)
        {
            this.curveCM = curveCM;
            this.polylineCM = polylineCM;
            this.rotationMatrix = bestRotate;

            ClothoidX.Mathc.SVDJacobiProgram.MatShow(bestRotate, 2, 4);
        }

        public ClothoidCurve AddSegment(ClothoidSegment newSegment)
        {
            if (segments.Count > 0) newSegment.ShiftStartArcLength(segments[^1].ArcLengthEnd);
            else newSegment.ShiftStartArcLength(0);
            this.segments.Add(newSegment);
            return this;
        }

        public ClothoidCurve AddSegments(params ClothoidSegment[] segments)
        {
            for (int i = 0; i < segments.Length; i++)
            {
                AddSegment(segments[i]);
            }
            return this;
        }

        public void Reset()
        {
            this.segments.Clear();
            this.inputPolyline.Clear();
        }

        /// <summary>
        /// Reflect each segment along its common axis.
        /// </summary>
        public void Reflect()
        {
            for (int i = 0; i < Count; i++)
            {
                this[i].Reflect();
            }

            AngleOffset *= -1;
        }

        /// <summary>
        /// Sample positions from this curve, without applying the input polyline transformations.
        /// In local space of the curve.
        /// </summary>
        /// <param name="arcLength"></param>
        /// <returns></returns>
        public Vector3 SampleCurveFromArcLength(double arcLength)
        {
            //if (arcLength == 0) return Vector3.Zero;
            Vector3 value = Vector3.Zero;
            Vector3 offset = Vector3.Zero;
            double rotation = 0; //radians
            for (int i = 0; i < segments.Count; i++)
            {
                ClothoidSegment segment = segments[i];
                if (arcLength >= segment.ArcLengthStart && arcLength <= segment.ArcLengthEnd)
                {
                    //Debug.Log($"Now sampling curve type: {segment.LineType}");
                    value = segment.SampleSegmentByTotalArcLength(arcLength);
                    value = ClothoidSegment.RotateAboutAxisRad(value, Vector3.UnitY, rotation);
                    //Debug.Log($"Rotate sampled segment by {rotation}");
                    value += offset;
                    break;
                }
                offset += ClothoidSegment.RotateAboutAxisRad(segment.Offset, Vector3.UnitY, rotation);
                rotation -= segment.Rotation; //apply rotation last since rotation is applied around the origin
                //Debug.Log($"New offset and rotation along curve: {offset}, {rotation}");
            }
            value = ClothoidSegment.RotateAboutAxisRad(value, Vector3.UnitY, -AngleOffset);
            //Debug.Log($"Rotate curve by {-AngleOffset}");
            value += Offset;
            return value;
        }

        /// <summary>
        /// Get samples along the curve, if it was generated from a polyline this will return points that approximate the input polyline as well. If not the points will start from the origin.
        /// </summary>
        /// <param name="numSamples"></param>
        /// <returns></returns>
        public List<Vector3> GetSamples(int numSamples)
        {
            List<Vector3> points = new List<Vector3>();
            Vector3 point;
            double increment = TotalArcLength / numSamples;
            //Debug.Log($"Increment size: {increment}, totalArcLength: {TotalArcLength}, numSamples: {numSamples}");
            for (double arcLength = 0; arcLength < TotalArcLength; arcLength += increment)
            {
                point = SampleCurveFromArcLength(arcLength); //untranslated, unrotated
                if (point == Vector3.Zero && arcLength != 0) continue;
                //UnityEngine.Debug.Log($"Curve CM before applying: {curveCM.X}, {curveCM.Y}, {curveCM.Z}");
                point -= curveCM;
                point = GetRotatedPoint(point);
                point += polylineCM;
                points.Add(point);
            }
            point = SampleCurveFromArcLength(TotalArcLength);
            point -= curveCM;
            point = GetRotatedPoint(point);
            point += polylineCM;
            points.Add(point);

            return points;
        }

        /// <summary>
        /// Rotates a point by the rotation matrix. Be sure to treat the unrotatedPoint and rotatedPoint as a column vector.
        /// Be sure to rotate only points that are centered on the origin, and then translate them by the polyline center of mass afterwards.
        /// </summary>
        /// <param name="unrotatedPoint"></param>
        /// <returns></returns>
        public Vector3 GetRotatedPoint(Vector3 unrotatedPoint)
        {
            if (rotationMatrix == null) return unrotatedPoint;
            bool isZero = true;
            for (int i = 0; i < rotationMatrix.Length; i++)
            {
                if (!rotationMatrix[i].All(x => x == 0))
                {
                    isZero = false;
                    break;
                }
            }
            if (isZero)
            {
                return unrotatedPoint;
            }
            double[][] point = Mathc.SVDJacobiProgram.MatMake(3, 1);
            point[0][0] = unrotatedPoint.X;
            point[1][0] = unrotatedPoint.Y;
            point[2][0] = unrotatedPoint.Z;
            double[][] rotatedPoint = Mathc.SVDJacobiProgram.MatProduct(rotationMatrix, point);
            return new Vector3((float)rotatedPoint[0][0], (float)rotatedPoint[1][0], (float)rotatedPoint[2][0]);
        }

        /// <summary>
        /// Create some test segments and display them.
        /// </summary>
        /// <returns></returns>
        public static IEnumerator<List<Vector3>> TestAllConnectionTypes(int a = 0, int b = 10, int c = 15, int d = 20, int e = 35, int f = 49)
        {
            ClothoidSegment lineSeg = new ClothoidSegment(0, b, 0, 0);
            ClothoidSegment lineSeg2 = new ClothoidSegment(b, c, 0, 0);
            ClothoidSegment cplus = new ClothoidSegment(0, b, 1, 1);
            ClothoidSegment cplus2 = new ClothoidSegment(b, c, 1, 1);
            ClothoidSegment cminus = new ClothoidSegment(0, b, -1, -1);
            ClothoidSegment cminus2 = new ClothoidSegment(b, c, -1, -1);
            ClothoidSegment clplus = new ClothoidSegment(0, b, 0, 1);
            ClothoidSegment clplus2 = new ClothoidSegment(b, c, 0, 1);
            ClothoidSegment clminus = new ClothoidSegment(0, b, 0, -1);
            ClothoidSegment clminus2 = new ClothoidSegment(b, c, 0, -1);
            ClothoidSegment clplusminus = new ClothoidSegment(0, b, 1, -1);
            ClothoidSegment clplusminus2 = new ClothoidSegment(b, c, 1, -1);
            ClothoidSegment clminusplus = new ClothoidSegment(0, b, -1, 1);
            ClothoidSegment clminusplus2 = new ClothoidSegment(b, c, -1, 1);

            ClothoidSegment[] segments = new ClothoidSegment[] { lineSeg, cplus, cminus, clplus, clminus, clplusminus, clminusplus };
            ClothoidSegment[] segments2 = new ClothoidSegment[] { lineSeg2, cplus2, cminus2, clplus2, clminus2, clplusminus2, clminusplus2 };

            for (int i = 0; i < segments.Length; i++)
            {
                for (int j = 0; j < segments2.Length; j++)
                {
                    ClothoidCurve curve = new ClothoidCurve(new List<ClothoidSegment>() { segments[i], segments2[j] }, new List<Vector3>());

                    for (int k = 0; k < curve.segments.Count; k++)
                    {
                        //Console.WriteLine($"({i},{j}) Segment{k + 1}: {curve.segments[k]} -> {curve.segments[k].Description()}");
                    }

                    yield return curve.GetSamples(100);
                }
            }
            yield break;
        }

        /// <summary>
        /// Add a random curve using the sharpness constructor
        /// </summary>
        /// <returns></returns>
        public ClothoidCurve AddRandomSegment2()
        {
            double sharpness;
            double startCurvature;
            double newArcLength = (new Random().NextDouble() * 4) + 5;
            double shape = 0.5f;// UnityEngine.Random.value;

            if (segments.Count > 0)
            {
                ClothoidSegment lastSegment = segments[^1];

                if (shape > .66f)
                {
                    //line
                    sharpness = 0;
                    startCurvature = 0;
                }
                else if (shape < .33f)
                {
                    //circle
                    sharpness = 0;
                    startCurvature = lastSegment.EndCurvature;
                }
                else
                {
                    //clothoid
                    sharpness = (new Random().NextDouble() * .06) - .03;
                    startCurvature = lastSegment.EndCurvature;
                }
                ClothoidSegment newSegment = new ClothoidSegment(startCurvature, sharpness, newArcLength);
                return AddSegment(newSegment);
            }
            else
            {
                if (shape > .66f)
                {
                    //line
                    sharpness = 0;
                    startCurvature = 0;
                }
                else if (shape < .33f)
                {
                    //circle
                    sharpness = 0;
                    startCurvature = new Random().NextDouble() - .5; //UnityEngine.Random.Range(-.5f, .5f);
                }
                else
                {
                    //clothoid
                    sharpness = (new Random().NextDouble() * .06) - .03;
                    startCurvature = (new Random().NextDouble() * .6) - .3;
                }
                ClothoidSegment newSegment = new ClothoidSegment(startCurvature, sharpness, newArcLength);
                //Console.WriteLine($"new parameters: sharpness: {sharpness}, arcLength: {newArcLength}, startcurvature: {startCurvature}");
                return AddSegment(newSegment);
            }
        }

        /// <summary>
        /// Add 3 random segments: a line segment followed by a clothoid transition to a circlar arc segment.
        /// </summary>
        /// <returns></returns>
        public ClothoidCurve AddRandomSegment3()
        {
            float lineLength = (new Random().Next() * 5f) + 5f;
            float clothoidLength = (new Random().Next() * 5f) + 5f;
            float arcLength = (new Random().Next() * 5f) + 5f;
            float sharpness = (new Random().Next() * .016f) - .008f;
            float curvature = (new Random().Next() * .16f) - .08f;

            ClothoidSegment lineSegment = new ClothoidSegment(0, 0, lineLength);
            ClothoidSegment clothoidSegment = new ClothoidSegment(0, sharpness, clothoidLength);
            ClothoidSegment circularSegment = new ClothoidSegment(clothoidSegment.EndCurvature, 0, arcLength);
            return AddSegments(lineSegment, clothoidSegment, circularSegment);
        }

        public ClothoidCurve AddRandomSegment4()
        {
            float arcLength = (new Random().Next() * 5f) + 5f;
            float curvature = (new Random().Next() * .16f) - .08f;

            ClothoidSegment circularSegment = new ClothoidSegment(curvature, 0, arcLength);
            return AddSegments(circularSegment);
        }

        public ClothoidCurve AddLine(float length)
        {
            return AddSegment(new ClothoidSegment(0, 0, length));
        }

        public static ClothoidCurve GetRandomCurve()
        {
            ClothoidCurve c = new ClothoidCurve();
            float sharpness;
            float startCurvature;
            float newArcLength = (new Random().Next() * 4) + 5;
            float shape = 0.5f;// UnityEngine.Random.value;
            if (shape > .66f)
            {
                //line
                sharpness = 0;
                startCurvature = 0;
            }
            else if (shape < .33f)
            {
                //circle
                sharpness = 0;
                startCurvature = new Random().Next() - .5f;
            }
            else
            {
                //clothoid
                sharpness = (new Random().Next() * .06f) - .03f;
                startCurvature = (new Random().Next() * .6f) - .3f;
            }
            ClothoidSegment newSegment = new ClothoidSegment(startCurvature, sharpness, newArcLength);
            //Console.WriteLine($"new parameters: sharpness: {sharpness}, arcLength: {newArcLength}, startcurvature: {startCurvature}");
            return c.AddSegment(newSegment);
        }

        /// <summary>
        /// Get the position of a polyline node by its index.
        /// </summary>
        /// <param name="polylineIndex"></param>
        /// <returns></returns>
        public Vector3 GetPositionFromPolyline(int polylineIndex)
        {
            if (polylineIndex >= this.inputPolyline.Count) polylineIndex = this.inputPolyline.Count - 1;
            return this.inputPolyline[polylineIndex];
        }

        /// <summary>
        /// Get the arc length of a polyline node based on its index.
        /// </summary>
        /// <param name="polylineIndex"></param>
        /// <returns></returns>
        public float GetStartingArcLengthFromPolylineIndex(int polylineIndex)
        {
            return EstimateArcLength(polylineIndex);
        }

        protected float EstimateArcLength(Vector3 node)
        {
            if (this.inputPolyline.Contains(node) == false) throw new ArgumentException("Node not contained in the input polyline.");
            return EstimateArcLength(this.inputPolyline.IndexOf(node));
        }

        protected float EstimateArcLength(int polylineIndex)
        {
            if (polylineIndex > this.inputPolyline.Count - 2) polylineIndex = this.inputPolyline.Count - 2;
            if (this.inputPolyline == null) throw new ArgumentNullException();
            if (this.inputPolyline.Count == 0) throw new ArgumentOutOfRangeException();
            float sum = 0;
            for (int i = 0; i + 1 <= polylineIndex; i++)
            {
                sum += Vector3.Distance(this.inputPolyline[i], this.inputPolyline[i + 1]);
            }
            return sum;
        }

        /// <summary>
        /// Estimate the arc length at a given node along the polyline index
        /// </summary>
        /// <param name="polyline"></param>
        /// <param name="index"></param>
        /// <returns></returns>
        public static float EstimateArcLength(List<Vector3> polyline, int index)
        {
            if (index < 0 || index >= polyline.Count) throw new ArgumentOutOfRangeException();
            float sum = 0;
            Vector3 prev = polyline[0];
            for (int i = 1; i <= index; i++)
            {
                sum += Vector3.Distance(prev, polyline[i]);
                prev = polyline[i];
            }
            return sum;
        }


        /// <summary>
        /// Add two clothoid curves together. This operation is not commutative, the order matters. The end curve will be added on to the start curve.
        /// </summary>
        /// <param name="start"></param>
        /// <param name="end"></param>
        /// <returns></returns>
        public static ClothoidCurve operator +(ClothoidCurve start, ClothoidCurve end)
        {
            ClothoidCurve result = new ClothoidCurve();
            //set the offsets the same as the start offset
            result.Offset = start.Offset;
            result.AngleOffset = start.AngleOffset;
            // Add the segments one at a time
            result.AddSegments(start.segments.ToArray())
            .AddSegments(end.segments.ToArray());/*
            for (int i = 0; i < start.segments.Count; i++) {
                result.AddSegment(start.segments[i]);
            }
            for (int i = 0; i < end.segments.Count; i++) {
                result.AddSegment(end.segments[i]);
            }*/

            return result;
        }

        public static ClothoidCurve operator +(ClothoidCurve start, ClothoidSegment newSegment)
        {
            ClothoidCurve result = new ClothoidCurve
            {
                Offset = start.Offset,
                AngleOffset = start.AngleOffset
            };
            result.AddSegments(start.segments.ToArray()).AddSegment(newSegment);
            return result;
        }

        /// <summary>
        /// Create a clothoid curve from a list of segments
        /// </summary>
        /// <param name="segments"></param>
        /// <returns></returns>
        public static ClothoidCurve FromSegments(params ClothoidSegment[] segments)
        {
            return new ClothoidCurve().AddSegments(segments);
        }

        /// <summary>
        /// Helper function to convert a list of N Vector3(arcLength, 0, curvature) into a ClothoidCurve with N-1 ClothoidSegments
        /// </summary>
        /// <param name="lkNodes"></param>
        /// <returns></returns>
        public static ClothoidCurve FromLKGraph(List<Vector3> lkNodes, List<Vector3> polyline)
        {
            ClothoidSegment[] segments = new ClothoidSegment[lkNodes.Count - 1];
            for (int segmentIndex = 0; segmentIndex + 1 < lkNodes.Count; segmentIndex++)
            {
                segments[segmentIndex] = new ClothoidSegment(
                    lkNodes[segmentIndex].X,
                    lkNodes[segmentIndex + 1].X,
                    lkNodes[segmentIndex].Z,
                    lkNodes[segmentIndex + 1].Z);
                //Console.WriteLine(segment.ToString());
            }
            ClothoidCurve c = new ClothoidCurve(polyline);
            return c.AddSegments(segments);
        }

        /// <summary>
        /// Copy the segments of a clothoid curve into a list
        /// </summary>
        /// <param name="curve"></param>
        /// <returns></returns>
        public static List<ClothoidSegment> CopySegments(ClothoidCurve curve)
        {
            List<ClothoidSegment> segments = new List<ClothoidSegment>();
            for (int i = 0; i < curve.Count; i++)
            {
                segments.Add(new ClothoidSegment(curve[i].StartCurvature, curve[i].Sharpness, curve[i].TotalArcLength));
            }
            return segments;
        }

        public override string ToString()
        {
            string s = "";
            for (int i = 0; i < segments.Count; i++)
            {
                s += segments[i].Description() + ", ";
            }
            return s;
        }

        /// <summary>
        /// Return a three segment curve in local space. Start point is at the origin and start tangent is along the +X-axis.
        /// </summary>
        /// <param name="sharpness"></param>
        /// <param name="arcLength1"></param>
        /// <param name="arcLength2"></param>
        /// <param name="arcLength3"></param>
        /// <param name="startCurvature"></param>
        /// <param name="endCurvature"></param>
        /// <param name="start"></param>
        /// <param name="end"></param>
        /// <returns></returns>
        public static ClothoidCurve ThreeSegmentsLocal(float startCurvature, float sharpness, float arcLength1, float arcLength2, float arcLength3)
        {
            ClothoidCurve c = new ClothoidCurve();
            ClothoidSegment s1 = new ClothoidSegment(startCurvature, sharpness, arcLength1);
            ClothoidSegment s2 = new ClothoidSegment(s1.EndCurvature, -sharpness, arcLength2);
            ClothoidSegment s3 = new ClothoidSegment(s2.EndCurvature, sharpness, arcLength3);
            return c.AddSegments(s1, s2, s3);
            //return c.AddSegment(s1);
        }
    }
}
