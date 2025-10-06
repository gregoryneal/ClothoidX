using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace ClothoidX
{
    /// <summary>
    /// Methods described by Singh and McCrae
    /// </summary>
    public static class ClothoidSolutionSinghMcCrae
    {
        private const double CURVE_FACTOR = 1.2533141373155; // sqrt(pi / 2)
        /// <summary>
        /// Approximate an input polyine with a ClothoidCurve composed of segments with G2 continuity.
        /// Uses the method described by Singh McCrae in Computers & Graphics 33 (2009) 452–461:
        /// 
        /// Sketching Piecewise Clothoid Curves
        /// DOI: 10.1016/j.cag.2009.05.006
        /// 
        /// Author: Greg Pritchard, April 2025
        /// 
        /// The generated curve will not pass through every input polyline nodes, however they will usually be very close. This solution
        /// would be good for approximating curves for drawing applications, where exact node values don't necessarily need to be followed.
        /// 
        /// Some planned features that didn't make it over from the Singh McCrae paper:
        /// - G1 (tangent) discontinuities: The LK graph can be used to detect G1 discontinuities as large spikes. We can filter those out and generate
        ///   curves with sharp corners, given some additional processing.
        /// - G3 fitting: The paper describes a method to smooth the LK graph in between segments to generate G3 segment continuity.
        /// 
        /// Some problems that need addressed: 
        /// - Some curves generate anomalous curvature values, this is visible when viewing LK node graph and its segmented counterpart. No work has been
        ///   done to address this yet.
        /// - This method might be better suited for batching the spline generation. With a large number of nodes, the method becomes unstable.
        ///   For instance, split the polyline into segments with some maximum arc length, and join all of those curves together.
        /// - As of 9/21/2025 these methods are broken and need to be reworked. This is on my immediate TODO list.
        /// </summary>
        public static ClothoidCurve G2Spline(List<Vector3> points, bool loop = false)
        {
            Console.WriteLine("Solving G2 Spline for points:");
            for (int i = 0; i < points.Count; i++)
            {
                Console.WriteLine(points[i]);
            }

            Console.WriteLine($"Loop: {loop}");

            if (points.Count < 1) return new ClothoidCurve();
            if (loop) points.Add(points[0]);
            if (points.Count < 3) return new ClothoidCurve();
            //this maps polyline nodes to their respective arc length curvature point (you can also use postures here but this is how SinghMcCrae do it)
            List<Vector3> lkNodes = new List<Vector3>();
            Vector3 lkNode;

            //for each polyline node we need to calculate its arc length and curvature, we store those in the xz values of a Vector3
            for (int i = 0; i < points.Count; i++)
            {
                lkNode = GetLKPositionAtPolylineIndex(lkNodes, points, i);
                if (float.IsNaN(lkNode.X) || float.IsNaN(lkNode.Z)) continue;
                lkNodes.Add(lkNode);
            }
            //run the segmenting algorithm given by Singh & McCrae to regress the LK graph into approximate linear segments.
            List<Vector3> segmentedNodes = SegmentedRegression3(lkNodes, .01f);
            for (int i = 0; i < segmentedNodes.Count; i++)
            {
                Console.WriteLine($"(arc length, 0, curvature) segment node {i}: {segmentedNodes[i]}");
            }
            //build the initial curve from the segments, this is in local space so it doesn't line up with our polyline yet.
            ClothoidCurve c = FromLKGraph(segmentedNodes, points);

            //sample the actual curve at the estimated node arc lengths
            List<Vector3> curveSamples = GetCurveSamples(c, lkNodes);

            for (int i = 0; i < curveSamples.Count; i++)
            {
                Console.WriteLine($"(x, y, z) curve sample {i}: {curveSamples[i]}");
            }
            //center of mass of the polyline and curve as they currently sit
            List<Vector3> centerOfMass = GetCenterOfMass(points, curveSamples, 1);

            //Vector3 startPosition = polyline[0];


            //get the rotation matrix required to align the two curves
            double[][] rotationMatrix = GetRotationMatrix(points, curveSamples, centerOfMass[0], centerOfMass[1]);
            //use the rotation matrix to find the rotation angle
            //double angle = Math.Acos(rotationMatrix[0][0]);
            c.AddBestFitTranslationRotation(centerOfMass[0].ToVD(), centerOfMass[1].ToVD(), rotationMatrix);
            return c;
        }

        /// <summary>
        /// Helper function to convert a list of N Vector3(arcLength, 0, curvature) into a ClothoidCurve with N-1 ClothoidSegments.
        /// </summary>
        /// <param name="lkNodes"></param>
        /// <returns></returns>
        public static ClothoidCurve FromLKGraph(List<Vector3> lkNodes, List<Vector3> polyline)
        {
            ClothoidSegment[] segments = new ClothoidSegment[lkNodes.Count - 1];
            Vector3 offset = Vector3.Zero;
            double angle = 0;
            for (int segmentIndex = 0; segmentIndex + 1 < lkNodes.Count; segmentIndex++)
            {
                segments[segmentIndex] = new ClothoidSegment(
                    offset,
                    angle,
                    lkNodes[segmentIndex].X,
                    lkNodes[segmentIndex + 1].X,
                    lkNodes[segmentIndex].Z,
                    lkNodes[segmentIndex + 1].Z,
                    SolutionType.SINGHMCCRAE);
                //Console.WriteLine(segment.ToString());
                offset += segments[segmentIndex].RelativeEndpoint;
                angle = segments[segmentIndex].AngleEnd;
            }
            ClothoidCurve c = new ClothoidCurve(polyline);
            return c.AddSegments(segments);
        }

        /// <summary>
        /// Rotate a point by a rotation matrix.
        /// </summary>
        /// <param name="unrotatedPoint"></param>
        /// <param name="rotationMatrix"></param>
        /// <returns></returns>
        private static Vector3 GetRotatedPoint(Vector3 unrotatedPoint, double[][] rotationMatrix)
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
        /// Generate samples of the curve in world space.
        /// </summary>
        /// <param name="c"></param>
        /// <param name="cmCurve"></param>
        /// <param name="cmPolyline"></param>
        /// <param name="rotationMatrix"></param>
        /// <param name="numPoints"></param>
        /// <returns></returns>
        public static List<Vector3> Eval(ClothoidCurve c, Vector3 cmCurve, Vector3 cmPolyline, double[][] rotationMatrix, int numPoints)
        {
            List<Vector3> samples = new List<Vector3>(numPoints);
            double increment = c.TotalArcLength / (numPoints - 1);
            for (double arcLength = 0; arcLength < c.TotalArcLength; arcLength += increment)
            {
                for (int i = 0; i < c.Count; i++)
                {
                    if (arcLength <= c[i].ArcLengthEnd && arcLength >= c[i].ArcLengthStart)
                    {
                        //sample the curve at the given arc length
                        samples.Add(EvalAtArcLength(arcLength, c[i], cmCurve, cmPolyline, rotationMatrix));
                        break;
                    }
                }
            }
            samples.Add(EvalAtArcLength(c.TotalArcLength, c[^1], cmCurve, cmPolyline, rotationMatrix)); //ensure the last point is included
            return samples;
        }

        public static List<Vector3> Eval(ClothoidSegment c, Vector3 cmCurve, Vector3 cmPolyline, double[][] rotationMatrix, int numPoints)
        {
            List<Vector3> samples = new List<Vector3>(numPoints);
            double increment = c.TotalArcLength / (numPoints - 1);
            for (double arcLength = 0; arcLength < c.TotalArcLength; arcLength += increment)
            {
                samples.Add(EvalAtArcLength(arcLength, c, cmCurve, cmPolyline, rotationMatrix));
            }
            samples.Add(EvalAtArcLength(c.TotalArcLength, c, cmCurve, cmPolyline, rotationMatrix)); //ensure the last point is included
            return samples;
        }

        /// <summary>
        /// Evaluate a SinghMcCrae generated ClothoidCurve at a given arc length.
        /// </summary>
        /// <param name="arcLength"></param>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Vector3 EvalAtArcLength(double arcLength, ClothoidSegment c, Vector3 cmCurve, Vector3 cmPolyline, double[][] rotationMatrix)
        {
            Console.WriteLine($"ROTMAT: ");
            for (int i = 0; i < rotationMatrix.Length; i++)
            {
                Console.WriteLine();
                Console.Write("[");
                for (int j = 0; j < rotationMatrix[i].Length; j++)
                {
                    Console.Write($"{rotationMatrix[i][j]}");
                    if (j < rotationMatrix[i].Length - 1) Console.Write(",");
                }
                Console.Write("]");
            }

            Console.WriteLine($"CMCURVE: {cmCurve} | CMPOYLINE: {cmPolyline}");
            Vector3 sample = SampleSegmentLocalSpace(arcLength, c);
            //apply the rotation and translation to the sample
            sample -= cmCurve;
            sample = GetRotatedPoint(sample, rotationMatrix);
            sample += cmPolyline;
            return sample;
        }

        private static Mathc.VectorDouble SampleClothoidLocal(double t1, double t2, double t)
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

        private static Mathc.VectorDouble SampleClothoidLocal(double interpolation, ClothoidSegment segment)
        {
            double t1 = segment.StartCurvature * segment.B * CURVE_FACTOR;
            double t2 = segment.EndCurvature * segment.B * CURVE_FACTOR;
            double t = t1 + (interpolation * (t2 - t1));
            Mathc.VectorDouble s = SampleClothoidLocal(t1, t2, t);
            return segment.B * Math.PI * s;
        }

        public static Mathc.VectorDouble SampleSegmentLocalSpace(double arcLength, ClothoidSegment segment)
        {
            //if (arcLength > segment.ArcLengthEnd || arcLength < segment.ArcLengthStart) throw new ArgumentOutOfRangeException();
            Mathc.VectorDouble v;
            double interp = (arcLength - segment.ArcLengthStart) / segment.TotalArcLength;
            switch (segment.LineType)
            {
                case LineType.LINE:
                    v = new Mathc.VectorDouble(arcLength - segment.ArcLengthStart, 0, 0);
                    break;
                case LineType.CIRCLE:
                    double radius = -2f / (segment.StartCurvature + segment.EndCurvature); //this might be positive or negative
                    double circumference = 2 * Math.PI * radius;
                    //double fullSweepAngle_deg = 360.0 * TotalArcLength / circumference;
                    //double rotationAngle = interp * fullSweepAngle_deg; //same here, value in radians
                    double rot = interp * segment.TotalArcLength / radius;
                    Mathc.VectorDouble vector = Mathc.VectorDouble.RotateAboutAxis(new Mathc.VectorDouble(0, 0, radius), Mathc.VectorDouble.UnitY, rot);
                    v = new Mathc.VectorDouble(vector.X, vector.Y, vector.Z - radius);
                    break;
                default:
                    v = SampleClothoidLocal(interp, segment);
                    break;
            }

            v = Mathc.VectorDouble.RotateAboutAxis(v, Mathc.VectorDouble.UnitY, -segment.AngleStart);
            v += new Mathc.VectorDouble(segment.Position);
            return v;
        }

        public static Mathc.VectorDouble SampleSegmentByRelativeLengthD(double value, ClothoidSegment segment)
        {
            if (value < 0) value = 0;
            if (value > 1) value = 1;
            double totalArcLength = segment.ArcLengthStart + (value * (segment.ArcLengthEnd - segment.ArcLengthStart));
            return SampleSegmentLocalSpace(totalArcLength, segment);
        }

        /// <summary>
        /// Uses singular value decomposition to factorize the product of the centered input polyline and the centered curve as a product of a rotation, scaling, and
        /// another rotation. Since the two curves are very close in size due to the generation algorithm, we only need the rotation portion of the factorization.
        /// </summary>
        private static double[][] GetRotationMatrix(List<Vector3> polyline, List<Vector3> curveSamples, Vector3 curveCM, Vector3 polylineCM)
        {
            //1. Translate both sets of points so they are centered around the origin (point - cmPolyline)
            //2. Create a matrix for each set of points. A is the points on the arc length samples and B is the points on the input polyline.
            //2b. The matrices A and B are organized with row vectors, where each row is a point on the line. here is an example with 3 points
            // | ax ay az |
            // | bx by bz |
            // | cx cy cz |
            //2c. Find the transpose of A (A') and multiply it with B => A'B = M
            //3. Find the SVD of M => SVD(M) = USV' where U and V' are orthogonal (real square) matrices and S is a diagonal. V' is the transpose of V.
            //4. The optimal rotation is given by R = UV'.
            //5. To rotate a point by this amount we must multiply the *column* vector by the rotation vector P' = RP where P => [x, y, z]^T, we must read the result as a column vector as well.

            //double[][] rotationMatrix = new double[][] { new double[] { 0 }, new double[] { 0 }, new double[] { 0 } };
            double[][] rotationMatrix = Mathc.SVDJacobiProgram.MatMake(3, 3);

            //1. & 2.
            double[][] cip = new double[polyline.Count][];
            double[][] ccp = new double[polyline.Count][];
            for (int i = 0; i < curveSamples.Count; i++)
            {
                //Debug.Log("Accessing " + i);
                //Vector3 centeredInputPoint = clippedPolyline[i] - this.cmPolyline;
                Vector3 centeredInputPoint = polyline[i] - polylineCM;
                Vector3 centeredCurvePoint = curveSamples[i] - curveCM;
                cip[i] = new double[] { centeredInputPoint.X, centeredInputPoint.Y, centeredInputPoint.Z };
                ccp[i] = new double[] { centeredCurvePoint.X, centeredCurvePoint.Y, centeredCurvePoint.Z };
            }

            //3. 
            double[][] M = Mathc.SVDJacobiProgram.MatProduct(Mathc.SVDJacobiProgram.MatTranspose(ccp), cip);

            //4. Vh is V transpose.
            if (Mathc.SVDJacobiProgram.SVD_Jacobi(M, out double[][] U, out double[][] Vh, out double[] S))
            {

                //5.
                rotationMatrix = Mathc.SVDJacobiProgram.MatProduct(U, Vh);
            }
            else
            {
                Console.WriteLine("Error: No rotation matrix found for spline");
            }

            return rotationMatrix;
        }

        /// <summary>
        /// This calculates the center of mass of the curve and polyline by sampling the arc length points
        /// The endpoint weight is configurable.
        /// </summary>
        /// <param name="endpointWeight"></param>
        /// <returns>Return a new List<Vector3>() { curveCM, polylineCM }</returns>
        private static List<Vector3> GetCenterOfMass(List<Vector3> polyline, List<Vector3> curveSamples, float endpointWeight)
        {
            //assign weights to all points in the input polyline
            float weight = 1f;
            Vector3 polylineCM = Vector3.Zero;
            Vector3 curveCM = Vector3.Zero;
            float totalPolylineWeight = (weight * (polyline.Count - 2)) + (endpointWeight * 2);
            float totalCurveWeight = (weight * (curveSamples.Count - 2)) + (endpointWeight * 2);
            float w;

            for (int i = 0; i < polyline.Count; i++)
            {
                w = i == 0 || i == polyline.Count - 1 ? endpointWeight : weight;
                polylineCM += polyline[i] * w;
            }

            polylineCM /= totalPolylineWeight;

            for (int i = 0; i < curveSamples.Count; i++)
            {
                w = i == 0 || i == curveSamples.Count - 1 ? endpointWeight : weight;
                curveCM += curveSamples[i] * w;
            }

            curveCM /= totalCurveWeight;

            return new List<Vector3>() { curveCM, polylineCM };
        }

        /// <summary>
        /// Get a sample of the curve in local space for each node in the input polyline.
        /// </summary>
        private static List<Vector3> GetCurveSamples(ClothoidCurve c, List<Vector3> lkNodes)
        {
            List<Vector3> curveSamples = new List<Vector3>();
            /*
            for (int i = 0; i < lkNodes.Count; i++)
            {
                Vector3 node = lkNodes[i];
                Vector3 sample;
                //hack - for some reason sampling 0 here gives matrix non conformity error. 
                //TODO: fix the above error
                
                if (node.X == 0) sample = c.SampleCurveFromArcLength(node.X);
                else sample = c.SampleCurveFromArcLength(1E-5);
                
                
                curveSamples.Add(sample);
            }
            */

            double l;
            double k;
            for (int i = 0; i < lkNodes.Count; i++)
            {
                l = lkNodes[i].X;
                k = lkNodes[i].Z;
                for (int j = 0; j < c.Count; j++)
                {
                    if (l >= c[j].ArcLengthStart && l <= c[j].ArcLengthEnd)
                    {
                        curveSamples.Add(c[j].SampleClothoidLocal((l - c[j].ArcLengthStart) / (c[j].TotalArcLength)));
                    }
                }
            }
            return curveSamples;
        }

        private static Vector3 GetLKPositionAtPolylineIndex(List<Vector3> lkNodes, List<Vector3> polyline, int polylineNodeIndex)
        {
            if (polylineNodeIndex == 0)
            {
                //approximate what the start curvature should be with a linear approximation from the next two node segment curvatures.
                if (polyline.Count >= 3)
                {
                    Vector3 pointa = GetLKPositionAtPolylineIndex(lkNodes, polyline, 1);
                    Vector3 pointb = GetLKPositionAtPolylineIndex(lkNodes, polyline, 2);
                    float newArcLength = 0;
                    float newCurvature = ((pointb.Z - pointa.Z) * (newArcLength - pointb.X) / (pointb.X - pointa.X)) + pointb.Z;
                    return new Vector3(newArcLength, 0, newCurvature);
                }
                else
                {
                    return new Vector3(0, 0, GetLKPositionAtPolylineIndex(lkNodes, polyline, 1).Z);
                }
            }
            if (polylineNodeIndex == polyline.Count - 1)
            {
                if (polyline.Count >= 3)
                {
                    Vector3 pointa = GetLKPositionAtPolylineIndex(lkNodes, polyline, polyline.Count - 3);
                    Vector3 pointb = GetLKPositionAtPolylineIndex(lkNodes, polyline, polyline.Count - 2);
                    float newArcLength = ClothoidCurve.EstimateArcLength(polyline, polylineNodeIndex);
                    float newCurvature = ((pointb.Z - pointa.Z) * (newArcLength - pointb.X) / (pointb.X - pointa.X)) + pointb.Z;
                    return new Vector3(newArcLength, 0, newCurvature);
                }
                else
                {
                    return new Vector3(ClothoidCurve.EstimateArcLength(polyline, polylineNodeIndex), 0, GetLKPositionAtPolylineIndex(lkNodes, polyline, polyline.Count - 2).Z);
                }
            }
            int n = lkNodes.IndexOf(polyline[polylineNodeIndex]);
            if (n >= 0 && n < lkNodes.Count)
            {
                return lkNodes[n];
            }

            Vector3 point1 = polyline[polylineNodeIndex - 1];
            Vector3 point2 = polyline[polylineNodeIndex];
            Vector3 point3 = polyline[polylineNodeIndex + 1];
            //negatize the curvature because in this solution the convention is that positive curvature is a right turn, opposite to most conventions
            return new Vector3(ClothoidCurve.EstimateArcLength(polyline, polylineNodeIndex), 0, (float)-Mathc.MoretonSequinCurvature(point1.ToVD(), point2.ToVD(), point3.ToVD()));
        }

        /// <summary>
        /// This function uses the walk matrix to figure out where the segments should start and end on the given polyline node. These values are marked in the segment end array.
        /// </summary>
        /// <param name="walkMatrix"></param>
        /// <param name="startIndex"></param>
        /// <param name="endIndex"></param>
        /// <param name="segmentEnd"></param>
        private static void RecurseThroughMatrix(int[,] walkMatrix, int startIndex, int endIndex, bool[] segmentEnd)
        {
            if (startIndex + 1 > endIndex)
            {
                segmentEnd[startIndex] = true;
                segmentEnd[endIndex] = true;
            }

            if (walkMatrix[startIndex, endIndex] == -1)
            {
                segmentEnd[startIndex] = true;
                segmentEnd[endIndex] = true;
            }
            else
            {
                RecurseThroughMatrix(walkMatrix, startIndex, walkMatrix[startIndex, endIndex], segmentEnd);
                RecurseThroughMatrix(walkMatrix, walkMatrix[startIndex, endIndex], endIndex, segmentEnd);
            }
        }

        /// <summary>
        /// This is the error of the linear regression approximating straight line segments in the LK node map (x -> l -> arc length, y -> k -> curvature) space.
        /// It is given in section 4.2 of the SinghMcCrae paper. A and B should be points in LK space. This function will approximate the error of a linear
        /// regression from A to B by summing the vertical error between the discrete points between A and B and the straight line connecting A and B. 
        /// If A is the point immediately preceding B, then the error is zero by definition.
        /// </summary>
        /// <returns></returns>
        private static float GetFitErrors(Vector3 start, Vector3 end, List<Vector3> dataSet, float slopeMatrixValue, float zIntMatrixValue)
        {
            float totalError = 0f;
            int a = dataSet.IndexOf(start);
            int b = dataSet.IndexOf(end);
            for (int i = a; i <= b; i++)
            {
                totalError += Math.Abs(zIntMatrixValue + (slopeMatrixValue * dataSet[i].X) - dataSet[i].Y);
            }
            return totalError;
        }

        private static bool HeightLineFit(Vector3 a, Vector3 b, List<Vector3> dataList, out float aMatrixValue, out float bMatrixValue)
        {
            int indA = dataList.IndexOf(a);
            int indB = dataList.IndexOf(b);
            if (indA < 0 || indB < 0) throw new ArgumentException("Index of A or index of B is less than zero.");
            if (indA > indB) throw new ArgumentException("Index of a is greater than index of b.");

            float sumX = 0;
            float sumZ = 0;
            float sumXX = 0;
            float sumXZ = 0;
            int n = indB - indA + 1;

            for (int i = indA; i <= indB; i++)
            {
                sumX += dataList[i].X;
                sumZ += dataList[i].Z;
                sumXX += dataList[i].X * dataList[i].X;
                sumXZ += dataList[i].X * dataList[i].Z;
            }

            float m = ((n * sumXZ) - (sumX * sumZ)) / ((n * sumXX) - (sumX * sumX));
            float bValue = (sumZ - (m * sumX)) / n;
            aMatrixValue = m;
            bMatrixValue = bValue;
            return true;
        }

        private static List<Vector3> SegmentedRegression3(List<Vector3> dataSet, float allowableError)
        {
            //The list of points in LK space that determine our clothoid segments
            List<Vector3> segmentedPointSet = new List<Vector3>();
            float[,] errorMatrix = new float[dataSet.Count, dataSet.Count];
            int[,] walkMatrix = new int[dataSet.Count, dataSet.Count];
            float[,] slopeMatrix = new float[dataSet.Count, dataSet.Count];
            float[,] zIntMatrix = new float[dataSet.Count, dataSet.Count];

            for (int i = 0; i < dataSet.Count; i++)
            {
                for (int j = 0; j < dataSet.Count; j++)
                {
                    errorMatrix[i, j] = 0f;
                    walkMatrix[i, j] = -1;
                }
            }

            float segmentCost = allowableError;
            for (int i = 0; i + 1 < dataSet.Count; i++)
            {
                errorMatrix[i, i + 1] = segmentCost;
                if (HeightLineFit(dataSet[i], dataSet[i + 1], dataSet, out float slope, out float zIntercept))
                {
                    slopeMatrix[i, i + 1] = slope;
                    zIntMatrix[i, i + 1] = zIntercept;
                }
            }

            //populate the first diagonal which is the cost of line segment between neighbouring points.
            for (int j = 2; j < dataSet.Count; j++)
            {
                for (int i = 0; i + j < dataSet.Count; i++)
                {
                    // do linear regression on segments between i and i+j inclusive
                    // if that error + penalty for segment is less than 
                    // sum of errors [i][i+j-1] and [i+1][i+j], then use
                    // that
                    // otherwise, use sum of [i][i+j-1] and [i+1][i+j]
                    Vector3 start = dataSet[i];
                    Vector3 end = dataSet[i + j];
                    if (HeightLineFit(start, end, dataSet, out float m, out float b))
                    {
                        slopeMatrix[i, i + j] = m;
                        zIntMatrix[i, i + j] = b;

                        float fitError = GetFitErrors(start, end, dataSet, m, b);
                        float minError = fitError + segmentCost;
                        int minIndex = -1;

                        for (int k = i + 1; k < i + j; k++)
                        {
                            if (errorMatrix[i, k] + errorMatrix[k, i + j] < minError)
                            {
                                minIndex = k;
                                minError = errorMatrix[i, k] + errorMatrix[k, i + j];
                            }
                        }

                        walkMatrix[i, i + j] = minIndex;
                        errorMatrix[i, i + j] = minError;
                    }
                }
            }

            bool[] segmentEnd = new bool[dataSet.Count];
            Array.Fill(segmentEnd, false);

            RecurseThroughMatrix(walkMatrix, 0, dataSet.Count - 1, segmentEnd);

            int endIndex = -1;
            int endNextIndex;

            for (int i = 1; i < dataSet.Count; i++)
            {
                if (segmentEnd[i])
                {
                    endIndex = i;
                    break;
                }
            }

            segmentedPointSet.Add(new Vector3(0, 0, zIntMatrix[0, endIndex]));
            for (int i = 0; i < dataSet.Count - 1; i++)
            {
                if (segmentEnd[i])
                {
                    for (int j = i + 1; j < dataSet.Count; j++)
                    {
                        if (segmentEnd[j])
                        {
                            endIndex = j;
                            break;
                        }
                    }

                    float s1 = slopeMatrix[i, endIndex];
                    float b1 = zIntMatrix[i, endIndex];

                    endNextIndex = endIndex;
                    for (int j = endIndex + 1; j < dataSet.Count; j++)
                    {
                        if (segmentEnd[j])
                        {
                            endNextIndex = j;
                            break;
                        }
                    }

                    if (endNextIndex > endIndex)
                    {
                        float s2 = slopeMatrix[endIndex, endNextIndex];
                        float b2 = zIntMatrix[endIndex, endNextIndex];

                        float xIntersect = (b2 - b1) / (s1 - s2);
                        float zIntersect = (s1 * xIntersect) + b1;

                        segmentedPointSet.Add(new Vector3(xIntersect, 0, zIntersect));
                    }
                    else
                    {
                        segmentedPointSet.Add(new Vector3(dataSet[^1].X, 0, (s1 * dataSet[^1].X) + b1));
                        break;
                    }
                }
            }

            bool foundOvershoot;
            int indexToRemove;
            do
            {
                foundOvershoot = false;
                indexToRemove = -1;
                for (int i = 1; i < segmentedPointSet.Count - 1; i++)
                {
                    if ((segmentedPointSet[i].X > segmentedPointSet[i + 1].X && segmentedPointSet[i].X > segmentedPointSet[i - 1].X) ||
                        (segmentedPointSet[i].X < segmentedPointSet[i - 1].X && segmentedPointSet[i].X < segmentedPointSet[i + 1].X))
                    {
                        indexToRemove = i;
                        foundOvershoot = true;
                        break;
                    }
                }
                if (foundOvershoot && indexToRemove >= 0)
                {
                    segmentedPointSet.RemoveAt(indexToRemove);
                }
            } while (foundOvershoot);

            return segmentedPointSet;
        }

    }
}