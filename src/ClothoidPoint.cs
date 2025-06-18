using System.Numerics;

namespace ClothoidX
{
    /// <summary>
    /// This class encodes information about each point on the input polyline to a clothoid curve. It isn't necessary for most of the solutions but we use it for storing information
    /// for the StandardFrame class and Bertolazzi solution methods.
    /// </summary>
    public class HermiteData
    {
        public double x;
        public double z;
        /// <summary>
        /// The tangent angle in radians from the positive X-axis.
        /// </summary>
        public double tangentAngle;
        public double curvature;

        public Vector3 Position => new Vector3((float)x, 0, (float)z);


        public static HermiteData[] FromPostures(Posture[] postures)
        {
            HermiteData[] data = new HermiteData[postures.Length];

            for (int i = 0; i < data.Length; i++)
            {
                HermiteData d = new HermiteData()
                {
                    curvature = postures[i].Curvature,
                    x = postures[i].X,
                    z = postures[i].Z,
                    tangentAngle = postures[i].Angle
                };

                data[i] = d;
            }

            return data;
        }
    }

}