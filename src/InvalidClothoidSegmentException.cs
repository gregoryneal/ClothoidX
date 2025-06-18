using System;

namespace ClothoidX
{
    [Serializable]
    public class InvalidClothoidSegmentException : Exception
    {
        public InvalidClothoidSegmentException()
        { }

        public InvalidClothoidSegmentException(string message)
            : base(message)
        { }

        public InvalidClothoidSegmentException(string message, Exception innerException)
            : base(message, innerException)
        { }
    }
}