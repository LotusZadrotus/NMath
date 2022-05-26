using System.Runtime.CompilerServices;

namespace NMath.MathObjects
{
    public interface IProductable
    {
        int Row { get; }
        int Column { get; }
        double this[int index1,int index2] { get; set; }
        
    }
}