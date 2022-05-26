namespace NMath.MathObjects
{
    public class IdentityMatrix: Matrix
    {
        public IdentityMatrix(int row) : base(row, row)
        {
            for (int i = 0; i < row; i++)
            {
                this[i, i] = 1;
            }
        }
    }
}