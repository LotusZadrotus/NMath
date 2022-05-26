namespace NMath.MathObjects
{
    public class HilbertMtrx : Matrix
    {
        public HilbertMtrx(int n):base(n,n)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    this[i, j] = 1.000 / (i + j + 1.000);
                }
            }
        }

        public override string ToString()
        {
            string ret = "";
            for (int i = 0; i < this.Row; i++)
            {
                for (int j = 0; j < this.Column; j++)
                {
                    ret+='1'+"\\"+(i+j+1).ToString()+' ';
                }

                ret += '\n';
            }
            return ret;
        }
    }
}