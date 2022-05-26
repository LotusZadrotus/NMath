using System;
using System.Globalization;
using NMath.MathObjects;

namespace NMath
{
    public static partial class NMath
    {
        public static class MatrixOperation
        {
            
            public static Matrix TriangleInverse(Matrix inputMatrix)
            {
                if (inputMatrix.IsTriangle("Upper"))
                {
                    return UpInv(inputMatrix);
                }

                return UndInv(inputMatrix);
            }

            private static Matrix UpInv(Matrix a)
            {
                double sum;
                Matrix ret = new Matrix(a);
                for (int i = 0; i < a.Row; i++)
                {
                    for (int j = i; j < a.Row; j++)
                    {
                        if (i==j)
                        {
                            ret[i, j] = 1 / a[i, j];
                        }
                        else
                        {
                            sum = 0;
                            for (int k = 0; k < j; k++)
                            {
                                sum += ret[i, k] * a[k, j];
                            }
                            ret[i, j] = - sum / a[i, i];
                        }
                    }
                }
                return ret;
            }
            private static Matrix UndInv(Matrix a)
            {
                double sum;
                Matrix ret = new Matrix(a);
                for (int i = 0; i < a.Row; i++)
                {
                    ret[i, i] = 1 / a[i, i];
                    for (int j = 0; j < i; j++)
                    {
                        sum = 0;
                        for (int k = 0; k < i; k++)
                        {
                            sum += ret[k, j] * a[i, k];
                        }
                        ret[i, j] = - sum / a[i, i];
                    }
                }
                return ret;
            }
            public static class MtrxDet
            {
               
                private static Matrix Gauss(Matrix inputMatrix, out int s)
                {
                    s = 0;
                    Matrix matrix = new Matrix(inputMatrix);
                    if (!matrix.IsSquare()) throw new Exception();
                    double[] nums = new Double[matrix.Row];
                    for (int i = 0; i < matrix.Row; i++)
                    {
                        int max = i;
                        for (int k = 0; k < matrix.Row; k++)
                        {
                            double a = 0;
                            if (System.Math.Abs(matrix[k, i]) > a)
                            {
                                a = matrix[k, i];
                                max = k;
                            }
                        }

                        if (matrix[i, i] < Math.Abs(matrix[max, i])) s++;
                        for (int k = 0; k < matrix.Column; k++)
                        {
                            var num = matrix[i, k];
                            matrix[i, k] = matrix[max, k];
                            matrix[max, k] = num;

                        }

                        for (int j = 0; j < matrix.Row; j++)
                        {
                            nums[j] = matrix[j, i] / matrix[i, i];
                        }

                        for (int j = i + 1; j < matrix.Row; j++)
                        {
                            for (int k = 0; k < matrix.Column; k++)
                            {
                                matrix[j, k] = matrix[j, k] - matrix[i, k] * nums[j];
                            }
                        }
                    }

                    return matrix;
                }

                public static double GaussDet(Matrix inputMatrix)
                {
                    var matrix = Gauss(inputMatrix, out int s);
                    double det = 1;
                    for (int i = 0; i < matrix.Row; i++)
                    {
                        det *= matrix[i, i];
                    }

                    Console.WriteLine(s);
                    return Math.Pow(-1, s) * det;
                }
            }

            public static class MtrxProduct
            {
                private static (Matrix, Matrix) reshape(Matrix a, Matrix b)
                {
                    int row, column = 0;
                    row = Math.Max(a.Row, b.Row);
                    column = Math.Max(a.Column, b.Column);
                    if (row % 2 != 0) row++;
                    if (column % 2 != 0) column++;
                    Matrix ret1 = new Matrix(row, column);
                    Matrix ret2 = new Matrix(row, column);
                    for (int i = 0; i < a.Row; i++)
                    {
                        for (int j = 0; j < a.Column; j++)
                        {
                            ret1[i, j] = a[i, j];
                        }
                    }

                    for (int i = 0; i < a.Row; i++)
                    {
                        for (int j = 0; j < a.Column; j++)
                        {
                            ret2[i, j] = b[i, j];
                        }
                    }

                    return (ret1, ret2);
                }

                private static (Matrix, Matrix) StrassenReshape(IProductable a, IProductable b, out int r, out int c)
                {
                    if (a.Row == b.Row && a.Column == b.Column && Math.Log(a.Row, 2) % 1 == 0 &&
                        Math.Log(a.Column, 2) % 1 == 0 && Math.Log(a.Row, 2)!=0 && Math.Log(a.Column, 2)!=0)
                    {
                        r = 0;
                        c = 0;
                        return ((Matrix) a, (Matrix) b);
                    }

                    int row, column = 0;
                    row = Math.Max(a.Row, b.Row);
                    column = Math.Max(a.Column, b.Column);
                    row = Convert.ToInt32(Math.Pow(2, Math.Ceiling(Math.Log(row, 2))));
                    column = Convert.ToInt32(Math.Pow(2, Math.Ceiling(Math.Log(column, 2))));
                    r = row - a.Row;
                    c = column - b.Column;
                    Matrix ret1 = new Matrix(row, column);
                    Matrix ret2 = new Matrix(row, column);
                    for (int i = 0; i < a.Row; i++)
                    {
                        for (int j = 0; j < a.Column; j++)
                        {
                            ret1[i, j] = a[i, j];
                        }
                    }

                    for (int i = 0; i < b.Row; i++)
                    {
                        for (int j = 0; j < b.Column; j++)
                        {
                            ret2[i, j] = b[i, j];
                        }
                    }
                    return (ret1, ret2);
                }

                private static (Matrix[,], Matrix[,]) StrassenReform(Matrix a, Matrix b)
                {
                    Matrix[,] stds1 = new Matrix[2, 2];
                    Matrix[,] stds2 = new Matrix[2, 2];
                    stds1[0, 0] = new Matrix(a.Row / 2, a.Column / 2);
                    stds1[0, 1] = new Matrix(a.Row / 2, a.Column / 2);
                    stds1[1, 0] = new Matrix(a.Row / 2, a.Column / 2);
                    stds1[1, 1] = new Matrix(a.Row / 2, a.Column / 2);
                    stds2[0, 0] = new Matrix(a.Row / 2, a.Column / 2);
                    stds2[0, 1] = new Matrix(a.Row / 2, a.Column / 2);
                    stds2[1, 0] = new Matrix(a.Row / 2, a.Column / 2);
                    stds2[1, 1] = new Matrix(a.Row / 2, a.Column / 2);
                    for (int i = 0; i < a.Row / 2; i++)
                    {
                        for (int j = 0; j < a.Column / 2; j++)
                        {
                            stds1[0, 0][i, j] = a[i, j];
                            stds2[0, 0][i, j] = b[i, j];
                        }

                        for (int j = a.Column / 2; j < a.Column; j++)
                        {
                            stds1[0, 1][i, j - a.Column / 2] = a[i, j];
                            stds2[0, 1][i, j - a.Column / 2] = b[i, j];
                        }
                    }

                    for (int i = a.Row / 2; i < a.Row; i++)
                    {
                        for (int j = 0; j < a.Column / 2; j++)
                        {
                            stds1[1, 0][i - a.Column / 2, j] = a[i, j];
                            stds2[1, 0][i - a.Column / 2, j] = b[i, j];
                        }

                        for (int j = a.Column / 2; j < a.Column; j++)
                        {
                            stds1[1, 1][i - a.Column / 2, j - a.Column / 2] = a[i, j];
                            stds2[1, 1][i - a.Column / 2, j - a.Column / 2] = b[i, j];
                        }
                    }

                    return (stds1, stds2);
                }

                public static Matrix StrassenProduct(IProductable a, IProductable b)
                {
                    var values = StrassenReshape(a, b, out int r, out int p);
                    Matrix ou1 = new Matrix(values.Item1.Row, values.Item2.Column);
                    if (values.Item1.Row == 2)
                    {
                        var p11 = (values.Item1[0, 0] + values.Item1[1, 1]) * (values.Item2[0, 0] + values.Item2[1, 1]);
                        var p12 = (values.Item1[1, 0] + values.Item1[1, 1]) * values.Item2[0, 0];
                        var p13 = (values.Item1[0, 0]) * (values.Item2[0, 1] - values.Item2[1, 1]);
                        var p14 = (values.Item1[1, 1]) * (values.Item2[1, 0] - values.Item2[0, 0]);
                        var p15 = (values.Item1[0, 0] + values.Item1[0, 1]) * values.Item2[1, 1];
                        var p16 = (values.Item1[1, 0] - values.Item1[0, 0]) * (values.Item2[0, 0] + values.Item2[0, 1]);
                        var p17 = (values.Item1[0, 1] - values.Item1[1, 1]) * (values.Item2[1, 0] + values.Item2[1, 1]);
                        Matrix c = new Matrix(2, 2);
                        c[0, 0] = p11 + p14 - p15 + p17;
                        c[0, 1] = p13 + p15;
                        c[1, 0] = p12 + p14;
                        c[1, 1] = p11 - p12 + p13 + p16;
                        Matrix k = new Matrix(c.Row-r,c.Column-p);
                        for (int i = 0; i < k.Row; i++)
                        {
                            for (int j = 0; j < k.Column; j++)
                            {
                                k[i, j] = c[i, j];
                            }
                        }
                        return k;
                    }

                    var values1 = StrassenReform(values.Item1, values.Item2);
                    var p1 = (values1.Item1[0, 0] + values1.Item1[1, 1]) * (values1.Item2[0, 0] + values1.Item2[1, 1]);
                    var p2 = (values1.Item1[1, 0] + values1.Item1[1, 1]) * values1.Item2[0, 0];
                    var p3 = (values1.Item1[0, 0]) * (values1.Item2[0, 1] - values1.Item2[1, 1]);
                    var p4 = (values1.Item1[1, 1]) * (values1.Item2[1, 0] - values1.Item2[0, 0]);
                    var p5 = (values1.Item1[0, 0] + values1.Item1[0, 1]) * values1.Item2[1, 1];
                    var p6 = (values1.Item1[1, 0] - values1.Item1[0, 0]) * (values1.Item2[0, 0] + values1.Item2[0, 1]);
                    var p7 = (values1.Item1[0, 1] - values1.Item1[1, 1]) * (values1.Item2[1, 0] + values1.Item2[1, 1]);
                    Matrix[,] c1 = new Matrix[2, 2];
                    c1[0, 0] = p1 + p4 - p5 + p7;
                    c1[0, 1] = p3 + p5;
                    c1[1, 0] = p2 + p4;
                    c1[1, 1] = p1 - p2 + p3 + p6;
                    for (int i = 0; i < ou1.Row / 2; i++)
                    {
                        for (int j = 0; j < ou1.Column / 2; j++)
                        {
                            ou1[i, j] = c1[0, 0][i, j];
                        }

                        for (int j = ou1.Column / 2; j < ou1.Column; j++)
                        {
                            ou1[i, j] = c1[0, 1][i, j - ou1.Column / 2];
                        }
                    }

                    for (int i = ou1.Row / 2; i < ou1.Row; i++)
                    {
                        for (int j = 0; j < ou1.Column / 2; j++)
                        {
                            ou1[i, j] = c1[1, 0][i - ou1.Column / 2, j];
                        }

                        for (int j = ou1.Column / 2; j < ou1.Column; j++)
                        {
                            ou1[i, j] = c1[1, 1][i - ou1.Column / 2, j - ou1.Column / 2];
                            //Console.WriteLine($"{i},{j}");
                        }
                    }

                    /*IProductable hoh;
                    if (/*ou1.Row - r == 0 || ou1.Column - p == 1)
                    {
                        hoh = new Vector(ou1.Column - p);
                        for (int i = 0; i < ou1.Column - p; i++)
                        {
                            hoh[i, 0] = ou1[i, 0];
                        }
                    }
                    if (ou1.Row - r == 1)
                    {
                        hoh = new Vector(ou1.Row - r);
                        for (int i = 0; i < ou1.Row - r; i++)
                        {
                            hoh[i, 0] = ou1[0, i];
                        }
                    }*/
                    Matrix iuk = new Matrix(ou1.Row - r, ou1.Column - p);
                    for (int i = 0; i < iuk.Row; i++)
                    {
                        for (int j = 0; j < iuk.Column; j++)
                        {
                            iuk[i, j] = ou1[i, j];
                        }
                    }

                    return iuk;
                }

                /* private static void reform(ref Matrix mtr, int a, int b/*,out int CToDelete,out int RToDelete*/
                /*{
                    
                    // CToDelete = 0;
                    // RToDelete = 0;
                    for (int i = mtr.Row; i < a; i++)
                    {
                        mtr.NewRow();
                        // RToDelete += 1;
                    }
                    for (int i = mtr.Column; i < b; i++)
                    {
                        mtr.NewColumn();
                        // CToDelete += 1;
                    }

                }

                public static Matrix. Stdr(Matrix. a, Matrix. b)
                {
                    if (a.Column != b.Row) throw new Exception();
                    int ok = System.Math.Max(System.Math.Max(a.Column, a.Row), System.Math.Max(b.Column, b.Row));

                    double sum;
                    var mtr0 = new Matrix.(a);
                    var mtr1 = new Matrix.(b);
                    reform(ref mtr0, ok, ok /*out columnToDelete);
                    reform(ref mtr1, ok, ok/*, out rowToDelete);

                    var matrx = new Matrix.(mtr0.Row, mtr1.Column);

                    try
                    {
                        for (int i = 0; i < mtr0.Row; i++)
                        {
                            for (int j = 0; j < mtr0.Column; j++)
                            {
                                sum = 0;
                                for (int k = 0; k < mtr1.Column; k++)
                                {
                                    sum += mtr0[i, k] * mtr1[k, j];
                                }

                                matrx[i, j] = sum;
                            }
                        }
                    }
                    catch (Exception e)
                    {
                        // ignored
                    }

                    for (int i = matrx.Column; i > b.Column; i--)
                    {
                        matrx.RemoveColumn();
                    }
                    for (int i = matrx.Row; i > a.Row; i--)
                    {
                        matrx.RemoveRow();
                    }
                    return matrx;
                }

                
            }*/
                public static Matrix Transpose(Matrix matrix)
                {

                    var matrix1 = new Matrix(matrix.Column, matrix.Row);
                    for (int i = 0; i < matrix.Row; i++)
                    {
                        for (int j = 0; j < matrix.Column; j++)
                        {
                            matrix1[j, i] = matrix[i, j];
                        }
                    }

                    return matrix1;
                }

                public static Matrix Invert(Matrix InputMatrix)
                {
                    var matrix = new Matrix(InputMatrix);
                    var matrix1 = new IdentityMatrix(matrix.Row);
                    double[] nums1 = new Double[matrix.Row];
                    double[] nums2 = new Double[matrix.Row];
                    for (int i = 0; i < matrix.Row; i++)
                    {
                        if (matrix[i, i] == 0)
                        {
                            if (i == matrix.Row) break;
                            for (int j = i + 1; j < matrix.Row; j++)
                            {
                                if (matrix[j, i] != 0)
                                {
                                    for (int k = 0; k < matrix.Column; k++)
                                    {
                                        var num = matrix[j, k];
                                        matrix[j, k] = matrix[i, k];
                                        matrix[i, k] = num;
                                    }
                                }
                            }
                        }
                        var a = matrix[i, i];
                        for (int j = matrix.Row - 1; j >= 0; j--)
                        {
                            matrix[i, j] = matrix[i, j] / a;
                            matrix1[i, j] = matrix1[i, j] / a;
                        }
                        for (int j = i + 1; j < matrix.Row; j++)
                        {
                            var koef = matrix[j, i] / matrix[i, i];
                            for (int k = 0; k < matrix.Row; k++)
                            {
                                matrix[j, k] = matrix[j, k] - matrix[i, k] * koef;
                                matrix1[j, k] = matrix1[j, k] - matrix1[i, k] * koef;
                            }
                        } 
                    }
                    for (var i = matrix.Row - 1; i >= 0 ; i--)
                    {
                        var a = matrix[i, i];
                        for (int j = matrix.Row - 1; j >= 0; j--)
                        {
                            matrix[i, j] = matrix[i, j] / a;
                            matrix1[i, j] = matrix1[i, j] / a;
                        }
                        for (int j = i-1;j>=0; j--)
                        {
                            var koef = matrix[j, i] / matrix[i, i];
                            for (int k=matrix.Row-1; k >= 0; k--)
                            {
                                
                                matrix[j, k] = matrix[j, k] - matrix[i, k] * koef;
                                matrix1[j, k] = matrix1[j, k] - matrix1[i, k] * koef;
                            }
                        }
                    }
                    return matrix1;
                }
            }

            public static class MtrxEigen
            {
                public static (double,Vector) Eigen_max(Matrix InputMatrix, double e)
                {
                    double numprev, numnow;
                    numnow = 0;
                    numprev = double.NaN;
                    Vector ret = new Vector(InputMatrix.Row);
                    ret.Transpose();
                    Vector ret1;
                    Vector ret2=new Vector(InputMatrix.Row);
                    Random rnd = new Random();
                    for (int i = 0; i < ret.Size; i++)
                        ret[i] = rnd.NextDouble();
                    //normi(ret);
                    
                    for (int i = 0; ; i++)
                    {
                        ret1 = (Vector)(InputMatrix * ret);
                        numnow = ret1 * ret;
                        //Console.WriteLine(ret1);
                        if(Math.Abs(numnow-numprev)<e) break;
                        Console.WriteLine($"Eigen value{i}: {numnow}\nEigen vector{i}: {ret}");
                        Console.WriteLine(ret1.Norm(2));
                        ret1.Transpose();
                        normi(ret1);
                        ret = new Vector(ret1);
                        ret.Transpose();
                        
                        
                        numprev = numnow;
                        numnow = 0;
                    }
                    Console.WriteLine($"Ax= {InputMatrix*ret}\nEigen value * Eigen vector: {numnow*ret}");
                    //Console.WriteLine(ret1.GetType());
                    return (numnow,ret);
                }

                private static void normi(Vector a)
                {
                    var ag = a.Norm(2);
                    for (int i = 0; i < a.Size; i++)
                    {
                        a[i] = a[i] / ag;
                    }
                }
            }
        }
    }
}