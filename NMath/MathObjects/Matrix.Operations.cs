using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace NMath.MathObjects
{
    public partial class Matrix
    {
        public double Norm1
        {
            get
            {
                double num1 = 0.0;
                for (int index2 = 0; index2 < this._column; ++index2)
                {
                    double num2 = 0.0;
                    for (int index1 = 0; index1 < this._row; ++index1)
                        num2 += System.Math.Abs(this[index1, index2]);
                    if (num1 < num2)
                        num1 = num2;
                }
                return num1;
            }
        }

        public double NormInf
        {
            get
            {
                double num1 = 0.0;
                for (int index1 = 0; index1 < this._column; ++index1)
                {
                    double num2 = 0.0;
                    for (int index2 = 0; index2 < this._row; ++index2)
                        num2 += System.Math.Abs(this[index1, index2]);
                    if (num1 < num2)
                        num1 = num2;
                }
                return num1;
            }
        }
        public static Matrix operator +(Matrix a, Matrix b)
        {
            if (a.Row != b.Row || a.Column != b.Column)
                throw new Exception();
            Matrix matrix = new Matrix(a.Row, a.Column);
            for (int index1 = 0; index1 < a.Row; ++index1)
            {
                for (int index2 = 0; index2 < a.Column; ++index2)
                    matrix[index1, index2] = a[index1, index2] + b[index1, index2];
            }
            return matrix;
        }

        public static Matrix operator -(Matrix a, Matrix b)
        {
            if (a.Row != b.Row || a.Column != b.Column)
                throw new Exception();
            Matrix matrix = new Matrix(a.Row, a.Column);
            for (int index1 = 0; index1 < a.Row; ++index1)
            {
                for (int index2 = 0; index2 < a.Column; ++index2)
                    matrix[index1, index2] = a[index1, index2] - b[index1, index2];
            }
            return matrix;
        }
        public static Matrix operator *(double a, Matrix matrix)
        {
            for (int index1 = 0; index1 < matrix.Row; ++index1)
            {
                for (int index2 = 0; index2 < matrix.Column; ++index2)
                    matrix[index1, index2] *= a;
            }
            return matrix;
        }
        public static Matrix operator *(IProductable a, Matrix b) => Matrix.NumFuncs["Strassen"](b, a);
        public static Matrix operator *(Matrix a, IProductable b) => Matrix.NumFuncs["Strassen"](b, a);
        public static Matrix operator *(Matrix a, Matrix b) => Matrix.NumFuncs["Strassen"](a, b);
        public static explicit operator Vector(Matrix matrix)
        {
            if (matrix.Row == 1)
            {
                Vector ret = new Vector(matrix.Column);
                ret.Transpose();
                for (int i = 0; i < ret.Size; i++)
                {
                    ret[i] = matrix[0, i];
                }
                //Console.WriteLine(matrix);
                //Console.WriteLine(ret);
                return ret;
            }

            if (matrix.Column == 1)
            {
                Vector ret = new Vector(matrix.Row);
                for (int i = 0; i < ret.Size; i++)
                {
                    ret[i] = matrix[i, 0];
                }
                //Console.WriteLine(matrix);
                //Console.WriteLine(ret);
                return ret;
            }
            throw new Exception("Cannot convert 2 dim matrix to vector");
        }

        public Matrix Transpose()
        {
            List<List<double>> list1 = new List<List<double>>();
            for (int i = 0; i < this.Column; i++)
            {
                List<double> list2 = new List<double>();
                for (int j = 0; j < this._row; j++)
                {
                    list2.Add(this[j,i]);
                }
                list1.Add(list2);
            }
            this._rowList = list1;
            return this;
        }
    }
}