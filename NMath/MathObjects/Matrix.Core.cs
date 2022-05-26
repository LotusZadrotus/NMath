using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace NMath.MathObjects
{
    public partial class Matrix : IEnumerable, IEnumerable<double>,IProductable
    {
        private int _row;
        private int _column;
        private List<List<double>> _rowList;
        public int Row => this._row;

        public int Column => this._column;
        
        public Matrix(int row, int column)
        {
            this._rowList = new List<List<double>>();
            this._row = row;
            this._column = column;
            for (int index1 = 0; index1 < row; ++index1)
            {
                List<double> doubleList = new List<double>();
                for (int index2 = 0; index2 < column; ++index2)
                    doubleList.Add(0.0);
                this._rowList.Add(doubleList);
            }
        }

        public Matrix(double[,] a)
        {
            this._column = a.GetLength(1);
            this._row = a.GetLength(0);
            this._rowList = new List<List<double>>();
            for (int index1 = 0; index1 < this._row; ++index1)
            {
                List<double> doubleList = new List<double>();
                for (int index2 = 0; index2 < this._column; ++index2)
                    doubleList.Add(a[index1, index2]);
                this._rowList.Add(doubleList);
            }
        }

        public Matrix(double[] a)
        {
            this._column = a.GetLength(0);
            this._row = a.GetLength(0);
            this._rowList = new List<List<double>>();
            for (int index = 0; index < this._row; ++index)
                this._rowList.Add(new List<double>() { a[index] });
        }
        public static Dictionary<string, Func<IProductable, IProductable, Matrix>> NumFuncs = new Dictionary<string, Func<IProductable, IProductable, Matrix>>()
        {
            {
                "Strassen",
                (Func<IProductable, IProductable, Matrix>) (NMath.MatrixOperation.MtrxProduct.StrassenProduct)
            }
        };

        public Vector this[int index]
        {
            get
            {
                Vector ret = new Vector(this._row);
                for (int i = 0; i < this.Column; i++)
                    ret[i] = this[index, i];
                return ret;
            }
        }
        public Matrix(Matrix a)
        {
            this._row = a.Row;
            this._column = a._column;
            this._rowList = new List<List<double>>();
            for (int index1 = 0; index1 < this._row; ++index1)
            {
                List<double> doubleList = new List<double>();
                for (int index2 = 0; index2 < this._column; ++index2)
                    doubleList.Add(a[index1, index2]);
                this._rowList.Add(doubleList);
            }
        }
        public double this[int index1, int index2]
        {
            get => this._rowList[index1][index2];
            set => this._rowList[index1][index2] = value;
        }
        public override string ToString()
        {
            int num1 = 0;
            string str1 = "";
            int num2 = this.Select<double, int>((Func<double, int>) (a => a.ToString().Length)).Prepend<int>(0).Max();
            foreach (double num3 in (IEnumerable<double>) this)
            {
                string str2 = num3.ToString();
                if (str2.Length != num2)
                {
                    int length = str2.Length;
                    for (int index = 0; index < num2 - length; ++index)
                        str2 = str2.Insert(0, " ");
                }
                string str3 = str2 + " ";
                ++num1;
                if (num1 % this._column == 0)
                    str3 += "\n";
                str1 += str3;
            }
            return str1;
        }

        IEnumerator IEnumerable.GetEnumerator() => ((IEnumerable) this).GetEnumerator();

        IEnumerator<double> IEnumerable<double>.GetEnumerator()
        {
            foreach (List<double> row in this._rowList)
            {
                foreach (double num in row)
                    yield return num;
            }
        }

        public override bool Equals(object obj)
        {
            if (!(obj is Matrix)) return false;
            Matrix matrix = (Matrix)obj;
            if (matrix.Row != this.Row || matrix.Column != this.Column) return false;
            for (int index1 = 0; index1 < this._row; ++index1)
            {
                for (int index2 = 0; index2 < this._column; ++index2)
                {
                    if (System.Math.Abs(matrix[index1, index2] - this[index1, index2]) > 0.0001)
                        return false;
                }
            }
            return true;
        }
        public bool IsSquare() => this._column == this._row;

        public bool IsTriangle(string type=null)
        {
            switch (type)
            {
                case "Upper":
                    return isUpperTriangle();
                case "Under":
                    return isUnderTriangle();
                default: return isUnderTriangle() || isUpperTriangle();
            }
        }

        private bool isUnderTriangle()
        {
            for (int i = 0; i < this._row; i++)
            {
                for (int j = i+1; j < this._column; j++)
                {
                    if (this[i, j] != 0) return false;
                }
            }
            return true;
        }
        private bool isUpperTriangle()
        {
            for (int i = 1; i < this._row; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    // Console.WriteLine(this[i,j]);
                    if (this[i, j] != 0) return false;
                }
            }
            return true;
        }
    }
}