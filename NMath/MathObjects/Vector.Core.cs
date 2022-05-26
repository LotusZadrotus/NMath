using System;
using System.Collections;
using System.Collections.Generic;

namespace NMath.MathObjects
{
    public class IdentityVector : Vector
    {
        public IdentityVector(int Size): base(Size)
        {
            for (int i = 0; i < Size; i++)
            {
                this[i] = 1;
            }
        }
    }
    public partial class Vector : IEnumerable, IEnumerable<double>,IProductable
    {
        public int Row
        {
            get
            {
                if (!isStroke)
                    return _size;
                else
                    return 1;
            }
        }

        public int Column
        {
            get
            {
                if (isStroke)
                    return _size;
                else
                    return 1;
            }
        }

        private int _size;
        private List<double> _values;
        private Dictionary<int, Func<double>> _norm;
        private bool isStroke=true; //Если вектор является строкой возвращает true
        public int Size => this._size;
        public double Norm(int p)
        {
            this._norm = new Dictionary<int, Func<double>>()
            {
                {
                    1,
                    (Func<double>) (() =>
                    {
                        double num1 = 0.0;
                        foreach (double num2 in (IEnumerable<double>) this)
                            num1 += System.Math.Abs(num2);
                        return num1;
                    })
                },
                {
                    2,
                    (Func<double>) (() =>
                    {
                        double d = 0.0;
                        foreach (double x in this._values)
                            d += System.Math.Pow(x, 2.0);
                        return System.Math.Sqrt(d);
                    })
                },
                {
                    3,
                    (Func<double>) (() =>
                    {
                        double num1 = 0.0;
                        foreach (double num2 in this._values)
                        {
                            if (System.Math.Abs(num2) > num1)
                                num1 = num2;
                        }
                        return num1;
                    })
                }
            };
            return this._norm[p]();
        }
        public Vector(int n)
        {
            this._values = new List<double>();
            this._size = n;
            this._values.Clear();
            for (int index = 0; index < n; ++index)
            {
                this._values.Add(0.0);
            }
        }
        public Vector(Vector vector)
        {
            this._values = new List<double>();
            this._size = vector.Size;
            for (int index = 0; index < vector.Size; ++index)
            {
                this._values.Add(vector[index]);
            }
        }
        public Vector(double[] aInts)
        {
            this._values = new List<double>();
            this._size = aInts.Length;
            this._values.Clear();
            foreach (double aInt in aInts)
                this._values.Add(aInt);
        }
        public double this[int index]
        {
            get => this._values[index];
            set => this._values[index] = value;
        }

        public double this[int index1, int index2]
        {
            get
            {
                if (isStroke && index1 == 0) return this._values[index2];
                else if (isStroke && index1 != 0) throw new Exception("Vector is stroke");
                if (!isStroke && index2 == 0) return this._values[index1];
                else if (!isStroke && index2 != 0) throw new Exception("Vector is column");
                throw new Exception();
            }
            set
            {
                if (isStroke && index1 == 0) this._values[index2] = value;
                else if (isStroke && index1 != 0) throw new Exception("Vector is stroke");
                if (!isStroke && index2 == 0) this._values[index1] = value;
                else if (!isStroke && index2 != 0) throw new Exception("Vector is column");
            }
        }
        public double[] ToDoubles()
        {
            double[] numArray = new double[this._size];
            for (int index = 0; index < this._size; ++index)
                numArray[index] = this[index];
            return numArray;
        }
        public override bool Equals(object obj)
        {
            if ((object) (obj as Vector) == null)
                return false;
            Vector vector = (Vector) obj;
            if (vector.Size != this.Size)
                return false;
            for (int index = 0; index < vector.Size; ++index)
            {
                if (System.Math.Abs(vector[index] - this[index]) > 0.0001)
                    return false;
            }
            return true;
        }
        IEnumerator IEnumerable.GetEnumerator() => ((IEnumerable) this).GetEnumerator();

        IEnumerator<double> IEnumerable<double>.GetEnumerator()
        {
            foreach (double num in this._values)
                yield return num;
        }
        public override string ToString()
        {
            string str = "";
            foreach (double num in this._values)
                str = str + num.ToString() + " ";
            return str;
        }
    }
}