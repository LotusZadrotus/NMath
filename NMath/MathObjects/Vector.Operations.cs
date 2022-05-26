using System;
using System.Runtime.CompilerServices;
using System.Text.RegularExpressions;

namespace NMath.MathObjects
{
    public partial class Vector
    {
        public static Vector operator +(Vector v1, Vector v2)
        {
            if (v1._size != v2._size)
                throw new Exception("Vectors have different size");
            Vector vector = new Vector(v1._size);
            for (int index = 0; index < v1._size; ++index)
                vector[index] = v1[index] + v2[index];
            return vector;
        }

        public static Vector operator -(Vector v1, Vector v2)
        {
            if (v1._size != v2._size)
                throw new Exception("Vectors have different size");
            Vector vector = new Vector(v1._size);
            for (int index = 0; index < v1._size; ++index)
                vector[index] = v1[index] - v2[index];
            return vector;
        }

        public static Vector operator *(double a, Vector v)
        {
            Vector vector = new Vector(v._size);
            for (int index = 0; index < v._size; ++index)
                vector[index] = a * v[index];
            return vector;
        }

        public static Vector operator *(Vector v, double a)
        {
            Vector vector = new Vector(v._size);
            for (int index = 0; index < v._size; ++index)
                vector[index] = a * v[index];
            return vector;
        }

        public static bool operator ==(Vector v1, Vector v2)
        {
            if (v1._size != v2._size)
                return false;
            for (int index = 0; index < v1._size; ++index)
            {
                if (System.Math.Abs(v1[index] - v2[index]) > 0.0001)
                    return false;
            }
            return true;
        }

        public static bool operator !=(Vector v1, Vector v2)
        {
            if (v1._size != v2._size)
                return true;
            for (int index = 0; index < v1._size; ++index)
            {
                if (System.Math.Abs(v1[index] - v2[index]) > 0.0001)
                    return true;
            }
            return false;
        }

        public static Vector operator %(Vector v1, Vector v2)
        {
            if(v1._size != v2._size) throw new Exception("Vector have different sizes");
            Vector ret =new Vector(v1._size);
            for (int i = 0; i < v1._size; i++)
                ret[i] = v1[i] * v2[i];
            return ret;
        }
        public static double operator *(Vector v1, Vector v2)
        {
            double num = 0.0;
            if (v1.Size != v2.Size)
                throw new Exception("Vectors have different sizes");
            for (int index = 0; index < v1.Size; ++index)
                num += v1[index] * v2[index];
            return num;
        }

        public static Vector operator /(Vector v, double a)
        {
            Vector ret = new Vector(v.Size);
            for (int i = 0; i < ret._size; i++)
            {
                ret[i] = v[i] / a;
            }
            return ret;
        }
        public void Add(double a)
        {
            this._values.Add(a);
            this._size++;
        }
        public void Transpose()
        {
            if (isStroke)
            {
                this.isStroke = false;
            }
            else isStroke = true;
        }
        
    }
}