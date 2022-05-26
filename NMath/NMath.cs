using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using NMath.MathObjects;

namespace NMath
{
    public static partial class NMath
    {
        public static class Equation
            {
                public static Vector SolveGauss1(Matrix inputMatrix, Vector inputVector)
                {
                    Matrix matrix = new Matrix(inputMatrix);
                    Vector vector = new Vector(inputVector);
                    Vector solve = new Vector(vector.Size);
                    if (!matrix.IsSquare()) throw new Exception();
                    if (matrix.Row != vector.Size) throw new Exception();
                    double[] nums = new Double[matrix.Row];
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

                            vector[j] = vector[j] - vector[i] * nums[j];
                        }
                    }

                    solve[solve.Size - 1] = vector[vector.Size - 1] / matrix[matrix.Row - 1, matrix.Column - 1];

                    for (int i = solve.Size - 2; i >= 0; i--)
                    {
                        double sum = 0;
                        for (int j = i; j < solve.Size; j++)
                        {
                            sum += solve[j] * matrix[i, j];
                        }

                        solve[i] = (vector[i] - sum) / matrix[i, i];
                    }

                    return solve;
                }

                public static Vector BCG(Matrix inputMatrix, Vector inputVector, Vector x0)
                {
                    int i = 0;
                    Vector r10 = inputVector - (Vector)(inputMatrix*x0);
                    Vector r20 = new Vector(r10);
                    Vector p10 = new Vector(r10);
                    Vector p20 = new Vector(r10);
                    do
                    {
                        Vector r11 = new Vector(r10);
                        Vector r21 = new Vector(r20);
                        double a = (r10 * r20) / ((Vector) (inputMatrix * p10) * p20);
                        x0 = x0 + a * p10;
                        
                        r10 = r11 - a * (Vector) (inputMatrix * p10);
                        r20 = r21 - a * (Vector) (inputMatrix.Transpose() * p20);
                        double b = (r10 * r20) / (r11 * r21);
                        if (b==0) break;
                        p10 = r10 + b * p10;
                        p20 = r20 + b * p20;
                        i++;
                    } while (r10.Norm(2)>.00001);

                    return x0;
                }

                public static Vector SolveGaussWithMainElementsChoosing(Matrix inputMatrix, Vector inputVector)
                {
                    Matrix matrix = new Matrix(inputMatrix);
                    Vector vector = new Vector(inputVector);
                    Vector solve = new Vector(vector.Size);
                    if (!matrix.IsSquare()) throw new Exception();
                    if (matrix.Row != vector.Size) throw new Exception();
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

                        for (int k = 0; k < matrix.Column; k++)
                        {
                            var num = matrix[i, k];
                            matrix[i, k] = matrix[max, k];
                            matrix[max, k] = num;

                        }

                        var num1 = vector[max];
                        vector[max] = vector[i];
                        vector[i] = num1;
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

                            vector[j] = vector[j] - vector[i] * nums[j];
                        }
                    }

                    solve[solve.Size - 1] = vector[vector.Size - 1] / matrix[matrix.Row - 1, matrix.Column - 1];

                    for (int i = solve.Size - 2; i >= 0; i--)
                    {
                        double sum = 0;
                        for (int j = i; j < solve.Size; j++)
                        {
                            sum += solve[j] * matrix[i, j];
                        }

                        solve[i] = (vector[i] - sum) / matrix[i, i];
                    }

                    return solve;
                }

                private static double Iterations(Matrix inputMatrix, Vector inputVector)
                {
                    return System.Math.Floor(System.Math.Log((.0001 * (1 - inputMatrix.NormInf) / inputVector.Norm(3)),
                        inputMatrix.NormInf));
                }

                public static Vector Rel(Matrix a, Vector b, double omega)
                {
                    double sum1;
                    double sum2;
                    Vector k1 = new Vector(b.Size);
                    Vector k2 = new Vector(b.Size);
                    for (int k = 0; k < 100; k++)
                    {
                        for (int i = 0; i < k1.Size; i++)
                        {
                            sum1 = 0;
                            sum2 = 0;
                            for (int j = 0; j < i; j++)
                            {
                                sum1 += a[i, j] * k2[j];
                            }

                            for (int j = i + 1; j < k1.Size; j++)
                            {
                                sum2 += a[i, j] * k1[j];
                            }

                            k2[i] = (b[i] - sum1 - sum2) / a[i, i];
                        }

                        k2 = omega * k2 + (1 - omega) * k1;
                        if ((k2 - k1).Norm(2) < 0.001) break;
                        k1=new Vector(k2);
                    }

                    return k2;
                }
                public static Vector SeidelIteration(Matrix inputMatrix, Vector inputVector)
                {
                    double sum1;
                    double sum2;
                    Vector k1 = new Vector(inputVector.Size);
                    Vector k2 = new Vector(inputVector.Size);
                    for (int k = 0; k < 100; k++)
                    {
                        for (int i = 0; i < k1.Size; i++)
                        {
                            sum1 = 0;
                            sum2 = 0;
                            for (int j = 0; j < i; j++)
                            {
                                sum1 += inputMatrix[i, j] * k2[j];
                            }

                            for (int j = i + 1; j < k1.Size; j++)
                            {
                                sum2 += inputMatrix[i, j] * k1[j];
                            }

                            k2[i] = (inputVector[i] - sum1 - sum2) / inputMatrix[i, i];
                        }
                        if ((k2 - k1).Norm(2) < 0.001) break;
                        k1=new Vector(k2);
                    }

                    return k2;
                }

                // public static Vector SimpleIteration(Matrix inputMatrix, Vector inputVector, out int iterations)
                // {
                //     iterations = 1;
                //     //List<Vector> solves1 = new List<Vector>();
                //     List<Vector> solves = new List<Vector>();
                //     Matrix jacobMatrix = JacobMatrix(inputMatrix);
                //     Vector jacobVector = JacobVector(inputVector, inputMatrix);
                //     if (jacobMatrix.NormInf > 1) throw new Exception("Matrix norm is higher then 1");
                //     // ReSharper disable once HeapView.BoxingAllocation
                //     Console.WriteLine($"Jacob matrix: \n{jacobMatrix} \nJacob vector: {jacobVector}");
                //     Console.WriteLine("Number of iterations: {0}", Iterations(jacobMatrix, jacobVector));
                //     solves.Add(new Vector(inputVector.Size));
                //     //solves1.Add(new Vector(inputVector.Size));
                //     for (int i = 1;; i++)
                //     {
                //         solves.Add(SolveVector(solves[i - 1], jacobVector, jacobMatrix));
                //         //solves1.Add(SolveVectorButDimaSucks(solves[i - 1], jacobVector, jacobMatrix));
                //         //Console.WriteLine(solves[i]);
                //         Console.WriteLine($"Iteration {i} vector: {solves[i]}");
                //         if ((solves[i] - solves[i - 1]).Norm(3) < 0.0001)
                //         {
                //             iterations = i;
                //             break;
                //         }
                //     }
                //
                //     return solves[iterations];
                // }
                //
                // public static Vector SeidelIteration(Matrix inputMatrix, Vector inputVector)
                // {
                //     int iterations = 1;
                //
                //     List<Vector> solves = new List<Vector>();
                //     Matrix jacobMatrix = JacobMatrix(inputMatrix);
                //     Vector jacobVector = JacobVector(inputVector, inputMatrix);
                //     // if (jacobMatrix.NormInf > 1) throw new Exception("Matrix norm is higher then 1");
                //     // Console.WriteLine($"Jacob matrix: \n{jacobMatrix} \nJacob vector: {jacobVector}");
                //     solves.Add(new Vector(inputVector.Size));
                //     for (int i = 1;; i++)
                //     {
                //         solves.Add(SolveSeidelVector(solves[i - 1], jacobVector, jacobMatrix));
                //         if ((solves[i] - solves[i - 1]).Norm(3) < 0.0001)
                //         {
                //             iterations = i;
                //             break;
                //         }
                //     }
                //
                //     return solves[iterations];
                // }

                private static Vector SolveVector(Vector solveVector, Vector jacobVector, Matrix jacobMatrix)
                {
                    Vector solve = new Vector(jacobVector.Size);
                    for (int i = 0; i < jacobVector.Size; i++)
                    {
                        double sum = 0;
                        for (int j = 0; j < jacobVector.Size; j++)
                            sum += solveVector[j] * jacobMatrix[i, j];
                        sum += jacobVector[i];
                        solve[i] = sum;
                    }
                    return solve;
                }

                private static double sum(params Vector[] a)
                {
                    double sum = 0;
                    foreach (var vector in a)
                    {
                        for (int i = 0; i < vector.Size; i++)
                        {
                            sum += vector[i];
                        }
                    }
                    return sum;
                }
                private static Vector SolveVectorButDimaSucks(Vector solveVector, Vector jacobVector,
                    Matrix jacobMatrix)
                {
                    Vector solve = new Vector(jacobVector.Size);
                    for (int i = 0; i < jacobVector.Size; i++)
                    {
                        solve[i] = sum(jacobMatrix[i] % solveVector) + jacobVector[i];
                    }
                    return solve;
                }
                private static Vector SolveSeidelVector(Vector solveVector, Vector jacobVector, Matrix jacobMatrix)
                {
                    Vector solve = new Vector(jacobVector.Size);
                    for (int i = 0; i < jacobVector.Size; i++)
                    {
                        double sum = 0;
                        for (int j = 0; j < jacobVector.Size; j++)
                        {
                            if (j < i - 1) sum += solve[j] * jacobMatrix[i, j];
                            else sum += solveVector[j] * jacobMatrix[i, j];
                        }

                        sum += jacobVector[i];
                        solve[i] = sum;
                    }

                    return solve;
                }

                private static Vector JacobVector(Vector inputVector, Matrix inputMatrix)
                {
                    Vector jacobVector = new Vector(inputVector.Size);
                    for (int i = 0; i < inputVector.Size; i++)
                        jacobVector[i] = inputVector[i] / inputMatrix[i, i];
                    return jacobVector;
                }

                private static Matrix JacobMatrix(Matrix inputMatrix)
                {
                    Matrix jacob = new Matrix(inputMatrix.Row, inputMatrix.Column);
                    for (int i = 0; i < jacob.Row; i++)
                    {
                        for (int j = 0; j < jacob.Column; j++)
                        {
                            if (i != j)
                                jacob[i, j] = -(inputMatrix[i, j] / inputMatrix[i, i]);
                            else jacob[i, j] = 0;
                        }
                    }

                    return jacob;
                }
            }

        }
    }
