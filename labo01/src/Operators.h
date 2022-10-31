#pragma once

/**
 * @file Operators.h
 *
 * @brief Implémentation de divers opérateurs arithmétiques pour les matrices et les vecteurs.
 *
 * Nom: Erick Santiago Chiappe Reyes
 * Code permanent : CHIE19079806
 * Email : erick-santiago.chiappe-reyes.1@ens.etsmtl.ca
 *
 */

#include "Matrix.h"
#include "Vector.h"

namespace gti320 {

    /** Multiplication : Matrice * Matrice (générique) */
    template<typename _Scalar, int RowsA, int ColsA, int StorageA, int RowsB, int ColsB, int StorageB>
    Matrix<_Scalar, RowsA, ColsB> operator*(const Matrix<_Scalar, RowsA, ColsA, StorageA>& A, const Matrix<_Scalar, RowsB, ColsB, StorageB>& B)
    {
        /** TODO: Implémenter. */
        assert(A.cols() == B.rows());
        Matrix<_Scalar, RowsA, ColsB> matrice = Matrix<_Scalar, RowsA, ColsB>(A.rows(), B.cols());
        for (int j = 0; j < B.rows(); j++)
        {
            for (int k = 0; k < A.cols(); k++)
            {
                for (int i = 0; i < A.rows(); i++)
                {
                    matrice(i, j) += A(i, k) * B(k, j);
                }
            }
        }
        return matrice;
    }

    /**
     * Multiplication : Matrice (colonne) * Matrice (ligne)
     *
     * Spécialisation de l'opérateur de multiplication pour le cas où les matrices
     * ont un stockage à taille dynamique et où la matrice de gauche utilise un
     * stockage par colonnes et celle de droite un stockage par lignes.
     */
    template<typename _Scalar>
    Matrix<_Scalar, Dynamic, Dynamic> operator*(const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& A, const Matrix<_Scalar, Dynamic, Dynamic, RowStorage>& B)
    {
        /** TODO: Implémenter. */
        assert(A.cols() == B.rows());
        Matrix<_Scalar, Dynamic, Dynamic> matrice = Matrix<_Scalar, Dynamic, Dynamic>(A.rows(), B.cols());
        for (int k = 0; k < A.cols(); k++)
        {
            for (int j = 0; j < B.rows(); j++)
            {
                for (int i = 0; i < A.rows(); i++)
                {
                    matrice(i, j) += A(i, k) * B(k, j);
                }
            }
        }
        return matrice;
    }

    /**
     * Multiplication : Matrice (ligne) * Matrice (colonne)
     *
     * Spécialisation de l'opérateur de multiplication pour le cas où les matrices
     * ont un stockage à taille dynamique et où la matrice de gauche utilise un
     * stockage par lignes et celle de droite un stockage par colonnes.
     */
    template<typename _Scalar>
    Matrix<_Scalar, Dynamic, Dynamic> operator*(const Matrix<_Scalar, Dynamic, Dynamic, RowStorage>& A, const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& B)
    {
        /** TODO: Implémenter. */
        assert(A.cols() == B.rows());
        Matrix<_Scalar, Dynamic, Dynamic> matrice = Matrix<_Scalar, Dynamic, Dynamic>(A.rows(), B.cols());
        for (int j = 0; j < B.rows(); j++)
        {
            for (int i = 0; i < A.rows(); i++)
            {
                for (int k = 0; k < A.cols(); k++)
                {
                    matrice(i, j) += A(i, k) * B(k, j);
                }
            }
        }
        return matrice;
    }

    /** Addition : Matrice + Matrice (générique) */
    template<typename _Scalar, int Rows, int Cols, int StorageA, int StorageB>
    Matrix<_Scalar, Rows, Cols> operator+(const Matrix<_Scalar, Rows, Cols, StorageA>& A, const Matrix<_Scalar, Rows, Cols, StorageB>& B)
    {
        /** TODO: Implémenter. */
        assert(A.cols() == B.cols() && A.rows() == B.rows());
        Matrix<_Scalar, Rows, Cols> matrice = Matrix<_Scalar, Rows, Cols>(A.rows(), A.cols());
        for (int j = 0; j < A.cols(); j++)
        {
            for (int i = 0; i < A.rows(); i++)
            {
                matrice(i, j) = A(i, j) + B(i, j);
            }
        }
        return matrice;
    }

    /**
     * Addition : Matrice (colonne) + Matrice (colonne)
     *
     * Spécialisation de l'opérateur d'addition pour le cas où les deux matrices
     * sont stockées par colonnes.
     */
    template<typename _Scalar>
    Matrix<_Scalar, Dynamic, Dynamic> operator+(const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& A, const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& B)
    {
        /** TODO: Implémenter. */
        assert(A.cols() == B.cols() && A.rows() == B.rows());
        Matrix<_Scalar, Dynamic, Dynamic> matrice = Matrix<_Scalar, Dynamic, Dynamic>(A.rows(), A.cols());
        for (int j = 0; j < A.cols(); j++)
        {
            for (int i = 0; i < A.rows(); i++)
            {
                matrice(i, j) = A(i, j) + B(i, j);
            }
        }
        return matrice;
    }

    /**
     * Addition : Matrice (ligne) + Matrice (ligne)
     *
     * Spécialisation de l'opérateur d'addition pour le cas où les deux matrices
     * sont stockées par lignes.
     */
    template<typename _Scalar>
    Matrix<_Scalar, Dynamic, Dynamic, RowStorage> operator+(const Matrix<_Scalar, Dynamic, Dynamic, RowStorage>& A, const Matrix<_Scalar, Dynamic, Dynamic, RowStorage>& B)
    {
        /** TODO: Implémenter. */
        assert(A.cols() == B.cols() && A.rows() == B.rows());
        Matrix<_Scalar, Dynamic, Dynamic, RowStorage> matrice = Matrix<_Scalar, Dynamic, Dynamic, RowStorage>(A.rows(), A.cols());
        for (int i = 0; i < A.rows(); i++) 
        {
            for (int j = 0; j < A.cols(); j++) 
            {
                matrice(i, j) = A(i, j) + B(i, j);
            }
        }
        return matrice;
    }

    /**
     * Multiplication  : Scalaire * Matrice (colonne)
     *
     * Spécialisation de l'opérateur de multiplication par un scalaire pour le
     * cas d'une matrice stockée par colonnes.
     */
    template<typename _Scalar, int _Rows, int _Cols>
    Matrix<_Scalar, _Rows, _Cols, ColumnStorage> operator*(const _Scalar& a, const Matrix<_Scalar, _Rows, _Cols, ColumnStorage>& A)
    {
        /** TODO: Implémenter. */
        assert(A.cols() > 0 && A.rows() > 0);
        Matrix<_Scalar, _Rows, _Cols, ColumnStorage> matrice = Matrix<_Scalar, _Rows, _Cols, ColumnStorage>(A.rows(), A.cols());
        for (int j = 0; j < A.cols(); j++)
        {
            for (int i = 0; i < A.rows(); i++)
            {
                matrice(i, j) = A(i, j) * a;
            }
        }
        return matrice;
    }

    /**
     * Multiplication  : Scalaire * Matrice (ligne)
     *
     * Spécialisation de l'opérateur de multiplication par un scalaire pour le
     * cas d'une matrice stockée par lignes.
     */
    template<typename _Scalar, int _Rows, int _Cols>
    Matrix<_Scalar, _Rows, _Cols, RowStorage> operator*(const _Scalar& a, const Matrix<_Scalar, _Rows, _Cols, RowStorage>& A)
    {
        /** TODO: Implémenter. */
        assert(A.cols() > 0 && A.rows() > 0);
        Matrix<_Scalar, _Rows, _Cols, RowStorage> matrice = Matrix<_Scalar, _Rows, _Cols, RowStorage>(A.rows(), A.cols());
        for (int i = 0; i < A.rows(); i++)
        {
            for (int j = 0; j < A.cols(); j++)
            {
                matrice(i, j) = A(i, j) * a;
            }
        }
        return matrice;
    }

    /**
     * Multiplication : Matrice (ligne) * Vecteur
     *
     * Spécialisation de l'opérateur de multiplication matrice*vecteur pour le
     * cas où la matrice est représentée par lignes.
     */
    template<typename _Scalar, int _Rows, int _Cols>
    Vector<_Scalar, _Rows> operator*(const Matrix<_Scalar, _Rows, _Cols, RowStorage>& A, const Vector<_Scalar, _Cols>& v)
    {
        /** TODO: Implémenter. */
        assert(A.cols() == v.size());
        Vector<_Scalar, _Rows> vector = Vector<_Scalar, _Rows>(A.rows());
        for (int i = 0; i < A.rows(); i++)
        {
            for (int j = 0; j < A.cols(); j++) 
            {
                vector(i) += A(i, j) * v(j);
            }
        }
        return vector;
    }

    /**
     * Multiplication : Matrice (colonne) * Vecteur
     *
     * Spécialisation de l'opérateur de multiplication matrice*vecteur pour le
     * cas où la matrice est représentée par colonnes.
     */
    template<typename _Scalar, int _Rows, int _Cols>
    Vector<_Scalar, _Rows> operator*(const Matrix<_Scalar, _Rows, _Cols, ColumnStorage>& A, const Vector<_Scalar, _Cols>& v)
    {
        /** TODO: Implémenter. */ 
        assert(A.cols() == v.size());
        Vector<_Scalar, _Rows> vector = Vector<_Scalar, _Rows>(A.rows());
        for (int j = 0; j < A.cols(); j++)
        {
            for (int i = 0; i < A.rows(); i++)
            {
                vector(i) += A(i, j) * v(j);
            }
        }
        return vector;
    }

    /** Multiplication : Scalaire * Vecteur */
    template<typename _Scalar, int _Rows>
    Vector<_Scalar, _Rows> operator*(const _Scalar& a, const Vector<_Scalar, _Rows>& v)
    {
        /** TODO: Implémenter. */
        assert(v.size() > 0);
        Vector<_Scalar, _Rows> vector = Vector<_Scalar, _Rows>(v.size());
        for (int i = 0; i < v.size(); i++)
        {
            vector(i) = v(i) * a;
        }
        return vector;
    }

    /** Addition : Vecteur + Vecteur */
    template<typename _Scalar, int _RowsA, int _RowsB>
    Vector<_Scalar, _RowsA> operator+(const Vector<_Scalar, _RowsA>& a, const Vector<_Scalar, _RowsB>& b)
    {
        /** TODO: Implémenter. */
        assert(a.size() == b.size());
        Vector<_Scalar, _RowsA> vector = Vector<_Scalar, _RowsA>(a.size());
        for (int i = 0; i < a.size(); i++)
        {
            vector(i) = a(i) + b(i);
        }
        return vector;
    }

    /** Soustraction : Vecteur - Vecteur */
    template<typename _Scalar, int _RowsA, int _RowsB>
    Vector<_Scalar, _RowsA> operator-(const Vector<_Scalar, _RowsA>& a, const Vector<_Scalar, _RowsB>& b)
    {
        /** TODO: Implémenter. */
        assert(a.size() == b.size());
        Vector<_Scalar, _RowsA> vector = Vector<_Scalar, _RowsA>(a.size());
        for (int i = 0; i < a.size(); i++)
        {
            vector(i) = a(i) - b(i);
        }
        return vector;
    }

    /** Multiplication : Vecteur * Vecteur (produit dyadique) */
    template<typename _Scalar, int _Rows>
    Matrix<_Scalar, _Rows, _Rows> operator*(const Vector<_Scalar, _Rows>& v1, const Vector<_Scalar, _Rows>& v2)
    {
        assert(v1.size() > 0 && v2.size() > 0 && v1.size() == v2.size());
        Matrix<_Scalar, _Rows, _Rows> matrix = Matrix<_Scalar, _Rows, _Rows>(v1.size(), v1.size());
        for (int j = 0; j < v1.size(); j++)
        {
            for (int i = 0; i < v1.size(); i++)
            {
                matrix(i, j) += v1(i) * v1(j);
            }
        }
        return matrix;
    }
} // namespace gti320
