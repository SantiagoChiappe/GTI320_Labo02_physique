#pragma once

/**
 * @file Math3D.h
 *
 * @brief Fonctions pour l'initialisation et la manipulation de matrices de
 * rotation, de matrices de transformations en coordonnées homogènes et de
 * vecteurs 3D.
 *
 * Nom: Erick Santiago Chiappe Reyes
 * Code permanent : CHIE19079806
 * Email : erick-santiago.chiappe-reyes.1@ens.etsmtl.ca
 *
 */

#include "Matrix.h"
#include "Operators.h"
#include "Vector.h"

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <math.h>

namespace gti320 {

    // Deux types de vecteurs 3D considérés ici
    typedef Vector<double, 3> Vector3d;
    typedef Vector<float, 3> Vector3f;

    // Dans le cadre de ce projet, nous considérons seulement deux
    // cas :
    //
    //  - les rotations
    //  - les translations
    //
    // Deux types de matrices en coordonnées homogèes :
    typedef Matrix<double, 4, 4, ColumnStorage> Matrix4d;
    typedef Matrix<float, 4, 4, ColumnStorage> Matrix4f;
    //
    // Deux types de matrices pour les rotations
    typedef Matrix<double, 3, 3, ColumnStorage> Matrix3d;
    typedef Matrix<float, 3, 3, ColumnStorage> Matrix3f;

    /** Initialise et retourne la matrice identité */
    template<>
    inline void Matrix4d::setIdentity()
    {
        /** TODO: Redéfinir la fonction membre de façon optimisée pour les Matrix4d */
        assert(this->rows() > 0 && this->cols() > 0);
        this->setZero();
        for (int i = 0; i < this->rows(); i++)
        {
            (*this)(i, i) = 1;
        }
    }

    /**
     * Calcul de la matrice inverse SPÉCIALISÉ pour le cas
     * d'une matrice de transformation en coordonnées homogènes.
     */
    template<>
    inline Matrix4d Matrix4d::inverse() const
    {
        /** TODO: Implémenter. */
        assert(this->rows() > 0 && this->cols() > 0);
        Matrix4d matrice = Matrix4d(*this);
        auto submatriceRotation = matrice.block(0, 0, 3, 3);
        submatriceRotation = submatriceRotation.transpose<double, 3, 3, ColumnStorage>();
        auto submatriceTranslate = matrice.block(0, 3, 3, 1);
        
        auto tx = submatriceTranslate(0, 0);
        auto ty = submatriceTranslate(1, 0);
        auto tz = submatriceTranslate(2, 0);

        submatriceTranslate(0, 0) = -(tx * submatriceRotation(0, 0) + ty * submatriceRotation(0, 1) + tz * submatriceRotation(0, 2));
        submatriceTranslate(1, 0) = -(tx * submatriceRotation(1, 0) + ty * submatriceRotation(1, 1) + tz * submatriceRotation(1, 2));
        submatriceTranslate(2, 0) = -(tx * submatriceRotation(2, 0) + ty * submatriceRotation(2, 1) + tz * submatriceRotation(2, 2));

        return matrice;
    }

    /**
     * Calcul de la matrice inverse SPÉCIALISÉ pour le cas
     * d'une matrice de rotation.
     */
    template<>
    inline Matrix3d Matrix3d::inverse() const
    {
        /** TODO: Implémenter. */
        assert(this->rows() > 0 && this->cols() > 0);
        return this->transpose<double, 3, 3, ColumnStorage>();
    }

    /**
     * Multiplication d'une matrice 4x4 avec un vecteur 3D où la matrice
     * représente une transformation en coordonnées homogène.
     */
    template<typename _Scalar>
    Vector<_Scalar, 3> operator*(const Matrix<_Scalar, 4, 4, ColumnStorage>& A, const Vector<_Scalar, 3>& v)
    {
        /** TODO: Implémenter. */
        assert(A.cols() == v.size());
        Vector<_Scalar, 4> vectorHomogene = Vector<_Scalar, 4>();
        for (int i = 0; i < v.size(); i++) 
        {
            vectorHomogene(i) = v(i);
        }
        vectorHomogene(3) = 1;
        vectorHomogene = A * vectorHomogene;
        Vector<_Scalar, 3> vectorTransformee = Vector<_Scalar, 3>();
        for (int i = 0; i < v.size(); i++)
        {
            vectorTransformee(i) = vectorHomogene(i);
        }
        return vectorTransformee;
    }

    /**
     * Créer et retourne une matrice de rotation définie par les angles
     * d'Euler XYZ exprimés en radians.
     *
     * La matrice doit correspondre au produit : Rz * Ry * Rx.
     */
    template<typename _Scalar>
    static Matrix<_Scalar, 3, 3> makeRotation(_Scalar x, _Scalar y, _Scalar z)
    {
        /** TODO: Implémenter. */
        double sinX = sin(x);
        double sinY = sin(y);
        double sinZ = sin(z);
        double cosX = cos(x);
        double cosY = cos(y);
        double cosZ = cos(z);

        Matrix<_Scalar, 3, 3> matrice = Matrix<_Scalar, 3, 3>();

        matrice(0, 0) = cosY * cosZ;
        matrice(1, 0) = cosY * sinZ;
        matrice(2, 0) = -sinY;

        matrice(0, 1) = -cosX * sinZ + sinX * sinY * cosZ;
        matrice(1, 1) = cosX * cosZ + sinX * sinY * sinZ;
        matrice(2, 1) = sinX * cosY;

        matrice(0, 2) = sinX * sinZ + cosX * sinY * cosZ;
        matrice(1, 2) = -sinX * cosZ + cosX * sinY * sinZ;
        matrice(2, 2) = cosX * cosY;
        return matrice;
    }
} // namespace gti320
