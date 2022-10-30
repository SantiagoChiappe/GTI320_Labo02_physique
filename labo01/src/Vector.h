#pragma once

/**
 * @file Vector.h
 *
 * @brief Implémentation de vecteurs simples
 *
 * Nom: Erick Santiago Chiappe Reyes
 * Code permanent : CHIE19079806
 * Email : erick-santiago.chiappe-reyes.1@ens.etsmtl.ca
 *
 */

#include "MatrixBase.h"

namespace gti320 {

    /**
     * Classe vecteur générique.
     *
     * Cette classe réutilise la classe `MatrixBase` et ses spécialisations de
     * templates pour les manipulation bas niveau.
     */
    template<typename _Scalar = double, int _Rows = Dynamic>
    class Vector : public MatrixBase<_Scalar, _Rows, 1> {
    public:
        /** Constructeur par défaut */
        Vector() : MatrixBase<_Scalar, _Rows, 1>() {}

        /** Contructeur à partir d'un taille (rows). */
        explicit Vector(int rows) : MatrixBase<_Scalar, _Rows, 1>(rows, 1) {}

        /** Constructeur de copie */
        Vector(const Vector& other) : MatrixBase<_Scalar, _Rows, 1>(other) {}

        /** Destructeur */
        ~Vector() {}

        /** Opérateur de copie */
        Vector& operator=(const Vector& other)
        {
            /** TODO: Implémenter. */
            MatrixBase<_Scalar, _Rows, 1>::operator=(other);
            return *this;
        }

        /** Accesseur à une entrée du vecteur (lecture seule) */
        _Scalar operator()(int i) const
        {
            /** TODO: Implémenter. */
            return this->data()[i];
        }

        /** Accesseur à une entrée du vecteur (lecture et écriture) */
        _Scalar& operator()(int i)
        {
            /** TODO: Implémenter. */
            return this->m_storage.data()[i];
        }

        /** Modifie le nombre de lignes du vecteur */
        void resize(int _rows)
        {
            MatrixBase<_Scalar, _Rows, 1>::resize(_rows, 1);
        }

        /** Produit scalaire de *this et other */
        inline _Scalar dot(const Vector& other) const
        {
            /** TODO: Implémenter. */
            assert(this->size() == other.size());
            _Scalar dot = 0;
            for (int i = 0; i < this->size(); i++)
            {
                dot += (*this)(i) * other(i);
            }
            return dot;
        }

        /** Retourne la norme euclidienne du vecteur */
        inline _Scalar norm() const
        {
            /** TODO: Implémenter. */
            assert(this->size() > 0);
            _Scalar norm = 0;
            for (int i = 0; i < this->size(); i++)
            {
                norm += pow((*this)(i), 2);
            }
            return sqrt(norm);
        }
    };
} // namespace gti320
