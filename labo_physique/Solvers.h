#pragma once

/**
 * @file Solvers.hpp
 *
 * @brief Implémentation de plusieurs algorihtmes de solveurs pour un système
 *        d'équations linéaires
 *
 * Nom:
 * Code permanent :
 * Email :
 *
 */

#include "Math3D.h"

namespace gti320
{
    // Identification des solveurs
    enum eSolverType { kNone, kGaussSeidel, kColorGaussSeidel, kCholesky };

    // Paramètres de convergences pour les algorithmes itératifs
    static const float eps = 1e-4f;
    static const float tau = 1e-5f;

    static bool testConvergence(const Matrix<float, Dynamic, Dynamic>& A,
        const Vector<float, Dynamic>& b,
        Vector<float, Dynamic>& x, Vector<float, Dynamic>& x_prec)
    {
        if ((A * x - b).norm() / b.norm() < eps)
        {
            return false;
        }
        auto deltaX = x - x_prec;
        if (deltaX.norm() / x.norm() < tau)
        {
            return false;
        }
        return true;
    }

    /**
     * Résout Ax = b avec la méthode Gauss-Seidel
     */
    static void gaussSeidel(const Matrix<float, Dynamic, Dynamic>& A,
        const Vector<float, Dynamic>& b,
        Vector<float, Dynamic>& x, int k_max)
    {
        // TODO 
        //
        // Implémenter la méthode de Gauss-Seidel

        for(int k = 0; k < k_max; k++)
        {
            auto x_prec = x;

            for (int i = 0; i < A.rows(); i++)
            {
                x(i) = b(i);
                for (int j = 0; i < i - 1; j++) 
                {
                    x(i) = x(i) - A(i, j) * x(j);
                }
                for (int j = i + 1; j < A.rows(); j++)
                {
                    x(i) = x(i) - A(i, j) * x(j);
                }
                x(i) = x(i) / A(i, i);
            }

            if (!testConvergence(A, b, x, x_prec))
            {
                return;
            }
        }
    }
    

    /**

     * Résout Ax = b avec la méthode Gauss-Seidel (coloration de graphe)
     */
    static void gaussSeidelColor(const Matrix<float, Dynamic, Dynamic>& A, const Vector<float, Dynamic>& b, Vector<float, Dynamic>& x, const Partitions& P, const int maxIter)
    {
        // TODO 
        //
        // Implémenter la méthode de Gauss-Seidel avec coloration de graphe.
        // Les partitions avec l'index de chaque particule sont stockées dans la table des tables, P.

    }

    /**
     * Résout Ax = b avec la méthode de Cholesky
     */
    static void cholesky(const Matrix<float, Dynamic, Dynamic>& A,
        const Vector<float, Dynamic>& b,
        Vector<float, Dynamic>& x)
    {
        // TODO 
        //
        // Calculer la matrice L de la factorisation de Cholesky



        // TODO
        //
        // Résoudre Ly = b



        // TODO
        //
        // Résoudre L^t x = y
        // 
        // Remarque : ne pas caculer la transposer de L, c'est inutilement
        // coûteux.



    }

}
