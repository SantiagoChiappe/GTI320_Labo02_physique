#pragma once

/**
 * @file Solvers.hpp
 *
 * @brief Implémentation de plusieurs algorihtmes de solveurs pour un système
 *        d'équations linéaires
 *
 * Nom: Erick Santiago Chiappe Reyes
 * Code permanent : CHIE19079806
 * Email : erick-santiago.chiappe-reyes.1@ens.etsmtl.ca
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
            return true;
        }
        auto deltaX = x - x_prec;
        if (deltaX.norm() / x.norm() < tau)
        {
            return true;
        }
        return false;
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

        for (int k = 0; k < k_max; k++)
        {
            auto x_prec = x;

            for (int i = 0; i < A.rows(); i++)
            {
                x(i) = b(i);
                for (int j = 0; j < i; j++)
                {
                    x(i) = x(i) - A(i, j) * x(j);
                }
                for (int j = i + 1; j < A.rows(); j++)
                {
                    x(i) = x(i) - A(i, j) * x(j);
                }
                x(i) = x(i) / A(i, i);
            }

            if (testConvergence(A, b, x, x_prec)) return;
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
        
        bool converged = false;

        while (!converged) 
        {
            for (int max = 0; max < maxIter; max++)
            {
                auto x_prec = x;

                for (int c = 0; c < P.size(); ++c)
                {
                    const std::vector<int>& inds = P[c];
                    #pragma omp parallel for

                    for (int k = 0; k < inds.size(); ++k)
                    {
                        // iteration de Gauss-Seidel

                        for (int ii = 0; ii < 2; ii++)
                        {
                            int i = inds[k] * 2 + ii;
                            x(i) = b(i);
                            for (int j = 0; j < i; j++)
                            {
                                x(i) = x(i) - A(i, j) * x(j);
                            }
                            for (int j = i + 1; j < A.rows(); j++)
                            {
                                x(i) = x(i) - A(i, j) * x(j);
                            }
                            x(i) = x(i) / A(i, i);
                        }

                        converged = testConvergence(A, b, x, x_prec);
                    }
                }
            }
        }
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

        Matrix<float, Dynamic, Dynamic> L(A.rows(), A.cols());
        L.setZero();

        for (int i = 0; i < A.rows(); i++) 
        {
            for (int k = 0; k <= i; k++) 
            {
                int s = 0;
                for (int j = 0; j < k; j++) 
                {
                    s = s + L(i, j) * L(k, j);
                }
                if (i == k) 
                {
                    L(i, k) = sqrt(A(i, i) - s);
                }
                else 
                {
                    L(i, k) = (A(i, k) - s) / L(k, k);
                }
            }
        }

        Vector<float, Dynamic> y(A.rows());
        y.setZero();

        // TODO
        //
        // Résoudre Ly = b

        for (int i = 0; i < A.rows(); i++) 
        {
            y(i) = b(i);
            for (int j = 0; j < i; j++) 
            {
                y(i) = y(i) - L(i, j) * y(j);
            }
            y(i) = y(i) / L(i, i);
        }

        // TODO
        //
        // Résoudre L^t x = y
        // 
        // Remarque : ne pas caculer la transposer de L, c'est inutilement
        // coûteux.

        for (int i = A.rows() - 1; i >= 0; i--) 
        {
            x(i) = y(i);
            for (int j = i + 1; j < A.rows(); j++) 
            {
                x(i) = x(i) - L(j, i) * x(j);
            }
            x(i) = x(i) / L(i, i);
        }
    }

}
