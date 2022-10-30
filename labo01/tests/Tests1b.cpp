/**
 * @file Tests1b.cpp
 *
 * @brief Tests unitaires de la partie 1b
 *
 * Nom:
 * Code permanent :
 * Email :
 *
 */

#include "Matrix.h"
#include "Vector.h"
#include "Math3D.h"
#include "Operators.h"

#include <gtest/gtest.h>
#include <chrono>

using namespace gti320; 

namespace
{
    /**
     * Multiplication  matrice * vecteur,  utilisant une implémentation naive
     */
    template<typename _Scalar>
    static inline Vector<_Scalar, Dynamic> naiveMatrixMult(const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& A, const Vector<_Scalar, Dynamic>& v)
    {
        assert(A.cols() == v.rows());

        Vector<_Scalar, Dynamic> b(A.rows());
        assert(b.rows() == A.rows());

        for (int i = 0; i < A.rows(); ++i) {
            b(i) = 0.0;
            for (int j = 0; j < A.cols(); ++j) {
                b(i) += A(i, j) * v(j);
            }
        }

        return b;
    }

    /**
     * Addition  matrice + matrice,  utilisant une implémentation naive
     */
    template<typename _Scalar>
    static inline Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage> naiveMatrixAddition(const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& A, const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& B)
    {
        assert(A.cols() == B.cols() && A.rows() == B.rows());

        Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage> C(A.rows(), A.cols());
        assert(C.rows() == A.rows() && C.cols() == A.cols());
        for (int i = 0; i < C.rows(); ++i) {
            for (int j = 0; j < C.cols(); ++j) {
                C(i, j) = A(i, j) + B(i, j);
            }
        }
        return C;
    }

    /**
     * Multiplication  matrice * matrice,  utilisant une implémentation naive.
     */
    template<typename _Scalar, int _Storage>
    static inline Matrix<_Scalar, Dynamic, Dynamic, _Storage> naiveMatrixMult(const Matrix<_Scalar, Dynamic, Dynamic, _Storage>& A, const Matrix<_Scalar, Dynamic, Dynamic, _Storage>& B)
    {
        assert(A.cols() == B.rows());
        Matrix<_Scalar, Dynamic, Dynamic> product(A.rows(), B.cols());
        for (int i = 0; i < A.rows(); ++i)
        {
            for (int j = 0; j < B.cols(); ++j)
            {
                for (int k = 0; k < A.cols(); ++k)
                {
                    product(i, j) += A(i, k) * B(k, j);
                }
            }
        }
        return product;
    }
}

/**
 * Test pour les opérateurs d'arithmétique matricielle.
 */
TEST(TestLabo1b, MatrixMatrixOperators)
{
    // Opérations arithmétiques avec matrices à taille dynamique
    {
        // Test : matrice identité
        Matrix<double> A(6, 6);
        A.setIdentity();
        EXPECT_DOUBLE_EQ(A(0, 0), 1.0);
        EXPECT_DOUBLE_EQ(A(1, 1), 1.0);
        EXPECT_DOUBLE_EQ(A(2, 2), 1.0);
        EXPECT_DOUBLE_EQ(A(3, 3), 1.0);
        EXPECT_DOUBLE_EQ(A(4, 4), 1.0);
        EXPECT_DOUBLE_EQ(A(5, 5), 1.0);
        EXPECT_DOUBLE_EQ(A(0, 1), 0.0);
        EXPECT_DOUBLE_EQ(A(1, 0), 0.0);

        // Test : produit  scalaire * matrice
        const double alpha = 2.5;
        Matrix<double> B = alpha * A;
        EXPECT_DOUBLE_EQ(B(0, 0), alpha);
        EXPECT_DOUBLE_EQ(B(1, 1), alpha);
        EXPECT_DOUBLE_EQ(B(2, 2), alpha);
        EXPECT_DOUBLE_EQ(B(3, 3), alpha);
        EXPECT_DOUBLE_EQ(B(4, 4), alpha);
        EXPECT_DOUBLE_EQ(B(5, 5), alpha);
        EXPECT_DOUBLE_EQ(B(0, 1), 0.0);
        EXPECT_DOUBLE_EQ(B(1, 0), 0.0);

        // Test : produit  matrice * matrice
        Matrix<double> C = A * B;
        EXPECT_DOUBLE_EQ(C(0, 0), A(0, 0) * B(0, 0));
        EXPECT_DOUBLE_EQ(C(1, 1), A(1, 1) * B(1, 1));
        EXPECT_DOUBLE_EQ(C(2, 2), A(2, 2) * B(2, 2));
        EXPECT_DOUBLE_EQ(C(3, 3), A(3, 3) * B(3, 3));
        EXPECT_DOUBLE_EQ(C(4, 4), A(4, 4) * B(4, 4));
        EXPECT_DOUBLE_EQ(C(5, 5), A(5, 5) * B(5, 5));
        EXPECT_DOUBLE_EQ(C(0, 1), 0.0);
        EXPECT_DOUBLE_EQ(C(2, 3), 0.0);

        // Test : addition  matrice * matrice
        Matrix<double> A_plus_B = A + B;
        EXPECT_DOUBLE_EQ(A_plus_B(0, 0), A(0, 0) + B(0, 0));
        EXPECT_DOUBLE_EQ(A_plus_B(1, 1), A(1, 1) + B(1, 1));
        EXPECT_DOUBLE_EQ(A_plus_B(2, 2), A(2, 2) + B(2, 2));
        EXPECT_DOUBLE_EQ(A_plus_B(3, 3), A(3, 3) + B(3, 3));
        EXPECT_DOUBLE_EQ(A_plus_B(4, 4), A(4, 4) + B(4, 4));
        EXPECT_DOUBLE_EQ(A_plus_B(5, 5), A(5, 5) + B(5, 5));
        EXPECT_DOUBLE_EQ(A_plus_B(0, 1), 0.0);
        EXPECT_DOUBLE_EQ(A_plus_B(2, 3), 0.0);
    }

    // Opérations arithmétique avec matrices à stockage par lignes et par
    // colonnes.
    {
        // Création d'un matrice à stockage par lignes
        Matrix<double, Dynamic, Dynamic, RowStorage> A(5, 5);
        A(0, 0) = 0.8147;    A(0, 1) = 0.0975;    A(0, 2) = 0.1576;    A(0, 3) = 0.1419;    A(0, 4) = 0.6557;
        A(1, 0) = 0.9058;    A(1, 1) = 0.2785;    A(1, 2) = 0.9706;    A(1, 3) = 0.4218;    A(1, 4) = 0.0357;
        A(2, 0) = 0.1270;    A(2, 1) = 0.5469;    A(2, 2) = 0.9572;    A(2, 3) = 0.9157;    A(2, 4) = 0.8491;
        A(3, 0) = 0.9134;    A(3, 1) = 0.9575;    A(3, 2) = 0.4854;    A(3, 3) = 0.7922;    A(3, 4) = 0.9340;
        A(4, 0) = 0.6324;    A(4, 1) = 0.9649;    A(4, 2) = 0.8003;    A(4, 3) = 0.9595;    A(4, 4) = 0.6787;

        // Test : transposée (le résultat est une matrice à stockage par
        //        colonnes)
        Matrix<double, Dynamic, Dynamic, ColumnStorage> B = A.transpose();

        // Test : multiplication  matrix(ligne) * matrice(colonne)
        // Note : teste seulement la première et la dernière colonne
        const auto C = A * B;
        EXPECT_NEAR(C(0, 0), 1.14815820000000, 1e-3); EXPECT_NEAR(C(0, 4), 1.31659795000000, 1e-3);
        EXPECT_NEAR(C(1, 0), 1.00133748000000, 1e-3); EXPECT_NEAR(C(1, 4), 2.04727044000000, 1e-3);
        EXPECT_NEAR(C(2, 0), 0.99433707000000, 1e-3); EXPECT_NEAR(C(2, 4), 2.82896409000000, 1e-3);
        EXPECT_NEAR(C(3, 0), 1.63883925000000, 1e-3); EXPECT_NEAR(C(3, 4), 3.28401323000000, 1e-3);
        EXPECT_NEAR(C(4, 0), 1.31659795000000, 1e-3); EXPECT_NEAR(C(4, 4), 3.35271580000000, 1e-3);


        // Test : multiplication  matrice(colonne) * matrice(ligne)
        // Note : teste seulement la première et la dernière colonne
        const auto C2 = B * A;
        EXPECT_NEAR(C2(0, 0), 2.73456805000000, 1e-3); EXPECT_NEAR(C2(0, 4), 1.95669703000000, 1e-3);
        EXPECT_NEAR(C2(1, 0), 1.88593811000000, 1e-3); EXPECT_NEAR(C2(1, 4), 2.08742862000000, 1e-3);
        EXPECT_NEAR(C2(2, 0), 2.07860468000000, 1e-3); EXPECT_NEAR(C2(2, 4), 1.94727447000000, 1e-3);
        EXPECT_NEAR(C2(3, 0), 1.94434955000000, 1e-3); EXPECT_NEAR(C2(3, 4), 2.27675041000000, 1e-3);
        EXPECT_NEAR(C2(4, 0), 1.95669703000000, 1e-3); EXPECT_NEAR(C2(4, 4), 2.48517748000000, 1e-3);

        // Test : addition  matrice(ligne) + matrice(ligne)
        // Note : teste seulement la première et la dernière colonne
        const auto A_plus_A = A + A;
        EXPECT_DOUBLE_EQ(A_plus_A(0, 0), A(0, 0) + A(0, 0)); EXPECT_DOUBLE_EQ(A_plus_A(0, 4), A(0, 4) + A(0, 4));
        EXPECT_DOUBLE_EQ(A_plus_A(1, 0), A(1, 0) + A(1, 0)); EXPECT_DOUBLE_EQ(A_plus_A(1, 4), A(1, 4) + A(1, 4));
        EXPECT_DOUBLE_EQ(A_plus_A(2, 0), A(2, 0) + A(2, 0)); EXPECT_DOUBLE_EQ(A_plus_A(2, 4), A(2, 4) + A(2, 4));
        EXPECT_DOUBLE_EQ(A_plus_A(3, 0), A(3, 0) + A(3, 0)); EXPECT_DOUBLE_EQ(A_plus_A(3, 4), A(3, 4) + A(3, 4));
        EXPECT_DOUBLE_EQ(A_plus_A(4, 0), A(4, 0) + A(4, 0)); EXPECT_DOUBLE_EQ(A_plus_A(4, 4), A(4, 4) + A(4, 4));

        // Test : addition  matrice(colonne) + matrice(colonne)
        // Note : teste seulement la première et la dernière colonne
        const auto B_plus_B = B + B;
        EXPECT_DOUBLE_EQ(B_plus_B(0, 0), B(0, 0) + B(0, 0)); EXPECT_DOUBLE_EQ(B_plus_B(0, 4), B(0, 4) + B(0, 4));
        EXPECT_DOUBLE_EQ(B_plus_B(1, 0), B(1, 0) + B(1, 0)); EXPECT_DOUBLE_EQ(B_plus_B(1, 4), B(1, 4) + B(1, 4));
        EXPECT_DOUBLE_EQ(B_plus_B(2, 0), B(2, 0) + B(2, 0)); EXPECT_DOUBLE_EQ(B_plus_B(2, 4), B(2, 4) + B(2, 4));
        EXPECT_DOUBLE_EQ(B_plus_B(3, 0), B(3, 0) + B(3, 0)); EXPECT_DOUBLE_EQ(B_plus_B(3, 4), B(3, 4) + B(3, 4));
        EXPECT_DOUBLE_EQ(B_plus_B(4, 0), B(4, 0) + B(4, 0)); EXPECT_DOUBLE_EQ(B_plus_B(4, 4), B(4, 4) + B(4, 4));

    }
}

/**
 * Test pour la multiplication  matrice * vecteur
 */
TEST(TestLabo1b, MatrixVectorOperators)
{
    // Vecteur à taille dynamique
    Vector<double> v(5);
    v(0) = 1.0;
    v(1) = 2.0;
    v(2) = 4.0;
    v(3) = 8.0;
    v(4) = 16.0;

    // Test : multiplication par la matrice identité
    {
        Matrix<double> M(5, 5);
        M.setIdentity();

        const auto b = M * v;
        EXPECT_DOUBLE_EQ(b(0), 1.0);
        EXPECT_DOUBLE_EQ(b(1), 2.0);
        EXPECT_DOUBLE_EQ(b(2), 4.0);
        EXPECT_DOUBLE_EQ(b(3), 8.0);
        EXPECT_DOUBLE_EQ(b(4), 16.0);
    }

    // Test : multiplication par une matrice à taille dynamique avec stockage par ligne.
    {
        Matrix<double, Dynamic, Dynamic, RowStorage> M(5, 5);
        M.setIdentity();
        M = 2.0 * M;

        Vector<double> b2 = M * v;
        EXPECT_DOUBLE_EQ(b2(0), 2.0);
        EXPECT_DOUBLE_EQ(b2(1), 4.0);
        EXPECT_DOUBLE_EQ(b2(2), 8.0);
        EXPECT_DOUBLE_EQ(b2(3), 16.0);
        EXPECT_DOUBLE_EQ(b2(4), 32.0);
    }
    // Test : autres
    {
        Vector<double, 3> u; u(0) = 89.63885274026418; u(1) = 33.98267508829284; u(2) = 58.18085749457168;
        Matrix<double, 3, 3, ColumnStorage> m33;
        m33(0, 0) = -37.761291293282206; m33(0, 1) = -58.77833883088219; m33(0, 2) = -80.31292383738193; m33(1, 0) = 82.42200051847567; m33(1, 1) = -21.616077216647483; m33(1, 2) = -81.41765352024458; m33(2, 0) = -1.826095889555674; m33(2, 1) = -72.65845670305742; m33(2, 2) = -68.43559055076148;
        Matrix<double, 4, 4, ColumnStorage> m44;
        m44(0, 0) = -0.6810804556978075; m44(0, 1) = -0.3717461369121619; m44(0, 2) = -0.6308202775413543; m44(0, 3) = 65.25663257716863; m44(1, 0) = -0.7322085856273928; m44(1, 1) = 0.3457881173505905; m44(1, 2) = 0.5867718168356992; m44(1, 3) = -72.82746334890318; m44(2, 0) = 0.0; m44(2, 1) = 0.8615308396047228; m44(2, 2) = -0.5077052416609282; m44(2, 3) = -50.951590508183095; m44(3, 0) = 0.0; m44(3, 1) = 0.0; m44(3, 2) = 0.0; m44(3, 3) = 1.0;
        Vector<double, 3> v = m33 * u;
        EXPECT_NEAR(v(0), -10054.998796995307, 1e-3);
        EXPECT_NEAR(v(1), 1916.692541293076, 1e-3);
        EXPECT_NEAR(v(2), -6614.459208480986, 1e-3);
        Vector<double, 3> w = m44 * u;
        EXPECT_NEAR(w(0), -45.129230953732865, 1e-3);
        EXPECT_NEAR(w(1), -92.57210823264782, 1e-3);
        EXPECT_NEAR(w(2), -51.21319422167322, 1e-3);
    }
    {
        Vector<double, 3> u; u(0) = 21.26102327359945; u(1) = -63.632003812151375; u(2) = 21.749676556807557;
        Matrix<double, 3, 3, ColumnStorage> m33;
        m33(0, 0) = -72.46587508929531; m33(0, 1) = 26.74413795088728; m33(0, 2) = 66.40155978725781; m33(1, 0) = 30.554539344667063; m33(1, 1) = 34.707428553000256; m33(1, 2) = -78.73225643210714; m33(2, 0) = 29.88002386285075; m33(2, 1) = 84.40743250651678; m33(2, 2) = 89.33116036637102;
        Matrix<double, 4, 4, ColumnStorage> m44;
        m44(0, 0) = -0.8156253236747055; m44(0, 1) = 0.0; m44(0, 2) = 0.5785804450381398; m44(0, 3) = -36.44550865077372; m44(1, 0) = -0.4886547783002762; m44(1, 1) = 0.5354367269787551; m44(1, 2) = -0.6888570381774334; m44(1, 3) = 26.24571159164249; m44(2, 0) = -0.3097932197851331; m44(2, 1) = -0.8445753438280554; m44(2, 2) = -0.43671575374937205; m44(2, 3) = 58.09911341828064; m44(3, 0) = 0.0; m44(3, 1) = 0.0; m44(3, 2) = 0.0; m44(3, 3) = 1.0;
        Vector<double, 3> v = m33 * u;
        EXPECT_NEAR(v(0), -1798.269296618443, 1e-3);
        EXPECT_NEAR(v(1), -3271.2835658593685, 1e-3);
        EXPECT_NEAR(v(2), -2792.8103398519333, 1e-3);
        Vector<double, 3> w = m44 * u;
        EXPECT_NEAR(w(0), -41.2026000982854, 1e-3);
        EXPECT_NEAR(w(1), -33.19691864907287, 1e-3);
        EXPECT_NEAR(w(2), 95.75618766923024, 1e-3);
    }
    {
        Vector<double, 3> u; u(0) = 60.89920397428068; u(1) = 41.99263628824008; u(2) = 2.8397811977043403;
        Matrix<double, 3, 3, ColumnStorage> m33;
        m33(0, 0) = 37.756340946034896; m33(0, 1) = 34.695096214807364; m33(0, 2) = 71.82563147437446; m33(1, 0) = 74.91874228380777; m33(1, 1) = 11.371524438829937; m33(1, 2) = 22.662878030996694; m33(2, 0) = -33.183135422386044; m33(2, 1) = -48.81237469926953; m33(2, 2) = 77.28631280104915;
        Matrix<double, 4, 4, ColumnStorage> m44;
        m44(0, 0) = 0.09071108841607319; m44(0, 1) = 0.15875989354807826; m44(0, 2) = -0.9831412892555038; m44(0, 3) = -20.28034322391565; m44(1, 0) = 0.8682639927927074; m44(1, 1) = -0.49610244790735053; m44(1, 2) = 0.0; m44(1, 3) = -65.91591665493024; m44(2, 0) = -0.487738800238444; m44(2, 1) = -0.8536261812883538; m44(2, 2) = -0.18284749208295362; m44(2, 3) = -28.08938117342305; m44(3, 0) = 0.0; m44(3, 1) = 0.0; m44(3, 2) = 0.0; m44(3, 3) = 1.0;
        Vector<double, 3> v = m33 * u;
        EXPECT_NEAR(v(0), 3960.2387427031354, 1e-3);
        EXPECT_NEAR(v(1), 5104.369672559086, 1e-3);
        EXPECT_NEAR(v(2), -3851.1106117734685, 1e-3);
        Vector<double, 3> w = m44 * u;
        EXPECT_NEAR(w(0), -10.881269828726621, 1e-3);
        EXPECT_NEAR(w(1), -33.87198031100272, 1e-3);
        EXPECT_NEAR(w(2), -94.15754648234068, 1e-3);
    }
    {
        Vector<double, 3> u; u(0) = 5.188467455547638; u(1) = 18.80451508110363; u(2) = -21.656727889135595;
        Matrix<double, 3, 3, ColumnStorage> m33;
        m33(0, 0) = 87.9457818703861; m33(0, 1) = 73.05874068692381; m33(0, 2) = 39.940251688233246; m33(1, 0) = 75.10237717805055; m33(1, 1) = 72.18489393198146; m33(1, 2) = 99.55423077694917; m33(2, 0) = -23.406607601688265; m33(2, 1) = 36.642344803461384; m33(2, 2) = 59.655570064408465;
        Matrix<double, 4, 4, ColumnStorage> m44;
        m44(0, 0) = 0.17067228277067772; m44(0, 1) = -0.6737686724718157; m44(0, 2) = -0.7189622715340582; m44(0, 3) = -18.335201879756795; m44(1, 0) = -0.9693828421385298; m44(1, 1) = -0.24555428191629253; m44(1, 2) = 0.0; m44(1, 3) = 29.845708309352148; m44(2, 0) = -0.1765442643114522; m44(2, 1) = 0.6969496901700587; m44(2, 2) = -0.6950491005034012; m44(2, 3) = -31.894796135192152; m44(3, 0) = 0.0; m44(3, 1) = 0.0; m44(3, 2) = 0.0; m44(3, 3) = 1.0;
        Vector<double, 3> v = m33 * u;
        EXPECT_NEAR(v(0), 965.1628555052321, 1e-3);
        EXPECT_NEAR(v(1), -408.9507197542059, 1e-3);
        EXPECT_NEAR(v(2), -724.347344278602, 1e-3);
        Vector<double, 3> w = m44 * u;
        EXPECT_NEAR(w(0), -14.549197180540624, 1e-3);
        EXPECT_NEAR(w(1), 20.198567783425617, 1e-3);
        EXPECT_NEAR(w(2), -4.652500106271733, 1e-3);
    }
    {
        Vector<double, 3> u; u(0) = 37.47349664672282; u(1) = -60.74514266373965; u(2) = -40.46471876369226;
        Matrix<double, 3, 3, ColumnStorage> m33;
        m33(0, 0) = -18.595918421354753; m33(0, 1) = -63.34971612297557; m33(0, 2) = 51.78315639657973; m33(1, 0) = 60.68830258069576; m33(1, 1) = -63.44130757420927; m33(1, 2) = 0.1227343978225548; m33(2, 0) = 69.14823593166864; m33(2, 1) = 49.45675635266332; m33(2, 2) = 81.76280170074256;
        Matrix<double, 4, 4, ColumnStorage> m44;
        m44(0, 0) = -0.250151688477416; m44(0, 1) = -0.2581347286114716; m44(0, 2) = -0.9331616122819132; m44(0, 3) = 15.902624286413442; m44(1, 0) = -0.7181238873483027; m44(1, 1) = 0.6959152839389017; m44(1, 2) = 0.0; m44(1, 3) = -12.358184768214258; m44(2, 0) = 0.6494014283720508; m44(2, 1) = 0.6701256445360971; m44(2, 2) = -0.3594570980846258; m44(2, 3) = -95.12313555065694; m44(3, 0) = 0.0; m44(3, 1) = 0.0; m44(3, 2) = 0.0; m44(3, 3) = 1.0;
        Vector<double, 3> v = m33 * u;
        EXPECT_NEAR(v(0), 1055.9425967082998, 1e-3);
        EXPECT_NEAR(v(1), 6122.987769732011, 1e-3);
        EXPECT_NEAR(v(2), -3721.5403091682488, 1e-3);
        Vector<double, 3> w = m44 * u;
        EXPECT_NEAR(w(0), 59.969118945113756, 1e-3);
        EXPECT_NEAR(w(1), -81.54227105743792, 1e-3);
        EXPECT_NEAR(w(2), -96.94934080054827, 1e-3);
    }
    {
        Vector<double, 3> u; u(0) = -18.242981648343232; u(1) = 56.21809486615007; u(2) = 68.42388784251386;
        Matrix<double, 3, 3, ColumnStorage> m33;
        m33(0, 0) = 72.32617531563085; m33(0, 1) = 16.921578626504015; m33(0, 2) = -27.34080821434577; m33(1, 0) = 75.31926937290174; m33(1, 1) = -81.4337114340538; m33(1, 2) = 49.516021610363254; m33(2, 0) = 40.651084247784695; m33(2, 1) = -19.615350275581747; m33(2, 2) = 96.61096468197644;
        Matrix<double, 4, 4, ColumnStorage> m44;
        m44(0, 0) = 0.8628144975844739; m44(0, 1) = -0.39604187777629973; m44(0, 2) = 0.3141687027784186; m44(0, 3) = -15.592609294116613; m44(1, 0) = -0.5055206650158348; m44(1, 1) = -0.6759578736217757; m44(1, 2) = 0.5362180622943203; m44(1, 3) = 10.857285613068584; m44(2, 0) = 0.0; m44(2, 1) = -0.6214754895699024; m44(2, 2) = -0.7834336065448367; m44(2, 3) = -95.09026591014809; m44(3, 0) = 0.0; m44(3, 1) = 0.0; m44(3, 2) = 0.0; m44(3, 3) = 1.0;
        Vector<double, 3> v = m33 * u;
        EXPECT_NEAR(v(0), -2238.9105712501655, 1e-3);
        EXPECT_NEAR(v(1), -2564.0174545638056, 1e-3);
        EXPECT_NEAR(v(2), 4766.163205213616, 1e-3);
        Vector<double, 3> w = m44 * u;
        EXPECT_NEAR(w(0), -32.1009941127323, 1e-3);
        EXPECT_NEAR(w(1), 18.768550516577342, 1e-3);
        EXPECT_NEAR(w(2), -183.634007166056, 1e-3);
    }
    {
        Vector<double, 3> u; u(0) = -59.27939726364973; u(1) = 96.57872437189005; u(2) = -50.71931054476342;
        Matrix<double, 3, 3, ColumnStorage> m33;
        m33(0, 0) = -53.51188145034156; m33(0, 1) = -39.773533094835024; m33(0, 2) = -81.91533947789623; m33(1, 0) = -48.98656664587613; m33(1, 1) = -30.35481126192778; m33(1, 2) = -24.41812181189411; m33(2, 0) = 50.987332445586844; m33(2, 1) = -15.614790715203753; m33(2, 2) = 1.953005394977339;
        Matrix<double, 4, 4, ColumnStorage> m44;
        m44(0, 0) = -0.9781360701985449; m44(0, 1) = 0.0; m44(0, 2) = -0.20796593032645286; m44(0, 3) = 31.079492180635583; m44(1, 0) = 0.17988007252338098; m44(1, 1) = 0.5018585333987298; m44(1, 2) = -0.846038516832336; m44(1, 3) = 89.95140758767585; m44(2, 0) = 0.10436947679053606; m44(2, 1) = -0.8649497167205016; m44(2, 2) = -0.4908859336542387; m44(2, 3) = -66.61272043333199; m44(3, 0) = 0.0; m44(3, 1) = 0.0; m44(3, 2) = 0.0; m44(3, 3) = 1.0;
        Vector<double, 3> v = m33 * u;
        EXPECT_NEAR(v(0), 3485.5645301169448, 1e-3);
        EXPECT_NEAR(v(1), 1210.7354976539962, 1e-3);
        EXPECT_NEAR(v(2), -4629.6099911875845, 1e-3);
        Vector<double, 3> w = m44 * u;
        EXPECT_NEAR(w(0), 99.6106974667983, 1e-3);
        EXPECT_NEAR(w(1), 170.66757254759543, 1e-3);
        EXPECT_NEAR(w(2), -131.4380242858257, 1e-3);
    }
    {
        Vector<double, 3> u; u(0) = -33.64679408120854; u(1) = 57.74585710730173; u(2) = -38.78100535372986;
        Matrix<double, 3, 3, ColumnStorage> m33;
        m33(0, 0) = -48.76215910463826; m33(0, 1) = -17.88064937607517; m33(0, 2) = -10.666843689761578; m33(1, 0) = 54.129177432133815; m33(1, 1) = -50.3409670888066; m33(1, 2) = 67.92880818886752; m33(2, 0) = -36.14536761608753; m33(2, 1) = 92.53432169416001; m33(2, 2) = 68.69365938680332;
        Matrix<double, 4, 4, ColumnStorage> m44;
        m44(0, 0) = 0.4751126538005149; m44(0, 1) = 0.5241078534879362; m44(0, 2) = -0.7068089728568109; m44(0, 3) = 7.21266115268952; m44(1, 0) = -0.7408885299063788; m44(1, 1) = 0.671628011813954; m44(1, 2) = 0.0; m44(1, 3) = -99.71424185587028; m44(2, 0) = 0.4747127051720829; m44(2, 1) = 0.5236666608245202; m44(2, 2) = 0.7074044641427562; m44(2, 3) = -11.057694172067542; m44(3, 0) = 0.0; m44(3, 1) = 0.0; m44(3, 2) = 0.0; m44(3, 3) = 1.0;
        Vector<double, 3> v = m33 * u;
        EXPECT_NEAR(v(0), 1021.8278247323343, 1e-3);
        EXPECT_NEAR(v(1), -7362.603053042772, 1e-3);
        EXPECT_NEAR(v(2), 3895.640286793221, 1e-3);
        Vector<double, 3> w = m44 * u;
        EXPECT_NEAR(w(0), 48.90246330164035, 1e-3);
        EXPECT_NEAR(w(1), -36.001982853511294, 1e-3);
        EXPECT_NEAR(w(2), -24.22453095406812, 1e-3);
    }
    {
        Vector<double, 3> u; u(0) = 15.282052878141911; u(1) = -2.259998336973041; u(2) = 91.08586195638563;
        Matrix<double, 3, 3, ColumnStorage> m33;
        m33(0, 0) = -88.39582865245377; m33(0, 1) = -64.58071224128338; m33(0, 2) = -4.675696828923634; m33(1, 0) = -86.1584474938788; m33(1, 1) = -34.543435043656444; m33(1, 2) = -62.495745585883625; m33(2, 0) = -59.183456388092395; m33(2, 1) = -57.14307430989349; m33(2, 2) = 44.41751022537625;
        Matrix<double, 4, 4, ColumnStorage> m44;
        m44(0, 0) = 0.44386855289227883; m44(0, 1) = 0.5950916475208161; m44(0, 2) = 0.6699601770286612; m44(0, 3) = -72.72024311218823; m44(1, 0) = -0.8960919080949868; m44(1, 1) = 0.29477162558569525; m44(1, 2) = 0.33185686823727545; m44(1, 3) = 15.476209050615353; m44(2, 0) = 0.0; m44(2, 1) = -0.7476467212531113; m44(2, 2) = 0.6640966648007446; m44(2, 3) = 99.8824543736142; m44(3, 0) = 0.0; m44(3, 1) = 0.0; m44(3, 2) = 0.0; m44(3, 3) = 1.0;
        Vector<double, 3> v = m33 * u;
        EXPECT_NEAR(v(0), -1630.8073013173835, 1e-3);
        EXPECT_NEAR(v(1), -6931.088700045262, 1e-3);
        EXPECT_NEAR(v(2), 3270.5057477108207, 1e-3);
        Vector<double, 3> w = m44 * u;
        EXPECT_NEAR(w(0), -6.258026348579378, 1e-3);
        EXPECT_NEAR(w(1), 31.343370633361097, 1e-3);
        EXPECT_NEAR(w(2), 162.0619518560263, 1e-3);
    }
    {
        Vector<double, 3> u; u(0) = -45.53190350588168; u(1) = 98.18869916536156; u(2) = 82.67567746788419;
        Matrix<double, 3, 3, ColumnStorage> m33;
        m33(0, 0) = -51.297389589304785; m33(0, 1) = -66.59459838930371; m33(0, 2) = 27.56805306675247; m33(1, 0) = 96.48690448918092; m33(1, 1) = 9.10079699500524; m33(1, 2) = 83.88487605323621; m33(2, 0) = 41.570057963368356; m33(2, 1) = 24.52264689542079; m33(2, 2) = -76.44940699537474;
        Matrix<double, 4, 4, ColumnStorage> m44;
        m44(0, 0) = 0.774965095070908; m44(0, 1) = 0.0; m44(0, 2) = 0.6320040359220331; m44(0, 3) = -94.78257020280037; m44(1, 0) = -0.5539089538751407; m44(1, 1) = 0.4815245168438749; m44(1, 2) = 0.6792046896887713; m44(1, 3) = 63.87247725975902; m44(2, 0) = -0.3043254380407359; m44(2, 1) = -0.8764326212997054; m44(2, 2) = 0.3731646929748866; m44(2, 3) = -95.85244653326883; m44(3, 0) = 0.0; m44(3, 1) = 0.0; m44(3, 2) = 0.0; m44(3, 3) = 1.0;
        Vector<double, 3> v = m33 * u;
        EXPECT_NEAR(v(0), -1923.9617306372265, 1e-3);
        EXPECT_NEAR(v(1), 3435.6019505357876, 1e-3);
        EXPECT_NEAR(v(2), -5805.423584529748, 1e-3);
        Vector<double, 3> w = m44 * u;
        EXPECT_NEAR(w(0), -77.81684429970429, 1e-3);
        EXPECT_NEAR(w(1), 192.52698008315832, 1e-3);
        EXPECT_NEAR(w(2), -137.20006524672516, 1e-3);
    }
}

/**
 * Opérateurs d'arithmétique vectorielle
 */
TEST(TestLabo1b, VectorOperators)
{
    Vector<double> v(5);
    v(0) = 0.1;
    v(1) = 0.2;
    v(2) = 0.4;
    v(3) = 0.8;
    v(4) = 1.6;

    // Test : multiplication  scalaire * vecteur
    const double alpha = 4.0;
    const auto v2 = alpha * v;
    EXPECT_DOUBLE_EQ(v2(0), alpha * v(0));
    EXPECT_DOUBLE_EQ(v2(1), alpha * v(1));
    EXPECT_DOUBLE_EQ(v2(2), alpha * v(2));
    EXPECT_DOUBLE_EQ(v2(3), alpha * v(3));
    EXPECT_DOUBLE_EQ(v2(4), alpha * v(4));

    // Test : addition  vecteur + vecteur
    const auto v3 = v + v2;
    EXPECT_DOUBLE_EQ(v3(0), v(0) + v2(0));
    EXPECT_DOUBLE_EQ(v3(1), v(1) + v2(1));
    EXPECT_DOUBLE_EQ(v3(2), v(2) + v2(2));
    EXPECT_DOUBLE_EQ(v3(3), v(3) + v2(3));
    EXPECT_DOUBLE_EQ(v3(4), v(4) + v2(4));
}

/**
 * Mathématiques 3D
 */
TEST(TestLabo1b, Math3D)
{
    // Test : norme d'un vecteur de dimension 3
    Vector3d v;
    v.setZero();
    v(1) = 2.0;
    EXPECT_EQ(v.rows(), 3);
    EXPECT_EQ(v.cols(), 1);
    EXPECT_DOUBLE_EQ(v(0), 0.0);
    EXPECT_DOUBLE_EQ(v(1), 2.0);
    EXPECT_DOUBLE_EQ(v(2), 0.0);
    EXPECT_DOUBLE_EQ(v.norm(), 2.0);

    // Test : calcul de la norme d'un deuxième vecteur 3D
    Vector3d v2;
    v2(0) = 4.0;
    v2(1) = 2.0;
    v2(2) = 5.0;
    EXPECT_EQ(v2.rows(), 3);
    EXPECT_EQ(v2.cols(), 1);
    EXPECT_DOUBLE_EQ(v2(0), 4.0);
    EXPECT_DOUBLE_EQ(v2(1), 2.0);
    EXPECT_DOUBLE_EQ(v2(2), 5.0);
    EXPECT_DOUBLE_EQ(v2.norm(), 6.7082039324993690892275210061938);

    // Test : produit scalaire 
    EXPECT_DOUBLE_EQ(v.dot(v2), 4.0);

    // Test : matrice identité 4x4
    Matrix4d M;
    M.setIdentity();
    EXPECT_DOUBLE_EQ(M(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(M(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(M(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(M(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(M(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(M(1, 2), 0.0);
    EXPECT_DOUBLE_EQ(M(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(M(2, 1), 0.0);
    EXPECT_DOUBLE_EQ(M(2, 2), 1.0);

    // Test : création d'une matrice de rotation de 45 degrés autour de l'axe des x
    const auto Rx = makeRotation<double>(M_PI / 4.0, 0, 0);
    EXPECT_NEAR(Rx(0, 0), 1, 1e-3); EXPECT_NEAR(Rx(0, 1), 0, 1e-3); EXPECT_NEAR(Rx(0, 2), 0, 1e-3);
    EXPECT_NEAR(Rx(1, 0), 0, 1e-3); EXPECT_NEAR(Rx(1, 1), 0.7071, 1e-3); EXPECT_NEAR(Rx(1, 2), -0.7071, 1e-3);
    EXPECT_NEAR(Rx(2, 0), 0, 1e-3); EXPECT_NEAR(Rx(2, 1), 0.7071, 1e-3); EXPECT_NEAR(Rx(2, 2), 0.7071, 1e-3);

    // Test : création d'une matrice de rotation de 45 degrés autour de l'axe des y
    const auto Ry = makeRotation<double>(0, M_PI / 4.0, 0);
    EXPECT_NEAR(Ry(0, 0), 0.7071, 1e-3); EXPECT_NEAR(Ry(0, 1), 0, 1e-3); EXPECT_NEAR(Ry(0, 2), 0.7071, 1e-3);
    EXPECT_NEAR(Ry(1, 0), 0, 1e-3); EXPECT_NEAR(Ry(1, 1), 1, 1e-3); EXPECT_NEAR(Ry(1, 2), 0, 1e-3);
    EXPECT_NEAR(Ry(2, 0), -0.7071, 1e-3); EXPECT_NEAR(Ry(2, 1), 0, 1e-3); EXPECT_NEAR(Ry(2, 2), 0.7071, 1e-3);

    // Test : création d'une matrice de rotation de 45 degrés autour de l'axe des z
    const auto Rz = makeRotation<double>(0, 0, M_PI / 4.0);
    EXPECT_NEAR(Rz(0, 0), 0.7071, 1e-3); EXPECT_NEAR(Rz(0, 1), -0.7071, 1e-3); EXPECT_NEAR(Rz(0, 2), 0, 1e-3);
    EXPECT_NEAR(Rz(1, 0), 0.7071, 1e-3); EXPECT_NEAR(Rz(1, 1), 0.7071, 1e-3); EXPECT_NEAR(Rz(1, 2), 0, 1e-3);
    EXPECT_NEAR(Rz(2, 0), 0, 1e-3); EXPECT_NEAR(Rz(2, 1), 0, 1e-3); EXPECT_NEAR(Rz(2, 2), 1, 1e-3);

    // Test : création d'une matrice de rotation quelconque.
    const auto Rxyz = makeRotation<double>(M_PI / 3.0, -M_PI / 6.0, M_PI / 4.0);
    EXPECT_NEAR(Rxyz(0, 0), 0.6124, 1e-3); EXPECT_NEAR(Rxyz(0, 1), -0.6597, 1e-3); EXPECT_NEAR(Rxyz(0, 2), 0.4356, 1e-3);
    EXPECT_NEAR(Rxyz(1, 0), 0.6124, 1e-3); EXPECT_NEAR(Rxyz(1, 1), 0.0474, 1e-3); EXPECT_NEAR(Rxyz(1, 2), -0.7891, 1e-3);
    EXPECT_NEAR(Rxyz(2, 0), 0.5, 1e-3); EXPECT_NEAR(Rxyz(2, 1), 0.75, 1e-3); EXPECT_NEAR(Rxyz(2, 2), 0.4330, 1e-3);

    // Test : création d'une transformation homogène via la sous-matrice 3x3 en
    // utilisant la fonction `block`
    M.block(0, 0, 3, 3) = Rxyz;
    M(0, 3) = -0.1;
    M(1, 3) = 1.0;
    M(2, 3) = 2.1;

    // Test : calcule l'inverse de la matrice M et vérifie que M^(-1) * M * v = v
    const Matrix4d Minv = M.inverse();
    const Vector3d v3 = Minv * (M * v2);
    EXPECT_DOUBLE_EQ(v3(0), v2(0));
    EXPECT_DOUBLE_EQ(v3(1), v2(1));
    EXPECT_DOUBLE_EQ(v3(2), v2(2));

    // Test : translation d'un vecteur 3D effectuée avec une matrice 4x4 en coordonnées homogènes
    Matrix4d T;
    T.setIdentity();
    T(0, 3) = 1.2;
    T(1, 3) = 2.5;
    T(2, 3) = -4.0;
    const Vector3d t = T * v3;
    EXPECT_DOUBLE_EQ(t(0), v3(0) + 1.2);
    EXPECT_DOUBLE_EQ(t(1), v3(1) + 2.5);
    EXPECT_DOUBLE_EQ(t(2), v3(2) - 4.0);

    // Test : inverse d'un matrice de rotation
    const Matrix3d Rinv = Rxyz.inverse();
    const Matrix3d RT = Rxyz.transpose<double, 3, 3, ColumnStorage>();
    EXPECT_DOUBLE_EQ(Rinv(0, 0), RT(0, 0));
    EXPECT_DOUBLE_EQ(Rinv(1, 1), RT(1, 1));
    EXPECT_DOUBLE_EQ(Rinv(0, 2), RT(0, 2));
}

TEST(TestLabo1b, RigidTransformInverse) {
    {
        Matrix<double, 4, 4, ColumnStorage> m;
        m(0, 0) = -0.7899179842354741; m(0, 1) = 0.5558875173441905; m(0, 2) = 0.2588795979645317; m(0, 3) = 63.53646575655944; m(1, 0) = -0.6132125065435026; m(1, 1) = -0.7160740240561843; m(1, 2) = -0.3334792555626483; m(1, 3) = 46.40415080808421; m(2, 0) = 0.0; m(2, 1) = -0.4221694684991984; m(2, 2) = 0.9065169275127214; m(2, 3) = 37.91755407005138; m(3, 0) = 0.0; m(3, 1) = 0.0; m(3, 2) = 0.0; m(3, 3) = 1.0;
        Matrix<double, 4, 4, ColumnStorage> mi = m.inverse();
        EXPECT_NEAR(mi(0, 0), -0.7899179842354742, 1e-3);
        EXPECT_NEAR(mi(0, 1), -0.6132125065435026, 1e-3);
        EXPECT_NEAR(mi(0, 2), 5.0321798256853854e-17, 1e-3);
        EXPECT_NEAR(mi(0, 3), 78.64420258691568, 1e-3);
        EXPECT_NEAR(mi(1, 0), 0.5558875173441905, 1e-3);
        EXPECT_NEAR(mi(1, 1), -0.7160740240561843, 1e-3);
        EXPECT_NEAR(mi(1, 2), -0.42216946849919834, 1e-3);
        EXPECT_NEAR(mi(1, 3), 13.917312440360101, 1e-3);
        EXPECT_NEAR(mi(2, 0), 0.25887959796453164, 1e-3);
        EXPECT_NEAR(mi(2, 1), -0.33347925556264824, 1e-3);
        EXPECT_NEAR(mi(2, 2), 0.9065169275127212, 1e-3);
        EXPECT_NEAR(mi(2, 3), -35.34637765902901, 1e-3);
        EXPECT_NEAR(mi(3, 0), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 1), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 2), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 3), 1.0, 1e-3);
    }
    {
        Matrix<double, 4, 4, ColumnStorage> m;
        m(0, 0) = 0.6589112893907124; m(0, 1) = -0.32720889915223983; m(0, 2) = 0.6773258071482646; m(0, 3) = -47.579657162583366; m(1, 0) = -0.7522206542720485; m(1, 1) = -0.2866202043457096; m(1, 2) = 0.5933067888937011; m(1, 3) = 30.029355381737275; m(2, 0) = 0.0; m(2, 1) = -0.9004350030826229; m(2, 2) = -0.4349905806147956; m(2, 3) = 80.66547089822612; m(3, 0) = 0.0; m(3, 1) = 0.0; m(3, 2) = 0.0; m(3, 3) = 1.0;
        Matrix<double, 4, 4, ColumnStorage> mi = m.inverse();
        EXPECT_NEAR(mi(0, 0), 0.6589112893907124, 1e-3);
        EXPECT_NEAR(mi(0, 1), -0.7522206542720486, 1e-3);
        EXPECT_NEAR(mi(0, 2), 5.551115123125783e-17, 1e-3);
        EXPECT_NEAR(mi(0, 3), 53.939474602384124, 1e-3);
        EXPECT_NEAR(mi(1, 0), -0.3272088991522398, 1e-3);
        EXPECT_NEAR(mi(1, 1), -0.28662020434570956, 1e-3);
        EXPECT_NEAR(mi(1, 2), -0.9004350030826228, 1e-3);
        EXPECT_NEAR(mi(1, 3), 65.67254627057903, 1e-3);
        EXPECT_NEAR(mi(2, 0), 0.6773258071482645, 1e-3);
        EXPECT_NEAR(mi(2, 1), 0.5933067888937011, 1e-3);
        EXPECT_NEAR(mi(2, 2), -0.43499058061479545, 1e-3);
        EXPECT_NEAR(mi(2, 3), 49.49902929898342, 1e-3);
        EXPECT_NEAR(mi(3, 0), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 1), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 2), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 3), 1.0, 1e-3);
    }
    {
        Matrix<double, 4, 4, ColumnStorage> m;
        m(0, 0) = -0.12921677297502387; m(0, 1) = 0.0; m(0, 2) = -0.9916163701663668; m(0, 3) = -85.67624239250289; m(1, 0) = 0.5536164401240612; m(1, 1) = -0.8296411663169171; m(1, 2) = -0.07214133611645579; m(1, 3) = -80.61895745285052; m(2, 0) = -0.8226857618837724; m(2, 1) = -0.5582969954713223; m(2, 2) = 0.1072035542387071; m(2, 3) = -59.16846832515195; m(3, 0) = 0.0; m(3, 1) = 0.0; m(3, 2) = 0.0; m(3, 3) = 1.0;
        Matrix<double, 4, 4, ColumnStorage> mi = m.inverse();
        EXPECT_NEAR(mi(0, 0), -0.12921677297502387, 1e-3);
        EXPECT_NEAR(mi(0, 1), 0.5536164401240611, 1e-3);
        EXPECT_NEAR(mi(0, 2), -0.8226857618837723, 1e-3);
        EXPECT_NEAR(mi(0, 3), -15.115883774598393, 1e-3);
        EXPECT_NEAR(mi(1, 0), -0.0, 1e-3);
        EXPECT_NEAR(mi(1, 1), -0.829641166316917, 1e-3);
        EXPECT_NEAR(mi(1, 2), -0.5582969954713222, 1e-3);
        EXPECT_NEAR(mi(1, 3), -99.91838398100924, 1e-3);
        EXPECT_NEAR(mi(2, 0), -0.991616370166367, 1e-3);
        EXPECT_NEAR(mi(2, 1), -0.07214133611645579, 1e-3);
        EXPECT_NEAR(mi(2, 2), 0.10720355423870709, 1e-3);
        EXPECT_NEAR(mi(2, 3), -84.43085369439521, 1e-3);
        EXPECT_NEAR(mi(3, 0), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 1), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 2), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 3), 1.0, 1e-3);
    }
    {
        Matrix<double, 4, 4, ColumnStorage> m;
        m(0, 0) = 0.6178488588834932; m(0, 1) = 0.656563779114608; m(0, 2) = 0.4326508887441579; m(0, 3) = -3.573389719025613; m(1, 0) = -0.7862968825935693; m(1, 1) = 0.5159084191865928; m(1, 2) = 0.33996428553015007; m(1, 3) = 9.98340634305292; m(2, 0) = 0.0; m(2, 1) = -0.5502385909468138; m(2, 2) = 0.8350074808244923; m(2, 3) = 57.03623931292819; m(3, 0) = 0.0; m(3, 1) = 0.0; m(3, 2) = 0.0; m(3, 3) = 1.0;
        Matrix<double, 4, 4, ColumnStorage> mi = m.inverse();
        EXPECT_NEAR(mi(0, 0), 0.6178488588834931, 1e-3);
        EXPECT_NEAR(mi(0, 1), -0.7862968825935694, 1e-3);
        EXPECT_NEAR(mi(0, 2), -0.0, 1e-3);
        EXPECT_NEAR(mi(0, 3), 10.057736045453359, 1e-3);
        EXPECT_NEAR(mi(1, 0), 0.656563779114608, 1e-3);
        EXPECT_NEAR(mi(1, 1), 0.5159084191865929, 1e-3);
        EXPECT_NEAR(mi(1, 2), -0.5502385909468138, 1e-3);
        EXPECT_NEAR(mi(1, 3), 28.57917482608178, 1e-3);
        EXPECT_NEAR(mi(2, 0), 0.4326508887441579, 1e-3);
        EXPECT_NEAR(mi(2, 1), 0.33996428553015007, 1e-3);
        EXPECT_NEAR(mi(2, 2), 0.8350074808244924, 1e-3);
        EXPECT_NEAR(mi(2, 3), -49.473657871198526, 1e-3);
        EXPECT_NEAR(mi(3, 0), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 1), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 2), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 3), 1.0, 1e-3);
    }
    {
        Matrix<double, 4, 4, ColumnStorage> m;
        m(0, 0) = -0.30532791727234754; m(0, 1) = 0.42993321248800975; m(0, 2) = -0.8496659906892062; m(0, 3) = -4.900019523419871; m(1, 0) = -0.9522472698486095; m(1, 1) = -0.13785349298616878; m(1, 2) = 0.272436325656807; m(1, 3) = 61.254180348498465; m(2, 0) = 0.0; m(2, 1) = 0.8922745358191345; m(2, 2) = 0.45149324771113447; m(2, 3) = 41.98722428204593; m(3, 0) = 0.0; m(3, 1) = 0.0; m(3, 2) = 0.0; m(3, 3) = 1.0;
        Matrix<double, 4, 4, ColumnStorage> mi = m.inverse();
        EXPECT_NEAR(mi(0, 0), -0.3053279172723474, 1e-3);
        EXPECT_NEAR(mi(0, 1), -0.9522472698486094, 1e-3);
        EXPECT_NEAR(mi(0, 2), 0.0, 1e-3);
        EXPECT_NEAR(mi(0, 3), 56.83301324799237, 1e-3);
        EXPECT_NEAR(mi(1, 0), 0.42993321248800975, 1e-3);
        EXPECT_NEAR(mi(1, 1), -0.13785349298616878, 1e-3);
        EXPECT_NEAR(mi(1, 2), 0.8922745358191346, 1e-3);
        EXPECT_NEAR(mi(1, 3), -26.91334720059331, 1e-3);
        EXPECT_NEAR(mi(2, 0), -0.8496659906892063, 1e-3);
        EXPECT_NEAR(mi(2, 1), 0.27243632565680703, 1e-3);
        EXPECT_NEAR(mi(2, 2), 0.4514932477111345, 1e-3);
        EXPECT_NEAR(mi(2, 3), -39.80819202150404, 1e-3);
        EXPECT_NEAR(mi(3, 0), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 1), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 2), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 3), 1.0, 1e-3);
    }
    {
        Matrix<double, 4, 4, ColumnStorage> m;
        m(0, 0) = -0.8426665329724763; m(0, 1) = -0.13329457405486883; m(0, 2) = 0.5216758291656588; m(0, 3) = 54.75958681320196; m(1, 0) = -0.5384358032376251; m(1, 1) = 0.20860959822407732; m(1, 2) = -0.81643672217793; m(1, 3) = 91.96965943681118; m(2, 0) = 0.0; m(2, 1) = -0.9688728461755547; m(2, 2) = -0.24755889792871488; m(2, 3) = 99.00199341303772; m(3, 0) = 0.0; m(3, 1) = 0.0; m(3, 2) = 0.0; m(3, 3) = 1.0;
        Matrix<double, 4, 4, ColumnStorage> mi = m.inverse();
        EXPECT_NEAR(mi(0, 0), -0.8426665329724762, 1e-3);
        EXPECT_NEAR(mi(0, 1), -0.5384358032376252, 1e-3);
        EXPECT_NEAR(mi(0, 2), -2.7755575615628914e-17, 1e-3);
        EXPECT_NEAR(mi(0, 3), 95.66382861923651, 1e-3);
        EXPECT_NEAR(mi(1, 0), -0.13329457405486886, 1e-3);
        EXPECT_NEAR(mi(1, 1), 0.20860959822407737, 1e-3);
        EXPECT_NEAR(mi(1, 2), -0.9688728461755547, 1e-3);
        EXPECT_NEAR(mi(1, 3), 84.03374523091135, 1e-3);
        EXPECT_NEAR(mi(2, 0), 0.5216758291656588, 1e-3);
        EXPECT_NEAR(mi(2, 1), -0.8164367221779302, 1e-3);
        EXPECT_NEAR(mi(2, 2), -0.2475588979287149, 1e-3);
        EXPECT_NEAR(mi(2, 3), 71.02947881694217, 1e-3);
        EXPECT_NEAR(mi(3, 0), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 1), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 2), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 3), 1.0, 1e-3);
    }
    {
        Matrix<double, 4, 4, ColumnStorage> m;
        m(0, 0) = -0.7135830808790772; m(0, 1) = 0.0; m(0, 2) = -0.7005706150582712; m(0, 3) = 46.98894352455778; m(1, 0) = 0.5815768734392242; m(1, 1) = 0.5575439866421585; m(1, 2) = -0.5923791380862643; m(1, 3) = -32.31315706340996; m(2, 0) = 0.39059893364393755; m(2, 1) = -0.8301473983330724; m(2, 2) = -0.39785395571371457; m(2, 3) = 7.702713412699751; m(3, 0) = 0.0; m(3, 1) = 0.0; m(3, 2) = 0.0; m(3, 3) = 1.0;
        Matrix<double, 4, 4, ColumnStorage> mi = m.inverse();
        EXPECT_NEAR(mi(0, 0), -0.7135830808790772, 1e-3);
        EXPECT_NEAR(mi(0, 1), 0.581576873439224, 1e-3);
        EXPECT_NEAR(mi(0, 2), 0.39059893364393744, 1e-3);
        EXPECT_NEAR(mi(0, 3), 49.314428298230084, 1e-3);
        EXPECT_NEAR(mi(1, 0), -5.551115123125783e-17, 1e-3);
        EXPECT_NEAR(mi(1, 1), 0.5575439866421584, 1e-3);
        EXPECT_NEAR(mi(1, 2), -0.8301473983330723, 1e-3);
        EXPECT_NEAR(mi(1, 3), 24.41039390978577, 1e-3);
        EXPECT_NEAR(mi(2, 0), -0.7005706150582711, 1e-3);
        EXPECT_NEAR(mi(2, 1), -0.5923791380862641, 1e-3);
        EXPECT_NEAR(mi(2, 2), -0.3978539557137144, 1e-3);
        EXPECT_NEAR(mi(2, 3), 16.841987936840617, 1e-3);
        EXPECT_NEAR(mi(3, 0), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 1), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 2), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 3), 1.0, 1e-3);
    }
    {
        Matrix<double, 4, 4, ColumnStorage> m;
        m(0, 0) = 0.655199358191179; m(0, 1) = -0.7015590749923188; m(0, 2) = 0.2802296653136303; m(0, 3) = -24.378695814492147; m(1, 0) = -0.7554560219005917; m(1, 1) = -0.6084550818878111; m(1, 2) = 0.24304035117451178; m(1, 3) = 36.943036725120294; m(2, 0) = 0.0; m(2, 1) = -0.3709410702804682; m(2, 2) = -0.928656407063011; m(2, 3) = -32.64950630748629; m(3, 0) = 0.0; m(3, 1) = 0.0; m(3, 2) = 0.0; m(3, 3) = 1.0;
        Matrix<double, 4, 4, ColumnStorage> mi = m.inverse();
        EXPECT_NEAR(mi(0, 0), 0.6551993581911789, 1e-3);
        EXPECT_NEAR(mi(0, 1), -0.7554560219005915, 1e-3);
        EXPECT_NEAR(mi(0, 2), 5.1550786254351327e-17, 1e-3);
        EXPECT_NEAR(mi(0, 3), 43.88174541248007, 1e-3);
        EXPECT_NEAR(mi(1, 0), -0.7015590749923187, 1e-3);
        EXPECT_NEAR(mi(1, 1), -0.6084550818878111, 1e-3);
        EXPECT_NEAR(mi(1, 2), -0.37094107028046824, 1e-3);
        EXPECT_NEAR(mi(1, 3), -6.735959663194609, 1e-3);
        EXPECT_NEAR(mi(2, 0), 0.28022966531363025, 1e-3);
        EXPECT_NEAR(mi(2, 1), 0.24304035117451175, 1e-3);
        EXPECT_NEAR(mi(2, 2), -0.9286564070630109, 1e-3);
        EXPECT_NEAR(mi(2, 3), -32.467188070139514, 1e-3);
        EXPECT_NEAR(mi(3, 0), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 1), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 2), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 3), 1.0, 1e-3);
    }
    {
        Matrix<double, 4, 4, ColumnStorage> m;
        m(0, 0) = -0.12316316348384228; m(0, 1) = 0.47851608243384774; m(0, 2) = -0.8693981792095123; m(0, 3) = -97.68194268706864; m(1, 0) = -0.9684362001788849; m(1, 1) = -0.24926156178416892; m(1, 2) = 0.0; m(1, 3) = -93.30408292578261; m(2, 0) = -0.2167075479620758; m(2, 1) = 0.8419566691161013; m(2, 2) = 0.4941121390809829; m(2, 3) = -48.17348057895294; m(3, 0) = 0.0; m(3, 1) = 0.0; m(3, 2) = 0.0; m(3, 3) = 1.0;
        Matrix<double, 4, 4, ColumnStorage> mi = m.inverse();
        EXPECT_NEAR(mi(0, 0), -0.12316316348384225, 1e-3);
        EXPECT_NEAR(mi(0, 1), -0.9684362001788847, 1e-3);
        EXPECT_NEAR(mi(0, 2), -0.2167075479620758, 1e-3);
        EXPECT_NEAR(mi(0, 3), -112.82942545947078, 1e-3);
        EXPECT_NEAR(mi(1, 0), 0.47851608243384763, 1e-3);
        EXPECT_NEAR(mi(1, 1), -0.2492615617841688, 1e-3);
        EXPECT_NEAR(mi(1, 2), 0.8419566691161011, 1e-3);
        EXPECT_NEAR(mi(1, 3), 64.04524235620795, 1e-3);
        EXPECT_NEAR(mi(2, 0), -0.8693981792095123, 1e-3);
        EXPECT_NEAR(mi(2, 1), 0.0, 1e-3);
        EXPECT_NEAR(mi(2, 2), 0.49411213908098284, 1e-3);
        EXPECT_NEAR(mi(2, 3), -61.12140157794279, 1e-3);
        EXPECT_NEAR(mi(3, 0), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 1), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 2), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 3), 1.0, 1e-3);
    }
    {
        Matrix<double, 4, 4, ColumnStorage> m;
        m(0, 0) = -0.8820615884895145; m(0, 1) = -0.11346792289034995; m(0, 2) = -0.45726620757093356; m(0, 3) = -48.8253349197771; m(1, 0) = -0.12758811649101956; m(1, 1) = 0.9918272392560481; m(1, 2) = 0.0; m(1, 3) = 50.033634014775515; m(2, 0) = 0.4535290802601621; m(2, 1) = 0.058341734158967; m(2, 2) = -0.8893298687290875; m(2, 3) = 92.22065810812657; m(3, 0) = 0.0; m(3, 1) = 0.0; m(3, 2) = 0.0; m(3, 3) = 1.0;
        Matrix<double, 4, 4, ColumnStorage> mi = m.inverse();
        EXPECT_NEAR(mi(0, 0), -0.8820615884895142, 1e-3);
        EXPECT_NEAR(mi(0, 1), -0.12758811649101953, 1e-3);
        EXPECT_NEAR(mi(0, 2), 0.453529080260162, 1e-3);
        EXPECT_NEAR(mi(0, 3), -78.50800560549041, 1e-3);
        EXPECT_NEAR(mi(1, 0), -0.11346792289034992, 1e-3);
        EXPECT_NEAR(mi(1, 1), 0.991827239256048, 1e-3);
        EXPECT_NEAR(mi(1, 2), 0.05834173415896697, 1e-3);
        EXPECT_NEAR(mi(1, 3), -60.545143551904374, 1e-3);
        EXPECT_NEAR(mi(2, 0), -0.45726620757093345, 1e-3);
        EXPECT_NEAR(mi(2, 1), 6.120531779240618e-18, 1e-3);
        EXPECT_NEAR(mi(2, 2), -0.8893298687290875, 1e-3);
        EXPECT_NEAR(mi(2, 3), 59.68841003726312, 1e-3);
        EXPECT_NEAR(mi(3, 0), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 1), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 2), 0.0, 1e-3);
        EXPECT_NEAR(mi(3, 3), 1.0, 1e-3);
    }
}

/**
 * Test des performance de la multiplication  matrice * vecteur
 * pour de grandes dimensions.
 */
TEST(TestLabo1b, PerformanceMatrixVector)
{
    Matrix<double> A(16384, 16384);     // grande matrice avec stockage colonne
    Vector<double> v(16384);            // grand vecteur

    using namespace std::chrono;
    // Test : multiplication avec l'algorithme naif.
    high_resolution_clock::time_point t = high_resolution_clock::now();
    naiveMatrixMult(A, v);
    const duration<double> naive_t = duration_cast<duration<double>>(high_resolution_clock::now() - t);

    // Test : multiplication avec l'implémentation spécifique pour les matrices avec
    // stockage par colonnes.
    t = high_resolution_clock::now();
    A* v;
    const duration<double> optimal_t = duration_cast<duration<double>>(high_resolution_clock::now() - t);

    EXPECT_TRUE(optimal_t < 0.4 * naive_t)
        << "Naive time: " << duration_cast<std::chrono::milliseconds>(naive_t).count() << " ms, "
        << "optimized time: " << duration_cast<std::chrono::milliseconds>(optimal_t).count() << " ms";
}

/**
 * Test des performances de l'addition  matrice + matrice
 * pour de grandes dimensions.
 */
TEST(TestLabo1b, PerformanceLargeMatrixMatrix)
{
    // deux grandes matrices à stockage par colonnes
    Matrix<double> A(16384, 16384);
    Matrix<double> B(16384, 16384);

    using namespace std::chrono;
    high_resolution_clock::time_point t = high_resolution_clock::now();
    // Test : addition avec l'algorithme naif 
    naiveMatrixAddition(A, B);
    const duration<double> naive_t = duration_cast<duration<double>>(high_resolution_clock::now() - t);

    // Test : addition avec l'implémentation spécifique pour les matrices à
    // stockage par colonnes.
    t = high_resolution_clock::now();
    A + B;
    const duration<double> optimal_t = duration_cast<duration<double>>(high_resolution_clock::now() - t);

    EXPECT_TRUE(optimal_t < 0.4 * naive_t);
}

/**
 * Test pour la matrice à taille fixe 4D
 */
TEST(TestLabo1b, Matrix4x4SizeTest)
{
    Matrix4d M;

    EXPECT_EQ(M.cols(), 4);
    EXPECT_EQ(M.rows(), 4);
}

/**
 * Test pour la matrice à taille fixe 3D
 */
TEST(TestLabo1b, Matrix3x3SizeTest)
{
    Matrix3d M;

    EXPECT_EQ(M.cols(), 3);
    EXPECT_EQ(M.rows(), 3);
}
