#include "ParticleSystem.h"

using namespace gti320;

/**
 * Calcule des forces qui affectent chacune des particules.
 *
 * Étant donné une particule p, la force est stockée dans le vecteur p.f
 * Les forces prisent en compte sont : la gravité et la force des ressorts.
 */
void ParticleSystem::computeForces()
{
    // TODO 
    //
    // Calcul de la force gravitationnelle sur chacune des particules

    float fGrav = -9.81f;

    for (Particle& p : m_particles)
    {
        p.f(0) = 0;
        p.f(1) = p.m * fGrav;
    }

    // TODO
    //
    // Pour chaque ressort, ajouter la force exercée à chacune des exptrémités sur
    // les particules appropriées. Pour un ressort s, les deux particules
    // auxquelles le ressort est attaché sont m_particles[s.index0] et
    // m_particles[s.index1]. On rappelle que les deux forces sont de même
    // magnitude mais dans des directions opposées.
    for (const Spring& s : m_springs)
    {
        Particle& p1 = m_particles[s.index0];
        Particle& p2 = m_particles[s.index1];

        Vector<float, 2> diffPosParticles = p2.x - p1.x;
        float alpha = this->CalculAlpha(diffPosParticles, s);
        Vector<float, 2> f1 = alpha * diffPosParticles;

        p1.f = p1.f + f1;
        p2.f = p2.f - f1;
    }
}

/**
 * Assemble les données du système dans les vecteurs trois vecteurs d'état _pos,
 * _vel et _force.
 */
void ParticleSystem::pack(Vector<float, Dynamic>& _pos,
    Vector<float, Dynamic>& _vel,
    Vector<float, Dynamic>& _force)
{
    // TODO 
    //
    // Copier les données des particules dans les vecteurs. Attention, si on a
    // changé de modèle, il est possible que les vecteurs n'aient pas la bonne
    // taille. Rappel : la taille de ces vecteurs doit être 2 fois le nombre de
    // particules.

    int numParticles = m_particles.size() * 2;

    _pos.resize(numParticles);
    _pos.setZero();

    _vel.resize(numParticles);
    _vel.setZero();

    _force.resize(numParticles);
    _force.setZero();
    
    for (int i = 0; i < m_particles.size(); i++) 
    {
        Particle& p = m_particles[i];
        int pos = i * 2;

        _pos(pos) = p.x(0);
        _pos(pos + 1) = p.x(1);

        _vel(pos) = p.v(0);
        _vel(pos + 1) = p.v(1);

        _force(pos) = p.f(0);
        _force(pos + 1) = p.f(1);
    }
}

/**
 * Copie les données des vecteurs d'états dans le particules du système.
 */
void ParticleSystem::unpack(const Vector<float, Dynamic>& _pos,
    const Vector<float, Dynamic>& _vel)
{
    // TODO 
    //
    // Mise à jour des propriétés de chacune des particules à partir des valeurs
    // contenues dans le vecteur d'état.

    for (int i = 0; i < m_particles.size(); i++)
    {
        Particle& p = m_particles[i];

        if (!p.fixed) 
        {
            int pos = i * 2;

            p.x(0) = _pos(pos);
            p.x(1) = _pos(pos + 1);

            p.v(0) = _vel(pos);
            p.v(1) = _vel(pos + 1);
        }
    }
}



/**
 * Construction de la matrice de masses.
 */
void ParticleSystem::buildMassMatrix(Matrix<float, Dynamic, Dynamic>& M)
{
    const int numParticles = m_particles.size();
    const int dim = 2 * numParticles;
    M.resize(dim, dim);
    M.setZero();

    // TODO 
    //
    // Inscrire la masse de chacune des particules dans la matrice de masses M
    //
    for (int i = 0; i < numParticles; ++i)
    {
        Particle& p = m_particles[i];
        float masse = p.fixed ? 1e9 : p.m;
        int pos = i * 2;
        M(pos, pos) = masse;
        M(pos + 1, pos + 1) = masse;
    }
}


/**
 * Construction de la matrice de rigidité.
 */
void ParticleSystem::buildDfDx(Matrix<float, Dynamic, Dynamic>& dfdx)
{
    const int numParticles = m_particles.size();
    const int numSprings = m_springs.size();
    const int dim = 2 * numParticles;
    dfdx.resize(dim, dim);
    dfdx.setZero();

    // Pour chaque ressort...
    for (const Spring& spring : m_springs)
    {
        // TODO
        //
        // Calculer le coefficients alpha et le produit dyadique tel que décrit au cours.
        // Utiliser les indices spring.index0 et spring.index1 pour calculer les coordonnées des endroits affectés.
        //
        // Astuce: créer une matrice de taille fixe 2 par 2 puis utiliser la classe SubMatrix pour accumuler 
        // les modifications sur la diagonale (2 endroits) et pour mettre à jour les blocs non diagonale (2 endroits).

        auto df_moins = this->MRigid(false, m_particles[spring.index0], m_particles[spring.index1], spring)(0, 0);
        auto df_plus = this->MRigid(true, m_particles[spring.index0], m_particles[spring.index1], spring)(0, 0);

        int pos1 = spring.index0 * 2;
        int pos2 = spring.index1 * 2;

        auto submatrice1 = dfdx.block(pos1, pos1, 2, 2);
        auto submatrice2 = dfdx.block(pos2, pos1, 2, 2);
        auto submatrice3 = dfdx.block(pos1, pos2, 2, 2);
        auto submatrice4 = dfdx.block(pos2, pos2, 2, 2);

        submatrice1.sumDiagonalValue(df_moins);
        submatrice2.sumDiagonalValue(df_plus);
        submatrice3.sumDiagonalValue(df_plus);
        submatrice4.sumDiagonalValue(df_moins);
    }
}

Matrix<float, 2, 2> ParticleSystem::MRigid(const bool positive, const Particle& p1, const Particle& p2, const Spring& spring)
{
    Vector<float, 2> diffPosParticles = p2.x - p1.x;
    float alpha = this->CalculAlpha(diffPosParticles, spring);

    Matrix<float, 2, 2> MAlpha(2, 2);

    auto otherCalc = spring.k * (spring.l0 / diffPosParticles.norm());

    // produit dyadique
    auto v = (1 / diffPosParticles.norm()) * diffPosParticles;
    Matrix<float, 2, 2> prodDyadique = v * v;

    if (positive) 
    {
        MAlpha.diagonalValue(alpha);
        return MAlpha + otherCalc * prodDyadique;
    }
    else 
    {
        MAlpha.diagonalValue(-alpha);
        return MAlpha + (-otherCalc * prodDyadique);
    }
}

float ParticleSystem::CalculAlpha(const Vector<float, 2>& diffPosParticles, const Spring& spring) 
{
    return spring.k * (1 - (spring.l0 / diffPosParticles.norm()));
}
