/**
 * @file GraphColoring.cpp
 *
 * @brief Algorithme glouton pour la coloration de graphe.
 *
 * Nom: Erick Santiago Chiappe Reyes
 * Code permanent : CHIE19079806
 * Email : erick-santiago.chiappe-reyes.1@ens.etsmtl.ca
 *
 */


#include "GraphColoring.h"
#include "ParticleSystem.h"

using namespace gti320;

void GraphColoring::color(ParticleSystem& particleSystem)
{
    // La palette de couleurs
    ColorList C;

    std::vector<Particle>& particles = particleSystem.getParticles();
    std::vector<Spring>& springs = particleSystem.getSprings();

    // Initialiser toutes les particules avec color = -1 
    for (Particle& p : particles)
    {
        p.color = -1;
    }

    // TODO Calculer les couleurs de toutes les particules du système. 
    //      Boucler sur chaque particule et appeler la fonction findColor.

    for (Particle& p : particles) 
    {
        p.color = this->findColor(p, particles, springs, C);
    }

    // TODO Construire les partitions qui correspond à chaque couleur. 
    //     Les partitions sont représentées par un tableau d'indices de particules, un pour chaque couleur. 
    //     Stocker les partitions dans m_partitions.
    m_partitions.resize(C.size());
    int i = 0;
    for (Particle& p : particles)
    {
        auto& partition = m_partitions[p.color];
        partition.push_back(i);
        i++;
    }
}

int GraphColoring::findColor(const Particle& p, const std::vector<Particle>& particles, const std::vector<Spring>& springs, ColorList& C) const
{
    // TODO Trouver la première couleur de la palette C qui n'est pas attribuée à une particule voisine. 
    //      Si une couleur est introuvable, ajouter une nouvelle couleur à la palette et retournez la couleur. 
    //      Utiliser la fonction findNeighbors pour assembler une liste de particules qui sont directement connectées à la particule p par un ressort (les voisines).
    int nbColors = C.size() > 0 ? C.size() : 1;
    Vector<int, Dynamic> count(nbColors);
    for (auto y : findNeighbors(p, particles, springs)) 
    {
        if(y->color >= 0) count(y->color)++;
    }

    for (int c = 0; c < nbColors; c++) 
    {
        if (count(c) == 0) 
        {
            if (C.size() == 0) C.push_back(c);
            return c;
        }
    }

    int c = nbColors;
    C.push_back(c);
    return c;
}

NeighborList GraphColoring::findNeighbors(const Particle& p, const std::vector<Particle>& particles, const std::vector<Spring>& springs) const
{
    NeighborList N;

    for (const Spring& s : springs)
    {
        if (&particles[s.index0] == &p)
        {
            N.push_back(&particles[s.index1]);
        }
        else if (&particles[s.index1] == &p)
        {
            N.push_back(&particles[s.index0]);
        }
    }

    return N;
}