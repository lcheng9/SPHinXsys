#include "particle_generator_mesh2.h"

#include "base_particle_dynamics.h"

namespace SPH
{
//=================================================================================================//
GeneratingMethod<UnstructuredMesh2>::GeneratingMethod(ANSYSMesh2 &ansys_mesh)
    : elements_centroids_(ansys_mesh.elements_centroids_),
      elements_volumes_(ansys_mesh.elements_volumes_) {}
//=================================================================================================//
ParticleGenerator<BaseParticles, UnstructuredMesh2>::
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles, ANSYSMesh2 &ansys_mesh)
    : ParticleGenerator<BaseParticles>(sph_body, base_particles),
      GeneratingMethod<UnstructuredMesh2>(ansys_mesh) {}
//=================================================================================================//
void ParticleGenerator<BaseParticles, UnstructuredMesh2>::prepareGeometricData()
{
    for (size_t i = 0; i != elements_centroids_.size(); ++i)
    {
        addPositionAndVolumetricMeasure(elements_centroids_[i], elements_volumes_[i]);
    }
}
//=================================================================================================//
} // namespace SPH
