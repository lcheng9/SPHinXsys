
#include "unstructured_mesh2.h"

#include "base_particle_dynamics.h"

namespace SPH
{
//=================================================================================================//
BaseInnerRelationInFVM2::BaseInnerRelationInFVM2(RealBody &real_body, ANSYSMesh2 &ansys_mesh)
    : BaseInnerRelation(real_body), real_body_(&real_body),
      node_coordinates_(ansys_mesh.node_coordinates_),
      mesh_topology_(ansys_mesh.mesh_topology_),
      pos_(base_particles_.getVariableDataByName<Vecd>("Position")),
      Vol_(base_particles_.getVariableDataByName<Real>("VolumetricMeasure"))
{
    subscribeToBody();
    inner_configuration_.resize(base_particles_.ParticlesBound(), Neighborhood());
};
//=============================================================================================//
} // namespace SPH
