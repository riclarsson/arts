#include <workspace.h>

#include "enumsAtmKey.h"
#include "enumsFieldComponent.h"

////////////////////////////////////////////////////////////////////////////////
// Retrieval code.  This wraps Jacobian and Covmat code.
////////////////////////////////////////////////////////////////////////////////

void RetrievalInit(JacobianTargets& jacobian_targets,
                   CovarianceMatrix& model_state_covariance_matrix,
                   JacobianTargetsDiagonalCovarianceMatrixMap&
                       covariance_matrix_diagonal_blocks) {
  model_state_covariance_matrixInit(model_state_covariance_matrix);
  jacobian_targetsInit(jacobian_targets);
  covariance_matrix_diagonal_blocks.clear();
}

void RetrievalAddSurface(JacobianTargets& jacobian_targets,
                         JacobianTargetsDiagonalCovarianceMatrixMap&
                             covariance_matrix_diagonal_blocks,
                         const SurfaceKey& key,
                         const Numeric& d,
                         const BlockMatrix& matrix,
                         const BlockMatrix& inverse) {
  jacobian_targetsAddSurface(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jacobian_targets.surf().back().type}] = {matrix, inverse};
}

void RetrievalAddSurface(JacobianTargets& jacobian_targets,
                         JacobianTargetsDiagonalCovarianceMatrixMap&
                             covariance_matrix_diagonal_blocks,
                         const SurfaceTypeTag& key,
                         const Numeric& d,
                         const BlockMatrix& matrix,
                         const BlockMatrix& inverse) {
  jacobian_targetsAddSurface(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jacobian_targets.surf().back().type}] = {matrix, inverse};
}

void RetrievalAddSurface(JacobianTargets& jacobian_targets,
                         JacobianTargetsDiagonalCovarianceMatrixMap&
                             covariance_matrix_diagonal_blocks,
                         const SurfacePropertyTag& key,
                         const Numeric& d,
                         const BlockMatrix& matrix,
                         const BlockMatrix& inverse) {
  jacobian_targetsAddSurface(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jacobian_targets.surf().back().type}] = {matrix, inverse};
}

void RetrievalAddAtmosphere(JacobianTargets& jacobian_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const AtmKey& key,
                            const Numeric& d,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse) {
  jacobian_targetsAddAtmosphere(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jacobian_targets.atm().back().type}] = {matrix, inverse};
}

void RetrievalAddAtmosphere(JacobianTargets& jacobian_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const SpeciesEnum& key,
                            const Numeric& d,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse) {
  jacobian_targetsAddAtmosphere(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jacobian_targets.atm().back().type}] = {matrix, inverse};
}

void RetrievalAddAtmosphere(JacobianTargets& jacobian_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const SpeciesIsotope& key,
                            const Numeric& d,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse) {
  jacobian_targetsAddAtmosphere(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jacobian_targets.atm().back().type}] = {matrix, inverse};
}

void RetrievalAddAtmosphere(JacobianTargets& jacobian_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const QuantumIdentifier& key,
                            const Numeric& d,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse) {
  jacobian_targetsAddAtmosphere(jacobian_targets, key, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jacobian_targets.atm().back().type}] = {matrix, inverse};
}

void RetrievalAddSpeciesVMR(JacobianTargets& jacobian_targets,
                            JacobianTargetsDiagonalCovarianceMatrixMap&
                                covariance_matrix_diagonal_blocks,
                            const SpeciesEnum& species,
                            const Numeric& d,
                            const BlockMatrix& matrix,
                            const BlockMatrix& inverse) {
  jacobian_targetsAddSpeciesVMR(jacobian_targets, species, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jacobian_targets.atm().back().type}] = {matrix, inverse};
}

void RetrievalAddSpeciesIsotopologueRatio(
    JacobianTargets& jacobian_targets,
    JacobianTargetsDiagonalCovarianceMatrixMap&
        covariance_matrix_diagonal_blocks,
    const SpeciesIsotope& species,
    const Numeric& d,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  jacobian_targetsAddSpeciesIsotopologueRatio(jacobian_targets, species, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jacobian_targets.atm().back().type}] = {matrix, inverse};
}

void RetrievalAddMagneticField(JacobianTargets& jacobian_targets,
                               JacobianTargetsDiagonalCovarianceMatrixMap&
                                   covariance_matrix_diagonal_blocks,
                               const String& component,
                               const Numeric& d,
                               const BlockMatrix& matrix,
                               const BlockMatrix& inverse) {
  jacobian_targetsAddMagneticField(jacobian_targets, component, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jacobian_targets.atm().back().type}] = {matrix, inverse};
}

void RetrievalAddWindField(JacobianTargets& jacobian_targets,
                           JacobianTargetsDiagonalCovarianceMatrixMap&
                               covariance_matrix_diagonal_blocks,
                           const String& component,
                           const Numeric& d,
                           const BlockMatrix& matrix,
                           const BlockMatrix& inverse) {
  jacobian_targetsAddWindField(jacobian_targets, component, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jacobian_targets.atm().back().type}] = {matrix, inverse};
}

void RetrievalAddTemperature(JacobianTargets& jacobian_targets,
                             JacobianTargetsDiagonalCovarianceMatrixMap&
                                 covariance_matrix_diagonal_blocks,
                             const Numeric& d,
                             const BlockMatrix& matrix,
                             const BlockMatrix& inverse) {
  jacobian_targetsAddTemperature(jacobian_targets, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jacobian_targets.atm().back().type}] = {matrix, inverse};
}

void RetrievalAddPressure(JacobianTargets& jacobian_targets,
                          JacobianTargetsDiagonalCovarianceMatrixMap&
                              covariance_matrix_diagonal_blocks,
                          const Numeric& d,
                          const BlockMatrix& matrix,
                          const BlockMatrix& inverse) {
  jacobian_targetsAddPressure(jacobian_targets, d);
  covariance_matrix_diagonal_blocks[JacobianTargetType{
      jacobian_targets.atm().back().type}] = {matrix, inverse};
}

void RetrievalAddSensorFrequencyPolyFit(
    JacobianTargets& jacobian_targets,
    JacobianTargetsDiagonalCovarianceMatrixMap&
        covariance_matrix_diagonal_blocks,
    const ArrayOfSensorObsel& measurement_sensor,
    const Numeric& d,
    const Index& elem,
    const Index& polyorder,
    const BlockMatrix& matrix,
    const BlockMatrix& inverse) {
  jacobian_targetsAddSensorFrequencyPolyFit(
      jacobian_targets, measurement_sensor, d, elem, polyorder);
      auto keyk = JacobianTargetType{
      jacobian_targets.sensor().back().type};
  covariance_matrix_diagonal_blocks[keyk] = {matrix, inverse};
}
