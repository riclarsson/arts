#ifndef montecarlo_h
#define montecarlo_h

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include "messages.h"
#include "arts.h"
#include "ppath.h"
#include "matpackI.h"
#include "special_interp.h"
#include "check_input.h"
#include <stdexcept>
#include <cmath>
#include "rte.h"
#include "lin_alg.h"
#include "auto_md.h"
#include "logic.h"
#include "physics_funcs.h"
#include "xml_io.h"
#include "rng.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;

void Cloudbox_ppathCalc(
        //  Output:
              Ppath&          ppath,
              Ppath&          ppath_step,
        //  Input:
        const Agenda&         ppath_step_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Matrix&         r_geoid,
        const Matrix&         z_surface,
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         rte_pos,
              const Vector&         rte_los);

void Cloudbox_ppath_rteCalc(
                             Ppath&                ppathcloud,
                             Ppath&                ppath,
                             Ppath&                ppath_step,
                             Vector&               rte_pos,
                             Vector&               rte_los,
                             Vector&               cum_l_step,
                             ArrayOfMatrix&        TArray,
                             ArrayOfMatrix&        ext_matArray,
                             ArrayOfVector&        abs_vecArray,
                             Vector&               t_ppath,
                             //Vector&               scat_za_grid,
                             //Vector&               scat_aa_grid,
                             Tensor3&              ext_mat,
                             Matrix&               abs_vec,
                             Numeric&              rte_pressure,
                             Numeric&              rte_temperature,
                             Vector&               rte_vmr_list,
                             Matrix&               iy,
                             GridPos&              rte_gp_p,
                             GridPos&              rte_gp_lat,
                             GridPos&              rte_gp_lon,
                             //Matrix&               i_space,
                             //Matrix&               ground_emission,
                             //Matrix&               ground_los, 
                             //Tensor4&              ground_refl_coeffs,
                             Index&                f_index,
                             Matrix&               pnd_ppath,
                             const Agenda&         ppath_step_agenda,
                             const Index&          atmosphere_dim,
                             const Vector&         p_grid,
                             const Vector&         lat_grid,
                             const Vector&         lon_grid,
                             const Tensor3&        z_field,
                             const Matrix&         r_geoid,
                             const Matrix&         z_surface,
                             const ArrayOfIndex&   cloudbox_limits,
                             const Index&          record_ppathcloud,
                             const Index&          record_ppath,
                             const Agenda&         opt_prop_gas_agenda,
                             const Agenda&         scalar_gas_absorption_agenda,
                             const Index&          stokes_dim,
                             const Tensor3&        t_field,
                             const Tensor4&        vmr_field,
                             const Agenda&         rte_agenda,
                             const Agenda&         iy_space_agenda,
                             const Agenda&         iy_surface_agenda,
                             const Agenda&         iy_cloudbox_agenda,
                             const Vector&         f_grid,
                             const Index&          photon_number,
                             const Index&          scattering_order,
                             const Tensor4&        pnd_field,
                             const ArrayOfSingleScatteringData& scat_data_mono
);


void cloudbox_ppath_start_stepping(
                                   Ppath&          ppath,
                                   const Index&          atmosphere_dim,
                                   ConstVectorView       p_grid,
                                   ConstVectorView       lat_grid,
                                   ConstVectorView       lon_grid,
                                   ConstTensor3View      z_field,
                                   ConstMatrixView       r_geoid,
                                   ConstMatrixView       z_surface,
                                   ConstVectorView       rte_pos,
                                   ConstVectorView       rte_los );
          

void cum_l_stepCalc(
                      Vector& cum_l_step,
                      const Ppath& ppath
                      );

void findZ11max(Vector& Z11maxvector,
        const  ArrayOfSingleScatteringData& scat_data_mono);

Matrix interp( ConstVectorView itw,
               ArrayOfMatrix a,    
               const GridPos&  tc );

Vector interp( ConstVectorView itw,
                ArrayOfVector a,    
               const GridPos&  tc );

void interp_scat_angle_temperature(//Output:
                                   VectorView pha_mat_int,
                                   Numeric& theta_rad,
                                   //Input:
                                   const SingleScatteringData& scat_data,
                                   const Numeric& za_sca,
                                   const Numeric& aa_sca,
                                   const Numeric& za_inc,
                                   const Numeric& aa_inc,
                                   const Numeric& rte_temperature
                                   );

void interpTArray(Matrix& T,
                  Vector& Kabs,
                  Numeric& temperature,
                  MatrixView&  K,
                  Vector& rte_pos,//maybe these should be VectorViews?
                  Vector& rte_los,
                  VectorView& pnd_vec,
                  ArrayOfGridPos& gp,
                  const ArrayOfMatrix& TArray,
                  const ArrayOfMatrix& ext_matArray,
                  const ArrayOfVector& abs_vecArray,
                  const Vector& t_ppath,
                  const Matrix& pnd_ppath,
                  const Vector& cum_l_step,
                  const Numeric& pathlength,
                  const Index& stokes_dim,
                  const Ppath& ppath
                 );

bool is_anyptype30(const ArrayOfSingleScatteringData& scat_data_mono);

void montecarloGetIncoming(
                           Matrix&               iy,
                           Vector&               rte_pos,
                           Vector&               rte_los,
                           GridPos&              rte_gp_p,
                           GridPos&              rte_gp_lat,
                           GridPos&              rte_gp_lon,
                           Ppath&                ppath,
                           Ppath&                ppath_step,
                           //Matrix&               i_space,
                           //Matrix&               ground_emission,
                           //Matrix&               ground_los, 
                           //Tensor4&              ground_refl_coeffs,
                           //Vector&               scat_za_grid,
                           //Vector&               scat_aa_grid,
                           const Agenda&         ppath_step_agenda,
                           const Agenda&         rte_agenda,
                           const Agenda&         iy_space_agenda,
                           const Agenda&         iy_surface_agenda,
                           const Agenda&         iy_cloudbox_agenda,
                           //const Tensor3&        t_field,
                           const Vector&         p_grid,
                           const Vector&         lat_grid,
                           const Vector&         lon_grid,
                           const Tensor3&        z_field,
                           const Matrix&         r_geoid,
                           const Matrix&         z_surface,
                           const ArrayOfIndex&   cloudbox_limits,
                           const Ppath&          ppathcloud,
                           const Index&          atmosphere_dim,
                           const Vector&         f_grid,
                           const Index&          stokes_dim
                           );

void opt_propCalc(
                  MatrixView& K,
                  VectorView& K_abs,
                  const Numeric za,
                  const Numeric aa,
                  const ArrayOfSingleScatteringData& scat_data_mono,
                  const Index&          stokes_dim,
                  const VectorView& pnd_vec,
                  const Numeric& rte_temperature
                  );

void opt_propExtract(
                     MatrixView& K_spt,
                     VectorView& K_abs_spt,
                     const SingleScatteringData& scat_data,
                     const Numeric& za,
                     const Numeric& aa,
                     const Numeric& rte_temperature,
                     const Index& stokes_dim
                     );

void pha_mat_singleCalc(
                        MatrixView& Z,                  
                        Numeric za_sca, 
                        Numeric aa_sca, 
                        Numeric za_inc, 
                        Numeric aa_inc,
                        const ArrayOfSingleScatteringData& scat_data_mono,
                        const Index&          stokes_dim,
                        const VectorView& pnd_vec,
                        const Numeric& rte_temperature
                        );

void pha_mat_singleExtract(
                           MatrixView& Z_spt,
                           const SingleScatteringData& scat_data,
                           const Numeric& za_sca,
                           const Numeric& aa_sca,
                           const Numeric& za_inc,
                           const Numeric& aa_inc,
                           const Numeric& rte_temperature,
                           const Index& stokes_dim
                           );

void ppathRecordMC(
                   const Ppath& ppath,
                   const String name,
                   const Index& photon_number,
                   const Index& scattering_order
                   );


void Sample_los (
                   VectorView& new_rte_los,
                   Numeric& g_los_csc_theta,
                   MatrixView& Z,
                   Rng& rng,
                   const VectorView& rte_los,
                   const ArrayOfSingleScatteringData& scat_data_mono,
                   const Index&          stokes_dim,
                   const VectorView& pnd_vec,
                   const bool& anyptype30,
                   const VectorView& Z11maxvector,
                   Numeric Csca,
                   const Numeric& rte_temperature
                   );

void Sample_ppathlength (
                         Numeric& pathlength, 
                         Numeric& g,
                         Rng& rng,
                         const ArrayOfMatrix& TArray,
                         const ConstVectorView& cum_l_step
                         );

void Sample_ppathlengthLOS (
                         Numeric& pathlength, 
                         Numeric& g,
                         Rng& rng,
                         const ArrayOfMatrix& TArray,
                         const ConstVectorView& cum_l_step
                         );


void TArrayCalc(
                //output
                ArrayOfMatrix& TArray,
                ArrayOfMatrix& ext_matArray,
                ArrayOfVector& abs_vecArray,
                Vector& t_ppath,
                Tensor3& ext_mat,
                Matrix& abs_vec,
                Numeric&   rte_pressure,
                Numeric&   rte_temperature,
                Vector&    rte_vmr_list,
                Matrix&  pnd_ppath,
                //input
                const Ppath& ppath,
                const Agenda& opt_prop_gas_agenda,
                const Agenda& scalar_gas_absorption_agenda,
                const Index& stokes_dim,
                const Vector&    p_grid,
                const Vector&    lat_grid,
                const Vector&    lon_grid,
                const Tensor3&   t_field,
                const Tensor4&   vmr_field,
                const Index&     atmosphere_dim,
                const Tensor4&   pnd_field,
                const ArrayOfSingleScatteringData& scat_data_mono
                );

#endif  // montecarlo_h
