#ifndef _Castro_F_H_
#define _Castro_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void ca_network_init();

  void burner_init();

  void ca_extern_init(const int* name, const int* namlen);

  void ca_get_num_spec(int* nspec);

  void ca_get_num_aux(int* naux);

  void ca_get_spec_names(int* spec_names, int* ispec, int* len);

  void ca_get_spec_az(int* ispec, amrex::Real* A, amrex::Real* Z);

  void ca_get_aux_names(int* aux_names, int* iaux, int* len);

  void ca_get_qvar(int* qvar);

  void ca_get_nqaux(int* nqaux);

  void ca_set_special_tagging_flag(const amrex::Real& max_den, const int* flag);

  void ca_get_method_params(int* HYP_GROW);

  void ca_set_method_params
    (const int& dm,
     const int& Density, const int& Xmom, 
     const int& Eden,    const int& Eint, 
     const int& Temp     , const int& FirstAdv, 
     const int& FirstSpec, const int& FirstAux, 
     const int& NumAdv);

  void ca_set_castro_method_params();

  void ca_init_godunov_indices();

  void ca_set_problem_params
    (const int& dm,
     const int* physbc_lo, const int* physbc_hi,
     const int& Interior_value, const int& Inflow_value,
     const int& Outflow_value, const int& Symmetry_value,
     const int& SlipWall_value, const int& NoSlipWall_value,
     const int& coord_type,
     const amrex::Real* problo, const amrex::Real* probhi,
     const amrex::Real* center);

  void ca_set_grid_info
    (const int& max_level, const amrex::Real* dx_level,
     const int* domlo_level, const int* domhi_level,
     const int* ref_ratio, const int* n_error_buf,
     const int* blocking_factor);

  void ca_get_tagging_params(const int* name, const int* namlen);

  void ca_ctoprim(const int* lo, const int* hi,
	          const amrex::Real* u, const int* u_lo, const int* u_hi,
		  const amrex::Real* q, const int* q_lo, const int* q_hi,
		  const amrex::Real* qaux, const int* qa_lo, const int* qa_hi, const int* idx);

  void ca_check_initial_species
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state), const int* idx);

  void ca_enforce_minimum_density
    (const amrex::Real* S_old, const int* s_old_lo, const int* s_old_hi,
     const amrex::Real* S_new, const int* s_new_lo, const int* s_new_hi,
     const amrex::Real* vol, const int* vol_lo, const int* vol_hi,
     const int* lo, const int* hi,
     const amrex::Real* frac_change, const int* verbose, const int* idx);

  void ca_normalize_species
    (BL_FORT_FAB_ARG_3D(S_new), const int* lo, const int* hi, const int* idx);

  void ca_compute_temp
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(state), const int* idx);

  void ca_enforce_consistent_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state), const int* idx);

  void ca_reset_internal_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(S_new),
     const int& verbose, const int* idx);

  void ca_hypfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_denfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_xmomfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_ymomfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_zmomfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_estdt
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, amrex::Real* dt, const int* idx);

  void ca_check_timestep
    (const BL_FORT_FAB_ARG_3D(state_old),
     const BL_FORT_FAB_ARG_3D(state_new),
     const int* lo, const int* hi,
     const amrex::Real* dx, const amrex::Real* dt_old, amrex::Real* dt);

  void ca_mol_single_stage
    (const amrex::Real* time,
     const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(stateout),
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(qaux),
     BL_FORT_FAB_ARG_3D(update),
     const amrex::Real dx[], const amrex::Real* dt,
     D_DECL(BL_FORT_FAB_ARG_3D(xflux),
            BL_FORT_FAB_ARG_3D(yflux),
            BL_FORT_FAB_ARG_3D(zflux)),
     D_DECL(BL_FORT_FAB_ARG_3D(xarea),
            BL_FORT_FAB_ARG_3D(yarea),
            BL_FORT_FAB_ARG_3D(zarea)),
     const BL_FORT_FAB_ARG_3D(volume),
     amrex::Real* cflLoc,
     const int&  verbose, const int* idx);

  void flush_output();

  void ca_summass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass, const int* idx);

  void ca_derpres
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

  void ca_dervel
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

#ifdef __cplusplus
}
#endif


BL_FORT_PROC_DECL(CA_INITDATA,ca_initdata)
    (const int& level, const amrex::Real& time,
     const int* lo, const int* hi,
     const int& num_state,
     BL_FORT_FAB_ARG(state),
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* xhi);

#endif
