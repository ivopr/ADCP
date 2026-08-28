/*
** Module with all the simulation parameters **
** and its default values
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>
#include<float.h>
#include"error.h"
#include"params.h"



/***********************************************************/
/****       INITIALISING AND FINALISING  ROUTINES       ****/
/***********************************************************/

/* initialising simulation parameters (including model parameters) */
/* default values should be set here */
void param_initialise(simulation_params *self) {
  /* general simulation */
  self->infile = stdin;
  self->outfile = stdout;
  self->infile_name = NULL;
  self->outfile_name = NULL;
  self->pace = 0;
  self->stretch = 16;
  self->tmask = 0x0;
  self->seed = 0;
  self->prm = NULL;
  self->acceptance_rate = 0.5;
  self->amplitude = -0.1;
  self->keep_amplitude_fixed = 0;
  self->acceptance = 1.0;
  self->accept_counter = 0;
  self->reject_counter = 0;
  self->acceptance_rate_tolerance = 0.03;
  self->amplitude_changing_factor = 0.9;
  self->number_initial_MC = 0;
  self->MC_lookup_table = NULL; //lookup table for random moves
  self->MC_lookup_table_n = NULL; //number of valid elements in the lookup table for random moves
  /* peptide */
  self->seq = NULL;
  self->sequence = NULL;
  self->NAA = 0;
  self->Nchains = 0;
  /* thermodynamic beta */
  self->thermobeta = 0;
  self->lowtemp = 0;
  self->beta1 = 1;
  self->beta2 = 0;
  self->bstp = 1;
  self->intrvl = 16384;
  self->nswap_per_try = 16384;
  self->energy_gradient.assign(36, 0);
  self->energy_probe_1_this.assign(36, 0);
  self->energy_probe_1_last.assign(36, 0);
  self->energy_probe_1_calc.assign(36, 1);
  /* nested sampling */
  self->NS = 0;
  self->iter = 0;
  self->iter_start = 0;
  self->iter_max = 1000;
  self->logLstar = DBL_MAX;
  self->logZ = -DBL_MAX;
  self->logfactor = 1;
  self->alpha = 1;
  self->logX = 0;
  self->logX_start = 0;
  self->Delta_logX = 0;
  self->log_DeltaX = 0;
  self->H = 0;
  self->N = 0;

  /* checkpointing */
  self->num_NS_per_checkpoint = 0;
  self->checkpoint_filename = NULL;
  self->checkpoint_file = NULL;
  self->checkpoint_counter = 0;
  self->restart_from_checkpoint = 0;
  self->checkpoint = 0;

  model_param_initialise(&(self->protein_model));
  flex_params_initialise(&(self->flex_params));

}

/* finalising simulation parameters (including model parameters) */
/* default values should be set here */
void param_finalise(simulation_params *self) {
  /* general simulation */

  if (self->infile && self->infile!=stdin) {
    fclose(self->infile);
    self->infile = NULL;
  }

  if (self->outfile && self->outfile!=stdout) {
    fclose(self->outfile);
    self->outfile = NULL;
  }

  if (self->infile_name) free(self->infile_name);
  if (self->outfile_name) free(self->outfile_name);

  self->pace = 0;
  self->stretch = 0;
  self->tmask = 0x0;
  self->seed = 0;

  if (self->prm) free(self->prm);

  self->acceptance_rate = 0.;
  self->amplitude = 0.;
  self->keep_amplitude_fixed = 0;
  self->acceptance = 0.;
  self->accept_counter = 0;
  self->reject_counter = 0;
  self->acceptance_rate_tolerance = 0;
  self->amplitude_changing_factor = 0;
  self->number_initial_MC = 0;
  if (self->MC_lookup_table) free(self->MC_lookup_table);
  if (self->MC_lookup_table_n) free(self->MC_lookup_table_n);


  /* peptide */

  if (self->seq) free(self->seq);
  if (self->sequence) free(self->sequence);
  self->NAA = 0;
  self->Nchains = 0;
  /* thermodynamic beta */
  self->thermobeta = 0;
  self->lowtemp = 0;
  self->beta1 = 0;
  self->beta2 = 0;
  self->bstp = 0;
  self->intrvl = 0;
  self->nswap_per_try = 0;
  self->energy_gradient.clear();
  self->energy_probe_1_calc.clear();
  self->energy_probe_1_this.clear();
  self->energy_probe_1_last.clear();
  /* nested sampling */
  self->NS = 0;
  self->iter = 0;
  self->iter_start = 0;
  self->iter_max = 0;
  self->logLstar = 0;
  self->logZ = 0;
  self->logfactor = 0;
  self->alpha = 0;
  self->logX = 0;
  self->logX_start = 0;
  self->Delta_logX = 0;
  self->log_DeltaX = 0;
  self->H = 0;
  self->N = 0;

  /* checkpointing */
  self->num_NS_per_checkpoint = 0;
  if (self->checkpoint_filename) free(self->checkpoint_filename);
  if (self->checkpoint_file) {
	fclose(self->checkpoint_file);
	self->checkpoint_file = NULL;
  }
  self->checkpoint_counter = 0;
  self->restart_from_checkpoint = 0;
  self->checkpoint = 0;


  model_param_finalise(&(self->protein_model));
  flex_params_finalise(&(self->flex_params));

}


/* initialise model parameters */
void model_param_initialise(model_params *self) {

  /* gamma atoms */
  self->fixit = 1;
  self->use_gamma_atoms = LINUS_GAMMA;
  self->use_original_gamma_atoms = 1;
  self->use_3_states = 1;
  self->fix_chi_angles = 0;
  self->fix_CA_atoms = 0;
  /* atomic radii */
  /* radii from Ward et al, 1999 */
  // rca = 1.75, rcb = 1.75, rc = 1.65, rn = 1.55, ro = 1.40; 
  /* radii from Hopfinger, 1973 */
  // rca = 1.57, rcb = 1.57, rc = 1.42, rn = 1.29, ro = 1.29
  // rs = 1.8 from wikipedia I'm sorry ,NB
  /* LINUS */
  // rca = 1.85, rcb = 2.0, rc = 1.85, rn = 1.75, ro = 1.6, rs = 2.0; 
  self->rca = 2.43;
  self->rcb = 1.97;
  self->rc = 1.82;
  self->rn = 1.74;
  self->ro = 1.98;
  self->rs = 3.10;
  self->rring = 2.00;
  self->vdw_depth_ca = 0.018;
  self->vdw_depth_cb = 0.018;
  self->vdw_depth_c = 0.018;
  self->vdw_depth_n = 0.018;
  self->vdw_depth_o = 0.018;
  self->vdw_depth_s = 0.018;
  self->vdw_depth_ring = 0.2;

  /* The vdW cutoff distances have to be calculated later, because
     the cutoff calculating routine depends on vdw.c */
  self->vdw_gamma_gamma_cutoff = NULL;
  self->vdw_gamma_nongamma_cutoff = NULL;
  self->vdw_backbone_cutoff = 500;
  self->vdw_use_extended_cutoff = 0;
  self->vdw_extended_cutoff = DEFAULT_EXTENDED_VDW_CUTOFF_GAMMA;
  self->vdw_potential = LJ_VDW_POTENTIAL;
  self->vdw_function = NULL;
  self->vdw_clash_energy_at_hard_cutoff = 30; //default value for LJ
  self->vdw_lj_neighbour_hard = 0;
  self->vdw_lj_hbonded_hard = 0;
  
  self->rel_vdw_cutoff = 2.0;
  self->vdw_uniform_depth = 0;
  vdw_param_zero(self);
  vdw_param_calculate(self);

  /* stress */
  self->stress_k = 98.;
  self->stress_angle = 1.20427718387608740808;

  /* hydrogen bond */
  self->hboh2 = 4.04;
  self->hbohn = 0.928;
  self->hbcoh = 0.772;
  //self->hboh_decay_width = 1.0;
  //self->hbohn_decay_width = 0.1;
  //self->hbcoh_decay_width = 0.1;
  self->hbs = 4.98;
  /* contact parameters */
  self->touch2 = 38.44;
  self->part = 1.0;
  self->split = -1.0;
  self->sts = 0.0;
  /* biasing force constants */
  self->contact_map_file.clear();
  self->bias_eta_beta = 3.7;
  self->bias_eta_alpha = 15.3;
  self->bias_kappa_alpha_3 = 0.0;
  self->bias_kappa_alpha_4 = 0.0;
  self->bias_kappa_beta = 0.85;
  self->prt = 1.;
  self->bias_r_alpha = 5.39;
  self->bias_r_beta = 5.39;
  /* hydrophobicity */
  self->kauzmann_param = 0.122;
  self->hydrophobic_cutoff_range = 2.8;
  self->hydrophobic_min_separation = 2;
  // 1/dist potential form
  //self->hydrophobic_min_cutoff = 2.0;
  //self->hydrophobic_max_cutoff = 8.0;
  //self->hydrophobic_max_Eshift = 0.125;
  // Spline potential form
  //self->hydrophobic_r = 6.0;
  //self->hydrophobic_half_delta = 2.0;
  /* electrostatics */
  self->recip_dielectric_param = 0.;
  self->debye_length_param = 0.;
  self->electrostatic_min_separation = 2;
  /* side chain hydrogen bond parameters */
  self->sidechain_hbond_strength_s2b = 0.;
  self->sidechain_hbond_strength_b2s = 0.;
  self->sidechain_hbond_strength_s2s = 0.;
  self->sidechain_hbond_cutoff = 0.;
  self->sidechain_hbond_decay_width = 0.;
  self->sidechain_hbond_min_separation = 2;
  self->sidechain_hbond_angle_cutoff = -1.;
  /* secondary radius of gyration */
  self->srgy_param = 0.0;
  self->srgy_offset = 0.0;
  self->hphobic_srgy_param = 0.0;
  self->hphobic_srgy_offset = 0.0;
  /*sbond */
  self->Sbond_strength = 0;
  self->Sbond_distance = 2.2;
  self->Sbond_cutoff = 0.2;
  self->Sbond_dihedral_cutoff = 0.35;

  /* fixed amino acids */
  self->fixed_aalist_file.clear();

  /*optimizing strategy*/
  self->opt = 0;
  self->opt_totE_weight = 1.0;
  self->opt_firstlastE_weight = 0.0;
  self->opt_extE_weight = 0.0;

  /* external potential */
  self->external_potential_type = 0;
  for (int i=0; i<3; i++) {
    self->external_direction[i] = EXTERNAL_NONE;
    self->external_k[i] = 0.0;
    self->external_r0[i] = 0.0;
  }
  self->external_ztip = 0.0;
  self->external_constrained_aalist_file.clear();
  self->external_potential_type2 = 0;
  for (int i=0; i<3; i++) {
    self->external_direction2[i] = EXTERNAL_NONE;
    self->external_k2[i] = 0.0;
    self->external_r02[i] = 0.0;
  }
  self->external_ztip2 = 0.0;
  self->external_constrained_aalist_file2.clear();

  //CAUTION!: aadict.c depends on params.c's model_params.  This means that
  //    initialize_sidechain_properties will have to be called after all updates
  //    of the vdW parameters; it can't be called from here, due to circular dependencies.
  self->sidechain_properties = (sidechain_properties_*)calloc( 31, sizeof(sidechain_properties_) );
  /* vdw parameters might have changed */
  //initialize_sidechain_properties(self);

}

/* finalise model parameters */
void model_param_finalise(model_params *self) {

  /* gamma atoms */
  self->use_gamma_atoms = 0;
  self->use_original_gamma_atoms = 0;
  self->use_3_states = 0;
  self->fix_chi_angles = 0;
  self->fix_CA_atoms = 0;
  /* atomic radii */
  self->rca = 0.;
  self->rcb = 0.;
  self->rc = 0.;
  self->rn = 0.;
  self->ro = 0.;
  self->rs = 0.;
  self->rring = 0.;
  self->vdw_depth_ca = 0.;
  self->vdw_depth_cb = 0.;
  self->vdw_depth_c = 0.;
  self->vdw_depth_n = 0.;
  self->vdw_depth_o = 0.;
  self->vdw_depth_s = 0.;
  self->vdw_depth_ring = 0.;
  vdw_param_zero(self);
  if (self->vdw_gamma_gamma_cutoff) free(self->vdw_gamma_gamma_cutoff);
  if (self->vdw_gamma_nongamma_cutoff) free(self->vdw_gamma_nongamma_cutoff);
  self->vdw_backbone_cutoff = 0;
  self->vdw_use_extended_cutoff = 0;
  self->vdw_extended_cutoff = 0;
  self->vdw_potential = 0;
  self->vdw_function = NULL;
  self->vdw_clash_energy_at_hard_cutoff = 0;
  self->vdw_lj_neighbour_hard = 0;
  self->vdw_lj_hbonded_hard = 0;

  /* stress */
  self->stress_k = 0.;
  self->stress_angle = 0.;

  /* hydrogen bond */
  self->hboh2 = 0.;
  self->hbohn = 0.;
  self->hbcoh = 0.;
  //self->hboh_decay_width = 0.;
  //self->hbohn_decay_width = 0.;
  //self->hbcoh_decay_width = 0.;
  self->hbs = 0.;
  /* contact parameters */
  self->touch2 = 0.;
  self->part = 0.;
  self->split = 0.;
  self->sts = 0.;
  /* biasing force constants */
  self->contact_map_file.clear();
  self->bias_eta_beta = 0.;
  self->bias_eta_alpha = 0.;
  self->bias_kappa_alpha_3 = 0.;
  self->bias_kappa_alpha_4 = 0.;
  self->bias_kappa_beta = 0.;
  self->prt = 0.;
  self->bias_r_alpha = 0.;
  self->bias_r_beta = 0.;
  /* hydrophobicity */
  self->kauzmann_param = 0.;
  self->hydrophobic_cutoff_range = 0.;
  self->hydrophobic_min_separation = 0;
  //self->hydrophobic_min_cutoff = 0.;
  //self->hydrophobic_max_cutoff = 0.;
  //self->hydrophobic_max_Eshift = 0.;
  //self->hydrophobic_r = 0.;
  //self->hydrophobic_half_delta = 0.;
  /* electrostatics */
  self->recip_dielectric_param = 0.;
  self->debye_length_param = 0.;
  self->electrostatic_min_separation = 0;
  /* side chain hydrogen bond parameters */
  self->sidechain_hbond_strength_s2b = 0.;
  self->sidechain_hbond_strength_b2s = 0.;
  self->sidechain_hbond_strength_s2s = 0.;
  self->sidechain_hbond_cutoff = 0.;
  self->sidechain_hbond_decay_width = 0.;
  self->sidechain_hbond_min_separation = 0;
  self->sidechain_hbond_angle_cutoff = -0.;
  /* secondary radius of gyration */
  self->srgy_param = 0.;
  self->srgy_offset = 0.;
  self->hphobic_srgy_param = 0.;
  self->hphobic_srgy_offset = 0.;

  /* fixed amino acids */
  self->fixed_aalist_file.clear();

  /*optimizing strategy*/
  self->opt = 0;
  self->opt_totE_weight = 1.0;
  self->opt_firstlastE_weight = 0.0;
  self->opt_extE_weight = 0.0;

  /* external potential */
  self->external_potential_type = 0;
  for (int i=0; i<3; i++) {
    self->external_direction[i] = EXTERNAL_NONE;
    self->external_k[i] = 0.0;
    self->external_r0[i] = 0.0;
  }
  self->external_ztip = 0.0;
  self->external_constrained_aalist_file.clear();
  /* external potential #2 */
  self->external_potential_type2 = 0;
  for (int i=0; i<3; i++) {
    self->external_direction2[i] = EXTERNAL_NONE;
    self->external_k2[i] = 0.0;
    self->external_r02[i] = 0.0;
  }
  self->external_ztip2 = 0.0;
  self->external_constrained_aalist_file2.clear();

  if (self->sidechain_properties) free(self->sidechain_properties);
}


void flex_params_initialise(FLEX_params *self){
  self->number_of_processors = 0;
  self->flex_cmd.clear();
  self->output_path = "/tmp/";
  self->outputpdb_filename = "out.pdb";
  self->fromnode = 7;
  self->tonode = 7;
  self->step = 0.01;
  self->freq = 100;
  self->totalconf = 100;
  self->flex_dir.clear();
  self->filenames_to_read_in.clear();
  self->size_of_filename_to_read_in = 0;
  self->hstrength_cutoff = 1.0;
  self->only_bias_hbonds = 0;

  self->acceptance_rate_aim = 0.5;
  self->acceptance_rate_tolerance = 0.03;
  self->flex_changing_factor = 0.8;

  self->MCiter = 5000;

}

void flex_params_finalise(FLEX_params *self){
  self->number_of_processors = 0;
  self->fromnode = 7;
  self->tonode = 7;
  self->output_path.clear();
  self->outputpdb_filename.clear();
  self->flex_cmd.clear();
  self->flex_dir.clear();
  self->step = 0.01;
  self->freq = 100;
  self->totalconf = 100;

  self->filenames_to_read_in.clear();
  self->size_of_filename_to_read_in = 0;

  self->hstrength_cutoff = 1.0;
  self->only_bias_hbonds = 0;

  self->acceptance_rate_aim = 0.1;
  self->acceptance_rate_tolerance = 0.03;
  self->flex_changing_factor = 0.8;

}



/* vdW depths and shifts at the cutoffs for the atomic pairs */
void vdw_param_zero(model_params *self) {

   /* depths */

   self->vdw_depth_ca_ca = 0;
   self->vdw_depth_cb_cb = 0;
   self->vdw_depth_c_c   = 0;
   self->vdw_depth_n_n   = 0;
   self->vdw_depth_o_o   = 0;

   self->vdw_depth_ca_cb = 0;
   self->vdw_depth_ca_c  = 0;
   self->vdw_depth_ca_n  = 0;
   self->vdw_depth_ca_o  = 0;

   self->vdw_depth_cb_c  = 0;
   self->vdw_depth_cb_n  = 0;
   self->vdw_depth_cb_o  = 0;

   self->vdw_depth_c_n   = 0;
   self->vdw_depth_c_o   = 0;

   self->vdw_depth_n_o   = 0;

   /* shifts */
   self->vdw_shift = 0;

   self->vdw_Eshift_ca_ca = 0;
   self->vdw_Eshift_cb_cb = 0;
   self->vdw_Eshift_c_c   = 0;
   self->vdw_Eshift_n_n   = 0;
   self->vdw_Eshift_o_o   = 0;

   self->vdw_Eshift_ca_cb = 0;
   self->vdw_Eshift_ca_c  = 0;
   self->vdw_Eshift_ca_n  = 0;
   self->vdw_Eshift_ca_o  = 0;

   self->vdw_Eshift_cb_c  = 0;
   self->vdw_Eshift_cb_n  = 0;
   self->vdw_Eshift_cb_o  = 0;

   self->vdw_Eshift_c_n   = 0;
   self->vdw_Eshift_c_o   = 0;

   self->vdw_Eshift_n_o   = 0;

}


/***********************************************************/
/****                DEFAULT  PARAMETERS                ****/
/***********************************************************/


/* Set CD learnt values as default for all parameters learnt
   when using the LJ vdW interaction potential */
/* TODO: update parameters */
void set_lj_default_params(model_params *self) {

  /* The vdW cutoff distances have to be calculated later, because
     the cutoff calculating routine depends on vdw.c */
  self->vdw_potential = LJ_VDW_POTENTIAL;

  /* vdW  */
  /* atomic radii */
  /* radii from Ward et al, 1999 */
  // rca = 1.75, rcb = 1.75, rc = 1.65, rn = 1.55, ro = 1.40; 
  /* radii from Hopfinger, 1973 */
  // rca = 1.57, rcb = 1.57, rc = 1.42, rn = 1.29, ro = 1.29
  // rs = 1.8 from wikipedia I'm sorry ,NB
  /* LINUS */
  // rca = 1.85, rcb = 2.0, rc = 1.85, rn = 1.75, ro = 1.6, rs = 2.0; 
  self->rca = 2.43;
  self->rcb = 1.97;
  self->rc = 1.82;
  self->rn = 1.74;
  self->ro = 1.98;
  self->rs = 3.10;
  self->vdw_depth_ca = 0.018;
  self->vdw_depth_cb = 0.018;
  self->vdw_depth_c = 0.018;
  self->vdw_depth_n = 0.018;
  self->vdw_depth_o = 0.018;
  self->vdw_depth_s = 0.018;

  vdw_param_zero(self);
  vdw_param_calculate(self);

  self->vdw_clash_energy_at_hard_cutoff = 30; //default value for LJ

  /* stress */
  self->stress_k = 98.;

  /* hydrogen bond */
  self->hboh2 = 4.04;
  self->hbohn = 0.928;
  self->hbcoh = 0.772;
  self->hbs = 4.98;

  /* biasing force constants */
  self->bias_eta_beta = 3.7;
  self->bias_eta_alpha = 15.3;
  self->bias_kappa_beta = 0.85;
  self->bias_r_beta = 5.39;

  /* hydrophobicity */
  self->kauzmann_param = 0.122;
  //self->hydrophobic_cutoff_range = 2.8;
  //self->hydrophobic_min_separation = 2;

  /* electrostatics */
  /* side chain hydrogen bond parameters */
  /* secondary radius of gyration */

}


/* Set CD learnt values as default for all parameters learnt
   when using the hard cutoff vdW interaction potential */
/* TODO: update parameters */
void set_hard_cutoff_default_params(model_params *self) {

  /* The vdW cutoff distances have to be calculated later, because
     the cutoff calculating routine depends on vdw.c */
  self->vdw_potential = HARD_CUTOFF_VDW_POTENTIAL;

  /* vdW  */
  /* atomic radii */
  /* radii from Ward et al, 1999 */
  // rca = 1.75, rcb = 1.75, rc = 1.65, rn = 1.55, ro = 1.40; 
  /* radii from Hopfinger, 1973 */
  // rca = 1.57, rcb = 1.57, rc = 1.42, rn = 1.29, ro = 1.29
  // rs = 1.8 from wikipedia I'm sorry ,NB
  /* LINUS */
  // rca = 1.85, rcb = 2.0, rc = 1.85, rn = 1.75, ro = 1.6, rs = 2.0; 
  self->rca = 1.57;
  self->rcb = 1.57;
  self->rc = 1.42;
  self->rn = 1.29;
  self->ro = 1.29;
  self->rs = 2.00;
  self->vdw_depth_ca = 0;
  self->vdw_depth_cb = 0;
  self->vdw_depth_c = 0;
  self->vdw_depth_n = 0;
  self->vdw_depth_o = 0;
  self->vdw_depth_s = 0;

  vdw_param_zero(self);
  vdw_param_calculate(self);

  self->vdw_clash_energy_at_hard_cutoff = 10.4; //default value for LJ

  /* stress */
  self->stress_k = 90.;

  /* hydrogen bond */
  self->hboh2 = 4.05;
  self->hbohn = 0.930;
  self->hbcoh = 0.770;
  self->hbs = 4.95;

  /* biasing force constants */
  self->bias_eta_beta = 4.5;
  self->bias_eta_alpha = 18;
  self->bias_kappa_beta = 0.8;
  self->prt = 1.;
  self->bias_r_beta = 5.65;

  /* hydrophobicity */
  self->kauzmann_param = 0.13;
  //self->hydrophobic_cutoff_range = 2.8;
  //self->hydrophobic_min_separation = 2;

  /* electrostatics */
  /* side chain hydrogen bond parameters */
  /* secondary radius of gyration */

}


/***********************************************************/
/****                 COPYING  ROUTINES                 ****/
/***********************************************************/


/* Copy a string, making sure the null character is copied, too */
/* This will allocate new memory for the new string, so it has to have been freed before if it was used. */
void copy_string(char** const to, const char* const from){
  if (from==NULL) {
    *to=NULL;
  }
  else {
    //if (*to) free(*to); //self would crash things in model_params_copy
    //fprintf(stderr,"string from: %s, of length %d.\n",from,(int)strlen(from));
    *to = (char*)malloc((strlen(from)+1)*sizeof(char));
    if (*to == NULL) stop("Unable to allocate memory in copy_string.");
    strcpy(*to, from);
    //fprintf(stderr,"string from: %s, to: %s.\n",from,*to);
  }
}

/* Copy model parameters */
void model_params_copy(model_params *to, model_params *from) {
  *to = *from;
  /* contact_map_file, fixed_aalist_file, external_constrained_aalist_file(2)
     are std::string, already deep-copied by "*to = *from" above. */
  /* sidechain_properties */
  sidechain_properties_ *temp;
  temp = (sidechain_properties_*)malloc(sizeof(sidechain_properties_) * 31);
  if (!temp) stop("Unable to allocate temp memory in model_params_copy.");
  to->sidechain_properties = temp;
  memcpy(to->sidechain_properties, from->sidechain_properties,
	31 * sizeof(sidechain_properties_));
  /* vdw cutoff matrices */
  if (from->vdw_gamma_gamma_cutoff) {
    double *temp1;
    temp1 = (double*)malloc(sizeof(double) * 702);
    if (!temp1) stop("Unable to allocate temp1 memory in model_params_copy.");
    to->vdw_gamma_gamma_cutoff = temp1;
    memcpy(to->vdw_gamma_gamma_cutoff, from->vdw_gamma_gamma_cutoff,
	702 * sizeof(double));
  }
  if (from->vdw_gamma_nongamma_cutoff) {
    double *temp2;
    temp2 = (double*)malloc(sizeof(double) * 702);
    if (!temp2) stop("Unable to allocate temp1 memory in model_params_copy.");
    to->vdw_gamma_nongamma_cutoff = temp2;
    memcpy(to->vdw_gamma_nongamma_cutoff, from->vdw_gamma_nongamma_cutoff,
	702 * sizeof(double));
  }
}


/* Copying simulation parameters (including model parameters) */
void sim_params_copy(simulation_params *to, simulation_params *from) {

  *to = *from;

  //TODO file pointers are copied, and files will be closed if the copy is finalised
  //to->infile = from->infile;
  //to->outfile = from->outfile;
  //to->checkpoint_file = from->checkpoint_file;

  //strings
  copy_string(&(to->prm), from->prm);
  copy_string(&(to->seq), from->seq);
  copy_string(&(to->sequence), from->sequence);
  copy_string(&(to->infile_name), from->infile_name);
  copy_string(&(to->outfile_name), from->outfile_name);
  copy_string(&(to->checkpoint_filename), from->checkpoint_filename);

  // energy_gradient, energy_probe_1_this, energy_probe_1_last, energy_probe_1_calc
  // are std::vector members, already deep-copied by the "*to = *from" above.

  int* temp2;
  if (from->seq != NULL && from->sequence != NULL && from->MC_lookup_table != NULL && from->MC_lookup_table_n != NULL) {
    int N = 4 * (from->NAA + strlen(from->sequence) - strlen(from->seq));
    temp2 = (int*)malloc(sizeof(int) * N);
    if (!temp2) stop("Unable to allocate temp memory2 in sim_params_copy.");
    to->MC_lookup_table = temp2;
    memcpy(to->MC_lookup_table, from->MC_lookup_table, sizeof(int) * N);
    //MC lookup table length
    temp2 = (int*)malloc(sizeof(int) * 4);
    if (!temp2) stop("Unable to allocate temp memory2 in sim_params_copy.");
    to->MC_lookup_table_n = temp2;
    memcpy(to->MC_lookup_table_n, from->MC_lookup_table_n, sizeof(int) * 4);
  } else {
    to->MC_lookup_table = NULL;
    to->MC_lookup_table_n = NULL;
  }

  //to->protein_model = malloc(sizeof(model_params));
  //if (!to->protein_model) stop("Unable to allocate protein_model memory in sim_params_copy."); 




  model_params_copy(&(to->protein_model), &(from->protein_model));
  flex_params_copy(&(to->flex_params), &(from->flex_params));
}


void flex_params_copy(FLEX_params *to, FLEX_params *from){
  *to = *from;
  // output_path, outputpdb_filename, flex_cmd, flex_dir are std::string,
  // already deep-copied by "*to = *from" above.

  // filenames_to_read_in is std::vector<std::string>, already deep-copied
  // by "*to = *from" above.
}


/***********************************************************/
/****                 READING  ROUTINES                 ****/
/***********************************************************/

void flex_setup_command(FLEX_params *self){
  self->filenames_to_read_in.clear();

  self->size_of_filename_to_read_in = 0;
  self->mode = self->fromnode + (int)(rand()/(double)RAND_MAX * (self->tonode-self->fromnode+1)) % (self->tonode-self->fromnode+1);
  //for(mode = self->fromnode; mode <= self->tonode; mode++){
    self->freq = self->totalconf/2 + (int)(rand()/(double)RAND_MAX * (self->totalconf/2)) % (self->totalconf/2);
    //for(freq = self->freq; freq <= self->totalconf; freq+=self->freq){
      //self->freq

      char pos_filename[DEFAULT_LONG_STRING_LENGTH];
      char neg_filename[DEFAULT_LONG_STRING_LENGTH];
      sprintf(pos_filename,"%sRuns/Mode%02d-pos/out_froda_%08d.pdb",self->output_path.c_str(),self->mode,self->freq);
      sprintf(neg_filename,"%sRuns/Mode%02d-neg/out_froda_%08d.pdb",self->output_path.c_str(),self->mode,self->freq);
      self->filenames_to_read_in.push_back(pos_filename);
      self->filenames_to_read_in.push_back(neg_filename);
      self->size_of_filename_to_read_in += 2;

    //}

  //}

  //fprintf(stdout,"Mode: %d Freq %d\n",mode,freq); fflush(stdout);

  char flex_cmd_buf[1000];
  sprintf(flex_cmd_buf,"%s%s %s %s %d %d %lf %d %d %s",self->flex_dir.c_str(),"flex_script.sh ",self->output_path.c_str(),self->outputpdb_filename.c_str(),self->mode,self->mode,self->step, self->freq,self->freq,self->flex_dir.c_str());
  self->flex_cmd = flex_cmd_buf;

}




/* read in model parameters */
int flex_param_read(char *prm, FLEX_params *self){
  char outdir_string[256]="";
  char flex_dir[256]="";
  int k = sscanf(prm,"FLEX=%d,%d,%d,%d,%255[^,],%255[^,],%lf,%lf,%lf",&(self->number_of_processors),&(self->totalconf),&(self->fromnode),&(self->tonode),outdir_string,flex_dir,&(self->acceptance_rate_aim),&(self->hstrength_cutoff),&(self->step));

  int returnval = 5;

  //if(self->number_of_processors == 0) return 4;
  if(self->hstrength_cutoff > 1.0 || self->hstrength_cutoff <= 0.0){
    stop("FLEX hbond strength cutoff must be in (0,1]\n");
  }
  if(self->fromnode < 7 || self->tonode < 7){
     stop("Error: node numbers for FLEX must be at least 7\n");
  }
  if(self->totalconf == 0) stop("Error totalconf = 0\n");

  if(k>=5){
    self->output_path = outdir_string;
    returnval += strlen(outdir_string);

    for (int digit=1; digit<=self->totalconf; digit*=10) {
      returnval ++;
    }
    for (int digit=1; digit<=self->fromnode; digit*=10) {
      returnval ++;
    }
    for (int digit=1; digit<=self->tonode; digit*=10) {
      returnval ++;
    }
    for (int digit=1; digit<=self->number_of_processors; digit*=10) {
          returnval ++;
    }

    returnval += 4;
  }
  if(k >= 6){
    self->flex_dir = flex_dir;
	returnval += strlen(flex_dir) + 1;
  }

  flex_setup_command(self);

  return returnval;
}

void model_param_read(char *prm, /* input line */
		model_params *self,FLEX_params *nma_params) {

    int k;
    int start;
    int explicit_contact_map_file = 0;
    int found_param;
    char error_string[DEFAULT_LONG_STRING_LENGTH]="";

    if (prm == NULL) {
	fprintf(stderr,"WARNING! No command line parameter line was given.\n");
	return;
    }
    //fprintf(stderr,"%s\n",prm);

    /* Read parameters one by one */
    while (strlen(prm)>0) {

	//fprintf(stderr,"%s\n",prm);

	start = 1;
	found_param = 0;

	/* backbone H-bond */
	k=sscanf(prm, "Hbond=%lf,%lf,%lf,%lf",
	//k=sscanf(prm, "Hbond=%lf,%lf,%lf,%lf,%lf,%lf,%lf",
		&(self->hboh2),
		&(self->hbohn),
		&(self->hbcoh),
		&(self->hbs) /*,
		&(self->hboh_decay_width),
		&(self->hbohn_decay_width),
		&(self->hbcoh_decay_width) */ );
	//fprintf(stderr,"Read in %d h-bond params.\n",k);
	if (k>0) {
		found_param += 1;
		start = 6;
	}

	/* optimizing strategy*/
	k = sscanf(prm, "Opt=%d,%lf,%lf,%lf",
		&(self->opt),
		&(self->opt_totE_weight),
		&(self->opt_extE_weight),
		&(self->opt_firstlastE_weight));
	//fprintf(stderr,"Read in %d optimizing strategy totE %g.\n", self->opt, self->opt_totE_weight);

	if (k>0) {
		found_param += 1;
		start = 4;
	}


	/* contact parameters */
	k=sscanf(prm, "Contact=%lf,%lf,%lf,%lf",
		&(self->touch2),
		&(self->part),
		&(self->split),
		&(self->sts));
	//fprintf(stderr,"Read in %d contact params.\n",k);
	if (k>0) {
		found_param += 1;
		start = 8;
	}

	char gamma_string[256];
	/* gamma atoms of the model */
	k=sscanf(prm, "Gamma=%255[^,],%d,%d,%d,%d",
		gamma_string,
		&(self->use_original_gamma_atoms),
		&(self->use_3_states),
		&(self->fix_chi_angles),
		&(self->fix_CA_atoms));
	//fprintf(stderr,"Read in %d gamma params.\n",k);
	if (k>0) {
//		fprintf(stderr,"The gamma string is %s.\n",gamma_string);
		if (strcmp(gamma_string,"NONE")==0) {
			self->use_gamma_atoms = NO_GAMMA;
		} else if (strcmp(gamma_string,"LINUS_GAMMA")==0) {
			self->use_gamma_atoms = LINUS_GAMMA;
		} else if (strcmp(gamma_string,"CORRECT_GAMMA")==0) {
			self->use_gamma_atoms = CORRECT_GAMMA;
		} else if (strcmp(gamma_string,"CORRECT_KMQR_GAMMA")==0) {
			self->use_gamma_atoms = CORRECT_KMQR_GAMMA;
		} else {
			sprintf(error_string,"Unknown value for gamma model (%s).  It must be one of NONE, LINUS_GAMMA, CORRECT_GAMMA and CORRECT_KMQR_GAMMA.",gamma_string);
			stop(error_string);
		}
		if (self->use_gamma_atoms==NO_GAMMA) {
			fprintf(stderr,"WARNING! Not using gamma atoms.\n");
			fprintf(stderr,"WARNING! If you want to set an extended cutoff value for self model to make it faster, do not forget to specify if in the VDW=... parameters, otherwise the default value for the gamma atom including model will be used.\n");
		}

		found_param += 1;
		start = strlen(gamma_string)+6; //B=gamma_string...
	}

	/* vdW radii */
	char potential_string[256];
	k=sscanf(prm, "VDW=%255[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d,%lf,%lf",
		   potential_string,
		   &(self->rca),
		   &(self->rcb),
		   &(self->rc),
		   &(self->rn),
		   &(self->ro),
		   &(self->rs),
		   &(self->rel_vdw_cutoff),
		   &(self->vdw_depth_ca),
		   &(self->vdw_depth_cb),
		   &(self->vdw_depth_c),
		   &(self->vdw_depth_n),
		   &(self->vdw_depth_o),
		   &(self->vdw_depth_s),
		   &(self->vdw_uniform_depth),
		   &(self->vdw_use_extended_cutoff),
		   &(self->vdw_extended_cutoff),
		   &(self->vdw_clash_energy_at_hard_cutoff));
	if (k>0) {
		//fprintf(stderr,"%s\n",potential_string);
		if (strcmp(potential_string,"hard_cutoff")==0) {
			self->vdw_potential = HARD_CUTOFF_VDW_POTENTIAL;
			//set defaults for the hard cutoff model
			if (k<18) self->vdw_clash_energy_at_hard_cutoff = 10.4;
			if (k<14) {
				self->vdw_depth_ca = 0;
				self->vdw_depth_cb = 0;
				self->vdw_depth_c = 0;
				self->vdw_depth_n = 0;
				self->vdw_depth_o = 0;
				self->vdw_depth_s = 0;
			}
		} else if (strcmp(potential_string,"lj")==0) {
			self->vdw_potential = LJ_VDW_POTENTIAL;
		} else {
			sprintf(error_string,"Unknown value for vdW potential model (%s).  It must be one of hard_cutoff and lj.",potential_string);
			stop(error_string);
		}
		found_param += 1;
		start = strlen(potential_string) + 4;
	}

	vdw_param_calculate(self);
	//CAUTION!: aadict.c depends on params.c's model_params.  This means that
	//    initialize_sidechain_properties will have to be called after all updates
	//    of the vdW parameters; it can't be called from here, due to circular dependencies.
	//initialize_sidechain_properties(self);

	/* LJ parameters */
	//k=sscanf(prm, "LJ=%d,%d",
	//	   potenti,
	//	   &(self->vdw_lj_neighbour_hard),
	//	   &(self->vdw_lj_hbonded_hard));
	//if (k>0) {
	//	found_param += 1;
	//	start = 3;
	//}

	/* stress */
	k=sscanf(prm, "Stress=%lf,%lf",
		   &(self->stress_k),
		   &(self->stress_angle));
	if (k>0) {
		found_param += 1;
		start = 7;
	}

	/* bias potential */
	char contact_map_file[256];
	if ((k=sscanf(prm, "Bias=%255[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
		   contact_map_file,
		   &(self->bias_eta_beta),
		   &(self->bias_eta_alpha),
		   &(self->bias_kappa_alpha_3),
		   &(self->bias_kappa_alpha_4),
		   &(self->bias_kappa_beta),
		   &(self->prt),
		   &(self->bias_r_alpha),
		   &(self->bias_r_beta))) > 0) {
		self->contact_map_file = contact_map_file;
		explicit_contact_map_file = 1;

		found_param += 1;
		start = strlen(contact_map_file)+5; //B=contact_map_file...
	}
	//fprintf(stderr,"Read in %d bias params.\n",k);

	/* hydrophobic parameters */
	// 1/dist potential form
	//k=sscanf(prm, "Hydrophobic=%lf,%d,%lf",
	//	&(self->kauzmann_param),
	//	&(self->hydrophobic_min_separation),
	//	&(self->hydrophobic_max_cutoff));
	//if (k==3) self->hydrophobic_max_Eshift = 1.0 / self->hydrophobic_max_cutoff;
	// Spline potential form
	//k=sscanf(prm, "Hydrophobic=%lf,%d,%lf,%lf",
	//	&(self->kauzmann_param),
	//	&(self->hydrophobic_min_separation),
	//	&(self->hydrophobic_r),
	//	&(self->hydrophobic_half_delta));
	// Linear decay potential form (default)
	k=sscanf(prm, "Hydrophobic=%lf,%d,%lf",
		&(self->kauzmann_param),
		&(self->hydrophobic_min_separation),
		&(self->hydrophobic_cutoff_range));
	//fprintf(stderr,"Read in %d hydrophobic params.\n",k);
	if (k>0) {
		found_param += 1;
		start = 12;
	}
	
	k=sscanf(prm, "Fixit=%d",
		&(self->fixit));
	if (k>0) {
		found_param += 1;
		start = 6;
	}

	/* electrostatic parameters */
	k=sscanf(prm, "Elec=%lf,%lf,%d",
		&(self->recip_dielectric_param),
		&(self->debye_length_param),
		&(self->electrostatic_min_separation));
	//fprintf(stderr,"Read in %d electrostatic params.\n",k);
	if (k>0) {
		found_param += 1;
		start = 5;
	}

	/* side chain hydrogen bond parameters */
	k=sscanf(prm, "SchHbond=%lf,%lf,%lf,%lf,%lf,%d",
		&(self->sidechain_hbond_strength_s2b),
		&(self->sidechain_hbond_strength_b2s),
		&(self->sidechain_hbond_strength_s2s),
		&(self->sidechain_hbond_angle_cutoff),
		&(self->sidechain_hbond_decay_width),
		&(self->sidechain_hbond_min_separation));
	//fprintf(stderr,"Read in %d sidechain h-bond params.\n",k);
	if (k>0) {
		found_param += 1;
		start = 9;
	}

	/* secondary radius of gyration parameters */
	k=sscanf(prm, "Rgyr=%lf,%lf,%lf,%lf",
		&(self->srgy_param),
		&(self->srgy_offset),
		&(self->hphobic_srgy_param),
		&(self->hphobic_srgy_offset));
	//fprintf(stderr,"Read in %d secondary radius of gyration params.\n",k);
	if (k>0) {
		found_param += 1;
		start = 5;

	}

	/* S-S bond */
	k=sscanf(prm,"SSbond=%lf,%lf,%lf,%lf",
		&(self->Sbond_strength),
		&(self->Sbond_distance),
		&(self->Sbond_cutoff),
		&(self->Sbond_dihedral_cutoff));
	if (k>0) {
		found_param += 1;
		start = 7;
	}

	/* fixed amino acids */
	char fixed_file[256];
	if ((k=sscanf(prm, "fixed=%255[^,]",fixed_file)) == 1) {
		fprintf(stderr,"setting fixed amino acid list file %s\n",fixed_file);
		self->fixed_aalist_file = fixed_file;
		found_param += 1;
		start = strlen(fixed_file) + strlen("fixed=");
	}


	/* external potential */
	int xdir = 0;
	int ydir = 0;
	int zdir = 0;
	char constraint_file[256];
	k=sscanf(prm,"external=%d",&(self->external_potential_type));
	if (k>0) {
		fprintf(stderr,"found external potential parameters %d \n", k);
		int l;
		if (self->external_potential_type == 1) {
			l=sscanf(prm+9,"%d,%255[^,],%lf,%lf",
				&(self->external_potential_type), constraint_file, &(self->external_k[0]), &(self->external_r0[0])); //save r0, k into 1st vector element
			fprintf(stderr, "found external potential parameters %d \n", l);
			//stop("unimplemented");
			if (l>1) {
				fprintf(stderr,"setting constrained list file %s\n",constraint_file);
				self->external_constrained_aalist_file = constraint_file;
			}

		} else if (self->external_potential_type == 2) { /* unimplemented */
			stop("unimplemented");
			l=sscanf(prm+9,"%d,%d,%lf,%lf,%d,%lf,%lf,%d,%lf,%lf",
				&(self->external_potential_type), &xdir, &(self->external_k[0]), &(self->external_r0[0]),
				&ydir, &(self->external_k[1]), &(self->external_r0[1]),
				&zdir, &(self->external_k[2]), &(self->external_r0[2]));
			if (l>0) {
				/* check x direction */
				if (xdir == -1) self->external_direction[0] = EXTERNAL_NEGATIVE;
				else if (xdir == 1) self->external_direction[0] = EXTERNAL_POSITIVE;
				else if (xdir == 2) self->external_direction[0] = EXTERNAL_POSNEG;
				else if (xdir == 0) self->external_direction[0] = EXTERNAL_NONE;
				else stop("external direction has to be one of -1(negative), 0(none), 1(positive) or 2 (positive and negative)");
				/* check y direction */
				if (ydir == -1) self->external_direction[1] = EXTERNAL_NEGATIVE;
				else if (ydir == 1) self->external_direction[1] = EXTERNAL_POSITIVE;
				else if (ydir == 2) self->external_direction[1] = EXTERNAL_POSNEG;
				else if (ydir == 0) self->external_direction[1] = EXTERNAL_NONE;
				else stop("external direction has to be one of -1(negative), 0(none), 1(positive) or 2 (positive and negative)");
				/* check z direction */
				if (zdir == -1) self->external_direction[2] = EXTERNAL_NEGATIVE;
				else if (zdir == 1) self->external_direction[2] = EXTERNAL_POSITIVE;
				else if (zdir == 2) self->external_direction[2] = EXTERNAL_POSNEG;
				else if (zdir == 0) self->external_direction[2] = EXTERNAL_NONE;
				else stop("external direction has to be one of -1(negative), 0(none), 1(positive) or 2 (positive and negative)");

			}
		} else if (self->external_potential_type == 3) { /* conical potential */
			l=sscanf(prm+9,"%d,%255[^,],%lf,%lf,%lf",
				&(self->external_potential_type), constraint_file, &(self->external_k[0]), &(self->external_r0[0]),&(self->external_ztip)); //save r0, k into 1st vector element
			if (l==5) {
				fprintf(stderr,"setting constrained list file %s\n",constraint_file);
				self->external_constrained_aalist_file = constraint_file;
			} else {
				stop("conical potential needs 4 parameters (constraint list file,k,r0,ztip)");
			}
		}
		else if (self->external_potential_type == 5) {
			l = sscanf(prm + 9, "%d,%255[^,],%lf,%lf",
				&(self->external_potential_type), constraint_file, &(self->external_k[0]), &(self->external_r0[0])); //save r0, k into 1st vector element
			fprintf(stderr, "found external potential parameters %d \n", l);
			//stop("unimplemented");
			if (l > 1) {
				fprintf(stderr, "setting constrained list file %s\n", constraint_file);
				self->external_constrained_aalist_file = constraint_file;
			}
		} else {
			stop("unknown external potential type, has to be 1 (cylindrical around z) or 3 (conical around z).");
		}
		found_param += 1;
		start = strlen(constraint_file) + 9 + 2; //2 comes from the type and the comma after
		//fprintf(stderr,"setting found_param to %d and start to %d\n",found_param,start);
	}
	/* external potential 2 */
	char constraint_file2[256];
	k=sscanf(prm,"external2=%d",&(self->external_potential_type2));
	if (k>0) {
		fprintf(stderr,"found external2 potential parameters %d\n", self->external_potential_type2);
		int l;
		if (self->external_potential_type2 == 1) {
			l=sscanf(prm+10,"%d,%255[^,],%lf,%lf",
				&(self->external_potential_type2), constraint_file2, &(self->external_k2[0]), &(self->external_r02[0])); //save r0, k into 1st vector element
			if (l>1) {
				fprintf(stderr,"setting constrained list file %s\n",constraint_file2);
				self->external_constrained_aalist_file2 = constraint_file2;
			}
		} else if (self->external_potential_type2 == 3) { /* conical potential */
			l=sscanf(prm+9,"%d,%255[^,],%lf,%lf,%lf",
				&(self->external_potential_type2), constraint_file2, &(self->external_k2[0]), &(self->external_r02[0]),&(self->external_ztip2)); //save r0, k into 1st vector element
			if (l==5) {
				fprintf(stderr,"setting constrained list file %s\n",constraint_file2);
				self->external_constrained_aalist_file2 = constraint_file2;
			} else {
				stop("conical potential needs 4 parameters (constraint list file,k,r0,ztip)");
			}
		} else if (self->external_potential_type2 == 4) {
			l = sscanf(prm + 10, "%d,%255[^,],%lf,%lf", &(self->external_potential_type2), constraint_file2, &(self->external_k2[0]), &(self->external_r02[0]));
			fprintf(stderr, "setting cyclic peptide %s \n", constraint_file2);
		} else {
			stop("unknown external2 potential type, has to be 1.");
		}
		found_param += 1;
		start = strlen(constraint_file2) + 10 + 2; //2 comes from the type and the comma after
		//fprintf(stderr,"setting found_param to %d and start to %d\n",found_param,start);
	}

    /*dummy position for FLEX */
	char temp[DEFAULT_LONG_STRING_LENGTH];
	k = sscanf(prm,"FLEX=%s",temp);
	if(k>0){
	  start = flex_param_read(prm, nma_params);
	  found_param += 1;

	}



	/* Check for more than one parameters found */
	if (found_param > 1) {
		sprintf(error_string,"The same section of the parameter line was used to read in more than one parameters.  Prm: %s",prm);
		stop(error_string);
	}

	/* Check for unknown parameters */
	if (found_param == 0 ) {
		sprintf(error_string,"Unknown parameter in prm string: %s",prm);
		stop(error_string);
	}

	/* Shift the prm to the next parameter entry */
	int next = 0;
	for (int i=start; (next==0 && i<=strlen(prm)); i++) {
		if (isalpha(prm[i])) next = i;
	}
	if (next > 0 && next < strlen(prm)) {
		/* the new prm is >=1 long */
		prm = prm+next;
	} else {
		//model_param_print(*self);
		if (!explicit_contact_map_file) fprintf(stderr,"WARNING: No contact map specified.");
		return;
	}

    }

}




/* calculating vdW depths and shifts at the cutoffs for the atomic pairs,
   not to have to recalculate on the fly */
void vdw_param_calculate(model_params *self) {

   /* depths */

   // quick check that nothing is negative
   if (self->vdw_depth_ca < 0) self->vdw_depth_ca = 0;
   if (self->vdw_depth_cb < 0) self->vdw_depth_cb = 0;
   if (self->vdw_depth_c  < 0) self->vdw_depth_c  = 0;
   if (self->vdw_depth_n  < 0) self->vdw_depth_n  = 0;
   if (self->vdw_depth_o  < 0) self->vdw_depth_o  = 0;
   if (self->vdw_depth_s  < 0) self->vdw_depth_s  = 0;

   if (self->vdw_uniform_depth) {
   self->vdw_depth_ca_sqrt = self->vdw_depth_cb_sqrt =
			     self->vdw_depth_c_sqrt  =
			     self->vdw_depth_n_sqrt  =
			     self->vdw_depth_o_sqrt  =
			     self->vdw_depth_s_sqrt  = sqrt(self->vdw_depth_ca);
   self->vdw_depth_ca_ca = self->vdw_depth_cb_cb =
			   self->vdw_depth_c_c   =
			   self->vdw_depth_n_n   =
			   self->vdw_depth_o_o   =
			   self->vdw_depth_ca_cb =
			   self->vdw_depth_ca_c  =
			   self->vdw_depth_ca_n  =
			   self->vdw_depth_ca_o  =

			   self->vdw_depth_cb_c  =
			   self->vdw_depth_cb_n  =
			   self->vdw_depth_cb_o  =

			   self->vdw_depth_c_n   =
			   self->vdw_depth_c_o   =

			   self->vdw_depth_n_o   = self->vdw_depth_ca;
   /* shifts */

   double shift = pow(1/self->rel_vdw_cutoff,12) - 2 * pow(1/self->rel_vdw_cutoff,6);
   self->vdw_shift = shift;

   self->vdw_Eshift_ca_ca = self->vdw_Eshift_cb_cb =
			    self->vdw_Eshift_c_c   =
			    self->vdw_Eshift_n_n   =
			    self->vdw_Eshift_o_o   =
			    self->vdw_Eshift_ca_cb =
			    self->vdw_Eshift_ca_c  =
			    self->vdw_Eshift_ca_n  =
			    self->vdw_Eshift_ca_o  =
			    self->vdw_Eshift_cb_c  =
			    self->vdw_Eshift_cb_n  =
			    self->vdw_Eshift_cb_o  =
			    self->vdw_Eshift_c_n   =
			    self->vdw_Eshift_c_o   =
			    self->vdw_Eshift_n_o   = self->vdw_depth_ca_ca * shift;
   } else {
   self->vdw_depth_ca_sqrt = sqrt(self->vdw_depth_ca);
   self->vdw_depth_cb_sqrt = sqrt(self->vdw_depth_cb);
   self->vdw_depth_c_sqrt = sqrt(self->vdw_depth_c);
   self->vdw_depth_n_sqrt = sqrt(self->vdw_depth_n);
   self->vdw_depth_o_sqrt = sqrt(self->vdw_depth_o);
   self->vdw_depth_s_sqrt = sqrt(self->vdw_depth_s);
  
   self->vdw_depth_ca_ca = self->vdw_depth_ca;
   self->vdw_depth_cb_cb = self->vdw_depth_cb;
   self->vdw_depth_c_c   = self->vdw_depth_c;
   self->vdw_depth_n_n   = self->vdw_depth_n;
   self->vdw_depth_o_o   = self->vdw_depth_o;

   self->vdw_depth_ca_cb = sqrt(self->vdw_depth_ca * self->vdw_depth_cb);
   self->vdw_depth_ca_c  = sqrt(self->vdw_depth_ca * self->vdw_depth_c );
   self->vdw_depth_ca_n  = sqrt(self->vdw_depth_ca * self->vdw_depth_n );
   self->vdw_depth_ca_o  = sqrt(self->vdw_depth_ca * self->vdw_depth_o );

   self->vdw_depth_cb_c  = sqrt(self->vdw_depth_cb * self->vdw_depth_c );
   self->vdw_depth_cb_n  = sqrt(self->vdw_depth_cb * self->vdw_depth_n );
   self->vdw_depth_cb_o  = sqrt(self->vdw_depth_cb * self->vdw_depth_o );

   self->vdw_depth_c_n   = sqrt(self->vdw_depth_c  * self->vdw_depth_n );
   self->vdw_depth_c_o   = sqrt(self->vdw_depth_c  * self->vdw_depth_o );

   self->vdw_depth_n_o   = sqrt(self->vdw_depth_n  * self->vdw_depth_o );

   /* shifts */

   double shift = pow(1/self->rel_vdw_cutoff,12) - 2 * pow(1/self->rel_vdw_cutoff,6);
   self->vdw_shift = shift;

   self->vdw_Eshift_ca_ca = self->vdw_depth_ca_ca * shift;
   self->vdw_Eshift_cb_cb = self->vdw_depth_cb_cb * shift;
   self->vdw_Eshift_c_c   = self->vdw_depth_c_c   * shift;
   self->vdw_Eshift_n_n   = self->vdw_depth_n_n   * shift;
   self->vdw_Eshift_o_o   = self->vdw_depth_o_o   * shift;

   self->vdw_Eshift_ca_cb = self->vdw_depth_ca_cb * shift;
   self->vdw_Eshift_ca_c  = self->vdw_depth_ca_c  * shift;
   self->vdw_Eshift_ca_n  = self->vdw_depth_ca_n  * shift;
   self->vdw_Eshift_ca_o  = self->vdw_depth_ca_o  * shift;

   self->vdw_Eshift_cb_c  = self->vdw_depth_cb_c  * shift;
   self->vdw_Eshift_cb_n  = self->vdw_depth_cb_n  * shift;
   self->vdw_Eshift_cb_o  = self->vdw_depth_cb_o  * shift;

   self->vdw_Eshift_c_n   = self->vdw_depth_c_n   * shift;
   self->vdw_Eshift_c_o   = self->vdw_depth_c_o   * shift;

   self->vdw_Eshift_n_o   = self->vdw_depth_n_o   * shift;

  }

}


/***********************************************************/
/****                 PRINTING ROUTINES                 ****/
/***********************************************************/


/* Printing the vdW backbone cutoff and the gamma - gamma and gamma - nongamma cutoff matrices */
void print_vdw_cutoff_distances(model_params *mod_params, FILE *outfile) {

	int l, m;

	fprintf(outfile,"vdW backbone cutoff: %g\n",mod_params->vdw_backbone_cutoff);

	if (mod_params->use_gamma_atoms == NO_GAMMA) {
	    fprintf(outfile,"vdW cutoff matrices are not used (no gamma atoms in the model).\n");
	    return;
	}

	if (mod_params->vdw_gamma_gamma_cutoff == NULL) {
	    fprintf(outfile,"vdW gamma - gamma cutoff was not allocated.\n");
	    return;
	}
	if (mod_params->vdw_gamma_nongamma_cutoff == NULL) {
	    fprintf(outfile,"vdW gamma - gamma cutoff was not allocated.\n");
	    return;
	}

	/* Print the cutoff matrices */
	/* The 26x26 array of gamma - gamma contact cutoffs */
	fprintf(outfile,"const double maxvdw_gamma_gamma[702] ={ ");
	for(l = 0; l < 27; l++){
	    for(m = 0; m < 26; m++){
		fprintf(outfile,"%f, ",mod_params->vdw_gamma_gamma_cutoff[l*26+m]);
	    }
	    fprintf(outfile,"\n");
	}
	fprintf(outfile," };\n");

	/* The 26x26 array of gamma - nongamma contact cutoffs */
	fprintf(outfile,"const double maxvdw_gamma_nongamma[702] ={ ");
	for(l = 0; l < 27; l++){
	    for(m = 0; m < 26; m++){
		fprintf(outfile,"%f, ",mod_params->vdw_gamma_nongamma_cutoff[l*26+m]);
	    }
	    fprintf(outfile,"\n");
	}
	fprintf(outfile," };\n");

}

/* print model parameters */
void model_param_print(model_params self, FILE *outfile) {

  fprintf(outfile,"================MODEL=PARAMETERS===================\n");
  /* gamma atoms */
  fprintf(outfile,"use gamma atoms?(no:%d,LINUS:%d,correct:%d,LINUS+correct_KMQR:%d) %d\n",NO_GAMMA,LINUS_GAMMA,CORRECT_GAMMA,CORRECT_KMQR_GAMMA,self.use_gamma_atoms);
  fprintf(outfile,"use original sidechain dihedrals when reading (1), or pick a random one(0)? %d\n",self.use_original_gamma_atoms);
  fprintf(outfile,"adjust to 3 states when reading in? %d\n",self.use_3_states);
  fprintf(outfile,"fix the side chain dihedral angles in the mc? %d\n",self.fix_chi_angles);
  fprintf(outfile,"fix CA atoms in the mc (only pivot peptide bonds)? %d\n",self.fix_CA_atoms);
  /* atomic radii */
  fprintf(outfile,"VDW RADII\n");
  fprintf(outfile,"r(CA_) %g\n",self.rca);
  fprintf(outfile,"r(CB_) %g\n",self.rcb);
  fprintf(outfile,"r(C__) %g\n",self.rc);
  fprintf(outfile,"r(O__) %g\n",self.ro);
  fprintf(outfile,"r(N__) %g\n",self.rn);
  fprintf(outfile,"r(S__) %g\n",self.rs);
  fprintf(outfile,"rel. cutoff %g\n",self.rel_vdw_cutoff);
  fprintf(outfile,"uniform depth? %d\n",self.vdw_uniform_depth);
  fprintf(outfile,"eps(CA_) %g\n",self.vdw_depth_ca);
  fprintf(outfile,"eps(CB_) %g\n",self.vdw_depth_cb);
  fprintf(outfile,"eps(C__) %g\n",self.vdw_depth_c);
  fprintf(outfile,"eps(N__) %g\n",self.vdw_depth_n);
  fprintf(outfile,"eps(O__) %g\n",self.vdw_depth_o);
  fprintf(outfile,"eps(S__) %g\n",self.vdw_depth_s);
  fprintf(outfile,"vdW potential model (%d: hard cutoff, %d: LJ) %d\n",HARD_CUTOFF_VDW_POTENTIAL,LJ_VDW_POTENTIAL,self.vdw_potential);
  fprintf(outfile,"Use hard cutoff part of LJ for neighbouring amino acid pairs? (0:no, 1:yes) %d\n",self.vdw_lj_neighbour_hard);
  fprintf(outfile,"Use hard cutoff part of LJ for H-bonded amino acid pairs (and their neighbours)? (0:no, 1:yes) %d\n",self.vdw_lj_hbonded_hard);
  fprintf(outfile,"Use and extended vdW cutoff? (0:no,1:yes): %d\n",self.vdw_use_extended_cutoff);
  fprintf(outfile,"vdW function pointer allocated? (0:no,1:yes) %d",(self.vdw_function!=NULL));
  fprintf(outfile,"vdW clash energy at hard cutoff %g",self.vdw_clash_energy_at_hard_cutoff);
  fprintf(outfile,"extended vdW cutoff value %g\n",self.vdw_extended_cutoff);
  print_vdw_cutoff_distances(&self,outfile);
  /* stress */
  fprintf(outfile,"STRESS\n");
  fprintf(outfile,"stress_k %g\n",self.stress_k);
  fprintf(outfile,"stress_angle %g\n",self.stress_angle);

  /* hydrogen bond */
  fprintf(outfile,"H-BOND\n");
  fprintf(outfile,"OH cutoff square %g\n",self.hboh2);
  fprintf(outfile,"OHN angle cutoff %g\n",self.hbohn);
  fprintf(outfile,"COH angle cutoff %g\n",self.hbcoh);
  //fprintf(outfile,"OH cutoff decay width %g\n",self.hboh_decay_width);
  //fprintf(outfile,"OHN angle decay width %g\n",self.hbohn_decay_width);
  //fprintf(outfile,"COH angle decay width %g\n",self.hbcoh_decay_width);
  fprintf(outfile,"strength %g\n",self.hbs);
  /* contact parameters */
  fprintf(outfile,"CONTACT\n");
  fprintf(outfile,"cutoff square %g\n",self.touch2);
  fprintf(outfile,"part %g\n",self.part);
  fprintf(outfile,"split %g\n",self.split);
  fprintf(outfile,"sts %g\n",self.sts);
  /* biasing force constants */
  fprintf(outfile,"BIAS\n");
  fprintf(outfile,"contact map file %s\n",self.contact_map_file.c_str());
  fprintf(outfile,"eta_beta %g\n",self.bias_eta_beta);
  fprintf(outfile,"eta_alpha %g\n",self.bias_eta_alpha);
  fprintf(outfile,"kappa_alpha_3 %g\n",self.bias_kappa_alpha_3);
  fprintf(outfile,"kappa_alpha_4 %g\n",self.bias_kappa_alpha_4);
  fprintf(outfile,"kappa_beta %g\n",self.bias_kappa_beta);
  fprintf(outfile,"CA(0) or CB(1) %g\n",self.prt);
  fprintf(outfile,"r0_a %g\n",self.bias_r_alpha);
  fprintf(outfile,"r0_b %g\n",self.bias_r_beta);
  /* hydrophobicity */
  fprintf(outfile,"HYDROPHOBICITY\n");
  fprintf(outfile,"k_h %g\n",self.kauzmann_param);
  fprintf(outfile,"delta_decay %g\n",self.hydrophobic_cutoff_range);
  fprintf(outfile,"from i,i+%d\n",self.hydrophobic_min_separation);
  //fprintf(outfile,"hydrophobic_min_cutoff (if 1/dist) %g\n",self.hydrophobic_min_cutoff);
  //fprintf(outfile,"hydrophobic_max_cutoff (if 1/dist) %g\n",self.hydrophobic_max_cutoff);
  //fprintf(outfile,"hydrophobic_max_Eshift (if 1/dist) %g\n",self.hydrophobic_max_Eshift);
  //fprintf(outfile,"hydrophobic_r (if spline) %g\n",self.hydrophobic_r);
  //fprintf(outfile,"hydrophobic_half_delta (if spline) %g\n",self.hydrophobic_half_delta);
  /* electrostatics */
  fprintf(outfile,"ELECTROSTATICS\n");
  fprintf(outfile,"permittivity %g\n",self.recip_dielectric_param);
  fprintf(outfile,"Debye length %g\n",self.debye_length_param);
  fprintf(outfile,"from i,i+%d\n",self.electrostatic_min_separation);
  /* side chain hydrogen bond parameters */
  fprintf(outfile,"SIDE CHAIN H-BOND\n");
  fprintf(outfile,"H_s2b %g\n",self.sidechain_hbond_strength_s2b);
  fprintf(outfile,"H_b2s %g\n",self.sidechain_hbond_strength_b2s);
  fprintf(outfile,"H_s2s %g\n",self.sidechain_hbond_strength_s2s);
  fprintf(outfile,"cutoff %g\n",self.sidechain_hbond_cutoff);
  fprintf(outfile,"cos(angle cutoff) %g\n",self.sidechain_hbond_angle_cutoff);
  fprintf(outfile,"delta_decay %g\n",self.sidechain_hbond_decay_width);
  fprintf(outfile,"from i,i+%d\n",self.sidechain_hbond_min_separation);
  /* secondary radius of gyration */
  fprintf(outfile,"srgy_param %g\n",self.srgy_param);
  fprintf(outfile,"srgy_offset %g\n",self.srgy_offset);
  fprintf(outfile,"hphobic_srgy_param %g\n",self.hphobic_srgy_param);
  fprintf(outfile,"hphobic_srgy_offset %g\n",self.hphobic_srgy_offset);
  /* SS bond */
  fprintf(outfile,"S-bond strength %g\n", self.Sbond_strength);
  fprintf(outfile,"S-bond distance %g\n", self.Sbond_distance);
  fprintf(outfile,"S-bond cutoff %g\n",	self.Sbond_cutoff);
  fprintf(outfile,"S-bond dihedral cutoff %g\n",	self.Sbond_dihedral_cutoff);
  fprintf(outfile,"C-S-S angle hard coded\n");  
  fprintf(outfile,"FIXED AMINO ACIDS\n");
  fprintf(outfile,"fixed amino acid list file: %s\n",self.fixed_aalist_file.c_str());
  fprintf(outfile,"EXTERNAL CONSTRAINTS\n");
  fprintf(outfile,"constrained amino acid list file: %s\n",self.external_constrained_aalist_file.c_str());
  fprintf(outfile,"external type: %d\n",self.external_potential_type);
  fprintf(outfile,"external_k[0]: %g\n",self.external_k[0]);
  fprintf(outfile,"external_r0[0]: %g\n",self.external_r0[0]);
  if (self.external_potential_type == 3) fprintf(outfile,"external_ztip: %g\n",self.external_ztip);
  fprintf(outfile,"constrained amino acid list file(2): %s\n",self.external_constrained_aalist_file2.c_str());
  fprintf(outfile,"external type(2): %d\n",self.external_potential_type2);
  fprintf(outfile,"external_k[0](2): %g\n",self.external_k2[0]);
  fprintf(outfile,"external_r0[0](2): %g\n",self.external_r02[0]);
  if (self.external_potential_type2 == 3) fprintf(outfile,"external_ztip(2): %g\n",self.external_ztip2);
  fprintf(outfile,"============END=MODEL=PARAMETERS===================\n");

}


void flex_param_print(FLEX_params self, FILE *outfile){
  fprintf(outfile,"================FLEX=PARAMETERS===================\n");
  fprintf(outfile,"Number of processors: %d\n",self.number_of_processors);
  if(self.number_of_processors == 0) return;
  fprintf(outfile,"Modes: %d to %d\n",self.fromnode,self.tonode);
  fprintf(outfile,"Hbond strength for FLEX: %lf\n",self.hstrength_cutoff);
  if(self.only_bias_hbonds == 1.0){fprintf(outfile,"Using Only Bias H-bonds\n");}
  else{fprintf(outfile,"Using All H-bonds\n");}
  fprintf(outfile,"Max Totconf: %d Step: %lf\n",self.totalconf, self.step);
  fprintf(outfile,"Output path: %s FLEX dir: %s\n",self.output_path.c_str(),self.flex_dir.c_str());
  fprintf(outfile, "Acceptance Rate Aim %lf\n",self.acceptance_rate_aim);
  fprintf(outfile,"FLEX command: %s\n",self.flex_cmd.c_str());
  fprintf(outfile,"============END=FLEX=PARAMETERS===================\n");
}

/* print simulation parameters (including model_parameters) */
void param_print(simulation_params self, FILE *outfile) {

  fprintf(outfile,"==============SIMULATION=PARAMETERS================\n");
  /* general simulation */
  fprintf(outfile,"--GENERAL--\n");
  fprintf(outfile,"infile_name %s\n",self.infile_name);
  fprintf(outfile,"outfile_name %s\n",self.outfile_name);
  fprintf(outfile,"pace %d\n",self.pace);
  fprintf(outfile,"stretch %d\n",self.stretch);
  fprintf(outfile,"test mask %x\n",self.tmask);
  fprintf(outfile,"random seed %d\n",self.seed);
  fprintf(outfile,"parameters %s\n",self.prm);
  fprintf(outfile,"acceptance ratio %g\n",self.acceptance_rate);
  fprintf(outfile,"amplitude %g\n",self.amplitude);
  fprintf(outfile,"keep_amplitude_fixed %d\n",self.keep_amplitude_fixed);
  fprintf(outfile,"acceptance %g\n",self.acceptance);
  fprintf(outfile,"accept_counter %d\n",self.accept_counter);
  fprintf(outfile,"reject_counter %d\n",self.reject_counter);
  fprintf(outfile,"acceptance_rate_tolerance %g\n",self.acceptance_rate_tolerance);
  fprintf(outfile,"amplitude_changing_factor %g\n",self.amplitude_changing_factor);
  if (self.MC_lookup_table != NULL && self.MC_lookup_table_n != NULL && self.Nchains > 0) {
    for (int i=0; i<4; i++) {
      for (int j=0; j<(self.NAA-1+self.Nchains);j++) {
        fprintf(outfile,"%d ",self.MC_lookup_table[i*(self.NAA+self.Nchains)+j]);
      }
      fprintf(outfile, " (%d valid elements)\n",self.MC_lookup_table_n[i]);
    }
  }
  /* peptide or multi-chain protein */
  fprintf(outfile,"--PEPTIDE--\n");
  fprintf(outfile,"number of chains: %d\n", self.Nchains);
  if (self.seq == NULL || self.sequence == NULL) {
    fprintf(outfile,"amino acid sequence %s\n",self.seq);
    fprintf(outfile,"amino acid sequence %s\n",self.sequence);
  } else {
    if (strlen(self.seq) == strlen(self.sequence)) { // 1 chain
      fprintf(outfile,"amino acid sequence %s\n",self.seq);
    } else {
      fprintf(outfile,"amino acid sequence %s\n",self.sequence);
      fprintf(outfile,"number of chains: %d\n", (int)strlen(self.sequence) - (int)strlen(self.seq));
    }
  }
  fprintf(outfile,"number of amino acids %d\n",self.NAA);
  /* thermodynamic beta */
  fprintf(outfile,"--TEMPERATURE--\n");
  fprintf(outfile,"thermobeta %g\n",self.thermobeta);
  fprintf(outfile,"lowtemp %d\n",self.lowtemp);
  fprintf(outfile,"beta1 %g\n",self.beta1);
  fprintf(outfile,"beta2 %g\n",self.beta2);
  fprintf(outfile,"bstp %g\n",self.bstp);
  fprintf(outfile,"intrvl %d\n",self.intrvl);
  fprintf(outfile, "nswap_per_try %d\n", self.nswap_per_try);
  /* nested sampling */
  fprintf(outfile,"--NESTED SAMPLING--\n");
  fprintf(outfile,"do nested sampling? %d\n",self.NS);
  fprintf(outfile,"NS iter %d\n",self.iter);
  fprintf(outfile,"NS iter_start %d\n",self.iter_start);
  fprintf(outfile,"NS iter_max %d\n",self.iter_max);
  fprintf(outfile,"NS logLstar %g\n",self.logLstar);
  fprintf(outfile,"NS logZ %g\n",self.logZ);
  fprintf(outfile,"NS logfactor %g\n",self.logfactor);
  fprintf(outfile,"NS alpha %g\n",self.alpha);
  fprintf(outfile,"NS logX %g\n",self.logX);
  fprintf(outfile,"NS logX_start %g\n",self.logX_start);
  fprintf(outfile,"NS Delta_logX %g\n",self.Delta_logX);
  fprintf(outfile,"NS log_DeltaX %g\n",self.log_DeltaX);
  fprintf(outfile,"NS H %g\n",self.H);
  fprintf(outfile,"NS N %d\n",self.N);
  /* checkpointing */
  fprintf(outfile,"--CHECKPOINTING--\n");
  fprintf(outfile,"num_NS_per_checkpoint %d\n",self.num_NS_per_checkpoint);
  fprintf(outfile,"checkpoint filename %s\n",self.checkpoint_filename);
  fprintf(outfile,"checkpoint_counter %d\n",self.checkpoint_counter);
  fprintf(outfile,"restart_from_checkpoint %d\n",self.restart_from_checkpoint);
  fprintf(outfile,"checkpoint %d\n",self.checkpoint);

  model_param_print(self.protein_model, outfile);
  flex_param_print(self.flex_params,outfile);
  fprintf(outfile,"==========END=SIMULATION=PARAMETERS================\n");

}
