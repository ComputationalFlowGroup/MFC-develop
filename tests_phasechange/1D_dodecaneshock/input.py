#!/usr/bin/python

# Dependencies and Logistics ===================================================

import math, json
eps = 1E-8
tfinal = 473E-6
step_stop = int(40000)
step_save = int(step_stop/50)
dt = float(tfinal/float(step_stop))
# print "dt=",dt

# Case Analysis Configuration ==================================================
# Configuring case dictionary
print(json.dumps({                                                              
                    # Logistics ================================================
                    'case_dir'                     : '\'.\'',                  
                    'run_time_info'                : 'T',                      
                    # ==========================================================
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 : 0.E+00,                   
                    'x_domain%end'                 : 1.E+00,                   
                    'm'                            : 999,                     
                    'n'                            : 0,                        
                    'p'                            : 0,                        
                    'dt'                           : dt,                       
                    't_step_start'                 : 0,                        
                    't_step_stop'                  : step_stop,                
                    't_step_save'                  : step_save,                
		    # ==========================================================
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 2,                        
                    'model_eqns'                   : 3,                        
                    'alt_soundspeed'               : 'F',                      
                    'num_fluids'                   : 2,                        
		    'adv_alphan'                   : 'T',                      
		    'mpp_lim'                      : 'T',                      
		    'mixture_err'                  : 'T',                      
		    'time_stepper'                 : 3,                        
                    'weno_vars'                    : 2,                        
                    'weno_order'                   : 5,                        
                    'weno_eps'                     : 1.E-016,                  
                    'char_decomp'                  : 'F',                      
                    'mapped_weno'                  : 'F',                      
                    'null_weights'                 : 'F',                      
                    'mp_weno'                      : 'F',                      
		    'riemann_solver'               : 2,                        
                    'wave_speeds'                  : 1,                        
                    'avg_state'                    : 2,                        
                    'commute_err'                  : 'F',                      
                    'split_err'                    : 'F',                      
                    'bc_x%beg'                     : -3,                       
                    'bc_x%end'                     : -3,                       
                    'relax_model'                  :  3,                      
                    'palpha_eps'                   : 1E-12,
                    'ptgalpha_eps'                 : 1E-4, 
                    # ==========================================================
                                                                               
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        
                    'precision'                    : 1,                        
                    'prim_vars_wrt'                :'T',                       
		    'parallel_io'                  :'F',                       
                    # ==========================================================
                                                                                
		    # Patch 1: High pressured water ============================
                    'patch_icpp(1)%geometry'       : 1,                        
                    'patch_icpp(1)%x_centroid'     : 0.5E+00,                  
                    'patch_icpp(1)%length_x'       : 1.E+00,                   
                    'patch_icpp(1)%vel(1)'         : 0.E+00,                   
                    'patch_icpp(1)%pres'           : 1.E+08,                   
                    'patch_icpp(1)%alpha_rho(1)'   : 500.E+00*(1.-eps),        
                    'patch_icpp(1)%alpha_rho(2)'   : 2.E+00*eps,               
                    'patch_icpp(1)%alpha(1)'       : 1.-eps,                   
                    'patch_icpp(1)%alpha(2)'       : eps,                      
                    # ==========================================================

                    # Patch 2: Air bubble ======================================
                    'patch_icpp(2)%geometry'       : 1,                        
                    'patch_icpp(2)%x_centroid'     : 0.875E+00,                
                    'patch_icpp(2)%length_x'       : 0.25E+00,                 
                    'patch_icpp(2)%alter_patch(1)' : 'T',                      
                    'patch_icpp(2)%vel(1)'         : 0.E+00,                   
                    'patch_icpp(2)%pres'           : 1.E+05,                   
                    'patch_icpp(2)%alpha_rho(1)'   : 500.E+00*eps,             
                    'patch_icpp(2)%alpha_rho(2)'   : 2.E+00*(1.-eps),          
                    'patch_icpp(2)%alpha(1)'       : eps,                      
                    'patch_icpp(2)%alpha(2)'       : 1.-eps,                   
                    # ==========================================================
 
		    # Fluids Physical Parameters ===============================
                    'fluid_pp(1)%gamma'            : 1.E+00/(2.35E+00-1.E+00), 
                    'fluid_pp(1)%pi_inf'           : 2.35E+00*4.E+08/(2.35E+00-1.E+00), 
                    'fluid_pp(1)%cv'               : 1077.7E+00,               
                    'fluid_pp(1)%qv'           	   : -775.269E+03,             
                    'fluid_pp(1)%qvp'              : 0.E+00,                   
                    'fluid_pp(2)%gamma'            : 1.E+00/(1.025E+00-1.E+00),
                    'fluid_pp(2)%pi_inf'           : 0.E+00,                   
                    'fluid_pp(2)%cv'               : 1956.45E+00,              
                    'fluid_pp(2)%qv'               : -237.547E+03,             
                    'fluid_pp(2)%qvp'              : -24.4E+03,                
	            # ==========================================================

}))
# ==============================================================================
