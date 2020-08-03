classdef VCA
  % VCA 
  
  properties    
  end
  
  methods (Static)	
	
	% Voltage sensitivities 
	[S_out] = makeS( mpw, mode, incr);
	[S_out] = cropS( S, SET_FROM, SET_TO)
	[dVldVg_out, dQgdVg_out, dVldQl_out, dVldQg_out, JAC_out, Ybus_out] = csvc_sen(mpw, debug)
	[plt_set, plt_qlt] = pilots_conejo1993(S_GL, S_LL, C_LL, Q_x, plt_ini, plt_max, debug)

	% Clustering
	[LBL, mean_rms_dst] = vcs(S, link, param);	
	[LBL, LBL_GEN, LBL_LOD, RST, CQ, eXp] = vcs_sc(S_GL, S_LL, param, ADJ, debug)
	[LBL, T] = vca_rest(adj, bus, LBL)	

	[mpw_out, dv_mean, qc_perf, svc_obj] = test_vca_plt(mpw0, mpw, SET_PQG, SET_PLT, param, debug);     
	[mpw] = set_gencosts_svc(mpw, idx_pqg);
	
  end
end

