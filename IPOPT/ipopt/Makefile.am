# Copyright (C) 2003, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 632 2004-07-20 23:52:15Z andreasw $

AUTOMAKE_OPTIONS = foreign

@SETPPF77COMPILE@

noinst_LIBRARIES = libipcore.a

IPINCLUDEDIR=@IPINCLUDEDIR@

preproc_sources = check_stop.F \
	constr.F \
	error.F \
	estimate_ws.F \
	eval_h_dummy.F \
	eval_hesscon_v_dummy.F \
	eval_hesslag_v_dummy.F \
	eval_hessobj_v_dummy.F \
	filter.F \
	get_scale.F \
	get_step_full.F \
	initparams.F \
	initpoint.F \
	ipopt.F \
	iter_out.F \
	ma27_call.F \
	ma28_call.F \
	ma47_call.F \
	ma48_call.F \
	ma48sym_call.F \
	ma57_call.F \
	mainloop.F \
	mc19_call.F \
	mc29_call.F \
	resto_filter.F \
	resto_kkt.F \
	resto_tron.F \
	solve_barrier.F \
	timer.F \
	wgsmp_call.F \
	wssmp_call.F

preproc_objects = constr.$(OBJEXT) \
	error.$(OBJEXT) \
	estimate_ws.$(OBJEXT) \
	eval_h_dummy.$(OBJEXT) \
	eval_hesscon_v_dummy.$(OBJEXT) \
	eval_hesslag_v_dummy.$(OBJEXT) \
	eval_hessobj_v_dummy.$(OBJEXT) \
	filter.$(OBJEXT) \
	get_scale.$(OBJEXT) \
	get_step_full.$(OBJEXT) \
	initparams.$(OBJEXT) \
	initpoint.$(OBJEXT) \
	ipopt.$(OBJEXT) \
	ma27_call.$(OBJEXT) \
	ma28_call.$(OBJEXT) \
	ma47_call.$(OBJEXT) \
	ma48_call.$(OBJEXT) \
	ma48sym_call.$(OBJEXT) \
	ma57_call.$(OBJEXT) \
	mainloop.$(OBJEXT) \
	mc19_call.$(OBJEXT) \
	mc29_call.$(OBJEXT) \
	resto_filter.$(OBJEXT) \
	resto_kkt.$(OBJEXT) \
	resto_tron.$(OBJEXT) \
	solve_barrier.$(OBJEXT) \
	timer.$(OBJEXT) \
	wgsmp_call.$(OBJEXT) \
	wssmp_call.$(OBJEXT)

ipmalloc_source = ip_malloc.c

wsmpextra_source = hsl2csc.c hsl2gcsc.c

libipcore_a_SOURCES = armijo.f \
	aug_lag.f \
	c_out.f \
	calc_bar.f \
	calc_nrm.f \
	check_basis.f \
	compute_y.f \
	cutalpha.f \
	dasv2f.f \
	filter_add.f \
	filter_check.f \
	find_param.f \
	get_bb.f \
	get_c.f \
	get_d.f \
	get_exactw.f \
	get_f.f \
	get_g.f \
	get_h.f \
	get_hv.f \
	get_ifix.f \
	get_lambda.f \
	get_pz.f \
	get_pz_cg.f \
	get_rg.f \
	get_rgb.f \
	get_rv.f \
	get_sigma.f \
	get_step_red.f \
	get_wcorr.f \
	get_ypy.f \
	get_zpz.f \
	get_zwzv.f \
	init_bnd.f \
	init_mem.f \
	ipopt_check.f \
	linesearch.f \
	mod_chol.f \
	new_basis.f \
	partition.f \
	read_entry.f \
	reorder_ib.f \
	reorder_x.f \
	res_hv.f  \
	slen.f \
	update_b.f \
	update_b_lm.f \
	update_mpec_eta.f \
	update_mu.f \
	update_nu.f \
	$(preproc_sources) \
	ffinite.c \
	Ipopt_Interface.c \
	IPOPT.INC \
	TIMER.INC

if USE_MALLOC
  libipcore_a_SOURCES += $(ipmalloc_source)
endif

if USE_WSMP
  libipcore_a_SOURCES += $(wsmpextra_source)
endif

$(preproc_objects): $(IPINCLUDEDIR)/config_f.h $(IPINCLUDEDIR)/Ipopt.h

armijo.$(OBJEXT) aug_lag.$(OBJEXT) c_out.$(OBJEXT) calc_bar.$(OBJEXT) \
  calc_nrm.$(OBJEXT) check_basis.$(OBJEXT) compute_y.$(OBJEXT) \
  constr.$(OBJEXT) cutalpha.$(OBJEXT) error.$(OBJEXT) estimate_ws.$(OBJEXT) \
  filter.$(OBJEXT) filter_add.$(OBJEXT) get_bb.$(OBJEXT) get_c.$(OBJEXT) \
  get_exactw.$(OBJEXT) get_f.$(OBJEXT) get_g.$(OBJEXT) get_h.$(OBJEXT) \
  get_hv.$(OBJEXT) get_ifix.$(OBJEXT) get_lambda.$(OBJEXT) get_pz.$(OBJEXT) \
  get_pz_cg.$(OBJEXT) get_rg.$(OBJEXT) get_rgb.$(OBJEXT) get_rv.$(OBJEXT) \
  get_scale.$(OBJEXT) get_sigma.$(OBJEXT) get_step_full.$(OBJEXT) \
  get_step_red.$(OBJEXT) get_wcorr.$(OBJEXT) get_ypy.$(OBJEXT) \
  get_zpz.$(OBJEXT) get_zwzv.$(OBJEXT) init_bnd.$(OBJEXT) \
  init_mem.$(OBJEXT) initparams.$(OBJEXT) initpoint.$(OBJEXT) \
  ipopt.$(OBJEXT) ipopt_check.$(OBJEXT) iter_out.$(OBJEXT) \
  linesearch.$(OBJEXT) ma27_call.$(OBJEXT) ma28_call.$(OBJEXT) \
  ma47_call.$(OBJEXT) ma48_call.$(OBJEXT) ma48sym_call.$(OBJEXT) \
  ma57_call.$(OBJEXT) mainloop.$(OBJEXT) mc19_call.$(OBJEXT) \
  mc29_call.$(OBJEXT) new_basis.$(OBJEXT) partition.$(OBJEXT) \
  res_hv.$(OBJEXT) resto_filter.$(OBJEXT) resto_kkt.$(OBJEXT) \
  resto_tron.$(OBJEXT) solve_barrier.$(OBJEXT) update_b.$(OBJEXT) \
  update_b_lm.$(OBJEXT) update_mpec_eta.$(OBJEXT) update_mu.$(OBJEXT) \
  update_nu.$(OBJEXT) wgsmp_call.$(OBJEXT) \
  wssmp_call.$(OBJEXT): $(srcdir)/IPOPT.INC

filter.$(OBJEXT) get_bb.$(OBJEXT) get_exactw.$(OBJEXT) get_hv.$(OBJEXT) \
  get_pz.$(OBJEXT) get_pz_cg.$(OBJEXT) get_step_full.$(OBJEXT) \
  get_ypy.$(OBJEXT) get_zwzv.$(OBJEXT) init_mem.$(OBJEXT) \
  initparams.$(OBJEXT) ipopt.$(OBJEXT) resto_filter.$(OBJEXT) \
  resto_tron.$(OBJEXT) solve_barrier.$(OBJEXT): $(srcdir)/TIMER.INC
