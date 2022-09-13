# Changelog

## Date: 05/03/2017

Author: J. Bock

Main updates and bugfix done in Mistra, as compared to the first version distributed
to Elise Droste in late December 2016/early January 2017:

- `activity.f`:
  - check that xo4 = xo**4 > tiny(0.)
  - plus cosmetic improvements, comments, etc...

- `nuc.f`:
  - missing initialisation of nn and nh in SR `ternucl`

- `kpp.f`:
  - bugfix SR `initc`: xiod was missing (used to initialise iodine concentration in particles)
  - updated SR `liq_parm`: update_now to get all aqueous rates of tot if a bin is activated
  - updated SRs `st_coeff_*`: implicit none and common factor calculated only once.
  - removed SRs `vmean_*`: were already no longer used in the first distributed version (automated version developed end 2016)
  - updated SRs `henry_*`: implicit none, temperature factor calculated only once
  - updated SR `cw_rc`: little inconsistency in .gt. and .lt. tests, plus mostly cosmetic work
  - updated SR `fast_k_mt_*`: print error message about the species list only once, at the beginning on one run
  - updated SR `equil_co_*`: commented rates for HCHO (unused)
  - updated SR `kpp_driver`: commented initial check sl1 & sion1 >0: at the moment, still in each mech
  - updated SR `ionbalance`: implicit none
  - updated SR `dry_cw_rc`: implicit none, cleaned
  - updated SR `gasdrydep`: potential division by 0 when f0=0 avoided
  - updated functions: `fhet_t` (implicit none), flsc6, uplim, uplip: handle unexpected [H+] (0 or very small)
  - removed unused function `fdhet_a`

- `str.f`:
  - cleaned multiple print* used for debug
  - added call `v_mean_init` in restart case (allocate arrays)
  - bugfix SR `initm`: rho(k) needs xm1(k) but was defined before xm1(k) itself
  - changed unit for file `profma.out`: was unit 6, reserved for stdout. Now: 26
  - updated function `vterm`
  - updated SR `claf`, with comments
  - updated SR `subkon`
  - added SR `diff_wat_vap` and SR `therm_conduct_air` (both new -- formerly handled in internal functions)
  - rewritten SR `advec` (change case u<0, and should be more efficient in case u=0) note also CB59 -> arguments
  - rewritten SR `oneD_dist` (but still needs improvements; at least, no particles are "lost" now)

## Date: 2022

Apply the GPL license

A range of minor updates, fixes mostly to improve usability, make everyting consistent, and clarify things.
