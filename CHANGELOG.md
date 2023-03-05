# Changelog

## 2000-2015

Changes to the Mistra codebase between 2000 (version _u6_) and 2009 (version _7.4.1_) are documented
in __Appendix I__ of the manual (`doc/manual_Mistra_v9.0.pdf`).

## 2015-2017

Main updates and bugfixes to Mistra (J. Bock, 05/03/2017) compared to the version shared
with E. Droste in late December 2016/early January 2017:

- `activity.f`:
  - check that `xo4 = xo**4 > tiny(0.)`
  - plus cosmetic improvements, comments, etc...

- `nuc.f`:
  - missing initialisation of `nn` and `nh` in SR `ternucl`

- `kpp.f`:
  - bugfix SR `initc`: switch `xiod` was missing (used to initialise iodine concentration in particles)
  - updated SR `liq_parm`: use `update_now` to get all aqueous rates of `tot` if a bin is activated
  - updated SRs `st_coeff_*`: implicit none and common factor calculated only once
  - removed SRs `vmean_*`: were already no longer used in the first distributed version (automated version developed end 2016)
  - updated SRs `henry_*`: implicit none, temperature factor calculated only once
  - updated SR `cw_rc`: little inconsistency in `.gt.` and `.lt.` tests, plus mostly cosmetic work
  - updated SR `fast_k_mt_*`: print error message about the species list only once, at the beginning on one run
  - updated SR `equil_co_*`: commented rates for HCHO (unused)
  - updated SR `kpp_driver`: commented initial check `sl1` & `sion1 >0`: at the moment, still in each mech
  - updated SR `ionbalance`: implicit none
  - updated SR `dry_cw_rc`: implicit none, cleaned
  - updated SR `gasdrydep`: potential division by zero when `f0=0` avoided
  - updated functions: `fhet_t` (implicit none), `flsc6`, `uplim`, `uplip`: handle unexpected [H+] (0 or very small)
  - removed unused function `fdhet_a`

- `str.f`:
  - cleaned multiple `print*` used for debug
  - added call `v_mean_init` in restart case (allocate arrays)
  - bugfix SR `initm`: variable `rho(k)` needs `xm1(k)` but was defined before `xm1(k)` itself
  - changed unit for file `profma.out`: was unit 6, reserved for stdout, now is unit 26
  - updated function `vterm`
  - updated SR `claf`, with comments
  - updated SR `subkon`
  - added SR `diff_wat_vap` and SR `therm_conduct_air` (both new -- formerly handled in internal functions)
  - rewritten SR `advec` (change case `u<0`, and should be more efficient in case `u=0`). Note also CB59 -> arguments
  - rewritten SR `oneD_dist` (but still needs improvements; at least, no particles are "lost" now)

## 2017-2023

Changes to the Mistra codebase between October 2017 and July 2022 are documented with the corresponding
Pull Requests in the __github repository__: https://github.com/MistraModel/Mistra. See also the notes related
to version _9.0_ in __Appendix I__ of the manual (`doc/manual_Mistra_v9.0.pdf`).

Additional minor changes, from August 2022 onwards, are documented with Pull Requests #23 and later:
- apply GPLv2 license (file `LICENSE`), and add corresponding copyright notes to all source files
- move the files of example runs to the `examples/` directory, and the Ferret plotting scripts to the `scripts/` directory
- several minor updates and fixes, mostly cosmetic, to improve usability, consistency, and to clarify some points
- update the information files: `CHANGELOG.md`, `CREDITS.md`, `README.md`
