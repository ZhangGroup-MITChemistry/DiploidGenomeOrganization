/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(tanhlr/cut/domainab,PairTanhlrCutDomainab)

#else

#ifndef LMP_PAIR_TANHLR_CUT_DOMAINAB_H
#define LMP_PAIR_TANHLR_CUT_DOMAINAB_H

#include "pair.h"

namespace LAMMPS_NS {

class PairTanhlrCutDomainab : public Pair {
 public:
  PairTanhlrCutDomainab(class LAMMPS *);
  ~PairTanhlrCutDomainab();

  virtual void compute(int, int);

  virtual double single(int, int, int, int, double, double, double, double &);

  virtual void settings(int, char **);
  virtual void coeff(int, char **);

  virtual double init_one(int, int);

  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);

  virtual double memory_usage();

 protected:
  double cut_global;
  double **cut;
  double ***htanh,**sigmah,**rmh, **rmh4;
  double ***ptanh,***offset,**ideal_potential;

  int nscale, *domain, *active;

  void allocate();
};

}

#endif
#endif
