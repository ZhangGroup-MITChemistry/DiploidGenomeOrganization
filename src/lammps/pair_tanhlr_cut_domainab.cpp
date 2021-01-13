/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Arben Jusufi, Axel Kohlmeyer (Temple U.)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_tanhlr_cut_domainab.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"
#include "integrate.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairTanhlrCutDomainab::PairTanhlrCutDomainab(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairTanhlrCutDomainab::~PairTanhlrCutDomainab()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(htanh);
    memory->destroy(sigmah);
    memory->destroy(rmh);
    memory->destroy(ptanh);
    memory->destroy(offset);
    memory->destroy(ideal_potential);
  }
}

/* ---------------------------------------------------------------------- */

void PairTanhlrCutDomainab::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,imol,jmol;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,rexp,utanh,factor_lj,tanr,alpha;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imol = atom->molecule[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      jmol = atom->molecule[j];

      int iglobal, jglobal;

      iglobal = atom->tag[i];
      jglobal = atom->tag[j];
      // we do not want to apply pair potential for beads in the same segment 
      
      if ((imol<47) && (jmol<47) && (rsq < cutsq[itype][jtype])) {
      //if (rsq < cutsq[itype][jtype]) {
        
        int iscale = 3;
        // **** add intra-chrom interaction to inter-homologs ****
        // if ( (imol==jmol) || (imol%2==1 && (jmol-imol)==1) || (imol%2==0 && (imol-jmol)==1 ) ) {
        // **** ****

        if (imol==jmol) {
          if (domain[iglobal-1] == domain[jglobal-1]) {
            iscale = 1;
          } else {
            iscale = 2;
          }
          // 	} else if (abs(jmol-imol)==1){
    	    	// iscale = 3;
    	  } 

    	// else {
    	//     iscale = 1;
  	  //   }

    	else {
    		if ((imol+1)/2 == (jmol+1)/2) {
    			iscale = 3;
    		}
    		else{
    			iscale = 1;
    		}
  	    }

        alpha = ptanh[iscale][itype][jtype];
	
        if ( abs(jmol-imol)>1 || ((jmol-imol) == 1 && (imol%2 == 0)) || ((imol-jmol) == 1 && (jmol%2 == 0)) ) {
          if ( active[iglobal-1] > 0 && active[jglobal-1] > 0 ) {
            alpha = ideal_potential[active[iglobal-1]-1][active[jglobal-1]-1];
          } else {
          	alpha = ptanh[3][itype][jtype];
          }
        } 

		// printf("%d %d %d %d %d %d %12.6f\n",imol, jmol, itype, jtype, active[iglobal-1],active[jglobal-1],alpha*2);
        // debug
        //printf("ever here i? %d %d %d %d %d \n", iglobal, jglobal, domain[iglobal], domain[jglobal], iscale);
        //printf("type: %d %d %12.6f %12.6f\n", itype,jtype, ptanh[1][itype][jtype], ptanh[2][itype][jtype]);

        r = sqrt(rsq);
        if (r <= rmh[itype][jtype]) { 
            rexp = (rmh[itype][jtype]-r)*sigmah[itype][jtype];
            // for us, p is 0.5 * h
            tanr = tanh(rexp);
            utanh = alpha*(1.0+ tanr);
            // the extra negative sign is taken care in f
            fpair = factor_lj/r*alpha*sigmah[itype][jtype]*(1-tanr*tanr);
        } else {
            utanh = alpha*rmh4[itype][jtype] / rsq / rsq;
            // there is an extra r here which is to normalize position vector
            fpair = factor_lj/r*alpha*rmh4[itype][jtype]* (4.0) /rsq/rsq/r;
        }


        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          //evdwl = utanh - offset[iscale][itype][jtype];
          //The energy with the correct offset is 
          evdwl = utanh - alpha*offset[iscale][itype][jtype];
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(htanh,nscale+1,n+1,n+1,"pair:htanh");
  memory->create(sigmah,n+1,n+1,"pair:sigmah");
  memory->create(rmh,n+1,n+1,"pair:rmh");
  memory->create(rmh4,n+1,n+1,"pair:rmh4");
  memory->create(ptanh,nscale+1,n+1,n+1,"pair:ptanh");
  memory->create(offset,nscale+1,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::settings(int narg, char **arg)
{
  // format: tanh/cut/multiscale cutoff nscale lowerB upperB
  //
  if (narg != 6) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  
  nscale = atoi(arg[1]);
  int nChrom = atoi(arg[2]);
    
  // pass through the domain file
  FILE *fp = fopen(arg[3],"r");
  if (fp == NULL)
    error->all(FLERR,"Cannot open the domain file");

  int ntotal;
  ntotal = atom->natoms;

  memory->create(domain,ntotal,"pair:domain");

  char line[1024];
  char *ptr;
  int idx;
  // Loop through the rest of the lines
  while(fgets(line,1024,fp)) {
    ptr = strtok(line," \t\n\r\f");

    // skip empty lines
    if (!ptr) continue;

    // skip comment lines starting with #
    if (*ptr == '#') continue;

    idx = atoi(ptr)-1;
    ptr = strtok(NULL," \t\n\r\f");

    // The second site
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted domain file");
    // printf("reading: %d %d\n", idx,atoi(ptr));
    domain[idx] = atoi(ptr); 

  }
  fclose(fp); 

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }

  // pass through the ab state file
  fp = fopen(arg[4],"r");
  if (fp == NULL)
    error->all(FLERR,"Cannot open the ideal active id file");

  memory->create(active,ntotal,"pair:active");

  // Loop through the rest of the lines
  while(fgets(line,102400,fp)) {
    ptr = strtok(line," \t\n\r\f");

    // skip empty lines
    if (!ptr) continue;

    // skip comment lines starting with #
    if (*ptr == '#') continue;

    idx = atoi(ptr)-1;
    ptr = strtok(NULL," \t\n\r\f");

    // The second site
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted ideal active id file");
    active[idx] = atoi(ptr); 
  }
  fclose(fp); 

  // 
  fp = fopen(arg[5],"r");
  if (fp == NULL)
    error->all(FLERR,"Cannot open the ideal chromosome potential file");
  
  memory->create(ideal_potential,nChrom,nChrom,"pair:ideal");

  int ichr, jchr;
  while(fgets(line,102400,fp)) {
    ptr = strtok(line," \t\n\r\f");
    if (!ptr) continue;
    if (*ptr == '#') continue;
    ichr = atoi(ptr);
    ptr = strtok(NULL," \t\n\r\f");
    jchr = atoi(ptr);
//    printf("%d\t%d\n",ichr,jchr);

    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr) error->all(FLERR,"Incorrectly formatted ideal chromosome potential file");
//    printf("%12.6f\n",force->numeric(FLERR, ptr) * 0.5);
    ideal_potential[ichr][jchr] = force->numeric(FLERR, ptr) * 0.5;
    ideal_potential[jchr][ichr] = ideal_potential[ichr][jchr];
  }
  fclose(fp);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::coeff(int narg, char **arg)
{
  if (narg < 4+nscale || narg > 5+nscale) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double htanh_one[nscale+1];
  for (int s=1; s <= nscale; s++)
    htanh_one[s] = force->numeric(FLERR,arg[s+1]);
  double rmh_one = force->numeric(FLERR,arg[nscale+2]);
  double sigmah_one = force->numeric(FLERR,arg[nscale+3]);

  double cut_one = cut_global;
  if (narg == nscale+5) cut_one = force->numeric(FLERR,arg[nscale+4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      for (int s = 1; s <= nscale; s++) {
        htanh[s][i][j] = htanh_one[s];
      }
      sigmah[i][j] = sigmah_one;
      rmh[i][j] = rmh_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTanhlrCutDomainab::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  for (int s=1; s <= nscale; s++)
    ptanh[s][i][j] = htanh[s][i][j] * 0.5;

  if (offset_flag) {
    double rexp = (rmh[i][j]-cut[i][j])*sigmah[i][j];
    for (int s=1; s <= nscale; s++) {
        //offset[s][i][j] = ptanh[s][i][j] * (1.0 + tanh(rexp));
        //offset[s][i][j] = ptanh[s][i][j]*pow(rmh[i][j],4.0)/pow(cut[i][j],4.0);
        //The above line is the old offset 
        //The below line does not have the prefactor. The prefactor will be accounted 
        //on a case by case basis within the compute function
        offset[s][i][j] = pow(rmh[i][j],4.0)/pow(cut[i][j],4.0);
    }
  } else {
    //offset[i][j] = 0.0;
    for (int s=1; s <= nscale; s++) {
        offset[s][i][j] = 0.0;
    }
  }

  for (int s=1; s <= nscale; s++) {
    htanh[s][j][i] = htanh[s][i][j];
    ptanh[s][j][i] = ptanh[s][i][j];
    offset[s][j][i] = offset[s][i][j];
  }
  sigmah[j][i] = sigmah[i][j];
  rmh[j][i] = rmh[i][j];
  rmh4[i][j] = pow(rmh[i][j],4.0);
  rmh4[j][i] = rmh4[i][j];
  cut[j][i] = cut[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        for (int s=1; s <= nscale; s++) {
            fwrite(&htanh[s][i][j],sizeof(double),1,fp);
        }
        fwrite(&rmh[i][j],sizeof(double),1,fp);
        fwrite(&sigmah[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          for (int s=1; s <= nscale; s++) {
            fread(&htanh[s][i][j],sizeof(double),1,fp);
          }
          fread(&rmh[i][j],sizeof(double),1,fp);
          fread(&sigmah[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        for (int s=1; s <= nscale; s++) {
            MPI_Bcast(&htanh[s][i][j],1,MPI_DOUBLE,0,world);
        }
        MPI_Bcast(&rmh[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigmah[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTanhlrCutDomainab::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairTanhlrCutDomainab::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r, rexp,utanh,phitanh, tanr;

  int iglobal, jglobal;
  int imol, jmol;
  imol = atom->molecule[i];
  jmol = atom->molecule[j];

  iglobal = atom->tag[i];
  jglobal = atom->tag[j];

  // ptanh index starts with 1
  int iscale = 3;
  if (imol==jmol) {
    if (domain[iglobal-1] == domain[jglobal-1]) {
        iscale = 1;
    } else {
        iscale = 2;
    }
  } 

  r=sqrt(rsq);
  if (r <= rmh[itype][jtype]) {
      rexp = (rmh[itype][jtype]-r)*sigmah[itype][jtype];
      tanr = tanh(rexp);
      utanh = ptanh[iscale][itype][jtype]*(1.0 + tanr);
      fforce = factor_lj/r*ptanh[iscale][itype][jtype]*(1.0-tanr*tanr)*sigmah[itype][jtype];
  } else {
      utanh = ptanh[iscale][itype][jtype]*rmh4[itype][jtype] / rsq / rsq;
      // there is an extra r here which is to normalize position vector
      fforce = factor_lj/r*ptanh[iscale][itype][jtype]*rmh4[itype][jtype]* (4.0) /rsq/rsq/r;
  }


  phitanh = utanh - offset[iscale][itype][jtype];
  return factor_lj*phitanh;
}

/* ---------------------------------------------------------------------- */
double PairTanhlrCutDomainab::memory_usage()
{
  const int n=atom->ntypes;

  double bytes = Pair::memory_usage();

  bytes += 7*((n+1)*(n+1) * sizeof(double) + (n+1)*sizeof(double *));
  bytes += 1*((n+1)*(n+1) * sizeof(int) + (n+1)*sizeof(int *));

  return bytes;
}

