/*********************************************************************/

/* AutoDock CrankPep, a peptide docking engine                       */
/* Copyright (C) 2019 MICHEL SANNER                                  */
/*                                                                   */ 
/* This library is free software; you can redistribute it and/or     */
/* modify it under the terms of the GNU Library General Public       */
/* License as published by the Free Software Foundation; either      */
/* version 2 of the License, or (at your option) any later version.  */
/*                                                                   */ 
/* This library is distributed in the hope that it will be useful,   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Library General Public License for more details.                  */
/*                                                                   */
/* You should have received a copy of the GNU Library General Public */
/* License along with this library; if not, see                      */
/* <https://www.gnu.org/licenses/>.                                  */
/*********************************************************************/

#ifndef CANONICALAA_H_
#define CANONICALAA_H_

#ifdef __cplusplus
extern "C" {
#endif

struct _ARG {
  int nbRot;
  int nbAtoms;
  int atypes[11];
  double charges[11];
  double coords[81][11][3];
};

struct _ASN {
  int nbRot;
  int nbAtoms;
  int atypes[5];
  double charges[5];
  double coords[18][5][3];
};

struct _ASP {
  int nbRot;
  int nbAtoms;
  int atypes[3];
  double charges[3];
  double coords[9][3][3];
};

struct _CYS {
  int nbRot;
  int nbAtoms;
  int atypes[2];
  double charges[2];
  double coords[3][2][3];
};

struct _GLN {
  int nbRot;
  int nbAtoms;
  int atypes[6];
  double charges[6];
  double coords[36][6][3];
};

struct _GLU {
  int nbRot;
  int nbAtoms;
  int atypes[4];
  double charges[4];
  double coords[27][4][3];
};

struct _HIS {
  int nbRot;
  int nbAtoms;
  int atypes[5];
  double charges[5];
  double coords[9][5][3];
};

struct _ILE {
  int nbRot;
  int nbAtoms;
  int atypes[3];
  double charges[3];
  double coords[9][3][3];
};

struct _LEU {
  int nbRot;
  int nbAtoms;
  int atypes[3];
  double charges[3];
  double coords[9][3][3];
};

struct _LYS {
  int nbRot;
  int nbAtoms;
  int atypes[7];
  double charges[7];
  double coords[81][7][3];
};

struct _MET {
  int nbRot;
  int nbAtoms;
  int atypes[3];
  double charges[3];
  double coords[27][3][3];
};

struct _PHE {
  int nbRot;
  int nbAtoms;
  int atypes[6];
  double charges[6];
  double coords[6][6][3];
};

struct _PRO {
  int nbRot;
  int nbAtoms;
  int atypes[2];
  double charges[2];
  double coords[2][2][3];
};

struct _SER {
  int nbRot;
  int nbAtoms;
  int atypes[2];
  double charges[2];
  double coords[3][2][3];
};

struct _THR {
  int nbRot;
  int nbAtoms;
  int atypes[3];
  double charges[3];
  double coords[3][3][3];
};

struct _TRP {
  int nbRot;
  int nbAtoms;
  int atypes[10];
  double charges[10];
  double coords[9][10][3];
};

struct _TYR {
  int nbRot;
  int nbAtoms;
  int atypes[8];
  double charges[8];
  double coords[6][8][3];
};

struct _VAL {
  int nbRot;
  int nbAtoms;
  int atypes[2];
  double charges[2];
  double coords[3][2][3];
};


#ifdef __cplusplus
}
#endif

#endif /* CANONICALAA_H_ */
