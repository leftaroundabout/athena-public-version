//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file exp-atm_cloud.c
//  \brief Problem generator for problem of a plasma interacting with a
//         static cloud of neutral gas
//
// The shock-cloud problem consists of a planar shock impacting a single spherical cloud
// Input parameters are:
//    - problem/velo_ambient    = velocity outside of the cloud (in x-direction)
//    - problem/density_plasma
//    - problem/b0              = initial magnetic field (in z-direction)
//    - problem/energy_intl     = internal energy density (without the magnetic component)
//    - problem/radius_physical = unnormalised radius of the cloud.
//
// The cloud radius is fixed at 1.0.  The center of the coordinate system defines the
// center of the cloud, and should be in the middle of the cloud. The shock is initially
// at x1=-2.0.  A typical grid domain should span x1 in [-3.0,7.0] , y and z in 
//[-2.5,2.5] (see input file in /tst).
//========================================================================================

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

// postshock flow variables are shared with IIB function
static Real uamb, rho, b0, e0;

static Real gmma1;

// fixes BCs on L-x1 (left edge) of grid to postshock flow.
void ShockCloudInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
// Set IIB value function pointer
  EnrollUserBoundaryFunction(INNER_X1, ShockCloudInnerX1);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the shock-cloud interaction test
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real gmma  = peos->GetGamma();
  gmma1 = gmma - 1.0;

  // Read input parameters
  Real rad    = pin->GetReal("problem","radius");
  uamb = pin->GetReal("problem","velo_ambient");
  rho = pin->GetReal("problem","density_plasma");
  e0 = pin->GetReal("problem","energy_intl");
  if (MAGNETIC_FIELDS_ENABLED) b0 = pin->GetReal("problem","b0");
  
  // Initialize the grid
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    phydro->u(IDN,k,j,i) = rho;
    phydro->u(IM1,k,j,i) = uamb*rho;
    phydro->u(IM2,k,j,i) = 0.0;
    phydro->u(IM3,k,j,i) = 0.0;
    phydro->u(IEN,k,j,i) = e0;

    // cloud interior
    Real diag = sqrt(SQR(pcoord->x1v(i)) + SQR(pcoord->x2v(j)) + SQR(pcoord->x3v(k)));
    if (diag < rad) {
      phydro->u(IM1,k,j,i) = 0;
    }
  }}}

  // initialize interface B, assuming longitudinal field only B=(1,0,0)
  if (MAGNETIC_FIELDS_ENABLED) {
    Real bx = 0.0;
    Real by = 0.0;
    Real bz = b0;

    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      pfield->b.x1f(k,j,i) = bx;
    }}}
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x2f(k,j,i) = by;
    }}}
    for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x3f(k,j,i) = bz;
    }}}

    // add magnetic component to total energy

    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      phydro->u(IEN,k,j,i) += 0.5*(bx*bx + by*by + bz*bz);
    }}}
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ShockCloudInnerX1()
//  \brief Sets boundary condition on left X boundary (iib) 
// Note quantities at this boundary are held fixed at the downstream state

void ShockCloudInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=1; i<=(NGHOST); ++i) {
      prim(IDN,k,j,is-i) = rho;
      prim(IVX,k,j,is-i) = uamb;
      prim(IVY,k,j,is-i) = 0.0;
      prim(IVZ,k,j,is-i) = 0.0;
      prim(IPR,k,j,is-i) = e0*gmma1;
    }
  }}
}
