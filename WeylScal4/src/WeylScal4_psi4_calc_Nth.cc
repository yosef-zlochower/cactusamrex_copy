/*  File produced by Kranc */

#define KRANC_C

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Kranc.hh"
#include "Differencing.h"

namespace WeylScal4 {

extern "C" void WeylScal4_psi4_calc_Nth_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_WeylScal4_psi4_calc_Nth_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(WeylScal4_psi4_calc_Nth_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % WeylScal4_psi4_calc_Nth_calc_every != WeylScal4_psi4_calc_Nth_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi4i_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi4i_group.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4::Psi4r_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4::Psi4r_group.");
  return;
}

static void WeylScal4_psi4_calc_Nth_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Include user-supplied include files */
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const ptrdiff_t cctkLbnd1 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[0];
  const ptrdiff_t cctkLbnd2 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[1];
  const ptrdiff_t cctkLbnd3 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[2];
  const CCTK_REAL t CCTK_ATTRIBUTE_UNUSED = cctk_time;
  const CCTK_REAL cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(0);
  const CCTK_REAL cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(1);
  const CCTK_REAL cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(2);
  const CCTK_REAL dt CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_TIME;
  const CCTK_REAL dx CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dxi CCTK_ATTRIBUTE_UNUSED = pow(dx,-1);
  const CCTK_REAL dyi CCTK_ATTRIBUTE_UNUSED = pow(dy,-1);
  const CCTK_REAL dzi CCTK_ATTRIBUTE_UNUSED = pow(dz,-1);
  const CCTK_REAL khalf CCTK_ATTRIBUTE_UNUSED = 0.5;
  const CCTK_REAL kthird CCTK_ATTRIBUTE_UNUSED = 
    0.333333333333333333333333333333;
  const CCTK_REAL ktwothird CCTK_ATTRIBUTE_UNUSED = 
    0.666666666666666666666666666667;
  const CCTK_REAL kfourthird CCTK_ATTRIBUTE_UNUSED = 
    1.33333333333333333333333333333;
  const CCTK_REAL hdxi CCTK_ATTRIBUTE_UNUSED = 0.5*dxi;
  const CCTK_REAL hdyi CCTK_ATTRIBUTE_UNUSED = 0.5*dyi;
  const CCTK_REAL hdzi CCTK_ATTRIBUTE_UNUSED = 0.5*dzi;
  /* Initialize predefined quantities */
  const CCTK_REAL p1o12dx CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL p1o12dy CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL p1o12dz CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL p1o144dxdy CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o144dxdz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o144dydz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o180dx2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dx,-2);
  const CCTK_REAL p1o180dy2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dy,-2);
  const CCTK_REAL p1o180dz2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dz,-2);
  const CCTK_REAL p1o2dx CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dx,-1);
  const CCTK_REAL p1o2dy CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dy,-1);
  const CCTK_REAL p1o2dz CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dz,-1);
  const CCTK_REAL p1o3600dxdy CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o3600dxdz CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o3600dydz CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dxdy CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o4dxdz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dydz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o5040dx2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dx,-2);
  const CCTK_REAL p1o5040dy2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dy,-2);
  const CCTK_REAL p1o5040dz2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dz,-2);
  const CCTK_REAL p1o60dx CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL p1o60dy CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL p1o60dz CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o705600dydz CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o840dx CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dx,-1);
  const CCTK_REAL p1o840dy CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dy,-1);
  const CCTK_REAL p1o840dz CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dz,-1);
  const CCTK_REAL p1odx2 CCTK_ATTRIBUTE_UNUSED = pow(dx,-2);
  const CCTK_REAL p1ody2 CCTK_ATTRIBUTE_UNUSED = pow(dy,-2);
  const CCTK_REAL p1odz2 CCTK_ATTRIBUTE_UNUSED = pow(dz,-2);
  const CCTK_REAL pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dx,-2);
  const CCTK_REAL pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dy,-2);
  const CCTK_REAL pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dz,-2);
  /* Assign local copies of arrays functions */
  
  
  /* Calculate temporaries and arrays functions */
  /* Copy local copies back to grid functions */
  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];
  #pragma omp parallel
  CCTK_LOOP3(WeylScal4_psi4_calc_Nth,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL gxxL CCTK_ATTRIBUTE_UNUSED = gxx[index];
    CCTK_REAL gxyL CCTK_ATTRIBUTE_UNUSED = gxy[index];
    CCTK_REAL gxzL CCTK_ATTRIBUTE_UNUSED = gxz[index];
    CCTK_REAL gyyL CCTK_ATTRIBUTE_UNUSED = gyy[index];
    CCTK_REAL gyzL CCTK_ATTRIBUTE_UNUSED = gyz[index];
    CCTK_REAL gzzL CCTK_ATTRIBUTE_UNUSED = gzz[index];
    CCTK_REAL kxxL CCTK_ATTRIBUTE_UNUSED = kxx[index];
    CCTK_REAL kxyL CCTK_ATTRIBUTE_UNUSED = kxy[index];
    CCTK_REAL kxzL CCTK_ATTRIBUTE_UNUSED = kxz[index];
    CCTK_REAL kyyL CCTK_ATTRIBUTE_UNUSED = kyy[index];
    CCTK_REAL kyzL CCTK_ATTRIBUTE_UNUSED = kyz[index];
    CCTK_REAL kzzL CCTK_ATTRIBUTE_UNUSED = kzz[index];
    CCTK_REAL xL CCTK_ATTRIBUTE_UNUSED = x[index];
    CCTK_REAL yL CCTK_ATTRIBUTE_UNUSED = y[index];
    CCTK_REAL zL CCTK_ATTRIBUTE_UNUSED = z[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    CCTK_REAL PDstandard1gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard2gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard3gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard22gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard33gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard23gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard1gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard2gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard3gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard33gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard12gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard13gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard23gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard1gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard2gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard3gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard22gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard12gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard13gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard23gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard1gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard2gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard3gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard11gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard33gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard13gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard1gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard2gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard3gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard11gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard12gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard13gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard23gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard1gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard2gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard3gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard11gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard22gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard12gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard2kxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard3kxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard1kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard2kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard3kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard1kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard2kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard3kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard1kyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard3kyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard1kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard2kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard3kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard1kzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandard2kzz CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDstandard1gxx = PDstandardfdOrder21(&gxx[index]);
        PDstandard2gxx = PDstandardfdOrder22(&gxx[index]);
        PDstandard3gxx = PDstandardfdOrder23(&gxx[index]);
        PDstandard22gxx = PDstandardfdOrder222(&gxx[index]);
        PDstandard33gxx = PDstandardfdOrder233(&gxx[index]);
        PDstandard23gxx = PDstandardfdOrder223(&gxx[index]);
        PDstandard1gxy = PDstandardfdOrder21(&gxy[index]);
        PDstandard2gxy = PDstandardfdOrder22(&gxy[index]);
        PDstandard3gxy = PDstandardfdOrder23(&gxy[index]);
        PDstandard33gxy = PDstandardfdOrder233(&gxy[index]);
        PDstandard12gxy = PDstandardfdOrder212(&gxy[index]);
        PDstandard13gxy = PDstandardfdOrder213(&gxy[index]);
        PDstandard23gxy = PDstandardfdOrder223(&gxy[index]);
        PDstandard1gxz = PDstandardfdOrder21(&gxz[index]);
        PDstandard2gxz = PDstandardfdOrder22(&gxz[index]);
        PDstandard3gxz = PDstandardfdOrder23(&gxz[index]);
        PDstandard22gxz = PDstandardfdOrder222(&gxz[index]);
        PDstandard12gxz = PDstandardfdOrder212(&gxz[index]);
        PDstandard13gxz = PDstandardfdOrder213(&gxz[index]);
        PDstandard23gxz = PDstandardfdOrder223(&gxz[index]);
        PDstandard1gyy = PDstandardfdOrder21(&gyy[index]);
        PDstandard2gyy = PDstandardfdOrder22(&gyy[index]);
        PDstandard3gyy = PDstandardfdOrder23(&gyy[index]);
        PDstandard11gyy = PDstandardfdOrder211(&gyy[index]);
        PDstandard33gyy = PDstandardfdOrder233(&gyy[index]);
        PDstandard13gyy = PDstandardfdOrder213(&gyy[index]);
        PDstandard1gyz = PDstandardfdOrder21(&gyz[index]);
        PDstandard2gyz = PDstandardfdOrder22(&gyz[index]);
        PDstandard3gyz = PDstandardfdOrder23(&gyz[index]);
        PDstandard11gyz = PDstandardfdOrder211(&gyz[index]);
        PDstandard12gyz = PDstandardfdOrder212(&gyz[index]);
        PDstandard13gyz = PDstandardfdOrder213(&gyz[index]);
        PDstandard23gyz = PDstandardfdOrder223(&gyz[index]);
        PDstandard1gzz = PDstandardfdOrder21(&gzz[index]);
        PDstandard2gzz = PDstandardfdOrder22(&gzz[index]);
        PDstandard3gzz = PDstandardfdOrder23(&gzz[index]);
        PDstandard11gzz = PDstandardfdOrder211(&gzz[index]);
        PDstandard22gzz = PDstandardfdOrder222(&gzz[index]);
        PDstandard12gzz = PDstandardfdOrder212(&gzz[index]);
        PDstandard2kxx = PDstandardfdOrder22(&kxx[index]);
        PDstandard3kxx = PDstandardfdOrder23(&kxx[index]);
        PDstandard1kxy = PDstandardfdOrder21(&kxy[index]);
        PDstandard2kxy = PDstandardfdOrder22(&kxy[index]);
        PDstandard3kxy = PDstandardfdOrder23(&kxy[index]);
        PDstandard1kxz = PDstandardfdOrder21(&kxz[index]);
        PDstandard2kxz = PDstandardfdOrder22(&kxz[index]);
        PDstandard3kxz = PDstandardfdOrder23(&kxz[index]);
        PDstandard1kyy = PDstandardfdOrder21(&kyy[index]);
        PDstandard3kyy = PDstandardfdOrder23(&kyy[index]);
        PDstandard1kyz = PDstandardfdOrder21(&kyz[index]);
        PDstandard2kyz = PDstandardfdOrder22(&kyz[index]);
        PDstandard3kyz = PDstandardfdOrder23(&kyz[index]);
        PDstandard1kzz = PDstandardfdOrder21(&kzz[index]);
        PDstandard2kzz = PDstandardfdOrder22(&kzz[index]);
        break;
      }
      
      case 4:
      {
        PDstandard1gxx = PDstandardfdOrder41(&gxx[index]);
        PDstandard2gxx = PDstandardfdOrder42(&gxx[index]);
        PDstandard3gxx = PDstandardfdOrder43(&gxx[index]);
        PDstandard22gxx = PDstandardfdOrder422(&gxx[index]);
        PDstandard33gxx = PDstandardfdOrder433(&gxx[index]);
        PDstandard23gxx = PDstandardfdOrder423(&gxx[index]);
        PDstandard1gxy = PDstandardfdOrder41(&gxy[index]);
        PDstandard2gxy = PDstandardfdOrder42(&gxy[index]);
        PDstandard3gxy = PDstandardfdOrder43(&gxy[index]);
        PDstandard33gxy = PDstandardfdOrder433(&gxy[index]);
        PDstandard12gxy = PDstandardfdOrder412(&gxy[index]);
        PDstandard13gxy = PDstandardfdOrder413(&gxy[index]);
        PDstandard23gxy = PDstandardfdOrder423(&gxy[index]);
        PDstandard1gxz = PDstandardfdOrder41(&gxz[index]);
        PDstandard2gxz = PDstandardfdOrder42(&gxz[index]);
        PDstandard3gxz = PDstandardfdOrder43(&gxz[index]);
        PDstandard22gxz = PDstandardfdOrder422(&gxz[index]);
        PDstandard12gxz = PDstandardfdOrder412(&gxz[index]);
        PDstandard13gxz = PDstandardfdOrder413(&gxz[index]);
        PDstandard23gxz = PDstandardfdOrder423(&gxz[index]);
        PDstandard1gyy = PDstandardfdOrder41(&gyy[index]);
        PDstandard2gyy = PDstandardfdOrder42(&gyy[index]);
        PDstandard3gyy = PDstandardfdOrder43(&gyy[index]);
        PDstandard11gyy = PDstandardfdOrder411(&gyy[index]);
        PDstandard33gyy = PDstandardfdOrder433(&gyy[index]);
        PDstandard13gyy = PDstandardfdOrder413(&gyy[index]);
        PDstandard1gyz = PDstandardfdOrder41(&gyz[index]);
        PDstandard2gyz = PDstandardfdOrder42(&gyz[index]);
        PDstandard3gyz = PDstandardfdOrder43(&gyz[index]);
        PDstandard11gyz = PDstandardfdOrder411(&gyz[index]);
        PDstandard12gyz = PDstandardfdOrder412(&gyz[index]);
        PDstandard13gyz = PDstandardfdOrder413(&gyz[index]);
        PDstandard23gyz = PDstandardfdOrder423(&gyz[index]);
        PDstandard1gzz = PDstandardfdOrder41(&gzz[index]);
        PDstandard2gzz = PDstandardfdOrder42(&gzz[index]);
        PDstandard3gzz = PDstandardfdOrder43(&gzz[index]);
        PDstandard11gzz = PDstandardfdOrder411(&gzz[index]);
        PDstandard22gzz = PDstandardfdOrder422(&gzz[index]);
        PDstandard12gzz = PDstandardfdOrder412(&gzz[index]);
        PDstandard2kxx = PDstandardfdOrder42(&kxx[index]);
        PDstandard3kxx = PDstandardfdOrder43(&kxx[index]);
        PDstandard1kxy = PDstandardfdOrder41(&kxy[index]);
        PDstandard2kxy = PDstandardfdOrder42(&kxy[index]);
        PDstandard3kxy = PDstandardfdOrder43(&kxy[index]);
        PDstandard1kxz = PDstandardfdOrder41(&kxz[index]);
        PDstandard2kxz = PDstandardfdOrder42(&kxz[index]);
        PDstandard3kxz = PDstandardfdOrder43(&kxz[index]);
        PDstandard1kyy = PDstandardfdOrder41(&kyy[index]);
        PDstandard3kyy = PDstandardfdOrder43(&kyy[index]);
        PDstandard1kyz = PDstandardfdOrder41(&kyz[index]);
        PDstandard2kyz = PDstandardfdOrder42(&kyz[index]);
        PDstandard3kyz = PDstandardfdOrder43(&kyz[index]);
        PDstandard1kzz = PDstandardfdOrder41(&kzz[index]);
        PDstandard2kzz = PDstandardfdOrder42(&kzz[index]);
        break;
      }
      
      case 6:
      {
        PDstandard1gxx = PDstandardfdOrder61(&gxx[index]);
        PDstandard2gxx = PDstandardfdOrder62(&gxx[index]);
        PDstandard3gxx = PDstandardfdOrder63(&gxx[index]);
        PDstandard22gxx = PDstandardfdOrder622(&gxx[index]);
        PDstandard33gxx = PDstandardfdOrder633(&gxx[index]);
        PDstandard23gxx = PDstandardfdOrder623(&gxx[index]);
        PDstandard1gxy = PDstandardfdOrder61(&gxy[index]);
        PDstandard2gxy = PDstandardfdOrder62(&gxy[index]);
        PDstandard3gxy = PDstandardfdOrder63(&gxy[index]);
        PDstandard33gxy = PDstandardfdOrder633(&gxy[index]);
        PDstandard12gxy = PDstandardfdOrder612(&gxy[index]);
        PDstandard13gxy = PDstandardfdOrder613(&gxy[index]);
        PDstandard23gxy = PDstandardfdOrder623(&gxy[index]);
        PDstandard1gxz = PDstandardfdOrder61(&gxz[index]);
        PDstandard2gxz = PDstandardfdOrder62(&gxz[index]);
        PDstandard3gxz = PDstandardfdOrder63(&gxz[index]);
        PDstandard22gxz = PDstandardfdOrder622(&gxz[index]);
        PDstandard12gxz = PDstandardfdOrder612(&gxz[index]);
        PDstandard13gxz = PDstandardfdOrder613(&gxz[index]);
        PDstandard23gxz = PDstandardfdOrder623(&gxz[index]);
        PDstandard1gyy = PDstandardfdOrder61(&gyy[index]);
        PDstandard2gyy = PDstandardfdOrder62(&gyy[index]);
        PDstandard3gyy = PDstandardfdOrder63(&gyy[index]);
        PDstandard11gyy = PDstandardfdOrder611(&gyy[index]);
        PDstandard33gyy = PDstandardfdOrder633(&gyy[index]);
        PDstandard13gyy = PDstandardfdOrder613(&gyy[index]);
        PDstandard1gyz = PDstandardfdOrder61(&gyz[index]);
        PDstandard2gyz = PDstandardfdOrder62(&gyz[index]);
        PDstandard3gyz = PDstandardfdOrder63(&gyz[index]);
        PDstandard11gyz = PDstandardfdOrder611(&gyz[index]);
        PDstandard12gyz = PDstandardfdOrder612(&gyz[index]);
        PDstandard13gyz = PDstandardfdOrder613(&gyz[index]);
        PDstandard23gyz = PDstandardfdOrder623(&gyz[index]);
        PDstandard1gzz = PDstandardfdOrder61(&gzz[index]);
        PDstandard2gzz = PDstandardfdOrder62(&gzz[index]);
        PDstandard3gzz = PDstandardfdOrder63(&gzz[index]);
        PDstandard11gzz = PDstandardfdOrder611(&gzz[index]);
        PDstandard22gzz = PDstandardfdOrder622(&gzz[index]);
        PDstandard12gzz = PDstandardfdOrder612(&gzz[index]);
        PDstandard2kxx = PDstandardfdOrder62(&kxx[index]);
        PDstandard3kxx = PDstandardfdOrder63(&kxx[index]);
        PDstandard1kxy = PDstandardfdOrder61(&kxy[index]);
        PDstandard2kxy = PDstandardfdOrder62(&kxy[index]);
        PDstandard3kxy = PDstandardfdOrder63(&kxy[index]);
        PDstandard1kxz = PDstandardfdOrder61(&kxz[index]);
        PDstandard2kxz = PDstandardfdOrder62(&kxz[index]);
        PDstandard3kxz = PDstandardfdOrder63(&kxz[index]);
        PDstandard1kyy = PDstandardfdOrder61(&kyy[index]);
        PDstandard3kyy = PDstandardfdOrder63(&kyy[index]);
        PDstandard1kyz = PDstandardfdOrder61(&kyz[index]);
        PDstandard2kyz = PDstandardfdOrder62(&kyz[index]);
        PDstandard3kyz = PDstandardfdOrder63(&kyz[index]);
        PDstandard1kzz = PDstandardfdOrder61(&kzz[index]);
        PDstandard2kzz = PDstandardfdOrder62(&kzz[index]);
        break;
      }
      
      case 8:
      {
        PDstandard1gxx = PDstandardfdOrder81(&gxx[index]);
        PDstandard2gxx = PDstandardfdOrder82(&gxx[index]);
        PDstandard3gxx = PDstandardfdOrder83(&gxx[index]);
        PDstandard22gxx = PDstandardfdOrder822(&gxx[index]);
        PDstandard33gxx = PDstandardfdOrder833(&gxx[index]);
        PDstandard23gxx = PDstandardfdOrder823(&gxx[index]);
        PDstandard1gxy = PDstandardfdOrder81(&gxy[index]);
        PDstandard2gxy = PDstandardfdOrder82(&gxy[index]);
        PDstandard3gxy = PDstandardfdOrder83(&gxy[index]);
        PDstandard33gxy = PDstandardfdOrder833(&gxy[index]);
        PDstandard12gxy = PDstandardfdOrder812(&gxy[index]);
        PDstandard13gxy = PDstandardfdOrder813(&gxy[index]);
        PDstandard23gxy = PDstandardfdOrder823(&gxy[index]);
        PDstandard1gxz = PDstandardfdOrder81(&gxz[index]);
        PDstandard2gxz = PDstandardfdOrder82(&gxz[index]);
        PDstandard3gxz = PDstandardfdOrder83(&gxz[index]);
        PDstandard22gxz = PDstandardfdOrder822(&gxz[index]);
        PDstandard12gxz = PDstandardfdOrder812(&gxz[index]);
        PDstandard13gxz = PDstandardfdOrder813(&gxz[index]);
        PDstandard23gxz = PDstandardfdOrder823(&gxz[index]);
        PDstandard1gyy = PDstandardfdOrder81(&gyy[index]);
        PDstandard2gyy = PDstandardfdOrder82(&gyy[index]);
        PDstandard3gyy = PDstandardfdOrder83(&gyy[index]);
        PDstandard11gyy = PDstandardfdOrder811(&gyy[index]);
        PDstandard33gyy = PDstandardfdOrder833(&gyy[index]);
        PDstandard13gyy = PDstandardfdOrder813(&gyy[index]);
        PDstandard1gyz = PDstandardfdOrder81(&gyz[index]);
        PDstandard2gyz = PDstandardfdOrder82(&gyz[index]);
        PDstandard3gyz = PDstandardfdOrder83(&gyz[index]);
        PDstandard11gyz = PDstandardfdOrder811(&gyz[index]);
        PDstandard12gyz = PDstandardfdOrder812(&gyz[index]);
        PDstandard13gyz = PDstandardfdOrder813(&gyz[index]);
        PDstandard23gyz = PDstandardfdOrder823(&gyz[index]);
        PDstandard1gzz = PDstandardfdOrder81(&gzz[index]);
        PDstandard2gzz = PDstandardfdOrder82(&gzz[index]);
        PDstandard3gzz = PDstandardfdOrder83(&gzz[index]);
        PDstandard11gzz = PDstandardfdOrder811(&gzz[index]);
        PDstandard22gzz = PDstandardfdOrder822(&gzz[index]);
        PDstandard12gzz = PDstandardfdOrder812(&gzz[index]);
        PDstandard2kxx = PDstandardfdOrder82(&kxx[index]);
        PDstandard3kxx = PDstandardfdOrder83(&kxx[index]);
        PDstandard1kxy = PDstandardfdOrder81(&kxy[index]);
        PDstandard2kxy = PDstandardfdOrder82(&kxy[index]);
        PDstandard3kxy = PDstandardfdOrder83(&kxy[index]);
        PDstandard1kxz = PDstandardfdOrder81(&kxz[index]);
        PDstandard2kxz = PDstandardfdOrder82(&kxz[index]);
        PDstandard3kxz = PDstandardfdOrder83(&kxz[index]);
        PDstandard1kyy = PDstandardfdOrder81(&kyy[index]);
        PDstandard3kyy = PDstandardfdOrder83(&kyy[index]);
        PDstandard1kyz = PDstandardfdOrder81(&kyz[index]);
        PDstandard2kyz = PDstandardfdOrder82(&kyz[index]);
        PDstandard3kyz = PDstandardfdOrder83(&kyz[index]);
        PDstandard1kzz = PDstandardfdOrder81(&kzz[index]);
        PDstandard2kzz = PDstandardfdOrder82(&kzz[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL detg CCTK_ATTRIBUTE_UNUSED = 2*gxyL*gxzL*gyzL - 
      gzzL*pow(gxyL,2) + gyyL*(gxxL*gzzL - pow(gxzL,2)) - gxxL*pow(gyzL,2);
    
    CCTK_REAL invdetg CCTK_ATTRIBUTE_UNUSED = pow(detg,-1);
    
    CCTK_REAL gInv11 CCTK_ATTRIBUTE_UNUSED = invdetg*(gyyL*gzzL - 
      pow(gyzL,2));
    
    CCTK_REAL gInv12 CCTK_ATTRIBUTE_UNUSED = (gxzL*gyzL - 
      gxyL*gzzL)*invdetg;
    
    CCTK_REAL gInv13 CCTK_ATTRIBUTE_UNUSED = (-(gxzL*gyyL) + 
      gxyL*gyzL)*invdetg;
    
    CCTK_REAL gInv21 CCTK_ATTRIBUTE_UNUSED = (gxzL*gyzL - 
      gxyL*gzzL)*invdetg;
    
    CCTK_REAL gInv22 CCTK_ATTRIBUTE_UNUSED = invdetg*(gxxL*gzzL - 
      pow(gxzL,2));
    
    CCTK_REAL gInv23 CCTK_ATTRIBUTE_UNUSED = (gxyL*gxzL - 
      gxxL*gyzL)*invdetg;
    
    CCTK_REAL gInv31 CCTK_ATTRIBUTE_UNUSED = (-(gxzL*gyyL) + 
      gxyL*gyzL)*invdetg;
    
    CCTK_REAL gInv32 CCTK_ATTRIBUTE_UNUSED = (gxyL*gxzL - 
      gxxL*gyzL)*invdetg;
    
    CCTK_REAL gInv33 CCTK_ATTRIBUTE_UNUSED = invdetg*(gxxL*gyyL - 
      pow(gxyL,2));
    
    CCTK_REAL gamma111 CCTK_ATTRIBUTE_UNUSED = 0.5*(gInv11*PDstandard1gxx 
      + 2*(gInv12*PDstandard1gxy + gInv13*PDstandard1gxz) - 
      gInv12*PDstandard2gxx - gInv13*PDstandard3gxx);
    
    CCTK_REAL gamma211 CCTK_ATTRIBUTE_UNUSED = 0.5*(gInv21*PDstandard1gxx 
      + 2*(gInv22*PDstandard1gxy + gInv23*PDstandard1gxz) - 
      gInv22*PDstandard2gxx - gInv23*PDstandard3gxx);
    
    CCTK_REAL gamma311 CCTK_ATTRIBUTE_UNUSED = 0.5*(gInv31*PDstandard1gxx 
      + 2*(gInv32*PDstandard1gxy + gInv33*PDstandard1gxz) - 
      gInv32*PDstandard2gxx - gInv33*PDstandard3gxx);
    
    CCTK_REAL gamma121 CCTK_ATTRIBUTE_UNUSED = 0.5*(gInv12*PDstandard1gyy 
      + gInv11*PDstandard2gxx + gInv13*(PDstandard1gyz + PDstandard2gxz - 
      PDstandard3gxy));
    
    CCTK_REAL gamma221 CCTK_ATTRIBUTE_UNUSED = 0.5*(gInv22*PDstandard1gyy 
      + gInv21*PDstandard2gxx + gInv23*(PDstandard1gyz + PDstandard2gxz - 
      PDstandard3gxy));
    
    CCTK_REAL gamma321 CCTK_ATTRIBUTE_UNUSED = 0.5*(gInv32*PDstandard1gyy 
      + gInv31*PDstandard2gxx + gInv33*(PDstandard1gyz + PDstandard2gxz - 
      PDstandard3gxy));
    
    CCTK_REAL gamma131 CCTK_ATTRIBUTE_UNUSED = 0.5*(gInv13*PDstandard1gzz 
      + gInv11*PDstandard3gxx + gInv12*(PDstandard1gyz - PDstandard2gxz + 
      PDstandard3gxy));
    
    CCTK_REAL gamma231 CCTK_ATTRIBUTE_UNUSED = 0.5*(gInv23*PDstandard1gzz 
      + gInv21*PDstandard3gxx + gInv22*(PDstandard1gyz - PDstandard2gxz + 
      PDstandard3gxy));
    
    CCTK_REAL gamma331 CCTK_ATTRIBUTE_UNUSED = 0.5*(gInv33*PDstandard1gzz 
      + gInv31*PDstandard3gxx + gInv32*(PDstandard1gyz - PDstandard2gxz + 
      PDstandard3gxy));
    
    CCTK_REAL gamma122 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv11*(-PDstandard1gyy + 2*PDstandard2gxy) + 
      gInv12*PDstandard2gyy + gInv13*(2*PDstandard2gyz - PDstandard3gyy));
    
    CCTK_REAL gamma222 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv21*(-PDstandard1gyy + 2*PDstandard2gxy) + 
      gInv22*PDstandard2gyy + gInv23*(2*PDstandard2gyz - PDstandard3gyy));
    
    CCTK_REAL gamma322 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv31*(-PDstandard1gyy + 2*PDstandard2gxy) + 
      gInv32*PDstandard2gyy + gInv33*(2*PDstandard2gyz - PDstandard3gyy));
    
    CCTK_REAL gamma132 CCTK_ATTRIBUTE_UNUSED = 0.5*(gInv13*PDstandard2gzz 
      + gInv11*(-PDstandard1gyz + PDstandard2gxz + PDstandard3gxy) + 
      gInv12*PDstandard3gyy);
    
    CCTK_REAL gamma232 CCTK_ATTRIBUTE_UNUSED = 0.5*(gInv23*PDstandard2gzz 
      + gInv21*(-PDstandard1gyz + PDstandard2gxz + PDstandard3gxy) + 
      gInv22*PDstandard3gyy);
    
    CCTK_REAL gamma332 CCTK_ATTRIBUTE_UNUSED = 0.5*(gInv33*PDstandard2gzz 
      + gInv31*(-PDstandard1gyz + PDstandard2gxz + PDstandard3gxy) + 
      gInv32*PDstandard3gyy);
    
    CCTK_REAL gamma133 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv11*(-PDstandard1gzz + 2*PDstandard3gxz) + 
      gInv12*(-PDstandard2gzz + 2*PDstandard3gyz) + gInv13*PDstandard3gzz);
    
    CCTK_REAL gamma233 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv21*(-PDstandard1gzz + 2*PDstandard3gxz) + 
      gInv22*(-PDstandard2gzz + 2*PDstandard3gyz) + gInv23*PDstandard3gzz);
    
    CCTK_REAL gamma333 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(gInv31*(-PDstandard1gzz + 2*PDstandard3gxz) + 
      gInv32*(-PDstandard2gzz + 2*PDstandard3gyz) + gInv33*PDstandard3gzz);
    
    CCTK_REAL xmoved CCTK_ATTRIBUTE_UNUSED = xL - xorig;
    
    CCTK_REAL ymoved CCTK_ATTRIBUTE_UNUSED = yL - yorig;
    
    CCTK_REAL zmoved CCTK_ATTRIBUTE_UNUSED = zL - zorig;
    
    CCTK_REAL va1 CCTK_ATTRIBUTE_UNUSED = -ymoved;
    
    CCTK_REAL va2 CCTK_ATTRIBUTE_UNUSED = offset + xmoved;
    
    CCTK_REAL va3 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL vb1 CCTK_ATTRIBUTE_UNUSED = offset + xmoved;
    
    CCTK_REAL vb2 CCTK_ATTRIBUTE_UNUSED = ymoved;
    
    CCTK_REAL vb3 CCTK_ATTRIBUTE_UNUSED = zmoved;
    
    CCTK_REAL vc1 CCTK_ATTRIBUTE_UNUSED = ((-(gInv13*va2) + 
      gInv12*va3)*vb1 + (gInv13*va1 - gInv11*va3)*vb2 + (-(gInv12*va1) + 
      gInv11*va2)*vb3)*pow(detg,0.5);
    
    CCTK_REAL vc2 CCTK_ATTRIBUTE_UNUSED = ((-(gInv23*va2) + 
      gInv22*va3)*vb1 + (gInv23*va1 - gInv21*va3)*vb2 + (-(gInv22*va1) + 
      gInv21*va2)*vb3)*pow(detg,0.5);
    
    CCTK_REAL vc3 CCTK_ATTRIBUTE_UNUSED = ((-(gInv33*va2) + 
      gInv32*va3)*vb1 + (gInv33*va1 - gInv31*va3)*vb2 + (-(gInv32*va1) + 
      gInv31*va2)*vb3)*pow(detg,0.5);
    
    CCTK_REAL wa1 CCTK_ATTRIBUTE_UNUSED = va1;
    
    CCTK_REAL wa2 CCTK_ATTRIBUTE_UNUSED = va2;
    
    CCTK_REAL wa3 CCTK_ATTRIBUTE_UNUSED = va3;
    
    CCTK_REAL omega11 CCTK_ATTRIBUTE_UNUSED = 2*(gyzL*wa2*wa3 + 
      wa1*(gxyL*wa2 + gxzL*wa3)) + gxxL*pow(wa1,2) + gyyL*pow(wa2,2) + 
      gzzL*pow(wa3,2);
    
    CCTK_REAL ea1 CCTK_ATTRIBUTE_UNUSED = wa1*pow(omega11,-0.5);
    
    CCTK_REAL ea2 CCTK_ATTRIBUTE_UNUSED = wa2*pow(omega11,-0.5);
    
    CCTK_REAL ea3 CCTK_ATTRIBUTE_UNUSED = wa3*pow(omega11,-0.5);
    
    CCTK_REAL omega12 CCTK_ATTRIBUTE_UNUSED = ea1*(gxxL*vb1 + gxyL*vb2 + 
      gxzL*vb3) + ea2*(gxyL*vb1 + gyyL*vb2 + gyzL*vb3) + ea3*(gxzL*vb1 + 
      gyzL*vb2 + gzzL*vb3);
    
    CCTK_REAL wb1 CCTK_ATTRIBUTE_UNUSED = -(ea1*omega12) + vb1;
    
    CCTK_REAL wb2 CCTK_ATTRIBUTE_UNUSED = -(ea2*omega12) + vb2;
    
    CCTK_REAL wb3 CCTK_ATTRIBUTE_UNUSED = -(ea3*omega12) + vb3;
    
    CCTK_REAL omega22 CCTK_ATTRIBUTE_UNUSED = 2*(gyzL*wb2*wb3 + 
      wb1*(gxyL*wb2 + gxzL*wb3)) + gxxL*pow(wb1,2) + gyyL*pow(wb2,2) + 
      gzzL*pow(wb3,2);
    
    CCTK_REAL eb1 CCTK_ATTRIBUTE_UNUSED = wb1*pow(omega22,-0.5);
    
    CCTK_REAL eb2 CCTK_ATTRIBUTE_UNUSED = wb2*pow(omega22,-0.5);
    
    CCTK_REAL eb3 CCTK_ATTRIBUTE_UNUSED = wb3*pow(omega22,-0.5);
    
    CCTK_REAL omega13 CCTK_ATTRIBUTE_UNUSED = ea1*(gxxL*vc1 + gxyL*vc2 + 
      gxzL*vc3) + ea2*(gxyL*vc1 + gyyL*vc2 + gyzL*vc3) + ea3*(gxzL*vc1 + 
      gyzL*vc2 + gzzL*vc3);
    
    CCTK_REAL omega23 CCTK_ATTRIBUTE_UNUSED = eb1*(gxxL*vc1 + gxyL*vc2 + 
      gxzL*vc3) + eb2*(gxyL*vc1 + gyyL*vc2 + gyzL*vc3) + eb3*(gxzL*vc1 + 
      gyzL*vc2 + gzzL*vc3);
    
    CCTK_REAL wc1 CCTK_ATTRIBUTE_UNUSED = -(ea1*omega13) - eb1*omega23 + 
      vc1;
    
    CCTK_REAL wc2 CCTK_ATTRIBUTE_UNUSED = -(ea2*omega13) - eb2*omega23 + 
      vc2;
    
    CCTK_REAL wc3 CCTK_ATTRIBUTE_UNUSED = -(ea3*omega13) - eb3*omega23 + 
      vc3;
    
    CCTK_REAL omega33 CCTK_ATTRIBUTE_UNUSED = 2*(gyzL*wc2*wc3 + 
      wc1*(gxyL*wc2 + gxzL*wc3)) + gxxL*pow(wc1,2) + gyyL*pow(wc2,2) + 
      gzzL*pow(wc3,2);
    
    CCTK_REAL ec1 CCTK_ATTRIBUTE_UNUSED = wc1*pow(omega33,-0.5);
    
    CCTK_REAL ec2 CCTK_ATTRIBUTE_UNUSED = wc2*pow(omega33,-0.5);
    
    CCTK_REAL ec3 CCTK_ATTRIBUTE_UNUSED = wc3*pow(omega33,-0.5);
    
    CCTK_REAL isqrt2 CCTK_ATTRIBUTE_UNUSED = 0.707106781186547524;
    
    CCTK_REAL n1 CCTK_ATTRIBUTE_UNUSED = -(eb1*isqrt2);
    
    CCTK_REAL n2 CCTK_ATTRIBUTE_UNUSED = -(eb2*isqrt2);
    
    CCTK_REAL n3 CCTK_ATTRIBUTE_UNUSED = -(eb3*isqrt2);
    
    CCTK_REAL rm1 CCTK_ATTRIBUTE_UNUSED = ec1*isqrt2;
    
    CCTK_REAL rm2 CCTK_ATTRIBUTE_UNUSED = ec2*isqrt2;
    
    CCTK_REAL rm3 CCTK_ATTRIBUTE_UNUSED = ec3*isqrt2;
    
    CCTK_REAL im1 CCTK_ATTRIBUTE_UNUSED = ea1*isqrt2;
    
    CCTK_REAL im2 CCTK_ATTRIBUTE_UNUSED = ea2*isqrt2;
    
    CCTK_REAL im3 CCTK_ATTRIBUTE_UNUSED = ea3*isqrt2;
    
    CCTK_REAL rmbar1 CCTK_ATTRIBUTE_UNUSED = ec1*isqrt2;
    
    CCTK_REAL rmbar2 CCTK_ATTRIBUTE_UNUSED = ec2*isqrt2;
    
    CCTK_REAL rmbar3 CCTK_ATTRIBUTE_UNUSED = ec3*isqrt2;
    
    CCTK_REAL imbar1 CCTK_ATTRIBUTE_UNUSED = -(ea1*isqrt2);
    
    CCTK_REAL imbar2 CCTK_ATTRIBUTE_UNUSED = -(ea2*isqrt2);
    
    CCTK_REAL imbar3 CCTK_ATTRIBUTE_UNUSED = -(ea3*isqrt2);
    
    CCTK_REAL nn CCTK_ATTRIBUTE_UNUSED = isqrt2;
    
    CCTK_REAL R1212 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(-2*(gamma122*(gxxL*gamma111 + gxyL*gamma211 + gxzL*gamma311) + 
      gamma222*(gxyL*gamma111 + gyyL*gamma211 + gyzL*gamma311) + 
      (gxzL*gamma111 + gyzL*gamma211 + gzzL*gamma311)*gamma322) - 
      PDstandard11gyy + 2*(gamma121*(gxxL*gamma121 + gxyL*gamma221 + 
      gxzL*gamma321) + gamma221*(gxyL*gamma121 + gyyL*gamma221 + 
      gyzL*gamma321) + gamma321*(gxzL*gamma121 + gyzL*gamma221 + 
      gzzL*gamma321) + PDstandard12gxy) - PDstandard22gxx);
    
    CCTK_REAL R1213 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(2*(gamma121*(gxxL*gamma131 + gxyL*gamma231 + gxzL*gamma331) + 
      gamma221*(gxyL*gamma131 + gyyL*gamma231 + gyzL*gamma331) + 
      gamma321*(gxzL*gamma131 + gyzL*gamma231 + gzzL*gamma331)) - 
      2*(gamma132*(gxxL*gamma111 + gxyL*gamma211 + gxzL*gamma311) + 
      gamma232*(gxyL*gamma111 + gyyL*gamma211 + gyzL*gamma311) + 
      (gxzL*gamma111 + gyzL*gamma211 + gzzL*gamma311)*gamma332) - 
      PDstandard11gyz + PDstandard12gxz + PDstandard13gxy - PDstandard23gxx);
    
    CCTK_REAL R1223 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(2*(gamma122*(gxxL*gamma131 + gxyL*gamma231 + gxzL*gamma331) + 
      gamma222*(gxyL*gamma131 + gyyL*gamma231 + gyzL*gamma331) + 
      gamma322*(gxzL*gamma131 + gyzL*gamma231 + gzzL*gamma331)) - 
      2*(gamma132*(gxxL*gamma121 + gxyL*gamma221 + gxzL*gamma321) + 
      gamma232*(gxyL*gamma121 + gyyL*gamma221 + gyzL*gamma321) + 
      (gxzL*gamma121 + gyzL*gamma221 + gzzL*gamma321)*gamma332) - 
      PDstandard12gyz + PDstandard13gyy + PDstandard22gxz - PDstandard23gxy);
    
    CCTK_REAL R1313 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(-2*(gamma133*(gxxL*gamma111 + gxyL*gamma211 + gxzL*gamma311) + 
      gamma233*(gxyL*gamma111 + gyyL*gamma211 + gyzL*gamma311) + 
      (gxzL*gamma111 + gyzL*gamma211 + gzzL*gamma311)*gamma333) - 
      PDstandard11gzz + 2*(gamma131*(gxxL*gamma131 + gxyL*gamma231 + 
      gxzL*gamma331) + gamma231*(gxyL*gamma131 + gyyL*gamma231 + 
      gyzL*gamma331) + gamma331*(gxzL*gamma131 + gyzL*gamma231 + 
      gzzL*gamma331) + PDstandard13gxz) - PDstandard33gxx);
    
    CCTK_REAL R1323 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(2*(gamma132*(gxxL*gamma131 + gxyL*gamma231 + gxzL*gamma331) + 
      gamma232*(gxyL*gamma131 + gyyL*gamma231 + gyzL*gamma331) + 
      (gxzL*gamma131 + gyzL*gamma231 + gzzL*gamma331)*gamma332) - 
      2*(gamma133*(gxxL*gamma121 + gxyL*gamma221 + gxzL*gamma321) + 
      gamma233*(gxyL*gamma121 + gyyL*gamma221 + gyzL*gamma321) + 
      (gxzL*gamma121 + gyzL*gamma221 + gzzL*gamma321)*gamma333) - 
      PDstandard12gzz + PDstandard13gyz + PDstandard23gxz - PDstandard33gxy);
    
    CCTK_REAL R2323 CCTK_ATTRIBUTE_UNUSED = 
      0.5*(-2*(gamma133*(gxxL*gamma122 + gxyL*gamma222 + gxzL*gamma322) + 
      gamma233*(gxyL*gamma122 + gyyL*gamma222 + gyzL*gamma322) + 
      (gxzL*gamma122 + gyzL*gamma222 + gzzL*gamma322)*gamma333) - 
      PDstandard22gzz + 2*(gamma132*(gxxL*gamma132 + gxyL*gamma232 + 
      gxzL*gamma332) + gamma232*(gxyL*gamma132 + gyyL*gamma232 + 
      gyzL*gamma332) + gamma332*(gxzL*gamma132 + gyzL*gamma232 + 
      gzzL*gamma332) + PDstandard23gyz) - PDstandard33gyy);
    
    CCTK_REAL R4p1212 CCTK_ATTRIBUTE_UNUSED = kxxL*kyyL + R1212 - 
      pow(kxyL,2);
    
    CCTK_REAL R4p1213 CCTK_ATTRIBUTE_UNUSED = -(kxyL*kxzL) + kxxL*kyzL + 
      R1213;
    
    CCTK_REAL R4p1223 CCTK_ATTRIBUTE_UNUSED = -(kxzL*kyyL) + kxyL*kyzL + 
      R1223;
    
    CCTK_REAL R4p1313 CCTK_ATTRIBUTE_UNUSED = kxxL*kzzL + R1313 - 
      pow(kxzL,2);
    
    CCTK_REAL R4p1323 CCTK_ATTRIBUTE_UNUSED = -(kxzL*kyzL) + kxyL*kzzL + 
      R1323;
    
    CCTK_REAL R4p2323 CCTK_ATTRIBUTE_UNUSED = kyyL*kzzL + R2323 - 
      pow(kyzL,2);
    
    CCTK_REAL Ro111 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro112 CCTK_ATTRIBUTE_UNUSED = kxxL*gamma121 - kyyL*gamma211 
      + kxyL*(-gamma111 + gamma221) - kyzL*gamma311 + kxzL*gamma321 + 
      PDstandard1kxy - PDstandard2kxx;
    
    CCTK_REAL Ro113 CCTK_ATTRIBUTE_UNUSED = kxxL*gamma131 - kyzL*gamma211 
      + kxyL*gamma231 - kzzL*gamma311 + kxzL*(-gamma111 + gamma331) + 
      PDstandard1kxz - PDstandard3kxx;
    
    CCTK_REAL Ro121 CCTK_ATTRIBUTE_UNUSED = -(kxxL*gamma121) + 
      kyyL*gamma211 + kxyL*(gamma111 - gamma221) + kyzL*gamma311 - 
      kxzL*gamma321 - PDstandard1kxy + PDstandard2kxx;
    
    CCTK_REAL Ro122 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro123 CCTK_ATTRIBUTE_UNUSED = -(kxzL*gamma121) + 
      kxyL*gamma131 + kyyL*gamma231 - kzzL*gamma321 + kyzL*(-gamma221 + 
      gamma331) + PDstandard2kxz - PDstandard3kxy;
    
    CCTK_REAL Ro131 CCTK_ATTRIBUTE_UNUSED = -(kxxL*gamma131) + 
      kyzL*gamma211 - kxyL*gamma231 + kzzL*gamma311 + kxzL*(gamma111 - 
      gamma331) - PDstandard1kxz + PDstandard3kxx;
    
    CCTK_REAL Ro132 CCTK_ATTRIBUTE_UNUSED = kxzL*gamma121 - kxyL*gamma131 
      - kyyL*gamma231 + kzzL*gamma321 + kyzL*(gamma221 - gamma331) - 
      PDstandard2kxz + PDstandard3kxy;
    
    CCTK_REAL Ro133 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro211 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro212 CCTK_ATTRIBUTE_UNUSED = kxxL*gamma122 - kyyL*gamma221 
      + kxyL*(-gamma121 + gamma222) - kyzL*gamma321 + kxzL*gamma322 + 
      PDstandard1kyy - PDstandard2kxy;
    
    CCTK_REAL Ro213 CCTK_ATTRIBUTE_UNUSED = kxxL*gamma132 - kyzL*gamma221 
      + kxyL*gamma232 - kzzL*gamma321 + kxzL*(-gamma121 + gamma332) + 
      PDstandard1kyz - PDstandard3kxy;
    
    CCTK_REAL Ro221 CCTK_ATTRIBUTE_UNUSED = -(kxxL*gamma122) + 
      kyyL*gamma221 + kxyL*(gamma121 - gamma222) + kyzL*gamma321 - 
      kxzL*gamma322 - PDstandard1kyy + PDstandard2kxy;
    
    CCTK_REAL Ro222 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro223 CCTK_ATTRIBUTE_UNUSED = -(kxzL*gamma122) + 
      kxyL*gamma132 + kyyL*gamma232 - kzzL*gamma322 + kyzL*(-gamma222 + 
      gamma332) + PDstandard2kyz - PDstandard3kyy;
    
    CCTK_REAL Ro231 CCTK_ATTRIBUTE_UNUSED = -(kxxL*gamma132) + 
      kyzL*gamma221 - kxyL*gamma232 + kzzL*gamma321 + kxzL*(gamma121 - 
      gamma332) - PDstandard1kyz + PDstandard3kxy;
    
    CCTK_REAL Ro232 CCTK_ATTRIBUTE_UNUSED = kxzL*gamma122 - kxyL*gamma132 
      - kyyL*gamma232 + kzzL*gamma322 + kyzL*(gamma222 - gamma332) - 
      PDstandard2kyz + PDstandard3kyy;
    
    CCTK_REAL Ro233 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro311 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro312 CCTK_ATTRIBUTE_UNUSED = kxxL*gamma132 - kyyL*gamma231 
      + kxyL*(-gamma131 + gamma232) - kyzL*gamma331 + kxzL*gamma332 + 
      PDstandard1kyz - PDstandard2kxz;
    
    CCTK_REAL Ro313 CCTK_ATTRIBUTE_UNUSED = kxxL*gamma133 - kyzL*gamma231 
      + kxyL*gamma233 - kzzL*gamma331 + kxzL*(-gamma131 + gamma333) + 
      PDstandard1kzz - PDstandard3kxz;
    
    CCTK_REAL Ro321 CCTK_ATTRIBUTE_UNUSED = -(kxxL*gamma132) + 
      kyyL*gamma231 + kxyL*(gamma131 - gamma232) + kyzL*gamma331 - 
      kxzL*gamma332 - PDstandard1kyz + PDstandard2kxz;
    
    CCTK_REAL Ro322 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Ro323 CCTK_ATTRIBUTE_UNUSED = -(kxzL*gamma132) + 
      kxyL*gamma133 + kyyL*gamma233 - kzzL*gamma332 + kyzL*(-gamma232 + 
      gamma333) + PDstandard2kzz - PDstandard3kyz;
    
    CCTK_REAL Ro331 CCTK_ATTRIBUTE_UNUSED = -(kxxL*gamma133) + 
      kyzL*gamma231 - kxyL*gamma233 + kzzL*gamma331 + kxzL*(gamma131 - 
      gamma333) - PDstandard1kzz + PDstandard3kxz;
    
    CCTK_REAL Ro332 CCTK_ATTRIBUTE_UNUSED = kxzL*gamma132 - kxyL*gamma133 
      - kyyL*gamma233 + kzzL*gamma332 + kyzL*(gamma232 - gamma333) - 
      PDstandard2kzz + PDstandard3kyz;
    
    CCTK_REAL Ro333 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Rojo11 CCTK_ATTRIBUTE_UNUSED = (gInv23 + 
      gInv32)*(-(kxyL*kxzL) + kxxL*kyzL + R1213) + gInv22*(kxxL*kyyL + R1212 
      - pow(kxyL,2)) + gInv33*(kxxL*kzzL + R1313 - pow(kxzL,2));
    
    CCTK_REAL Rojo12 CCTK_ATTRIBUTE_UNUSED = (kxyL*kxzL - 
      kxxL*kyzL)*gInv13 + (-(kxzL*kyyL) + kxyL*kyzL)*gInv32 - gInv21*R1212 - 
      gInv31*R1213 + gInv23*R1223 + gInv33*(-(kxzL*kyzL) + kxyL*kzzL + R1323) 
      + gInv12*(-(kxxL*kyyL) + pow(kxyL,2));
    
    CCTK_REAL Rojo13 CCTK_ATTRIBUTE_UNUSED = (kxyL*kxzL - 
      kxxL*kyzL)*gInv12 + (kxzL*kyzL - kxyL*kzzL)*gInv23 - gInv21*R1213 + 
      gInv22*(kxzL*kyyL - kxyL*kyzL - R1223) - gInv31*R1313 - gInv32*R1323 + 
      gInv13*(-(kxxL*kzzL) + pow(kxzL,2));
    
    CCTK_REAL Rojo21 CCTK_ATTRIBUTE_UNUSED = (-(kxzL*kyyL) + 
      kxyL*kyzL)*gInv23 + (kxyL*kxzL - kxxL*kyzL)*gInv31 - gInv12*R1212 - 
      gInv13*R1213 + gInv32*R1223 + gInv33*(-(kxzL*kyzL) + kxyL*kzzL + R1323) 
      + gInv21*(-(kxxL*kyyL) + pow(kxyL,2));
    
    CCTK_REAL Rojo22 CCTK_ATTRIBUTE_UNUSED = (gInv13 + gInv31)*(kxzL*kyyL 
      - kxyL*kyzL - R1223) + gInv11*(kxxL*kyyL + R1212 - pow(kxyL,2)) + 
      gInv33*(kyyL*kzzL + R2323 - pow(kyzL,2));
    
    CCTK_REAL Rojo23 CCTK_ATTRIBUTE_UNUSED = (kxzL*kyzL - 
      kxyL*kzzL)*gInv13 + (-(kxzL*kyyL) + kxyL*kyzL)*gInv21 + 
      gInv11*(-(kxyL*kxzL) + kxxL*kyzL + R1213) + gInv12*R1223 - gInv31*R1323 
      - gInv32*R2323 + gInv23*(-(kyyL*kzzL) + pow(kyzL,2));
    
    CCTK_REAL Rojo31 CCTK_ATTRIBUTE_UNUSED = (kxyL*kxzL - 
      kxxL*kyzL)*gInv21 + (kxzL*kyzL - kxyL*kzzL)*gInv32 - gInv12*R1213 + 
      gInv22*(kxzL*kyyL - kxyL*kyzL - R1223) - gInv13*R1313 - gInv23*R1323 + 
      gInv31*(-(kxxL*kzzL) + pow(kxzL,2));
    
    CCTK_REAL Rojo32 CCTK_ATTRIBUTE_UNUSED = (-(kxzL*kyyL) + 
      kxyL*kyzL)*gInv12 + (kxzL*kyzL - kxyL*kzzL)*gInv31 + 
      gInv11*(-(kxyL*kxzL) + kxxL*kyzL + R1213) + gInv21*R1223 - gInv13*R1323 
      - gInv23*R2323 + gInv32*(-(kyyL*kzzL) + pow(kyzL,2));
    
    CCTK_REAL Rojo33 CCTK_ATTRIBUTE_UNUSED = (gInv12 + 
      gInv21)*(-(kxzL*kyzL) + kxyL*kzzL + R1323) + gInv11*(kxxL*kzzL + R1313 
      - pow(kxzL,2)) + gInv22*(kyyL*kzzL + R2323 - pow(kyzL,2));
    
    CCTK_REAL Psi4rL CCTK_ATTRIBUTE_UNUSED = 2*(n1*(n2*R4p1212 + 
      n3*R4p1213) - n3*(n2*R4p1223 + n3*R4p1323))*(imbar1*imbar2 - 
      rmbar1*rmbar2) + 2*(-(imbar2*imbar3) + rmbar2*rmbar3)*(n1*(n2*R4p1223 - 
      n3*R4p1323) - n2*n3*R4p2323 + R4p1213*pow(n1,2)) + 2*(imbar1*imbar3 - 
      rmbar1*rmbar3)*(n1*n2*R4p1213 + n1*n3*R4p1313 + n2*n3*R4p1323 + 
      R4p1223*pow(n2,2)) - (2*n2*n3*R4p1213 + R4p1212*pow(n2,2) + 
      R4p1313*pow(n3,2))*(pow(imbar1,2) - pow(rmbar1,2)) - (-2*n1*n3*R4p1223 
      + R4p1212*pow(n1,2) + R4p2323*pow(n3,2))*(pow(imbar2,2) - 
      pow(rmbar2,2)) - pow(nn,2)*((imbar2*imbar3 - rmbar2*rmbar3)*Rojo23 + 
      imbar1*(imbar2*(Rojo12 + Rojo21) + imbar3*(Rojo13 + Rojo31)) - 
      rmbar1*(rmbar2*(Rojo12 + Rojo21) + rmbar3*(Rojo13 + Rojo31)) + 
      (imbar2*imbar3 - rmbar2*rmbar3)*Rojo32 + Rojo11*(pow(imbar1,2) - 
      pow(rmbar1,2)) + Rojo22*(pow(imbar2,2) - pow(rmbar2,2)) + 
      Rojo33*(pow(imbar3,2) - pow(rmbar3,2))) - (2*n1*n2*R4p1323 + 
      R4p1313*pow(n1,2) + R4p2323*pow(n2,2))*(pow(imbar3,2) - pow(rmbar3,2)) 
      + 2*nn*((-(imbar1*imbar2) + rmbar1*rmbar2)*(n1*Ro112 + n2*Ro122 + 
      n3*Ro132) + (-(imbar1*imbar3) + rmbar1*rmbar3)*(n1*Ro113 + n2*Ro123 + 
      n3*Ro133) + (-(imbar1*imbar2) + rmbar1*rmbar2)*(n1*Ro211 + n2*Ro221 + 
      n3*Ro231) + (-(imbar2*imbar3) + rmbar2*rmbar3)*(n1*Ro213 + n2*Ro223 + 
      n3*Ro233) + (-(imbar1*imbar3) + rmbar1*rmbar3)*(n1*Ro311 + n2*Ro321 + 
      n3*Ro331) + (-(imbar2*imbar3) + rmbar2*rmbar3)*(n1*Ro312 + n2*Ro322 + 
      n3*Ro332) + (n1*Ro111 + n2*Ro121 + n3*Ro131)*(-pow(imbar1,2) + 
      pow(rmbar1,2)) + (n1*Ro212 + n2*Ro222 + n3*Ro232)*(-pow(imbar2,2) + 
      pow(rmbar2,2)) + (n1*Ro313 + n2*Ro323 + n3*Ro333)*(-pow(imbar3,2) + 
      pow(rmbar3,2)));
    
    CCTK_REAL Psi4iL CCTK_ATTRIBUTE_UNUSED = 2*((n1*(n2*R4p1212 + 
      n3*R4p1213) - n3*(n2*R4p1223 + n3*R4p1323))*(im2*rm1 + im1*rm2) + 
      nn*((im2*rm1 + im1*rm2)*(n1*(-Ro112 - Ro211) + n2*(-Ro122 - Ro221) + 
      n3*(-Ro132 - Ro231)) - 2*(im1*rm1*(n1*Ro111 + n2*Ro121 + n3*Ro131) + 
      im2*rm2*(n1*Ro212 + n2*Ro222 + n3*Ro232)) + (im3*rm1 + 
      im1*rm3)*(n1*(-Ro113 - Ro311) + n2*(-Ro123 - Ro321) + n3*(-Ro133 - 
      Ro331)) + (im3*rm2 + im2*rm3)*(n1*(-Ro213 - Ro312) + n2*(-Ro223 - 
      Ro322) + n3*(-Ro233 - Ro332)) - 2*im3*rm3*(n1*Ro313 + n2*Ro323 + 
      n3*Ro333)) + (im3*rm1 + im1*rm3)*(n1*(n2*R4p1213 + n3*R4p1313) + 
      n2*n3*R4p1323 + R4p1223*pow(n2,2))) - 2*((im3*rm2 + 
      im2*rm3)*(n1*(n2*R4p1223 - n3*R4p1323) - n2*n3*R4p2323 + 
      R4p1213*pow(n1,2)) + im3*rm3*(2*n1*n2*R4p1323 + R4p1313*pow(n1,2) + 
      R4p2323*pow(n2,2)) + im1*rm1*(2*n2*n3*R4p1213 + R4p1212*pow(n2,2) + 
      R4p1313*pow(n3,2)) + im2*rm2*(-2*n1*n3*R4p1223 + R4p1212*pow(n1,2) + 
      R4p2323*pow(n3,2))) - (im1*(2*rm1*Rojo11 + rm2*(Rojo12 + Rojo21) + 
      rm3*(Rojo13 + Rojo31)) + im2*(rm1*(Rojo12 + Rojo21) + 2*rm2*Rojo22 + 
      rm3*(Rojo23 + Rojo32)) + im3*(rm1*(Rojo13 + Rojo31) + rm2*(Rojo23 + 
      Rojo32) + 2*rm3*Rojo33))*pow(nn,2);
    /* Copy local copies back to grid functions */
    Psi4i[index] = Psi4iL;
    Psi4r[index] = Psi4rL;
  }
  CCTK_ENDLOOP3(WeylScal4_psi4_calc_Nth);
}
extern "C" void WeylScal4_psi4_calc_Nth(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_WeylScal4_psi4_calc_Nth
  DECLARE_CCTK_ARGUMENTS_CHECKED(WeylScal4_psi4_calc_Nth);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering WeylScal4_psi4_calc_Nth_Body");
  }
  if (cctk_iteration % WeylScal4_psi4_calc_Nth_calc_every != WeylScal4_psi4_calc_Nth_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "admbase::curv",
    "admbase::metric",
    "grid::coordinates",
    "WeylScal4::Psi4i_group",
    "WeylScal4::Psi4r_group"};
  AssertGroupStorage(cctkGH, "WeylScal4_psi4_calc_Nth", 5, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "WeylScal4_psi4_calc_Nth", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "WeylScal4_psi4_calc_Nth", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "WeylScal4_psi4_calc_Nth", 3, 3, 3);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "WeylScal4_psi4_calc_Nth", 4, 4, 4);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, WeylScal4_psi4_calc_Nth_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving WeylScal4_psi4_calc_Nth_Body");
  }
}

} // namespace WeylScal4
