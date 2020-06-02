MODULE VEGE

! Level set model of fire spread across terrain.

USE COMP_FUNCTIONS
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE PART
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
IMPLICIT NONE
PRIVATE
PUBLIC INITIALIZE_LEVEL_SET_FIRESPREAD_1,INITIALIZE_LEVEL_SET_FIRESPREAD_2,LEVEL_SET_FIRESPREAD, &
       BNDRY_VEG_MASS_ENERGY_TRANSFER,RAISED_VEG_MASS_ENERGY_TRANSFER
INTEGER :: IZERO
INTEGER  :: LIMITER_LS
REAL(EB) :: B_ROTH,BETA_OP_ROTH,C_ROTH,E_ROTH
REAL(EB), POINTER, DIMENSION(:,:) :: PHI_LS_P
REAL(EB), PARAMETER :: PHI_LS_MIN=-1._EB, PHI_LS_MAX=1._EB

TYPE(LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP         !WFDS6
TYPE(LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC  !WFDS6

CONTAINS

SUBROUTINE BNDRY_VEG_MASS_ENERGY_TRANSFER(T,NM)
!
! Issues:
! 1. Are SF%VEG_FUEL_FLUX_L and SF%VEG_MOIST_FLUX_L needed in linear degradation model?
USE PHYSICAL_FUNCTIONS, ONLY : DRAG,GET_MASS_FRACTION,GET_SPECIFIC_HEAT,GET_VISCOSITY,GET_CONDUCTIVITY
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES)
REAL(EB) :: DT_BC,RDT_BC
REAL(EB), INTENT(IN) ::T
INTEGER, INTENT(IN) :: NM
INTEGER  ::  IW
INTEGER  ::  I,IIG,JJG,KKG,KKG_L,KGRID,KLOC_GAS
REAL(EB) :: CP_GAS,CP_MOIST_AND_VEG,DZVEG_L,ETAVEG_H,H_CONV_L, &
            KAPPA_VEG,K_GAS,MU_GAS,QRADM_INC,QRADP_INC,RHO_GAS, &
            TMP_BOIL,TMP_CHAR_MAX,TMP_FILM,TMP_G,DTMP_L,RE_VEG_PART,U2,V2,RE_D,Y_O2,ZVEG,ZLOC_GAS_B,ZLOC_GAS_T
!REAL(EB) :: H_CONV_FDS_WALL,DTMP_FDS_WALL,QCONF_FDS_WALL,LAMBDA_AIR,TMPG_A
INTEGER  IIVEG_L,IVEG_L,J,LBURN,LBURN_NEW,NVEG_L,I_FUEL
!REAL(EB), ALLOCATABLE, DIMENSION(:) :: VEG_DIV_QRNET_EMISS,VEG_DIV_QRNET_INC,
!         VEG_QRNET_EMISS,VEG_QRNET_INC,VEG_QRM_EMISS,VEG_QRP_EMISS, VEG_QRM_INC,VEG_QRP_INC
REAL(EB) :: VEG_DIV_QRNET_EMISS(60),VEG_DIV_QRNET_INC(60),VEG_QRNET_EMISS(0:60),VEG_QRNET_INC(0:60), &
            VEG_QRM_EMISS(0:60),VEG_QRP_EMISS(0:60), VEG_QRM_INC(0:60),VEG_QRP_INC(0:60)
REAL(EB) :: H_H2O_VEG,A_H2O_VEG,E_H2O_VEG,H_PYR_VEG,A_PYR_VEG,E_PYR_VEG,RH_PYR_VEG,                  &
            H_CHAR_VEG,A_CHAR_VEG,E_CHAR_VEG,BETA_CHAR_VEG,NU_CHAR_VEG,NU_ASH_VEG,NU_O2_CHAR_VEG
REAL(EB) :: CP_ASH,CP_CHAR,CP_H2O,CP_VEG,CP_TOTAL,DTMP_VEG,Q_VEG_CHAR,TMP_VEG,TMP_VEG_NEW, &
            CHAR_ENTHALPY_FRACTION_VEG
REAL(EB) :: CHAR_FCTR,CHAR_FCTR2,MPA_MOIST,MPA_MOIST_LOSS,MPA_MOIST_LOSS_MAX,MPA_MOIST_MIN,DMPA_VEG, &
            MPA_CHAR,MPA_VEG,MPA_CHAR_MIN,MPA_VEG_MIN,MPA_VOLIT,MPA_VOLIT_LOSS_MAX,MPA_CHAR_LOSS,MPA_ASH
REAL(EB) :: DETA_VEG,ETA_H,ETAFM_VEG,ETAFP_VEG,VEG_TMP_FACE
REAL(EB) :: QCONF_L,Q_FOR_DRYING,Q_UPTO_DRYING,Q_VEG_MOIST,Q_VEG_VOLIT,QNET_VEG,Q_FOR_VOLIT,Q_VOLIT,Q_UPTO_VOLIT
REAL(EB) :: C_DRAG,CM,CN,NUSS_HILPERT_CYL_FORCEDCONV,NUSS_MORGAN_CYL_FREECONV,HCON_VEG_FORCED,HCON_VEG_FREE, &
            LENGTH_SCALE,RAYLEIGH_NUM,ZGRIDCELL,ZGRIDCELL0,VEG_DRAG_RAMP_FCTR,VEG_DRAG_MIN
!LOGICAL  :: H_VERT_CYLINDER_LAMINAR,H_CYLINDER_RE

!INTEGER  :: IC,II,IOR,JJ,KK,IW_CELL

TYPE (WALL_TYPE),    POINTER :: WC =>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF =>NULL()

!TYPE (WALL_TYPE),    POINTER :: WC1 =>NULL() !to handle qrad on slopes
!TYPE (SURFACE_TYPE), POINTER :: SF1 =>NULL() !to handle qrad on slopes

CALL POINT_TO_MESH(NM)
!IF (VEG_LEVEL_SET_COUPLED .OR. VEG_LEVEL_SET_UNCOUPLED) RETURN

TMP_BOIL           = 373._EB
TMP_CHAR_MAX       = 1300._EB
CP_ASH             = 800._EB !J/kg/K specific heat of ash
CP_H2O             = 4190._EB !J/kg/K specific heat of water
DT_BC              = T - VEG_CLOCK_BC
RDT_BC             = 1.0_EB/DT_BC
VEG_DRAG(:,:,1:10) = 0.0_EB
!VEG_DRAG(:,:,0) = -1.0_EB !default value when no veg is present (set in init.f90)

I_FUEL = 0
IF (N_REACTIONS>0) I_FUEL = REACTION(1)%FUEL_SMIX_INDEX

! Loop through vegetation wall cells and burn
!
VEG_WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
  WC  => WALL(IW)
  IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE VEG_WALL_CELL_LOOP

  SF  => SURFACE(WC%SURF_INDEX)
!
  IF (.NOT. SF%VEGETATION) CYCLE VEG_WALL_CELL_LOOP

  H_H2O_VEG = SF%VEG_H_H2O !J/kg
  H_PYR_VEG = SF%VEG_H_PYR !J/kg
  RH_PYR_VEG = 1._EB/H_PYR_VEG
  CHAR_FCTR  = 1._EB - SF%VEG_CHAR_FRACTION
  CHAR_FCTR2 = 1._EB/CHAR_FCTR

!Gas quantities 
  IIG = WC%ONE_D%IIG
  JJG = WC%ONE_D%JJG
  KKG = WC%ONE_D%KKG
  IF(SF%VEG_NO_BURN .OR. T <= DT_BC) WC%VEG_HEIGHT = SF%VEG_HEIGHT
  VEG_DRAG(IIG,JJG,0) = REAL(KKG,EB) !for terrain location in drag calc in velo.f90
!if (nm==1 .and. iig==9 .and. jjg==5) print '(A,1x,1ES15.7,1I3,4ES15.7)','t,kkg,uvel k=1,2,3,4', &
!                                                          t,kkg,u(9,5,1),u(9,5,2),u(9,5,3),u(9,5,4)
!

!-- Simple Drag implementation, assumes veg height is <= grid cell height.
!   No Reynolds number dependence
!VEG_DRAG(IIG,JJG,1) = SF%VEG_DRAG_INI*(SF%VEG_CHAR_FRACTION + CHAR_FCTR*WC%VEG_HEIGHT/SF%VEG_HEIGHT)
!VEG_DRAG(IIG,JJG,1) = VEG_DRAG(IIG,JJG,1)*SF%VEG_HEIGHT/(Z(KKG)-Z(KKG-1))

!-- Drag varies with height above the terrain according to the fraction of the grid cell occupied by veg
!   veg height can be < or >= than grid cell height, drage is Reynolds number dependent
!   Implemented in velo.f90 
!   KKG is the grid cell in the gas phase bordering the terrain (wall). For no terrain, KKG=1 along the "ground" 
!   The Z() array is the height of the gas-phase cell. Z(0) = zmin for the current mesh 

  BF_DRAG: IF (WC%VEG_HEIGHT > 0.0_EB) THEN
 
    VEG_DRAG_RAMP_FCTR = 1.0_EB
!   IF (T-T_BEGIN <= 5.0_EB) VEG_DRAG_RAMP_FCTR = 0.20_EB*(T-T_BEGIN)

    DO KGRID=0,5
      KLOC_GAS   = KKG + KGRID            !gas-phase grid index
      ZLOC_GAS_T = Z(KLOC_GAS)  -Z(KKG-1) !height above terrain of gas-phase grid cell top
      ZLOC_GAS_B = Z(KLOC_GAS-1)-Z(KKG-1) !height above terrain of gas-phase grid cell bottom

      IF (ZLOC_GAS_T <= WC%VEG_HEIGHT) THEN !grid cell filled with veg
        IF (.NOT. SF%VEG_UNIT_DRAG_COEFF) THEN
          TMP_G = TMP(IIG,JJG,KLOC_GAS)
          RHO_GAS  = RHO(IIG,JJG,KLOC_GAS)
          ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,ZZ(IIG,JJG,KLOC_GAS,1:N_TRACKED_SPECIES))
          CALL GET_VISCOSITY(ZZ_GET,MU_GAS,TMP_G)
          U2 = 0.25*(U(IIG,JJG,KLOC_GAS)+U(IIG-1,JJG,KLOC_GAS))**2
          V2 = 0.25*(V(IIG,JJG,KLOC_GAS)+V(IIG,JJG-1,KLOC_GAS))**2
          RE_VEG_PART = 4._EB*RHO_GAS*SQRT(U2 + V2 + W(IIG,JJG,KLOC_GAS)**2)/SF%VEG_SV/MU_GAS !for cylinder particle
          C_DRAG = 0.0_EB
          IF (RE_VEG_PART > 0.0_EB) C_DRAG = DRAG(RE_VEG_PART,2) !2 is for cylinder, 1 is for sphere
        ELSE
          C_DRAG = 1.0_EB
        ENDIF
        VEG_DRAG(IIG,JJG,KGRID+1)= C_DRAG*SF%VEG_DRAG_INI*VEG_DRAG_RAMP_FCTR

      ENDIF

      IF (ZLOC_GAS_T >  WC%VEG_HEIGHT .AND. ZLOC_GAS_B < WC%VEG_HEIGHT) THEN !grid cell is partially filled with veg
        IF (.NOT. SF%VEG_UNIT_DRAG_COEFF) THEN
          TMP_G = TMP(IIG,JJG,KLOC_GAS)
          RHO_GAS  = RHO(IIG,JJG,KLOC_GAS)
          ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,ZZ(IIG,JJG,KLOC_GAS,1:N_TRACKED_SPECIES))
          CALL GET_VISCOSITY(ZZ_GET,MU_GAS,TMP_G)
          U2 = 0.25*(U(IIG,JJG,KLOC_GAS)+U(IIG-1,JJG,KLOC_GAS))**2
          V2 = 0.25*(V(IIG,JJG,KLOC_GAS)+V(IIG,JJG-1,KLOC_GAS))**2
          RE_VEG_PART = 4._EB*RHO_GAS*SQRT(U2 + V2 + W(IIG,JJG,KLOC_GAS)**2)/SF%VEG_SV/MU_GAS !for cylinder particle
          C_DRAG = 0.0_EB
          IF (RE_VEG_PART > 0.0_EB) C_DRAG = DRAG(RE_VEG_PART,2) !2 is for cylinder, 1 is for sphere
        ELSE
          C_DRAG = 1.0_EB
        ENDIF
        VEG_DRAG(IIG,JJG,KGRID+1)= &
                   C_DRAG*SF%VEG_DRAG_INI*(WC%VEG_HEIGHT-ZLOC_GAS_B)*VEG_DRAG_RAMP_FCTR/(ZLOC_GAS_T-ZLOC_GAS_B)

        IF (KGRID == 0) THEN !compute minimum drag based on user input
         VEG_DRAG_MIN = C_DRAG*SF%VEG_DRAG_INI*SF%VEG_POSTFIRE_DRAG_FCTR*VEG_DRAG_RAMP_FCTR* &
                          SF%VEG_HEIGHT/(ZLOC_GAS_T-ZLOC_GAS_B)
         VEG_DRAG(IIG,JJG,1) = MAX(VEG_DRAG(IIG,JJG,1),VEG_DRAG_MIN)
        ENDIF
!if(iig==20.and.jjg==20)print '(A,1x,4ES12.3)','C_DRAG,DRAG_INI,RAMP_FACTR,VEG_DRAG', &
!   C_DRAG,SF%VEG_DRAG_INI,VEG_DRAG_RAMP_FCTR,veg_drag(iig,jjg,1)
      ENDIF

    ENDDO

! ELSE IF (WC%VEG_HEIGHT == 0.0_EB) THEN !veg is burned away, approx drag as SF%VEG_POSTFIRE_DRAG_FCTR*original
!   IF (.NOT. SF%VEG_UNIT_DRAG_COEFF) THEN
!     TMP_G = TMP(IIG,JJG,KKG)
!     RHO_GAS  = RHO(IIG,JJG,KKG)
!     ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(IIG,JJG,KKG,1:N_TRACKED_SPECIES)
!     CALL GET_VISCOSITY(ZZ_GET,MU_GAS,TMP_G)
!     U2 = 0.25*(U(IIG,JJG,KKG)+U(IIG-1,JJG,KKG))**2
!     V2 = 0.25*(V(IIG,JJG,KKG)+V(IIG,JJG-1,KKG))**2
!     RE_VEG_PART = 4._EB*RHO_GAS*SQRT(U2 + V2 + W(IIG,JJG,KKG)**2)/SF%VEG_SV/MU_GAS !for cylinder particle
!     C_DRAG = 0.0_EB
!     IF (RE_VEG_PART > 0.0_EB) C_DRAG = DRAG(RE_VEG_PART,2) !2 is for cylinder, 1 is for sphere
!   ELSE
!     C_DRAG = 1.0_EB
!   ENDIF
!   VEG_DRAG(IIG,JJG,1)= C_DRAG*SF%VEG_DRAG_INI*SF%VEG_POSTFIRE_DRAG_FCTR*VEG_DRAG_RAMP_FCTR

  ENDIF BF_DRAG

  IF(SF%VEG_NO_BURN) CYCLE VEG_WALL_CELL_LOOP

! Initialize quantities
  Q_VEG_MOIST     = 0.0_EB
  Q_VEG_VOLIT     = 0.0_EB
  Q_UPTO_VOLIT    = 0.0_EB
  Q_VOLIT         = 0.0_EB
  Q_VEG_CHAR      = 0.0_EB
  MPA_MOIST_LOSS  = 0.0_EB
  MPA_VOLIT       = 0.0_EB
  MPA_CHAR_LOSS   = 0.0_EB
  SF%VEG_DIVQNET_L          = 0.0_EB
  SF%VEG_MOIST_FLUX_L       = 0.0_EB
  SF%VEG_FUEL_FLUX_L        = 0.0_EB
  WC%ONE_D%MASSFLUX(I_FUEL) = 0.0_EB 
  WC%ONE_D%MASSFLUX_SPEC(I_FUEL) = 0.0_EB

  IF (I_WATER > 0) WC%ONE_D%MASSFLUX(I_WATER) = 0.0_EB

  IF(WC%VEG_HEIGHT == 0.0_EB) THEN
    WC%ONE_D%TMP_F = MAX(TMP(IIG,JJG,KKG),TMPA) !Tveg=Tgas if veg is completely burned
    CYCLE VEG_WALL_CELL_LOOP
  ENDIF

! Vegetation variables and minimum bounds
  NVEG_L = SF%NVEG_L
  LBURN  = 0

!  Minimum bound on dry veg.Newer, linear or Arrhenius degradation and char
  MPA_VEG_MIN   = 0.001_EB*SF%VEG_LOAD/REAL(NVEG_L,EB) !kg/m^2

  MPA_CHAR_MIN  = SF%VEG_CHAR_FRACTION*MPA_VEG_MIN !kg/m^2
  MPA_MOIST_MIN = 0.0001_EB*SF%VEG_MOISTURE*SF%VEG_LOAD/REAL(NVEG_L,EB) !ks/m^2

  IF (SF%VEG_MOISTURE == 0.0_EB) MPA_MOIST_MIN = MPA_VEG_MIN
  DZVEG_L   = SF%VEG_HEIGHT/REAL(NVEG_L,EB)
  KAPPA_VEG = SF%VEG_KAPPA
  DETA_VEG  = DZVEG_L*KAPPA_VEG

! Find the number of computational grids cells, in the grid for the veg, with burned veg 
! and the resulting height of the unburned veg. 
! Vegetation burns downward from the top. Array index, IVEG_L for WC% quantities, starts at top of veg top.
! LBURN is the number of computational cells with burned veg. 
! LBURN=0 when no burning has occurred. 
! LBURN=2, for example, means the top two veg grid cells have burned
! LBURN = NVEG_L when veg is completely burned away

  IF (SF%VEG_CHAR_OXIDATION) THEN
    DO IVEG_L = 1,NVEG_L 
      IF(WC%ONE_D%VEG_CHARMASS_L(IVEG_L) <= MPA_CHAR_MIN .AND. WC%ONE_D%VEG_FUELMASS_L(IVEG_L) <= MPA_VEG_MIN ) LBURN = IVEG_L
    ENDDO
  ELSE
    DO IVEG_L = 1,NVEG_L 
      IF(WC%ONE_D%VEG_FUELMASS_L(IVEG_L) <= MPA_VEG_MIN) LBURN = IVEG_L
    ENDDO
  ENDIF

  LBURN_NEW          = LBURN
  WC%VEG_HEIGHT      = REAL(NVEG_L-LBURN,EB)*DZVEG_L
  MPA_VOLIT_LOSS_MAX = SF%FIRELINE_MLR_MAX*DT_BC*DZVEG_L 
  MPA_MOIST_LOSS_MAX = MPA_VOLIT_LOSS_MAX

! Determine the gas-phase vertical grid cell index, SF%VEG_KGAS_L, for each cell in the vegetation grid. 
! This is needed for cases in which the vegetation height is larger than the height of the first gas-phase grid cell
! The WC% and SF% indices are related. As the WC% index goes from LBURN+1 to NVEG_L the SF% index goes 
! from 1 to NVEG_L - LBURN.
! Also, with increasing index value in WC% and SF% we pass from the top of the vegetation to the bottom

  DO IVEG_L = 1, NVEG_L - LBURN
   SF%VEG_KGAS_L(NVEG_L-LBURN-IVEG_L+1) = KKG 
   ZVEG = REAL(IVEG_L,EB)*DZVEG_L 
   ZGRIDCELL0 = 0.0_EB
   DO KGRID = 0,5
     ZGRIDCELL = ZGRIDCELL0 + Z(KKG+KGRID) - Z(KKG+KGRID-1)
     IF (ZVEG > ZGRIDCELL0 .AND. ZVEG <= ZGRIDCELL) SF%VEG_KGAS_L(NVEG_L-LBURN-IVEG_L+1) = KKG + KGRID
     ZGRIDCELL0 = ZGRIDCELL
   ENDDO
  ENDDO

! Factors for computing divergence of incident and self emission radiant fluxes
! in vegetation fuel bed. These need to be recomputed as the height of the
! vegetation surface layer decreases with burning

! Factors for computing decay of +/- incident fluxes
  SF%VEG_FINCM_RADFCT_L(:) =  0.0_EB
  SF%VEG_FINCP_RADFCT_L(:) =  0.0_EB
  ETA_H = KAPPA_VEG*WC%VEG_HEIGHT

  DO IVEG_L = 0,NVEG_L - LBURN
    ETAFM_VEG = REAL(IVEG_L,EB)*DETA_VEG
    ETAFP_VEG = ETA_H - ETAFM_VEG
    SF%VEG_FINCM_RADFCT_L(IVEG_L) = EXP(-ETAFM_VEG)
    SF%VEG_FINCP_RADFCT_L(IVEG_L) = EXP(-ETAFP_VEG)
  ENDDO

!  Integrand for computing +/- self emission fluxes
  SF%VEG_SEMISSP_RADFCT_L(:,:) = 0.0_EB
  SF%VEG_SEMISSM_RADFCT_L(:,:) = 0.0_EB
! q+
  DO IIVEG_L = 0,NVEG_L-LBURN !veg grid coordinate
    DO IVEG_L = IIVEG_L,NVEG_L-1-LBURN !integrand index
     ETAFM_VEG = REAL((IVEG_L-IIVEG_L),EB)*DETA_VEG
     ETAFP_VEG = ETAFM_VEG + DETA_VEG
     SF%VEG_SEMISSP_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-ETAFM_VEG)*(1.0_EB - EXP(-DETA_VEG))
    ENDDO
  ENDDO
! q-
  DO IIVEG_L = 0,NVEG_L-LBURN
    DO IVEG_L = 1,IIVEG_L
     ETAFM_VEG = REAL((IIVEG_L-IVEG_L),EB)*DETA_VEG
     ETAFP_VEG = ETAFM_VEG + DETA_VEG
     SF%VEG_SEMISSM_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-ETAFM_VEG)*(1.0_EB - EXP(-DETA_VEG))
    ENDDO
  ENDDO
!
! -----------------------------------------------
! compute CONVECTIVE HEAT FLUX on vegetation
! -----------------------------------------------
! Divergence of convective and radiative heat fluxes

  DO I=1,NVEG_L-LBURN
    KKG_L  = SF%VEG_KGAS_L(I)
    TMP_G  = TMP(IIG,JJG,KKG_L)
    DTMP_L = TMP_G - WC%ONE_D%VEG_TMP_L(I+LBURN)

!Convective heat correlation for laminar flow (Holman see ref above) 
    IF (SF%VEG_HCONV_CYLLAM) H_CONV_L = 1.42_EB*(ABS(DTMP_L)/DZVEG_L)**0.25

!Convective heat correlation that accounts for air flow using forced convection correlation for
!a cylinder in a cross flow, Hilpert Correlation; Incropera & Dewitt Forth Edition p. 370
    IF(SF%VEG_HCONV_CYLRE) THEN 
     RHO_GAS  = RHO(IIG,JJG,KKG_L)
     TMP_FILM = 0.5_EB*(TMP_G + WC%ONE_D%VEG_TMP_L(I+LBURN))
     ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,ZZ(IIG,JJG,KKG,1:N_TRACKED_SPECIES))
     CALL GET_VISCOSITY(ZZ_GET,MU_GAS,TMP_FILM)
     CALL GET_CONDUCTIVITY(ZZ_GET,K_GAS,TMP_FILM) !W/m/K
     U2 = 0.25*(U(IIG,JJG,KKG_L)+U(IIG-1,JJG,KKG_L))**2
     V2 = 0.25*(V(IIG,JJG,KKG_L)+V(IIG,JJG-1,KKG_L))**2
     RE_VEG_PART = 4._EB*RHO_GAS*SQRT(U2 + V2 + W(IIG,JJG,KKG_L)**2)/SF%VEG_SV/MU_GAS

     IF(RE_VEG_PART < 4._EB) THEN
       CN = 0.989_EB
       CM = 0.330_EB
     ELSE IF (RE_VEG_PART >= 4._EB .AND. RE_VEG_PART < 40._EB) THEN
       CN = 0.911_EB
       CM = 0.385_EB
     ELSE
       CN = 0.683_EB
       CM = 0.466_EB
     ENDIF
     H_CONV_L = 0.25_EB*SF%VEG_SV*K_GAS*CN*(RE_VEG_PART**CM)*PR_ONTH !W/K/m^2
    ENDIF
!
! Use largest of natural and forced convective heat transfer
   
    HCONV_CYLMAX: IF(SF%VEG_HCONV_CYLMAX) THEN 
      RHO_GAS  = RHO(IIG,JJG,KKG_L)
      TMP_FILM = 0.5_EB*(TMP_G + WC%ONE_D%VEG_TMP_L(I+LBURN))
      ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,ZZ(IIG,JJG,KKG_L,1:N_TRACKED_SPECIES))
      CALL GET_VISCOSITY(ZZ_GET,MU_GAS,TMP_FILM)
      CALL GET_SPECIFIC_HEAT(ZZ_GET,CP_GAS,TMP_FILM)
      CALL GET_CONDUCTIVITY(ZZ_GET,K_GAS,TMP_FILM)
      U2 = 0.25*(U(IIG,JJG,KKG_L)+U(IIG-1,JJG,KKG_L))**2
      V2 = 0.25*(V(IIG,JJG,KKG_L)+V(IIG,JJG-1,KKG_L))**2
      RE_VEG_PART = 4._EB*RHO_GAS*SQRT(U2 + V2 + W(IIG,JJG,KKG_L)**2)/SF%VEG_SV/MU_GAS !for cylinder SV

! - Forced convection heat transfer coefficients in a layer
!
! Hilpert Correlation (Incropera & DeWitt Fourth Edition, p. 370) for cylinder in crossflow,
! forced convection
      IF(RE_VEG_PART < 4._EB) THEN
        CN = 0.989_EB
        CM = 0.330_EB
      ELSE IF (RE_VEG_PART >= 4._EB .AND. RE_VEG_PART < 40._EB) THEN
        CN = 0.911_EB
        CM = 0.385_EB
      ELSE
        CN = 0.683_EB
        CM = 0.466_EB
      ENDIF
      NUSS_HILPERT_CYL_FORCEDCONV = CN*(RE_VEG_PART**CM)*PR_ONTH !Nusselt number
      HCON_VEG_FORCED = 0.25_EB*SF%VEG_SV*K_GAS*NUSS_HILPERT_CYL_FORCEDCONV !W/m^2 from Hilpert (cylinder)

! - Free convection heat transfer coefficients
      LENGTH_SCALE = 4._EB/SF%VEG_SV !horizontal cylinder diameter
      RAYLEIGH_NUM = 9.8_EB*ABS(DTMP_L)*LENGTH_SCALE**3*RHO_GAS**2*CP_GAS/(TMP_FILM*MU_GAS*K_GAS)

! Morgan correlation (Incropera & DeWitt, 4th Edition, p. 501-502) for horizontal cylinder, free convection
      IF (RAYLEIGH_NUM < 0.01_EB) THEN
        CN = 0.675_EB
        CM = 0.058_EB
      ELSE IF (RAYLEIGH_NUM >= 0.01_EB .AND. RAYLEIGH_NUM < 100._EB) THEN
        CN = 1.02_EB
        CM = 0.148_EB
      ELSE IF (RAYLEIGH_NUM >= 100._EB .AND. RAYLEIGH_NUM < 10**4._EB) THEN
        CN = 0.85_EB
        CM = 0.188_EB
      ELSE IF (RAYLEIGH_NUM >= 10**4._EB .AND. RAYLEIGH_NUM < 10**7._EB) THEN
        CN = 0.48_EB
        CM = 0.25_EB
      ELSE IF (RAYLEIGH_NUM >= 10**7._EB .AND. RAYLEIGH_NUM < 10**12._EB) THEN
        CN = 0.125_EB
        CM = 0.333_EB
      ENDIF

      NUSS_MORGAN_CYL_FREECONV = CN*RAYLEIGH_NUM**CM
      HCON_VEG_FREE = 0.25_EB*SF%VEG_SV*K_GAS*NUSS_MORGAN_CYL_FREECONV !W/K/m^2

      H_CONV_L = MAX(HCON_VEG_FORCED,HCON_VEG_FREE)
    ENDIF HCONV_CYLMAX

    QCONF_L  = H_CONV_L*DTMP_L
    SF%VEG_DIVQNET_L(I) = SF%VEG_PACKING*SF%VEG_SV*QCONF_L*DZVEG_L !W/m^2 see Mell et al. 2007 IJWF accessory pub

  ENDDO
!
! -----------------------------------------------
! Compute +/- radiation fluxes and their divergence due to self emission within vegetation
! -----------------------------------------------
  LAYER_RAD_FLUXES: IF (LBURN < NVEG_L) THEN
    VEG_QRP_EMISS   = 0.0_EB ; VEG_QRM_EMISS = 0.0_EB 
    VEG_QRNET_EMISS = 0.0_EB ; VEG_DIV_QRNET_EMISS = 0.0_EB
! qe+
    DO J=0,NVEG_L-LBURN !veg grid coordinate loop
      DO I=J,NVEG_L-LBURN !integrand loop 
!        VEG_QRP_EMISS(J) =  VEG_QRP_EMISS(J) + SF%VEG_SEMISSP_RADFCT_L(I,J)*WC%VEG_TMP_L(I+LBURN)**4
         IF (I==0) THEN
           KKG_L = SF%VEG_KGAS_L(1)
         ELSE
           KKG_L = SF%VEG_KGAS_L(I)
         ENDIF
         TMP_G = TMP(IIG,JJG,KKG_L)
         WC%ONE_D%VEG_TMP_L(LBURN)    = TMP_G !for top of fuel bed
         WC%ONE_D%VEG_TMP_L(NVEG_L+1) = WC%ONE_D%VEG_TMP_L(NVEG_L) !for bottom of fuel bed
         VEG_TMP_FACE = 0.5_EB*(WC%ONE_D%VEG_TMP_L(I+LBURN)+WC%ONE_D%VEG_TMP_L(I+LBURN+1))
         VEG_QRP_EMISS(J) =  VEG_QRP_EMISS(J) + SF%VEG_SEMISSP_RADFCT_L(I,J)*VEG_TMP_FACE**4

      ENDDO
    ENDDO
! qe-
    DO J=0,NVEG_L-LBURN  !veg grid coordinate
      DO I=0,J           !integrand for q-
!        VEG_QRM_EMISS(J) = VEG_QRM_EMISS(J) + SF%VEG_SEMISSM_RADFCT_L(I,J)*WC%VEG_TMP_L(I+LBURN)**4
         IF (I==0) THEN
           KKG_L = SF%VEG_KGAS_L(1)
         ELSE
           KKG_L = SF%VEG_KGAS_L(I)
         ENDIF
         TMP_G = TMP(IIG,JJG,KKG_L)
         WC%ONE_D%VEG_TMP_L(LBURN)    = TMP_G
         WC%ONE_D%VEG_TMP_L(NVEG_L+1) = WC%ONE_D%VEG_TMP_L(NVEG_L) 
         VEG_TMP_FACE = 0.5_EB*(WC%ONE_D%VEG_TMP_L(I+LBURN)+WC%ONE_D%VEG_TMP_L(I+LBURN+1))
         VEG_QRM_EMISS(J) =  VEG_QRM_EMISS(J) + SF%VEG_SEMISSM_RADFCT_L(I,J)*VEG_TMP_FACE**4

      ENDDO
    ENDDO
    VEG_QRP_EMISS =  VEG_QRP_EMISS*SIGMA
    VEG_QRM_EMISS =  VEG_QRM_EMISS*SIGMA
!
    DO I=0,NVEG_L-LBURN
      VEG_QRNET_EMISS(I) = VEG_QRP_EMISS(I)-VEG_QRM_EMISS(I)
    ENDDO
!    DO I=1,NVEG_L-LBURN
!      VEG_QRNET_EMISS(I)  = VEG_QRNET_EMISS(I) - VEG_QRM_EMISS(I)
!    ENDDO
!
    DO I=1,NVEG_L-LBURN
      VEG_DIV_QRNET_EMISS(I) = VEG_QRNET_EMISS(I-1) - VEG_QRNET_EMISS(I)
    ENDDO
!
! Compute +/- radiation fluxes and their divergence due to incident fluxes on boundaries
    QRADM_INC = WC%ONE_D%Q_RAD_IN/SF%EMISSIVITY !sigma*Ta^4 + flame
!   QRADM_INC = QRADIN(IW)/E_WALL(IW) + SIGMA*TMP_F(IW)**4 ! as done in FDS4
!   print*,'vege: QRADIN(IW)',qradin(iw)

! Adjust incident radiant flux to account for sloped terrain
! assumes user put VEG_NO_BURN=.TRUE. for vertical faces
! sets qrad on cell downspread of vertical face = qrad on cell face upspread of vertical face
!   QRADM_INC = QRADM_INC*1.0038_EB !adjustment for horizontal faces assuming 5 degree slope
!   II = WC%II
!   JJ = WC%JJ
!   KK = WC%KK
!   IC = CELL_INDEX(II-1,JJ,KK)
!   IOR = 1
!   IW_CELL = WALL_INDEX(IC,IOR) 
!   WC1 => WALL(IW_CELL)
!   SF1 => SURFACE(WC1%SURF_INDEX)
!print*,'vege: i,j,k,iw,sf',ii,jj,kk,iw,sf1%veg_no_burn
!   IF(SF1%VEG_NO_BURN) THEN
!print*,'vege: in vertical face qrad determination'
!!   QRADM_INC_SLOPE_VERTFACE = QRADM_INC_SLOPE_VERTFACE + WALL(IW_CELL)%RADIN/WALL(IW_CELL)%E_WALL
!!   QRADM_INC_SLOPE_VERTFACE = QRADM_INC_SLOPE_VERTFACE*0.0872_EB !assumes 5 degree slope
!!   QRADM_INC = QRADM_INC + QRADM_INC_SLOPE_VERTFACE !adjustment for adjacent vertical faces

!   IOR = -3
!   IW_CELL = WALL_INDEX(IC,IOR)
!adjustment for horizontal faces downspread of vertical face
!set flux = to max of flux up or downspread 
!print*,'vege: i,j,k,iw,qr',ii,jj,kk,wall(iw_cell)%qradin,wall(iw)%qradin
!   WALL(IW)%QRADIN = MAX(WALL(IW_CELL)%QRADIN,WALL(IW)%QRADIN) 
!   QRADM_INC = 1.0038_EB*WALL(IW)%QRADIN/WALL(IW_CELL)%E_WALL !assumes 5 degree slope!!!
!print*,'vege: qradm_inc,wallqrad',qradm_inc,wall(iw)%qradin
!   ENDIF

    ETAVEG_H  = (NVEG_L - LBURN)*DETA_VEG
    !this QRADP_INC ensures zero net radiant fluxes at bottom of vegetation (Albini)
    IF(SF%VEG_GROUND_ZERO_RAD) QRADP_INC = QRADM_INC*SF%VEG_FINCM_RADFCT_L(NVEG_L-LBURN) + VEG_QRM_EMISS(NVEG_L-LBURN)
    !this QRADP_INC assumes the ground stays at user specified temperature
    IF(.NOT. SF%VEG_GROUND_ZERO_RAD) QRADP_INC = SIGMA*SF%VEG_GROUND_TEMP**4
!   QRADP_INC = SIGMA*WC%VEG_TMP_L(NVEG_L)**4 
!   IF(.NOT. SF%VEG_GROUND_ZERO_RAD) QRADP_INC = SIGMA*TMP_G**4
!   QRADP_INC = SIGMA*WC%VEG_TMP_L(NVEG_L)**4*EXP(-ETAVEG_H) + VEG_QRM_EMISS(NVEG_L-LBURN) !fds4
    VEG_QRM_INC   = 0.0_EB ; VEG_QRP_INC = 0.0_EB 
    VEG_QRNET_INC = 0.0_EB ; VEG_DIV_QRNET_INC = 0.0_EB
    DO I=0,NVEG_L-LBURN
      VEG_QRM_INC(I)   = QRADM_INC*SF%VEG_FINCM_RADFCT_L(I)
      VEG_QRP_INC(I)   = QRADP_INC*SF%VEG_FINCP_RADFCT_L(I)
      VEG_QRNET_INC(I) = VEG_QRP_INC(I)-VEG_QRM_INC(I)
    ENDDO
    DO I=1,NVEG_L-LBURN
      VEG_DIV_QRNET_INC(I) = VEG_QRNET_INC(I-1) - VEG_QRNET_INC(I)
    ENDDO
  ENDIF LAYER_RAD_FLUXES
!
! Add divergence of net radiation flux to divergence of convection flux
  DO I=1,NVEG_L-LBURN
!if(nm==2 .and. iig==26 .and. jjg==18) print '(A,1x,I3,7ES13.3)','I,time,hveg,Tveg,Divqc,Divqr,mh2o,mfuel', &
!   i,t,wc%veg_height,wc%veg_tmp_l(i),sf%veg_divqnet_l(i), -(veg_div_qrnet_inc(i) + veg_div_qrnet_emiss(i)), & 
!   wc%veg_moistmass_l(i),wc%veg_fuelmass_l(i)
    SF%VEG_DIVQNET_L(I)= SF%VEG_DIVQNET_L(I) - (VEG_DIV_QRNET_INC(I) + VEG_DIV_QRNET_EMISS(I)) !includes self emiss
!   SF%VEG_DIVQNET_L(I)= SF%VEG_DIVQNET_L(I) - VEG_DIV_QRNET_INC(I)                            !no self emiss
  ENDDO
!
!
!      ************** Boundary Fuel Non-Arrehnius (Linear in temp) Degradation model *************************
! Drying occurs if qnet > 0 with Tveg held at 100 c
! Pyrolysis occurs according to Morvan & Dupuy empirical formula. Linear
! temperature dependence with qnet factor
!

  IF_VEG_DEGRADATION_LINEAR: IF (SF%VEG_DEGRADATION == 'LINEAR') THEN

    LAYER_LOOP1: DO IVEG_L = LBURN+1,NVEG_L
!
! Compute temperature of vegetation
!
      MPA_CHAR    = WC%ONE_D%VEG_CHARMASS_L(IVEG_L)
      MPA_VEG     = WC%ONE_D%VEG_FUELMASS_L(IVEG_L)
      MPA_MOIST   = WC%ONE_D%VEG_MOISTMASS_L(IVEG_L)
      TMP_VEG     = WC%ONE_D%VEG_TMP_L(IVEG_L)
      QNET_VEG    = SF%VEG_DIVQNET_L(IVEG_L-LBURN)
      CP_VEG      = (0.01_EB + 0.0037_EB*TMP_VEG)*1000._EB !J/kg/K
      CP_CHAR     = 420._EB + 2.09_EB*TMP_VEG + 6.85E-4_EB*TMP_VEG**2 !J/kg/K Park etal. C&F 2010 147:481-494
      CP_TOTAL    = CP_H2O*MPA_MOIST +  CP_VEG*MPA_VEG + CP_CHAR*MPA_CHAR
      DTMP_VEG    = DT_BC*QNET_VEG/CP_TOTAL
      TMP_VEG_NEW = TMP_VEG + DTMP_VEG

      IF_DIVQ_L_GE_0: IF(QNET_VEG > 0._EB) THEN 
!if(nm==2 .and. iig==26 .and. jjg==18) print '(A,1x,I3,2ES13.3)','I,time,Tveg_new',iveg_l,t,tmp_veg_new 

! -- drying of veg layer (Linear)
      IF(MPA_MOIST > MPA_MOIST_MIN .AND. TMP_VEG_NEW >= TMP_BOIL) THEN
        Q_UPTO_DRYING  = MAX(CP_TOTAL*(TMP_BOIL-TMP_VEG),0.0_EB)
        Q_FOR_DRYING   = DT_BC*QNET_VEG - Q_UPTO_DRYING
        MPA_MOIST_LOSS = MIN(Q_FOR_DRYING/H_H2O_VEG,MPA_MOIST_LOSS_MAX)
!       Q_FOR_DRYING   = (TMP_VEG_NEW - TMP_BOIL)/DTMP_VEG * QNET_VEG
!       MPA_MOIST_LOSS = MIN(DT_BC*Q_FOR_DRYING/H_H2O_VEG,MPA_MOIST_LOSS_MAX)
        MPA_MOIST_LOSS = MIN(MPA_MOIST_LOSS,MPA_MOIST-MPA_MOIST_MIN)
        TMP_VEG_NEW    = TMP_BOIL
        WC%ONE_D%VEG_MOISTMASS_L(IVEG_L) = MPA_MOIST - MPA_MOIST_LOSS !kg/m^2
        IF( WC%ONE_D%VEG_MOISTMASS_L(IVEG_L) <= MPA_MOIST_MIN ) WC%ONE_D%VEG_MOISTMASS_L(IVEG_L) = 0.0_EB
        IF (I_WATER > 0) WC%ONE_D%MASSFLUX(I_WATER) = WC%ONE_D%MASSFLUX(I_WATER) + RDT_BC*MPA_MOIST_LOSS
      ENDIF

! -- pyrolysis (Linear)
      IF_VOLITIZATION: IF (MPA_MOIST <= MPA_MOIST_MIN) THEN

        IF(TMP_VEG_NEW >= 400._EB .AND. MPA_VEG > MPA_VEG_MIN) THEN
          Q_UPTO_VOLIT = MAX(CP_TOTAL*(400._EB-TMP_VEG),0.0_EB)
!         Q_UPTO_VOLIT = CP_VEG*MPA_VEG*(400._EB-TMP_VEG)
          Q_FOR_VOLIT  = DT_BC*QNET_VEG - Q_UPTO_VOLIT
          Q_VOLIT      = Q_FOR_VOLIT*0.01_EB*(TMP_VEG-400._EB)

          MPA_VOLIT    = CHAR_FCTR*Q_VOLIT*RH_PYR_VEG
          MPA_VOLIT    = MAX(MPA_VOLIT,0._EB)
          MPA_VOLIT    = MIN(MPA_VOLIT,MPA_VOLIT_LOSS_MAX) !user specified max

          DMPA_VEG     = CHAR_FCTR2*MPA_VOLIT
          DMPA_VEG     = MIN(DMPA_VEG,(MPA_VEG-MPA_VEG_MIN))
          MPA_VEG      = MPA_VEG - DMPA_VEG

          MPA_VOLIT    = CHAR_FCTR*DMPA_VEG
          MPA_CHAR     = MPA_CHAR + SF%VEG_CHAR_FRACTION*DMPA_VEG
          Q_VOLIT      = MPA_VOLIT*H_PYR_VEG 

          TMP_VEG_NEW  = TMP_VEG + (Q_FOR_VOLIT-Q_VOLIT)/(MPA_VEG*CP_VEG + MPA_CHAR*CP_CHAR)
          TMP_VEG_NEW  = MIN(TMP_VEG_NEW,500._EB)
          WC%ONE_D%VEG_CHARMASS_L(IVEG_L) = MPA_CHAR
          WC%ONE_D%VEG_FUELMASS_L(IVEG_L) = MPA_VEG
          IF( WC%ONE_D%VEG_FUELMASS_L(IVEG_L) <= MPA_VEG_MIN ) WC%ONE_D%VEG_FUELMASS_L(IVEG_L) = 0.0_EB !**
          WC%ONE_D%MASSFLUX(I_FUEL)= WC%ONE_D%MASSFLUX(I_FUEL) + RDT_BC*MPA_VOLIT
        ENDIF        

      ENDIF IF_VOLITIZATION

      ENDIF IF_DIVQ_L_GE_0
      
      IF(MPA_VEG <= MPA_VEG_MIN) LBURN_NEW = MIN(LBURN_NEW+1,NVEG_L)
      IF(TMP_VEG_NEW < TMPA) TMP_VEG_NEW = TMPA !clip
      WC%ONE_D%VEG_TMP_L(IVEG_L) = TMP_VEG_NEW
      WC%VEG_HEIGHT        = REAL(NVEG_L-LBURN_NEW,EB)*DZVEG_L

    ENDDO LAYER_LOOP1

  ENDIF  IF_VEG_DEGRADATION_LINEAR

!      ************** Boundary Fuel Arrehnius Degradation model *************************
! Drying and pyrolysis occur according to Arrehnius expressions obtained 
! from the literature (Porterie et al., Num. Heat Transfer, 47:571-591, 2005
! Predicting wildland fire behavior and emissions using a fine-scale physical
! model

  IF_VEG_DEGRADATION_ARRHENIUS: IF(SF%VEG_DEGRADATION == 'ARRHENIUS') THEN
! Defaults:
!   A_H2O_VEG      = 600000._EB !1/s sqrt(K)
!   E_H2O_VEG      = 5800._EB !K

!   A_PYR_VEG      = 36300._EB !1/s
!   E_PYR_VEG      = 7250._EB !K

!   A_CHAR_VEG     = 430._EB !m/s
!   E_CHAR_VEG     = 9000._EB !K
!   H_CHAR_VEG     = -12.0E+6_EB !J/kg

!   BETA_CHAR_VEG  = 0.2_EB
!   NU_CHAR_VEG    = SF%VEG_CHAR_FRACTION
!   NU_ASH_VEG     = 0.1_EB
!   NU_O2_CHAR_VEG = 1.65_EB
!   CHAR_ENTHALPY_FRACTION_VEG = 0.5_EB

    A_H2O_VEG      = SF%VEG_A_H2O !1/2 sqrt(K)
    E_H2O_VEG      = SF%VEG_E_H2O !K

    A_PYR_VEG      = SF%VEG_A_PYR !1/s
    E_PYR_VEG      = SF%VEG_E_PYR !K

    A_CHAR_VEG     = SF%VEG_A_CHAR !m/s
    E_CHAR_VEG     = SF%VEG_E_CHAR !K
    H_CHAR_VEG     = SF%VEG_H_CHAR !J/kg

    BETA_CHAR_VEG  = SF%VEG_BETA_CHAR
    NU_CHAR_VEG    = SF%VEG_CHAR_FRACTION
    NU_ASH_VEG     = SF%VEG_ASH_FRACTION/SF%VEG_CHAR_FRACTION !fraction of char that can become ash
    NU_O2_CHAR_VEG = SF%VEG_NU_O2_CHAR
    CHAR_ENTHALPY_FRACTION_VEG = SF%VEG_CHAR_ENTHALPY_FRACTION

    LAYER_LOOP2: DO IVEG_L = LBURN+1,NVEG_L

      MPA_MOIST = WC%ONE_D%VEG_MOISTMASS_L(IVEG_L)
      MPA_VEG   = WC%ONE_D%VEG_FUELMASS_L(IVEG_L)
      MPA_CHAR  = WC%ONE_D%VEG_CHARMASS_L(IVEG_L)
      MPA_ASH   = WC%ONE_D%VEG_ASHMASS_L(IVEG_L)
      TMP_VEG   = WC%ONE_D%VEG_TMP_L(IVEG_L)

      TEMP_THRESEHOLD: IF (WC%ONE_D%VEG_TMP_L(IVEG_L) > 323._EB) THEN
              !arbitrary thresehold to prevent low-temp hrr reaction
              !added for drainage runs

! Drying of vegetation (Arrhenius)
      IF_DEHYDRATION_2: IF (MPA_MOIST > MPA_MOIST_MIN) THEN
        MPA_MOIST_LOSS = MIN(DT_BC*MPA_MOIST*A_H2O_VEG*EXP(-E_H2O_VEG/TMP_VEG)/SQRT(TMP_VEG), &
                         MPA_MOIST-MPA_MOIST_MIN)
        MPA_MOIST_LOSS = MIN(MPA_MOIST_LOSS,MPA_MOIST_LOSS_MAX) !user specified max
        MPA_MOIST      = MPA_MOIST - MPA_MOIST_LOSS
        WC%ONE_D%VEG_MOISTMASS_L(IVEG_L) = MPA_MOIST !kg/m^2
        IF (MPA_MOIST <= MPA_MOIST_MIN) WC%ONE_D%VEG_MOISTMASS_L(IVEG_L) = 0.0_EB
      ENDIF IF_DEHYDRATION_2

! Volitalization of vegetation(Arrhenius)
      IF_VOLITALIZATION_2: IF(MPA_VEG > MPA_VEG_MIN) THEN
        MPA_VOLIT = MAX(CHAR_FCTR*DT_BC*MPA_VEG*A_PYR_VEG*EXP(-E_PYR_VEG/TMP_VEG),0._EB)
        MPA_VOLIT = MIN(MPA_VOLIT,MPA_VOLIT_LOSS_MAX) !user specified max

        DMPA_VEG = CHAR_FCTR2*MPA_VOLIT
        DMPA_VEG = MIN(DMPA_VEG,(MPA_VEG - MPA_VEG_MIN))
        MPA_VEG  = MPA_VEG - DMPA_VEG

        MPA_VOLIT = CHAR_FCTR*DMPA_VEG
        MPA_CHAR  = MPA_CHAR + SF%VEG_CHAR_FRACTION*DMPA_VEG !kg/m^2

      ENDIF IF_VOLITALIZATION_2

      WC%ONE_D%VEG_FUELMASS_L(IVEG_L) = MPA_VEG
      WC%ONE_D%VEG_CHARMASS_L(IVEG_L) = MPA_CHAR

      WC%ONE_D%MASSFLUX(I_FUEL)= WC%ONE_D%MASSFLUX(I_FUEL) + MPA_VOLIT*RDT_BC
      IF (I_WATER > 0) WC%ONE_D%MASSFLUX(I_WATER) = WC%ONE_D%MASSFLUX(I_WATER) + MPA_MOIST_LOSS*RDT_BC

!Char oxidation oF Vegetation Layer within the Arrhenius pyrolysis model
!(note that this can be handled only approximately with the conserved
!scalar based gas-phase combustion model - no gas phase oxygen is consumed by
!the char oxidation reaction since it would be inconsistent with the state
!relation for oxygen based on the conserved scalar approach for gas phase
!combustion)
      IF_CHAR_OXIDATION: IF (SF%VEG_CHAR_OXIDATION .AND. MPA_CHAR > 0.0_EB) THEN
         KKG_L = SF%VEG_KGAS_L(IVEG_L-LBURN)
         ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,ZZ(IIG,JJG,KKG_L,1:N_TRACKED_SPECIES))
         CALL GET_MASS_FRACTION(ZZ_GET,O2_INDEX,Y_O2)
         TMP_G = TMP(IIG,JJG,KKG_L)
         CALL GET_VISCOSITY(ZZ_GET,MU_GAS,TMP_G)
         RE_D = RHO_GAS*SQRT(U2 + V2 + W(IIG,JJG,1)**2)*4._EB/SF%VEG_SV/MU_GAS 
         MPA_CHAR_LOSS = DT_BC*RHO_GAS*Y_O2*A_CHAR_VEG/NU_O2_CHAR_VEG*SF%VEG_SV*  &
                         SF%VEG_PACKING*EXP(-E_CHAR_VEG/WC%ONE_D%VEG_TMP_L(IVEG_L))*  &
                         (1+BETA_CHAR_VEG*SQRT(RE_D))
         MPA_CHAR_LOSS = MIN(MPA_CHAR,MPA_CHAR_LOSS)
         MPA_CHAR      = MPA_CHAR - MPA_CHAR_LOSS
         MPA_ASH       = MPA_ASH + NU_ASH_VEG*MPA_CHAR_LOSS
!        MPA_CHAR_CO2  = (1._EB + NU_O2_CHAR_VEG - NU_ASH_VEG)*MPA_CHAR_LOSS
         WC%ONE_D%VEG_CHARMASS_L(IVEG_L) = MPA_CHAR !kg/m^3
         WC%ONE_D%VEG_ASHMASS_L(IVEG_L)  = MPA_ASH

         IF (MPA_CHAR <= MPA_CHAR_MIN .AND. MPA_VEG <= MPA_VEG_MIN) WC%ONE_D%VEG_CHARMASS_L(IVEG_L) = 0.0_EB
       ENDIF IF_CHAR_OXIDATION

      ENDIF TEMP_THRESEHOLD

! Vegetation temperature (Arrhenius)
      CP_VEG = (0.01_EB + 0.0037_EB*TMP_VEG)*1000._EB !W/kg/K
      CP_CHAR= 420._EB + 2.09_EB*TMP_VEG + 6.85E-4_EB*TMP_VEG**2 !J/kg/K Park etal. C&F 2010 147:481-494
      Q_VEG_CHAR       = MPA_CHAR_LOSS*H_CHAR_VEG
      CP_MOIST_AND_VEG = CP_H2O*WC%ONE_D%VEG_MOISTMASS_L(IVEG_L) + CP_VEG*WC%ONE_D%VEG_FUELMASS_L(IVEG_L) + &
                         CP_CHAR*WC%ONE_D%VEG_CHARMASS_L(IVEG_L) + CP_ASH*WC%ONE_D%VEG_ASHMASS_L(IVEG_L)

      WC%ONE_D%VEG_TMP_L(IVEG_L) = WC%ONE_D%VEG_TMP_L(IVEG_L) + (DT_BC*SF%VEG_DIVQNET_L(IVEG_L-LBURN) - &
                             (MPA_MOIST_LOSS*H_H2O_VEG + MPA_VOLIT*H_PYR_VEG) + CHAR_ENTHALPY_FRACTION_VEG*Q_VEG_CHAR ) &
                             /CP_MOIST_AND_VEG
      WC%ONE_D%VEG_TMP_L(IVEG_L) = MAX( WC%ONE_D%VEG_TMP_L(IVEG_L), TMPA)
      WC%ONE_D%VEG_TMP_L(IVEG_L) = MIN( WC%ONE_D%VEG_TMP_L(IVEG_L), TMP_CHAR_MAX)

    ENDDO LAYER_LOOP2

  ENDIF IF_VEG_DEGRADATION_ARRHENIUS
  
  WC%ONE_D%VEG_TMP_L(LBURN) = MAX(TMP_G,TMPA)
  WC%ONE_D%MASSFLUX_SPEC(I_FUEL) = WC%ONE_D%MASSFLUX(I_FUEL)
  IF (I_WATER > 0) WC%ONE_D%MASSFLUX_SPEC(I_WATER) = WC%ONE_D%MASSFLUX(I_WATER)
 
! Temperature boundary condtions 
! Mass boundary conditions are determine in subroutine SPECIES_BC in wall.f90 for case SPECIFIED_MASS_FLUX
! TMP_F(IW) = WC%VEG_TMP_L(NVEG_L)
! IF (LBURN < NVEG_L)  TMP_F(IW) = WC%VEG_TMP_L(1+LBURN)

  IF (LBURN_NEW < NVEG_L) THEN
    WC%ONE_D%TMP_F = WC%ONE_D%VEG_TMP_L(1+LBURN_NEW)
!   WC%ONE_D%TMP_F = ((VEG_QRP_INC(0)+VEG_QRP_EMISS(0))/SIGMA)**.25 !as done in FDS4
  ELSE
    KKG_L = SF%VEG_KGAS_L(1)
    TMP_G = TMP(IIG,JJG,KKG_L)
    WC%ONE_D%TMP_F = MAX(TMP_G,TMPA) !Tveg=Tgas if veg is completely burned
  ENDIF

!if(wc%one_d%tmp_f > 500._EB) then 
! print '(A,1x,5I3,7ES16.8)','++++ nm,lburn,lburn_new+1,iig,jjg,tmpv_l,tmp_f,tmp_g,vegmass,moistmass,charmass,qnet_veg', &
!    nm,lburn,lburn_new+1,iig,jjg, &
!    wc%veg_tmp_l(lburn),wc%one_d%tmp_f,tmp_g,wc%veg_fuelmass_l(lburn),wc%veg_moistmass_l(lburn), &
!    wc%veg_charmass_l(lburn),qnet_veg
!endif

ENDDO VEG_WALL_CELL_LOOP

VEG_CLOCK_BC = T

END SUBROUTINE BNDRY_VEG_MASS_ENERGY_TRANSFER

SUBROUTINE INITIALIZE_LEVEL_SET_FIRESPREAD_1(NM)

! Set up the major arrays, like the level set value PHI_LS, and determine terrain height on each 2D mesh, Z_LS(I,J).
! After this routine, go back to main, exchange Z_LS, and return for more initialization.

USE COMPLEX_GEOMETRY, ONLY : IBM_IDCF
INTEGER, INTENT(IN) :: NM
INTEGER :: ICF,IW,I,J,SURF_INDEX
TYPE (MESH_TYPE),    POINTER :: M
TYPE (WALL_TYPE),    POINTER :: WC
TYPE (CFACE_TYPE),   POINTER :: CFA
TYPE (SURFACE_TYPE), POINTER :: SF

CALL POINT_TO_MESH(NM)

M => MESHES(NM)

! Loop through all SURFace types and find level set cases that need a calculated RoS

DO SURF_INDEX=0,N_SURF
   SF => SURFACE(SURF_INDEX)
   IF (SF%VEG_LSET_SPREAD .AND. SF%VEG_LSET_FUEL_INDEX>0) THEN
      SF%VEG_LSET_ROS = ROS_NO_WIND_NO_SLOPE(SF%VEG_LSET_FUEL_INDEX,SURF_INDEX)
   ENDIF
ENDDO

! Level set values (Phi). PHI1_LS is the first-order accurate estimate at the next time step.

ALLOCATE(M%PHI_LS(0:IBP1,0:JBP1)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_LS',IZERO)   ; PHI_LS  => M%PHI_LS  ; PHI_LS = PHI_LS_MIN
ALLOCATE(M%PHI1_LS(0:IBP1,0:JBP1)); CALL ChkMemErr('VEGE:LEVEL SET','PHI1_LS',IZERO)  ; PHI1_LS => M%PHI1_LS ; PHI1_LS = PHI_LS_MIN

! Wind speed components in the center of the first gas phsae cell above the ground.

ALLOCATE(M%U_LS(0:IBP1,0:JBP1)) ; CALL ChkMemErr('VEGE:LEVEL SET','U_LS',IZERO) ; U_LS => M%U_LS ; U_LS = 0._EB 
ALLOCATE(M%V_LS(0:IBP1,0:JBP1)) ; CALL ChkMemErr('VEGE:LEVEL SET','V_LS',IZERO) ; V_LS => M%V_LS ; V_LS = 0._EB

! Terrain height, Z_LS, and z index of the first gas cell above terrain, K_LS

ALLOCATE(M%Z_LS(0:IBP1,0:JBP1),STAT=IZERO) ; CALL ChkMemErr('READ','Z_LS',IZERO) ; Z_LS => M%Z_LS ; Z_LS = 0._EB
ALLOCATE(M%K_LS(0:IBP1,0:JBP1),STAT=IZERO) ; CALL ChkMemErr('READ','K_LS',IZERO) ; K_LS => M%K_LS ; K_LS = 0
ALLOCATE(M%LS_SURF_INDEX(0:IBP1,0:JBP1),STAT=IZERO) ; CALL ChkMemErr('READ','LS_SURF_INDEX',IZERO)  
LS_SURF_INDEX => M%LS_SURF_INDEX ; LS_SURF_INDEX = 0

IF (CC_IBM) THEN

   ALLOCATE(M%LS_KLO_TERRAIN(0:IBP1,0:JBP1),STAT=IZERO) ; CALL ChkMemErr('READ','LS_KLO_TERRAIN',IZERO)
   LS_KLO_TERRAIN => M%LS_KLO_TERRAIN ; LS_KLO_TERRAIN = 2*KBP1+1 ! Number larger that KBP1.
   DO ICF=1,M%N_CUTFACE_MESH
      IF (CUT_FACE(ICF)%STATUS /= 2) CYCLE ! IBM_INBOUNDARY == 2
      ! Location of CFACE with largest AREA, to define SURF_INDEX:
      IW  = MAXLOC(CUT_FACE(ICF)%AREA(1:CUT_FACE(ICF)%NFACE),DIM=1)
      CFA => CFACE( CUT_FACE(ICF)%CFACE_INDEX(IW) )
      IF (CFA%NVEC(KAXIS)>-TWO_EPSILON_EB .AND. CFA%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         ! Area averaged Z height of CFACES within this cut-cell (containing IBM_INBOUNDARY CFACES):
         Z_LS(CFA%ONE_D%IIG,CFA%ONE_D%JJG) = DOT_PRODUCT(CUT_FACE(ICF)%XYZCEN(KAXIS,1:CUT_FACE(ICF)%NFACE), &
                                                                   CUT_FACE(ICF)%  AREA(1:CUT_FACE(ICF)%NFACE))     / &
                                                               SUM(CUT_FACE(ICF)% AREA(1:CUT_FACE(ICF)%NFACE))
         IF (CFA%ONE_D%KKG < LS_KLO_TERRAIN(CFA%ONE_D%IIG,CFA%ONE_D%JJG)) & 
            LS_KLO_TERRAIN(CFA%ONE_D%IIG,CFA%ONE_D%JJG) = CFA%ONE_D%KKG
         IF (CFA%ONE_D%KKG > K_LS(CFA%ONE_D%IIG,CFA%ONE_D%JJG)) K_LS(CFA%ONE_D%IIG,CFA%ONE_D%JJG) = CFA%ONE_D%KKG
         LS_SURF_INDEX(CFA%ONE_D%IIG,CFA%ONE_D%JJG) = CFA%SURF_INDEX
      ENDIF
   ENDDO
   DO J=1,JBAR
      DO I=1,IBAR
         IF (K_LS(I,J)==KBAR .AND. FCVAR(I,J,K_LS(I,J),IBM_IDCF,KAXIS)>0) LS_SURF_INDEX(I,J) = 0
      ENDDO
   ENDDO

ELSE

   DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC => WALL(IW)
      IF (WC%ONE_D%IOR==3 .AND. WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
         Z_LS(WC%ONE_D%IIG,WC%ONE_D%JJG) = Z(WC%ONE_D%KKG-1)
         K_LS(WC%ONE_D%IIG,WC%ONE_D%JJG) = WC%ONE_D%KKG
         LS_SURF_INDEX(WC%ONE_D%IIG,WC%ONE_D%JJG)= WC%SURF_INDEX
      ENDIF
   ENDDO

ENDIF

Z_LS(1:IBAR,   0) = 2._EB*Z_LS(1:IBAR,   1) - Z_LS(1:IBAR,   2)
Z_LS(1:IBAR,JBP1) = 2._EB*Z_LS(1:IBAR,JBAR) - Z_LS(1:IBAR,JBM1)
Z_LS(   0,1:JBAR) = 2._EB*Z_LS(   1,1:JBAR) - Z_LS(   2,1:JBAR)
Z_LS(IBP1,1:JBAR) = 2._EB*Z_LS(IBAR,1:JBAR) - Z_LS(IBM1,1:JBAR)

Z_LS(   0,   0) = Z_LS(   1,   1)
Z_LS(IBP1,   0) = Z_LS(IBAR,   1)
Z_LS(   0,JBP1) = Z_LS(   1,JBAR)
Z_LS(IBP1,JBP1) = Z_LS(IBAR,JBAR)

END SUBROUTINE INITIALIZE_LEVEL_SET_FIRESPREAD_1


SUBROUTINE INITIALIZE_LEVEL_SET_FIRESPREAD_2(NM)

! Continuation of set up routine. First, retrieve terrain height, Z_LS, from other meshes. Then do various other set up chores.

INTEGER, INTENT(IN) :: NM
INTEGER :: I,IM1,IM2,IIG,IP1,IP2,J,JJG,JM1,JP1,KDUM,KWIND
REAL(EB) :: DZT_DUM
REAL(EB) :: G_EAST,G_WEST,G_SOUTH,G_NORTH
REAL(EB) :: PHX,PHY,MAG_PHI,UMF_TMP
REAL(EB) :: PHI_W_X,PHI_W_Y,MAG_PHI_S,UMF_X,UMF_Y
TYPE (MESH_TYPE),    POINTER :: M
TYPE (SURFACE_TYPE), POINTER :: SF

CALL POINT_TO_MESH(NM)

M => MESHES(NM)

! Retrieve terrain height, Z_LS, from other meshes above and below current mesh.

CALL GET_BOUNDARY_VALUES

! Allocate some work arrays

ALLOCATE(M%LS_WORK1(0:IBAR,0:JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','LS_WORK1',IZERO)
ALLOCATE(M%LS_WORK2(0:IBAR,0:JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','LS_WORK2',IZERO)

! Define spread rate across domain (including no burn areas)

ALLOCATE(M%ROS_HEAD(IBAR,JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_HEAD',IZERO)  ; ROS_HEAD => M%ROS_HEAD
ALLOCATE(M%ROS_FLANK(IBAR,JBAR))   ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_FLANK',IZERO) ; ROS_FLANK => M%ROS_FLANK
ALLOCATE(M%ROS_BACKU(IBAR,JBAR))   ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_BACKU',IZERO) ; ROS_BACKU => M%ROS_BACKU
ALLOCATE(M%WIND_EXP(IBAR,JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','WIND_EXP',IZERO)  ; WIND_EXP => M%WIND_EXP

! Assign spread rates (i.e., vegetation types) to locations on terrain

ROS_HEAD  = 0.0_EB
ROS_FLANK = 0.0_EB
ROS_BACKU = 0.0_EB
WIND_EXP  = 1.0_EB

LSET_ELLIPSE = .FALSE. ! Flag for the elliptical spread model
LSET_TAN2    = .FALSE. ! Flag for ROS proportional to Tan(slope)^2

! Flux limiters
! LIMITER_LS=1 MINMOD
! LIMITER_LS=2 SUPERBEE
! LIMITER_LS=3 First order upwinding

LIMITER_LS = I_FLUX_LIMITER
IF (LIMITER_LS > 3) LIMITER_LS = 1

! Flux terms

ALLOCATE(M%FLUX0_LS(IBAR,JBAR)); CALL ChkMemErr('VEGE:LEVEL SET','FLUX0_LS',IZERO) ; FLUX0_LS => M%FLUX0_LS
ALLOCATE(M%FLUX1_LS(IBAR,JBAR)); CALL ChkMemErr('VEGE:LEVEL SET','FLUX1_LS',IZERO) ; FLUX1_LS => M%FLUX1_LS

! Slopes (gradients)

ALLOCATE(M%DZTDX(IBAR,JBAR)); CALL ChkMemErr('VEGE:LEVEL SET','DZDTX',IZERO)   ; DZTDX => M%DZTDX
ALLOCATE(M%DZTDY(IBAR,JBAR)); CALL ChkMemErr('VEGE:LEVEL SET','DZDTY',IZERO)   ; DZTDY => M%DZTDY
ALLOCATE(M%MAG_ZT(IBAR,JBAR)); CALL ChkMemErr('VEGE:LEVEL SET','MAG_ZT',IZERO) ; MAG_ZT => M%MAG_ZT

! Rothermel 'Phi' factors for effects of Wind and Slope on ROS

ALLOCATE(M%PHI_WS(IBAR,JBAR))   ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_W',IZERO)   ; PHI_WS => M%PHI_WS    ; PHI_WS = 0.0_EB
ALLOCATE(M%PHI_S(IBAR,JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_S',IZERO)   ; PHI_S => M%PHI_S
ALLOCATE(M%PHI_S_X(IBAR,JBAR))  ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_S_X',IZERO) ; PHI_S_X => M%PHI_S_X
ALLOCATE(M%PHI_S_Y(IBAR,JBAR))  ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_S_Y',IZERO) ; PHI_S_Y => M%PHI_S_Y
ALLOCATE(M%PHI_W(IBAR,JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_W',IZERO)   ; PHI_W => M%PHI_W

! UMF = wind speed at mean flame heights

ALLOCATE(M%UMF(IBAR,JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','UMF',IZERO) ; M%UMF = 0._EB ; UMF => M%UMF
ALLOCATE(M%THETA_ELPS(IBAR,JBAR))    ; CALL ChkMemErr('VEGE:LEVEL SET','THETA_ELPS',IZERO) ; THETA_ELPS => M%THETA_ELPS
THETA_ELPS = 0.0_EB ! Normal to fireline

! ROS in X and Y directions

ALLOCATE(M%SR_X_LS(IBAR,JBAR)) ; CALL ChkMemErr('VEGE:LEVEL SET','SR_X_LS',IZERO) ; SR_X_LS => M%SR_X_LS
ALLOCATE(M%SR_Y_LS(IBAR,JBAR)) ; CALL ChkMemErr('VEGE:LEVEL SET','SR_Y_LS',IZERO) ; SR_Y_LS => M%SR_Y_LS

! Compute components of terrain slope gradient and magnitude of gradient

GRADIENT_ILOOP: DO I = 1,IBAR

   IM1=I-1 ; IM2=I-2
   IP1=I+1 ; IP2=I+2
   IF (I==1) IM1 = I
   IF (I==IBAR) IP1 = I

   DO J = 1,JBAR

      JM1=J-1
      JP1=J+1
      IF (J==1) JM1 = J
      IF (J==IBAR) JP1 = J

      G_EAST  = 0.5_EB*( Z_LS(I,J) + Z_LS(IP1,J) )
      G_WEST  = 0.5_EB*( Z_LS(I,J) + Z_LS(IM1,J) )
      G_NORTH = 0.5_EB*( Z_LS(I,J) + Z_LS(I,JP1) )
      G_SOUTH = 0.5_EB*( Z_LS(I,J) + Z_LS(I,JM1) )

      DZTDX(I,J) = (G_EAST-G_WEST)   * RDX(I)
      DZTDY(I,J) = (G_NORTH-G_SOUTH) * RDY(J)
      MAG_ZT(I,J) = SQRT(DZTDX(I,J)**2 + DZTDY(I,J)**2)

   ENDDO

ENDDO GRADIENT_ILOOP

! Initialize arrays for head, flank, and back fire spread rates with values
! explicitly declared in the input file or from FARSITE head fire and ellipse
! based flank and back fires.
! Fill arrays for the horizontal component of the velocity arrays.
! Initialize level set scalar array PHI

DO JJG=1,JBAR
   DO IIG=1,IBAR

      SF => SURFACE(LS_SURF_INDEX(IIG,JJG))

      IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE

      ! Ignite landscape at user specified location if ignition is at time zero

      IF (SF%VEG_LSET_IGNITE_T == 0.0_EB) PHI_LS(IIG,JJG) = PHI_LS_MAX

      ! Wind field

      U_LS(IIG,JJG) = 0.5_EB*(U(IIG-1,JJG,K_LS(IIG,JJG))+U(IIG,JJG,K_LS(IIG,JJG)))
      V_LS(IIG,JJG) = 0.5_EB*(V(IIG,JJG-1,K_LS(IIG,JJG))+V(IIG,JJG,K_LS(IIG,JJG)))

      ROS_HEAD(IIG,JJG)  = SF%VEG_LSET_ROS_HEAD
      ROS_FLANK(IIG,JJG) = SF%VEG_LSET_ROS_FLANK
      ROS_BACKU(IIG,JJG) = SF%VEG_LSET_ROS_BACK
      WIND_EXP(IIG,JJG)  = SF%VEG_LSET_WIND_EXP

      ! If any surfaces uses tan^2 function for slope, tan^2 will be used throughout simulation

      IF (SF%VEG_LSET_TAN2) LSET_TAN2=.TRUE.

      ! Use assumed ellipse shape of fireline as in Farsite

      IF_ELLIPSE_UNCOUPLED: IF (SF%VEG_LSET_ELLIPSE) THEN

         ROS_HEAD(IIG,JJG) = SF%VEG_LSET_ROS
         SF%VEG_LSET_HT = MAX(0.001_EB,SF%VEG_LSET_HT)

         ! If any surfaces set to ellipse, then elliptical model used for all surfaces

         IF (.NOT. LSET_ELLIPSE) LSET_ELLIPSE = .TRUE.

         ! Find wind at ~6.1 m height for Farsite
         KWIND = 0
         DO KDUM = K_LS(IIG,JJG),KBAR
            IF (ZC(KDUM)-ZC(K_LS(IIG,JJG)) >= 6.1_EB) KWIND = KDUM
         ENDDO
         IF (ZC(KBAR) < 6.1_EB) KWIND=1

         U_LS(IIG,JJG) = 0.5_EB*(U(IIG-1,JJG,KWIND)+U(IIG,JJG,KWIND))
         V_LS(IIG,JJG) = 0.5_EB*(V(IIG,JJG-1,KWIND)+V(IIG,JJG,KWIND))

         ! Wind at midflame height (UMF). From Andrews 2012, USDA FS Gen Tech Rep. RMRS-GTR-266 (with added SI conversion)

         UMF_TMP = 1.83_EB / LOG((20.0_EB + 1.18_EB * SF%VEG_LSET_HT) /(0.43_EB * SF%VEG_LSET_HT))

         ! Factor 60 converts U from m/s to m/min which is used in elliptical model.
         UMF_X = UMF_TMP * U_LS(IIG,JJG) * 60.0_EB
         UMF_Y = UMF_TMP * V_LS(IIG,JJG) * 60.0_EB

         ! Variables used in Phi_W formulas below (Rothermel model)
         B_ROTH = 0.15988_EB * (SF%VEG_LSET_SIGMA**0.54_EB)
         C_ROTH = 7.47_EB * EXP(-0.8711_EB * (SF%VEG_LSET_SIGMA**0.55_EB))
         E_ROTH = 0.715_EB * EXP(-0.01094_EB * SF%VEG_LSET_SIGMA)
         BETA_OP_ROTH = 0.20395_EB * (SF%VEG_LSET_SIGMA**(-0.8189_EB))! Optimum packing ratio

         ! Components of wind factor - affects spread rate
         PHI_W_X = C_ROTH * ((3.281_EB * ABS(UMF_X))**B_ROTH) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**(-E_ROTH)
         PHI_W_X = SIGN(PHI_W_X,UMF_X)

         PHI_W_Y = C_ROTH * ((3.281_EB * ABS(UMF_Y))**B_ROTH) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**(-E_ROTH)
         PHI_W_Y = SIGN(PHI_W_Y,UMF_Y)

         PHI_W(IIG,JJG) =  SQRT(PHI_W_X**2 + PHI_W_Y**2)

         ! Limit effect to slope lte 80 degrees. Phi_s_x,y are slope factors
         DZT_DUM = MIN(5.67_EB,ABS(DZTDX(IIG,JJG))) ! 5.67 ~ tan 80 deg
         PHI_S_X(IIG,JJG) = 5.275_EB * ((SF%VEG_LSET_BETA)**(-0.3_EB)) * DZT_DUM**2
         PHI_S_X(IIG,JJG) = SIGN(PHI_S_X(IIG,JJG),DZTDX(IIG,JJG))
         DZT_DUM = MIN(1.73_EB,ABS(DZTDY(IIG,JJG))) ! 1.73 ~ tan 60 deg
         PHI_S_Y(IIG,JJG) = 5.275_EB * ((SF%VEG_LSET_BETA)**(-0.3_EB)) * DZT_DUM**2
         PHI_S_Y(IIG,JJG) = SIGN(PHI_S_Y(IIG,JJG),DZTDY(IIG,JJG))

         MAG_PHI_S = SQRT(PHI_S_X(IIG,JJG)**2 + PHI_S_Y(IIG,JJG)**2)

         PHI_S(IIG,JJG) = MAG_PHI_S  !5.275 * MAG_ZT(I,J)**2 * (SF%VEG_LSET_BETA)**-0.3

         ! Slope factor

         IF (MAG_PHI_S > 0.0_EB) THEN

            PHX = PHI_W_X + PHI_S_X(IIG,JJG)
            PHY = PHI_W_Y + PHI_S_Y(IIG,JJG)
            MAG_PHI = SQRT(PHX**2 + PHY**2)

            ! Total phi (phi_w + phi_s) for use in spread rate section
            PHI_WS(IIG,JJG) = MAG_PHI

            ! Theta_elps is angle of direction (0 to 2pi) of highest spread rate
            ! 0<=theta_elps<=2pi as measured clockwise from Y-axis
            THETA_ELPS(IIG,JJG) = ATAN2(PHY,PHX)

            !"Effective midflame windspeed" used in length-to-breadth ratio calculation (spread rate routine)
            ! is the wind + slope effect obtained by solving Phi_w eqs. above for UMF
            ! 8/8/13 - Changed phi_ws to Phi_s below to match Farsite, i.e., instead of adding phi_w and phi_s
            ! and then calculating effective wind speed, phi_s is converted to an effected wind speed and added
            ! to UMF calculated from the wind. Effective U has units of m/min in Wilson formula.
            ! 0.3048 ~= 1/3.281
            ! if phi_s < 0 then a complex value (NaN) results. Using abs(phi_s) and sign function to correct.

            UMF_TMP = (((ABS(PHI_S_X(IIG,JJG)) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**E_ROTH)/C_ROTH)**(1/B_ROTH))*0.3048
            UMF_TMP = SIGN(UMF_TMP,PHI_S_X(IIG,JJG))
            UMF_X = UMF_X + UMF_TMP

            UMF_TMP = (((ABS(PHI_S_Y(IIG,JJG)) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**E_ROTH)/C_ROTH)**(1/B_ROTH))*0.3048
            UMF_TMP = SIGN(UMF_TMP,PHI_S_Y(IIG,JJG))
            UMF_Y = UMF_Y + UMF_TMP

         ELSE

            PHI_WS(IIG,JJG) = SQRT (PHI_W_X**2 + PHI_W_Y**2)
            THETA_ELPS(IIG,JJG) = ATAN2(PHI_W_Y,PHI_W_X)

         ENDIF

         UMF(IIG,JJG) = SQRT(UMF_X**2 + UMF_Y**2)

         ! The following two lines convert ATAN2 output to compass system (0 to 2 pi CW from +Y-axis)
         THETA_ELPS(IIG,JJG) = PIO2 - THETA_ELPS(IIG,JJG)
         IF (THETA_ELPS(IIG,JJG) < 0.0_EB) THETA_ELPS(IIG,JJG) = 2.0_EB*PI + THETA_ELPS(IIG,JJG)

      ENDIF IF_ELLIPSE_UNCOUPLED

   ENDDO
ENDDO

END SUBROUTINE INITIALIZE_LEVEL_SET_FIRESPREAD_2


SUBROUTINE LEVEL_SET_FIRESPREAD(T,DT,NM)

! Predictor: Estimate PHI_LS at next time step. Estimated value is called PHI1_LS.
! Corrector: Correct PHI_LS at next time step.

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT
INTEGER :: IIG,IW,JJG,IC
INTEGER :: KDUM,KWIND,ICF,IKT
REAL(EB) :: UMF_TMP,PHX,PHY,MAG_PHI,PHI_W_X,PHI_W_Y,UMF_X,UMF_Y,UMAG
TYPE (ONE_D_M_AND_E_XFER_TYPE), POINTER :: ONE_D
TYPE (SURFACE_TYPE), POINTER :: SF

CALL POINT_TO_MESH(NM)

CALL GET_BOUNDARY_VALUES

! Loop over terrain surface cells and update level set field

DO JJG=1,JBAR
   DO IIG=1,IBAR

      SF => SURFACE(LS_SURF_INDEX(IIG,JJG))

      IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE

      ! Ignite landscape at user specified location(s) and time(s)

      IF (SF%VEG_LSET_IGNITE_T > 0.0_EB .AND. SF%VEG_LSET_IGNITE_T < DT)  PHI_LS(IIG,JJG) = PHI_LS_MAX
      IF (SF%VEG_LSET_IGNITE_T >= T .AND. SF%VEG_LSET_IGNITE_T <= T + DT) PHI_LS(IIG,JJG) = PHI_LS_MAX

      ! Variable update when level set is coupled to CFD computation

      IF_CFD_COUPLED: IF (VEG_LEVEL_SET_COUPLED) THEN

         U_LS(IIG,JJG) = 0.5_EB*(U(IIG-1,JJG,K_LS(IIG,JJG))+U(IIG,JJG,K_LS(IIG,JJG)))
         V_LS(IIG,JJG) = 0.5_EB*(V(IIG,JJG-1,K_LS(IIG,JJG))+V(IIG,JJG,K_LS(IIG,JJG)))

         ! AU grassland ROS for infinite head and 6% moisutre

         UMAG     = SQRT(U_LS(IIG,JJG)**2 + V_LS(IIG,JJG)**2)
         ROS_HEAD(IIG,JJG)  = SF%VEG_LSET_ROS_HEAD*(0.165_EB + 0.534_EB*UMAG)*0.523_EB

         ! Use assumed ellipse shape of fireline as in Farsite

         IF_ELLIPSE_COUPLED: IF (SF%VEG_LSET_ELLIPSE) THEN

            ROS_HEAD(IIG,JJG) = SF%VEG_LSET_ROS

            ! Find wind at ~6.1 m height for Farsite
            KWIND = 0
            DO KDUM = K_LS(IIG,JJG),KBAR
              IF(ZC(KDUM)-ZC(K_LS(IIG,JJG)) >= 6.1_EB) THEN
               KWIND = KDUM
               ENDIF
            ENDDO

            U_LS(IIG,JJG) = 0.5_EB*(U(IIG-1,JJG,KWIND)+U(IIG,JJG,KWIND))
            V_LS(IIG,JJG) = 0.5_EB*(V(IIG,JJG-1,KWIND)+V(IIG,JJG,KWIND))

            ! Wind at midflame height (UMF). From Andrews 2012, USDA FS Gen Tech Rep. RMRS-GTR-266 (with added SI conversion)
            UMF_TMP = 1.83_EB / LOG((20.0_EB + 1.18_EB * SF%VEG_LSET_HT) /(0.43_EB * SF%VEG_LSET_HT))

            ! Factor 60 converts U from m/s to m/min which is used in elliptical model.
            UMF_X = UMF_TMP * U_LS(IIG,JJG) * 60.0_EB
            UMF_Y = UMF_TMP * V_LS(IIG,JJG) * 60.0_EB

            ! Components of wind factor - affects spread rate
            PHI_W_X = C_ROTH * ((3.281_EB * ABS(UMF_X))**B_ROTH) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**(-E_ROTH)
            PHI_W_X = SIGN(PHI_W_X,UMF_X)

            PHI_W_Y = C_ROTH * ((3.281_EB * ABS(UMF_Y))**B_ROTH) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**(-E_ROTH)
            PHI_W_Y = SIGN(PHI_W_Y,UMF_Y)

            PHI_W(IIG,JJG) =  SQRT(PHI_W_X**2 + PHI_W_Y**2)

            ! Slope factor

            IF (PHI_S(IIG,JJG) > 0.0_EB) THEN

               PHX = PHI_W_X + PHI_S_X(IIG,JJG)
               PHY = PHI_W_Y + PHI_S_Y(IIG,JJG)
               MAG_PHI = SQRT(PHX**2 + PHY**2)

               ! Total phi (phi_w + phi_s) for use in spread rate section
               PHI_WS(IIG,JJG) = MAG_PHI

               ! Theta_elps is angle of direction (0 to 2pi) of highest spread rate
               ! 0<=theta_elps<=2pi as measured clockwise from Y-axis
               THETA_ELPS(IIG,JJG) = ATAN2(PHY,PHX)

               !"Effective midflame windspeed" used in length-to-breadth ratio calculation (spread rate routine)
               ! is the wind + slope effect obtained by solving Phi_w eqs. above for UMF
               ! 8/8/13 - Changed phi_ws to Phi_s below to match Farsite, i.e., instead of adding phi_w and phi_s
               ! and then calculating effective wind speed, phi_s is converted to an effected wind speed and added
               ! to UMF calculated from the wind. Effective U has units of m/min in Wilson formula.
               ! 0.3048 ~= 1/3.281
               !if phi_s < 0 then a complex value (NaN) results. Using abs(phi_s) and sign function to correct.

               UMF_TMP = (((ABS(PHI_S_X(IIG,JJG)) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**E_ROTH)/C_ROTH)**(1/B_ROTH))*0.3048
               UMF_TMP = SIGN(UMF_TMP,PHI_S_X(IIG,JJG))
               UMF_X = UMF_X + UMF_TMP

               UMF_TMP = (((ABS(PHI_S_Y(IIG,JJG)) * (SF%VEG_LSET_BETA / BETA_OP_ROTH)**E_ROTH)/C_ROTH)**(1/B_ROTH))*0.3048
               UMF_TMP = SIGN(UMF_TMP,PHI_S_Y(IIG,JJG))
               UMF_Y = UMF_Y + UMF_TMP

            ELSE

               PHI_WS(IIG,JJG) = SQRT(PHI_W_X**2 + PHI_W_Y**2)
               THETA_ELPS(IIG,JJG) = ATAN2(PHI_W_Y,PHI_W_X)

            ENDIF

            UMF(IIG,JJG) = SQRT(UMF_X**2 + UMF_Y**2)

            ! The following two lines convert ATAN2 output to compass system (0 to 2 pi CW from +Y-axis)

            THETA_ELPS(IIG,JJG) = PIO2 - THETA_ELPS(IIG,JJG)
            IF (THETA_ELPS(IIG,JJG) < 0.0_EB) THETA_ELPS(IIG,JJG) = 2.0_EB*PI + THETA_ELPS(IIG,JJG)

         ENDIF IF_ELLIPSE_COUPLED

      ENDIF IF_CFD_COUPLED

   ENDDO
ENDDO

! Runge-Kutta Scheme

CALL LEVEL_SET_SPREAD_RATE
CALL LEVEL_SET_ADVECT_FLUX

IF (PREDICTOR) THEN
   PHI1_LS(1:IBAR,1:JBAR) = PHI_LS(1:IBAR,1:JBAR) - DT*FLUX0_LS(1:IBAR,1:JBAR)
   PHI1_LS = MAX(PHI_LS_MIN,MIN(PHI_LS_MAX,PHI1_LS))
ELSE
   PHI_LS(1:IBAR,1:JBAR) = PHI_LS(1:IBAR,1:JBAR) - 0.5_EB*DT*(FLUX0_LS(1:IBAR,1:JBAR) + FLUX1_LS(1:IBAR,1:JBAR))
   PHI_LS = MAX(PHI_LS_MIN,MIN(PHI_LS_MAX,PHI_LS))
ENDIF

! Loop over all cells and assign the value of PHI_LS to the appropriate WALL or
! CFACE cells. Also, if PHI_LS increases above 0, set the ignition time T_IGN.

IF (.NOT.PREDICTOR) THEN
   IF (CC_IBM) THEN
      DO JJG=1,JBAR
         DO IIG=1,IBAR
            DO IKT=LS_KLO_TERRAIN(IIG,JJG),K_LS(IIG,JJG)
               ! Loop over all CFACEs corresponding to IIG,JJG and set ONE_D%T_IGN and ONE_D%PHI_LS as below
               ICF = CCVAR(IIG,JJG,IKT,3); IF(ICF<1) CYCLE  ! IBM_IDCF = 3 CUT_FCE container for this cell.
               DO IW=1,CUT_FACE(ICF)%NFACE ! All IBM_INBOUNDARY CFACES on this cell.
                  ONE_D => CFACE( CUT_FACE(ICF)%CFACE_INDEX(IW) ) % ONE_D
                  IF (PHI_LS(IIG,JJG)>=0._EB .AND. ONE_D%T_IGN>1.E5_EB) ONE_D%T_IGN = T
                  ONE_D%PHI_LS = PHI_LS(IIG,JJG)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO JJG=1,JBAR
         DO IIG=1,IBAR
            IF (K_LS(IIG,JJG)<1) CYCLE
            IC = CELL_INDEX(IIG,JJG,K_LS(IIG,JJG))
            IW = WALL_INDEX(IC,-3)
            ONE_D => WALL(IW)%ONE_D
            IF (PHI_LS(IIG,JJG)>=0._EB .AND. ONE_D%T_IGN>1.E5_EB) ONE_D%T_IGN = T
            ONE_D%PHI_LS = PHI_LS(IIG,JJG)
         ENDDO
      ENDDO
   ENDIF
ENDIF

END SUBROUTINE LEVEL_SET_FIRESPREAD


SUBROUTINE GET_BOUNDARY_VALUES

! Retrieve various quantities from neighboring meshes.

INTEGER :: IIG,JJG,II,JJ,IOR

IF (PREDICTOR) THEN
   PHI_LS_P => PHI_LS
ELSE
   PHI_LS_P => PHI1_LS
ENDIF

! Fetch values of various quantities at six faces of current mesh.

DO II=1,IBAR
   IIG=II ; JJ=0 ; JJG=1 ; IOR=-2
   CALL FILL_BOUNDARY_VALUES
   IIG=II ; JJ=JBP1 ; JJG=JBAR ; IOR=2
   CALL FILL_BOUNDARY_VALUES
ENDDO
DO JJ=1,JBAR
   JJG=JJ ; II=0 ; IIG=1 ; IOR=-1
   CALL FILL_BOUNDARY_VALUES
   JJG=JJ ; II=IBP1 ; IIG=IBAR ; IOR=1
   CALL FILL_BOUNDARY_VALUES
ENDDO
DO JJ=1,JBAR
   JJG=JJ
   DO II=1,IBAR
      IIG=II
      IOR = -3
      CALL FILL_BOUNDARY_VALUES
      IOR =  3
      CALL FILL_BOUNDARY_VALUES
   ENDDO
ENDDO

CONTAINS

SUBROUTINE FILL_BOUNDARY_VALUES

USE COMPLEX_GEOMETRY, ONLY : IBM_CGSC,IBM_SOLID,IBM_CUTCFE
INTEGER :: IW,IIO,JJO,N_INT_CELLS,NOM,IC
REAL(EB) :: PHI_LS_OTHER,U_LS_OTHER,V_LS_OTHER,Z_LS_OTHER
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC
LOGICAL :: SOLID_CELL

! Grab boundary values of PHI_LS from MPI storage arrays

IF (IOR==3) THEN  ! get the CELL_INDEX of the grid cell adjacent to the exterior boundary of the current mesh
   IC = CELL_INDEX(IIG,JJG,KBAR)
ELSEIF (IOR==-3) THEN
   IC = CELL_INDEX(IIG,JJG,1)
ELSE
   IC = CELL_INDEX(IIG,JJG,1)
ENDIF
IW = WALL_INDEX(IC,IOR)
EWC=>EXTERNAL_WALL(IW)
NOM = EWC%NOM
IF (NOM==0) RETURN  ! there is no other mesh adjacent to the boundary

PHI_LS_OTHER = 0._EB
U_LS_OTHER = 0._EB
V_LS_OTHER = 0._EB
Z_LS_OTHER = 0._EB
DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
   DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
      IF (PREDICTOR) THEN
         PHI_LS_OTHER = PHI_LS_OTHER + OMESH(NOM)%PHI_LS(IIO,JJO)
      ELSE
         PHI_LS_OTHER = PHI_LS_OTHER + OMESH(NOM)%PHI1_LS(IIO,JJO)
      ENDIF
      U_LS_OTHER = U_LS_OTHER + OMESH(NOM)%U_LS(IIO,JJO)
      V_LS_OTHER = V_LS_OTHER + OMESH(NOM)%V_LS(IIO,JJO)
      Z_LS_OTHER = Z_LS_OTHER + OMESH(NOM)%Z_LS(IIO,JJO)
   ENDDO
ENDDO
N_INT_CELLS = (EWC%IIO_MAX-EWC%IIO_MIN+1) * (EWC%JJO_MAX-EWC%JJO_MIN+1)

SELECT CASE(IOR)
   CASE(-2:2) 
      PHI_LS_P(II,JJ) = PHI_LS_OTHER/REAL(N_INT_CELLS,EB)
      U_LS(II,JJ)     = U_LS_OTHER/REAL(N_INT_CELLS,EB)
      V_LS(II,JJ)     = V_LS_OTHER/REAL(N_INT_CELLS,EB)
      Z_LS(II,JJ)     = Z_LS_OTHER/REAL(N_INT_CELLS,EB)
   CASE(3)  ! only grab a PHI_LS value from the other mesh if the (II,JJ) cell of the current mesh has no terrain surface
      SOLID_CELL = .FALSE.
      IF (CC_IBM) THEN
         IF (CCVAR(II,JJ,KBAR,IBM_CGSC)==IBM_SOLID .OR. CCVAR(II,JJ,KBAR,IBM_CGSC)==IBM_CUTCFE) SOLID_CELL = .TRUE.
      ELSE
         IF (SOLID(CELL_INDEX(II,JJ,KBAR))) SOLID_CELL = .TRUE.
      ENDIF
      IF (.NOT.SURFACE(LS_SURF_INDEX(II,JJ))%VEG_LSET_SPREAD .AND. SOLID_CELL) THEN
         PHI_LS_P(II,JJ) = PHI_LS_OTHER/REAL(N_INT_CELLS,EB)
         U_LS(II,JJ) = U_LS_OTHER/REAL(N_INT_CELLS,EB)
         V_LS(II,JJ) = V_LS_OTHER/REAL(N_INT_CELLS,EB)
         Z_LS(II,JJ) = Z_LS_OTHER/REAL(N_INT_CELLS,EB)
      ENDIF
   CASE(-3)  ! only grab a PHI_LS value from the other mesh if the (II,JJ) cell of the current mesh has no terrain surface
      SOLID_CELL = .FALSE.
      IF (CC_IBM) THEN
         IF (CCVAR(II,JJ,1,IBM_CGSC)==IBM_SOLID .OR. CCVAR(II,JJ,1,IBM_CGSC)==IBM_CUTCFE) SOLID_CELL = .TRUE.
      ELSE
         IF (SOLID(CELL_INDEX(II,JJ,1))) SOLID_CELL = .TRUE.
      ENDIF
      IF (.NOT.SURFACE(LS_SURF_INDEX(II,JJ))%VEG_LSET_SPREAD .AND. .NOT.SOLID_CELL) THEN
         PHI_LS_P(II,JJ) = PHI_LS_OTHER/REAL(N_INT_CELLS,EB)
         U_LS(II,JJ) = U_LS_OTHER/REAL(N_INT_CELLS,EB)
         V_LS(II,JJ) = V_LS_OTHER/REAL(N_INT_CELLS,EB)
         Z_LS(II,JJ) = Z_LS_OTHER/REAL(N_INT_CELLS,EB)
      ENDIF
END SELECT

END SUBROUTINE FILL_BOUNDARY_VALUES

END SUBROUTINE GET_BOUNDARY_VALUES


SUBROUTINE LEVEL_SET_SPREAD_RATE

! Compute components of spread rate vector

INTEGER :: I,J,IM1,IP1,JM1,JP1
REAL(EB) :: COS_THETA_WIND,COS_THETA_SLOPE,COS_THETA_WIND_H,COS_THETA_WIND_B, &
            COS_THETA_SLOPE_H,COS_THETA_SLOPE_B,DPHIDX,DPHIDY,F_EAST,F_WEST,F_NORTH,F_SOUTH, &
            GRAD_SLOPE_DOT_NORMAL_FIRELINE,MAG_F,MAG_SR,MAG_U,WIND_DOT_NORMAL_FIRELINE,NEXP_WIND
REAL(EB) :: ROS_BACKS,ROS_HEADS
REAL(EB) :: RAD_TO_DEGREE,DEGREES_SLOPE,SLOPE_FACTOR
REAL(EB) :: COS_THETA,SIN_THETA,XSF,YSF,UMF_DUM
REAL(EB) :: AROS,A_ELPS,A_ELPS2,BROS,B_ELPS2,B_ELPS,C_ELPS,DENOM,ROS_TMP,LB,LBD,HB
REAL(EB), DIMENSION(:) :: NORMAL_FIRELINE(2)

RAD_TO_DEGREE = 90._EB/ASIN(1._EB)

IF (PREDICTOR) THEN
   PHI_LS_P => PHI_LS
ELSE
   PHI_LS_P => PHI1_LS
ENDIF

SR_X_LS = 0.0_EB ; SR_Y_LS = 0.0_EB

FLUX_ILOOP: DO J=1,JBAR

   JM1 = J-1
   JP1 = J+1

   DO I=1,IBAR

      IM1 = I-1
      IP1 = I+1

      F_EAST  = 0.5_EB*( PHI_LS_P(I,J) + PHI_LS_P(IP1,J) )
      F_WEST  = 0.5_EB*( PHI_LS_P(I,J) + PHI_LS_P(IM1,J) )
      F_NORTH = 0.5_EB*( PHI_LS_P(I,J) + PHI_LS_P(I,JP1) )
      F_SOUTH = 0.5_EB*( PHI_LS_P(I,J) + PHI_LS_P(I,JM1) )

      DPHIDX = (F_EAST-F_WEST)   * RDX(I)
      DPHIDY = (F_NORTH-F_SOUTH) * RDY(J)

      MAG_F = SQRT(DPHIDX**2 + DPHIDY**2)

      IF (MAG_F > 0._EB) THEN   !components of unit vector normal to PHI contours
         NORMAL_FIRELINE(1) = -DPHIDX/MAG_F
         NORMAL_FIRELINE(2) = -DPHIDY/MAG_F
         XSF =  DPHIDY
         YSF = -DPHIDX
         GRAD_SLOPE_DOT_NORMAL_FIRELINE = DZTDX(I,J)*(DPHIDY/MAG_F) + DZTDY(I,J)*(-DPHIDY/MAG_F)
      ELSE
        NORMAL_FIRELINE = 0._EB
        GRAD_SLOPE_DOT_NORMAL_FIRELINE = 0._EB
        XSF=0._EB
        YSF=0._EB
      ENDIF

      COS_THETA_SLOPE = 0.0_EB ; COS_THETA_SLOPE_H = 0.0_EB ; COS_THETA_SLOPE_B = 0.0_EB

      IF (MAG_ZT(I,J) > 0.0_EB) COS_THETA_SLOPE = GRAD_SLOPE_DOT_NORMAL_FIRELINE/MAG_ZT(I,J)

      DEGREES_SLOPE = ATAN(MAG_ZT(I,J))*RAD_TO_DEGREE

      IF (LSET_ELLIPSE) THEN

         ! Effective wind direction (theta) is clockwise from y-axis (Richards 1990)
         COS_THETA = COS(THETA_ELPS(I,J)) !V_LS(I,J) / MAG_U
         SIN_THETA = SIN(THETA_ELPS(I,J)) !U_LS(I,J) / MAG_U

         ROS_TMP = ROS_HEAD(I,J) * (1.0_EB + PHI_WS(I,J))

         ! Mag of wind speed at midflame ht must be in units of m/s here
         UMF_DUM = UMF(I,J)/60.0_EB

         ! Length to breadth ratio of ellipse based on effective UMF
         LB = 0.936_EB * EXP(0.2566_EB * UMF_DUM) + 0.461_EB * EXP(-0.1548_EB * UMF_DUM) - 0.397_EB

         ! Constraint LB max = 8 from Finney 2004
         LB = MAX(1.0_EB,MIN(LB,8.0_EB))
         LBD = SQRT(LB**2 - 1.0_EB)

         ! Head to back ratio based on LB
         HB = (LB + LBD) / (LB - LBD)

         ! A_ELPS and B_ELPS notation is consistent with Farsite and Richards 
         B_ELPS =  0.5_EB * (ROS_TMP + ROS_TMP/HB)
         B_ELPS2 = B_ELPS**2
         A_ELPS =  B_ELPS / LB
         A_ELPS2=  A_ELPS**2
         C_ELPS =  B_ELPS - (ROS_TMP/HB)

         ! Denominator used in spread rate equation from Richards, Intnl. J. Num. Methods Eng. 1990 
         ! and in LS vs Farsite paper, Bova et al., Intnl. J. Wildland Fire, 25(2):229-241, 2015  
         AROS  = XSF*COS_THETA - YSF*SIN_THETA
         BROS  = XSF*SIN_THETA + YSF*COS_THETA
         DENOM = A_ELPS2*BROS**2 + B_ELPS2*AROS**2

         IF (DENOM > 0._EB) THEN
            DENOM = 1._EB / SQRT(DENOM)
         ELSE
            DENOM = 0._EB
         ENDIF

!        This is with A_ELPS2 and B_ELPS2 notation consistent with Finney and Richards and in 
!        Bova et al. 2015 IJWF 2015
         SR_X_LS(I,J) = DENOM * ( A_ELPS2*COS_THETA*BROS - B_ELPS2*SIN_THETA*AROS) + C_ELPS*SIN_THETA
         SR_Y_LS(I,J) = DENOM * (-A_ELPS2*SIN_THETA*BROS - B_ELPS2*COS_THETA*AROS) + C_ELPS*COS_THETA

         ! Project spread rates from slope to horizontal plane

         IF (ABS(DZTDX(I,J)) > 0._EB) SR_X_LS(I,J) = SR_X_LS(I,J) * ABS(COS(ATAN(DZTDX(I,J))))
         IF (ABS(DZTDY(I,J)) > 0._EB) SR_Y_LS(I,J) = SR_Y_LS(I,J) * ABS(COS(ATAN(DZTDY(I,J))))

         MAG_SR = SQRT(SR_X_LS(I,J)**2 + SR_Y_LS(I,J)**2)

      ELSE ! McArthur Spread Model

         WIND_DOT_NORMAL_FIRELINE = U_LS(I,J)*NORMAL_FIRELINE(1) + V_LS(I,J)*NORMAL_FIRELINE(2)
         MAG_U  = SQRT(U_LS(I,J)**2 + V_LS(I,J)**2)

         COS_THETA_WIND = 0.0_EB ; COS_THETA_WIND_H = 0.0_EB ; COS_THETA_WIND_B = 0.0_EB
         IF(MAG_U > 0.0_EB) COS_THETA_WIND = WIND_DOT_NORMAL_FIRELINE/MAG_U

         GRAD_SLOPE_DOT_NORMAL_FIRELINE = DZTDX(I,J)*NORMAL_FIRELINE(1) + DZTDY(I,J)*NORMAL_FIRELINE(2)
         COS_THETA_SLOPE = 0.0_EB ; COS_THETA_SLOPE_H = 0.0_EB ; COS_THETA_SLOPE_B = 0.0_EB

         IF (MAG_ZT(I,J) > 0.0_EB) COS_THETA_SLOPE = GRAD_SLOPE_DOT_NORMAL_FIRELINE/MAG_ZT(I,J)

         DEGREES_SLOPE = ATAN(MAG_ZT(I,J))*RAD_TO_DEGREE

         SLOPE_FACTOR  = MAG_ZT(I,J)**2
         IF (SLOPE_FACTOR > 3._EB) SLOPE_FACTOR = 3._EB

         ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)

         IF (DEGREES_SLOPE >= 5._EB  .AND. DEGREES_SLOPE < 10._EB) ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)
         IF (DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =         ROS_HEAD(I,J)
         IF (DEGREES_SLOPE >= 20._EB)                              ROS_HEADS = 3._EB*  ROS_HEAD(I,J)

         MAG_SR    = 0.0_EB
         ROS_HEADS = 0.0_EB
         ROS_BACKS = 0.0_EB

         NEXP_WIND = WIND_EXP(I,J)

         ! Spread with the wind and upslope

         IF (COS_THETA_WIND >= 0._EB .AND. COS_THETA_SLOPE >= 0._EB) THEN
            IF (.NOT. LSET_TAN2) THEN
                IF (DEGREES_SLOPE >= 5._EB  .AND. DEGREES_SLOPE < 10._EB) ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)
                IF (DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =         ROS_HEAD(I,J)
                IF (DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  3._EB*ROS_HEAD(I,J)
            ELSEIF (DEGREES_SLOPE > 0._EB) THEN
                ROS_HEADS = ROS_HEAD(I,J) * SLOPE_FACTOR !Dependence on TAN(slope)^2
            ENDIF
            MAG_SR = ROS_FLANK(I,J)*(1._EB + COS_THETA_WIND**NEXP_WIND*COS_THETA_SLOPE) + &
                     (ROS_HEAD(I,J) - ROS_FLANK(I,J))*COS_THETA_WIND**NEXP_WIND + &
                     (ROS_HEADS     - ROS_FLANK(I,J))*COS_THETA_SLOPE  !magnitude of spread rate
         ENDIF

         ! Spread with the wind and downslope

         IF (COS_THETA_WIND >= 0._EB .AND. COS_THETA_SLOPE < 0._EB) THEN
            IF (DEGREES_SLOPE >= 5._EB  .AND. DEGREES_SLOPE < 10._EB) ROS_HEADS =  0.33_EB*ROS_HEAD(I,J)
            IF (DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =  0.50_EB*ROS_HEAD(I,J)
            IF (DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  0.75_EB*ROS_HEAD(I,J)
            MAG_SR = ROS_FLANK(I,J)*(1._EB + COS_THETA_WIND*COS_THETA_SLOPE) + &
                     (ROS_HEAD(I,J) - ROS_FLANK(I,J))*COS_THETA_WIND**NEXP_WIND + &
                     (ROS_HEADS     - ROS_FLANK(I,J))*COS_THETA_SLOPE  !magnitude of spread rate
         ENDIF

         ! Spread against the wind and upslope

         IF (COS_THETA_WIND <  0._EB .AND. COS_THETA_SLOPE >= 0._EB) THEN
            IF (.NOT. LSET_TAN2) THEN
                IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_BACKS = -0.33_EB*ROS_BACKU(I,J)
                IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_BACKS =         -ROS_BACKU(I,J)
                IF(DEGREES_SLOPE >= 20._EB)                              ROS_BACKS = -3.0_EB*ROS_BACKU(I,J)
            ELSEIF (DEGREES_SLOPE > 0._EB) THEN
                ROS_HEADS = ROS_HEAD(I,J) * SLOPE_FACTOR !Dependence on TAN(slope)^2
            ENDIF
            MAG_SR = ROS_FLANK(I,J)*(1._EB - ABS(COS_THETA_WIND)**NEXP_WIND*COS_THETA_SLOPE) + &
                     (ROS_FLANK(I,J) - ROS_BACKU(I,J))*(-ABS(COS_THETA_WIND)**NEXP_WIND) + &
                     (ROS_FLANK(I,J) - ROS_BACKS)*COS_THETA_SLOPE  !magnitude of spread rate
         ENDIF

         ! Spread against the wind and downslope

         IF (COS_THETA_WIND <  0._EB .AND. COS_THETA_SLOPE < 0._EB) THEN
            IF (DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_BACKS = 0.33_EB*ROS_BACKU(I,J)
            IF (DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_BACKS = 0.50_EB*ROS_BACKU(I,J)
            IF (DEGREES_SLOPE >= 20._EB)                              ROS_BACKS = 0.75_EB*ROS_BACKU(I,J)
            MAG_SR = ROS_FLANK(I,J)*(1._EB - ABS(COS_THETA_WIND)**NEXP_WIND*COS_THETA_SLOPE) + &
                     (ROS_FLANK(I,J) - ROS_BACKU(I,J))*(-ABS(COS_THETA_WIND)**NEXP_WIND) + &
                     (ROS_FLANK(I,J) - ROS_BACKS)*COS_THETA_SLOPE  !magnitude of spread rate
         ENDIF

         SR_X_LS(I,J) = MAG_SR*NORMAL_FIRELINE(1) !spread rate components
         SR_Y_LS(I,J) = MAG_SR*NORMAL_FIRELINE(2)

      ENDIF !Ellipse or McArthur Spread

   ENDDO

ENDDO FLUX_ILOOP

END SUBROUTINE LEVEL_SET_SPREAD_RATE


SUBROUTINE LEVEL_SET_ADVECT_FLUX

! Use the spread rate [SR_X_LS,SR_Y_LS] to compute the limited scalar gradient
! and take dot product with spread rate vector to get advective flux

INTEGER :: I,IM1,IP1,IP2,J,JM1,JP1,JP2
REAL(EB), DIMENSION(:) :: Z(4)
REAL(EB), POINTER, DIMENSION(:,:) :: FLUX_LS_P,F_X,F_Y
REAL(EB) :: DPHIDX,DPHIDY,SR_X_AVG,SR_Y_AVG

F_X => LS_WORK1
F_Y => LS_WORK2

IF (PREDICTOR) THEN
   PHI_LS_P => PHI_LS
   FLUX_LS_P => FLUX0_LS
ELSE
   PHI_LS_P => PHI1_LS
   FLUX_LS_P => FLUX1_LS
ENDIF

DO J=1,JBAR
   DO I=0,IBAR
      IM1 = I-1 ; IF (IM1<0) IM1 = I
      IP1 = I+1
      IP2 = I+2 ; IF (IP2>IBP1) IP2 = IP1
      Z(1) = PHI_LS_P(IM1,J)
      Z(2) = PHI_LS_P(I,J)
      Z(3) = PHI_LS_P(IP1,J)
      Z(4) = PHI_LS_P(IP2,J)
      SR_X_AVG = 0.5_EB*(SR_X_LS(MIN(IP1,IBAR),J)+SR_X_LS(MAX(1,I),J))
      F_X(I,J) = SCALAR_FACE_VALUE_LS(SR_X_AVG,Z,LIMITER_LS)
   ENDDO
ENDDO

DO J=0,JBAR
   DO I=1,IBAR
      JM1 = J-1 ; IF (JM1<0) JM1 = J
      JP1 = J+1
      JP2 = J+2 ; IF (JP2>JBP1) JP2 = JP1
      Z(1) = PHI_LS_P(I,JM1)
      Z(2) = PHI_LS_P(I,J)
      Z(3) = PHI_LS_P(I,JP1)
      Z(4) = PHI_LS_P(I,JP2)
      SR_Y_AVG = 0.5_EB*(SR_Y_LS(I,MIN(JP1,JBAR))+SR_Y_LS(I,MAX(1,J)))
      F_Y(I,J) = SCALAR_FACE_VALUE_LS(SR_Y_AVG,Z,LIMITER_LS)
   ENDDO
ENDDO

DO J=1,JBAR
   DO I=1,IBAR
      DPHIDX = (F_X(I,J)-F_X(I-1,J))*RDX(I)
      DPHIDY = (F_Y(I,J)-F_Y(I,J-1))*RDY(J)
      FLUX_LS_P(I,J) = SR_X_LS(I,J)*DPHIDX + SR_Y_LS(I,J)*DPHIDY
   ENDDO
ENDDO

END SUBROUTINE LEVEL_SET_ADVECT_FLUX


REAL(EB) FUNCTION SCALAR_FACE_VALUE_LS(SR_XY,Z,LIMITER)

! This function computes the scalar value on a face.
! The scalar is denoted Z, and the velocity is denoted U.
! The gradient (computed elsewhere) is a central difference across
! the face subject to a flux limiter.  The flux limiter choices are:
!
! LIMITER = 1 implements the MINMOD limiter
! LIMITER = 2 implements the SUPERBEE limiter of Roe
! LIMITER = 3 implements first-order upwinding (monotone)
!
!                    location of face
!
!                            f
!    |     o     |     o     |     o     |     o     |
!                     SRXY        SRXY
!                 (if f_east)  (if f_west)
!         Z(1)        Z(2)        Z(3)        Z(4)
!
INTEGER, INTENT(IN) :: LIMITER
REAL(EB) :: SR_XY
REAL(EB), INTENT(IN), DIMENSION(4) :: Z
REAL(EB) :: B,DZLOC,DZUP,R,ZUP,ZDWN

IF (SR_XY > 0._EB) THEN  ! flow is left to right

   DZLOC = Z(3)-Z(2)
   DZUP  = Z(2)-Z(1)
   IF (ABS(DZLOC) > 0._EB) THEN
      R = DZUP/DZLOC
   ELSE
      R = 0._EB
   ENDIF
   ZUP  = Z(2)
   ZDWN = Z(3)

ELSE  ! flow is right to left

   DZLOC = Z(3)-Z(2)
   DZUP  = Z(4)-Z(3)
   IF (ABS(DZLOC) > 0._EB) THEN
      R = DZUP/DZLOC
   ELSE
      R = 0._EB
   ENDIF
   ZUP  = Z(3)
   ZDWN = Z(2)
ENDIF

IF (LIMITER==1) THEN
   B = MAX(0._EB,MIN(1._EB,R))
ELSEIF (LIMITER==2) THEN
   B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
ELSEIF (LIMITER==3) THEN
   B = 0._EB
ENDIF

SCALAR_FACE_VALUE_LS = ZUP + 0.5_EB * B * ( ZDWN - ZUP )

END FUNCTION SCALAR_FACE_VALUE_LS

SUBROUTINE RAISED_VEG_MASS_ENERGY_TRANSFER(T,NM)
    
! Mass and energy transfer between gas and raised vegetation fuel elements 
!
USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION,GET_SPECIFIC_HEAT,GET_CONDUCTIVITY,GET_VISCOSITY
USE MATH_FUNCTIONS, ONLY : AFILL2
USE TRAN, ONLY: GET_IJK
!arrays for debugging
!REAL(EB), POINTER, DIMENSION(:,:,:) :: HOLD1,HOLD2,HOLD3,HOLD4
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW !,RHOP

REAL(EB), INTENT(IN) :: T
REAL(EB) :: RE_D,RCP_GAS,CP_GAS
REAL(EB) :: DT_FE,RDT_FE,V_CELL,V_VEG
REAL(EB) :: CP_ASH,CP_H2O,CP_CHAR,H_VAP_H2O,TMP_H2O_BOIL
REAL(EB) :: K_GAS,MU_GAS,RHO_GAS,RRHO_GAS_NEW,TMP_FILM,TMP_GAS,UBAR,VBAR,WBAR,UREL,VREL,WREL
REAL(EB) :: CHAR_FCTR,CHAR_FCTR2,CP_VEG,DTMP_VEG,MPV_MOIST,MPV_MOIST_MIN,DMPV_VEG,MPV_VEG,MPV_VEG_MIN, &
            SV_VEG,TMP_VEG,TMP_VEG_NEW
!REAL(EB) :: TMP_IGNITOR
REAL(EB) :: MPV_ADDED,MPV_MOIST_LOSS,MPV_VOLIT,MPV_CHAR_LOSS_MAX,MPV_MOIST_LOSS_MAX,MPV_VOLIT_MAX
REAL(EB) :: QCON_VEG,QNET_VEG,QRAD_VEG,QREL,TMP_GMV,Q_FOR_DRYING,Q_VOLIT,Q_FOR_VOLIT, &
            Q_UPTO_VOLIT
REAL(EB) :: H_SENS_VEG_VOLIT,Q_ENTHALPY,Q_VEG_MOIST,Q_VEG_VOLIT,Q_VEG_CHAR
REAL(EB) :: MW_AVERAGE,MW_VEG_MOIST_TERM,MW_VEG_VOLIT_TERM
REAL(EB) :: XI,YJ,ZK
REAL(EB) :: A_H2O_VEG,E_H2O_VEG,A_PYR_VEG,E_PYR_VEG,H_PYR_VEG,R_H_PYR_VEG
REAL(EB) :: A_CHAR_VEG,E_CHAR_VEG,BETA_CHAR_VEG,NU_CHAR_VEG,NU_ASH_VEG,NU_O2_CHAR_VEG, &
            MPV_ASH,MPV_ASH_MAX,MPV_CHAR,MPV_CHAR_LOSS,MPV_CHAR_MIN,MPV_CHAR_CO2,MPV_CHAR_O2,Y_O2, &
            H_CHAR_VEG ,ORIG_PACKING_RATIO,CP_VEG_FUEL_AND_CHAR_MASS,CP_MASS_VEG_SOLID,     &
            TMP_CHAR_MAX
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES)
INTEGER :: I,II,JJ,KK,IIX,JJY,KKZ,IPC,I_FUEL
INTEGER, INTENT(IN) :: NM
LOGICAL :: VEG_DEGRADATION_LINEAR,VEG_DEGRADATION_ARRHENIUS
!INTEGER :: IDT
REAL(EB) :: Q_VEG_CHAR_TOTAL,MPV_CHAR_CO2_TOTAL,MPV_CHAR_O2_TOTAL,MPV_CHAR_LOSS_TOTAL, &
            MPV_MOIST_LOSS_TOTAL,MPV_VOLIT_TOTAL,VEG_VF
REAL(EB) :: VEG_CRITICAL_MASSFLUX,VEG_CRITICAL_MASSSOURCE
REAL(EB) :: CM,CN
REAL(EB) :: HCON_VEG_FORCED,HCON_VEG_FREE,LENGTH_SCALE,NUSS_HILPERT_CYL_FORCEDCONV,NUSS_MORGAN_CYL_FREECONV,RAYLEIGH_NUM, &
            R_VEG_CYL_DIAM,HC_VERT_CYL,HC_HORI_CYL

!place holder
REAL(EB) :: RCP_TEMPORARY

!Debug
REAL(EB)TOTAL_BULKDENS_MOIST,TOTAL_BULKDENS_DRY_FUEL,TOTAL_MASS_DRY_FUEL,TOTAL_MASS_MOIST


!IF (.NOT. TREE_MESH(NM)) RETURN !Exit if raised veg is not present in mesh
CALL POINT_TO_MESH(NM)

UU => U
VV => V
WW => W

! Initializations
DT_FE  =  T - VEG_CLOCK_FE
RDT_FE = 1._EB/DT_FE
RCP_TEMPORARY = 1._EB/1010._EB

!Critical mass flux (kg/(s m^2)
VEG_CRITICAL_MASSFLUX = 0.0025_EB !kg/s/m^2 for qradinc=50 kW/m^2, M=4% measured by McAllister Fire Safety J., 61:200-206 2013
!VEG_CRITICAL_MASSFLUX = 0.0035_EB !kg/s/m^2 largest measured by McAllister Fire Safety J., 61:200-206 2013
!VEG_CRITICAL_MASSFLUX = 999999._EB !kg/s/m^2 for testing

!Constants for Arrhenius pyrolyis and Arrhenius char oxidation models
!are from the literature (Porterie et al., Num. Heat Transfer, 47:571-591, 2005)
CP_H2O       = 4190._EB !J/kg/K specific heat of water
TMP_H2O_BOIL = 373.15_EB
TMP_CHAR_MAX = 1300._EB !K

!Kinetic constants used by multiple investigators from Porterie or Morvan papers
!VEG_A_H2O      = 600000._EB !1/s sqrt(K)
!VEG_E_H2O      = 5800._EB !K
!VEG_A_PYR      = 36300._EB !1/s
!VEG_E_PYR      = 7250._EB !K
!VEG_E_CHAR     = 9000._EB !K
!VEG_BETA_CHAR  = 0.2_EB
!VEG_NU_O2_CHAR = 1.65_EB

!CP_ASH         = 800._EB !J/kg/K

!Kinetic constants used by Morvan and Porterie mostly obtained from Grishin
!VEG_H_PYR      = 418000._EB !J/kg 
!VEG_A_CHAR     = 430._EB !m/s 
!VEG_H_CHAR     = -12.0E+6_EB ! J/kg

!Kinetic constants used by Yolanda and Paul
!VEG_H_PYR      = 418000._EB !J/kg 
!VEG_A_CHAR     = 215._EB !m/s Yolanda, adjusted from Morvan, Porterie values based on HRR exp
!VEG_H_CHAR     = -32.74E+6_EB !J/kg via Susott

!Kinetic constants used by Shankar
!VEG_H_PYR      = 418._EB !J/kg Shankar
!VEG_A_CHAR     = 430._EB !m/s Porterie, Morvan
!VEG_H_CHAR     = -32.74E+6_EB !J/kg Shankar via Susott

!Kinetic constants used by me for ROS vs Slope excelsior experiments
!VEG_H_PYR      = 711000._EB !J/kg excelsior Catchpole et al. (via Susott)
!VEG_A_CHAR     = 430._EB !m/s Porterie, Morvan
!VEG_H_CHAR     = -32.74E+6_EB !J/kg via Susott

!R_H_PYR_VEG    = 1._EB/H_PYR_VEG

!D_AIR  = 2.6E-5_EB  ! Water Vapor - Air binary diffusion (m2/s at 25 C, Incropera & DeWitt, Table A.8) 
!SC_AIR = 0.6_EB     ! NU_AIR/D_AIR (Incropera & DeWitt, Chap 7, External Flow)
!PR_AIR = 0.7_EB     

! Working arrays
!IF(N_TREES_OUT > 0) TREE_OUTPUT_DATA(:,:,NM) = 0._EB !for output of veg data
!DMPVDT_FM_VEG  = 0.0_EB

!Clear arrays and scalars
!HOLD1 => WORK4 ; WORK4 = 0._EB
!HOLD2 => WORK5 ; WORK5 = 0._EB
!HOLD3 => WORK6 ; WORK6 = 0._EB
!HOLD4 => WORK7 ; WORK7 = 0._EB
TOTAL_BULKDENS_MOIST    = 0.0_EB
TOTAL_BULKDENS_DRY_FUEL = 0.0_EB
TOTAL_MASS_MOIST    = 0.0_EB
TOTAL_MASS_DRY_FUEL = 0.0_EB
V_VEG               = 0.0_EB

!print*,'vege h-m transfer: NM, NLP',nm,nlp

PARTICLE_LOOP: DO I=1,NLP

 LP  => LAGRANGIAN_PARTICLE(I)
 IPC =  LP%CLASS_INDEX
 LPC => LAGRANGIAN_PARTICLE_CLASS(IPC)
 IF (.NOT. LPC%VEG_WFDS_FE) CYCLE PARTICLE_LOOP !Ensure WFDS FE vegetation model is to be used
!IF (LPC%MASSLESS) CYCLE PARTICLE_LOOP   !Skip PARTICLE if massless

 THERMAL_CALC: IF (.NOT. LPC%VEG_STEM) THEN   !compute heat transfer, etc if thermally thin

! Intialize quantities
 LP%VEG_MLR      = 0.0_EB
 LP%VEG_Q_CHAROX = 0.0_EB
 Q_VEG_CHAR      = 0.0_EB
 Q_VEG_MOIST     = 0.0_EB
 Q_VEG_VOLIT     = 0.0_EB
 Q_UPTO_VOLIT    = 0.0_EB
 Q_VOLIT         = 0.0_EB
 MPV_MOIST_LOSS  = 0.0_EB
 MPV_CHAR_LOSS   = 0.0_EB
 MPV_CHAR_CO2    = 0.0_EB
 MPV_CHAR_O2     = 0.0_EB
 MPV_VOLIT       = 0.0_EB
 MPV_ADDED       = 0.0_EB
 MW_VEG_MOIST_TERM = 0.0_EB
 MW_VEG_VOLIT_TERM = 0.0_EB
 CP_VEG_FUEL_AND_CHAR_MASS = 0.0_EB
 CP_MASS_VEG_SOLID         = 0.0_EB
 VEG_DEGRADATION_LINEAR    = .FALSE.
 VEG_DEGRADATION_ARRHENIUS = .FALSE.
 MPV_CHAR_CO2_TOTAL   = 0.0_EB
 MPV_CHAR_O2_TOTAL   = 0.0_EB 
 MPV_CHAR_LOSS_TOTAL  = 0.0_EB 
 MPV_MOIST_LOSS_TOTAL = 0.0_EB 
 MPV_VOLIT_TOTAL  = 0.0_EB 
 Q_VEG_CHAR_TOTAL = 0.0_EB

! Vegetation variables
 VEG_VF             = LPC%VEG_VOLUME_FRACTION !volume fraction of vegetation in cell
 NU_CHAR_VEG        = LPC%VEG_CHAR_FRACTION
 NU_ASH_VEG         = LPC%VEG_ASH_FRACTION/LPC%VEG_CHAR_FRACTION !fraction of char that can become ash
 CHAR_FCTR          = 1._EB - LPC%VEG_CHAR_FRACTION !factor used to determine volatile mass
 CHAR_FCTR2         = 1._EB/CHAR_FCTR !factor used to determine char mass
 SV_VEG             = LPC%VEG_SV !surface-to-volume ration 1/m
 TMP_VEG            = LP%VEG_TMP
 MPV_VEG            = LP%VEG_FUEL_MASS !bulk density of dry veg kg/m^3
 MPV_CHAR           = LP%VEG_CHAR_MASS !bulk density of char
 MPV_ASH            = LP%VEG_ASH_MASS  !bulk density of ash 
 MPV_MOIST          = LP%VEG_MOIST_MASS !bulk density of moisture in veg
 MPV_VEG_MIN        = LPC%VEG_FUEL_MPV_MIN 
 MPV_CHAR_MIN       = MPV_VEG_MIN*LPC%VEG_CHAR_FRACTION
 MPV_MOIST_MIN      = LPC%VEG_MOIST_MPV_MIN
 MPV_ASH_MAX        = LPC%VEG_ASH_MPV_MAX   !maxium ash bulk density
 MPV_MOIST_LOSS_MAX = LPC%VEG_DEHYDRATION_RATE_MAX*DT_FE
 MPV_VOLIT_MAX      = LPC%VEG_BURNING_RATE_MAX*DT_FE
 MPV_CHAR_LOSS_MAX  = LPC%VEG_CHAROX_RATE_MAX*DT_FE
 ORIG_PACKING_RATIO = LPC%VEG_BULK_DENSITY/LPC%VEG_DENSITY 
 H_VAP_H2O          = LPC%VEG_H_H2O !J/kg/K heat of vaporization of water
 A_H2O_VEG          = LPC%VEG_A_H2O !1/s sqrt(K)
 E_H2O_VEG          = LPC%VEG_E_H2O !K
 H_PYR_VEG          = LPC%VEG_H_PYR !J/kg 
 A_PYR_VEG          = LPC%VEG_A_PYR !1/s
 E_PYR_VEG          = LPC%VEG_E_PYR !K
 H_CHAR_VEG         = LPC%VEG_H_CHAR ! J/kg
 A_CHAR_VEG         = LPC%VEG_A_CHAR !m/s 
 E_CHAR_VEG         = LPC%VEG_E_CHAR !K
 BETA_CHAR_VEG      = LPC%VEG_BETA_CHAR
 NU_O2_CHAR_VEG     = LPC%VEG_NU_O2_CHAR

! Thermal degradation approach parameters
 IF(LPC%VEG_DEGRADATION == 'LINEAR') VEG_DEGRADATION_LINEAR = .TRUE.
 IF(LPC%VEG_DEGRADATION == 'ARRHENIUS') VEG_DEGRADATION_ARRHENIUS = .TRUE.

 R_H_PYR_VEG    = 1._EB/H_PYR_VEG

!Bound on volumetric mass flux
 VEG_CRITICAL_MASSSOURCE = VEG_CRITICAL_MASSFLUX*SV_VEG*LP%VEG_PACKING_RATIO

! Determine grid cell quantities of the vegetation fuel element
 CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,II,JJ,KK)
 IIX = FLOOR(XI+0.5_EB)
 JJY = FLOOR(YJ+0.5_EB)
 KKZ = FLOOR(ZK+0.5_EB)
 V_CELL = DX(II)*DY(JJ)*DZ(KK)

! Gas velocities in vegetation grid cell
 UBAR = AFILL2(UU,II-1,JJY,KKZ,XI-II+1,YJ-JJY+.5_EB,ZK-KKZ+.5_EB)
 VBAR = AFILL2(VV,IIX,JJ-1,KKZ,XI-IIX+.5_EB,YJ-JJ+1,ZK-KKZ+.5_EB)
 WBAR = AFILL2(WW,IIX,JJY,KK-1,XI-IIX+.5_EB,YJ-JJY+.5_EB,ZK-KK+1)
 UREL = LP%U - UBAR
 VREL = LP%V - VBAR
 WREL = LP%W - WBAR
 QREL = MAX(1.E-6_EB,SQRT(UREL*UREL + VREL*VREL + WREL*WREL))


! Gas thermophysical quantities
 RHO_GAS  = RHO(II,JJ,KK)
 TMP_GAS  = TMP(II,JJ,KK)
 TMP_FILM = 0.5_EB*(TMP_GAS + TMP_VEG)
!print '(A,6ES12.3)','vege:tmp_veg,mpv_veg,mpv_char,mpv_ash,mpv_moist,tmp_gas',tmp_veg,mpv_veg,mpv_char,mpv_ash,mpv_moist,tmp_gas

! Assuming gas is air
!RHO_AIR  = 101325./(287.05*TMP_FILM) !rho_air = standard pressure / (ideal gas constant*gas temp)
!MU_AIR   =  (0.000001458_EB*TMP_FILM**1.5_EB)/(TMP_FILM+110.4_EB) !kg/m/s
!K_AIR    = (0.002495_EB*TMP_FILM**1.5_EB)/(TMP_FILM+194._EB) !W/m.K

!Use full gas composition
 ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,ZZ(II,JJ,KK,1:N_TRACKED_SPECIES))
 CALL GET_VISCOSITY(ZZ_GET,MU_GAS,TMP_FILM) 
 CALL GET_CONDUCTIVITY(ZZ_GET,K_GAS,TMP_FILM) !W/m/K
 CALL GET_SPECIFIC_HEAT(ZZ_GET,CP_GAS,TMP_FILM)

! Veg thermophysical properties
 TMP_GMV  = TMP_GAS - TMP_VEG
 CP_VEG   = (0.01_EB + 0.0037_EB*TMP_VEG)*1000._EB !J/kg/K Ritchie IAFSS 1997:177-188
 CP_CHAR  = 420._EB + 2.09_EB*TMP_VEG + 6.85E-4_EB*TMP_VEG**2 !J/kg/K Park etal. C&F 2010 147:481-494
 CP_ASH   = 1244._EB*(TMP_VEG/TMPA)**0.315 !J/kg/K Lautenberger & Fernandez-Pell, C&F 2009 156:1503-1513
 R_VEG_CYL_DIAM = 0.25_EB*SV_VEG

! Convective heat flux on thermal elements

 IF_QCONV_MAX: IF (LPC%VEG_HCONV_CYLMAX) THEN !Max of forced and free

   RE_D     = RHO_GAS*QREL*4._EB/(SV_VEG*MU_GAS)

! - Forced convection heat transfer coefficients on veg particles
!
! Hilpert Correlation (Incropera & DeWitt Fourth Edition, p. 370) for cylinder in crossflow,
! forced convection
   IF(RE_D < 4._EB) THEN
     CN = 0.989_EB
     CM = 0.330_EB
   ELSE IF (RE_D >= 4._EB .AND. RE_D < 40._EB) THEN
     CN = 0.911_EB
     CM = 0.385_EB
   ELSE
     CN = 0.683_EB
     CM = 0.466_EB
   ENDIF
   NUSS_HILPERT_CYL_FORCEDCONV = CN*(RE_D**CM)*PR_ONTH !Nusselt number
!print '(A,2x,2ES12.4)','nuss Hilpert,Re', NUSS_HILPERT_CYL_FORCEDCONV,re_d
   HCON_VEG_FORCED = 0.25_EB*SV_VEG*K_GAS*NUSS_HILPERT_CYL_FORCEDCONV !W/m^2 from Hilpert (cylinder)

! - Free convection heat transfer coefficients
   LENGTH_SCALE = 4._EB/SV_VEG !horizontal cylinder diameter
   RAYLEIGH_NUM = 9.8_EB*ABS(TMP_GMV)*LENGTH_SCALE**3*RHO_GAS**2*CP_GAS/(TMP_FILM*MU_GAS*K_GAS)
!print*,'ZZ_GET',ZZ_GET(:)

! Morgan correlation free convection (Incropera & DeWitt, 4th Edition, p. 501-502) for horizontal cylinder of diameter
! 4/SV_VEG, free convection
   IF (RAYLEIGH_NUM < 0.01_EB) THEN
     CN = 0.675_EB
     CM = 0.058_EB
   ELSE IF (RAYLEIGH_NUM >= 0.01_EB .AND. RAYLEIGH_NUM < 100._EB) THEN
     CN = 1.02_EB
     CM = 0.148_EB
   ELSE IF (RAYLEIGH_NUM >= 100._EB .AND. RAYLEIGH_NUM < 10**4._EB) THEN
     CN = 0.85_EB
     CM = 0.188_EB
   ELSE IF (RAYLEIGH_NUM >= 10**4._EB .AND. RAYLEIGH_NUM < 10**7._EB) THEN
     CN = 0.48_EB
     CM = 0.25_EB
   ELSE IF (RAYLEIGH_NUM >= 10**7._EB .AND. RAYLEIGH_NUM < 10**12._EB) THEN
     CN = 0.125_EB
     CM = 0.333_EB
   ENDIF
   NUSS_MORGAN_CYL_FREECONV = CN*RAYLEIGH_NUM**CM
   HCON_VEG_FREE = 0.25_EB*SV_VEG*K_GAS*NUSS_MORGAN_CYL_FREECONV !W/m^2

! Holman correlation for free convection (used in FDS) 
! In single particle test gave same result as Morgan free convection above
!  HCON_VEG_FREE = 1.31_EB*ABS(TMP_GMV)**ONTH

   QCON_VEG = MAX(HCON_VEG_FORCED,HCON_VEG_FREE)*TMP_GMV !W/m^2

!Weighting of Nusselt numbers following Morvan et al. "A 3D physical model ...", 101:39-52 2018
   QCON_VEG = 0.25_EB*SV_VEG*K_GAS*SQRT(NUSS_HILPERT_CYL_FORCEDCONV**2._EB + NUSS_MORGAN_CYL_FREECONV**2._EB)

 ENDIF IF_QCONV_MAX

 IF (LPC%VEG_HCONV_CYLLAM) THEN  !Laminar flow, cyl of diameter 4/sv_veg
  HC_VERT_CYL = 1.42_EB*(ABS(TMP_GMV)*R_VEG_CYL_DIAM)**0.25_EB !Holman vertical cylinder
  HC_HORI_CYL = 1.32_EB*(ABS(TMP_GMV)*R_VEG_CYL_DIAM)**0.25_EB !Holman horizontal cylinder
  QCON_VEG    = TMP_GMV*0.5_EB*(HC_VERT_CYL+HC_HORI_CYL) !average of vertical and horizontal cylinder 
 ENDIF 

 
! IF (TMP_VEG >= TMP_GAS )QCON_VEG = SV_VEG*(0.5_EB*K_AIR*0.683_EB*RE_D**0.466_EB)*0.5_EB*TMP_GMV !W/m^2 from Porterie
! IF (TMP_VEG <  TMP_GAS ) QCON_VEG = TMP_GMV*1.42_EB*(ABS(TMP_GMV)/DZ(KK))**0.25_EB !Holman
!RE_D     = RHO_GAS*QREL*2._EB/(SV_VEG*MU_GAS)
!QCON_VEG = SV_VEG*(0.5_EB*K_GAS*0.683_EB*RE_D**0.466_EB)*0.5_EB*TMP_GMV !W/m^2 from Porterie (cylinder)
! QCON_VEG = TMP_GMV*1.42_EB*(ABS(TMP_GMV)/DZ(KK))**0.25_EB !Holman vertical cylinders of length dz, Laminar flow
! QCON_VEG = TMP_GMV*1.32_EB*(ABS(TMP_GMV)*SV_VEG*0.25_EB)**0.25_EB !Holman horizontal cylinders of diameter 4/sv_veg, Laminar flow

 QCON_VEG = SV_VEG*LP%VEG_PACKING_RATIO*QCON_VEG !W/m^3
 LP%VEG_DIVQC = QCON_VEG
 QRAD_VEG     = LP%VEG_DIVQR
!print '(A,1x,5E13.5)','vege:tmp_veg,tmp_gas,tmp_gmv,qcon,qrad',tmp_veg-273,tmp_gas-273,tmp_gmv,qcon_veg,qrad_veg
! Divergence of net heat flux
 QNET_VEG = QCON_VEG + QRAD_VEG !W/m^3

! Update temperature of vegetation
!CP_VEG_FUEL_AND_CHAR_MASS = CP_VEG*MPV_VEG + CP_CHAR*MPV_CHAR
!DTMP_VEG    = DT_FE*QNET_VEG/(CP_VEG_FUEL_AND_CHAR_MASS + CP_H2O*MPV_MOIST)
 CP_MASS_VEG_SOLID = CP_VEG*MPV_VEG + CP_CHAR*MPV_CHAR + CP_ASH*MPV_ASH
!print '(A,4ES12.3)','mpv_veg,mpv_char,mpv_ash,mpv_moist', mpv_veg,mpv_char,mpv_ash,mpv_moist
!print '(A,5ES12.3)','vege:dt,tmpveg,qnet_veg,cp_mass_veg_solid,cph20,mpv_moist',tmp_veg,qnet_veg,cp_mass_veg_solid,cp_h2o,mpv_moist
 DTMP_VEG    = DT_FE*QNET_VEG/(CP_MASS_VEG_SOLID + CP_H2O*MPV_MOIST)
 TMP_VEG_NEW = TMP_VEG + DTMP_VEG
 IF (TMP_VEG_NEW < TMPA) TMP_VEG_NEW = TMP_GAS
!print*,'---------------------------------------------------------'
!print 1113,ii,jj,kk,idt
!1113 format(2x,4(I3))
!print 1112,tmp_veg_new,tmp_veg,qnet_veg,cp_mass_veg_solid,cp_h2o,dtmp_veg
!1112 format(2x,6(e15.5))

! Set temperature of inert ignitor elements
!IF(LP%IGNITOR) THEN
! TMP_IGNITOR = LPC%VEG_INITIAL_TEMPERATURE
! TMP_VEG_NEW = TMP_GAS
! IF(T>=LP%VEG_IGN_TON .AND. T<=LP%VEG_IGN_TON+LP%VEG_IGN_TRAMPON) THEN
!   TMP_VEG_NEW = &
!     TMPA + (TMP_IGNITOR-TMPA)*(T-LP%VEG_IGN_TON)/LP%VEG_IGN_TRAMPON
! ENDIF  
! IF(T>LP%VEG_IGN_TON+LP%VEG_IGN_TRAMPON) TMP_VEG_NEW = TMP_IGNITOR
! IF(T>=LP%VEG_IGN_TOFF .AND. T<=LP%VEG_IGN_TOFF+LP%VEG_IGN_TRAMPOFF)THEN 
!   TMP_VEG_NEW = &
!     TMP_IGNITOR - (TMP_IGNITOR-TMP_GAS)*(T-LP%VEG_IGN_TOFF)/LP%VEG_IGN_TRAMPOFF
! ENDIF
! IF(T > LP%VEG_IGN_TOFF+LP%VEG_IGN_TRAMPOFF) THEN
!  LP%R = 0.0001_EB*LPC%KILL_RADIUS !remove ignitor element
!  TMP_VEG_NEW = TMP_GAS
! ENDIF
!ENDIF

!      ************** Fuel Element Linear Pyrolysis Degradation model *************************
! Drying occurs if qnet > 0 with Tveg held at 100 c
! Pyrolysis occurs if qnet > 0 according to Morvan & Dupuy empirical formula. Linear
! temperature dependence with qnet factor. 
! Char oxidation occurs if qnet > 0 (user must request char ox) after pyrolysis is completed.
!
 IF_VEG_DEGRADATION_LINEAR: IF(VEG_DEGRADATION_LINEAR) THEN
!  IF_NET_HEAT_INFLUX: IF (QNET_VEG > 0.0_EB .AND. .NOT. LP%IGNITOR) THEN !dehydrate or pyrolyze 
   IF_NET_HEAT_INFLUX: IF (QNET_VEG > 0.0_EB) THEN !dehydrate or pyrolyze 

! Drying of fuel element vegetation 
     IF_DEHYDRATION: IF (MPV_MOIST > MPV_MOIST_MIN .AND. TMP_VEG_NEW > TMP_H2O_BOIL) THEN
       Q_FOR_DRYING   = (TMP_VEG_NEW - TMP_H2O_BOIL)/DTMP_VEG * QNET_VEG
       MPV_MOIST_LOSS = MIN(DT_FE*Q_FOR_DRYING/H_VAP_H2O,MPV_MOIST-MPV_MOIST_MIN)
       MPV_MOIST_LOSS = MIN(MPV_MOIST_LOSS,MPV_MOIST_LOSS_MAX) !use specified max
       TMP_VEG_NEW       = TMP_H2O_BOIL
       LP%VEG_MOIST_MASS = MPV_MOIST - MPV_MOIST_LOSS !kg/m^3
       IF (LP%VEG_MOIST_MASS <= MPV_MOIST_MIN) LP%VEG_MOIST_MASS = 0.0_EB
       Q_VEG_MOIST       = MPV_MOIST_LOSS*CP_H2O*(TMP_VEG_NEW - TMPA)
       MW_VEG_MOIST_TERM = MPV_MOIST_LOSS/MW_H2O
!      IF (I == 1) print*,MPV_MOIST,MPV_MOIST_LOSS
!Print '(A,1x,3ES13.3)','vege:tmp_veg_new,mpv_moist_loss,mass h2o',tmp_veg_new-273,mpv_moist_loss,lp%veg_moist_mass
     ENDIF IF_DEHYDRATION

! Volitalization of fuel element vegetation
     IF_VOLITALIZATION: IF(MPV_MOIST <= MPV_MOIST_MIN) THEN

       IF_MD_VOLIT: IF(MPV_VEG > MPV_VEG_MIN .AND. TMP_VEG_NEW >= 400._EB) THEN !Morvan & Dupuy volitalization
         Q_UPTO_VOLIT = CP_MASS_VEG_SOLID*MAX((400._EB-TMP_VEG),0._EB)
         Q_FOR_VOLIT  = DT_FE*QNET_VEG - Q_UPTO_VOLIT
         Q_VOLIT      = Q_FOR_VOLIT*0.01_EB*(MIN(500._EB,TMP_VEG)-400._EB)

!        MPV_VOLIT    = Q_VOLIT*R_H_PYR_VEG
         MPV_VOLIT    = CHAR_FCTR*Q_VOLIT*R_H_PYR_VEG
         MPV_VOLIT    = MAX(MPV_VOLIT,0._EB)
         MPV_VOLIT    = MIN(MPV_VOLIT,MPV_VOLIT_MAX) !user specified max

         DMPV_VEG     = CHAR_FCTR2*MPV_VOLIT
         DMPV_VEG     = MIN(DMPV_VEG,(MPV_VEG - MPV_VEG_MIN))
         MPV_VEG      = MPV_VEG - DMPV_VEG
!Print '(A,1x,1ES13.3)','vege: dmpv_veg',dmpv_veg

         MPV_VOLIT    = CHAR_FCTR*DMPV_VEG
!        MPV_CHAR     = MPV_CHAR + NU_CHAR_VEG*MPV_VOLIT !kg/m^3
!        MPV_CHAR     = MPV_CHAR + LPC%VEG_CHAR_FRACTION*MPV_VOLIT !kg/m^3
         MPV_CHAR     = MPV_CHAR + LPC%VEG_CHAR_FRACTION*DMPV_VEG !kg/m^3
         Q_VOLIT      = MPV_VOLIT*H_PYR_VEG
         CP_MASS_VEG_SOLID = CP_VEG*MPV_VEG + CP_CHAR*MPV_CHAR 
         TMP_VEG_NEW  = TMP_VEG + (Q_FOR_VOLIT-Q_VOLIT)/CP_MASS_VEG_SOLID
         TMP_VEG_NEW  = MIN(TMP_VEG_NEW,500._EB) !set to high pyrol temp if too hot

!Handle veg. fuel elements if element mass <= prescribed minimum
         IF (MPV_VEG <= MPV_VEG_MIN) THEN
           MPV_VEG = 0.0_EB
           IF(LPC%VEG_REMOVE_CHARRED .AND. .NOT. LPC%VEG_CHAR_OXIDATION) LP%ONE_D%BURNAWAY = .TRUE.
!           LP%R = 0.0001_EB*LPC%KILL_RADIUS !fuel element will be removed
         ENDIF
!Enthalpy of fuel element volatiles using Cp,volatiles(T) from Ritchie
         H_SENS_VEG_VOLIT = 0.0445_EB*(TMP_VEG**1.5_EB - TMP_GAS**1.5_EB) - 0.136_EB*(TMP_VEG - TMP_GAS)
         H_SENS_VEG_VOLIT = H_SENS_VEG_VOLIT*1000._EB !J/kg
         Q_VEG_VOLIT      = CHAR_FCTR*MPV_VOLIT*H_SENS_VEG_VOLIT !J/m^3
         MW_VEG_VOLIT_TERM= MPV_VOLIT/SPECIES(FUEL_INDEX)%MW
        ENDIF IF_MD_VOLIT

      LP%VEG_FUEL_MASS = MPV_VEG
      LP%VEG_CHAR_MASS = MPV_CHAR !kg/m^3

    ENDIF IF_VOLITALIZATION

   ENDIF IF_NET_HEAT_INFLUX

!Char oxidation of fuel element with the Linear pyrolysis model from Morvan and Dupuy, Comb.
!Flame, 138:199-210 (2004)
!(note that this can be handled only approximately with the conserved
!scalar based gas-phase combustion model - the oxygen is consumed by
!the char oxidation reaction is not accounted for since it would be inconsistent with the state
!relation for oxygen that is based on the conserved scalar approach used for gas phase
!combustion)
   IF_CHAR_OXIDATION_LIN: IF (LPC%VEG_CHAR_OXIDATION .AND. MPV_VEG <= MPV_VEG_MIN) THEN

     ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,ZZ(II,JJ,KK,1:N_TRACKED_SPECIES))
     CALL GET_MASS_FRACTION(ZZ_GET,O2_INDEX,Y_O2)
     MPV_CHAR_LOSS = DT_FE*RHO_GAS*Y_O2*A_CHAR_VEG/NU_O2_CHAR_VEG*SV_VEG*LP%VEG_PACKING_RATIO*  &
                      EXP(-E_CHAR_VEG/TMP_VEG)*(1+BETA_CHAR_VEG*SQRT(RE_D))
     MPV_CHAR_LOSS = MIN(MPV_CHAR_LOSS,MPV_CHAR_LOSS_MAX) !user bound
     MPV_CHAR_LOSS = MIN(MPV_CHAR,MPV_CHAR_LOSS)
     MPV_CHAR      = MPV_CHAR - MPV_CHAR_LOSS
     MPV_ASH       = MPV_ASH + NU_ASH_VEG*MPV_CHAR_LOSS
     MPV_CHAR_CO2  = (1._EB + NU_O2_CHAR_VEG - NU_ASH_VEG)*MPV_CHAR_LOSS
     MPV_CHAR_O2   = NU_O2_CHAR_VEG*MPV_CHAR_LOSS
     CP_MASS_VEG_SOLID = CP_VEG*MPV_VEG + CP_CHAR*MPV_CHAR + CP_ASH*MPV_ASH
     LP%VEG_CHAR_MASS  = MPV_CHAR !kg/m^3
     LP%VEG_ASH_MASS   = MPV_ASH

! Reduce fuel element size based on char consumption
!!    IF (MPV_VEG <= MPV_VEG_MIN) THEN !charring reduce veg elem size
!      LP%VEG_PACKING_RATIO = LP%VEG_PACKING_RATIO - MPV_CHAR_LOSS/(LPC%VEG_DENSITY*LPC%VEG_CHAR_FRACTION)
!      LP%VEG_SV     = LPC%VEG_SV*(ORIG_PACKING_RATIO/LP%VEG_PACKING_RATIO)**0.333_EB 
!      LP%VEG_KAPPA  = 0.25_EB*LP%VEG_SV*LP%VEG_PACKING_RATIO
!!    ENDIF

!remove fuel element if char ox is complete
      IF (MPV_CHAR <= MPV_CHAR_MIN) THEN 
        CP_MASS_VEG_SOLID = CP_ASH*MPV_ASH
        LP%VEG_CHAR_MASS = 0.0_EB
        IF(LPC%VEG_REMOVE_CHARRED) LP%ONE_D%BURNAWAY = .TRUE. !fuel element will be removed
      ENDIF

      Q_VEG_CHAR       = MPV_CHAR_LOSS*H_CHAR_VEG 
      LP%VEG_Q_CHAROX  = -Q_VEG_CHAR*RDT_FE
      Q_VEG_CHAR_TOTAL = Q_VEG_CHAR_TOTAL + Q_VEG_CHAR
      TMP_VEG_NEW  = TMP_VEG_NEW - LPC%VEG_CHAR_ENTHALPY_FRACTION*Q_VEG_CHAR/CP_MASS_VEG_SOLID
!     TMP_VEG_NEW  = MIN(TMP_CHAR_MAX,TMP_VEG_NEW)
!     print*,'vege: q_veg_char,temp_veg_new,',q_veg_char,tmp_veg_new
!          print*,'------------------'
!    ENDIF IF_CHAR_OXIDATION_LIN_2

   ENDIF IF_CHAR_OXIDATION_LIN
  
 ENDIF IF_VEG_DEGRADATION_LINEAR

!      ************** Fuel Element Arrehnius Degradation model *************************
! Drying and pyrolysis of fuel element occur according to Arrehnius expressions obtained 
! from the literature (Porterie et al., Num. Heat Transfer, 47:571-591, 2005
! Predicting wildland fire behavior and emissions using a fine-scale physical
! model
!
 IF_VEG_DEGRADATION_ARRHENIUS: IF(VEG_DEGRADATION_ARRHENIUS) THEN

!  TMP_VEG = TMPA + 5._EB/60._EB*T ; TMP_VEG_NEW = TMP_VEG !mimic TGA with 5 C/min heating rate

!  IF_NOT_IGNITOR1: IF (.NOT. LP%IGNITOR) THEN !dehydrate or pyrolyze 

! Drying of fuel element vegetation 
     IF_DEHYDRATION_2: IF (MPV_MOIST > MPV_MOIST_MIN) THEN
       MPV_MOIST_LOSS = MIN(DT_FE*MPV_MOIST*A_H2O_VEG*EXP(-E_H2O_VEG/TMP_VEG)/SQRT(TMP_VEG), &
                            MPV_MOIST-MPV_MOIST_MIN)
       MPV_MOIST_LOSS = MIN(MPV_MOIST_LOSS,MPV_MOIST_LOSS_MAX) !use specified max
       MPV_MOIST      = MPV_MOIST - MPV_MOIST_LOSS
       LP%VEG_MOIST_MASS = MPV_MOIST !kg/m^3
       IF (MPV_MOIST <= MPV_MOIST_MIN) LP%VEG_MOIST_MASS = 0.0_EB
       MW_VEG_MOIST_TERM = MPV_MOIST_LOSS/MW_H2O
       Q_VEG_MOIST  = MPV_MOIST_LOSS*CP_H2O*(TMP_VEG - TMPA)
!      IF (I == 1) print*,MPV_MOIST,MPV_MOIST_LOSS
     ENDIF IF_DEHYDRATION_2

! Volitalization of fuel element vegetation
     IF_VOLITALIZATION_2: IF(MPV_VEG > MPV_VEG_MIN) THEN
       MPV_VOLIT    = DT_FE*CHAR_FCTR*MPV_VEG*A_PYR_VEG*EXP(-E_PYR_VEG/TMP_VEG)
       MPV_VOLIT    = MIN(MPV_VOLIT,MPV_VOLIT_MAX) !user specified max

       DMPV_VEG     = CHAR_FCTR2*MPV_VOLIT
       DMPV_VEG     = MIN(DMPV_VEG,(MPV_VEG - MPV_VEG_MIN))
       MPV_VEG      = MPV_VEG - DMPV_VEG

       MPV_VOLIT    = CHAR_FCTR*DMPV_VEG 
       MPV_CHAR     = MPV_CHAR + LPC%VEG_CHAR_FRACTION*DMPV_VEG !kg/m^3
!      MPV_CHAR     = MPV_CHAR + LPC%VEG_CHAR_FRACTION*MPV_VOLIT !kg/m^3
       CP_MASS_VEG_SOLID = CP_VEG*MPV_VEG + CP_CHAR*MPV_CHAR + CP_ASH*MPV_ASH

!Handle veg. fuel elements if original element mass <= prescribed minimum
       IF (MPV_VEG <= MPV_VEG_MIN) THEN
!        MPV_VEG = MPV_VEG_MIN
         MPV_VEG = 0.0_EB
         CP_MASS_VEG_SOLID = CP_CHAR*MPV_CHAR + CP_ASH*MPV_ASH
!        IF(LPC%VEG_REMOVE_CHARRED .AND. .NOT. LPC%VEG_CHAR_OXIDATION) LP%R = 0.0001_EB*LPC%KILL_RADIUS !remove part
         IF(LPC%VEG_REMOVE_CHARRED .AND. .NOT. LPC%VEG_CHAR_OXIDATION) LP%ONE_D%BURNAWAY = .TRUE. !remove part
       ENDIF
!Enthalpy of fuel element volatiles using Cp,volatiles(T) from Ritchie
       H_SENS_VEG_VOLIT = 0.0445_EB*(TMP_VEG**1.5_EB - TMP_GAS**1.5_EB) - 0.136_EB*(TMP_VEG - TMP_GAS)
       H_SENS_VEG_VOLIT = H_SENS_VEG_VOLIT*1000._EB !J/kg
       Q_VEG_VOLIT      = MPV_VOLIT*H_SENS_VEG_VOLIT !J
       MW_VEG_VOLIT_TERM= MPV_VOLIT/SPECIES(FUEL_INDEX)%MW
     ENDIF IF_VOLITALIZATION_2

     LP%VEG_FUEL_MASS = MPV_VEG
     LP%VEG_CHAR_MASS = MPV_CHAR

!Char oxidation of fuel element within the Arrhenius pyrolysis model
!(note that this can be handled only approximately with the conserved
!scalar based gas-phase combustion model - no gas phase oxygen is consumed by
!the char oxidation reaction since it would be inconsistent with the state
!relation for oxygen based on the conserved scalar approach for gas phase
!combustion)
     IF_CHAR_OXIDATION: IF (LPC%VEG_CHAR_OXIDATION) THEN
       ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
       CALL GET_MASS_FRACTION(ZZ_GET,O2_INDEX,Y_O2)
       MPV_CHAR_LOSS = DT_FE*RHO_GAS*Y_O2*A_CHAR_VEG/NU_O2_CHAR_VEG*SV_VEG*LP%VEG_PACKING_RATIO*  &
                        EXP(-E_CHAR_VEG/TMP_VEG)*(1._EB+BETA_CHAR_VEG*SQRT(RE_D))
       MPV_CHAR_LOSS = MIN(MPV_CHAR_LOSS,MPV_CHAR_LOSS_MAX) !user bound
       MPV_CHAR_LOSS = MIN(MPV_CHAR,MPV_CHAR_LOSS)
       MPV_CHAR      = MPV_CHAR - MPV_CHAR_LOSS
       MPV_ASH       = MPV_ASH + NU_ASH_VEG*MPV_CHAR_LOSS
       MPV_CHAR_CO2  = (1._EB + NU_O2_CHAR_VEG - NU_ASH_VEG)*MPV_CHAR_LOSS
       MPV_CHAR_O2   = NU_O2_CHAR_VEG*MPV_CHAR_LOSS
       CP_MASS_VEG_SOLID = CP_VEG*MPV_VEG + CP_CHAR*MPV_CHAR + CP_ASH*MPV_ASH
       LP%VEG_CHAR_MASS = MPV_CHAR !kg/m^3
       LP%VEG_ASH_MASS  = MPV_ASH

! Reduce veg element size based on char consumption
!      LP%VEG_PACKING_RATIO = LP%VEG_PACKING_RATIO - MPV_CHAR_LOSS/(LPC%VEG_DENSITY*LPC%VEG_CHAR_FRACTION)
!      LP%VEG_SV     = LPC%VEG_SV*(ORIG_PACKING_RATIO/LP%VEG_PACKING_RATIO)**0.333_EB 
!      LP%VEG_KAPPA  = 0.25_EB*LP%VEG_SV*LP%VEG_PACKING_RATIO

! Remove partical if char is fully consumed
       IF (MPV_CHAR <= MPV_CHAR_MIN .AND. MPV_VEG <= MPV_VEG_MIN) THEN 
!        IF (MPV_ASH >= MPV_ASH_MAX .AND. MPV_VEG <= MPV_VEG_MIN) THEN 
!        CP_MASS_VEG_SOLID = CP_CHAR*MPV_CHAR_MIN
         CP_MASS_VEG_SOLID = CP_ASH*MPV_ASH
         LP%VEG_CHAR_MASS = 0.0_EB
!        IF(LPC%VEG_REMOVE_CHARRED) LP%R = 0.0001_EB*LPC%KILL_RADIUS !fuel element will be removed
         IF(LPC%VEG_REMOVE_CHARRED) LP%ONE_D%BURNAWAY = .TRUE. !fuel element will be removed
       ENDIF
!  ENDIF IF_CHAR_OXIDATION_2

     ENDIF IF_CHAR_OXIDATION

     Q_VEG_CHAR        = MPV_CHAR_LOSS*H_CHAR_VEG 
     LP%VEG_Q_CHAROX   = -Q_VEG_CHAR*RDT_FE
     Q_VEG_CHAR_TOTAL  =  Q_VEG_CHAR_TOTAL + Q_VEG_CHAR
     TMP_VEG_NEW  = TMP_VEG_NEW - (MPV_MOIST_LOSS*H_VAP_H2O + MPV_VOLIT*H_PYR_VEG + & 
                                  LPC%VEG_CHAR_ENTHALPY_FRACTION*Q_VEG_CHAR) / &
                                 (LP%VEG_MOIST_MASS*CP_H2O + CP_MASS_VEG_SOLID)
     TMP_VEG_NEW  = MIN(TMP_CHAR_MAX,TMP_VEG_NEW)
!print 1111,tmp_veg_new,mpv_moist,mpv_volit,q_veg_char,mpv_char_loss
!1111 format(2x,5(e15.5))
!    IF (MPV_VEG <= MPV_VEG_MIN) MPV_VOLIT = 0.0_EB

!  ENDIF IF_NOT_IGNITOR1
 ENDIF IF_VEG_DEGRADATION_ARRHENIUS

 LP%ONE_D%TMP_F = TMP_VEG_NEW
 LP%VEG_TMP     = TMP_VEG_NEW

 LP%VEG_EMISS = 4.*SIGMA*LP%VEG_KAPPA*TMP_VEG_NEW**4 !used in RTE solver

 MPV_CHAR_LOSS_TOTAL  = MPV_CHAR_LOSS_TOTAL  + MPV_CHAR_LOSS !needed for subcycling
 MPV_MOIST_LOSS_TOTAL = MPV_MOIST_LOSS_TOTAL + MPV_MOIST_LOSS !needed for subcycling
 MPV_VOLIT_TOTAL      = MPV_VOLIT_TOTAL      + MPV_VOLIT !needed for subcycling
 MPV_CHAR_CO2_TOTAL   = MPV_CHAR_CO2_TOTAL   + MPV_CHAR_CO2
 MPV_CHAR_O2_TOTAL    = MPV_CHAR_O2_TOTAL    + MPV_CHAR_O2
 MPV_ADDED = MPV_ADDED + MPV_MOIST_LOSS + MPV_VOLIT + MPV_CHAR_CO2
 !MPV_ADDED = MPV_ADDED + MPV_MOIST_LOSS + MPV_VOLIT + MPV_CHAR_CO2 - MPV_CHAR_O2

! Check if critical mass flux condition is met
!IF (MPV_ADDED*RDT_FE < VEG_CRITICAL_MASSSOURCE .AND. .NOT. LP%VEG_IGNITED) THEN
! MPV_ADDED      = 0.0_EB
! MW_AVERAGE     = 0.0_EB
! MPV_MOIST_LOSS = 0.0_EB
! MPV_VOLIT      = 0.0_EB
! Q_VEG_MOIST    = 0.0_EB
! Q_VEG_VOLIT    = 0.0_EB
!ELSE
! LP%VEG_IGNITED = .TRUE.
!ENDIF

! Add affects of fuel element thermal degradation of vegetation to velocity divergence

 CALL GET_SPECIFIC_HEAT(ZZ_GET,CP_GAS,TMP_GAS)
 RCP_GAS    = 1._EB/CP_GAS
 !MW_TERM    = MW_VEG_MOIST_TERM + MW_VEG_VOLIT_TERM
 MW_AVERAGE = R0/RSUM(II,JJ,KK)/RHO_GAS*(MW_VEG_MOIST_TERM + MW_VEG_VOLIT_TERM)
 Q_ENTHALPY = Q_VEG_MOIST + Q_VEG_VOLIT - (1.0_EB - LPC%VEG_CHAR_ENTHALPY_FRACTION)*Q_VEG_CHAR

 !D_LAGRANGIAN(II,JJ,KK) = D_LAGRANGIAN(II,JJ,KK)  +           & 
 !                         (-QCON_VEG*RCP_GAS + Q_ENTHALPY*RCP_GAS)/(RHO_GAS*TMP_GAS) + &
 !                         RDT_FE*MW_AVERAGE 
 D_SOURCE(II,JJ,KK) = D_SOURCE(II,JJ,KK)  + (RDT_FE*Q_ENTHALPY*RCP_GAS/(RHO_GAS*TMP_GAS) + RDT_FE*MW_AVERAGE) 


 TMP_VEG   = TMP_VEG_NEW
 D_SOURCE(II,JJ,KK) = D_SOURCE(II,JJ,KK) + (-QCON_VEG*RCP_GAS)/(RHO_GAS*TMP_GAS)
!M_DOT_PPP(II,JJ,KK,2) = MPV_VOLIT*RDT_FE
!M_DOT(2,NM) = MPV_VOLIT*V_CELL
!mpv_volit_sum = mpv_volit_sum + mpv_volit
!print '(A,1x,1ES13.4)','mpv_volit_sum = ',mpv_volit_sum

 IF (MPV_MOIST <= MPV_MOIST_MIN) THEN !for time sub cycling
   MPV_MOIST = 0.0_EB
   MW_VEG_MOIST_TERM = 0.0_EB
  Q_VEG_MOIST = 0.0_EB
 ENDIF
 IF (MPV_VEG <= MPV_VEG_MIN) THEN
   MPV_VOLIT = 0.0_EB
   MPV_VEG   = 0.0_EB
   MW_VEG_VOLIT_TERM = 0.0_EB
   Q_VEG_VOLIT = 0.0_EB
 ENDIF

!mpv_volit_total = 0.0_EB
!if (t<=10._EB) mpv_volit_total = 0.1_EB*dt_fe
!mpv_added = mpv_volit_total
!zz(ii,jj,kk,5) = 0.1
!print '(A,1x,2ES13.4)','vege:mpv_added,mpv_volit_total',mpv_added,mpv_volit_total

! Add water vapor, fuel vapor, and CO2 mass to total density
! MPV_ADDED     = MPV_MOIST_LOSS + MPV_VOLIT + MPV_CHAR_CO2
  LP%VEG_MLR    = MPV_ADDED*RDT_FE !kg/m^3/s used in FVX,FVY,FVZ along with drag in part.f90
  RHO(II,JJ,KK) = RHO_GAS + MPV_ADDED
  RRHO_GAS_NEW  = 1._EB/RHO(II,JJ,KK)
! print*,'NM =',NM
! print*,'** ',rho(ii,jj,kk)

!Need to have the following in the input file in order to inject water vapor from drying and CO2 from
!char oxidation, then n=I_WATER=4, n=I_CO2=5,n=3 are lumped combustion products CO2 and H2O, n=I_FUEL=2 is
!pyrolsis gas
!&SPEC ID='WATER VAPOR' / which is n=4
!&SPEC ID='CARBON DIOXIDE' / which is n=5

! Add gas species created by degradation of vegetation Yi_new = (Yi_old*rho_old + change in rho_i)/rho_new
! Add water vapor mass from drying to water vapor mass fraction
  IF (I_WATER > 0) THEN 
!  ZZ(II,JJ,KK,I_WATER) = ZZ(II,JJ,KK,I_WATER) +  MPV_MOIST_LOSS*RRHO_GAS_NEW
   ZZ(II,JJ,KK,I_WATER) = ZZ(II,JJ,KK,I_WATER) + (MPV_MOIST_LOSS_TOTAL - MPV_ADDED*ZZ(II,JJ,KK,I_WATER))*RRHO_GAS_NEW
!  ZZ(II,JJ,KK,I_WATER) = MIN(1._EB,ZZ(II,JJ,KK,I_WATER))
!  DMPVDT_FM_VEG(II,JJ,KK,I_WATER) = DMPVDT_FM_VEG(II,JJ,KK,I_WATER) + RDT_FE*MPV_MOIST_LOSS
  ENDIF

! Add fuel vapor mass from pyrolysis to fuel mass fraction
  I_FUEL = REACTION(1)%FUEL_SMIX_INDEX
  IF (I_FUEL /= 0) THEN 
!  ZZ(II,JJ,KK,I_FUEL) = ZZ(II,JJ,KK,I_FUEL) + MPV_VOLIT*RRHO_GAS_NEW
   ZZ(II,JJ,KK,I_FUEL) = ZZ(II,JJ,KK,I_FUEL) + (MPV_VOLIT_TOTAL - MPV_ADDED*ZZ(II,JJ,KK,I_FUEL))*RRHO_GAS_NEW
!  ZZ(II,JJ,KK,I_FUEL) = MIN(1._EB,ZZ(II,JJ,KK,I_FUEL))
!  DMPVDT_FM_VEG(II,JJ,KK,I_FUEL) = DMPVDT_FM_VEG(II,JJ,KK,I_FUEL) + RDT_FE*MPV_VOLIT
  ENDIF

! Add CO2 mass, due to production during char oxidation, to CO2 mass fraction
  IF (I_CO2 /= -1 .AND. LPC%VEG_CHAR_OXIDATION) THEN 
   ZZ(II,JJ,KK,I_CO2) = ZZ(II,JJ,KK,I_CO2) + (MPV_CHAR_CO2_TOTAL - MPV_ADDED*ZZ(II,JJ,KK,I_CO2))*RRHO_GAS_NEW
  ENDIF

! Remove O2 from gas due to char oxidation ***  this was incorrect and cannot be done because oygen is combined with
! CO2 and H2O to make the lumped AIR species
!IF (LPC%VEG_CHAR_OXIDATION) THEN 
! ZZ(II,JJ,KK,O2_INDEX) = ZZ(II,JJ,KK,O2_INDEX) - (MPV_CHAR_O2_TOTAL - MPV_ADDED*ZZ(II,JJ,KK,I_CO2))*RRHO_GAS_NEW
! ZZ(II,JJ,KK,O2_INDEX) = MAX(0.0_EB,ZZ(II,JJ,KK,O2_INDEX))
!ENDIF

!print '(A,1x,6I3)','vege:n_tracked_species,i_water,h2o_index,i_fuel,i_co2,os_index',n_tracked_species,i_water, &
!                            h2o_index,i_fuel,i_co2,o2_index

! WRITE(9998,'(A)')'T,TMP_VEG,QCON_VEG,QRAD_VEG'
!IF (II==0.5*IBAR .AND. JJ==0.5*JBAR .AND. KK==0.333*KBAR) THEN
!IF (II==12 .AND. JJ==12 .AND. KK==4) THEN 
!IF (II==20 .AND. JJ==20 .AND. KK==25) THEN !M=14% and 49% element burnout
!IF (II==27 .AND. JJ==20 .AND. KK==7) THEN !M=49% not full element burnout
! WRITE(9998,'(9(ES12.4))')T,TMP_GAS,TMP_VEG,QCON_VEG,QRAD_VEG,LP%VEG_MOIST_MASS,LP%VEG_FUEL_MASS, &
!                          MPV_MOIST_LOSS_MAX*RDT_FE,MPV_VOLIT_MAX*RDT_FE
!ENDIF

! V_VEG               = V_VEG + V_CELL
! TOTAL_MASS_MOIST    = TOTAL_MASS_MOIST + LP%VEG_MOIST_MASS*V_CELL
! TOTAL_MASS_DRY_FUEL = TOTAL_MASS_DRY_FUEL + LP%VEG_FUEL_MASS*V_CELL

 ENDIF THERMAL_CALC  ! end of thermally thin heat transfer, etc. calculations

! Fill arrays for outputting vegetation variables when OUTPUT_TREE=.TRUE.
! N_TREE = LP%VEG_N_TREE_OUTPUT
! IF (N_TREE /= 0) THEN
!  TREE_OUTPUT_DATA(N_TREE,1,NM) = TREE_OUTPUT_DATA(N_TREE,1,NM) + LP%TMP - 273._EB !C
!  TREE_OUTPUT_DATA(N_TREE,2,NM) = TREE_OUTPUT_DATA(N_TREE,2,NM) + TMP_GAS - 273._EB !C
!  TREE_OUTPUT_DATA(N_TREE,3,NM) = TREE_OUTPUT_DATA(N_TREE,3,NM) + LP%VEG_FUEL_MASS*V_CELL !kg
!  TREE_OUTPUT_DATA(N_TREE,4,NM) = TREE_OUTPUT_DATA(N_TREE,4,NM) + LP%VEG_MOIST_MASS*V_CELL !kg
!  TREE_OUTPUT_DATA(N_TREE,5,NM) = TREE_OUTPUT_DATA(N_TREE,5,NM) + LP%VEG_CHAR_MASS*V_CELL !kg
!  TREE_OUTPUT_DATA(N_TREE,6,NM) = TREE_OUTPUT_DATA(N_TREE,6,NM) + LP%VEG_ASH_MASS*V_CELL !kg
!  TREE_OUTPUT_DATA(N_TREE,7,NM) = TREE_OUTPUT_DATA(N_TREE,7,NM) + LP%VEG_DIVQC*V_CELL*0.001_EB !kW
!  TREE_OUTPUT_DATA(N_TREE,8,NM) = TREE_OUTPUT_DATA(N_TREE,8,NM) + LP%VEG_DIVQR*V_CELL*0.001_EB !kW
!  TREE_OUTPUT_DATA(N_TREE,9,NM) = TREE_OUTPUT_DATA(N_TREE,9,NM) + 1._EB !number of particles
!  TREE_OUTPUT_DATA(N_TREE,10,NM) = TREE_OUTPUT_DATA(N_TREE,10,NM) + MPV_CHAR_LOSS_TOTAL*V_CELL !kg 
!  TREE_OUTPUT_DATA(N_TREE,11,NM) = TREE_OUTPUT_DATA(N_TREE,11,NM) - Q_VEG_CHAR_TOTAL*V_CELL*RDT_FE*0.001_EB !kW

!! TREE_OUTPUT_DATA(N_TREE,10,NM) = TREE_OUTPUT_DATA(N_TREE,10,NM) + NUSS_HILPERT_CYL_FORCEDCONV
!! TREE_OUTPUT_DATA(N_TREE,11,NM) = TREE_OUTPUT_DATA(N_TREE,11,NM) + NUSS_MORGAN_CYL_FREECONV 

!! TREE_OUTPUT_DATA(N_TREE,4,NM) = TREE_OUTPUT_DATA(N_TREE,4,NM) + LP%VEG_PACKING_RATIO
!! TREE_OUTPUT_DATA(N_TREE,5,NM) = TREE_OUTPUT_DATA(N_TREE,5,NM) + LP%VEG_SV

! ENDIF

ENDDO PARTICLE_LOOP

!print*,'--------------------------------'
!print '(A,1x,I2,1x,ES12.4)','vege:nm,tree_output divqc ',nm,tree_output_data(1,7,nm)

! Write out total bulk
!TOTAL_BULKDENS_MOIST = TOTAL_MASS_MOIST/V_VEG
!TOTAL_BULKDENS_DRY_FUEL = TOTAL_MASS_DRY_FUEL/V_VEG
!WRITE(9999,'(5(ES12.4))')T,TOTAL_BULKDENS_DRY_FUEL,TOTAL_BULKDENS_MOIST,TOTAL_MASS_DRY_FUEL,TOTAL_MASS_MOIST

!VEG_TOTAL_DRY_MASS(NM)   = TOTAL_MASS_DRY_FUEL
!VEG_TOTAL_MOIST_MASS(NM) = TOTAL_MASS_MOIST

! Remove vegetation that has completely burned (i.e., LP%R has been set equal to zero)
CALL REMOVE_PARTICLES(T,NM)

VEG_CLOCK_FE = T
 
END SUBROUTINE RAISED_VEG_MASS_ENERGY_TRANSFER

!> \brief Calculate the Rothermel no-wind, no-slope rate of spread.
!> 
!>
!> \details The Rothermel model as described in Bachmann's thesis.

REAL(EB) FUNCTION ROS_NO_WIND_NO_SLOPE(ROTHERMEL_FUEL_INDEX,SURF_INDEX)

INTEGER, INTENT(IN) :: ROTHERMEL_FUEL_INDEX,SURF_INDEX
REAL(EB) :: w0d1, w0d2, w0d3, w0lh, w0lw, md1, md2, md3, mlh, mlw, svd1, svd2, svd3, svlh, svlw, depth, rhop, heat, st, se, mx
REAL(EB) :: swd1, swd2, swd3, swlh, swlw, swd, swl, swt, s2wt, sw2d, sw2l, swmd, swml, sigma, rhob, beta, &
            betaOpt, wnd, wnl, hnd1, hnd2, hnd3, hnlh, hnlw, hnd, hnl, bigW, hnmd, mfdead, mxlive, rml, rmd, etaMd, etaMl, etaM, &
            etas, gammaMax, bigA, gamma, bigIr, xi, epsd1, epsd2, epsd3, epslh, epslw, bigQd1, bigQd2, &
            bigQd3, bigQlh, bigQlw, hskz, hsk
TYPE(SURFACE_TYPE), POINTER :: SF

SF => SURFACE(SURF_INDEX)

md1 = SF%VEG_LSET_M1
md2 = SF%VEG_LSET_M10
md3 = SF%VEG_LSET_M100
mlw = SF%VEG_LSET_MLW
mlh = SF%VEG_LSET_MLH
               
SELECT CASE(ROTHERMEL_FUEL_INDEX)
   CASE(1)  ! 'Short Grass'
      w0d1=0.1659     ; w0d2=0.        ; w0d3=0.        ; w0lh=0.        ; w0lw=0. 
      svd1=11483.     ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.12         ; depth=0.3048   ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(2)  ! 'Timbergrass'
      w0d1=0.448      ; w0d2=0.224     ; w0d3=0.112     ; w0lh=0.112     ; w0lw=0. 
      svd1=9842.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.15         ; depth=0.3048   ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(3)  ! 'Tall Grass'
      w0d1=0.675      ; w0d2=0.        ; w0d3=0.        ; w0lh=0.        ; w0lw=0. 
      svd1=4921.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.25         ; depth=0.762    ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(4)  ! 'Chaparral'
      w0d1=1.123      ; w0d2=0.899     ; w0d3=0.448     ; w0lh=1.123     ; w0lw=0. 
      svd1=6562.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.20         ; depth=1.829    ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(5)  ! 'Brush'
      w0d1=0.224      ; w0d2=0.112     ; w0d3=0.        ; w0lh=0.        ; w0lw=0.448
      svd1=6562.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.20         ; depth=0.6096   ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(6)  ! 'Dormant Brush'
      w0d1=0.336      ; w0d2=0.56      ; w0d3=0.448     ; w0lh=0.        ; w0lw=0.   
      svd1=5741.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.25         ; depth=0.762    ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(7)  ! 'Southern Rough'
      w0d1=0.255      ; w0d2=0.419     ; w0d3=0.336     ; w0lh=0.        ; w0lw=0.083
      svd1=5741.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.40         ; depth=0.762    ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(8)  ! 'Closed Timber Litter'
      w0d1=0.336      ; w0d2=0.224     ; w0d3=0.56      ; w0lh=0.        ; w0lw=0.   
      svd1=6562.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.30         ; depth=0.06096  ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(9)  ! ID='Hardwood Litter'
      w0d1=0.655      ; w0d2=0.092     ; w0d3=0.034     ; w0lh=0.        ; w0lw=0.   
      svd1=8202.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.25         ; depth=0.06096  ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(10)  ! 'Timber'
      w0d1=0.675      ; w0d2=0.448     ; w0d3=1.123     ; w0lh=0.        ; w0lw=0.448
      svd1=6562.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.25         ; depth=0.3048   ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(11)  ! 'Light Slash'
      w0d1=0.336      ; w0d2=1.011     ; w0d3=1.235     ; w0lh=0.        ; w0lw=0.   
      svd1=4921.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.15         ; depth=0.3048   ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(12)  ! ID='Medium Slash'
      w0d1=0.899      ; w0d2=3.145     ; w0d3=3.706     ; w0lh=0.        ; w0lw=0.   
      svd1=4921.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.20         ; depth=0.70104  ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
   CASE(13)  ! 'Heavy Slash'
      w0d1=1.571      ; w0d2=5.165     ; w0d3=6.288     ; w0lh=0.        ; w0lw=0.   
      svd1=4921.      ; svd2=358.      ; svd3=98.       ; svlh=4921.     ; svlw=4921. 
      mx=0.25         ; depth=0.9144   ; rhop=512.      ; heat=18607.    ; st=0.0555      ; se=0.01
END SELECT

SF%VEG_LSET_HT = depth
      
! Auxiliary functions
   
swd1 = svd1*w0d1  
swd2 = svd2*w0d2 
swd3 = svd3*w0d3  
swlh = svlh*w0lh  
swlw = svlw*w0lw 
swd  = swd1 + swd2 + swd3  
swl  = swlh + swlw  
swt  = swd + swl  
s2wt = svd1**2*w0d1 + svd2**2*w0d2 + svd3**2*w0d3 + svlh**2*w0lh + svlw**2*w0lw  
sw2d = svd1*w0d1**2 + svd2*w0d2**2 + svd3*w0d3**2
sw2l = svlh*w0lh**2 + svlw*w0lw**2  
swmd = swd1*md1 + swd2*md2 + swd3*md3  
swml = swlh*mlh + swlw*mlw
   
! Characteristic surface-to-volume ratio [R(71,72)]

sigma = s2wt/swt
SF%VEG_LSET_SIGMA = sigma*0.01  ! Convert from 1/m to 1/cm
   
! Mean bulk density [R(74)]
   
rhob = (w0d1 + w0d2 + w0d3 + w0lh + w0lw)/depth
  
! Mean packing ratio [R(31,73)]
   
beta = rhob/rhop
SF%VEG_LSET_BETA = beta
   
! Optimal packing ratio [R(37)]
   
betaOpt = 8.8578*sigma**(-0.8189)
   
! Net fuel loading [R(60), adjusted by A.(p.88) and R(59)]
   
if (swd==0._eb) then
   wnd = 0._eb
else
   wnd = (sw2d/swd)*(1. - st)
endif

if (swl==0._eb) then
   wnl = 0._eb
else
   wnl = (sw2l/swl)*(1. - st)
endif
   
! Mineral damping coefficient [R(62)]
   
etas = 0.174*se**(-0.19)
   
! Ratio of "fine" fuel loadings,dead/living [Albini,p.89]
   
hnd1 = 0.20482*w0d1*Exp(-452.76/svd1) 
hnd2 = 0.20482*w0d2*exp(-452.76/svd2) 
hnd3 = 0.20482*w0d3*exp(-452.76/svd3) 
hnlh = 0.20482*w0lh*exp(-1640.42/svlh) 
hnlw = 0.20482*w0lw*exp(-1640.42/svlw) 
hnd = hnd1 + hnd2 + hnd3 
hnl = hnlh + hnlw
if (swl==0._eb) then
   bigW = 0._eb
else
   bigW = hnd/hnl
endif
   
! Moisture content of "fine" dead fuel [Albini,p.89]
   
hnmd   = hnd1*md1 + hnd2*md2 + hnd3*md3
mfdead = hnmd/hnd
   
! Moisture of extinction of living fuel [R(88),Albini,p.89]
   
mxlive = 2.9*bigW*(1.0 - (mfdead/mx)) - 0.226
   
! Moisture ratios [R(65,66)]
   
if (swl==0._eb) then
   rml = 0._eb
else
   rml = swml/(swl*mxlive)
endif

rmd = swmd/(swd*mx)
   
! Moisture damping coefficients [R(64)]
   
etaMd = 1.0 - (2.59*rmd) + (5.11*rmd**2) - (3.52*rmd**3) 
etaMl = 1.0 - (2.59*rml) + (5.11*rml**2) - (3.52*rml**3) 
etaM  = wnd*etaMd + wnl*etaMl
   
! Maximum reaction velocity [R(36,68)]
   
gammaMax = (0.16828*sigma**(1.5))/(29700 + 0.5997*sigma**(1.5))
   
! A [R(70),Albini p.88]
   
bigA = 340.53*sigma**(-0.7913)
   
! Potential reaction velocity [R(38)]
   
gamma = gammaMax*(beta/betaOpt)**(bigA)*exp(bigA*(1.0 - (beta/betaOpt)))
   
! Propagating flux ratio [R(42)]
   
xi = exp((0.792 + 0.37597*sqrt(sigma))*(beta + 0.1))/(192.0 + 0.0791*sigma)
   
! Effective heating number [R(14,77)]
   
epsd1 = exp(-452.76/svd1) 
epsd2 = exp(-452.76/svd2) 
epsd3 = exp(-452.76/svd3)
epslh = exp(-452.76/svlh) 
epslw = exp(-452.76/svlw)
   
! Heat of pre-ignition [R(12,78)]
   
bigQd1 = 581.5 + 2595.7*md1 
bigQd2 = 581.5 + 2595.7*md2 
bigQd3 = 581.5 + 2595.7*md3 
bigQlh = 581.5 + 2595.7*mlh 
bigQlw = 581.5 + 2595.7*mlw
   
! Heat sink [R(77)]
   
hskz = svd1*w0d1*epsd1*bigQd1 + svd2*w0d2*epsd2*bigQd2 + svd3*w0d3*epsd3*bigQd3 + svlh*w0lh*epslh*bigQlh + svlw*w0lw*epslw*bigQlw
hsk  = rhob*hskz/swt
   
! Reaction intensity [R(27,58),Albini,p.89]
   
bigIr = gamma*heat*etas*etaM
   
! Rate of spread [R(52)] and the rate of spread in the absence of wind and with no slope.
   
ROS_NO_WIND_NO_SLOPE = (bigIr*xi)/hsk

END FUNCTION ROS_NO_WIND_NO_SLOPE

END MODULE VEGE
