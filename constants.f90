!===============================================================================
! Модуль констант и параметров
!===============================================================================
MODULE CONSTANTS
  IMPLICIT NONE
  REAL(8), PARAMETER :: PI = 4.0_8 * DATAN(1.0_8)
  
  ! Параметры сетки
  INTEGER, PARAMETER :: N1cv_core = 20
  INTEGER, PARAMETER :: N1cv_shell = 50
  INTEGER, PARAMETER :: N1 = N1cv_shell + N1cv_core + 2
  INTEGER, PARAMETER :: N1_core = N1cv_core + 2
  INTEGER, PARAMETER :: N2 = 72
  INTEGER, PARAMETER :: N3 = 142
  INTEGER, PARAMETER :: Nmax = 142
  
  ! Параметры для PSI_OUT
  INTEGER, PARAMETER :: N1_PSI_OUT = 52
  
  ! Параметры для PSI_IN
  INTEGER, PARAMETER :: N1_PSI_IN = 32
  
  ! Численные параметры
  REAL(8), PARAMETER :: TOL_SMAX = 1.0D-3
  REAL(8), PARAMETER :: TOL_SMAX_CORRECT = 1.0D-7
  REAL(8), PARAMETER :: TOL_P = 1.0D-6
  REAL(8), PARAMETER :: TOL_T = 1.0D-6
  REAL(8), PARAMETER :: TOL_UVW = 1.0D-6
  
END MODULE CONSTANTS
