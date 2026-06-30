!===============================================================================
! Модуль переменных для DVMH (распределенные массивы)
! Заменяет VAR, VAR_d, VAR_PSI_OUT, VAR_PSI_OUT_d, VAR_PSI_IN, VAR_PSI_IN_d
!===============================================================================
MODULE VAR_DVMH
  USE CONSTANTS
  USE DVMPARA
  IMPLICIT NONE
  
  ! Параметры задачи
  INTEGER :: LAST, IWRITE, TIME_ITER
  REAL(8) :: TIME, DT
  INTEGER :: NTIMES_P, NTIMES_PC, NTIMES_T, NTIMES_UVW, NTIMES_F1, NTIMES_F2, NTIMES_F3
  LOGICAL :: FirstCallMP
  
  ! Физические параметры
  REAL(8) :: XL1, XL2, XL3, X1_0, X2_0, X3_0
  REAL(8) :: viscm, tolPSI, tolPSI_C, SMAX_B, sum_sq_f
  REAL(8) :: Pm, Br_out_avg, Br_in_avg, thc, rho_c, r_i, r_o
  REAL(8) :: Ranum, Pr, E, visk, dens, SMAX, TOL_P, TOL_T, TOL_UVW
  REAL(8) :: Ha, B0, Umax, omega_ic, resF1, resF2, resF3
  
  !$DVM NODES(NODES_COUNT)
  
  ! Распределенные одномерные массивы
  REAL(8), DIMENSION(:), ALLOCATABLE :: Xs1, X1, Xdif1, Kx1
  REAL(8), DIMENSION(:), ALLOCATABLE :: Xs2, X2, Xdif2, Kx2
  REAL(8), DIMENSION(:), ALLOCATABLE :: Xs3, X3, Xdif3, Kx3, X1_out
  REAL(8), DIMENSION(:), ALLOCATABLE :: cosXs2, cosX2, sinXs2, sinX2
  REAL(8), DIMENSION(:), ALLOCATABLE :: X1cv, X2cv, X3cv
  REAL(8), DIMENSION(:), ALLOCATABLE :: VVY1, VVY2, VVX1, VVX2, VYX1, VYX2, APX
  REAL(8), DIMENSION(:), ALLOCATABLE :: WWY1, WWY2, UYY, WZY1, WZY2, Spv1J, Spvy1, Spvy2
  REAL(8), DIMENSION(:), ALLOCATABLE :: VZY1, VZY2, Spwy1, Spwy2
  REAL(8), DIMENSION(:), ALLOCATABLE :: VVX1S, VVX2S, XsDifP4, XDifP4
  REAL(8), DIMENSION(:), ALLOCATABLE :: WWY1S, WWY2S, VZY1S, VZY2S, E1_Oz, X_3
  
  ! Распределенные трехмерные массивы полей
  !$DVM ARRAY U1(N1,N2,N3), U2(N1,N2,N3), U3(N1,N2,N3)
  !$DVM ARRAY T(N1,N2,N3), T0(N1,N2,N3), P(N1,N2,N3), PC(N1,N2,N3)
  !$DVM ARRAY F1(N1,N2,N3), F2(N1,N2,N3), F3(N1,N2,N3)
  !$DVM ARRAY S1(N1,N2,N3), S2(N1,N2,N3), S3(N1,N2,N3)
  !$DVM ARRAY B1(N1,N2,N3), B2(N1,N2,N3), B3(N1,N2,N3), PSI(N1,N2,N3)
  !$DVM ARRAY AIM(N1,N2,N3), AIP(N1,N2,N3), AJM(N1,N2,N3), AJP(N1,N2,N3)
  !$DVM ARRAY AKM(N1,N2,N3), AKP(N1,N2,N3), AP(N1,N2,N3), CON(N1,N2,N3)
  
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: &
    U1, U2, U3, T, T0, P, PC, &
    F1, F2, F3, S1, S2, S3, &
    U1_0, U2_0, U3_0, U1hat, U2hat, U3hat, DU1, DU2, DU3, &
    S1u1, S2u1, S2u1_1, S3u1, S1u2, S2u2, S3u2, S1u3, S2u3, &
    B1, B2, B3, PSI, A1, A2, A3, &
    VCV, VCV1, VCV2, VCV3, &
    F1_0, F2_0, F3_0, F1a, F2a, F3a, &
    U1a, U2a, U3a, B1a, B2a, B3a, TEMPS, &
    AIM, AIP, AJM, AJP, AKM, AKP, AP, CON, &
    D213, D231, D312, D321, D123, D132, &
    F_213, F_231, F_312, F_321, F_123, F_132, &
    J213, J231, J312, J321, J123, J132, &
    J1, J2, J3, PSI_C
  
  ! Двумерные граничные массивы
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: B2_o, B2_i, B3_o, B3_i
  
CONTAINS

  !------------------------------------------------------
  SUBROUTINE ALLOC_DVMH_MEMS
    IMPLICIT NONE
    INTEGER :: IERR
    
    ! Выделение памяти для распределенных массивов DVMH
    ! DVM автоматически управляет размещением на узлах
    ALLOCATE(Xs1(N1), X1(N1), Xdif1(N1), Kx1(N1), X_3(0:N3+1), &
             Xs2(N2), X2(N2), Xdif2(N2), Kx2(N2), &
             Xs3(N3), X3(N3), Xdif3(N3), Kx3(N3), &
             cosXs2(N2), cosX2(N2), sinXs2(N2), sinX2(N2), &
             X1cv(N1), X2cv(N2), X3cv(N3), E1_Oz(N1), &
             VVY1(N2), VVY2(N2), VVX1(N1), VVX2(N1), VYX1(N1), &
             VYX2(N1), APX(N1), WWY1(N2), WWY2(N2), UYY(N2), &
             WZY1(N2), WZY2(N2), Spv1J(N2), Spvy1(N2), Spvy2(N2), &
             VZY1(N2), VZY2(N2), Spwy1(N2), Spwy2(N2), &
             VVX1S(N1), VVX2S(N1), XsDifP4(N1), XDifP4(N1), &
             WWY1S(N2), WWY2S(N2), VZY1S(N2), VZY2S(N2), STAT = IERR)
    
    ! Трехмерные массивы
    ALLOCATE(U1(N1,N2,N3), U2(N1,N2,N3), U3(N1,N2,N3), &
             T(N1,N2,N3), T0(N1,N2,N3), P(N1,N2,N3), PC(N1,N2,N3), &
             F1(N1,N2,N3), F2(N1,N2,N3), F3(N1,N2,N3), &
             S1(N1,N2,N3), S2(N1,N2,N3), S3(N1,N2,N3), STAT = IERR)
    
    ! Остальные массивы...
    
  END SUBROUTINE ALLOC_DVMH_MEMS
  
  !------------------------------------------------------
  SUBROUTINE DEALLOC_DVMH_MEMS
    IMPLICIT NONE
    ! DVM автоматически освобождает распределенную память
    DEALLOCATE(Xs1, X1, Xdif1, Kx1, X_3)
    DEALLOCATE(U1, U2, U3, T, T0, P, PC, F1, F2, F3, S1, S2, S3)
  END SUBROUTINE DEALLOC_DVMH_MEMS
  
END MODULE VAR_DVMH
