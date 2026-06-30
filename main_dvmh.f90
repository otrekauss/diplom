!===============================================================================
! Главная программа CVMHD с использованием DVMH
!===============================================================================
PROGRAM MAIN_DVMH
  USE CONSTANTS
  USE VAR_DVMH
  USE DVMPARA
  USE FD_DVMH, ONLY: SIMPLER
  USE MD_DVMH, ONLY: CorrectB, MHD, ComputeSmaxB, ComputeMagPotential, &
                     ComputeBoundCondsF, FtoB
  USE GEOM_DVMH
  USE USER_DVMH
  USE OMP_LIB
  IMPLICIT NONE
  
  REAL(8) :: startTime, endTime
  INTEGER :: J, K
  
  ! Инициализация DVM
  CALL DVM_INIT()
  
  startTime = omp_get_wtime()
  
  ! Настройка задачи
  CALL SET_TASK_PARAMS
  CALL BUILD_GRID
  
  ! Выделение памяти (DVM автоматически распределяет по узлам)
  CALL ALLOC_DVMH_MEMS
  
  ! Основной цикл времени (раскомментировать для полноценной работы)
  ! DO TIME_ITER = 1, LAST
  !   TIME = TIME + DT
  !   CALL SIMPLER
  !   
  !   IF (SMAX_B >= TOL_SMAX_CORRECT) THEN
  !     CALL CorrectB 
  !   END IF
  !   
  !   CALL F0_F
  !   CALL SAVE_DATA
  !   CALL OUTPUT
  ! END DO
  
  endTime = omp_get_wtime()
  WRITE(*,'(A, 1P1E16.4)') 'Elapsed time:', endTime - startTime
  
  ! Вывод результатов
  CALL OUTPUT_GRID
  
  ! Освобождение памяти
  CALL DEALLOC_DVMH_MEMS
  
  ! Завершение DVM
  CALL DVM_FINALIZE()
  
END PROGRAM MAIN_DVMH

!------------------------------------------------------
SUBROUTINE F0_F
  USE VAR_DVMH
  IMPLICIT NONE
  
  !$DVM TASK SHARED
  F1_0 = F1
  F2_0 = F2
  F3_0 = F3
  U1_0 = U1
  U2_0 = U2
  U3_0 = U3
  T0 = T
  !$DVM END TASK
  
END SUBROUTINE F0_F

!------------------------------------------------------
SUBROUTINE FF1N1(F, N1) 
  USE CONSTANTS
  USE VAR_DVMH
  IMPLICIT NONE
  INTEGER :: I, J, N1
  
  !$DVM TASK SHARED PRIVATE(I,J)
  DO J = 1, N2
    DO I = 1, N1
      F(I,J,1)  = F(I,J,N3-1)
      F(I,J,N3) = F(I,J,2)  
    ENDDO
  ENDDO
  !$DVM END TASK
  
END SUBROUTINE FF1N1
