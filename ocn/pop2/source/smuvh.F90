#define LOGMSG()
!write(mytid+600,'(a,3i4)')"SMUV",__LINE__,k,j
!  CVS: $Id: smuvh.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ========================
      SUBROUTINE SMUV (X,Z,KK,fil_lat)
!     ========================
!     1-D zonal smoother

#include <def-undef.h>
use precision_mod
use param_mod
use msg_mod
use pconst_mod, only: sinu,pi,torad
use domain

      IMPLICIT NONE
      INTEGER :: JFS1,JFS2,JFN1,JFN2,KK, NCY
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT,KK,max_blocks_clinic),XS (IMT),Z (IMT,JMT,KM,max_blocks_clinic)
      INTEGER     :: NN(JMT), MAX_NN, iblock

!lhl      fil_lat=55.D0
!
      MAX_NN = 0
      do j =jst,jmt
         if (sinu(j).le.cos(fil_lat*torad)) then
            NN(j) = int(cos(fil_lat*torad)/sinu(j)*1.2D0)
            if (NN(j) .gt. MAX_NN) MAX_NN = NN(J)
         else 
            NN(j) = 0
         endif
      enddo

!
      DO NCY = 1,MAX_NN
!
!$OMP PARALLEL DO PRIVATE (iblock,k,j)
   do iblock = 1, nblocks_clinic
         do j =jst,jmt
            if (NN(j) .ge. NCY) then       
               DO K = 1,KK
                  DO I = 1,IMT
                     XS (I)= X (I,J,K,iblock)* Z (I,J,K,iblock)
                  END DO
                  DO I = 2,IMM
                     X (I,J,K,iblock)= (0.5D0*XS(I)+0.25D0*(XS(I-1)+XS(I+1)))*Z(I,J,K,iblock)
                  END DO
               ENDDO
            endif
         END DO
     END DO
!
   END DO

      RETURN
      END SUBROUTINE SMUV



!     ========================
      SUBROUTINE SMTS (X,Z,KK,fil_lat)
!     ========================
!     1-D zonal smoother

#include <def-undef.h>
use precision_mod
use param_mod
use msg_mod
use pconst_mod,only: sint,PI,torad
use domain
      IMPLICIT NONE

      INTEGER :: JFS1,JFS2,JFN1,JFN2,KK, NCY
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT,KK,max_blocks_clinic),XS (IMT),Z (IMT,JMT,KM,max_blocks_clinic)
      INTEGER     :: NN(JMT), MAX_NN, iblock


!     fil_lat=66.D0
!lhl      fil_lat=56.D0
!

      MAX_NN = 0
      do j =jst,jmt
         if (sint(j).le.cos(fil_lat*torad)) then
            NN(j) = int(cos(fil_lat*torad)/sint(j)*1.2D0)
            if (NN(j) .gt. MAX_NN) MAX_NN = NN(J)
         else
            NN(j) = 0
         endif
      enddo

!!$OMP PARALLEL DO PRIVATE (iblocks,K,J,xs)
   do iblock = 1, nblocks_clinic
      DO NCY = 1,MAX_NN
         do j =jst,jmt
            if (NN(j) .ge. NCY) then
               DO K = 1,KK
                  DO I = 1,IMT
                     XS (I)= X (I,J,K,iblock)* Z (I,J,K,iblock)
                  END DO
                  DO I = 2,IMM
                     X(I,J,K,iblock)=(XS(I)*(1.0D0-0.25D0*Z(I-1,J,K,iblock)-0.25D0*Z(I+1,J,K,iblock)) &
                               +0.25D0*(XS(I-1)+XS(I+1)))*Z(I,J,K,iblock)
                  END DO
               ENDDO
!              if (nx_proc == 1) then
!                 DO K = 1,KK
!                    X (1,J,K)= X (IMM,J,K)
!                    X (IMT,J,K)= X (2,J,K)
!                 ENDDO
!              endif
            endif
         enddo
         if (nx_proc .ne. 1) then
!           call exchange_2D_boundary(x,kk,NN,NCY,0)
!           call exchange_2D_boundary(x,kk,NN,NCY,1)
         end if
     END DO
  end do


      RETURN
      END SUBROUTINE SMTS



!     ========================
      SUBROUTINE SMZ0 (X,Z,fil_lat)
!     ========================
!     1-D zonal smoother

#include <def-undef.h>
use precision_mod
use param_mod
use msg_mod
use domain
use pconst_mod,only:sint,pi,torad
      IMPLICIT NONE

      INTEGER :: JFS1,JFS2,JFN1,JFN2,NCY
      REAL(r8)    :: fil_lat
      REAL(r8)    :: X (IMT,JMT,max_blocks_clinic),XS (IMT),Z (IMT,JMT,KM,max_blocks_clinic)
      INTEGER    :: NN(JMT), MAX_NN, iblock


!lhl      fil_lat=54.D0
!
         
      MAX_NN = 0
      do j =jst,jmt
         if (sint(j).le.cos(fil_lat*torad)) then
            NN(j) = int(cos(fil_lat*torad)/sint(j)*1.2D0)
            if (NN(j) .gt. MAX_NN) MAX_NN = NN(J)
         else
            NN(j) = 0
         endif
      enddo

!!$OMP PARALLEL DO PRIVATE (iblock,J,nn,xs)
      DO NCY = 1,MAX_NN
         do j =jst,jmt
            if (NN(j) .ge. NCY) then
               DO I = 1,IMT
                  XS (I)= X (I,J,IBLOCK)* Z (I,J,1,iblock)
               END DO
               DO I = 2,IMM
                  X(I,J,iblock)=(XS(I)*(1.0D0-0.25D0*Z(I-1,J,1,iblock)-0.25D0*Z(I+1,J,1,iblock)) &
                            +0.25D0*(XS(I-1)+XS(I+1)))*Z(I,J,1,iblock)
               END DO
!              if (nx_proc == 1) then
!                 X (1,J)= X (IMM,J)
!                 X (IMT,J)= X (2,J)
!              endif
            endif
         enddo
!        if (nx_proc .ne. 1) then
!           call exchange_2D_boundary(x,1,NN,NCY,0)
!           call exchange_2D_boundary(x,1,NN,NCY,1)
!        end if
      END DO


   RETURN
   END SUBROUTINE SMZ0

