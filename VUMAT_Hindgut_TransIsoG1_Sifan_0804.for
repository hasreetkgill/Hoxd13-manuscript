************************************************************************
!                                        ABAQUS VUMAT
!              Hindgut Growth--Endoderm and Mesenchyme
!     Mesenchyme--  Purely elastic Neo-Hookean material with no growth 
!     Endoderm  --  Area growth coupled to Neo-Hookean material
!                 VUMAT Hindgut_TransverseIso_Growth_Sifan_0804
!                                           Sifan Yin 
!                                        08/04/2022
************************************************************************
!     State Variables
!     --------------------------------------------------------------
!      stateNew(km,1) = growth factor  -- thetag_tau
!
!     Material Properties Vector
!     --------------------------------------------------------------
!       mu_m  = props(1)  --- Shear modulus for mesenchyme
!       lam_m = props(2)  --- Bulk modulus for mesenchyme
!       mu_e   = props(3) --- Shear modulus for endoderm
!       lam_e  = props(4)  --- Bulk modulus for endoderm
!
!**********************************************************************
      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C
      if (cmname(1:4) .eq. 'MESE') then
            call VUMAT_mese(
     +     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
     +     stressNew, stateNew, enerInternNew, enerInelasNew)
      else if (cmname(1:4) .eq. 'ENDO') then
            call vumat_endo(
     +     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
     +     stressNew, stateNew, enerInternNew, enerInelasNew)
      end if
      

      return
      end subroutine vumat
      
    
*******************************************************************************************
!             Endoderm -- Area growht coupled to elastic Neo-Hookean material             !
*******************************************************************************************  
      
      subroutine vumat_endo(
     +     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
     +     stressNew, stateNew, enerInternNew, enerInelasNew)

      include 'vaba_param.inc'

      dimension props(nprops), density(nblock), coordMp(nblock,*),
     +     charLength(nblock), strainInc(nblock,ndir+nshr),
     +     relSpinInc(nblock,nshr), tempOld(nblock),
     +     stretchOld(nblock,ndir+nshr),
     +     defgradOld(nblock,ndir+nshr+nshr),
     +     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     +     stateOld(nblock,nstatev), enerInternOld(nblock),
     +     enerInelasOld(nblock), tempNew(nblock),
     +     stretchNew(nblock,ndir+nshr),
     +     defgradNew(nblock,ndir+nshr+nshr),
     +     fieldNew(nblock,nfieldv),
     +     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     +     enerInternNew(nblock), enerInelasNew(nblock)

      character*80 cmname

      integer i,j,l,i1,j1,ii,jj,kk,ll,km,ifail

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),U_t(3,3),U_tau(3,3),Fp_t(3,3)
      real*8 Fp_tau(3,3),Me_t(3,3),Me_tau(3,3),nuP_t,nuP_tau,S_t,S_tau
      real*8 gBarP_t,gBarP_tau,T_tau(3,3),R_tau(3,3),U_inv(3,3),detF
      real*8 Fp_inv(3,3),Ee_tau(3,3),Re_tau(3,3),Ue_tau(3,3)
      real*8 pnu0,damage_t,damage_tau,mag_Dp_tau,pwrinct,stress_power
      real*8 Be_tau(3,3), Fg_tau(3,3), Fe_tau(3,3), Je, Jg, Fginv(3,3)
      ! Material Props
      real*8 mu_e, lambda_e, growth_tau
      real*8 coordx, coordy, coordz
      real*8 a0(3,1), tmp
      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0)


      ! Identity matrix for later use.
      !
      call onem(Iden)

      ! Read material properties needed here
      !
      mu_e = props(3)
      lambda_e  = props(4)
      !
      ! START LOOP OVER MATERIAL POINTS:
      !
      do km=1,nblock

           
         ! Copy old and new deformation gradients
         !
         F_t(1,1) = defgradOld(km,1)
         F_t(2,2) = defgradOld(km,2)
         F_t(3,3) = defgradOld(km,3)
         F_t(1,2) = defgradOld(km,4)
         F_tau(1,1) = defgradNew(km,1)
         F_tau(2,2) = defgradNew(km,2)
         F_tau(3,3) = defgradNew(km,3)
         F_tau(1,2) = defgradNew(km,4)
         U_tau(1,1) = stretchNew(km,1)
         U_tau(2,2) = stretchNew(km,2)
         U_tau(3,3) = stretchNew(km,3)
         U_tau(1,2) = stretchNew(km,4)
         if(nshr .lt. 2) then
            ! 2D case
            F_t(2,1) = defgradOld(km,5)
            F_t(1,3) = zero
            F_t(2,3) = zero
            F_t(3,1) = zero
            F_t(3,2) = zero
            F_tau(2,1) = defgradNew(km,5)
            F_tau(1,3) = zero
            F_tau(2,3) = zero
            F_tau(3,1) = zero
            F_tau(3,2) = zero
            U_tau(2,1) = U_tau(1,2)
            U_tau(1,3) = zero
            U_tau(2,3) = zero
            U_tau(3,1) = zero
            U_tau(3,2) = zero
         else
            ! 3D case
            F_t(2,3) = defgradOld(km,5)
            F_t(3,1) = defgradOld(km,6)
            F_t(2,1) = defgradOld(km,7)
            F_t(3,2) = defgradOld(km,8)
            F_t(1,3) = defgradOld(km,9)
            F_tau(2,3) = defgradNew(km,5)
            F_tau(3,1) = defgradNew(km,6)
            F_tau(2,1) = defgradNew(km,7)
            F_tau(3,2) = defgradNew(km,8)
            F_tau(1,3) = defgradNew(km,9)
            U_tau(2,3) = stretchNew(km,5)
            U_tau(3,1) = stretchNew(km,6)
            U_tau(2,1) = U_tau(1,2)
            U_tau(3,2) = U_tau(2,3)
            U_tau(1,3) = U_tau(3,1)
         end if


         ! Compute the relative volume change
         !
             call mdet(F_tau,detF)
!******************************************************
!              if (totalTime.eq.zero) then
!                    stateOld(km,1)  = one
!              endif
!********************************************************              
              
!             growth_t = stateOld(km,1) ! Growth parameter at time t 
              
              
             growth_tau = one + TotalTime         
             coordx = coordMp(km,1)
             coordy = coordMp(km,2)
             coordz = coordMp(km,3)
             
             

!      if (totalTime.eq.zero) then
!          thetag_tau = thetag_t
!      else
!          thetag_tau = thetag_t+dt*thetag_dot
!      endif

         stateNew(km,1) = growth_tau         
          !a0(1,1) = 0.0
          !a0(2,1) = 0.0
          !a0(3,1) = 1.0
         ! Set outer unit normal vector  
             a0(1,1) = -coordx
             a0(2,1) = -coordy
             a0(3,1) = 0.0
             tmp = dsqrt(a0(1,1)**2.0 + a0(2,1)**2.0 + a0(3,1)**2.0)
             a0 = a0/tmp
          
      
          Fg_tau = dsqrt(growth_tau)*Iden
     +        +(1.0 - dsqrt(growth_tau))*matmul(a0,transpose(a0)) 
         
          call m3inv(Fg_tau,Fginv)
         
          Fe_tau = matmul(F_tau, Fginv)
          Be_tau = matmul(Fe_tau,transpose(Fe_tau))
          call mdet(Fg_tau,Jg)
          Je = detF/Jg
          T_tau = ((lambda_e*dlog(Je) - mu_e)*Iden  + mu_e*Be_tau)/Je
         
          
!           open(unit=502, file='F:\Abaqus_VUMAT\OneEle_Initial_0802.txt')
 !          write(502,*)"________________________________"
 !          write(502,"(f16.8)" ) TotalTime
 !          write(502,*)"________________________________"      
     
C      open(unit=501, file='F:\Abaqus_VUMAT\OneEle_thetag0802.txt')
C     write(501,*)"________________________________"
C     write(501,"(f16.8,f16.8)" ) TotalTime,thetag_tau
C      write(501,"(f16.8,f16.8,f16.8)") Fg_tau(1,1),Fg_tau(1,2),Fg_tau(1,3)
C     write(501,"(f16.8,f16.8,f16.8)") Fg_tau(2,1),Fg_tau(2,2),Fg_tau(2,3)
C     write(501,"(f16.8,f16.8,f16.8)") Fg_tau(3,1),Fg_tau(3,2),Fg_tau(3,3)
C      write(501,*)"________________________________"

         ! Update state variables
         
         ! ABAQUS/Explicit uses stress measure (transpose(R) T R)
         !
         call m3inv(U_tau,U_inv)
         R_tau = matmul(F_tau,U_inv)
         T_tau = matmul(transpose(R_tau),matmul(T_tau,R_tau))

         do i=1,ndir
            stressNew(km,i) = T_tau(i,i)
         end do
         if(nshr.ne.0) then
            stressNew(km,ndir+1) = T_tau(1,2)
            if(nshr.ne.1) then
               stressNew(km, ndir+2) = T_tau(2,3)
               if(nshr.ne.2) then
                  stressNew(km,ndir+3) = T_tau(1,3)
               endif
            endif
         endif


         ! Update the specific internal energy
         !
         stress_power = 0.d0
         do i = 1,ndir
            stress_power = stress_power +
     +           0.5*((StressOld(km,i)+StressNew(km,i))*
     +           StrainInc(km,i))
         enddo
         
         select case (nshr)
         case(1)
            stress_power = stress_power + 
     +           0.5*((StressOld(km,ndir+1)+StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1))
         case(3)
            stress_power = stress_power + 
     +           0.5*(((StressOld(km,ndir+1) + StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1)) +
     +           ((StressOld(km,ndir+2)+ StressNew(km,ndir+2)) *
     +           StrainInc(km,ndir+2))+
     +           ((StressOld(km,ndir+3) + StressNew(km,ndir+3))*
     +           StrainInc(km,ndir+3)))
         end select
           
         enerInternNew(km) = enerInternOld(km) + 
     +        stress_power/density(km)
           
         enerInelasNew(km) = enerInelasOld(km) + 
     +        pwrinct/density(km)
           
           
      enddo ! end loop over material points

      end subroutine vumat_endo
      
      
      
!**************************************************************************!
!**************************************************************************!
!**************************************************************************!
!             white matter -- Purely elastic Neo-Hookean material                !
!**************************************************************************!
!**************************************************************************!
      subroutine vumat_mese(
     +     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
     +     stressNew, stateNew, enerInternNew, enerInelasNew)

      include 'vaba_param.inc'

      dimension props(nprops), density(nblock), coordMp(nblock,*),
     +     charLength(nblock), strainInc(nblock,ndir+nshr),
     +     relSpinInc(nblock,nshr), tempOld(nblock),
     +     stretchOld(nblock,ndir+nshr),
     +     defgradOld(nblock,ndir+nshr+nshr),
     +     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     +     stateOld(nblock,nstatev), enerInternOld(nblock),
     +     enerInelasOld(nblock), tempNew(nblock),
     +     stretchNew(nblock,ndir+nshr),
     +     defgradNew(nblock,ndir+nshr+nshr),
     +     fieldNew(nblock,nfieldv),
     +     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     +     enerInternNew(nblock), enerInelasNew(nblock)

      character*80 cmname

      integer i,j,l,i1,j1,ii,jj,kk,ll,km,ifail

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),U_t(3,3),U_tau(3,3),Fp_t(3,3)
      real*8 Fp_tau(3,3),Me_t(3,3),Me_tau(3,3),nuP_t,nuP_tau,S_t,S_tau
      real*8 gBarP_t,gBarP_tau,T_tau(3,3),R_tau(3,3),U_inv(3,3),detF
      real*8 Fp_inv(3,3),Ee_tau(3,3),Re_tau(3,3),Ue_tau(3,3),Fe_tau(3,3)
      real*8 pnu0,damage_t,damage_tau,mag_Dp_tau,pwrinct,stress_power
      real*8 nu1,nu3,nu5,effStr,Bdis(3,3),Bdis0(3,3),trBdis
      real*8 B_tau(3,3)
      ! Material Props
      !       mu_m  = props(1)  --- Shear modulus for mesenchyme
      !       lambda_m = props(2)  --- Bulk modulus for mesenchyme
      
      
      real*8 mu_m,lambda_m

      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0)


      ! Identity matrix for later use.
      !
      call onem(Iden)

      ! Read material properties needed here
      !
      mu_m = props(1)
      lambda_m = props(2)

      !
      ! START LOOP OVER MATERIAL POINTS:
      !
      do km=1,nblock

           
         ! Copy old and new deformation gradients
         !
         F_t(1,1) = defgradOld(km,1)
         F_t(2,2) = defgradOld(km,2)
         F_t(3,3) = defgradOld(km,3)
         F_t(1,2) = defgradOld(km,4)
         F_tau(1,1) = defgradNew(km,1)
         F_tau(2,2) = defgradNew(km,2)
         F_tau(3,3) = defgradNew(km,3)
         F_tau(1,2) = defgradNew(km,4)
         U_tau(1,1) = stretchNew(km,1)
         U_tau(2,2) = stretchNew(km,2)
         U_tau(3,3) = stretchNew(km,3)
         U_tau(1,2) = stretchNew(km,4)
         if(nshr .lt. 2) then
            ! 2D case
            F_t(2,1) = defgradOld(km,5)
            F_t(1,3) = zero
            F_t(2,3) = zero
            F_t(3,1) = zero
            F_t(3,2) = zero
            F_tau(2,1) = defgradNew(km,5)
            F_tau(1,3) = zero
            F_tau(2,3) = zero
            F_tau(3,1) = zero
            F_tau(3,2) = zero
            U_tau(2,1) = U_tau(1,2)
            U_tau(1,3) = zero
            U_tau(2,3) = zero
            U_tau(3,1) = zero
            U_tau(3,2) = zero
         else
            ! 3D case
            F_t(2,3) = defgradOld(km,5)
            F_t(3,1) = defgradOld(km,6)
            F_t(2,1) = defgradOld(km,7)
            F_t(3,2) = defgradOld(km,8)
            F_t(1,3) = defgradOld(km,9)
            F_tau(2,3) = defgradNew(km,5)
            F_tau(3,1) = defgradNew(km,6)
            F_tau(2,1) = defgradNew(km,7)
            F_tau(3,2) = defgradNew(km,8)
            F_tau(1,3) = defgradNew(km,9)
            U_tau(2,3) = stretchNew(km,5)
            U_tau(3,1) = stretchNew(km,6)
            U_tau(2,1) = U_tau(1,2)
            U_tau(3,2) = U_tau(2,3)
            U_tau(1,3) = U_tau(3,1)
         end if


          call mdet(F_tau,detF)
          B_tau = matmul(F_tau, transpose(F_tau))
          T_tau = ((lambda_m*dlog(detF) - mu_m)*Iden  + mu_m*B_tau)/detF
          call m3inv(U_tau,U_inv)
                 R_tau = matmul(F_tau,U_inv)
                 T_tau = matmul(transpose(R_tau),matmul(T_tau,R_tau))

         do i=1,ndir
               stressNew(km,i) = T_tau(i,i)
         end do
         if(nshr.ne.0) then
              stressNew(km,ndir+1) = T_tau(1,2)
               if(nshr.ne.1) then
                  stressNew(km, ndir+2) = T_tau(2,3)
               if(nshr.ne.2) then
                   stressNew(km,ndir+3) = T_tau(1,3)
               endif
            endif
         endif


         ! Update the specific internal energy
         !
         stress_power = 0.d0
         do i = 1,ndir
            stress_power = stress_power +
     +           0.5*((StressOld(km,i)+StressNew(km,i))*
     +           StrainInc(km,i))
         enddo
         
         select case (nshr)
         case(1)
            stress_power = stress_power + 
     +           0.5*((StressOld(km,ndir+1)+StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1))
         case(3)
            stress_power = stress_power + 
     +           0.5*(((StressOld(km,ndir+1) + StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1)) +
     +           ((StressOld(km,ndir+2)+ StressNew(km,ndir+2)) *
     +           StrainInc(km,ndir+2))+
     +           ((StressOld(km,ndir+3) + StressNew(km,ndir+3))*
     +           StrainInc(km,ndir+3)))
         end select
           
         enerInternNew(km) = enerInternOld(km) + 
     +        stress_power/density(km)
           
         enerInelasNew(km) = enerInelasOld(km) + 
     +        pwrinct/density(km)
           
           
      enddo ! end loop over material points

      end subroutine vumat_mese

      
***********************************************************************
c
c
c  The following are all utility routines used in fortran codes
c
c
c
C**********************************************************************
	SUBROUTINE ONEM(A)

C	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C	3 BY 3 MATRIX [A]
C**********************************************************************

        REAL*8 A(3,3)
        DATA ZERO/0.D0/
        DATA ONE/1.D0/

	DO 1 I=1,3
	  DO 1 J=1,3
	    IF (I .EQ. J) THEN
              A(I,J) = 1.0
            ELSE
              A(I,J) = 0.0
            ENDIF
1       CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE MTRANS(A,ATRANS)
 
C	THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 
C	MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C**********************************************************************

	REAL*8 A(3,3),ATRANS(3,3)

	DO 1 I=1,3
 	  DO 1 J=1,3
	    ATRANS(J,I) = A(I,J)
1	CONTINUE

	RETURN
	END


C**********************************************************************
	SUBROUTINE MDET(A,DET)
 
C 	THIS SUBROUTINE CALCULATES THE DETERMINANT
C 	OF A 3 BY 3 MATRIX [A].
C**********************************************************************

	REAL*8  A(3,3), DET

	DET =	  A(1,1)*A(2,2)*A(3,3) 
     +	        + A(1,2)*A(2,3)*A(3,1)
     +	        + A(1,3)*A(2,1)*A(3,2)
     +		- A(3,1)*A(2,2)*A(1,3)
     +		- A(3,2)*A(2,3)*A(1,1)
     +		- A(3,3)*A(2,1)*A(1,2)

	RETURN
	END

C**********************************************************************
	SUBROUTINE M3INV(A,AINV)

C 	THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX
C	[A] AND PLACES THE RESULT IN [AINV]. 
C 	IF DET(A) IS ZERO, THE CALCULATION
C 	IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
C**********************************************************************

	REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)

C	A(3,3)	        -- THE MATRIX WHOSE INVERSE IS DESIRED.
C	DET		-- THE COMPUTED DETERMINANT OF [A].
C	ACOFAC(3,3)	-- THE MATRIX OF COFACTORS OF A(I,J).
C			   THE SIGNED MINOR (-1)**(I+J)*M_IJ
C			   IS CALLED THE COFACTOR OF A(I,J).
C	AADJ(3,3)	-- THE ADJOINT OF [A]. IT IS THE MATRIX
C			   OBTAINED BY REPLACING EACH ELEMENT OF
C			   [A] BY ITS COFACTOR, AND THEN TAKING
C			   TRANSPOSE OF THE RESULTING MATRIX.
C	AINV(3,3)	-- RETURNED AS INVERSE OF [A].
C			   [AINV] = [AADJ]/DET.
C----------------------------------------------------------------------

	CALL MDET(A,DET)
	IF ( DET .EQ. 0.D0 ) THEN
	  write(*,10)
	  STOP
	ENDIF
	CALL MCOFAC(A,ACOFAC)
	CALL MTRANS(ACOFAC,AADJ)
	DO 1 I = 1,3
	DO 1 J = 1,3
	     AINV(I,J) = AADJ(I,J)/DET
1	CONTINUE
10	FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
     +         10X,'PROGRAM TERMINATED')

	RETURN
	END

C**********************************************************************
	SUBROUTINE MCOFAC(A,ACOFAC)
 
C 	THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C 	AND PLACES THE RESULT IN [ACOFAC]. 
C**********************************************************************

	REAL*8  A(3,3), ACOFAC(3,3)

	ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
	ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
	ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
	ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
	ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
	ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
	ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
	ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
	ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

	RETURN
	END

