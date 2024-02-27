*********************************************************
!    A large deformation umat for 
!      compressible NeoHookean material 
      
!    *Orientation is open, all the tensors are 
!       transformed to corotational coordinates  
      
!                   Sifan Yin
!                   01/04/2021
********************************************************

      
      

      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
     
      include 'aba_param.inc'

      DIMENSION stress(ntens),statev(nstatv),
     + ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     + stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     + props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      CHARACTER*80 cmname,file1
      CHARACTER*256 jobName,outDir,fileName
      
      if (CMNAME(1:4).eq.'ENDO') then
   
          call UMAT_ENDO(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
          
       elseif(CMNAME(1:4).eq.'MESE' )then
          call UMAT_MESE(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
      endif
      return 
      end subroutine UMAT
        
**********************************************************     
**********************************************************
      subroutine UMAT_ENDO(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
     
      include 'aba_param.inc'
      
      DIMENSION stress(ntens),statev(nstatv),
     + ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     + stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     + props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      CHARACTER*80 cmname,file1

      integer i,j,k,l

      real*8 Iden(3,3),Fal(3,3),theta_t,theta_tau,Tal(3,3)
      real*8 mu,Kbulk,effStr,detF,trB
      real*8 Baldis(3,3),trBdis,Baldis0(3,3),Bal(3,3)
      
      real*8 R(3,3),U(3,3),RT(3,3),Aal(3,3)
      real*8 Ginv(3,3),Ginval(3,3),SpTanModal(3,3,3,3)
      real*8 g1,g2,g3,detG,detA
     
      ! Parameters
      real*8 zero,one,two,three,four,half,third,eight,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0, three=3.d0,
     1      four=4.d0,third=1.d0/3.d0,eight=8.d0,nine=9.d0)
      
      ! Identity matrix
      call onem(Iden)
      
      ! Obtain old and new deformation gradients
      
      Fal = DFGRD1

      mu = props(1)
      Kbulk = props(2)
      g1 = one
      g2 = one-time(2)*time(2)+two*time(2)
      g3 = one+two*time(2)*time(2)
      
 
      
      call mdet(Fal,detF)
      call kpolarDecomp(Fal,U,R)
      RT = transpose(R)
      Ginv = zero
      Ginv(1,1) = one/g1
      Ginv(2,2) = one/g2
      Ginv(3,3) = one/g3
      detG = g1*g2*g3
      detA = detF/detG
      Ginval = matmul(RT,matmul(Ginv,R))
      Aal = matmul(Fal,Ginval)
      Bal = matmul(Aal,transpose(Aal))   
      trB = Bal(1,1)+Bal(2,2)+Bal(3,3)
      Baldis = detA**(-two/three)*Bal
      trBdis = Baldis(1,1)+Baldis(2,2)+Baldis(3,3)
      Baldis0 = Baldis - third*trBdis*Iden
      
      ! Compute the Cauchy stress
      !TEE = mu/detA*BEEdis0 + Kbulk*(one-one/detA)*Iden
      Tal = mu/detA*Baldis0 + Kbulk*(detA-one)*Iden
           
      
      ! Compute the material jacobian
      !call TangentModulus(Baldis,detA,props,SpTanMod)

      SpTanModal = zero
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      SpTanModal(i,j,k,l) = SpTanModal(i,j,k,l)
     1                   + mu/detA*
     1       (
     1       half*Iden(i,k)*Baldis(j,l)
     1       +half*Iden(j,k)*Baldis(i,l)
     1       +half*Iden(i,l)*Baldis(j,k)
     1       +half*Iden(j,l)*Baldis(i,k)
     1       -(two/three)*Iden(i,j)*Baldis(k,l)
     1       -(two/three)*Iden(k,l)*Baldis(i,j)
     1       +(two/nine)*trBdis*Iden(i,j)*Iden(k,l)
     1        )
     1       + Kbulk*(two*detA-one)*Iden(i,j)*Iden(k,l)
                  end do
              end do
          end do
      end do
      
      if (ntens.eq.6) then
          stress(1) = Tal(1,1)
          stress(2) = Tal(2,2)
          stress(3) = Tal(3,3)
          stress(4) = Tal(1,2)
          stress(5) = Tal(1,3)
          stress(6) = Tal(2,3)
      elseif(ntens.eq.4) then
          stress(1) = Tal(1,1)
          stress(2) = Tal(2,2)
          stress(3) = Tal(3,3)
          stress(4) = Tal(1,2)
      endif
      if (ntens.eq.6) then
          call jac3D(SpTanModal,ddsdde)
      elseif(ntens.eq.4) then
          call jac2D(SpTanModal,ddsdde)
      endif
      
      return
      end subroutine UMAT_ENDO
      
      INCLUDE 'kpolarDecomp.for'
*****************************************************      
******************************************************
          subroutine UMAT_mese(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
     
      include 'aba_param.inc'
      
      DIMENSION stress(ntens),statev(nstatv),
     + ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     + stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     + props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      CHARACTER*80 cmname,file1

      integer i,j,k,l

      real*8 Iden(3,3),Fal(3,3),Tal(3,3)
      real*8 mu,Kbulk,effStr,detF,trB
      real*8 Baldis(3,3),trBdis,Baldis0(3,3),Bal(3,3)
      
      real*8 R(3,3),U(3,3),RT(3,3),Aal(3,3)
      real*8 Ginv(3,3),Ginval(3,3),SpTanModal(3,3,3,3)
      real*8 g1,g2,g3,detG,detA
     
      ! Parameters
      real*8 zero,one,two,three,four,half,third,eight,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0, three=3.d0,
     1      four=4.d0,third=1.d0/3.d0,eight=8.d0,nine=9.d0)
      
      ! Identity matrix
      call onem(Iden)
      
      ! Obtain old and new deformation gradients
      
      Fal = DFGRD1

      mu = props(1)
      Kbulk = props(2)
      g1 = one 
      g2 = one 
      g3 = one
 
      
     
     
      call mdet(Fal,detF)
      call kpolarDecomp(Fal,U,R)
      RT = transpose(R)
      Ginv = zero
      Ginv(1,1) = one/g1
      Ginv(2,2) = one/g2
      Ginv(3,3) = one/g3
      detG = g1*g2*g3
      detA = detF/detG
      Ginval = matmul(RT,matmul(Ginv,R))
      Aal = matmul(Fal,Ginval)
      Bal = matmul(Aal,transpose(Aal))   
      trB = Bal(1,1)+Bal(2,2)+Bal(3,3)
      Baldis = detA**(-two/three)*Bal
      trBdis = Baldis(1,1)+Baldis(2,2)+Baldis(3,3)
      Baldis0 = Baldis - third*trBdis*Iden
      
      ! Compute the Cauchy stress
      !TEE = mu/detA*BEEdis0 + Kbulk*(one-one/detA)*Iden
      Tal = mu/detA*Baldis0 + Kbulk*(one-one/detA)*Iden
           
      
      ! Compute the material jacobian
      !call TangentModulus(Baldis,detA,props,SpTanMod)

      SpTanModal = zero
      do i=1,3
          do j=1,3
              do k=1,3
                  do l=1,3
                      SpTanModal(i,j,k,l) = SpTanModal(i,j,k,l)
     1                   + mu/detA*
     1       (
     1       half*Iden(i,k)*Baldis(j,l)
     1       +half*Iden(j,k)*Baldis(i,l)
     1       +half*Iden(i,l)*Baldis(j,k)
     1       +half*Iden(j,l)*Baldis(i,k)
     1       -(two/three)*Iden(i,j)*Baldis(k,l)
     1       -(two/three)*Iden(k,l)*Baldis(i,j)
     1       +(two/nine)*trBdis*Iden(i,j)*Iden(k,l)
     1        )
     1       + Kbulk*Iden(i,j)*Iden(k,l)
                  end do
              end do
          end do
      end do
      
      if (ntens.eq.6) then
          stress(1) = Tal(1,1)
          stress(2) = Tal(2,2)
          stress(3) = Tal(3,3)
          stress(4) = Tal(1,2)
          stress(5) = Tal(1,3)
          stress(6) = Tal(2,3)
      elseif(ntens.eq.4) then
          stress(1) = Tal(1,1)
          stress(2) = Tal(2,2)
          stress(3) = Tal(3,3)
          stress(4) = Tal(1,2)
      endif
      if (ntens.eq.6) then
          call jac3D(SpTanModal,ddsdde)
      elseif(ntens.eq.4) then
          call jac2D(SpTanModal,ddsdde)
      endif
      
      return
      end subroutine UMAT_mese 
      
      
      
      
******************************************************    
******************************************************     
   
C            
      	SUBROUTINE MDET(A,DET)
 
C 	THIS SUBROUTINE CALCULATES THE DETERMINANT
C 	OF A 3 BY 3 MATRIX [A].


	REAL*8  A(3,3), DET

	DET =	  A(1,1)*A(2,2)*A(3,3) 
     +	        + A(1,2)*A(2,3)*A(3,1)
     +	        + A(1,3)*A(2,1)*A(3,2)
     +		- A(3,1)*A(2,2)*A(1,3)
     +		- A(3,2)*A(2,3)*A(1,1)
     +		- A(3,3)*A(2,1)*A(1,2)

	RETURN
	END subroutine mdet
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
      
      
      subroutine jac2D(SpTanMod,ddsdde)

      real*8 SpTanMod(3,3,3,3),ddsdde(4,4)

      ddsdde(1,1) = SpTanMod(1,1,1,1)
      ddsdde(1,2) = SpTanMod(1,1,2,2)
      ddsdde(1,3) = SpTanMod(1,1,3,3)
      ddsdde(1,4) = SpTanMod(1,1,1,2)

      ddsdde(2,1) = SpTanMod(2,2,1,1)
      ddsdde(2,2) = SpTanMod(2,2,2,2)
      ddsdde(2,3) = SpTanMod(2,2,3,3)
      ddsdde(2,4) = SpTanMod(2,2,1,2)

      ddsdde(3,1) = SpTanMod(3,3,1,1)
      ddsdde(3,2) = SpTanMod(3,3,2,2)
      ddsdde(3,3) = SpTanMod(3,3,3,3)
      ddsdde(3,4) = SpTanMod(3,3,1,2)

      ddsdde(4,1) = SpTanMod(1,2,1,1)
      ddsdde(4,2) = SpTanMod(1,2,2,2)
      ddsdde(4,3) = SpTanMod(1,2,3,3)
      ddsdde(4,4) = SpTanMod(1,2,1,2)

      end subroutine jac2D

***********************************************************************

      subroutine jac3D(SpTanMod,ddsdde)

      real*8 SpTanMod(3,3,3,3),ddsdde(6,6)

      ddsdde(1,1) = SpTanMod(1,1,1,1)
      ddsdde(1,2) = SpTanMod(1,1,2,2)
      ddsdde(1,3) = SpTanMod(1,1,3,3)
      ddsdde(1,4) = SpTanMod(1,1,1,2)
      ddsdde(1,5) = SpTanMod(1,1,1,3)
      ddsdde(1,6) = SpTanmod(1,1,2,3)

      ddsdde(2,1) = SpTanMod(2,2,1,1)
      ddsdde(2,2) = SpTanMod(2,2,2,2)
      ddsdde(2,3) = SpTanMod(2,2,3,3)
      ddsdde(2,4) = SpTanMod(2,2,1,2)
      ddsdde(2,5) = SpTanMod(2,2,1,3)
      ddsdde(2,6) = SpTanmod(2,2,2,3)

      ddsdde(3,1) = SpTanMod(3,3,1,1)
      ddsdde(3,2) = SpTanMod(3,3,2,2)
      ddsdde(3,3) = SpTanMod(3,3,3,3)
      ddsdde(3,4) = SpTanMod(3,3,1,2)
      ddsdde(3,5) = SpTanMod(3,3,1,3)
      ddsdde(3,6) = SpTanmod(3,3,2,3)

      ddsdde(4,1) = SpTanMod(1,2,1,1)
      ddsdde(4,2) = SpTanMod(1,2,2,2)
      ddsdde(4,3) = SpTanMod(1,2,3,3)
      ddsdde(4,4) = SpTanMod(1,2,1,2)
      ddsdde(4,5) = SpTanMod(1,2,1,3)
      ddsdde(4,6) = SpTanmod(1,2,2,3)

      ddsdde(5,1) = SpTanMod(1,3,1,1)
      ddsdde(5,2) = SpTanMod(1,3,2,2)
      ddsdde(5,3) = SpTanMod(1,3,3,3)
      ddsdde(5,4) = SpTanMod(1,3,1,2)
      ddsdde(5,5) = SpTanMod(1,3,1,3)
      ddsdde(5,6) = SpTanmod(1,3,2,3)

      ddsdde(6,1) = SpTanMod(2,3,1,1)
      ddsdde(6,2) = SpTanMod(2,3,2,2)
      ddsdde(6,3) = SpTanMod(2,3,3,3)
      ddsdde(6,4) = SpTanMod(2,3,1,2)
      ddsdde(6,5) = SpTanMod(2,3,1,3)
      ddsdde(6,6) = SpTanmod(2,3,2,3)

      end subroutine jac3D     
      
      
*************************************************
************************************************
     