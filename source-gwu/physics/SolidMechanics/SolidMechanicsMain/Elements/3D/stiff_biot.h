!! IO: (ke,qs,XYZ,YoungsModulus,PoissonRatio,Qie)
    
#if SM_ELEM_BUILDTYPE == EIGHT_NODE_HEXAHEDRON
    USE quadrule_C3_GL4
    USE sm_element_interface, only: el05_ShapeFunc, ShapeFunc_Expand
#elif SM_ELEM_BUILDTYPE == TWSEVEN_NODE_HEXAHEDRON
    USE quadrule_C3_GL5
    USE sm_element_interface, only: el12_ShapeFunc, ShapeFunc_Expand
#endif  
    USE sm_Misc_interface, only: compute_Jacobian3d
    
    implicit none
#include "Flash.h"
#if SM_ELEM_BUILDTYPE == EIGHT_NODE_HEXAHEDRON
    integer, parameter :: nen_e = 8
#elif SM_ELEM_BUILDTYPE == TWSEVEN_NODE_HEXAHEDRON
    integer, parameter :: nen_e = 27
#endif  
	
    integer, parameter :: nee = NDIM*nen_e
    real, external :: DDOT
    
    real, dimension(nee,nee), intent(out) :: ke
    real, dimension(nee), intent(out)     :: Qs
	
    real, dimension(nen_e,NDIM), intent(in)   :: XYZ, Qie
    real, intent(in) :: YoungsModulus,PoissonRatio
    
    integer :: a, i1, i2
    real    :: xi, eta, zeta, detJ, dV0, mu, lambda
    real, dimension(3,3) :: Jinv
    real, dimension(nen_e)   :: NN, NNxi, NNeta, NNzeta, NNx, NNy, NNz
    real, dimension(nen_e,nen_e) :: dXX, dYY, dZZ, dYZ, dXZ, dXY, dYZ2, dXZ2, dXY2, tempD
    real, dimension(nen_e,3) :: DN
    real, dimension(6,nee)   :: DS
    real, dimension(nee,6)   :: Bt
    real, dimension(3,3)     :: S, G, L, F, R, U, Yy, Uinv, dF, Theta, tempB
    real                     :: detYy, TrU, TrL, TrTheta
    real, dimension(6)       :: Svec
    real, dimension(nen_e,3) :: XYZQe, tempA
    
    ! SVD parameters and such:
    integer :: INFO
    integer, parameter :: LWORK = 201
    real :: Sigma(3), Psi(3,3), LambdaT(3,3), WORK(LWORK)
    
    ! Init ke and qs
    ke = 0.
    qs = 0.
    
    ! get Xe and qe:
    !call get_Nodal_XYZ(e, nen_e, XYZ)
    !call get_Nodal_UVW_qi(e, nen_e, Qie)  
    XYZQe  = XYZ+Qie
        
    ! Get the material-constant matrix [C].
    mu = YoungsModulus/(2.*(1.+PoissonRatio))
    lambda = YoungsModulus*PoissonRatio/(1.+PoissonRatio)/(1.-2.*PoissonRatio)   

    ! loop over the integral
    do i1 = 1,quad_nt
        xi = quad_xi(i1)
        eta = quad_eta(i1)
        zeta = quad_zeta(i1)
    
        ! Get Shape Functions, and local derivatives
#if SM_ELEM_BUILDTYPE == EIGHT_NODE_HEXAHEDRON
        call el05_ShapeFunc(NN,NNxi,NNeta,NNzeta,xi,eta,zeta)
#elif SM_ELEM_BUILDTYPE == TWSEVEN_NODE_HEXAHEDRON
        call el12_ShapeFunc(NN,NNxi,NNeta,NNzeta,xi,eta,zeta)
#endif    
    
        ! Compute Jacobian
        !J = [Nxi*Xe, Neta*Xe, Nzeta*Xe]
        !detJ = det(J)
        call compute_Jacobian3d(nen_e, XYZ, NNxi, NNeta, NNzeta, Jinv, detJ)
    
        ! Derivative matrix DN
        !DN = [NNxi', NNeta', NNzeta']/J
        tempA(1:nen_e,1) = NNxi
        tempA(1:nen_e,2) = NNeta
        tempA(1:nen_e,3) = NNzeta
        call DGEMM('N','N', nen_e, 3, 3, 1./detJ, tempA, nen_e, Jinv, 3, 0.d0, DN, nen_e)                
    
        ! Build Nx, Ny, Nz
        NNx = DN(1:nen_e,1)
        NNy = DN(1:nen_e,2)
        NNz = DN(1:nen_e,3)
        
        ! Build dXX, dYY, dZZ, dYZ, dXZ, dXY
        !dXX = NNx'*NNx
        !dYY = NNy'*NNy
        !dZZ = NNz'*NNz
        !dYZ = NNy'*NNz
        !dXZ = NNx'*NNz
        !dXY = NNx'*NNy
        dXX = 0.
        call DGER( nen_e, nen_e, 1., NNx, 1, NNx, 1, dXX, nen_e)
        dYY = 0.
        call DGER( nen_e, nen_e, 1., NNy, 1, NNy, 1, dYY, nen_e)
        dZZ = 0.
        call DGER( nen_e, nen_e, 1., NNz, 1, NNz, 1, dZZ, nen_e)
        dYZ = 0.
        call DGER( nen_e, nen_e, 1., NNy, 1, NNz, 1, dYZ, nen_e)
        dXZ = 0.
        call DGER( nen_e, nen_e, 1., NNx, 1, NNz, 1, dXZ, nen_e)
        dXY = 0.
        call DGER( nen_e, nen_e, 1., NNx, 1, NNy, 1, dXY, nen_e)
        dYZ2 = Transpose(dYZ) + dYZ
        dXZ2 = Transpose(dXZ) + dXZ
        dXY2 = Transpose(dXY) + dXY
                       
        ! Build B
        !B = [ (X + qe)'*DXX
        !    (X + qe)'*DYY
        !    (X + qe)'*DZZ
        !    (X + qe)'*(DYZ' + DYZ)
        !    (X + qe)'*(DXZ' + DXZ)
        !    (X + qe)'*(DXY' + DXY)]
        tempA = 0.
        ! Part 1: (X + qe)'*DXX
        call DSYMM('L', 'L', nen_e, 3, 1., dXX, nen_e, XYZQe, nen_e, 0., tempA, nen_e)
        Bt(1:nee,1) = reshape(tempA, (/ nee /) )       
        ! Part 2: (X + qe)'*DYY
        call DSYMM('L', 'L', nen_e, 3, 1., dYY, nen_e, XYZQe, nen_e, 0., tempA, nen_e)
        Bt(1:nee,2) = reshape(tempA, (/ nee /) )
        ! Part 3: (X + qe)'*DZZ
        call DSYMM('L', 'L', nen_e, 3, 1., dZZ, nen_e, XYZQe, nen_e, 0., tempA, nen_e)
        Bt(1:nee,3) = reshape(tempA, (/ nee /) )
        ! Part 4: (X + qe)'*(DYZ' + DYZ)
        call DSYMM('L', 'L', nen_e, 3, 1., dYZ2, nen_e, XYZQe, nen_e, 0., tempA, nen_e)
        Bt(1:nee,4) = reshape(tempA, (/ nee /) )
        ! Part 5: (X + qe)'*(DXZ' + DXZ)
        call DSYMM('L', 'L', nen_e, 3, 1., dXZ2, nen_e, XYZQe, nen_e, 0., tempA, nen_e)
        Bt(1:nee,5) = reshape(tempA, (/ nee /) )
        ! Part 6: (X + qe)'*(DXY' + DXY)
        call DSYMM('L', 'L', nen_e, 3, 1., dXY2, nen_e, XYZQe, nen_e, 0., tempA, nen_e)
        Bt(1:nee,6) = reshape(tempA, (/ nee /) )
        
        !
        ! Build F (deformation gradient)
        ! 
        ! F = [ Nx*qe, Ny*qe, Nz*qe ] + eye(3)
        F = 0.
        do a = 1,3
            F(a,a) = 1.
        enddo
        call DGEMM('T', 'N', 3, 3, nen_e, 1.d0, Qie, nen_e, DN, nen_e, 1.d0, F, 3)

        !
        ! Compute the SVD:
        !
        ! F = Psi*diag(Sigma)*LambdaT
        ! Remember that DGESVD destroys F on its exit.
        call DGESVD( 'A', 'A', 3, 3, F, 3, Sigma, Psi, 3, LambdaT, 3, WORK, LWORK, INFO)
        
        ! Compute U = Lambda*diag(Sigma)*Lambda'
        do a = 1,3
            tempB(a,1:3) = Sigma(a)*LambdaT(a,1:3)
        enddo
        U = matmul(Transpose(LambdaT),tempB)
        
        ! Compute R = Psi*Lambda'
        R = matmul(Psi,LambdaT)
                
        ! Compute Yy = Tr(U)*eye(3)-U
        TrU = Sigma(1) + Sigma(2) + Sigma(3)
        Yy = -U
        do a = 1,3
            Yy(a,a) = Yy(a,a) + TrU
        enddo                
        call det33(Yy, detYy)
        
        ! Compute L = U - eye(3)
        L = U
        do a = 1,3
            L(a,a) = L(a,a) - 1.
        enddo
                
        ! Compute G = lambda*Tr_L*eye(3) + 2*mu*L
        TrL = TrU - 3.
        G = (2.*mu)*L
        do a = 1,3
            G(a,a) = G(a,a) + lambda*TrL
        enddo
                     
        ! Compute Uinv = Lambda*diag(1/Sigma)*Lambda'
        do a = 1,3
            tempB(a,1:3) = LambdaT(a,1:3)/Sigma(a)
        enddo
        Uinv = matmul(Transpose(LambdaT),tempB)
        
        ! Compute S = inv(U)*G   
        S = matmul(Uinv,G)
        
        ! Convert from a tensor to the Voigt vector
        Svec = (/ S(1,1), S(2,2), S(3,3), S(2,3), S(1,3), S(1,2) /)
        
        !
        ! Compute DS
        !
        do i2 = 1,nee
            ! dF = [Nx(:,k),Ny(:,k),Nz(:,k)]
            ! take advantage: since local DOF numbering is by direction, then only 1 row of dF is non-zero
            dF = 0.
            a = ceiling(real(i2,8)/real(nen_e,8))
            dF(a,1:3) = DN(mod(i2-1,nen_e)+1,1:3)
            
            ! Theta = R'*dF-1/detY*Y*(R'*dF-dF'*R)*Y*U
            tempB = matmul(Transpose(R),dF)
            Theta = tempB - matmul(Yy, matmul(tempB-Transpose(tempB), matmul(Yy, U) ) )/detYy
            TrTheta = Theta(1,1) + Theta(2,2) + Theta(3,3)
            
            ! dS = (lambda*Tr_theta*eye(3)+2*mu*theta-S*theta)*inv(U)
            tempB = 2.*mu*Theta - matmul(S,Theta)
            do a = 1,3
                tempB(a,a) = tempB(a,a) + lambda*TrTheta
            enddo
            tempB = matmul(tempB, Uinv)
            
            ! Load into DS
            DS(1:6,i2) = (/ tempB(1,1), tempB(2,2), tempB(3,3), tempB(2,3), tempB(1,3), tempB(1,2) /)           
            
        enddo    
        
        ! compute dV0
        dV0 = detJ*quad_w(i1)
    
        !
        ! Build ke = ks + kb:
        !
        !kb = kb + B'*DS*dV0
        call DGEMM( 'N', 'N', nee, nee, 6, dV0, Bt, nee, DS, 6, 1., ke, nee)
        
        ! ks = ks + (DXX*S(1) + DYY*S(2) + DZZ*S(3) + (DYZ'+DYZ)*S(4) + (DXZ'+DXZ)*S(5) + (DXY'+DXY)*S(6))*dV0
        ! Load S(i)*dII into tempD, and then add by blocks down the diagonal of ke
        tempD = 0.
        call DAXPY(nen_e*nen_e, Svec(1)*dV0, dXX,  1, tempD, 1)
        call DAXPY(nen_e*nen_e, Svec(2)*dV0, dYY,  1, tempD, 1)
        call DAXPY(nen_e*nen_e, Svec(3)*dV0, dZZ,  1, tempD, 1)
        call DAXPY(nen_e*nen_e, Svec(4)*dV0, dYZ2, 1, tempD, 1)
        call DAXPY(nen_e*nen_e, Svec(5)*dV0, dXZ2, 1, tempD, 1)
        call DAXPY(nen_e*nen_e, Svec(6)*dV0, dXY2, 1, tempD, 1)
        ! Load into ke:
        ke(1:nen_e, 1:nen_e) = ke(1:nen_e, 1:nen_e) + tempD
        ke((nen_e+1):(2*nen_e), (nen_e+1):(2*nen_e)) = ke((nen_e+1):(2*nen_e), (nen_e+1):(2*nen_e)) + tempD
        ke((2*nen_e+1):nee, (2*nen_e+1):nee) = ke((2*nen_e+1):nee, (2*nen_e+1):nee) + tempD
                
        !
        ! Build Qs
        !
        ! Qs = Qs + B'*S*dV0
        call DGEMV( 'N', nee, 6, dV0, Bt, nee, Svec, 1, 1., qs, 1)
    
    enddo    
