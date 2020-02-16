! --------------------------------------------------------------------------
! gglasso.f90: the BMD algorithm for group-lasso penalized learning.
! --------------------------------------------------------------------------
! 
! USAGE:
! 
! SUBROUTINE ls_f (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
!                     eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
! 
! SUBROUTINE log_f (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
!                     eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
! 
! SUBROUTINE hsvm_f (delta,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
!                     eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
! 
! SUBROUTINE sqsvm_f (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
!                     eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
! 
! INPUT ARGUMENTS:
!    delta = delta parameter in Huberized hinge loss, only available in HSVM case (hsvm_f).
!    bn = number of groups
!    bs(bn) = size of each group
!    ix(bn) = first index for each group
!    iy(bn) = last index for each group
!    gam(bn) = upper bound gamma_k in MM algorithm
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an observation vector.
!    y(nobs) = response variable. This argument should be in {-inf, inf} for regression. 
!                and should be a two-level factor {-1, 1} for classification.
!    pf(bn) = relative penalties for each group
!                pf(k) = 0 => kth group unpenalized
!    dfmax = limit the maximum number of variables in the model.
!            (one of the stopping criterion)
!    pmax = limit the maximum number of variables ever to be nonzero. 
!           For example once beta enters the model, no matter how many 
!           times it exits or re-enters model through the path, it will 
!           be counted only once. 
!    nlam = the number of lambda values
!    flmin = user control of lambda values (>=0)
!            flmin < 1.0 => minimum lambda = flmin*(largest lambda value)
!            flmin >= 1.0 => use supplied lambda values (see below)
!    ulam(nlam) = user supplied lambda values (ignored if flmin < 1.0)
!    eps = convergence threshold for coordinate majorization descent. 
!          Each inner coordinate majorization descent loop continues 
!          until the relative change in any coefficient is less than eps.
!    maxit = maximum number of outer-loop iterations allowed at fixed lambda value. 
!            (suggested values, maxit = 100000)
!    intr = whether to include the intercept in the model
! 
!
! OUTPUT:
! 
!    nalam = actual number of lambda values (solutions)
!    b0(nvars) = intercept values for each solution
!    beta(nvars, nlam) = compressed coefficient values for each solution
!    idx(pmax) = pointers to compressed coefficients
!    nbeta(nlam) = number of compressed coefficients for each solution
!    alam(nlam) = lambda values corresponding to each solution
!    npass = actual number of passes over the data for all lambda values
!    jerr = error flag:
!           jerr  = 0 => no error
!           jerr > 0 => fatal error - no output returned
!                    jerr < 7777 => memory allocation error
!                    jerr = 10000 => maxval(vp) <= 0.0
!           jerr < 0 => non fatal error - partial output:
!                    Solutions for larger lambdas (1:(k-1)) returned.
!                    jerr = -k => convergence for kth lambda value not reached
!                           after maxit (see above) iterations.
!                    jerr = -10000-k => number of non zero coefficients along path
!                           exceeds pmax (see above) at kth lambda value.
! 
! LICENSE: GNU GPL (version 2 or later)
! 
! AUTHORS:
!    * Yi Yang (yi.yang6@mcgill.ca) and + Hui Zou (hzou@stat.umn.edu), 
!    * Department of Mathematics and Statistics, McGill University
!    + School of Statistics, University of Minnesota.
! 
! REFERENCES:
!    Yang, Y. and Zou, H. (2015). 
!    A Fast Unified Algorithm for Computing Group-Lasso Penalized Learning Problems
!    Statistics and Computing.
!    25(6), 1129-1141.

! --------------------------------------------------
SUBROUTINE log_f (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                    eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    DOUBLE PRECISION, PARAMETER :: big=9.9E30
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
    INTEGER, PARAMETER :: mnlam = 6
    INTEGER:: mnl
    INTEGER:: bn
    INTEGER::bs(bn)
    INTEGER::ix(bn)
    INTEGER::iy(bn)
    INTEGER:: nobs
    INTEGER::nvars
    INTEGER::dfmax
    INTEGER::pmax
    INTEGER::nlam
    INTEGER::nalam
    INTEGER::npass
    INTEGER::jerr
    INTEGER::maxit
    INTEGER::intr
    INTEGER:: idx(pmax)
    INTEGER::nbeta(nlam)
    DOUBLE PRECISION:: flmin
    DOUBLE PRECISION::eps
    DOUBLE PRECISION:: x(nobs,nvars)
    DOUBLE PRECISION::y(nobs)
    DOUBLE PRECISION::pf(bn)
    DOUBLE PRECISION::ulam(nlam)
    DOUBLE PRECISION::gam(bn)
    DOUBLE PRECISION:: b0(nlam)
    DOUBLE PRECISION::beta(nvars,nlam)
    DOUBLE PRECISION::alam(nlam)
    ! - - - local declarations - - -
    DOUBLE PRECISION:: max_gam
    DOUBLE PRECISION::d
    DOUBLE PRECISION::t
    DOUBLE PRECISION::dif
    DOUBLE PRECISION::unorm
    DOUBLE PRECISION::al
    DOUBLE PRECISION::alf
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
    INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
    INTEGER:: g
    INTEGER::j
    INTEGER::l
    INTEGER::ni
    INTEGER::me
    INTEGER::start
    INTEGER::end
    ! - - - local declarations - - -
    DOUBLE PRECISION:: tlam
    INTEGER:: jx
    INTEGER:: jxx(bn)
    DOUBLE PRECISION:: ga(bn)
    DOUBLE PRECISION:: vl(nvars)
    DOUBLE PRECISION:: al0
! - - - begin - - -
! - - - allocate variables - - -
    ALLOCATE(b(0:nvars))
    ALLOCATE(oldbeta(0:nvars))
    ALLOCATE(r(1:nobs))
    ALLOCATE(oidx(1:bn))
! - - - checking pf - - -
    IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
    ENDIF
    pf=max(0.0D0,pf)
! - - - some initial setup - - -
    jxx = 0
    al = 0.0D0
    mnl = Min (mnlam, nlam)
    r = 0.0D0
    b=0.0D0
    oldbeta=0.0D0
    idx=0
    oidx=0
    npass=0
    ni=npass
    alf = 0.0D0
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf=flmin**(1.0D0/(nlam-1.0D0))
    ENDIF
    vl = matmul(y/(1.0D0+exp(r)), x) / nobs
    DO g = 1,bn
            ALLOCATE(u(bs(g)))
            u = vl(ix(g):iy(g))
            ga(g) = sqrt(dot_product(u,u))
            DEALLOCATE(u)
    END DO
    DO l=1, nlam
        al0 = al
        IF(flmin>=1.0D0) THEN
            al=ulam(l)
        ELSE
            IF(l > 2) THEN
                al=al*alf
            ELSE IF(l==1) THEN
                al=big
            ELSE IF(l==2) THEN
                al0 = 0.0D0
                DO g = 1,bn
                    IF(pf(g)>0.0D0) THEN
                        al0 = max(al0, ga(g) / pf(g))
                    ENDIF
                END DO
                al = al0 * alf
            ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
            IF(jxx(g) == 1) CYCLE
            IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
! --------- outer loop ----------------------------
        DO
            oldbeta(0)=b(0)
            IF(ni>0) THEN
                DO j=1,ni
                    g=idx(j)
                    oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
                ENDDO
            ENDIF
! --middle loop-------------------------------------
            DO
                npass=npass+1
                dif=0.0D0
                DO g=1,bn
                    IF(jxx(g) == 0) CYCLE
                    start=ix(g)
                    end=iy(g)
                    ALLOCATE(u(bs(g)))
                    ALLOCATE(dd(bs(g)))
                    ALLOCATE(oldb(bs(g)))
                    oldb=b(start:end)
                    u=matmul(y/(1.0D0+exp(r)),x(:,start:end))/nobs
                    u=gam(g)*b(start:end)+u
                    unorm=sqrt(dot_product(u,u))
                    t=unorm-pf(g)*al
                    IF(t>0.0D0) THEN
                        b(start:end)=u*t/(gam(g)*unorm)
                    ELSE
                        b(start:end)=0.0D0
                    ENDIF
                    dd=b(start:end)-oldb
                    IF(any(abs(dd)>0.0D0)) THEN
                        dif=max(dif,gam(g)**2*dot_product(dd,dd))
                        r=r+y*matmul(x(:,start:end),dd)
                        IF(oidx(g)==0) THEN
                            ni=ni+1
                            IF(ni>pmax) EXIT
                            oidx(g)=ni
                            idx(ni)=g
                        ENDIF
                    ENDIF
                    DEALLOCATE(u,dd,oldb)
                ENDDO
                IF(intr /= 0) THEN
                    d = sum(y/(1.0D0+exp(r)))
                    d = 4.0D0*d/nobs
                    IF(d /= 0.0D0) THEN
                        b(0)=b(0)+d
                        r=r+y*d
                        dif=max(dif,d**2)
                    ENDIF
                ENDIF
                IF (ni > pmax) EXIT
                IF (dif < eps) EXIT
                IF(npass > maxit) THEN
                    jerr=-l
                    RETURN
               ENDIF
! --inner loop----------------------
                DO
                    npass=npass+1
                    dif=0.0D0
                    DO j=1,ni
                        g=idx(j)
                        start=ix(g)
                        end=iy(g)
                        ALLOCATE(u(bs(g)))
                        ALLOCATE(dd(bs(g)))
                        ALLOCATE(oldb(bs(g)))
                        oldb=b(start:end)
                        u=matmul(y/(1.0D0+exp(r)),x(:,start:end))/nobs
                        u=gam(g)*b(start:end)+u
                        unorm=sqrt(dot_product(u,u))
                        t=unorm-pf(g)*al
                        IF(t>0.0D0) THEN
                            b(start:end)=u*t/(gam(g)*unorm)
                        ELSE
                            b(start:end)=0.0D0
                        ENDIF
                        dd=b(start:end)-oldb
                        IF(any(abs(dd)>0.0D0)) THEN
                            dif=max(dif,gam(g)**2*dot_product(dd,dd))
                            r=r+y*matmul(x(:,start:end),dd)
                        ENDIF
                        DEALLOCATE(u,dd,oldb)
                    ENDDO
                    IF(intr /= 0) THEN
                        d = sum(y/(1.0D0+exp(r)))
                        d = 4.0D0*d/nobs
                        IF(d/=0.0D0) THEN
                            b(0)=b(0)+d
                            r=r+y*d
                            dif=max(dif, d**2)
                        ENDIF
                    ENDIF
                    IF(dif<eps) EXIT
                    IF(npass > maxit) THEN
                        jerr=-l
                        RETURN
                    ENDIF
                ENDDO
            ENDDO
            IF(ni>pmax) EXIT
!--- final check ------------------------
            jx = 0
            max_gam = maxval(gam)
            IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
            IF (jx /= 0) CYCLE
            vl = matmul(y/(1.0D0+exp(r)), x) / nobs
            DO g = 1, bn
                IF(jxx(g) == 1) CYCLE
                ALLOCATE(u(bs(g)))
                u = vl(ix(g):iy(g))
                ga(g) = sqrt(dot_product(u,u))
                IF(ga(g) > al*pf(g))THEN
                    jxx(g) = 1
                    jx = 1
                ENDIF
                DEALLOCATE(u)
            ENDDO
            IF(jx == 1) CYCLE
            EXIT
        ENDDO
!---------- final update variable and save results------------
        IF(ni>pmax) THEN
            jerr=-10000-l
            EXIT
        ENDIF
        IF(ni>0) THEN
            DO j=1,ni
                g=idx(j)
                beta(ix(g):iy(g),l)=b(ix(g):iy(g))
            ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
            g=idx(j)
            IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
    ENDDO
    DEALLOCATE(b,oldbeta,r,oidx)
    RETURN
END SUBROUTINE log_f
