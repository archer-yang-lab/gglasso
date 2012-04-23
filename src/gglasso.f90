! --------------------------------------------------
SUBROUTINE ls_f (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                    eps,maxit,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
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
    INTEGER::ierr
    INTEGER::ni
    INTEGER::me
    INTEGER::start
    INTEGER::end
    ! - - - begin - - -
    ! - - - local declarations - - -
    DOUBLE PRECISION:: tlam
    INTEGER:: jx
    INTEGER:: jxx(bn)
    DOUBLE PRECISION:: ga(bn)
    DOUBLE PRECISION:: vl(nvars)
    DOUBLE PRECISION:: al0
! - - - allocate variables - - -
    ALLOCATE(b(0:nvars),STAT=jerr)
    ALLOCATE(oldbeta(0:nvars),STAT=ierr)
    jerr=jerr+ierr
    ALLOCATE(r(1:nobs),STAT=ierr)
    jerr=jerr+ierr
    ALLOCATE(oidx(1:bn),STAT=ierr)
    jerr=jerr+ierr
    IF(jerr/=0) RETURN
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
    r = y
    b=0.0D0
    oldbeta=0.0D0
    idx=0
    oidx=0
    npass=0
    ni=npass
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf=flmin**(1.0D0/(nlam-1.0D0))
    ENDIF
    vl = matmul(r, x)/nobs
    DO g = 1,bn
            ALLOCATE(u(bs(g)),STAT=ierr)
            IF(ierr/=0) RETURN
            u = vl(ix(g):iy(g))
            ga(g) = sqrt(dot_product(u,u))
            DEALLOCATE(u)
    END DO
    DO l=1,nlam
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
                    ALLOCATE(u(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    ALLOCATE(dd(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    ALLOCATE(oldb(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    IF(jerr/=0) RETURN
                    oldb=b(start:end)
                    u=matmul(r,x(:,start:end))/nobs
                    u=gam(g)*b(start:end)+u
                    unorm=sqrt(dot_product(u,u))
                    t=unorm-pf(g)*al
                    IF(t>0.0D0) THEN
                        b(start:end)=u*t/(gam(g)*unorm)
                    ELSE
                        b(start:end)=0.0D0
                    ENDIF
                    dd=b(start:end)-oldb
                    IF(any(dd/=0.0D0)) THEN
                        dif=max(dif,gam(g)**2*dot_product(dd,dd))
                        r=r-matmul(x(:,start:end),dd)
                        IF(oidx(g)==0) THEN
                            ni=ni+1
                            IF(ni>pmax) EXIT
                            oidx(g)=ni
                            idx(ni)=g
                        ENDIF
                    ENDIF
                    DEALLOCATE(u,dd,oldb)
                ENDDO
                d=sum(r)/nobs
                IF(d/=0.0D0) THEN
                    b(0)=b(0)+d
                    r=r-d
                    dif=max(dif,d**2)
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
                        ALLOCATE(u(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        ALLOCATE(dd(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        ALLOCATE(oldb(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        IF(jerr/=0) RETURN
                        oldb=b(start:end)
                        u=matmul(r,x(:,start:end))/nobs
                        u=gam(g)*b(start:end)+u
                        unorm=sqrt(dot_product(u,u))
                        t=unorm-pf(g)*al
                        IF(t>0.0D0) THEN
                            b(start:end)=u*t/(gam(g)*unorm)
                        ELSE
                            b(start:end)=0.0D0
                        ENDIF
                        dd=b(start:end)-oldb
                        IF(any(dd/=0.0D0)) THEN
                            dif=max(dif,gam(g)**2*dot_product(dd,dd))
                            r=r-matmul(x(:,start:end),dd)
                        ENDIF
                        DEALLOCATE(u,dd,oldb)
                    ENDDO
                    d=sum(r)/nobs
                    IF(d/=0.0D0) THEN
                        b(0)=b(0)+d
                        r=r-d
                        dif=max(dif,d**2)
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
            vl = matmul(r, x)/nobs
            DO g = 1, bn
                IF(jxx(g) == 1) CYCLE
                ALLOCATE(u(bs(g)),STAT=ierr)
                IF(ierr/=0) RETURN
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
END SUBROUTINE ls_f



! --------------------------------------------------
SUBROUTINE log_f (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                    eps,maxit,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
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
    INTEGER::ierr
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
    ALLOCATE(b(0:nvars),STAT=jerr)
    ALLOCATE(oldbeta(0:nvars),STAT=ierr)
    jerr=jerr+ierr
    ALLOCATE(r(1:nobs),STAT=ierr)
    jerr=jerr+ierr
    ALLOCATE(oidx(1:bn),STAT=ierr)
    jerr=jerr+ierr
    IF(jerr/=0) RETURN
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
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf=flmin**(1.0D0/(nlam-1.0D0))
    ENDIF
    vl = matmul(y/(1.0D0+exp(r)), x) / nobs
    DO g = 1,bn
            ALLOCATE(u(bs(g)),STAT=ierr)
            IF(ierr/=0) RETURN
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
                    ALLOCATE(u(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    ALLOCATE(dd(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    ALLOCATE(oldb(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    IF(jerr/=0) RETURN
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
                d = sum(y/(1.0D0+exp(r)))
                d = 4.0D0*d/nobs
                IF(d /= 0.0D0) THEN
                    b(0)=b(0)+d
                    r=r+y*d
                    dif=max(dif,d**2)
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
                        ALLOCATE(u(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        ALLOCATE(dd(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        ALLOCATE(oldb(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        IF(jerr/=0) RETURN
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
                    d = sum(y/(1.0D0+exp(r)))
                    d = 4.0D0*d/nobs
                    IF(d/=0.0D0) THEN
                        b(0)=b(0)+d
                        r=r+y*d
                        dif=max(dif, d**2)
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
                ALLOCATE(u(bs(g)),STAT=ierr)
                IF(ierr/=0) RETURN
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

! --------------------------------------------------
SUBROUTINE hsvm_f (delta,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                    eps,maxit,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
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
    INTEGER:: idx(pmax)
    INTEGER::nbeta(nlam)
    DOUBLE PRECISION:: flmin
    DOUBLE PRECISION::eps
    DOUBLE PRECISION::delta
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
    DOUBLE PRECISION::dl(nobs)
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
    INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
    INTEGER:: i
    INTEGER::g
    INTEGER::j
    INTEGER::l
    INTEGER::ierr
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
    ALLOCATE(b(0:nvars), STAT = jerr)
    ALLOCATE(oldbeta(0:nvars), STAT = ierr)
    jerr=jerr + ierr
    ALLOCATE(r(1:nobs), STAT = ierr)
    jerr=jerr + ierr
    ALLOCATE(oidx(1:bn), STAT = ierr)
    jerr=jerr + ierr
    IF(jerr /= 0) RETURN
! - - - checking pf - - -
    IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
    ENDIF
    pf = max(0.0D0, pf)
! - - - some initial setup - - -
    jxx = 0
    al = 0.0D0
    mnl = Min (mnlam, nlam)
    r = 0.0D0
    b = 0.0D0
    oldbeta = 0.0D0
    idx = 0
    oidx = 0
    npass = 0
    ni = npass
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
    ENDIF
    vl = 0.0
    CALL hsvmdrv(delta,nobs,nvars,x,y,r,vl)
    DO g = 1,bn
            ALLOCATE(u(bs(g)),STAT=ierr)
            IF(ierr/=0) RETURN
            u = vl(ix(g):iy(g))
            ga(g) = sqrt(dot_product(u,u))
            DEALLOCATE(u)
    END DO
    DO l=1,nlam
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
                    ALLOCATE(u(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    ALLOCATE(dd(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    ALLOCATE(oldb(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    IF(jerr/=0) RETURN
                    oldb=b(start:end)
                    u = 0.0D0
                    DO i = 1,nobs
                        IF (r(i) > 1.0D0) THEN
                            dl(i) = 0.0D0
                        ELSEIF (r(i) <= (1-delta)) THEN
                            dl(i) = 1.0D0
                        ELSE
                            dl(i) = (1.0D0 - r(i)) / delta
                        ENDIF
                        u = u + dl(i)*y(i)*x(i,start:end)/nobs
                    ENDDO
                    u=gam(g)*b(start:end) + u
                    unorm=sqrt(dot_product(u,u))
                    t=unorm-pf(g)*al
                    IF(t>0.0D0) THEN
                        b(start:end)=u*t/(gam(g)*unorm)
                    ELSE
                        b(start:end)=0.0D0
                    ENDIF
                    dd=b(start:end)-oldb
                    IF(any(dd/=0.0D0)) THEN
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
                d = 0.0D0
                DO i = 1,nobs
                    IF (r(i) > 1.0D0) THEN
                        dl(i) = 0.0D0
                    ELSEIF (r(i) <= (1-delta)) THEN
                        dl(i) = 1.0D0
                    ELSE
                        dl(i) = (1.0D0 - r(i)) / delta
                    ENDIF
                    d = d + dl(i)*y(i)
                ENDDO
                d = 0.5 * delta * d / nobs
                IF(d/=0.0D0) THEN
                    b(0)=b(0)+d
                    r=r+y*d
                    dif=max(dif, d**2)
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
                        ALLOCATE(u(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        ALLOCATE(dd(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        ALLOCATE(oldb(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        IF(jerr/=0) RETURN
                        oldb=b(start:end)
                        u = 0.0D0
                        DO i = 1,nobs
                            IF (r(i) > 1.0D0) THEN
                                dl(i) = 0.0D0
                            ELSEIF (r(i) <= (1-delta)) THEN
                                dl(i) = 1.0D0
                            ELSE
                                dl(i) = (1.0D0 - r(i)) / delta
                            ENDIF
                            u = u + dl(i)*y(i)*x(i,start:end)/nobs
                        ENDDO
                        u=gam(g)*b(start:end) + u
                        unorm=sqrt(dot_product(u,u))
                        t=unorm-pf(g)*al
                        IF(t>0.0D0) THEN
                            b(start:end)=u*t/(gam(g)*unorm)
                        ELSE
                            b(start:end)=0.0D0
                        ENDIF
                        dd=b(start:end)-oldb
                        IF(any(dd/=0.0D0)) THEN
                            dif=max(dif,gam(g)**2*dot_product(dd,dd))
                            r=r+y*matmul(x(:,start:end),dd)
                        ENDIF
                        DEALLOCATE(u,dd,oldb)
                    ENDDO
                    d = 0.0D0
                    DO i = 1,nobs
                        IF (r(i) > 1.0D0) THEN
                            dl(i) = 0.0D0
                        ELSEIF (r(i) <= (1-delta)) THEN
                            dl(i) = 1.0D0
                        ELSE
                            dl(i) = (1.0D0 - r(i)) / delta
                        ENDIF
                        d = d + dl(i)*y(i)
                    ENDDO
                    d = 0.5 * delta * d / nobs
                    IF(d/=0.0D0) THEN
                        b(0)=b(0)+d
                        r=r+y*d
                        dif=max(dif, d**2)
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
            CALL hsvmdrv(delta,nobs,nvars,x,y,r,vl)
            DO g = 1, bn
                IF(jxx(g) == 1) CYCLE
                ALLOCATE(u(bs(g)),STAT=ierr)
                IF(ierr/=0) RETURN
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
END SUBROUTINE hsvm_f





SUBROUTINE hsvmdrv(delta,nobs,nvars,x,y,r,vl)
IMPLICIT NONE
INTEGER:: nobs
INTEGER:: nvars
INTEGER:: i
DOUBLE PRECISION:: delta
DOUBLE PRECISION:: dl(nobs)
DOUBLE PRECISION:: y(nobs)
DOUBLE PRECISION:: r(nobs)
DOUBLE PRECISION:: x(nobs,nvars)
DOUBLE PRECISION:: vl(nvars)
 vl = 0.0
 DO i = 1, nobs
     IF (r(i) > 1.0D0) THEN
        dl (i) = 0.0D0
     ELSEIF (r(i) <= (1-delta)) THEN
        dl (i) = 1.0D0
     ELSE
        dl (i) = (1.0D0 - r(i)) / delta
     ENDIF
 ENDDO
 vl = matmul(dl*y, x) / nobs
END SUBROUTINE hsvmdrv

! --------------------------------------------------
SUBROUTINE sqsvm_f (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                    eps,maxit,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
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
    DOUBLE PRECISION::dl(nobs)
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
    INTEGER::ierr
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
    ALLOCATE(b(0:nvars),STAT=jerr)
    ALLOCATE(oldbeta(0:nvars),STAT=ierr)
    jerr=jerr+ierr
    ALLOCATE(r(1:nobs),STAT=ierr)
    jerr=jerr+ierr
    ALLOCATE(oidx(1:bn),STAT=ierr)
    jerr=jerr+ierr
    IF(jerr/=0) RETURN
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
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf=flmin**(1.0D0/(nlam-1.0D0))
    ENDIF
    dl = 2.0D0 * dim (1.0D0, r)
    vl = matmul(dl*y, x) / nobs
    DO g = 1,bn
            ALLOCATE(u(bs(g)),STAT=ierr)
            IF(ierr/=0) RETURN
            u = vl(ix(g):iy(g))
            ga(g) = sqrt(dot_product(u,u))
            DEALLOCATE(u)
    END DO
    DO l=1,nlam
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
                    ALLOCATE(u(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    ALLOCATE(dd(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    ALLOCATE(oldb(bs(g)),STAT=ierr)
                    jerr=jerr+ierr
                    IF(jerr/=0) RETURN
                    oldb=b(start:end)
                    dl = 2.0 * dim(1.0, r)
                    u=matmul(y*dl,x(:,start:end))/nobs
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
                dl = 2.0 * dim(1.0, r)
                d = dot_product(y,dl)
                d = 0.25*d/nobs
                IF(d /= 0.0D0) THEN
                    b(0)=b(0)+d
                    r=r+y*d
                    dif=max(dif,d**2)
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
                        ALLOCATE(u(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        ALLOCATE(dd(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        ALLOCATE(oldb(bs(g)),STAT=ierr)
                        jerr=jerr+ierr
                        IF(jerr/=0) RETURN
                        oldb=b(start:end)
                        dl = 2.0 * dim(1.0, r)
                        u=matmul(y*dl,x(:,start:end))/nobs
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
                    dl = 2.0 * dim(1.0, r)
                    d = dot_product(y,dl)
                    d = 0.25*d/nobs
                    IF(d/=0.0D0) THEN
                        b(0)=b(0)+d
                        r=r+y*d
                        dif=max(dif,d**2)
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
            dl = 2.0D0 * dim (1.0D0, r)
            vl = matmul(dl*y, x) / nobs
            DO g = 1, bn
                IF(jxx(g) == 1) CYCLE
                ALLOCATE(u(bs(g)),STAT=ierr)
                IF(ierr/=0) RETURN
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
END SUBROUTINE sqsvm_f
