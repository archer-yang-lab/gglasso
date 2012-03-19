c
c                          newGLMnet (6/1/11)
c
c
c                 Elastic net with squared-error loss
c
c dense predictor matrix:
c
c call elnet(ka,parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,maxit,
c            lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call spelnet(ka,parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,
c             isd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse row format
c
c
c other inputs:
c
c   ka = algorithm flag
c      ka=1 => covariance updating algorithm
c      ka=2 => naive algorithm
c   parm = penalty member index (0 <= parm <= 1)
c        = 0.0 => ridge
c        = 1.0 => lasso
c   no = number of observations
c   ni = number of predictor variables
c   y(no) = response vector (overwritten)
c   w(no)= observation weights (overwritten)
c   jd(jd(1)+1) = predictor variable deletion flag
c      jd(1) = 0  => use all variables
c      jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
c   vp(ni) = relative penalties for each predictor variable
c      vp(j) = 0 => jth variable unpenalized
c   ne = maximum number of variables allowed to enter largest model
c        (stopping criterion)
c   nx = maximum number of variables allowed to enter all models
c        along path (memory allocation, nx > ne).
c   nlam = (maximum) number of lamda values
c   flmin = user control of lamda values (>=0)
c      flmin < 1.0 => minimum lamda = flmin*(largest lamda value)
c      flmin >= 1.0 => use supplied lamda values (see below)
c   ulam(nlam) = user supplied lamda values (ignored if flmin < 1.0)
c   thr = convergence threshold for each lamda solution.
c      iterations stop when the maximum reduction in the criterion value
c      as a result of each parameter update over a single pass
c      is less than thr times the null criterion value.
c      (suggested value, thr=1.0e-5)
c   isd = standarization flag:
c      isd = 0 => regression on original predictor variables
c      isd = 1 => regression on standardized predictor variables
c      Note: output solutions always reference original
c            variables locations and scales.
c   maxit = maximum allowed number of passes over the data for all lambda
c      values (suggested values, maxit = 100000)
c
c output:
c
c   lmu = actual number of lamda values (solutions)
c   a0(lmu) = intercept values for each solution
c   ca(nx,lmu) = compressed coefficient values for each solution
c   ia(nx) = pointers to compressed coefficients
c   nin(lmu) = number of compressed coefficients for each solution
c   rsq(lmu) = R**2 values for each solution
c   alm(lmu) = lamda values corresponding to each solution
c   nlp = actual number of passes over the data for all lamda values
c   jerr = error flag:
c      jerr  = 0 => no error
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 10000 => maxval(vp) <= 0.0
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
c
c
c
c least-squares utility routines:
c
c
c uncompress coefficient vectors for all solutions:
c
c call solns(ni,nx,lmu,ca,ia,nin,b)
c
c input:
c
c    ni,nx = input to elnet
c    lmu,ca,ia,nin = output from elnet
c
c output:
c
c    b(ni,lmu) = all elnet returned solutions in uncompressed format
c
c
c uncompress coefficient vector for particular solution:
c
c call uncomp(ni,ca,ia,nin,a)
c
c input:
c
c    ni = total number of predictor variables
c    ca(nx) = compressed coefficient values for the solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for the solution
c
c output:
c
c    a(ni) =  uncompressed coefficient vector
c             referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor matrix:
c
c call modval(a0,ca,ia,nin,n,x,f);
c
c input:
c
c    a0 = intercept
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    n = number of predictor vectors (observations)
c    x(n,ni) = full (uncompressed) predictor matrix
c
c output:
c
c    f(n) = model predictions
c
c
c evaluate linear model from compressed coefficients and
c compressed predictor matrix:
c
c call cmodval(a0,ca,ia,nin,x,ix,jx,n,f);
c
c input:
c
c    a0 = intercept
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    x, ix, jx = predictor matrix in compressed sparse row format
c    n = number of predictor vectors (observations)
c
c output:
c
c    f(n) = model predictions
c
c
c
c
c          Symmetric binomial/multinomial logistic elastic net
c
c
c dense predictor matrix:
c
c call lognet (parm,no,ni,nc,x,y,o,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,
c              maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call splognet (parm,no,ni,nc,x,ix,jx,y,o,jd,vp,ne,nx,nlam,flmin,
c             ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse row format
c
c
c other inputs:
c
c   parm,no,ni,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,maxit
c    = same as elnet above.
c
c   nc = number of classes (distinct outcome values)
c        nc=1 => binomial two-class logistic regression
c            (all output references class 1)
c   y(no,max(2,nc)) = number of each class at each design point(overwritten)
c   o(no,nc) = observation off-sets for each class
c   kopt = optimization flag
c      kopt = 0 => Newton-Raphson (recommended)
c      kpot = 1 => modified Newton-Raphson (sometimes faster)
c
c
c output:
c
c   lmu,ia,nin,alm,nlp = same as elent above
c
c   a0(nc,lmu) = intercept values for each class at each solution
c   ca(nx,nc,lmu) = compressed coefficient values for each class at
c                each solution
c   dev0 = null deviance (intercept only model)
c   fdev(lmu) = fraction of devience explained by each solution
c   jerr = error flag
c      jerr = 0  => no error
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 8000 + k => null probability < 1.0e-5 for class k
c         jerr = 9000 + k => null probability for class k
c                            > 1.0 - 1.0e-5
c         jerr = 10000 => maxval(vp) <= 0.0
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
c    o(no,nc) = training data values for last (lmu_th) solution linear
c               combination.
c
c
c
c logistic/multinomial utilitity routines:
c
c
c uncompress coefficient vectors for all solutions:
c
c call lsolns(ni,nx,nc,lmu,ca,ia,nin,b)
c
c input:
c
c    ni,nx,nc = input to lognet
c    lmu,ca,ia,nin = output from lognet
c
c output:
c
c    b(ni,nc,lmu) = all lognet returned solutions in uncompressed format
c
c
c uncompress coefficient vector for particular solution:
c
c call luncomp(ni,nx,nc,ca,ia,nin,a)
c
c input:
c
c    ni, nx, nc = same as above
c    ca(nx,nc) = compressed coefficient values (for each class)
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients
c
c output:
c
c    a(ni,nc) =  uncompressed coefficient vectors
c                 referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor vectors:
c
c call lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans);
c
c input:
c
c    nt = number of observations
c    x(nt,ni) = full (uncompressed) predictor vectors
c    nc, nx = same as above
c    a0(nc) = intercepts
c    ca(nx,nc) = compressed coefficient values (for each class)
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients
c
c output:
c
c ans(nc,nt) = model predictions
c
c
c evaluate linear model from compressed coefficients and
c compressed predictor matrix:
c
c call lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f);
c
c input:
c
c    nc, nx = same as above
c    a0(nc) = intercept
c    ca(nx,nc) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    x, ix, jx = predictor matrix in compressed sparse row format
c    n = number of predictor vectors (observations)
c
c output:
c
c    f(nc,n) = model predictions
c
c
c
c
c                        Poisson elastic net
c
c
c dense predictor matrix:
c
c call fishnet (parm,no,ni,x,y,o,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,
c               isd,maxit,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c sparse predictor matrix:
c
c call spfishnet (parm,no,ni,x,ix,jx,y,o,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,
c               isd,maxit,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c    x, ix, jx = predictor data matrix in compressed sparse row format
c
c other inputs:
c
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   parm,no,ni,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,maxit
c    = same as elnet above
c
c output:
c
c   lmu,a0,ca,ia,nin,alm = same as elnet above
c   dev0,fdev = same as lognet above
c   nlp = total number of passes over predictor variables
c   jerr = error flag
c      jerr = 0  => no error
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 8888 => negative response count y values
c         jerr = 9999 => no positive observations weights
c         jerr = 10000 => maxval(vp) <= 0.0
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
c    o(no) = training data values for last (lmu_th) solution linear
c            combination.
c
c
c Poisson utility routines:
c
c
c same as elnet above:
c
c    call solns(ni,nx,lmu,ca,ia,nin,b)
c    call uncomp(ni,ca,ia,nin,a)
c    call modval(a0,ca,ia,nin,n,x,f);
c    call cmodval(a0,ca,ia,nin,x,ix,jx,n,f);
c
c compute deviance for given uncompressed data and set of uncompressed
c solutions
c
c call deviance(no,ni,x,y,o,w,nsol,a0,a,flog,jerr)
c
c input:
c
c   no = number of observations
c   ni = number of predictor variables
c   x(no,ni) = predictor data matrix flat file
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   w(no)= observation weights
c   nsol = number of solutions
c   a0(nsol) = intercept for each solution
c   a(ni,nsol) = solution coefficient vectors (uncompressed)
c
c output:
c
c   flog(nsol) = respective deviance values minus null deviance
c   jerr = error flag - see above
c
c
c compute deviance for given compressed data and set of uncompressed solutions
c
c call spdeviance(no,ni,x,ix,jx,y,o,w,nsol,a0,a,flog,jerr)
c
c input:
c
c   no = number of observations
c   ni = number of predictor variables
c   x, ix, jx = predictor data matrix in compressed sparse row format
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   w(no)= observation weights
c   nsol = number of solutions
c   a0(nsol) = intercept for each solution
c   a(ni,nsol) = solution coefficient vectors (uncompressed)
c
c output
c
c   flog(nsol) = respective deviance values minus null deviance
c   jerr = error flag - see above
c
c
c compute deviance for given compressed data and compressed solutions
c
c call cspdeviance(no,x,ix,jx,y,o,w,nx,lmu,a0,ca,ia,nin,flog,jerr)
c
c input:
c
c   no = number of observations
c   x, ix, jx = predictor data matrix in compressed sparse row format
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   w(no)= observation weights
c   nx = input to spfishnet
c   lmu,a0(lmu),ca(nx,lmu),ia(nx),nin(lmu) = output from spfishnet
c
c output
c
c   flog(lmu) = respective deviance values minus null deviance
c   jerr = error flag - see above
c
c
c
c          Elastic net with Cox proportional hazards model
c
c
c dense predictor matrix:
c
c call coxnet (parm,no,ni,x,y,d,o,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,
c              maxit,isd,lmu,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c input:
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c   y(no) = observation times
c   d(no) = died/censored indicator
c       d(i)=0.0 => y(i) = censoring time
c       d(i)=1.0 => y(i) = death time
c   o(no) = observation off-sets
c   parm,no,ni,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,maxit
c                = same as fishnet above
c
c output:
c
c   lmu,ca,ia,nin,dev0,fdev,alm,nlp = same as fishnet above
c   jerr = error flag
c      jerr = 0  => no error - output returned
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 8888 => all observations censored (d(i)=0.0)
c         jerr = 9999 => no positive observations weights
c         jerr = 10000 => maxval(vp) <= 0.0
c         jerr = 20000, 30000 => initialization numerical error
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
c         jerr = -30000-k => numerical error at kth lambda value
c    o(no) = training data values for last (lmu_th) solution linear
c            combination.
c
c
c
c coxnet utility routines:
c
c
c same as elnet above:
c
c    call solns(ni,nx,lmu,ca,ia,nin,b)
c    call uncomp(ni,ca,ia,nin,a)
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor matrix:
c
c call cxmodval(ca,ia,nin,n,x,f);
c
c input:
c
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    n = number of predictor vectors (observations)
c    x(n,ni) = full (uncompressed) predictor matrix
c
c output:
c
c    f(n) = model predictions
c
c
c compute log-likelihood for given data set and vectors of coefficients
c
c call loglike(no,ni,x,y,d,o,w,nvec,a,flog,jerr)
c
c input:
c
c   no = number of observations
c   ni = number of predictor variables
c   x(no,ni) = predictor data matrix flat file
c   y(no) = observation times
c   d(no) = died/censored indicator
c       d(i)=0.0 => y(i) = censoring time
c       d(i)=1.0 => y(i) = death time
c   o(no) = observation off-sets
c   w(no)= observation weights
c   nvec = number of coefficient vectors
c   a(ni,nvec) = coefficient vectors (uncompressed)
c
c output
c
c   flog(nvec) = respective log-likelihood values
c   jerr = error flag - see coxnet above
c
c              
      subroutine elnet  (ka,parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam    609 
     *,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam)                          610
      real ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)                          611
      integer jd(*),ia(nx),nin(nlam)                                        612
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10021                                     615
      jerr=10000                                                            615
      return                                                                615
10021 continue                                                              616
      allocate(vq(1:ni),stat=jerr)                                          616
      if(jerr.ne.0) return                                                  617
      vq=max(0.0,vp)                                                        617
      vq=vq*ni/sum(vq)                                                      618
      if(ka .ne. 1)goto 10041                                               619
      call elnetu  (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd    622 
     *,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                            623
10041 continue                                                              624
      call elnetn (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd,    627 
     *maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                              628
10031 continue                                                              628
      deallocate(vq)                                                        629
      return                                                                630
      end                                                                   631
      subroutine elnetu  (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,t    634 
     *hr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ulam(nlam)                           635
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         636
      integer jd(*),ia(nx),nin(nlam)                                        637
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           642
      allocate(xm(1:ni),stat=ierr)                                          642
      jerr=jerr+ierr                                                        643
      allocate(xs(1:ni),stat=ierr)                                          643
      jerr=jerr+ierr                                                        644
      allocate(ju(1:ni),stat=ierr)                                          644
      jerr=jerr+ierr                                                        645
      allocate(xv(1:ni),stat=ierr)                                          645
      jerr=jerr+ierr                                                        646
      allocate(vlam(1:nlam),stat=ierr)                                      646
      jerr=jerr+ierr                                                        647
      if(jerr.ne.0) return                                                  648
      call chkvars(no,ni,x,ju)                                              649
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  650
      if(maxval(ju) .gt. 0)goto 10071                                       650
      jerr=7777                                                             650
      return                                                                650
10071 continue                                                              651
      call standard(no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)               652
      if(jerr.ne.0) return                                                  653
      if(flmin.ge.1.0) vlam=ulam/ys                                         654
      call elnet1(parm,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxit,x    656 
     *v,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  657
10080 do 10081 k=1,lmu                                                      657
      alm(k)=ys*alm(k)                                                      657
      nk=nin(k)                                                             658
10090 do 10091 l=1,nk                                                       658
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          658
10091 continue                                                              659
10092 continue                                                              659
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         660
10081 continue                                                              661
10082 continue                                                              661
      deallocate(xm,xs,g,ju,xv,vlam)                                        662
      return                                                                663
      end                                                                   664
      subroutine standard (no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)        665
      real x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                  665
      integer ju(ni)                                                        666
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           669
      if(jerr.ne.0) return                                                  670
      w=w/sum(w)                                                            670
      v=sqrt(w)                                                             671
10100 do 10101 j=1,ni                                                       671
      if(ju(j).eq.0)goto 10101                                              672
      xm(j)=dot_product(w,x(:,j))                                           672
      x(:,j)=v*(x(:,j)-xm(j))                                               673
      xv(j)=dot_product(x(:,j),x(:,j))                                      673
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        674
10101 continue                                                              675
10102 continue                                                              675
      if(isd .ne. 0)goto 10121                                              675
      xs=1.0                                                                675
      goto 10131                                                            676
10121 continue                                                              677
10140 do 10141 j=1,ni                                                       677
      if(ju(j).eq.0)goto 10141                                              677
      x(:,j)=x(:,j)/xs(j)                                                   677
10141 continue                                                              678
10142 continue                                                              678
      xv=1.0                                                                679
10131 continue                                                              680
10111 continue                                                              680
      ym=dot_product(w,y)                                                   680
      y=v*(y-ym)                                                            680
      ys=sqrt(dot_product(y,y))                                             680
      y=y/ys                                                                680
      g=0.0                                                                 681
10150 do 10151 j=1,ni                                                       681
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                             681
10151 continue                                                              682
10152 continue                                                              682
      deallocate(v)                                                         683
      return                                                                684
      end                                                                   685
      subroutine elnet1 (beta,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,ulam,thr,    687 
     *maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    688 
     *9)
      real vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    689 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       690
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                           696
      jerr=jerr+ierr                                                        697
      allocate(mm(1:ni),stat=ierr)                                          697
      jerr=jerr+ierr                                                        698
      allocate(da(1:ni),stat=ierr)                                          698
      jerr=jerr+ierr                                                        699
      if(jerr.ne.0) return                                                  700
      bta=beta                                                              700
      omb=1.0-bta                                                           701
      if(flmin .ge. 1.0)goto 10171                                          701
      eqs=max(eps,flmin)                                                    701
      alf=eqs**(1.0/(nlam-1))                                               701
10171 continue                                                              702
      rsq=0.0                                                               702
      a=0.0                                                                 702
      mm=0                                                                  702
      nlp=0                                                                 702
      nin=nlp                                                               702
      iz=0                                                                  702
      mnl=min(mnlam,nlam)                                                   703
10180 do 10181 m=1,nlam                                                     704
      if(flmin .lt. 1.0)goto 10201                                          704
      alm=ulam(m)                                                           704
      goto 10191                                                            705
10201 if(m .le. 2)goto 10211                                                705
      alm=alm*alf                                                           705
      goto 10191                                                            706
10211 if(m .ne. 1)goto 10221                                                706
      alm=big                                                               706
      goto 10231                                                            707
10221 continue                                                              707
      alm=0.0                                                               708
10240 do 10241 j=1,ni                                                       708
      if(ju(j).eq.0)goto 10241                                              708
      if(vp(j).le.0.0)goto 10241                                            709
      alm=max(alm,abs(g(j))/vp(j))                                          710
10241 continue                                                              711
10242 continue                                                              711
      alm=alf*alm/max(bta,1.0e-3)                                           712
10231 continue                                                              713
10191 continue                                                              713
      dem=alm*omb                                                           713
      ab=alm*bta                                                            713
      rsq0=rsq                                                              713
      jz=1                                                                  714
10250 continue                                                              714
10251 continue                                                              714
      if(iz*jz.ne.0) go to 10260                                            714
      nlp=nlp+1                                                             714
      dlx=0.0                                                               715
10270 do 10271 k=1,ni                                                       715
      if(ju(k).eq.0)goto 10271                                              716
      ak=a(k)                                                               716
      u=g(k)+ak*xv(k)                                                       716
      v=abs(u)-vp(k)*ab                                                     716
      a(k)=0.0                                                              717
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         718
      if(a(k).eq.ak)goto 10271                                              719
      if(mm(k) .ne. 0)goto 10291                                            719
      nin=nin+1                                                             719
      if(nin.gt.nx)goto 10272                                               720
10300 do 10301 j=1,ni                                                       720
      if(ju(j).eq.0)goto 10301                                              721
      if(mm(j) .eq. 0)goto 10321                                            721
      c(j,nin)=c(k,mm(j))                                                   721
      goto 10301                                                            721
10321 continue                                                              722
      if(j .ne. k)goto 10341                                                722
      c(j,nin)=xv(j)                                                        722
      goto 10301                                                            722
10341 continue                                                              723
      c(j,nin)=dot_product(x(:,j),x(:,k))                                   724
10301 continue                                                              725
10302 continue                                                              725
      mm(k)=nin                                                             725
      ia(nin)=k                                                             726
10291 continue                                                              727
      del=a(k)-ak                                                           727
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      728
      dlx=max(xv(k)*del**2,dlx)                                             729
10350 do 10351 j=1,ni                                                       729
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                               729
10351 continue                                                              730
10352 continue                                                              730
10271 continue                                                              731
10272 continue                                                              731
      if(dlx.lt.thr)goto 10252                                              731
      if(nin.gt.nx)goto 10252                                               732
      if(nlp .le. maxit)goto 10371                                          732
      jerr=-m                                                               732
      return                                                                732
10371 continue                                                              733
10260 continue                                                              733
      iz=1                                                                  733
      da(1:nin)=a(ia(1:nin))                                                734
10380 continue                                                              734
10381 continue                                                              734
      nlp=nlp+1                                                             734
      dlx=0.0                                                               735
10390 do 10391 l=1,nin                                                      735
      k=ia(l)                                                               735
      ak=a(k)                                                               735
      u=g(k)+ak*xv(k)                                                       735
      v=abs(u)-vp(k)*ab                                                     736
      a(k)=0.0                                                              737
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         738
      if(a(k).eq.ak)goto 10391                                              739
      del=a(k)-ak                                                           739
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      740
      dlx=max(xv(k)*del**2,dlx)                                             741
10400 do 10401 j=1,nin                                                      741
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                  741
10401 continue                                                              742
10402 continue                                                              742
10391 continue                                                              743
10392 continue                                                              743
      if(dlx.lt.thr)goto 10382                                              743
      if(nlp .le. maxit)goto 10421                                          743
      jerr=-m                                                               743
      return                                                                743
10421 continue                                                              744
      goto 10381                                                            745
10382 continue                                                              745
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                      746
10430 do 10431 j=1,ni                                                       746
      if(mm(j).ne.0)goto 10431                                              747
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))            748
10431 continue                                                              749
10432 continue                                                              749
      jz=0                                                                  750
      goto 10251                                                            751
10252 continue                                                              751
      if(nin .le. nx)goto 10451                                             751
      jerr=-10000-m                                                         751
      goto 10182                                                            751
10451 continue                                                              752
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 752
      kin(m)=nin                                                            753
      rsqo(m)=rsq                                                           753
      almo(m)=alm                                                           753
      lmu=m                                                                 754
      if(m.lt.mnl)goto 10181                                                754
      if(flmin.ge.1.0)goto 10181                                            755
      me=0                                                                  755
10460 do 10461 j=1,nin                                                      755
      if(ao(j,m).ne.0.0) me=me+1                                            755
10461 continue                                                              755
10462 continue                                                              755
      if(me.gt.ne)goto 10182                                                756
      if(rsq-rsq0.lt.sml*rsq)goto 10182                                     756
      if(rsq.gt.rsqmax)goto 10182                                           757
10181 continue                                                              758
10182 continue                                                              758
      deallocate(a,mm,c,da)                                                 759
      return                                                                760
      end                                                                   761
      subroutine elnetn (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,th    763 
     *r,isd,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no),w(no),ulam(nlam)                           764
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         765
      integer jd(*),ia(nx),nin(nlam)                                        766
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          771
      allocate(xs(1:ni),stat=ierr)                                          771
      jerr=jerr+ierr                                                        772
      allocate(ju(1:ni),stat=ierr)                                          772
      jerr=jerr+ierr                                                        773
      allocate(xv(1:ni),stat=ierr)                                          773
      jerr=jerr+ierr                                                        774
      allocate(vlam(1:nlam),stat=ierr)                                      774
      jerr=jerr+ierr                                                        775
      if(jerr.ne.0) return                                                  776
      call chkvars(no,ni,x,ju)                                              777
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  778
      if(maxval(ju) .gt. 0)goto 10481                                       778
      jerr=7777                                                             778
      return                                                                778
10481 continue                                                              779
      call standard1(no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)                780
      if(jerr.ne.0) return                                                  781
      if(flmin.ge.1.0) vlam=ulam/ys                                         782
      call elnet2(parm,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,vlam,thr,maxit,x    784 
     *v,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  785
10490 do 10491 k=1,lmu                                                      785
      alm(k)=ys*alm(k)                                                      785
      nk=nin(k)                                                             786
10500 do 10501 l=1,nk                                                       786
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          786
10501 continue                                                              787
10502 continue                                                              787
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         788
10491 continue                                                              789
10492 continue                                                              789
      deallocate(xm,xs,ju,xv,vlam)                                          790
      return                                                                791
      end                                                                   792
      subroutine standard1 (no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)         793
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)                        793
      integer ju(ni)                                                        794
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           797
      if(jerr.ne.0) return                                                  798
      w=w/sum(w)                                                            798
      v=sqrt(w)                                                             799
10510 do 10511 j=1,ni                                                       799
      if(ju(j).eq.0)goto 10511                                              800
      xm(j)=dot_product(w,x(:,j))                                           800
      x(:,j)=v*(x(:,j)-xm(j))                                               801
      xv(j)=dot_product(x(:,j),x(:,j))                                      801
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        802
10511 continue                                                              803
10512 continue                                                              803
      if(isd .ne. 0)goto 10531                                              803
      xs=1.0                                                                803
      goto 10541                                                            804
10531 continue                                                              804
10550 do 10551 j=1,ni                                                       804
      if(ju(j).eq.0)goto 10551                                              804
      x(:,j)=x(:,j)/xs(j)                                                   804
10551 continue                                                              805
10552 continue                                                              805
      xv=1.0                                                                806
10541 continue                                                              807
10521 continue                                                              807
      ym=dot_product(w,y)                                                   807
      y=v*(y-ym)                                                            807
      ys=sqrt(dot_product(y,y))                                             807
      y=y/ys                                                                808
      deallocate(v)                                                         809
      return                                                                810
      end                                                                   811
      subroutine elnet2(beta,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,ulam,thr,m    813 
     *axit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    814 
     *9)
      real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    815 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       816
      real, dimension (:), allocatable :: a,g                                   
      integer, dimension (:), allocatable :: mm,ix                              
      allocate(a(1:ni),stat=jerr)                                           821
      allocate(mm(1:ni),stat=ierr)                                          821
      jerr=jerr+ierr                                                        822
      allocate(g(1:ni),stat=ierr)                                           822
      jerr=jerr+ierr                                                        823
      allocate(ix(1:ni),stat=ierr)                                          823
      jerr=jerr+ierr                                                        824
      if(jerr.ne.0) return                                                  825
      bta=beta                                                              825
      omb=1.0-bta                                                           825
      ix=0                                                                  826
      if(flmin .ge. 1.0)goto 10571                                          826
      eqs=max(eps,flmin)                                                    826
      alf=eqs**(1.0/(nlam-1))                                               826
10571 continue                                                              827
      rsq=0.0                                                               827
      a=0.0                                                                 827
      mm=0                                                                  827
      nlp=0                                                                 827
      nin=nlp                                                               827
      iz=0                                                                  827
      mnl=min(mnlam,nlam)                                                   827
      alm=0.0                                                               828
10580 do 10581 j=1,ni                                                       828
      if(ju(j).eq.0)goto 10581                                              828
      g(j)=abs(dot_product(y,x(:,j)))                                       828
10581 continue                                                              829
10582 continue                                                              829
10590 do 10591 m=1,nlam                                                     829
      alm0=alm                                                              830
      if(flmin .lt. 1.0)goto 10611                                          830
      alm=ulam(m)                                                           830
      goto 10601                                                            831
10611 if(m .le. 2)goto 10621                                                831
      alm=alm*alf                                                           831
      goto 10601                                                            832
10621 if(m .ne. 1)goto 10631                                                832
      alm=big                                                               832
      goto 10641                                                            833
10631 continue                                                              833
      alm0=0.0                                                              834
10650 do 10651 j=1,ni                                                       834
      if(ju(j).eq.0)goto 10651                                              834
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                            834
10651 continue                                                              835
10652 continue                                                              835
      alm0=alm0/max(bta,1.0e-3)                                             835
      alm=alf*alm0                                                          836
10641 continue                                                              837
10601 continue                                                              837
      dem=alm*omb                                                           837
      ab=alm*bta                                                            837
      rsq0=rsq                                                              837
      jz=1                                                                  838
      tlam=bta*(2.0*alm-alm0)                                               839
10660 do 10661 k=1,ni                                                       839
      if(ix(k).eq.1)goto 10661                                              839
      if(ju(k).eq.0)goto 10661                                              840
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                        841
10661 continue                                                              842
10662 continue                                                              842
10670 continue                                                              842
10671 continue                                                              842
      if(iz*jz.ne.0) go to 10260                                            843
10680 continue                                                              843
      nlp=nlp+1                                                             843
      dlx=0.0                                                               844
10690 do 10691 k=1,ni                                                       844
      if(ix(k).eq.0)goto 10691                                              844
      gk=dot_product(y,x(:,k))                                              845
      ak=a(k)                                                               845
      u=gk+ak*xv(k)                                                         845
      v=abs(u)-vp(k)*ab                                                     845
      a(k)=0.0                                                              846
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         847
      if(a(k).eq.ak)goto 10691                                              848
      if(mm(k) .ne. 0)goto 10711                                            848
      nin=nin+1                                                             848
      if(nin.gt.nx)goto 10692                                               849
      mm(k)=nin                                                             849
      ia(nin)=k                                                             850
10711 continue                                                              851
      del=a(k)-ak                                                           851
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        852
      y=y-del*x(:,k)                                                        852
      dlx=max(xv(k)*del**2,dlx)                                             853
10691 continue                                                              854
10692 continue                                                              854
      if(nin.gt.nx)goto 10672                                               855
      if(dlx .ge. thr)goto 10731                                            855
      ixx=0                                                                 856
10740 do 10741 k=1,ni                                                       856
      if(ix(k).eq.1)goto 10741                                              856
      if(ju(k).eq.0)goto 10741                                              857
      g(k)=abs(dot_product(y,x(:,k)))                                       858
      if(g(k) .le. ab*vp(k))goto 10761                                      858
      ix(k)=1                                                               858
      ixx=1                                                                 858
10761 continue                                                              859
10741 continue                                                              860
10742 continue                                                              860
      if(ixx.eq.1) go to 10680                                              861
      goto 10672                                                            862
10731 continue                                                              863
      if(nlp .le. maxit)goto 10781                                          863
      jerr=-m                                                               863
      return                                                                863
10781 continue                                                              864
10260 continue                                                              864
      iz=1                                                                  865
10790 continue                                                              865
10791 continue                                                              865
      nlp=nlp+1                                                             865
      dlx=0.0                                                               866
10800 do 10801 l=1,nin                                                      866
      k=ia(l)                                                               866
      gk=dot_product(y,x(:,k))                                              867
      ak=a(k)                                                               867
      u=gk+ak*xv(k)                                                         867
      v=abs(u)-vp(k)*ab                                                     867
      a(k)=0.0                                                              868
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         869
      if(a(k).eq.ak)goto 10801                                              870
      del=a(k)-ak                                                           870
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        871
      y=y-del*x(:,k)                                                        871
      dlx=max(xv(k)*del**2,dlx)                                             872
10801 continue                                                              873
10802 continue                                                              873
      if(dlx.lt.thr)goto 10792                                              873
      if(nlp .le. maxit)goto 10821                                          873
      jerr=-m                                                               873
      return                                                                873
10821 continue                                                              874
      goto 10791                                                            875
10792 continue                                                              875
      jz=0                                                                  876
      goto 10671                                                            877
10672 continue                                                              877
      if(nin .le. nx)goto 10841                                             877
      jerr=-10000-m                                                         877
      goto 10592                                                            877
10841 continue                                                              878
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 878
      kin(m)=nin                                                            879
      rsqo(m)=rsq                                                           879
      almo(m)=alm                                                           879
      lmu=m                                                                 880
      if(m.lt.mnl)goto 10591                                                880
      if(flmin.ge.1.0)goto 10591                                            881
      me=0                                                                  881
10850 do 10851 j=1,nin                                                      881
      if(ao(j,m).ne.0.0) me=me+1                                            881
10851 continue                                                              881
10852 continue                                                              881
      if(me.gt.ne)goto 10592                                                882
      if(rsq-rsq0.lt.sml*rsq)goto 10592                                     882
      if(rsq.gt.rsqmax)goto 10592                                           883
10591 continue                                                              884
10592 continue                                                              884
      deallocate(a,mm,g,ix)                                                 885
      return                                                                886
      end                                                                   887
      subroutine chkvars(no,ni,x,ju)                                        888
      real x(no,ni)                                                         888
      integer ju(ni)                                                        889
10860 do 10861 j=1,ni                                                       889
      ju(j)=0                                                               889
      t=x(1,j)                                                              890
10870 do 10871 i=2,no                                                       890
      if(x(i,j).eq.t)goto 10871                                             890
      ju(j)=1                                                               890
      goto 10872                                                            890
10871 continue                                                              891
10872 continue                                                              891
10861 continue                                                              892
10862 continue                                                              892
      return                                                                893
      end                                                                   894
      subroutine uncomp(ni,ca,ia,nin,a)                                     895
      real ca(*),a(ni)                                                      895
      integer ia(*)                                                         896
      a=0.0                                                                 896
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                   897
      return                                                                898
      end                                                                   899
      subroutine modval(a0,ca,ia,nin,n,x,f)                                 900
      real ca(nin),x(n,*),f(n)                                              900
      integer ia(nin)                                                       901
      f=a0                                                                  901
      if(nin.le.0) return                                                   902
10880 do 10881 i=1,n                                                        902
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                       902
10881 continue                                                              903
10882 continue                                                              903
      return                                                                904
      end                                                                   905
