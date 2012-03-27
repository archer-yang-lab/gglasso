! --------------------------------------------------
subroutine hreg_f (delta,bn,bs,ix,iy,gam,nobs,nvars,x,y,w,pf,dfmax,pmax,nlam,flmin,ulam,&
					eps,maxit,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
! --------------------------------------------------
	implicit none
! - - - arg types - - -
	double precision, parameter :: big=9.9e30
	DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
    INTEGER, PARAMETER :: mnlam = 6
      INTEGER :: mnl
	integer :: bn,bs(bn),ix(bn),iy(bn)                      
	integer :: nobs,nvars,dfmax,pmax,nlam,nalam,npass,jerr,maxit
	integer :: idx(pmax),nbeta(nlam)                   
	double precision :: flmin,eps,delta
    double precision :: x(nobs,nvars),y(nobs),w(nobs),pf(bn),ulam(nlam),gam(bn)           
    double precision :: b0(nlam),beta(nvars,nlam),alam(nlam) 
! - - - local declarations - - -
    double precision :: max_gam,d,t,dif,unorm,al,alf,dl(nobs)
 	double precision, dimension (:), allocatable :: b,oldbeta,r,oldb,u,dd
 	integer, dimension (:), allocatable :: oidx
    integer :: g,j,l,ierr,ni,me,start,end
! - - - begin - - -
! - - - local declarations - - -                    
	DOUBLE PRECISION :: tlam
	INTEGER :: jx
	INTEGER :: jxx(bn)
	double precision :: ga(bn)
	double precision :: vl(nvars)
	DOUBLE PRECISION :: al0           	    
! - - - allocate variables - - -
	allocate(b(0:nvars),stat=jerr)                                                                                                
	allocate(oldbeta(0:nvars),stat=ierr)                                          
	jerr=jerr+ierr   
	allocate(r(1:nobs),stat=ierr)                                           
	jerr=jerr+ierr                                                     
	allocate(oidx(1:bn),stat=ierr)                                          
	jerr=jerr+ierr                                                                                                                                                                                                    
	if(jerr/=0) return
! - - - checking pf - - -	
	if(maxval(pf) <= 0.0D0) then
		jerr=10000
		return
	endif
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
	if(flmin < 1.0D0) then
		flmin = Max (mfl, flmin)             
		alf=flmin**(1.0D0/(nlam-1.0D0))                                                                                              
	endif
	dl = 2.0 * sign(min(abs(r) , delta) , r)
    vl = matmul(dl,x) / nobs
    do g = 1,bn
	      	allocate(u(bs(g)),stat=ierr)  
		    if(ierr/=0) return
			u = vl(ix(g):iy(g))
    		ga(g) = sqrt(dot_product(u,u))
			deallocate(u)
	end do
	do l=1,nlam
		al0 = al   
		if(flmin>=1.0D0) then                             
	    	al=ulam(l)                                                         
	    else 
			if(l > 2) then                                            
	    		al=al*alf
	        else if(l==1) then
				al=big
			else if(l==2) then  
				al0 = 0.0D0
				do g = 1,bn
					if(pf(g)>0.0D0) then
			    		al0 = max(al0, ga(g) / pf(g))
					endif
				end do
				al = al0 * alf
			endif
		endif
		tlam = (2.0*al-al0)                                   
        do g = 1, bn
	        if(jxx(g) == 1) cycle
	        if(ga(g) > pf(g) * tlam) jxx(g) = 1
        enddo
! 		call intpr("jxx",-1,jxx,nvars)
! --------- outer loop ----------------------------                                                                                                        
		do  
		    oldbeta(0)=b(0)                                                           
		    if(ni>0) then
				do j=1,ni   
				  	g=idx(j)
					oldbeta(ix(g):iy(g))=b(ix(g):iy(g))    			                       
				enddo
			endif
! --middle loop-------------------------------------
			do
			    npass=npass+1                       
			    dif=0.0D0                                                              
			 	do g=1,bn
					if(jxx(g) == 0) cycle    
					start=ix(g)
					end=iy(g)
			      	allocate(u(bs(g)),stat=ierr)  
				    jerr=jerr+ierr 
				    allocate(dd(bs(g)),stat=ierr)                                           
				    jerr=jerr+ierr                                        
				    allocate(oldb(bs(g)),stat=ierr)                                           
				    jerr=jerr+ierr				                                                                                                                                                              
				    if(jerr/=0) return                
		      		oldb=b(start:end)      
                    dl = 2.0 * sign(min(abs(r) , delta) , r)
					u=matmul(dl,x(:,start:end))/nobs
					u=gam(g)*b(start:end)+u
					unorm=sqrt(dot_product(u,u)) 
					t=unorm-pf(g)*al
					if(t>0.0D0) then
		      			b(start:end)=u*t/(gam(g)*unorm)                                                                                                                    
					else           
						b(start:end)=0.0D0
					endif
                    dd=b(start:end)-oldb     		      		                                                                                           
		      		if(any(dd/=0.0D0)) then                                         
		      			dif=max(dif,gam(g)**2*dot_product(dd,dd))                                                  
						r=r-matmul(x(:,start:end),dd)
		      			if(oidx(g)==0) then                                           
		      				ni=ni+1      
		      				if(ni>pmax) exit                                             
		      				oidx(g)=ni                                                            
		      				idx(ni)=g      
		                endif
					endif
					deallocate(u,dd,oldb)
				enddo  
                dl = 2.0 * sign(min(abs(r) , delta) , r)
				d = sum(dl)                                         
				d = 0.25*d/nobs                                                   
			    if(d/=0.0D0) then                                            
			      	b(0)=b(0)+d    
			   		r=r-d                                                                                                                    
			      	dif=max(dif,d**2)                                                  
				endif  
				IF (ni > pmax) EXIT  
				if (dif < eps) exit
				if(npass > maxit) then                               
		      		jerr=-l                                                  
		      		return
		        endif
! --inner loop----------------------                                                       	
				do                                   
				    npass=npass+1
				    dif=0.0D0       
				 	do j=1,ni                                                     
				      	g=idx(j)  
				        start=ix(g)
				        end=iy(g)                                    
				      	allocate(u(bs(g)),stat=ierr)  
					    jerr=jerr+ierr                                        
					    allocate(dd(bs(g)),stat=ierr)                                           
					    jerr=jerr+ierr                                        
					    allocate(oldb(bs(g)),stat=ierr)                                           
					    jerr=jerr+ierr				                                                                                                                                                              
					    if(jerr/=0) return                                              
			      		oldb=b(start:end)      
	                    dl = 2.0 * sign(min(abs(r) , delta) , r)
						u=matmul(dl,x(:,start:end))/nobs
						u=gam(g)*b(start:end)+u
						unorm=sqrt(dot_product(u,u)) 
						t=unorm-pf(g)*al
						if(t>0.0D0) then
			      			b(start:end)=u*t/(gam(g)*unorm)                                                                                                                    
						else           
							b(start:end)=0.0D0
						endif
	                    dd=b(start:end)-oldb     		      		                                                                                           
			      		if(any(dd/=0.0D0)) then                                         
			      			dif=max(dif,gam(g)**2*dot_product(dd,dd))                                                  
							r=r-matmul(x(:,start:end),dd)
						endif
						deallocate(u,dd,oldb)                          
					enddo                           
	                dl = 2.0 * sign(min(abs(r) , delta) , r)
					d = sum(dl)                                         
					d = 0.25*d/nobs                                                   
				    if(d/=0.0D0) then                                            
				      	b(0)=b(0)+d    
				   		r=r-d                                                                                                                    
				      	dif=max(dif,d**2)                                                  
					endif  
					if(dif<eps) exit
					if(npass > maxit) then                               
			      		jerr=-l                                                  
			      		return
			        endif
				enddo
			enddo                                                      
		    if(ni>pmax) exit  
!--- final check ------------------------
		    jx = 0
		    max_gam = maxval(gam)
		    if(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1     
            IF (jx /= 0) cycle
            dl = 2.0 * sign(min(abs(r) , delta) , r)
		    vl = matmul(dl,x) / nobs
			do g = 1, bn                                            
	            if(jxx(g) == 1) cycle
		      	allocate(u(bs(g)),stat=ierr)  
			    if(ierr/=0) return
				u = vl(ix(g):iy(g))
	    		ga(g) = sqrt(dot_product(u,u))
	            if(ga(g) > al*pf(g))then                          
	      	        jxx(g) = 1
		            jx = 1
!  					call intpr("jx",-1,jx,1)
	            endif              
				deallocate(u)
            enddo 
            if(jx == 1) cycle
            exit                                                                                       
		    enddo     
!---------- final update variable and save results------------                             	                                             
	 	if(ni>pmax) then                                            
	    	jerr=-10000-l                                                      
	      	exit
	    endif
	    if(ni>0) then
		    do j=1,ni   
			  	g=idx(j) 
				beta(ix(g):iy(g),l)=b(ix(g):iy(g)) 
			enddo
		endif
		nbeta(l)=ni                                                         
	   	b0(l)=b(0)                                                         
	   	alam(l)=al                                                          
	   	nalam=l
	    IF (l < mnl) CYCLE   
	    me=0
		do j=1,ni   
		  	g=idx(j) 
			if(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
		enddo
		if(me>dfmax) exit                                                                                                           
	enddo    
	deallocate(b,oldbeta,r,oidx)                                         
	return                                                               
end subroutine hreg_f