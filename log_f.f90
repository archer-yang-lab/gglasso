! --------------------------------------------------
subroutine log_f (bn,bs,ix,iy,gam,nobs,nvars,x,y,w,pf,dfmax,pmax,nlam,flmin,ulam,&
					eps,maxit,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
! --------------------------------------------------
	implicit none
! - - - arg types - - -
	double precision, parameter :: big=9.9e30
	integer :: bn,bs(bn),ix(bn),iy(bn)                      
	integer :: nobs,nvars,dfmax,pmax,nlam,nalam,npass,jerr,maxit
	integer :: idx(pmax),nbeta(nlam)                   
	double precision :: flmin,eps
    double precision :: x(nobs,nvars),y(nobs),w(nobs),pf(bn),ulam(nlam),gam(bn)           
    double precision :: b0(nlam),beta(nvars,nlam),alam(nlam) 
! - - - local declarations - - -
! - - - local declarations - - -                    
    double precision :: d,t,dif,unorm,al,alf
 	double precision, dimension (:), allocatable :: b,oldbeta,r,oldb,u,dd
 	integer, dimension (:), allocatable :: oidx
    integer :: g,j,l,ctr,ierr,ni,me,start,end
! - - - begin - - -           	    
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
	pf=pf*bn/sum(pf)
! - - - some initial setup - - -   
	r = 0.0D0
	b=0.0D0                                                           
	oldbeta=0.0D0
	idx=0                                                                  
	oidx=0                                                                 
	npass=0                                                                
	ni=npass  
! --------- lambda loop ----------------------------                                                                                                                                                                 
	if(flmin < 1.0D0) then             
		alf=flmin**(1.0D0/(nlam-1.0D0))                                                                                              
	endif
	do l=1,nlam   
		if(flmin>=1.0D0) then                             
	    	al=ulam(l)                                                         
	    else 
			if(l > 2) then                                            
	    		al=al*alf
	        else if(l==1) then
				al=big
			else if(l==2) then  
				al=0.0D0
				do g = 1,bn
					if(pf(g)>0.0D0) then
				      	allocate(u(bs(g)),stat=ierr)  
					    if(ierr/=0) return
						u=matmul(y/(1.0D0+exp(r)),x(:,ix(g):iy(g)))/nobs
			    		al=max(al,sqrt(dot_product(u,u))/pf(g))
						deallocate(u)
					endif
				end do
				al=al*alf
			endif
		endif
		ctr=0     
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
					u=matmul(y/(1.0D0+exp(r)),x(:,start:end))/nobs
					u=gam(g)*b(start:end)+u
					unorm=sqrt(dot_product(u,u)) 
					t=unorm-pf(g)*al
					if(t>0.0D0) then
		      			b(start:end)=u*t/(gam(g)*unorm)                                                                                                                    
					else           
						b(start:end)=0.0D0
					endif
                    dd=b(start:end)-oldb     		      		                                                                                           
		      		if(any(abs(dd)>0.0D0)) then                                         
		      			dif=max(dif,maxval(abs(dd)))                                                  
						r=r+y*matmul(x(:,start:end),dd)
		      			if(oidx(g)==0) then                                           
		      				ni=ni+1      
		      				if(ni>pmax) exit                                             
		      				oidx(g)=ni                                                            
		      				idx(ni)=g      
		                endif
					endif
					deallocate(u,dd,oldb)
				enddo  
				if(ni>pmax) exit	                                            
			    d = sum(y/(1.0D0+exp(r)))  
				d = 4.0D0*d/nobs                                                     
			    if(d /= 0.0D0) then                                            
			      	b(0)=b(0)+d    
			   		r=r+y*d                                                                                                                    
			      	dif=max(dif,abs(d))                                                  
				endif
			    if(dif<eps) exit    
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
						u=matmul(y/(1.0D0+exp(r)),x(:,start:end))/nobs
						u=gam(g)*b(start:end)+u
						unorm=sqrt(dot_product(u,u)) 
						t=unorm-pf(g)*al
						if(t>0.0D0) then
			      			b(start:end)=u*t/(gam(g)*unorm)                                                                                                                    
						else           
							b(start:end)=0.0D0
						endif
	                    dd=b(start:end)-oldb     		      		                                                                                           
			      		if(any(abs(dd)>0.0D0)) then                                         
			      			dif=max(dif,maxval(abs(dd)))                                                  
							r=r+y*matmul(x(:,start:end),dd)
						endif
						deallocate(u,dd,oldb)                          
					enddo                           
				    d = sum(y/(1.0D0+exp(r)))  
					d = 4.0D0*d/nobs                                                        
				    if(d/=0.0D0) then                                            
				      	b(0)=b(0)+d    
				   		r=r+y*d                                                                                                                    
				      	dif=max(dif,abs(d))                                                  
					endif  
					if(dif<eps) exit   
				enddo
			enddo                                                      
		    if(ni>pmax) exit  
!--- final check ------------------------
		    if(all(abs(b-oldbeta)<eps)) exit
			ctr=ctr+1                                                                                                                
		  	if(ctr > maxit) then                                         
		    	jerr=-l                                                            
		      	return        
		    endif                
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
	    me=0
		do j=1,ni   
		  	g=idx(j) 
			if(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
		enddo
		if(me>dfmax) exit                                                                                                           
	enddo    
	deallocate(b,oldbeta,r,oidx)                                         
	return                                                               
end subroutine log_f