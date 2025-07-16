      subroutine r4fftsc  (a,n,st,ct,ifail)
c-----------------------------------------------------------------------
c  DMC -- real precision IMSL replacement routine
c         locally allocated, saved workspace objects for each problem
c         size (n) that is seen.  This may not be thread safe;
c         this is only intended to replace copied IMSL routines in some
c         serial TRANSP-related data processing applications.  For codes
c         requiring high performance FFT, use FFTW 2.1.5.
c
c--------------------------------
c   some IMSL comments retained...
c--------------------------------
c                                                                       
c   purpose             - compute the sine and cosine transforms of     
c                           a real valued sequence                      
c                                                                       
c   usage               - call fftsc (a,n,st,ct,iwk,wk,cwk)             
c                                                                       
c   arguments    a      - input real vector of length n which           
c                           contains the data to be transformed.        
c                n      - input number of data points to be transformed.
c                           n must be a positive even integer.          
c                st     - output real vector of length n/2+1            
c                           containing the coefficients of the          
c                           sine transform.                             
c                ct     - output real vector of length n/2+1            
c                           containing the coefficients of the          
c                           cosine transform.                           
c                                         
c              ifail    - 0 on exit: normal
c                       - 1 on exit: could not allocate workspace
c
c   notation            - information on special notation and           
c                           conventions is available in the manual      
c                           introduction or through imsl routine uhelp  
c                                                                       
c   remarks  1.  fftsc computes the sine transform, st, according       
c                to the following formula;                              
c                                                                       
c                  st(k+1) = 2.0 * sum from j = 0 to n-1 of             
c                            a(j+1)*sin(2.0*pi*j*k/n)                   
c                  for k=0,1,...,n/2 and pi=3.1415...                   
c                                                                       
c                fftsc computes the cosine transform, ct, according     
c                to the following formula;                              
c                                                                       
c                  ct(k+1) = 2.0 * sum from j = 0 to n-1 of             
c                            a(j+1)*cos(2.0*pi*j*k/n)                   
c                  for k=0,1,...,n/2 and pi=3.1415...                   
c            2.  the following relationship exists between the data     
c                and the coefficients of the sine and cosine transform  
c                                                                       
c                  a(j+1) = ct(1)/(2*n) + ct(n/2+1)/(2*n)*(-1)**j +     
c                           sum from k = 1 to n/2-1 of                  
c                            (ct(k+1)/n*cos((2.0*pi*j*k)/n) +           
c                             st(k+1)/n*sin((2.0*pi*j*k)/n))            
c                  for j=0,1,...,n-1 and pi=3.1415...                   
c                                                                       
c                                                                       
c-----------------------------------------------------------------------
c
      implicit NONE
c
c                                  specifications for arguments         
      integer            n
      real               a(n),st(*),ct(*)
      integer            ifail
c
c-----------------------------------------------------------------------
c
      type :: fftwk
         integer :: n_size
         real, dimension(:), pointer :: wka => NULL()
      end type fftwk
c
      real, dimension(:), pointer :: wka => NULL()
c
      integer, SAVE :: n_save = 0
c
c  collection of work arrays:
c
      integer, parameter :: nmax=50
      type (fftwk), dimension(nmax), SAVE :: fw
c
c-----------------------------------------------------------------------
      integer :: isize,in2p1,ino2,k1,k2,i,imatch
      real, dimension(:), allocatable :: awk
c-----------------------------------------------------------------------
c
      ifail = 0
c
      if(n.le.0) return
c
      in2p1 = n/2 + 1

      st(1:in2p1) = 0.0
      ct(1:in2p1) = 0.0
c
c-------------------------------
c
      imatch=0
      do i=1,n_save
         if(fw(i)%n_size.eq.n) then
            imatch=i
            exit
         endif
      enddo

      if(imatch.eq.0) then
         if(n_save.eq.nmax) then
            write(6,*) ' ?r4fftsc: nmax exceeded. '
            ifail=1
            return
         endif

         n_save = n_save + 1
         fw(n_save)%n_size = n

         isize = 2*n + 30
         allocate(fw(n_save)%wka(isize)); fw(n_save)%wka = 0.0

         call r4rffti(n,fw(n_save)%wka)

         imatch = n_save
      endif

      wka => fw(imatch)%wka
c
c-------------------------------
c
      allocate(awk(n))
      awk = a
c
      call r4rfftf(n,awk,wka)
c
      ct(1)=awk(1)*2.0
      st(1)=0.0
c
      k1=0
      k2=1
c
      ino2 = (n+1)/2
c
      do i=2,ino2
         k1=k1+2
         k2=k2+2

         ct(i) = awk(k1)*2.0
         st(i) = -awk(k2)*2.0
      enddo
c
      if(ino2.lt.in2p1) then
         k1=k1+2
         ct(i)=awk(k1)*2.0
      endif

      return
      end
