  MODULE dims
!    dims.inc
!      definition of cliches for pest  2._r8  feb. 23, 1983 jm
!       standard dimensions are used with a letter code eg. go24ij
!       where i can take values s,t,y,z and j can be a,b,c. these set
!       the dimensions of nsf,nths etc. the code table follows:
!            nsf0               nths0      nfm      nfn
!                           e     32         11       11
!         s   51            a     64         21       31
!         t   101           b     128        41       63
!         y   201           c     256        81       127
!         z   401           d     512        161      255
!         w   601     
!         h   1201
!         k   2401
!..............................................................
 INTEGER, PARAMETER :: ipolp=1,iasym = 2, iasym2 = iasym*iasym &
 &         , nmdv = 9, nsg=30,nsg1=nsg+1,nsg2=2*nsg
!
!..............................................................
!      thus
!          nsg      sets the max number of inter-rational surface regions
!          nsf      sets the max. number of psi surfaces + 1
!          nbig     is the larger of nsf and ngs. set through dif
!          nths     sets the max. number of theta points + 5
!          nfe      sets the max. number of finite elements + 3
!          nfm      sets the max. number of fourier modes
!          nfn      sets the max dimension in inverting ' gee '
!          nfv0     equals nfm * nfm
!          nfm21    equals 2 * nfm + 1
!          nfm2     equals nfm21 * nfm21
!          nmt      is the rank of the resulting matrix
!          nmt2     equals nmt * nmt
!          nmdv     sets the max. number of 2*mdiv + 1
!          nkb      equals 3 * nfm
!          nkc      equals 6 * nfm
!              these two dimension pot and kin
!
!     note: a synonymous set exists for grub23
!.
!
!234567890123456789012345678901234567890123456789012345678901234567890123
!
! end cliche
 END MODULE dims

 MODULE pstcom
!
!      definition of cliches for pest package (tearing version) ...
!      version p33lst for cray,general jacobian...rg 2/1/ 82._r8 
!      also contains cliches for inner region and matching...
!
!     cliche pstcom
!
!....................................................................
!
!      the dimensions are set through the use of the parameter statement
!      in dims.inc
!
 USE dims
!
 INTEGER, PARAMETER :: nrun = 20
!
!
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 dpsi
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 dtent
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 psinod
 INTEGER, DIMENSION(:), ALLOCATABLE :: 	 msub
 INTEGER, DIMENSION(:), ALLOCATABLE :: 	 jsub
 INTEGER, DIMENSION(:), ALLOCATABLE :: 	 jtot
 INTEGER, DIMENSION(:), ALLOCATABLE :: 	 lmax
 INTEGER, DIMENSION(:), ALLOCATABLE :: 	 lmin
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 fa
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 fpa
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 ga
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 gpa
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 pa
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 ppa
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 psia
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 psibig
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 qa
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 qpa
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 qppa
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 rhoa
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 alpha
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 psisin
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 sigavr
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 sigavp
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 di
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 dr
 REAL*8, DIMENSION(:), ALLOCATABLE ::      rgradpsi
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 ajphi
 REAL*8, dimension(:), allocatable :: betapol, aiaspect, alqolp
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 psinew
 COMPLEX*16, DIMENSION(nsg,nsg) 	 :: apr
 COMPLEX*16, DIMENSION(nsg,nsg) 	 :: bpr
 COMPLEX*16, DIMENSION(nsg,nsg) 	 :: delpr
 COMPLEX*16, DIMENSION(nsg,nsg) 	 :: gampr
 real*8, DIMENSION(nsg,nsg) 	 :: error_apr
 real*8, DIMENSION(nsg,nsg) 	 :: error_bpr
 real*8, DIMENSION(nsg,nsg) 	 :: error_delpr
 real*8, DIMENSION(nsg,nsg) 	 :: error_gampr
 LOGICAL	 	 	 :: lcont
 LOGICAL	 	 	 :: lsymz
 LOGICAL	 	 	 :: check1
 LOGICAL	 	 	 :: LADDMO
 LOGICAL	 	 	 :: fast
 LOGICAL	 	 	 :: ke
 LOGICAL	 	 	 :: lpless
 LOGICAL	 	 	 :: wall
 LOGICAL	 	 	 :: infwal
 LOGICAL	 	 	 :: lcirc
 LOGICAL	 	 	 :: lcub
 LOGICAL	 	 	 :: lmarg
 LOGICAL	 	 	 :: lsing
 LOGICAL	 	 	 :: lmsin
 LOGICAL	 	 	 :: lsub
 LOGICAL	 	 	 :: lpstps
 LOGICAL	 	 	 :: leqdsk
 LOGICAL	 	 	 :: liter
 LOGICAL	 	 	 :: lrego
 LOGICAL	 	 	 :: lextra
 LOGICAL	 	 	 :: lso
 LOGICAL	 	 	 :: pkbdry
 real*8                  :: wlambda ! delta W eigenvalue
 REAL*8	 	 	 :: n
 REAL*8	 	 	 :: kin
 REAL*8	 	 	 :: no2pi
 REAL*8	 	 	 :: no2pi2
 REAL*8	 	 	 :: zero
 REAL*8	 	 	 :: half
 REAL*8	 	 	 :: one
 REAL*8	 	 	 :: two
 REAL*8	 	 	 :: three
 REAL*8	 	 	 :: four
 REAL*8	 	 	 :: five
 REAL*8	 	 	 :: seven
 REAL*8	 	 	 :: pye
 REAL*8	 	 	 :: twopi
 REAL*8	 	 	 :: twopi2
 REAL*8	 	 	 :: bit
 REAL*8	 	 	 :: pt1
 REAL*8	 	 	 :: amu0
 REAL*8	 	 	 :: b0
 REAL*8	 	 	 :: gamma
 REAL*8	 	 	 :: hdtent
 REAL*8	 	 	 :: dth
 REAL*8	 	 	 :: dthinc
 REAL*8	 	 	 :: dlay
 REAL*8	 	 	 :: blay
 REAL*8	 	 	 :: frs
 REAL*8, DIMENSION(5) 	 :: dlays
 REAL*8	 	 	 :: r2
 REAL*8	 	 	 :: r4
 REAL*8	 	 	 :: r6
 REAL*8	 	 	 :: zma
 REAL*8	 	 	 :: ensq
 REAL*8	 	 	 :: effmax
 REAL*8	 	 	 :: themax
 REAL*8	 	 	 :: psiv
 REAL*8	 	 	 :: psidiv
 REAL*8	 	 	 :: psidif
 REAL*8	 	 	 :: f0
 REAL*8	 	 	 :: f02
 REAL*8	 	 	 :: xpress
 REAL*8	 	 	 :: rbth
 REAL*8, DIMENSION(20) 	 :: rbt
 REAL*8	 	 	 :: alx
 REAL*8	 	 	 :: alz
 REAL*8	 	 	 :: xzero
 REAL*8	 	 	 :: xma
 REAL*8	 	 	 :: r
 REAL*8	 	 	 :: p0
 REAL*8	 	 	 :: gp
 REAL*8	 	 	 :: psimin
 REAL*8	 	 	 :: psilim
 REAL*8	 	 	 :: psivax
 REAL*8	 	 	 :: psipls
 REAL*8	 	 	 :: betag
 REAL*8	 	 	 :: betap
 REAL*8	 	 	 :: wsq
 REAL*8	 	 	 :: wsqold
 REAL*8	 	 	 :: cb1
 REAL*8	 	 	 :: cb2
 REAL*8	 	 	 :: pm
 REAL*8	 	 	 :: upsiln
 REAL*8, DIMENSION(10) 	 :: phi

 REAL*8                  :: totalToroidalCurrent, b0SquareCentre
 REAL*8                  :: betaPoloidal, betaToroidal
 REAL*8                  :: betaStar, TroyonG, gStar, volumeAveragePressure
 REAL*8                  :: volumeAverageBSquare, volumeAverageBeta
 REAL*8                  :: pressurePeakingFactor
 REAL*8                  :: fourPiInductance, majorRadius, minorRadius
 REAL*8                  :: inverseAspectRatio, rMagnetic, rCentre, elongation
 REAL*8                  :: triangularity, plasmaSurface, plasmaVolume

 integer :: nIdealInstabilities, idealInstabilityNode(100)
 integer :: inputFormat
 character*80 :: inputFile

 integer                         :: isolver
 INTEGER	 	 	 :: outpst
 INTEGER	 	 	 :: inkin
 INTEGER	 	 	 :: outmap
 INTEGER	 	 	 :: outmod
 INTEGER	 	 	 :: outilc
 INTEGER	 	 	 :: outd2
 INTEGER	 	 	 :: outmp1
 INTEGER	 	 	 :: over
 INTEGER	 	 	 :: outvec
 INTEGER	 	 	 :: surf
 INTEGER	 	 	 :: scount
 INTEGER	 	 	 :: scl
 INTEGER	 	 	 :: scr
 INTEGER	 	 	 :: nsf0
 INTEGER	 	 	 :: nsf
 INTEGER	 	 	 :: nsf2
 INTEGER	 	 	 :: nfe
 INTEGER	 	 	 :: nke
 INTEGER	 	 	 :: nmt
 INTEGER	 	 	 :: nfe1
 INTEGER	 	 	 :: nmtg
 INTEGER	 	 	 :: nmt2
 INTEGER	 	 	 :: nfn
 INTEGER	 	 	 :: nths0
 INTEGER	 	 	 :: nths
 INTEGER	 	 	 :: nthsh
 INTEGER	 	 	 :: nthss
 INTEGER	 	 	 :: nthsfm
 INTEGER	 	 	 :: nths2
 INTEGER	 	 	 :: nthsq
 INTEGER	 	 	 :: nfm
 INTEGER	 	 	 :: nsg3
 INTEGER	 	 	 :: nkd
 INTEGER	 	 	 :: nfv0
 INTEGER	 	 	 :: nfmgbx
 INTEGER	 	 	 :: nfm21
 INTEGER	 	 	 :: nfm8
 INTEGER	 	 	 :: nkb
 INTEGER	 	 	 :: nkc
 INTEGER	 	 	 :: nad
 INTEGER	 	 	 :: nkc2
 INTEGER, DIMENSION(nrun) 	 :: mm
 INTEGER	 	 	 :: iequil
 INTEGER	 	 	 :: imap
 INTEGER	 	 	 :: imode
 INTEGER	 	 	 :: isymz
 INTEGER	 	 	 :: inpest
 INTEGER	 	 	 :: inpot
 INTEGER	 	 	 :: iomode
 INTEGER	 	 	 :: inmode
 INTEGER	 	 	 :: itty
 INTEGER	 	 	 :: neqdsk
 INTEGER	 	 	 :: nmpdsk
 INTEGER	 	 	 :: nmpout
 INTEGER	 	 	 :: nmod2
 INTEGER	 	 	 :: il
 INTEGER	 	 	 :: n1surf
 INTEGER	 	 	 :: npsurf
 INTEGER	 	 	 :: mth1
 INTEGER	 	 	 :: mth2
 INTEGER	 	 	 :: mthh
 INTEGER	 	 	 :: minc
 INTEGER	 	 	 :: mp
 INTEGER	 	 	 :: ibas
 INTEGER	 	 	 :: ipol
 INTEGER	 	 	 :: mdiv
 INTEGER	 	 	 :: m
 INTEGER	 	 	 :: mat
 INTEGER	 	 	 :: nx
 INTEGER	 	 	 :: nz
 INTEGER	 	 	 :: nosurf
 INTEGER	 	 	 :: mth
 INTEGER	 	 	 :: nj
 INTEGER, DIMENSION(nsg1) 	 :: msing
 INTEGER, DIMENSION(nsg) 	 :: msin
 INTEGER, DIMENSION(nsg1) 	 :: mbm
 INTEGER, DIMENSION(nsg1) 	 :: mbp
 INTEGER	 	 	 :: mbeg
 INTEGER	 	 	 :: mend
 INTEGER	 	 	 :: ms
 INTEGER	 	 	 :: nosing
 INTEGER	 	 	 :: nsing1
 INTEGER	 	 	 :: nsing2
 INTEGER	 	 	 :: nusurf
 INTEGER	 	 	 :: m1
 INTEGER	 	 	 :: m2
 INTEGER	 	 	 :: jmax1
 INTEGER	 	 	 :: jmax2
 INTEGER	 	 	 :: im1
 INTEGER	 	 	 :: im2
 INTEGER	 	 	 :: lpmax
!     endcliche
 END MODULE pstcom

 MODULE comggf
!.....................................................................
!
!     cliche aplet
!
! cliche for generalised green's function method...
!
 USE dims
! 
!
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 alxi1e
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 alxi1o
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 qdchie
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 qdchio
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 ddchie
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 ddchio
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 gm1dco
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 gm1dce
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 x1frbo
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 x1frbe
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 c1frbo
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 c1frbe
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 x1fpro
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 x1fpre
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 c1fpro
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 c1fpre
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 w0l1ou
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 w0l1eu
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 xisolo
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 xisole
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 xix0
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 chi0
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 xilog
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 chlog
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 xicon
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 chcon
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 xix1
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 chi1
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 xix2
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 chi2
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 xsml0
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 csml0
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 xsml1
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 csml1
 COMPLEX*16, DIMENSION(nsg) :: XiLogTerm
 COMPLEX*16, DIMENSION(nsg,nsg) 	 :: w1l1oo
 COMPLEX*16, DIMENSION(nsg,nsg) 	 :: w1l1oe
 COMPLEX*16, DIMENSION(nsg,nsg) 	 :: w1l1eo
 COMPLEX*16, DIMENSION(nsg,nsg) 	 :: w1l1ee
 COMPLEX*16, DIMENSION(2*nsg,2*nsg) 	 :: cjump
 COMPLEX*16, DIMENSION(2*nsg) 	 :: eigold
 LOGICAL	 	 	 :: laplet
 LOGICAL	 	 	 :: lchkap
 LOGICAL	 	 	 :: lsmoth
 LOGICAL	 	 	 :: lplot
 LOGICAL	 	 	 :: lplot2
 LOGICAL	 	 	 :: lzerop
 LOGICAL	 	 	 :: lsymhi
 LOGICAL	 	 	 :: lsmnus
 LOGICAL	 	 	 :: lpack
 LOGICAL	 	 	 :: ltransf
 INTEGER	 	 	 :: packmth
 INTEGER, DIMENSION(nsg1) 	 :: packpts
 REAL*8, DIMENSION(nsg1) 	 :: alfmsh
 REAL*8, DIMENSION(nsg1) 	 :: widmsh
 REAL*8	 	 	 :: betmsh
 REAL*8	 	 	 :: dlayb
 REAL*8	 	 	 :: pquad
 REAL*8, DIMENSION(nsg) 	 :: qslay, qpslay
 REAL*8, DIMENSION(nsg) 	 :: cmatch
 REAL*8, DIMENSION(nsg) 	 :: xsplus
 REAL*8, DIMENSION(nsg) 	 :: xsmnus
 REAL*8, DIMENSION(nsg) 	 :: alay
 REAL*8, DIMENSION(nsg) 	 :: q0lay
 REAL*8, DIMENSION(nsg) 	 :: hlay
 REAL*8, DIMENSION(nsg) 	 :: drlay
 REAL*8, DIMENSION(nsg) 	 :: xmu
 REAL*8, DIMENSION(nsg) 	 :: eflay
 REAL*8, DIMENSION(nsg) 	 :: gelay
 REAL*8, DIMENSION(nsg) 	 :: aklay
 REAL*8, DIMENSION(nsg) 	 :: etalay
 REAL*8, DIMENSION(nsg) 	 :: rholay
 REAL*8, dimension(nsg)    :: rsnorm, beplay, epslay, aqplay
 REAL*8	 	 	 :: slundq
 REAL*8, DIMENSION(nsg) 	 :: alphb
 REAL*8, DIMENSION(nsg) 	 :: alphs
 INTEGER	 	 	 :: outdel
 INTEGER	 	 	 :: outequ
 END MODULE comggf


 MODULE combla
!************************************************************
!c     cliches for inverse iteration subroutines...
!************************************************************
!
!
!
 COMPLEX*16, DIMENSION(:), ALLOCATABLE, target :: 	 a
 COMPLEX*16, DIMENSION(:),pointer :: 	 b
 COMPLEX*16, DIMENSION(:), ALLOCATABLE, target :: 	 x
 COMPLEX*16, DIMENSION(:),pointer :: 	 vb
 COMPLEX*16, DIMENSION(:),pointer  :: 	 c
 COMPLEX*16, DIMENSION(:), ALLOCATABLE :: 	 ut
 COMPLEX*16, DIMENSION(:), ALLOCATABLE :: 	 vt
 COMPLEX*16, DIMENSION(:), ALLOCATABLE :: 	 xt
 COMPLEX*16, DIMENSION(:), ALLOCATABLE :: 	 u
 COMPLEX*16, DIMENSION(:), ALLOCATABLE :: 	 va
!  23.1_r8 .97
!! to sort out  EQUIVALENCE (b,a),(x,vb), (c,x)
!
 INTEGER	 	 	 :: nad2
 INTEGER	 	 	 :: nmtg2
 INTEGER	 	 	 :: nfmgb2
 INTEGER	 	 	 :: nfmga
 INTEGER	 	 	 :: nfmgb
 INTEGER	 	 	 :: ipolp2
 INTEGER	 	 	 :: nfmg
 INTEGER	 	 	 :: nsfg
 INTEGER	 	 	 :: nfeg
 INTEGER	 	 	 :: nfmgc
!     endcliche
 END MODULE combla


 MODULE comivi
!     .........
!
!     cliche comivi
!     macro combla      (pstcom)
 REAL*8	 	 	 :: al0
 REAL*8	 	 	 :: xnorm
 REAL*8	 	 	 :: alam
 INTEGER	 	 	 :: ng
 INTEGER	 	 	 :: mf
 INTEGER	 	 	 :: ml
 INTEGER	 	 	 :: nlong
 INTEGER	 	 	 :: nsing
 INTEGER	 	 	 :: ndlt
 INTEGER	 	 	 :: nda
 INTEGER	 	 	 :: ndb
 INTEGER	 	 	 :: nds
 INTEGER	 	 	 :: ndscrc
 INTEGER	 	 	 :: neg
 INTEGER	 	 	 :: m12
 INTEGER	 	 	 :: m11
 INTEGER	 	 	 :: nit
 INTEGER	 	 	 :: nitmax
 INTEGER	 	 	 :: ndx
 INTEGER	 	 	 :: ndu
 INTEGER	 	 	 :: ndv
 INTEGER	 	 	 :: nconv
 INTEGER	 	 	 :: nlong2
!     endcliche
 END MODULE comivi


 MODULE l22com
!........................
!     cliche l22com
!aplet  27.4_r8 .93
!aplet  27.4_r8 .93
!     common / temstr / vacmat(nfm,nfm)
!
!
!
!.......................................................................
!
!
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 vacmat
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 vacmti
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 scqdp1
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 scqdp2
 LOGICAL	 	 	 :: kinkstab
 LOGICAL	 	 	 :: lfunin
 LOGICAL	 	 	 :: limag
 LOGICAL	 	 	 :: checki
 LOGICAL	 	 	 :: checkd
 LOGICAL	 	 	 :: checkv
 LOGICAL	 	 	 :: symvac
 LOGICAL	 	 	 :: lebc
 LOGICAL	 	 	 :: lscale
 REAL*8	 	 	 :: a
 REAL*8	 	 	 :: b
 REAL*8	 	 	 :: aw
 REAL*8	 	 	 :: bw
 REAL*8                    :: dw
 REAL*8                    :: sw
 REAL*8	 	 	 :: epsq
 REAL*8	 	 	 :: delg
 REAL*8	 	 	 :: xt
 REAL*8	 	 	 :: zt
 REAL*8	 	 	 :: xs
 REAL*8	 	 	 :: zs
 REAL*8	 	 	 :: xtp
 REAL*8	 	 	 :: ztp
 REAL*8	 	 	 :: aval
 REAL*8	 	 	 :: bval
 REAL*8	 	 	 :: a02pyr
 REAL*8	 	 	 :: var
 REAL*8	 	 	 :: varmin
 REAL*8	 	 	 :: varmax
 REAL*8	 	 	 :: alphar
 REAL*8	 	 	 :: betar
 REAL*8	 	 	 :: delr
 REAL*8	 	 	 :: psi0r
 REAL*8	 	 	 :: alphap
 REAL*8	 	 	 :: dlp
 REAL*8	 	 	 :: psi0p
 REAL*8	 	 	 :: gext
 REAL*8	 	 	 :: gext2
 REAL*8	 	 	 :: scale
 INTEGER, DIMENSION(11) 	 :: lstab
 INTEGER	 	 	 :: jtot1
 INTEGER	 	 	 :: jtot2
 INTEGER	 	 	 :: llp
 INTEGER	 	 	 :: nsing
 INTEGER	 	 	 :: nout
 INTEGER	 	 	 :: ieps
 INTEGER	 	 	 :: lmax1
 INTEGER	 	 	 :: nvar
!     endcliche
 END MODULE l22com




 MODULE mtrik1
!     cliche mtrik1
!
!
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 xinf
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 zinf

 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 xwal
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 zwal

!     endcliche
 END MODULE mtrik1
 MODULE newmet
!..........................
!     cliche newmet
!
!      use same common area as cliche metcom in mapping segment
!      so cannot CALL pstboth cliches in same subroutine...
!
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 xa
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 za
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 xdth
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 zdth
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 xdps
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 zdps
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 xjacob
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 xjprym
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 grpssq
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 grpsth
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 xsq
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 xsqdps
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 qdelp
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 delta
!     endcliche
 END MODULE newmet
 MODULE comfft
!
!
 INTEGER, DIMENSION(:), ALLOCATABLE :: 	 invfft
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 sfft
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 omfft
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 onfft
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 gm1
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 ttq
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 ttr
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 geem
! end cliche
 END MODULE comfft
 MODULE l21com
!.............................
!     cliche l21com
!      note: see note below w. r. t. splcof
!
!
 REAL*8, DIMENSION(:,:), ALLOCATABLE :: 	 psi
 REAL*8	 	 	 :: psimax
 REAL*8	 	 	 :: delpsi
 REAL*8	 	 	 :: psival
 REAL*8	 	 	 :: psint
 REAL*8	 	 	 :: psirat
 REAL*8	 	 	 :: psi0
 REAL*8	 	 	 :: tcuro
 REAL*8	 	 	 :: xsc
 REAL*8	 	 	 :: dx
 REAL*8	 	 	 :: dz
 INTEGER	 	 	 :: ncx
 INTEGER	 	 	 :: ncz
 INTEGER	 	 	 :: nxm
 INTEGER	 	 	 :: nzh
 INTEGER	 	 	 :: ndim
 INTEGER	 	 	 :: ndim2
!     endcliche
 END MODULE l21com


 MODULE r33com
!............................
!
!
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 wreg
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 wkin
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 rw1
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 ry2
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 rz3
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 rw4
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 ry5
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 rz6
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 rw7
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 ry9
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 rz10
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 rk1
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 rk4
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 rk7
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 ry8
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 tw
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 ty
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 tz
 COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: 	 td
!     endcliche
 END MODULE r33com


 MODULE l34com
!
!
 COMPLEX*16, DIMENSION(:), ALLOCATABLE :: 	 spot
 COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: 	 amat
 COMPLEX*16, DIMENSION(:), ALLOCATABLE :: 	 skin
!
 COMPLEX*16, DIMENSION(:), ALLOCATABLE :: akin, wpot, ashift
!
!...........................
 END MODULE l34com

 MODULE temps
!     cliche temps
!
!
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 temp1
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 temp2
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 temp3
 REAL*8, DIMENSION(:), ALLOCATABLE :: 	 temp4
!     endcliche
 END MODULE temps

module eigens
  real*8 alam, eigens_epsmac, eigens_epscon, dtry
  integer nda,ndb,ntype,nitmax, nsteps
  logical eigens_lquit, eigens_liter
end module eigens

