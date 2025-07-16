subroutine pstWresult
  !
  ! write results into pest3.nc file
  !
  USE pstcom
  USE newmet
  USE comggf

  USE ezcdf
  implicit none
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
  CHARACTER*10, parameter :: FILENAME='pest3.nc'
  integer cdfid
  character*3, parameter :: xint = 'INT'
  character*2, parameter :: xr8 = 'R8'
  integer dimlens(3)
  real(r8) zza(nosing, nosing)
  real(r8) zzb(nusurf,jmax1,nosing)

  CALL cdfOpn(cdfid, filename, 'w')

  dimlens = (/ 1, 1, 1 /)
  ! mapper
  call cdfDefVar(cdfid, "nsf0", dimlens, xint)
  call cdfDefVar(cdfid, "nths0", dimlens, xint)

  ! global
  call cdfDefVar(cdfid, "nosing", dimlens, xint)
  call cdfDefVar(cdfid, "jmax1", dimlens, xint)
  call cdfDefVar(cdfid, "nusurf", dimlens, xint)
  call cdfDefVar(cdfid, "mp", dimlens, xint)
  call cdfDefVar(cdfid, "mth1", dimlens, xint)
  call cdfDefVar(cdfid, "nIdealInstabilities", dimlens, xint)

  ! geometry
  call cdfDefVar(cdfid, "inverseAspectRatio", dimlens, xr8)
  call cdfDefVar(cdfid, "rMagnetic", dimlens, xr8)
  call cdfDefVar(cdfid, "rCentre", dimlens, xr8)
  call cdfDefVar(cdfid, "elongation", dimlens, xr8)
  call cdfDefVar(cdfid, "triangularity", dimlens, xr8)
  call cdfDefVar(cdfid, "plasmaSurface", dimlens, xr8)
  call cdfDefVar(cdfid, "plasmaVolume", dimlens, xr8)

  ! current, B, betas etc
  call cdfDefVar(cdfid, "totalToroidalCurrent", dimlens, xr8)
  call cdfDefVar(cdfid, "b0SquareCentre", dimlens, xr8)
  call cdfDefVar(cdfid, "betaPoloidal", dimlens, xr8)
  call cdfDefVar(cdfid, "betaToroidal", dimlens, xr8)
  call cdfDefVar(cdfid, "betaStar", dimlens, xr8)
  call cdfDefVar(cdfid, "TroyonG", dimlens, xr8)
  call cdfDefVar(cdfid, "gStar", dimlens, xr8)
  call cdfDefVar(cdfid, "volumeAveragePressure", dimlens, xr8)
  call cdfDefVar(cdfid, "volumeAverageBSquare", dimlens, xr8)
  call cdfDefVar(cdfid, "volumeAverageBeta", dimlens, xr8)
  call cdfDefVar(cdfid, "fourPiInductance", dimlens, xr8)
  call cdfDefVar(cdfid, "pressurePeakingFactor", dimlens, xr8)


  call cdfDefVar(cdfid, "psimin", dimlens, xr8)
  call cdfDefVar(cdfid, "psimax", dimlens, xr8)
  call cdfDefVar(cdfid, "psi0" , dimlens, xr8)
  call cdfDefVar(cdfid, "psia" , dimlens, xr8)


  ! layer
  call cdfDefVar(cdfid, "n", dimlens, xr8)
  dimlens = (/ nosing, 1, 1 /)
  call cdfDefVar(cdfid, "qslay", dimlens, xr8)
  call cdfDefVar(cdfid, "qpslay", dimlens, xr8)
  call cdfDefVar(cdfid, "psisin", dimlens, xr8)
  call cdfDefVar(cdfid, "drlay", dimlens, xr8)
  call cdfDefVar(cdfid, "eflay", dimlens, xr8)
  call cdfDefVar(cdfid, "gelay", dimlens, xr8)
  call cdfDefVar(cdfid, "hlay", dimlens, xr8)
  call cdfDefVar(cdfid, "aklay", dimlens, xr8)
  call cdfDefVar(cdfid, "xmu", dimlens, xr8)
  call cdfDefVar(cdfid, "cmatch", dimlens, xr8)
  call cdfDefVar(cdfid, "rsnorm", dimlens, xr8)
  call cdfDefVar(cdfid, "beplay", dimlens, xr8)
  call cdfDefVar(cdfid, "epslay", dimlens, xr8)
  call cdfDefVar(cdfid, "aqplay", dimlens, xr8)
  call cdfDefVar(cdfid, "xsmnus", dimlens, xr8)
  call cdfDefVar(cdfid, "xsplus", dimlens, xr8)

  dimlens = (/ nosing, nosing, 1 /)
  call cdfDefVar(cdfid, "aprim_re", dimlens, xr8)
  call cdfDefVar(cdfid, "bprim_re", dimlens, xr8)
  call cdfDefVar(cdfid, "gprim_re", dimlens, xr8)
  call cdfDefVar(cdfid, "dprim_re", dimlens, xr8)
  call cdfDefVar(cdfid, "aprim_im", dimlens, xr8)
  call cdfDefVar(cdfid, "bprim_im", dimlens, xr8)
  call cdfDefVar(cdfid, "gprim_im", dimlens, xr8)
  call cdfDefVar(cdfid, "dprim_im", dimlens, xr8)
  call cdfDefVar(cdfid, "error_aprim", dimlens, xr8)
  call cdfDefVar(cdfid, "error_bprim", dimlens, xr8)
  call cdfDefVar(cdfid, "error_gprim", dimlens, xr8)
  call cdfDefVar(cdfid, "error_dprim", dimlens, xr8)

  ! mesh & profiles
  dimlens = (/ mp, 1, 1 /)
  call cdfDefVar(cdfid, "psinod", dimlens, xr8)
  dimlens = (/ nusurf, 1, 1 /)
  call cdfDefVar(cdfid, "pa", dimlens, xr8)
  call cdfDefVar(cdfid, "ppa", dimlens, xr8)
  call cdfDefVar(cdfid, "qa", dimlens, xr8)
  call cdfDefVar(cdfid, "qpa", dimlens, xr8)
  call cdfDefVar(cdfid, "ga", dimlens, xr8)
  call cdfDefVar(cdfid, "gpa", dimlens, xr8)

  !grid
  dimlens = (/ nusurf, 1, 1 /)
  call cdfDefVar(cdfid, "psinew", dimlens, xr8)
  call cdfDefVAr(cdfid, "di", dimlens, xr8) ! ideal Mercier
  call cdfDefVAr(cdfid, "dr", dimlens, xr8) ! resistive Mercier

  dimlens = (/ mth1, nusurf, 1 /)
  call cdfDefVar(cdfid, "xa", dimlens, xr8)
  call cdfDefVar(cdfid, "za", dimlens, xr8)
  call cdfDefVar(cdfid, "xjacob", dimlens, xr8)
  call cdfDefVar(cdfid, "grpssq", dimlens, xr8)
  call cdfDefVar(cdfid, "delta", dimlens, xr8)

  ! big solution
  dimlens = (/ nusurf, jmax1, nosing /)
  call cdfDefVar(cdfid, "x1frbo_re", dimlens, xr8)
  call cdfDefVar(cdfid, "x1frbe_re", dimlens, xr8)
  call cdfDefVar(cdfid, "x1frbo_im", dimlens, xr8)
  call cdfDefVar(cdfid, "x1frbe_im", dimlens, xr8)

  ! small solution
  dimlens = (/ mp, jmax1, nosing /)
  call cdfDefVar(cdfid, "xisolo_re", dimlens, xr8)
  call cdfDefVar(cdfid, "xisole_re", dimlens, xr8)
  call cdfDefVar(cdfid, "xisolo_im", dimlens, xr8)
  call cdfDefVar(cdfid, "xisole_im", dimlens, xr8)

!!!!!!!!!!!!!! dump data !!!!!!!!!!!!!!!!!!!!!!

  call cdfPutVar(cdfid, "nsf0", nsf0)
  call cdfPutVar(cdfid, "nths0", nths0)
  call cdfPutVar(cdfid, "nosing", nosing)
  call cdfPutVar(cdfid, "jmax1", jmax1)
  call cdfPutVar(cdfid, "nusurf", nusurf)
  call cdfPutVar(cdfid, "mp", mp)
  call cdfPutVar(cdfid, "mth1", mth1)
  call cdfPutVar(cdfid, "nIdealInstabilities", nIdealInstabilities)

  call cdfPutVar(cdfid, "inverseAspectRatio", inverseAspectRatio)
  call cdfPutVar(cdfid, "rMagnetic", rMagnetic)
  call cdfPutVar(cdfid, "rCentre", rCentre)
  call cdfPutVar(cdfid, "elongation",elongation )
  call cdfPutVar(cdfid, "triangularity", triangularity)
  call cdfPutVar(cdfid, "plasmaSurface",plasmaSurface )
  call cdfPutVar(cdfid, "plasmaVolume",plasmaVolume )

  ! current, B, betas etc
  call cdfPutVar(cdfid, "totalToroidalCurrent", totalToroidalCurrent)
  call cdfPutVar(cdfid, "b0SquareCentre",b0SquareCentre )
  call cdfPutVar(cdfid, "betaPoloidal", betaPoloidal)
  call cdfPutVar(cdfid, "betaToroidal",betaToroidal )
  call cdfPutVar(cdfid, "betaStar",betaStar )
  call cdfPutVar(cdfid, "TroyonG", TroyonG)
  call cdfPutVar(cdfid, "gStar",gStar )
  call cdfPutVar(cdfid, "volumeAveragePressure", volumeAveragePressure)
  call cdfPutVar(cdfid, "volumeAverageBSquare",volumeAverageBSquare )
  call cdfPutVar(cdfid, "volumeAverageBeta",volumeAverageBeta )
  call cdfPutVar(cdfid, "fourPiInductance",fourPiInductance )
  call cdfPutVar(cdfid, "pressurePeakingFactor", pressurePeakingFactor)

  call cdfPutVar(cdfid, "psimin", psinod(1))
  call cdfPutVar(cdfid, "psimax", psinod(mp))
  call cdfPutVar(cdfid, "psi0", psisin(1))
  call cdfPutVar(cdfid, "psia", psisin(nosing+2))
  call cdfPutVar(cdfid, "n", n)

  if (nosing > 1) then
     call cdfPutVar(cdfid, "qslay", qslay)
     call cdfPutVar(cdfid, "qpslay", qpslay)
     call cdfPutVar(cdfid, "psisin", psisin(2:nsing1))
     call cdfPutVar(cdfid, "drlay", drlay)
     call cdfPutVar(cdfid, "eflay", eflay)
     call cdfPutVar(cdfid, "gelay", gelay)
     call cdfPutVar(cdfid, "hlay", hlay)
     call cdfPutVar(cdfid, "aklay", aklay)
     call cdfPutVar(cdfid, "xmu", xmu)
     call cdfPutVar(cdfid, "cmatch", cmatch)
     call cdfPutVar(cdfid, "rsnorm", rsnorm)
     call cdfPutVar(cdfid, "beplay", beplay)
     call cdfPutVar(cdfid, "epslay", epslay)
     call cdfPutVar(cdfid, "aqplay", aqplay)
     call cdfPutVar(cdfid, "xsmnus", xsmnus)
     call cdfPutVar(cdfid, "xsplus", xsplus)
     zza = real( apr(1:nosing,1:nosing) , r8)
     call cdfPutVar(cdfid, "aprim_re", zza)
     zza = real( bpr(1:nosing,1:nosing) , r8)
     call cdfPutVar(cdfid, "bprim_re", zza)
     zza = real( gampr(1:nosing,1:nosing) , r8)
     call cdfPutVar(cdfid, "gprim_re", zza)
     zza = real( delpr(1:nosing,1:nosing) , r8)
     call cdfPutVar(cdfid, "dprim_re", zza)
     zza = aimag( apr(1:nosing,1:nosing))
     call cdfPutVar(cdfid, "aprim_im", zza)
     zza = aimag( bpr(1:nosing,1:nosing))
     call cdfPutVar(cdfid, "bprim_im", zza)
     zza = aimag( gampr(1:nosing,1:nosing))
     call cdfPutVar(cdfid, "gprim_im", zza)
     zza = aimag( delpr(1:nosing,1:nosing))
     call cdfPutVar(cdfid, "dprim_im", zza)

     zza = error_apr(1:nosing,1:nosing)
     call cdfPutVar(cdfid, "error_aprim", zza)
     zza = error_bpr(1:nosing,1:nosing)
     call cdfPutVar(cdfid, "error_bprim", zza)
     zza = error_gampr(1:nosing,1:nosing)
     call cdfPutVar(cdfid, "error_gprim", zza)
     zza = error_delpr(1:nosing,1:nosing)
     call cdfPutVar(cdfid, "error_dprim", zza)
  else
     call cdfPutVar(cdfid, "qslay", qslay(1))
     call cdfPutVar(cdfid, "qpslay", qpslay(1))
     call cdfPutVar(cdfid, "psisin", psisin(2))
     call cdfPutVar(cdfid, "drlay", drlay(1))
     call cdfPutVar(cdfid, "eflay", eflay(1))
     call cdfPutVar(cdfid, "gelay", gelay(1))
     call cdfPutVar(cdfid, "hlay", hlay(1))
     call cdfPutVar(cdfid, "aklay", aklay(1))
     call cdfPutVar(cdfid, "xmu", xmu(1))
     call cdfPutVar(cdfid, "cmatch", cmatch(1))
     call cdfPutVar(cdfid, "rsnorm", rsnorm(1))
     call cdfPutVar(cdfid, "beplay", beplay(1))
     call cdfPutVar(cdfid, "epslay", epslay(1))
     call cdfPutVar(cdfid, "aqplay", aqplay(1))
     call cdfPutVar(cdfid, "xsmnus", xsmnus(1))
     call cdfPutVar(cdfid, "xsplus", xsplus(1))
     zza = real( apr(1:nosing,1:nosing) , r8)
     call cdfPutVar(cdfid, "aprim_re", zza(1,1))
     zza = real( bpr(1:nosing,1:nosing) , r8)
     call cdfPutVar(cdfid, "bprim_re", zza(1,1))
     zza = real( gampr(1:nosing,1:nosing) , r8)
     call cdfPutVar(cdfid, "gprim_re", zza(1,1))
     zza = real( delpr(1:nosing,1:nosing) , r8)
     call cdfPutVar(cdfid, "dprim_re", zza(1,1))
     zza = aimag( apr(1:nosing,1:nosing))
     call cdfPutVar(cdfid, "aprim_im", zza(1,1))
     zza = aimag( bpr(1:nosing,1:nosing))
     call cdfPutVar(cdfid, "bprim_im", zza(1,1))
     zza = aimag( gampr(1:nosing,1:nosing))
     call cdfPutVar(cdfid, "gprim_im", zza(1,1))
     zza = aimag( delpr(1:nosing,1:nosing))
     call cdfPutVar(cdfid, "dprim_im", zza(1,1))

     zza = error_apr(1:nosing,1:nosing)
     call cdfPutVar(cdfid, "error_aprim", zza(1,1))
     zza = error_bpr(1:nosing,1:nosing)
     call cdfPutVar(cdfid, "error_bprim", zza(1,1))
     zza = error_gampr(1:nosing,1:nosing)
     call cdfPutVar(cdfid, "error_gprim", zza(1,1))
     zza = error_delpr(1:nosing,1:nosing)
     call cdfPutVar(cdfid, "error_dprim", zza(1,1))
  endif

  call cdfPutVar(cdfid, "psinod", psinod(1:mp))
  call cdfPutVar(cdfid, "pa",  pa(1:nusurf))
  call cdfPutVar(cdfid, "ppa", ppa(1:nusurf))
  call cdfPutVar(cdfid, "qa",  qa(1:nusurf))
  call cdfPutVar(cdfid, "qpa", qpa(1:nusurf))
  call cdfPutVar(cdfid, "ga",  ga(1:nusurf))
  call cdfPutVar(cdfid, "gpa", gpa(1:nusurf))
  call cdfPutVar(cdfid, "psinew", psinew(1:nusurf))
  call cdfPutVar(cdfid, "di", di(1:nusurf))
  call cdfPutVar(cdfid, "dr", dr(1:nusurf))
  call cdfPutVar(cdfid, "xa", xa(1:mth1,1:nusurf))
  call cdfPutVar(cdfid, "za", za(1:mth1,1:nusurf))
  call cdfPutVar(cdfid, "xjacob", xjacob(1:mth1,1:nusurf))
  call cdfPutVar(cdfid, "grpssq", grpssq(1:mth1,1:nusurf))
  call cdfPutVar(cdfid, "delta", delta(1:mth1,1:nusurf))
  if(nosing == 1) then ! potential problem if jmax1==1
     zzb =   real(x1frbo(1:nusurf,1:jmax1,1:nosing),r8)
     call cdfPutVar(cdfid, "x1frbo_re",zzb(1:nusurf,1:jmax1,1))
     zzb =   real(x1frbe(1:nusurf,1:jmax1,1:nosing),r8)
     call cdfPutVar(cdfid, "x1frbe_re",zzb(1:nusurf,1:jmax1,1))
     zzb =  aimag(x1frbo(1:nusurf,1:jmax1,1:nosing))
     call cdfPutVar(cdfid, "x1frbo_im",zzb(1:nusurf,1:jmax1,1))
     zzb =  aimag(x1frbe(1:nusurf,1:jmax1,1:nosing))
     call cdfPutVar(cdfid, "x1frbe_im",zzb(1:nusurf,1:jmax1,1))
     zzb(1:mp,1:jmax1,1:nosing) =   real(xisolo(1:mp,1:jmax1,1:nosing),r8)
     call cdfPutVar(cdfid, "xisolo_re",zzb(1:mp,1:jmax1,1))
     zzb(1:mp,1:jmax1,1:nosing) =   real(xisole(1:mp,1:jmax1,1:nosing),r8)
     call cdfPutVar(cdfid, "xisole_re",zzb(1:mp,1:jmax1,1))
     zzb(1:mp,1:jmax1,1:nosing) =  aimag(xisolo(1:mp,1:jmax1,1:nosing))
     call cdfPutVar(cdfid, "xisolo_im",zzb(1:mp,1:jmax1,1))
     zzb(1:mp,1:jmax1,1:nosing) =  aimag(xisole(1:mp,1:jmax1,1:nosing))
     call cdfPutVar(cdfid, "xisole_im",zzb(1:mp,1:jmax1,1))
  else
     zzb =   real(x1frbo(1:nusurf,1:jmax1,1:nosing),r8)
     call cdfPutVar(cdfid, "x1frbo_re",zzb(1:nusurf,1:jmax1,1:nosing))
     zzb =   real(x1frbe(1:nusurf,1:jmax1,1:nosing),r8)
     call cdfPutVar(cdfid, "x1frbe_re",zzb(1:nusurf,1:jmax1,1:nosing))
     zzb =  aimag(x1frbo(1:nusurf,1:jmax1,1:nosing))
     call cdfPutVar(cdfid, "x1frbo_im",zzb(1:nusurf,1:jmax1,1:nosing))
     zzb =  aimag(x1frbe(1:nusurf,1:jmax1,1:nosing))
     call cdfPutVar(cdfid, "x1frbe_im",zzb(1:nusurf,1:jmax1,1:nosing))
     zzb(1:mp,1:jmax1,1:nosing) =   real(xisolo(1:mp,1:jmax1,1:nosing),r8)
     call cdfPutVar(cdfid, "xisolo_re",zzb(1:mp,1:jmax1,1:nosing))
     zzb(1:mp,1:jmax1,1:nosing) =   real(xisole(1:mp,1:jmax1,1:nosing),r8)
     call cdfPutVar(cdfid, "xisole_re",zzb(1:mp,1:jmax1,1:nosing))
     zzb(1:mp,1:jmax1,1:nosing) =  aimag(xisolo(1:mp,1:jmax1,1:nosing))
     call cdfPutVar(cdfid, "xisolo_im",zzb(1:mp,1:jmax1,1:nosing))
     zzb(1:mp,1:jmax1,1:nosing) =  aimag(xisole(1:mp,1:jmax1,1:nosing))
     call cdfPutVar(cdfid, "xisole_im",zzb(1:mp,1:jmax1,1:nosing))
  endif



  call cdfCls(cdfid)

  return
end subroutine pstWresult
