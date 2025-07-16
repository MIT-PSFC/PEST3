module cgeq_paths
  save
  !  (not needed in xplasma state file -- eq_save/eq_restore)
  !
  !  MDSplus path settings for EFIT trees
  !  this fragment used by geqdsk_mod and geqdsk_subs in xplasma
  !  ...used for resolution of incompatibilities of EFIT trees maintained
  !  on various experimental devices.
  !
  !  dmc 25 Sep 2003
 
  character*50 :: cgeq_apath = ' '
  character*50 :: cgeq_gpath = ' '
  character*10 :: cgeq_nbdry = ' '
 
end module cgeq_paths
