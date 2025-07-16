subroutine geqdsk_lunmsg_set(ilun)

  use geqdsk_aux

  lunmsg = ilun

end subroutine geqdsk_lunmsg_set


subroutine geqdsk_lunmsg_get(ilun)

  use geqdsk_aux

  ilun = lunmsg

end subroutine geqdsk_lunmsg_get
