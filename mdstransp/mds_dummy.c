#include <stdio.h>

#ifdef __NOMDSPLUS

int mdsput()
{MdsPut();
return(-1);}
int mdsput_()
{MdsPut();
return(-1);}
int MdsPut()
{
  printf(" __mdsdummy.c:  dummy MdsPut call.\n");
  return(-1);
}
int mdsdisconnect()
{MdsDisconnect();
return(-1);}
int mdsdisconnect_()
{MdsDisconnect();
return(-1);}
int MdsDisconnect()
{
  printf(" __mdsdummy.c:  dummy MdsDisconnect call.\n");
  return(-1);
}
 
int mdsconnect()
{MdsConnect();
return(-1);}
int mdsconnect_()
{MdsConnect();
return(-1);}
int MdsConnect()
{
  printf(" __mdsdummy.c:  dummy MdsConnect call.\n");
  return(-1);
}
int mdsvalue()
{MdsValue();
return(-1);}
int mdsvalue_()
{MdsValue();
return(-1);}
int MdsValue()
{
  printf(" __mdsdummy.c:  dummy MdsValue call.\n");
  return(-1);
}

int mdsopen()
{MdsOpen();
return(-1);}
int mdsopen_()
{MdsOpen();
return(-1);}

int MdsOpen()
{
  printf(" __mdsdummy.c:  dummy MdsOpen call.\n");
  return(-1);
}

int mdsclose()
{MdsClose();
return(-1);}
int mdsclose_()
{MdsClose();
return(-1);}

int MdsClose()
{
  printf(" __mdsdummy.c:  dummy MdsClose call.\n");
  return(-1);
}
#endif
