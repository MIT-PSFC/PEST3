/*
PEST3 API calls.
================
File generated on  Wed Feb 28 08:37:25 2001  script  api.py   
*/


/** pstGetError_Gampr **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgeterror_gampr_ pstgeterror_gampr
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgeterror_gampr_ PSTGETERROR_GAMPR
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgeterror_gampr_ (int *kdim,int *ksize,double *err_gampr);


/** pstGetMinorRadius **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetminorradius_ pstgetminorradius
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetminorradius_ PSTGETMINORRADIUS
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetminorradius_ (double *p);


/** pstSetIsolver **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstsetisolver_ pstsetisolver
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstsetisolver_ PSTSETISOLVER
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstsetisolver_ (int *k);


/** pstSetBetmsh **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstsetbetmsh_ pstsetbetmsh
#define pstsetalfmsh_ pstsetalfmsh
#define pstsetwidmsh_ pstsetwidmsh
#define pstsetpackmth_ pstsetpackmth
#define pstsetpackpts_ pstsetpackpts
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstsetbetmsh_ PSTSETBETMSH
#define pstsetalfmsh_ PSTSETALFMSH
#define pstsetwidmsh_ PSTSETWIDMSH
#define pstsetpackmth_ PSTSETPACKMTH
#define pstsetpackpts_ PSTSETPACKPTS
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstsetbetmsh_ (double *p);
void  pstsetalfmsh_ (int *j,double *p);
void  pstsetwidmsh_ (int *j,double *p);
void  pstsetpackmth_ (int *j);
void  pstsetpackpts_ (int *j,int *p);


/** pstSetMSIN **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstsetmsin_ pstsetmsin
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstsetmsin_ PSTSETMSIN
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstsetmsin_ (int *jdim,int *j);


/** pstSetDlayb **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstsetdlayb_ pstsetdlayb
#define pstsetdlay_ pstsetdlay
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstsetdlayb_ PSTSETDLAYB
#define pstsetdlay_ PSTSETDLAY
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstsetdlayb_ (double *p);
void  pstsetdlay_ (double *p);


/** pstGetMth **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetmth_ pstgetmth
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetmth_ PSTGETMTH
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetmth_ (int *j);


/** pstGetNosurf **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetnosurf_ pstgetnosurf
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetnosurf_ PSTGETNOSURF
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetnosurf_ (int *j);


/** pstGetXiLogTerm **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetxilogterm_ pstgetxilogterm
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetxilogterm_ PSTGETXILOGTERM
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetxilogterm_ (int *ksize,double *re_p,double *im_p);


/** pstGetRsnorm **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetrsnorm_ pstgetrsnorm
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetrsnorm_ PSTGETRSNORM
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetrsnorm_ (int *ksize,double *pnorm);


/** pstInit **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstinit_ pstinit
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstinit_ PSTINIT
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstinit_ (int *kmax1,int *kmeshMax,int *kinputFormat,char *kinputPath,
		int *kin_len,double *ptime,double *newQ,int *akima,
		double *x2axis,double *x2edge,int *np1,int *nr,double *separatrix, 
		int *ier, long long_len8);


/** pstGetfourPiInductance **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetfourpiinductance_ pstgetfourpiinductance
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetfourpiinductance_ PSTGETFOURPIINDUCTANCE
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetfourpiinductance_ (double *p);


/** pstGetTriangularity **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgettriangularity_ pstgettriangularity
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgettriangularity_ PSTGETTRIANGULARITY
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgettriangularity_ (double *p);


/** pstGetHlay **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgethlay_ pstgethlay
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgethlay_ PSTGETHLAY
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgethlay_ (int *ksize,double *p);


/** pstGetGampr **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetgampr_ pstgetgampr
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetgampr_ PSTGETGAMPR
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetgampr_ (int *kdim,int *ksize,double *re_gampr,double *im_gampr);


/** pstGetbetaPoloidal **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetbetapoloidal_ pstgetbetapoloidal
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetbetapoloidal_ PSTGETBETAPOLOIDAL
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetbetapoloidal_ (double *p);


/** pstGetCmatch **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetcmatch_ pstgetcmatch
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetcmatch_ PSTGETCMATCH
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetcmatch_ (int *ksize,double *pmatch);


/** pstGetEflay **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgeteflay_ pstgeteflay
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgeteflay_ PSTGETEFLAY
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgeteflay_ (int *ksize,double *p);


/** pstGetPsisin **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetpsisin_ pstgetpsisin
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetpsisin_ PSTGETPSISIN
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetpsisin_ (int *ksize,double *ppsi0,double *ppsis,double *ppsia);


/** pstSetWALL **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstsetwall_ pstsetwall
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstsetwall_ PSTSETWALL
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstsetwall_ (int *j);


/** pstSetLsymhi **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstsetlsymhi_ pstsetlsymhi
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstsetlsymhi_ PSTSETLSYMHI
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstsetlsymhi_ (int *k);


/** pstGetrMagnetic **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetrmagnetic_ pstgetrmagnetic
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetrmagnetic_ PSTGETRMAGNETIC
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetrmagnetic_ (double *p);


/** pstGetElongation **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetelongation_ pstgetelongation
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetelongation_ PSTGETELONGATION
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetelongation_ (double *p);


/** pstGetQslay **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetqslay_ pstgetqslay
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetqslay_ PSTGETQSLAY
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetqslay_ (int *ksize,double *p);


/** pstGetQBounds **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetqbounds_ pstgetqbounds
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetqbounds_ PSTGETQBOUNDS
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetqbounds_ (double *q0,double *qa,double *qmin,double *qmax);


/** pstGetTroyonG **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgettroyong_ pstgettroyong
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgettroyong_ PSTGETTROYONG
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgettroyong_ (double *p);


/** pstGetNosing **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetnosing_ pstgetnosing
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetnosing_ PSTGETNOSING
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetnosing_ (int *josing);


/** pstGetXmu **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetxmu_ pstgetxmu
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetxmu_ PSTGETXMU
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetxmu_ (int *ksize,double *pmu);


/** pstGetGSAverageError **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetgsaverageerror_ pstgetgsaverageerror
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetgsaverageerror_ PSTGETGSAVERAGEERROR
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetgsaverageerror_ (double *res);


/** pstSetUniformMesh **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstsetuniformmesh_ pstsetuniformmesh
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstsetuniformmesh_ PSTSETUNIFORMMESH
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstsetuniformmesh_ (int *k);


/** pstSetB **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstsetb_ pstsetb
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstsetb_ PSTSETB
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstsetb_ (double *pb);


/** pstGetXsmnus **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetxsmnus_ pstgetxsmnus
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetxsmnus_ PSTGETXSMNUS
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetxsmnus_ (int *ksize,double *psmnus);


/** pstSetMM **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstsetmm_ pstsetmm
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstsetmm_ PSTSETMM
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstsetmm_ (int *jdim,int *j);


/** pstSetN **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstsetn_ pstsetn
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstsetn_ PSTSETN
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstsetn_ (int *j);


/** pstSetINFWALL **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstsetinfwall_ pstsetinfwall
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstsetinfwall_ PSTSETINFWALL
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstsetinfwall_ (int *j);


/** pstGetMajorRadius **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetmajorradius_ pstgetmajorradius
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetmajorradius_ PSTGETMAJORRADIUS
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetmajorradius_ (double *p);


/** pstGetb0SquareCentre **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetb0squarecentre_ pstgetb0squarecentre
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetb0squarecentre_ PSTGETB0SQUARECENTRE
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetb0squarecentre_ (double *p);


/** pstGetIneg **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetineg_ pstgetineg
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetineg_ PSTGETINEG
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetineg_ (int *kneg,int *nodes,int *nodes_size);


/** pstGetError_Delpr **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgeterror_delpr_ pstgeterror_delpr
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgeterror_delpr_ PSTGETERROR_DELPR
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgeterror_delpr_ (int *kdim,int *ksize,double *err_delpr);


/** pstGetbetaToroidal **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetbetatoroidal_ pstgetbetatoroidal
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetbetatoroidal_ PSTGETBETATOROIDAL
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetbetatoroidal_ (double *p);


/** pstFree **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstfree_ pstfree
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstfree_ PSTFREE
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstfree_ (void);


/** pstGetWlambda **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetwlambda_ pstgetwlambda
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetwlambda_ PSTGETWLAMBDA
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetwlambda_ (double *p);


/** pstGetDelpr **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetdelpr_ pstgetdelpr
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetdelpr_ PSTGETDELPR
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetdelpr_ (int *kdim,int *ksize,double *re_delpr,double *im_delpr);


/** pstExec **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstexec_ pstexec
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstexec_ PSTEXEC
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstexec_ (void);


/** pstGettotalToroidalCurrent **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgettotaltoroidalcurrent_ pstgettotaltoroidalcurrent
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgettotaltoroidalcurrent_ PSTGETTOTALTOROIDALCURRENT
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgettotaltoroidalcurrent_ (double *p);


/** pstGetDrlay **/


#if __VMS || __IBM || __RS6000 || __HP || __ABS
#define pstgetdrlay_ pstgetdrlay
#endif /* __VMS || __IBM || __RS6000 || __HP || __ABS */

#ifdef __CRAY
#define pstgetdrlay_ PSTGETDRLAY
#endif /* __CRAY */

/* Default assume underscore appended */
void  pstgetdrlay_ (int *ksize,double *p);
