      real function rmaxvec(vec,n)
c returns the maximum value in VEC(N)
      real vec(n)
      rmaxvec=0.0
      do 10 i=1,n
         rmaxvec=max(rmaxvec,abs(vec(i)))
 10   continue
      return
      end
