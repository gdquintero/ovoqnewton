c Main code to test subroutine sysredneu
c June 6, 2023

      implicit none
      integer nnx, nnobs, nrowsa, ndepth,nnAb,nnobs2, sizeb
      integer dimaat
      parameter(nnx=10, nnobs=3200,nrowsa=100, ndepth=10)
c     parameter (nnAb=(nnx+1)*nrowsa + (ndepth-1)*nrowsa*nrowsa)
      parameter (nnAb=11410)
      parameter (dimaat = min0(nnobs,nnAb))
      parameter(nnobs2 = nnobs*2)
      parameter(sizeb = nnobs) 
      integer nobs
      integer nx, depth, ncomps(ndepth), npsi
      integer nAb
      integer ncoan,nada
      double precision yaux1(nrowsa), yaux2(nrowsa), zaux(nrowsa,nrowsa)
      double precision Ab(nnAb), y(nnobs), f(nnobs),vaux(2*nnobs)
      double precision jacob(nnobs*nnAb),jacosave(nnobs*nnAb)
      double precision fsave(nnobs),Abtrial(nnAb),ftrial(nnobs)      
      double precision jaca(nrowsa,nnAb), jacx(nnAb,nnAb)
      double precision aat(dimaat*dimaat), aatsave(dimaat*dimaat)
      double precision jacaux(nrowsa,nnAb)
      double precision a(nrowsa, nrowsa, ndepth), b(nrowsa, ndepth)
      double precision x(nnx, nnobs)
      double precision mataux(nnobs2)
      double precision seed
      integer jqua,nback
      integer iperm(nnobs)
      integer i, j, k, iobs, idep, kprob
      integer kon, max, nvar, nref
      integer idata, iout
      integer iadd, mgsyes
      real day1, day2
      double precision eps, delta,rmsdtrial, perout 
      double precision zmedio
      double precision z, under
      integer sorte, iqn, maxrow, ipri, inewran, option
      integer ia, ja, ka, ndatos
      double precision alfa, beta, gama
      integer esprin
      double precision rmsd 
      logical ceroum, lovo, elimi
      integer plovo, nout
      real start, finish
c  Data that will be read from fort.106
      double precision datos(4,6060)

      write(*,*)' If you have a LOVO problem type 1'
      i = 1 
      if(i.eq.1) then
      lovo = .true.
      else
      lovo = .false.
      endif

      if(lovo) then
      write(*,*)' Number of observations:'
      nobs = 100 
      write(*,*)' Number of outliers:'
      nout = 0
      write(*,*)' Percentage of error of each outlier:'
      perout = 0
      write(*,*)' p parameter of LOVO:'
      plovo = 50

      endif

      write(*,*)' If want to eliminate auxiliar printints type 1'
      write (*,*) ' Default = 1'
      i=1
      if(i.eq.1) then
      elimi = .true.
      else
      elimi = .false.
      endif      

      write(*, *)' Training Neural Networks with sysredneu'
      write(*, *) 
      write(*, *)' Memory is reserved for supporting ', nnAb, 
     *   ' coefficients '
      


 


      write(*, *)' If the independent terms y_i are only 0 or 1 type 1'
      write(*,*)' Default: 0'
      i = 0
      ceroum = .false.
      
      if(i.eq.1) ceroum = .true.
 
      seed = 28383723771717.d0
      write(*, *)' Number of function psi :'
      write(*,*)' psi = 3 is RELU except last level'
      write(*,*)' Default=3'
      npsi=3

      write(*, *)' If want detailed printing type 1, else type 0'
      write(*,*)' Default = 0'
      ipri = 0
      if(ipri.eq.1) esprin=1
      if(ipri.eq.0) then
      write(*, *)' Essential printing every how many iterations?'
      write(*,*)' Default = 1'
      esprin=1
      endif

      write(*, *)' If you want to use and detect rank with MGS, type 1'
      write(*, *)' otherwise type 0. Default: 0'
      mgsyes = 0       


      if(lovo) then
      maxrow = nobs
      else

      write(*, *)' Type Maximal number of observations handled'
      write(*, *)' at each iteration: Default=number of observations'
      read(*, *) maxrow
      write(*, *)' Maximal number of observations handled'
      write(*, *)' at each iteration', maxrow
      write(*, *)
      write(*, *)' (If ', maxrow,' is bigger than the number of'
      write(*, *)' observations, it will be decreased by the code)' 
      endif
c  This endif corresponds to "if(lovo) then maxrow = nobs"

c      write(*, *)' Warning!'
c      write(*, *)' In this code parameter sizeb is ', sizeb
c      write(*, *)' Optimal time execution is obtained if sizeb =',maxrow
c      write(*, *)

      write(*, *)
      write(*, *)' Under-relaxation parameter:'
      write(*, *)' Default = 1'
      under = 1
      write(*, *)
      write(*, *)' Under-relaxation parameter:', under

      write(*, *)
      write(*, *)' Number refinements (recall that 0 is "pure newton"):'
      write(*, *)' Default=0'
      nref=0
      write(*, *)
      write(*, *)' Refinements following each Newton iteration:', nref
      write(*, *)

      if(nref.gt.0) then
      iqn = 0
      write(*, *)' If want quasi-Newton refinements type 1, else type 0'
      write(*,*)' Default=0'
      iqn=0
      write(*, *)
      if(iqn.eq.1) then
      write(*, *)' Refinements will be of Quasi-Newton (secant) type'
      else
      write(*, *)' Refinements will repeat the previous Jacobian'
      endif
      endif
      write(*, *)


      write(*, *)' Maximal number of iterations, Default=100:'
      max=100
      write(*, *)
      write(*, *)' Maximal number of iterations allowed:', max
      write(*, *)      


      write(*, *)' Stopping delta regarding the norm of increment:'
      write(*,*)' Default = 0'
      delta = 0.
      write(*, *)' "Small" delta used for stopping by increment norm:',
     *              delta
      write(*, *)' Stopping epsilon regarding RMSD of residuals:'
      write(*,*)' Default = 0.0001'
      eps = 1.d-4
      write(*, *)
      write(*, *)' "Small" epsilon used for RMSD residual stopping:',eps

      write(*, *)
  
 

      write(*, *)' Data generation:'
      write(*, *)' If you want Saint Venant given in fort.106, type 1'
      write(*, *)' fort.106 contains 3060 observations, corresponding'
      write(*, *)' to 30 days'
      write(*, *)' If you want Saint Venant given in fort.107, type 2'
      write(*, *)' fort.107 contains 1442 observations, corresponding'
      write(*, *)' to 30 days'
 

      write(*, *)' Otherwise, type 3'
      write(*,*)' Default = 3'

      kprob=3
      if(kprob.eq.1.or.kprob.eq.2) then

c  Saint Venant synthetic from fort.106
      if(kprob.eq.1) then
      ndatos = 6060
      nvar = 3
     
      do iobs = 1, ndatos
      read(106, *) (datos(i, iobs),i=1,4)
      if(ceroum) then
      if(datos(4,iobs).gt.7.5d0) then
      datos(4,iobs) = 1.d0
      else
      datos(4,iobs) = 0.d0
      endif
      endif
      end do
      endif

      if(kprob.eq.2) then
      ndatos = 1442
      nvar = 3
     
      do iobs = 1, ndatos
      read(107, *) (datos(i, iobs),i=1,4)
      if(ceroum) then
      if(datos(4,iobs).gt.7.5d0) then
      datos(4,iobs) = 1.d0
      else
      datos(4,iobs) = 0.d0
      endif
      endif
      end do
      endif
 




      write(*, *)
      write(*, *)' Observations given ' 
      write(*, *)' in the following order:'
      write(*, *)' 1:Inlet discharge(m^3/s),2:station x(km)',
     * ' 3:time(days)' 
      write(*, *)
  
      if(.not.lovo) then
      write(*, *)' Number of observations nobs:'

      read(*, *) nobs 
      endif      
      if(nx.gt.nnx.or.nobs.gt.nnobs) then
      write(*, *)' Abort. Maximal nx and nobs are:', nnx, nnobs
      stop
      endif
      if(2*nobs.gt.nnobs) then
      write(*, *)' Aborted because 2*nobs must be < nnobs'
      stop
      endif 


      write(*, *)' Training data will go from day: '
      read(*, *) day1
      write(*, *)' to day: '
      read(*, *) day2
      write(*, *)' Training data are between days', day1,' and', day2
      do iobs = 1, nobs
3     call randin(seed, 1, ndatos, j)
      if((datos(3,j).lt.day1).or.(datos(3,j).gt.day2))
     *   go to 3
      x(1, iobs) = datos(1, j)
      x(2, iobs) = datos(2, j)
      x(3, iobs) = datos(3, j)
      y(iobs) = datos(4, j)     
      end do

c      write(*, *)' Test data will go from day: '
c      read(*, *) day1
c      write(*, *)' to day: '
c      read(*, *) day2
      day1 = 10.9
      day2 = 11.1
c  Estas filas no sirven para nada:
      do iobs = nobs+1, 2*nobs
5     call randin(seed, 1, ndatos,j)
      if((datos(3,j).lt.day1).or.(datos(3,j).gt.day2)) go to 5 
      x(1, iobs) = datos(1, j)
      x(2, iobs) = datos(2, j)
      x(3, iobs) = datos(3, j)
      y(iobs) = datos(4, j)
      end do 
c  Hasta acá


      
c      do iobs = 1, 2*nobs
c      call randin(seed, 1, ndatos, j) 
c      x(1, iobs) = datos(1, j)
c      x(2, iobs) = datos(2, j)
c      x(3, iobs) = datos(3, j)
c      y(iobs) = datos(4, j)
c      end do

      if(ipri.eq.1) then
      write(*, *)
      write(*, *)' Data for Saint Venant (Qinlet, x, days, z):'
      write(*, *)' Observed data for training:'
      write(*, *)

      do iobs = 1, nobs
      write(*, *) iobs, (x(i, iobs),i=1,3), y(iobs)
      end do
      write(*, *)
      write(*, *)' Data that will be used for test:'
      do iobs = nobs+1, 2*nobs
      write(*, *) iobs, (x(i, iobs),i=1,3), y(iobs)
      end do    
      endif


      





      write(*, *)' Observation data have been read'
      write(*, *)' If want to use x = (Qinlet, position, time) type 1'
      write(*, *)' If want to use x = (Qinlet, position) type 2'
      read(*, *) i
      if(i.eq.2) then
      nvar = 2
      write(*, *)' If you want to include column Q^2 type 1'
      read(*, *) jqua
      if(jqua.eq.1) then
      nvar = nvar + 1
      do iobs = 1, 2*nobs
      x(nvar, iobs) = x(1, iobs)**2
      end do 
      endif
      write(*, *)' If you also want to include column x^2 type 2'
      read(*, *) jqua
      if(jqua.eq.2) then
      nvar = nvar+1
      do iobs = 1, 2*nobs
      x(nvar, iobs) = x(2, iobs)**2
      end do
      endif
      write(*, *)' If you also want to include column x*Q type 3'
      read(*, *) jqua
      if(jqua.eq.3) then
      nvar = nvar+1
      do iobs = 1, 2*nobs
      x(nvar, iobs) = x(1, iobs)*x(2, iobs)
      end do
      endif
      write(*, *)' If you also want to include column x^2*Q type 4'
      read(*, *) jqua
      if(jqua.eq.4) then
      nvar = nvar+1
      do iobs = 1, 2*nobs
      x(nvar, iobs) = x(1, iobs)**2*x(2, iobs)
      end do
      endif
      write(*, *)' If you also want to include column x*Q^2 type 5'
      read(*, *) jqua
      if(jqua.eq.5) then
      nvar = nvar+1
      do iobs = 1, 2*nobs
      x(nvar, iobs) = x(1, iobs)*x(2, iobs)**2
      end do
      endif      

      write(*, *)' If you also want to include column x^2*Q^2 type 6'
      read(*, *) jqua
      if(jqua.eq.6) then
      nvar = nvar+1
      do iobs = 1, 2*nobs
      x(nvar, iobs) = x(1, iobs)**2*x(2, iobs)**2
      end do
      endif       

      endif
c  This endif corresponds to if(i.eq.2) 


      nx = nvar 

      write(*, *)' Reading Saint Venant data is complete'
      endif          


c      write(*, *)' If you want to add quadratic dependence type 1'
c      write(*, *)' Otherwise type 2'
c      read(*, *) i
c      if(i.ne.1) then
c      nx = nvar
c      else
c      write(*, *)' Variables x include quadratic terms.'
c      do iobs=1, 2*nobs
c      x(nvar+1, iobs) = x(1, iobs)**2
c      x(nvar+2, iobs) = x(2, iobs)**2
c      x(nvar+3, iobs) = x(3, iobs)**3
c      x(nvar+4, iobs) = x(1, iobs)**2
c      x(nvar+5, iobs) = x(2, iobs)**2
c      x(nvar+6, iobs) = x(3, iobs)**2 
c      x(nvar+7, iobs) = x(1, iobs)**2 * x(2, iobs)
c      x(nvar+8, iobs) = x(1, iobs)**2 * x(3, iobs)
c      x(nvar+9, iobs) = x(2, iobs)**2 * x(3, iobs)
c      x(nvar+10, iobs) = x(1, iobs) * x(2, iobs)**2
c      x(nvar+11, iobs) = x(1, iobs) * x(3, iobs)**2
c      x(nvar+12, iobs) = x(2, iobs) * x(3, iobs)**2      
c      x(nvar+13, iobs) = dsin(x(1, iobs))
c      x(nvar+14, iobs) = dsin(x(2, iobs))
c      x(nvar+15, iobs) = dcos(x(3, iobs))
c      x(nvar+16, iobs) = dcos(x(1, iobs))
c      x(nvar+17, iobs) = dcos(x(2, iobs))
c      x(nvar+18, iobs) = dcos(x(3, iobs))
c      end do 
c      write(*, *)' How many additional input xs you want?'
c      write(*, *)' The first one is Q-inlet squared'
c      read(*, *) iadd
c      nx = nvar+iadd
c      write(*, *)' nx including quadratic dependences:', nx
c      endif 
             

c  End of if(kprob.eq.1.or.kprob.ne.2) .......************************************************


     
      if(kprob.ne.1.and.kprob.ne.2) then
      write(*, *)' Dimension of x :'
      write(*,*)' Default = 10'
      nx=10 
      if(.not.lovo) then
      write(*, *)' Number of observations nobs:'
      read(*, *) nobs
      endif
      if(nx.gt.nnx.or.nobs.gt.nnobs) then
      write(*, *)' Abort. Maximal nx and nobs are:', nnx, nnobs
      stop
      endif

      write(*, *)' nx = ', nx,' nobs = ', nobs

      write(*, *)' Data generation'
      write(*, *)' Type 1 for all-random generation'
      write(*, *)' Type 2 for function ||x||^2'
      write(*, *)' Type 3 for function sum sin(x_i)'
      write(*, *)' Type 4 for characteristic function of y > ||x||^2'
      write(*,*)' Default: 2 '
      idata=2



      if(idata.eq.1) then
c  All random generation
c  We generate 2 nobs observations in order to use
c  the last nobs as Test Set

      write(*, *)' All-random data generation'      
      do iobs = 1, nobs+nobs
      do i = 1, nx
      call rondo(seed, x(i,iobs))
      end do
      call rondo(seed, y(iobs))
      end do
      endif

      if(idata.eq.4) then
      do iobs = 1, nobs+nobs
      z = 0.d0
      do i = 1, nx-1
      call rondo(seed, x(i, iobs))
      z = z + x(i, iobs)**2
      end do
      call rondo(seed, x(nx, iobs))
      if(x(nx,iobs).ge.z) then
      y(iobs) = 1.d0
      else
      y(iobs) = 0.d0
      endif
      end do
      endif


      if(idata.eq.2.or.idata.eq.3) then
      do iobs = 1, nobs+nobs
       do i = 1, nx
       call rondo(seed, x(i,iobs))
       end do
       y(iobs) = 0.d0
       do i = 1, nx
       if(idata.eq.2) then
       y(iobs) = y(iobs) + x(i, iobs)**2
       endif
       if(idata.eq.3) then
       y(iobs) = y(iobs) + dsin(x(i, iobs))
       endif
       end do
       end do
c  Introducing outliers in the training set
      do i = 1, nout
      call rando(seed, z)
      iout = nobs*z + 1
      if(iout.gt.nobs) iout=nobs
      call rondo(seed,z)
      z = z*perout/100.d0
      write(*,*)' y(',iout,')',' is an outlier'
      write(*,*)' Original y=', y(iout)
      y(iout) = y(iout)+z*y(iout)
      write(*,*)' Outlier y =', y(iout)
      end do



      endif
      endif
c  This endif corresponds to "if(kprob.ne.1.and.kprob.ne.2)"


      write(*, *)
      write(*, *) ' Depth of the Neural Network desired = '
      write(*,*)' Default = 10'
      depth=10
      if(depth.gt.ndepth) then
      write(*, *)' Abort. Maximal allowed depth is:', ndepth
      stop
      endif
      write(*, *)
      write(*, *)' Depth = ', depth
      write(*, *)
      write(*,*)' The number of components of level', depth,' is set'
      write(*,*)' equal to 1 by this code'
      write(*,*)
      write(*,*)' Type the number of components of the first',depth-1,
     * ' levels'
      write(*,*)' Default = (9,...,9)'
      write(*,*)' If you want default type 0'
      nada=0
      if(nada.eq.0) then
      do i = 1, depth-1
      ncomps(i) = 9
      end do      
      else
      write(*,*) ' Components of each level 1,...,depth-1:'
      read(*,*)(ncomps(i),i=1,depth-1)
      endif
      ncomps(depth) = 1

      do i = 1, depth
      if(nx.gt.nrowsa.or.(ncomps(i).gt.nrowsa)) then
      write(*, *)' Abort, Maximal number of components at each level:',
     *  nrowsa
      stop
      endif
      end do
      write(*, *)
      do i = 1,depth
      write(*, *)' Components of level', i,'  (Rows of A number', i,'):'
     *  , ncomps(i)   
      end do
      

      nAb = (nx+1)*ncomps(1)
      do i = 2, depth
      nAb = nAb + (ncomps(i-1)+1)*ncomps(i)
      end do



      write(*, *)' Total number of coefficients (A, b) :', nAb
      if(nAb.gt.nnAb) then
      write(*, *)' Abort. Total number of coeffs (A, b) exceeds', nnAb
      stop
      endif

c      write(*, *)' To continue type 1'
c      read(*, *) i
 

      if(ipri.eq.1)then
      write(*, *)
      write(*, *) ' Data:'
      write(*, *)' iobs   x_1,...x_{nx}   y'
      do iobs = 1, nobs
      write(*, *) iobs, (x(i, iobs),i=1,nx), y(iobs)
      end do
      endif

      write(*, *)

6     write(*, *)' To continue type 1'
      i=1

      write(*, *)' Generation of Initial NN for sysredneu'
      write(*, *)' If you want a random generation type 1'
      write(*, *)
      write(*, *)' If you want to use the NN that comes '
      write(*, *)' from previous run, type 2'
      write(*, *)' Default: 1'

      inewran=1

7     if(inewran.eq.2) then
      write(*, *) ' Initial NN comes from previous run'
 
      else
      write(*, *)' Initial NN is random'
      endif 
 

c Generation of initial coefficients (A, b)
      if(inewran.eq.1) then
      ncoan = nx

      do k = 1, depth
      do i = 1, ncomps(k)
      do j = 1, ncoan 
      call rondo(seed, a(i, j, k))
      end do
      call rondo(seed, b(i, k))
      b(i, k) = b(i, k)+1.d0
      end do
      ncoan = ncomps(k)
      end do
      endif
c  End of if(inewran.eq.1)*************************************************

      if(inewran.eq.2) then
      write(*, *) ' The new run adds a coordinate to x'
      nx = nx+1
      nAb = (nx+1)*ncomps(1)
      do i = 2, depth
      nAb = nAb + (ncomps(i-1)+1)*ncomps(i)
      end do



      write(*, *)' Total number of coefficients (A, b) :', nAb
      if(nAb.gt.nnAb) then
      write(*, *)' Abort. Total number of coeffs (A, b) exceeds', nnAb
      stop
      endif
 



      do i = 1, ncomps(1)
      a(i, nx, 1) = 0.d0
      end do
      write(*, *)' Decide which should be the new coordinate of x'
      write(*, *)' Option 1:', ' Previous x(1) squared'
      write(*, *)' Option 2:', ' Previous x(2) squared'
      write(*, *)' Option 3:', ' Product x(1)*x(2)'
      write(*, *)' Option 4:', ' |abs(x(i)|^alfa * |abs(x(j)|^beta'
      write(*, *)' Option 5:',
     *   ' |abs(x(i)|^alfa * |abs(x(j)|^beta * |abs(x(k)|^gama'        

      write(*, *)
      write(*, *)' Option 6:', ' x(1)^2 * x(2)'
      write(*, *)' Option 7:', ' x(1) * x(2)^2'
      write(*, *)' Option 8:', ' x(1)^2 * x(2)^2'
   
      read(*, *) option

      if(option.eq.6) then 
      do iobs=1,2*nobs
      x(nx, iobs) = x(1, iobs)**2*x(2, iobs)
      end do
      endif

      if(option.eq.7) then 
      do iobs=1,2*nobs
      x(nx, iobs) = x(1, iobs)*x(2, iobs)**2
      end do
      endif

      if(option.eq.8) then 
      do iobs=1,2*nobs
      x(nx, iobs) = x(1, iobs)**2*x(2, iobs)**2
      end do
      endif
 


      if(option.eq.5) then

      write(*, *)' Type i, j, k, alfa, beta and gama'
      read(*, *) ia, ja, ka, alfa, beta, gama 

      if(ia.lt.1.or.ia.gt.nx-1.or.ja.lt.1.or.ja.gt.nx-1.or.
     * ka.lt.1.or.ka.gt.nx-1) then
      write(*, *)' Abort.Indices i,j,k given for Option 5 out of range'
      stop
      endif
           
      do iobs = 1, 2*nobs
      if(x(ia,iobs).eq.0.d0.or.x(ja,iobs).eq.0.d0.or.x(ka,iobs).eq.0.d0)
     *  then
      x(nx, iobs) = 0.d0
      else
      x(nx, iobs) = dabs(x(ia,iobs))**alfa * dabs(x(ja,iobs))**beta
     *  *dabs(x(ka,iobs))**gama
      endif
      end do
      endif
 


      if(option.eq.4) then

      write(*, *)' Type i, j, alfa and beta'
      read(*, *) ia, ja, alfa, beta

      if(ia.lt.1.or.ia.gt.nx-1.or.ja.lt.1.or.ja.gt.nx-1) then
      write(*, *)' Abort. Indices i, j given for Option 4 out of range'
      stop
      endif
           
      do iobs = 1, 2*nobs
      if(x(ia, iobs).eq.0.d0.or.x(ja, iobs).eq.0.d0) then
      x(nx, iobs) = 0.d0
      else
      x(nx, iobs) = dabs(x(ia,iobs))**alfa * dabs(x(ja,iobs))**beta
      endif
      end do
      endif


      if(option.eq.1) then
      do iobs = 1, 2*nobs
      x(nx, iobs) = x(1, iobs)**2
      end do
      endif
      if(option.eq.2) then
      do iobs = 1, 2*nobs
      x(nx, iobs) = x(2, iobs)**2
      end do
      endif
      if(option.eq.3) then
      do iobs = 1, 2*nobs
      x(nx, iobs) = x(1, iobs) * x(2, iobs)
      end do
      endif       
      endif

c endif corresponding to if(inewran.eq.2)*******************************************
     
      if(ipri.eq.1) then
      ncoan = nx   
      do k = 1, depth
      
      write(*, *)' Matrix A number', k
      do i = 1, ncomps(k)
      write(*, *)(a(i, j, k),j=1,ncoan)
      end do
      ncoan = ncomps(k)
      end do
      
      do k = 1, depth
      write(*, *)' Vector b number', k
      write(*, *)(b(i, k),i=1,ncomps(k))
      end do
      endif
c end of printing***********************************************************************

      k = 1
      ncoan = nx
      do idep = 1, depth
      do i = 1, ncomps(idep)
      do j = 1, ncoan
      Ab(k) = a(i, j, idep) 
      k = k + 1
      end do
      Ab(k) = b(i, idep) 
      k = k + 1
      end do
      ncoan = ncomps(idep)
      end do
      if(ipri.eq.1) then
      write(*, *)' Vector Ab:'
      write(*, *)(Ab(i),i=1,nAb)
      endif

      maxrow = min0(maxrow, nobs)
      call cpu_time(start)
2     call sysredneu(x,nx,y,nobs,depth,ncomps,nAb,Ab,npsi,f,jacob,
     *jacosave,fsave,nnobs,nnx,nrowsa,jaca,jacx,jacaux,a,b,kon,nback,
     *max,seed,
     *under,eps,delta,maxrow,nref,iqn,ipri,aat,aatsave,vaux,iperm,yaux1,
     *yaux2,zaux,mgsyes,esprin,sizeb,rmsd,ceroum,lovo,plovo,mataux,
     *Abtrial,ftrial,elimi)
      call cpu_time(finish)

      write(*,*)' CPU-time=', finish-start
      write(*, *)' Number of iterations performed:', kon
      write(*, *)' Number of backtrackings:', nback
      zmedio=0.d0
      do j = 1, nobs
      zmedio = zmedio + y(j)
      end do
      zmedio = zmedio/dfloat(nobs)
      write(*, *)' Average observed y"s:', zmedio, ' 100 rmsd/ave:',
     *   100*rmsd/zmedio         


      if(nAb.le.10) then
      write(*, *)' Solution Ab found:'
      write(*, *)(Ab(i),i=1,nAb)
      endif
 
 
      if(max.eq.0) then 
      write(*, *)' Average observed y"s:', zmedio, ' 100 rmsd/ave:',
     *   100*rmsd/zmedio
      endif

      write(*, *)' Recall that:'

      write(*, *)' nx = ', nx,' nobs = ', nobs
      write(*, *)' depth = ', depth
      write(*, *)' Number of function psi:', npsi

      write(*, *)' Total number of coefficients (A, b) :', nAb         

      do i = 1,depth
      write(*, *)' Components of level', i,'  (Rows of A number', i,'):'
     *  , ncomps(i)   
      end do
     
c  Transform the vector Ab into the matrices a(i,j,idep) and b(i,idep) 
      k = 1
      ncoan = nx
      do idep = 1, depth
      do i = 1, ncomps(idep)
      do j = 1, ncoan
      a(i, j, idep) = Ab(k) 
      k = k + 1
      end do
      b(i, idep) = Ab(k) 
      k = k + 1
      end do
      ncoan = ncomps(idep)
      end do

      if(max.ne.0) then
      write(*, *)
      write(*, *)' If you want to expand the vector x using the'
      write(*, *)' the actual NN as initial point, type 2'
      write(*, *)' If you want to finish the optimization procedure'
      write(*, *)' going to test set verifications, type 1'  
      write(*,*)' If you want to stop this execution type 0'
      read(*, *) inewran
      write(*, *)
      if(inewran.eq.0) stop
      if(inewran.eq.2) then
      write(*, *)' Provide new values for epsilon and delta:'
      write(*, *)' Recall previous ones: ', eps, delta
      read(*, *) eps, delta
      write(*, *)' We will use the NN output obtained to initialize'
      write(*, *)' new optimization process'
      go to 7
      endif
      endif
 





      write(*, *)' To continue with Test Set, type 1'
      read(*, *) i
      if(i.ne.1) stop

      maxrow = nobs


      max = 0


      write(*, *)' Minimal and maximal time (in days) for test set:'
      read(*, *) day1, day2

      nobs = 2000
      maxrow = nobs
      do iobs = 1,nobs 
8     call randin(seed, 1, ndatos,j)
      if((datos(3,j).lt.day1).or.(datos(3,j).gt.day2)) go to 8  
      x(1, iobs) = datos(1, j)
      x(2, iobs) = datos(2, j)
      x(3, iobs) = datos(3, j)
      if(jqua.ge.1) x(3, iobs) = x(1, iobs)**2
      if(jqua.ge.2) x(4, iobs) = x(2, iobs)**2
      if(jqua.ge.3) x(5, iobs) = x(1, iobs)*x(2, iobs)
      if(jqua.ge.4) x(6, iobs) = x(1, iobs)**2*x(2, iobs)    
      if(jqua.ge.5) x(7, iobs) = x(1, iobs)*x(2, iobs)**2   
      if(jqua.ge.6) x(8, iobs) = x(1, iobs)**2*x(2, iobs)**2   
      y(iobs) = datos(4, j)
      end do 


      go to 2          


      stop
      end


 


      subroutine sysredneu(x,nx,y,nobs,depth,ncomps,nAb,Ab,npsi,f,jacob,
     *jacosave,fsave,nnobs,nnx,nrowsa,jaca,jacx,jacaux,a,b,kon,nback,
     *max,seed,
     *under,eps,delta,maxrow,nref,iqn,ipri,aat,aatsave,aux,iperm,yaux1,
     *yaux2,zaux,mgsyes,esprin,sizeb,rmsd,ceroum,lovo,plovo,mataux,
     *Abtrial,ftrial,elimi)

      implicit none
      integer nnobs, nrowsa, nnx, sizeb 
      integer nx, npsi, nref, kon, nback
      integer nobs, maxrow, nAb
c  positions for aux must be at least 2*nobs      
      double precision x(nnx, nobs), y(nobs), aux(2*nobs)
      double precision jaca(nrowsa, *),jacx(nrowsa, *),jacaux(nrowsa,*)
      double precision aat(min0(nAb,nobs),min0(nAb,nobs))
      double precision  aatsave(min0(nAb,nobs),min0(nAb,nobs))
      double precision yaux1(nrowsa), yaux2(nrowsa), zaux(nrowsa,nrowsa)
      double precision mataux(nobs*nobs)
      integer iperm(nobs)
      logical lovo, elimi
      logical ceroum 
      integer correc
      integer plovo
      double precision predicted
      logical posdef
      
      double precision ff(1)

      integer ndats, one, mgsyes
      double precision Ab(nAb), Abtrial(nAb), f(nobs), jacob(maxrow,nAb)
      double precision ftrial(nobs)
      double precision jacosave(maxrow,nAb), fsave(nobs)
      integer k, idep, depth, ncoan
      integer ncomps(depth)
      integer i, j
      parameter(one=1)
      double precision a(nrowsa, nrowsa, depth), b(nrowsa, depth)
      integer iobs
      integer modo
      double precision under
      integer max, mlin, iermgs
      integer sorte
      integer iqn, ipri, esprin, nada
      integer imax, ithis, nabs
      double precision fmax, dife, gnor
 

      double precision s(nAb)
      double precision snor, z, seed, add

      double precision error, err2, rmsd,rmsdtrial,eps, delta
      double precision sts
      double precision diago
 
      integer discar

c The observations are:
c     x(1, 1), ..., x(nx, 1),    y(1)
c     x(1, 2), ..., x(nx, 2),    y(2)
c ................................
c     x(1, nobs)...,x(nx, nobs), y(nobs)
c
c (This table explains the meaning of the inputs x, y, nx, and nobs.
c
c  Attention: Notice that each component of y is a single number. 
c  This means that we are restricted to Neural Networks that are scalar functions
c  In our notation ncomps(depth) is always equal to 1. 


c nobs is the number of observations
c nAb is the number of parameters (A, b) of the Neural Network
c Ab is the vector of parameters (A, b) in the form 
c ( row corresponding to A_1, b_1, ..., row corresponding to A_depth b_depth)
c where each "row corresponding to A_k, b_k 
c is is displayed in the form (row 1 of A_k, entry 1 of b_k, ..., row ncomp(k) of A_k, ..., entry ncomp(k) of b_k)
c 
      kon = 0
      nback = 0

      if(lovo) then
      nabs = plovo
      maxrow = plovo
      else
      nabs = nobs
      endif
     


      if(maxrow.gt.nabs) maxrow=nabs


      if(ncomps(depth).ne.1) then
      write(*, *)' Error.In this implementation ncomps(depth) must be 1'
      stop
      endif

      kon = 0


c Evaluation of the system, obtaining F(Ab) and Jacobian(Ab)   
c Transformation of vector Ab into matrices A and vectors b

      write(*, *)' In this problem n_z = ', nx,'  nobs=', nobs
      write(*, *)' Number of coefficients to be estimated:', nAb

      if(nAb.lt.nobs) then
      write(*, *)' nAb = ', nAb,' nobs =', nobs
      write(*, *)' Warning! This method was firstly developed ',
     * ' for nAb .ge. nobs!'
      endif
 
      if(.not.elimi) then
      write(*, *)' To continue type 1'
      read(*, *) nada 
      endif

     


c  Transform the vector Ab into the matrices a and the vectors b

1     k = 1
      ncoan = nx
      do idep = 1, depth
      do i = 1, ncomps(idep)
      do j = 1, ncoan
      a(i, j, idep) = Ab(k)
      k = k + 1
      end do
      b(i, idep) = Ab(k)
      k = k + 1
      end do
      ncoan = ncomps(idep)
      end do
c***********************************************************************
      if(.not.elimi) then
      if(mod(kon,esprin).eq.0) then
      write(*, *)
      write(*, *)' Iteration:', kon

c      write(*, *)' Ab=', (Ab(i),i=1,nAb)

      write(*, *)' Ab(1)=', Ab(1),' ... Ab(',nAb,') =', Ab(nAb)
      write(*, *)
      endif
      endif

      if(.not.elimi) then
      write(*,*)' To continue type 1'
      read(*,*) nada
      endif



c Effective evaluation of the system , obtaining F(Ab) and Jacobian(Ab)

      if(mod(kon, nref+1).eq.0) then
      modo = 2
      else
      modo = 1
      endif

      if(maxrow.lt.nobs) then
      do iobs = 1, maxrow
      call randin(seed,iobs,nobs,sorte)
      i = iperm(iobs)
      iperm(iobs) = iperm(sorte)
      iperm(sorte) = i
      end do      
      endif
 


c  System and perhaps Jacobian evaluation

      if(lovo) then
      do iobs=1,nobs
      iperm(iobs)=iobs
      end do
      endif

      do iobs = 1, nobs 

      
      call redneu(a, b, x(1, iperm(iobs)), nx, one, npsi, depth, 
     *   ncomps,ff(1), jaca, jacx, jacaux, modo,nrowsa,yaux1,yaux2,zaux)

      if(.not.elimi) then
      if(iobs.eq.1) then
      write(*,*)' ff(1) after first call redneu to compute f(1)', ff(1)
      write(*,*)' To continue type 1'
      read(*,*) nada
      endif
      endif

c      if(ipri.eq.1) then
c      write(*, *)' Observation ', iperm(iobs),' Predicted:', ff(1),
c     *     ' Observed:', y(iperm(iobs))
c      endif

      f(iobs) = ff(1) - y(iperm(iobs))

      if(modo.eq.2) then
      do j = 1, nAb
      jacob(iobs, j) = jaca(ncomps(depth), j)
      end do   
      endif     

      if(f(iobs).lt.0.d0) then
      f(iobs) = - f(iobs)
      if(modo.eq.2) then
      do j = 1, nAb
      jacob(iobs,j) = - jacob(iobs, j)
      end do
      endif
      endif
     
      end do
c  This end do corresponds to "do iobs = 1, nobs"       

c********************************************************************************
      if(lovo) then

c      write(*,*)' Modula of Residuals before Lovization:'
c      do iobs = 1, nobs
c      write(*,*)' iobs = ', iobs,' f(iobs)=', f(iobs)
c      end do
      

c  Compute the p functions that will be considered at this iteration of LOVO
      do j = 1, nobs-plovo
      fmax = -1.d0      
        do iobs = 1, nobs
          if(f(iobs).ge.fmax)then
          imax = iobs
          fmax = dabs(f(iobs))
          endif
        end do
      f(imax) = -1.d99
      end do
c   This end do corresponds to "do j = 1, nobs-plovo"      
       

c      write(*,*)' Modulos of Residuos preprocessados:'
c      do iobs = 1, nobs
c      write(*,*)' iobs=', iobs,' f(iobs)=', f(iobs)
c      end do

c      do iobs=1,nobs
c      write(*,*)' Jacob before Lovization:'
c      write(*,*)(jacob(iobs,j),j=1,nAb)
c      end do

c      write(*,*)' To continue type 1'
c      read(*,*) nada
 

      ithis = 0
      do iobs = 1, nobs
      if(f(iobs).gt.-1.d99) then
      ithis = ithis+1
      f(ithis) = f(iobs)
      if(modo.eq.2) then
      do j = 1, nAb
      jacob(ithis, j) = jacob(iobs,j)
      end do
      endif
c  This endif corresponds to "if(modo.eq.2)"
      endif
c  This endif corresponds to "if(f(iobs).gt.-1.d99)"
      end do
c  This end do corresponds to "do iobs=1,nobs"

c      write(*,*)' Modula of Residuals after Lovization:'
c      do iobs = 1, plovo
c      write(*,*)' iobs = ', iobs,' f(iobs)=', f(iobs)
c      end do
c      write(*,*)
c      do iobs = plovo+1,nobs
c      write(*,*)' iobs =', iobs,' f(iobs)=', f(iobs)
c      end do

c      write(*,*)' Jacob after Lovization:'
c      write(*,*)' ithis=', ithis,' must be ', plovo

c      Compute the gradient, which is not necessary in the case that nAb >0 plovo
      do j = 1, nAb
      aux(j) = 0.d0
      do iobs = 1, plovo
      aux(j) = aux(j) + jacob(iobs,j)*f(iobs)
      end do
      end do
      gnor = 0.d0
      do j = 1, nAb
      gnor = dmax1(gnor, dabs(aux(j)))
      end do


 


      endif
c  This endif corresponds to "if(lovo) then" 
c*************************************************************************************************************
 
      if(kon.ge.1.and.lovo) then
      write(*,*)' Comparison of Ab versus previous Abtrial (=?):'
      dife = 0.d0
      do i = 1,nAb
      dife = dmax1(dife, dabs(Ab(i)-Abtrial(i)))
c      write(*,*) i, Ab(i), Abtrial(i)
      end do

      if(.not.elimi) then
      write(*,*)' Difference between Ab and previous Abtrial:', dife
      write(*,*)' Comparison of f versus previous ftrial (should be =)'

      write(*,*)' f(1)=', f(1),' ftrial(1)=', ftrial(1)
      write(*,*)' To continue type 1'
      read(*,*) nada
      endif

c      do iobs = 1, plovo
c      write(*,*) iobs, f(iobs), ftrial(iobs)
c      end do
      endif

       

c End of system and (perhaps)  Jacobian evaluation


c  Recall that ncomps(depth) = 1 in this implementation

      if(lovo) maxrow = plovo

c      do iobs = 1, maxrow 
c      if(ipri.eq.1) then
c      write(*, *)' f(', iperm(iobs),') = ', f(iobs),
c     *   ' for data =', y(iperm(iobs))
c      endif
c      end do


c********************************************************************************
c  Test stopping criterion


      if(maxrow.eq.nobs.or.lovo) then
      nabs = nobs
      if(lovo) then
      nabs = plovo
      maxrow = plovo
      endif
      error = 0.d0
      err2 = 0.d0
      do iobs = 1, plovo 
      error = dmax1(error, dabs(f(iobs)))
      err2 = err2 + f(iobs)**2
      end do

      if(ceroum) then
      correc = 0
      do iobs = 1, nobs
      predicted = f(iobs) + y(iperm(iobs))
      if(predicted.ge.0.5d0) then
      if(y(iperm(iobs)).ge.0.5d0) correc = correc+1 
      else
      if(y(iperm(iobs)).lt.0.5d0) correc = correc +1
      endif
      end do
      write(*, *)' Number of corrected predictions:', correc,
     * ' Percentage:', 100*dfloat(correc)/dfloat(nobs)

      endif
     

      rmsd = dsqrt(err2/dfloat(plovo))
      if(mod(kon,esprin).eq.0) then
      write(*,*)
      write(*,*)
      write(*,*)' Iteration ', kon
      write(*, *)' Error (max diff between predicted & observed):',error
      write(*, *)' Sum of squares:', err2,' RMSD =', rmsd
      endif
      write(*,*)' Sup norm of the gradient:', gnor
      endif
      
      if(max.eq.0) then
      write(*, *)' Above resuts have been obtained for conference test'
      write(*, *)' or for test set, as maximum of iterations was 0'
      return
      endif  

      if(mod(kon,esprin).eq.0) then
      write(*, *)' Number of iterations up to now:', kon             
      endif

      if(lovo.or.(maxrow.eq.nobs)) then
      if(rmsd.le.eps) then
      write(*, *)
      write(*, *)'************************************'
      write(*, *)' Return by small RMSD residual error, kon=',kon
      return
      endif
      endif

      if(kon.ge.max) then
      write(*, *)' Number of iterations up to now:', kon     
      write(*, *)' Return by maximum of iterations', max,' reached'
      return
      endif 

c***********************************************************************
      if(nabs.le.nAb) then

      write(*,*)' nabs =', nabs,' maxrow =', maxrow


c      write(*, *)' Parameters entry of mgs:'
c      write(*, *)' nobs=', nobs,' nAb=', nAb

      do iobs = 1, maxrow 
      do j = 1, nAb
      jacosave(iobs, j) = jacob(iobs, j)
      end do
      fsave(iobs) = f(iobs)
      end do

c      write(*, *)' Independent term that enters to mgs:',(f(i),i=1,maxrow)
c      write(*, *)' Matrix jacob that enters to mgs:'
      
c      do i = 1, nobs 
c      write(*, *)(jacob(i,j),j=1,nAb)
c      end do 



      
      z = 0.d0
      do j = 1, nAb
      do i = 1, maxrow
      z = dmax1(z, dabs(jacob(i,j)))
      end do
      end do
      if(mod(kon,esprin).eq.1) then
      write(*, *)' Maximal element of ', maxrow,' rows of Jacob:', z
      endif

c**************************** mgysyes=1 ***********************************************
      if(mgsyes.eq.1) then

c  If mgsyes=1 and iermgs=0 (full row-rank of Jacobian) the increment
c  is computed by mgs
c  Observar que el parametro nobs en call mgs no cumple ningun papel. Tocar por maxrow. Trocado.
      call mgs(jacob, f, s,maxrow,nAb,iermgs,discar)

      if(mod(kon,esprin).eq.0) then
      write(*, *)' Discarded rows at mgs:', discar,' Rank of Jacobian:', 
     *   maxrow-discar
      write(*, *)' Error code ier of mgs:', iermgs
      if(iermgs.eq.0) write(*, *)' Increment computed by iermgs'
      endif
      endif 
c****************************************************************************************

      if(iermgs.ne.0.or.mgsyes.ne.1) then
c  Use Cholesky, either by decision of not using mgs or by rank-defficiency

      if(nref.eq.0.or.modo.ne.1.or.iqn.eq.1) then

      do iobs = 1, nabs
      do j = 1, nAb
      jacob(iobs, j) = jacosave(iobs, j)
      end do
      f(iobs) = fsave(iobs)
      end do

c      do j = 1, maxrow
c      do i = 1, maxrow
c      aat(i,j)=0.d0
c        do k = 1, nAb
c        aat(i,j)=aat(i,j) + jacosave(i,k)*jacosave(j,k)
c        end do
c        aatsave(i,j) = aat(i,j)
c      end do
c      end do

      do j = 1, maxrow
      do i = 1, maxrow
      aat(i,j)=0.d0
      end do
      end do

        do k = 1, nAb
        do j = 1, maxrow
        do i = 1, maxrow
        aat(i,j)=aat(i,j) + jacosave(i,k)*jacosave(j,k)
        end do
        end do
        end do

      do j = 1, maxrow
      do i = 1, maxrow
      aatsave(i,j) = aat(i,j)
      end do
      end do

c      write(*,*)' Matrix A A^t:'
c      do i = 1, maxrow
c      write(*,*)(aat(i,j),j=1,maxrow)
c      end do

c      write(*,*)' To continue type 1'
c      read(*,*) nada


       diago = 0.d0
       do i = 1, maxrow
       diago = dmax1(diago, dabs(aat(i,i)))
       end do
      if(mod(kon,esprin).eq.0) then
      write(*, *)' Maximal diagonal of aat:', diago           
      endif

      if(.not.elimi) then
      write(*,*)' To continue type 1'
      read(*,*) nada     
      endif

c  Decide the first augmentation of the diagonal of aat
c  If mgsyes=0 the augmentation should be 0, because
c  perhaps aat is nonsingular (so, posdef)
      if(mgsyes.eq.0) then
      z = 0.d0
      else
c  When mgsyes=1 this means that the Jacobian is rank-defficient
c  So, aat is singular, so we begin adding a positive element
c  in the diagonal 



      if(z.eq.0.d0) then
c      if(mod(kon,esprin).eq.0) then
c      write(*, *)' Diagonal of AA^t is null'  
c      write(*,*)' To continue type 1'
c      read(*,*) nada
c      endif
      z = 0.1d0
      else
      z = 0.1d0
      endif
      endif
c  Here finished the decision on the first diagonal augmentation z

3     do j  = 1, maxrow
      do i = 1, maxrow
      aat(i,j)=aatsave(i,j) 
      end do 
      aat(j, j) = aatsave(j,j)+z
      end do

c      write(*,*)' Matrix A A^t after possible correction z=', z
c      do i = 1, maxrow
c      write(*,*)(aat(i,j),j=1,maxrow)
c      end do
c      write(*,*)' To continue type 1'
c      read(*,*) nada

      do j = 1, maxrow
      do i = 1, maxrow
      mataux((j-1)*maxrow+i) = aat(i,j)
      end do
      end do


      call chole (maxrow, mataux, posdef, add)    
      if(posdef) then
      if(mod(kon,esprin).eq.0) then
      write(*, *)' For augmentation z=',z,' Matrix aat is posdef'
      endif
      else 
      if(mod(kon,esprin).eq.0) then
      write(*, *)' For augmentation z=',z,' Matrix aat isnt posdef yet'
      write(*,*)' To continue type 1'
      read(*,*) nada
      endif
c      z = dmax1(1.d-6*diago, 100.d0*z)
      z = dmax1(0.1d0, 10.d0*z) 
c      if(z.eq.0.d0) z=1.d0
      goto 3 
      endif

      endif
c  This endif corresponds to "  if(nref.eq.0.or.modo.ne.1.or.iqn.eq..1)   "   

      call sicho (maxrow, mataux, aux, f, aux(nabs+1)) 

      do j = 1, nAb
      s(j) = 0.d0
      do iobs=1,maxrow
      s(j) = s(j)+jacob(iobs,j)*aux(iobs)
      end do
      end do  

      endif
c  This endif corresponds to "if(iermgs.ne.0.or.mgsyes.ne.1)"
c************************************************************************************

c  Up to now the basic increment s(nAb components) was computed
c  perhaps by mgs, perhaps by Cholesky

      do iobs = 1, nabs
      do j = 1, nAb
      jacob(iobs, j) = jacosave(iobs, j)
      end do
      f(iobs) = fsave(iobs)
      end do          

      endif
c  This endif corresponds to "if(nobs.le.nAb)"
c  


c*************************************************************************

      if(nobs.gt.nAb) then
      do iobs = 1, maxrow 
      do j = 1, nAb
      jacosave(iobs, j) = jacob(iobs, j)
      end do
      fsave(iobs) = f(iobs)
      end do

      z = 0.d0
      do j = 1, nAb
      do i = 1, maxrow
      z = dmax1(z, dabs(jacob(i,j)))
      end do
      end do      

c  Use Cholesky, either by decision of not using mgs or by rank-defficiency

      if(nref.eq.0.or.modo.ne.1.or.iqn.eq.1) then

      do iobs = 1, maxrow
      do j = 1, nAb
      jacob(iobs, j) = jacosave(iobs, j)
      end do
      f(iobs) = fsave(iobs)
      end do

      do j = 1, nAb
      do i = 1, nAb
      aat(i,j)=0.d0
      do iobs = 1, maxrow
      aat(i,j) = aat(i,j) + jacob(iobs, i)*jacob(iobs,j)
      end do
      end do
      end do

      do j = 1, nAb
      do i = 1, nAb
      aatsave(i,j) = aat(i,j)
      end do
      end do

       diago = 0.d0
       do i = 1, nAb
       diago = dmax1(diago, dabs(aat(i,i)))
       end do
      if(mod(kon,esprin).eq.0) then
      write(*, *)' Maximal diagonal of A^t A :', diago           
      endif

      if(diago.eq.0.d0) then
      diago = 1.d0
      do i = 1, nAb
      aat(i,i)=1.d0
      end do
      endif


      z = 0.d0

4     do j  = 1, nAb 
      do i = 1, nAb
      aat(i,j)=aatsave(i,j) 
      end do 
      aat(j, j) = aatsave(j,j)+z
      end do

      call chole (nAb, aat, posdef, add)    
      if(posdef) then
          if(mod(kon,esprin).eq.0) then
          write(*, *)' For augmentation z=',z,' Matrix A^t A  is posdef'
          endif
      else 
      if(mod(kon,esprin).eq.0) then
      write(*, *)' For augmentation z=',z,' Matrix A^tA isnt posdef yet'
      endif
      z = dmax1(1.d-6*diago, 100.d0*z) 
      goto 4                  
      endif

c  Compute A^t F
      do j = 1, nAb
      aux(j) = 0.d0
      do i = 1, maxrow
      aux(j) = aux(j)+jacob(i,j)*f(i)
      end do
      end do

      call sicho (nAb, aat, s, aux, aux(nobs+1)) 





      endif
c  This endif corresponds to "  if(nref.eq.0.or.modo.ne.1.or.iqn.eq..1)   "   
      




      endif

c  This endif corresponds to if(nAb.lt.nobs)

c*************************************************************************


      snor = 0.d0
      do i = 1, nAb
      snor = dmax1(snor, dabs(s(i)))
      end do
      if(mod(kon,esprin).eq.0) then
      write(*, *)' Sup norm of the increment s:', snor
      endif
      if(snor.le.delta*dfloat(maxrow)/dfloat(nabs)) then
      write(*, *)' Finish because small increment'
      write(*, *)' Iterations performed:', kon
      return
      endif

c Compute Abtrial  

5     do i = 1, nAb
      s(i) = - under * s(i)
      Abtrial(i) = Ab(i) +  s(i)
      end do

      write(*,*)' The increment has been computed but we want to know'
      write(*,*)' whether the objective function decreased'
      if(.not.elimi) then
      write(*,*)' To continue type 1'
      read(*,*) nada
      endif
 
      if(lovo) then
      k = 1
      ncoan = nx
      do idep = 1, depth
      do i = 1, ncomps(idep)
      do j = 1, ncoan
      a(i, j, idep) = Abtrial(k)
      k = k + 1
      end do
      b(i, idep) = Abtrial(k)
      k = k + 1
      end do
      ncoan = ncomps(idep)
      end do      
      maxrow = nobs
      do iobs=1,nobs
      iperm(iobs)=iobs
      end do
      endif

      do iobs = 1, nobs  

      modo = 1 
      call redneu(a, b, x(1, iperm(iobs)), nx, one, npsi, depth, 
     *   ncomps,ff(1), jaca, jacx, jacaux, modo,nrowsa,yaux1,yaux2,zaux)

      

      if(ipri.eq.1) then
      write(*, *)' Observation ', iperm(iobs),' Predicted:', ff(1),
     *     ' Observed:', y(iperm(iobs))
      endif

      ftrial(iobs) = ff(1) - y(iperm(iobs))

c      if(iobs.eq.1) then
c      write(*,*)' ff(1) after call redneu:', ff(1)
c      write(*,*)' ftrial(1) after call redneu:', ftrial(1)
c      write(*,*)' To continue type 1'
c      read(*,*) nada
c      endif

      modo=1
      if(modo.eq.2) then
      do j = 1, nAb
      jacob(iobs, j) = jaca(ncomps(depth), j)
      end do   
      endif     

      if(ftrial(iobs).lt.0.d0) then
      ftrial(iobs) = - ftrial(iobs)
      if(modo.eq.2) then
      do j = 1, nAb
      jacob(iobs,j) = - jacob(iobs, j)
      end do
      endif
      endif
     
      end do
c  This end do corresponds to "do iobs = 1, maxrow"       

c********************************************************************************
      if(lovo) then

c      write(*,*)' Modula of Residuals before Lovization:'
c      do iobs = 1, nobs
c      write(*,*)' iobs = ', iobs,' ftrial(iobs)=', ftrial(iobs)
c      end do
      

c  Compute the p functions that will be considered at this iteration of LOVO
      do j = 1, nobs-plovo
      fmax = -1.d0      
        do iobs = 1, nobs
          if(ftrial(iobs).ge.fmax)then
          imax = iobs
          fmax = dabs(ftrial(iobs))
          endif
        end do
      ftrial(imax) = -1.d99
      end do
c   This end do corresponds to "do j = 1, nobs-plovo"      
       

c      write(*,*)' Modulos of Residuos preprocessados:'
c      do iobs = 1, nobs
c      write(*,*)' iobs=', iobs,' ftrial(iobs)=', ftrial(iobs)
c      end do


      ithis = 0
      do iobs = 1, nobs
      if(ftrial(iobs).gt.-1.d99) then
      ithis = ithis+1
      ftrial(ithis) = ftrial(iobs)
      if(modo.eq.2) then
      do j = 1, nAb
      jacob(ithis, j) = jacob(iobs,j)
      end do
      endif
c  This endif corresponds to "if(modo.eq.2)"
      endif
c  This endif corresponds to "if(f(iobs).gt.-1.d99)"
      end do
c  This end do corresponds to "do iobs=1,nobs"

c      write(*,*)' Modula of Residuals after Lovization:'
c      do iobs = 1, plovo
c      write(*,*)' iobs = ', iobs,' ftrial(iobs)=', ftrial(iobs)
c      end do
c      write(*,*)



      endif
c  This endif corresponds to "if(lovo) then" 
c*************************************************************************************************************
      if(maxrow.eq.nobs.or.lovo) then
      nabs = nobs
      if(lovo) then
      nabs = plovo
      maxrow = plovo
      endif
      error = 0.d0
      err2 = 0.d0
      do iobs = 1, plovo 
      error = dmax1(error, dabs(ftrial(iobs)))
      err2 = err2 + ftrial(iobs)**2
      end do

      rmsdtrial = dsqrt(err2/dfloat(plovo))
      endif     

      write(*,*)' rmsd trial =', rmsdtrial,' rmsd current=', rmsd

c      write(*,*)' To continue type 1'     
c      read(*,*) nada



      if(rmsdtrial.le.rmsd) then
      modo = 2
      kon = kon + 1
      do i = 1, nAb
      Ab(i) = Abtrial(i)
      end do
      write(*,*)' ftrial(1) before go to 1:', ftrial(1)
      go to 1
      endif

c  We increase the regularization parameter
      nback = nback+1
c     z = dmax1(1.d0, 2.d0*z)
      z = dmax1(0.1d0, 10.d0*z)
      write(*,*)' rmsd trial was bigger than rmsd current'
      write(*,*)' New regularization:', z

      do j  = 1, maxrow
      do i = 1, maxrow
      aat(i,j)=aatsave(i,j) 
      end do 
      aat(j, j) = aatsave(j,j)+z
      end do

c      write(*,*)' Matrix A A^t after possible correction z=', z
c      do i = 1, maxrow
c      write(*,*)(aat(i,j),j=1,maxrow)
c      end do
c      write(*,*)' To continue type 1'
c      read(*,*) nada

      do j = 1, maxrow
      do i = 1, maxrow
      mataux((j-1)*maxrow+i) = aat(i,j)
      end do
      end do


      call chole (maxrow, mataux, posdef, add)    
      if(.not.posdef) then
      write(*, *)' For augmentation z=',z,' Matrix aat isnt posdef yet'
      write(*,*)' This is not possible. Something is wrong'
      stop
      endif

c  This endif corresponds to "  if(nref.eq.0.or.modo.ne.1.or.iqn.eq..1)   "   

      call sicho (maxrow, mataux, aux, f, aux(nabs+1)) 

      do j = 1, nAb
      s(j) = 0.d0
      do iobs=1,maxrow
      s(j) = s(j)+jacob(iobs,j)*aux(iobs)
      end do
      end do 
      go to 5      

   




      end
      

   	subroutine mgs (a, b, x, maxrow, n, ier, discar)
	implicit none 
      integer m, n, ier  
      double precision a(maxrow, n), b(maxrow),  x(n)
      integer ii, i, j
      double precision anor, z, tol, tolanor
      double precision pesca
      double precision seed
      integer discar, maxrow
      double precision snor

c  Modified  Gram-Schmidt subroutine for solving A x = b
c  with A being m x n and m <= n.
c  The matrix and the independent term are  destroyed.
c  Computes minimum norm solution when the system is 
c  compatible. 
c  If the system is not compatible
c  the equations whose rows are dependent with respect
c  to other ones are eliminated.

      discar = 0

	ier = 0
	if(maxrow.eq.0.or.n.eq.0) then
	ier = 1
	return
	endif

c  Compute the maximal norm of the rows of A
	
	anor = 0.d0
	do j=1,n
	z = 0.d0
	do i=1,maxrow  
	z = z + a(i,j)*a(i,j)
	end do
	z = dsqrt(z)
	anor = dmax1(anor, z)
	end do

	if(anor.eq.0.d0) then
      write(*, *)' Matrix entering mgs is null'
      discar = maxrow
	ier = 2
      do j = 1, n
      x(j) = 0.d0
      end do
	return
	endif

	tol = 1.d-12

	tolanor = tol*anor

	do i = 1, maxrow
c  Normalize row i
	z = 0.d0
	do j = 1, n
	z = z + a(i, j)*a(i, j)
	end do
      z = dsqrt(z)

	if(z.lt.tolanor) then
c  Row i of A is a linear combination of rows 1,...,i-1
c  We replace this row by (0, ..., 0) and we set, 
c  accordingly, b(i) equal to 0.
c  Warning is represented by ier = 3
      z = 0.d0
      discar = discar+1
      ier = 3
      do j = 1, n 
      a(i, j) = 0.d0
      end do
      b(i) = 0.d0
      else 
	
	do j = 1, n
	a(i, j) = a(i, j)/z
	end do
	
c   Same multiplication in independent term
	b(i) = b(i)/z
      endif
	

	if(i.eq.maxrow) go to 1

      if(z.ne.0.d0) then
	do ii = i+1,  maxrow
c   Orthogonalize row ii with respect to row i
	pesca = 0.d0
	do j =1, n
	pesca = pesca + a(ii, j) * a(i, j)
	end do
	do j = 1, n
	a(ii, j) = a(ii, j) - pesca * a(i, j)
	end do
c   Same in independent term
	b(ii) = b(ii) - pesca * b(i)
	end do
      endif

	end do


c   Compute solution = Q(trasp) b
1	do j = 1, n
	x(j)= 0.d0
	do i = 1, maxrow
	x(j) = x(j) + a(i, j) * b(i)
	end do
	end do


	return
	end



 
 



      subroutine redneu(a, b, x, n, m, npsi, depth, ncomps, f, jaca,
     *   jacx, jacaux, modo, nrowsa, y, ynew, z) 

      implicit none
      integer nrowsa

      double precision a(nrowsa, nrowsa,depth), b(nrowsa,depth)
      double precision jaca(nrowsa, *), jacx(nrowsa,*), jacaux(nrowsa,*)
      double precision x(n), f(m)
      integer depth
      integer ncomps(depth)
      double precision y(nrowsa)
      integer mmax, nmax, i, nivel, j, k  
      integer n, m
      integer ncolult, ncpan
      double precision z(nrowsa,nrowsa), ynew(nrowsa)
      integer nada
      integer npsi
      integer modo
      double precision dpsi

c      if(modo.eq.1) then
c      do i = 1, n
c      y(i) = x(i)
c      end do
c      ncpan = n
c  Compute New y
c      do nivel = 1, depth
c      do i = 1, ncomps(nivel)
c      ynew(i) = b(i, nivel)
c      do j = 1, ncpan
c      ynew(i) = ynew(i) + a(i, j, nivel)*y(j)
c      end do
c      end do
c      ncpan = ncomps(nivel)
c      do i = 1, ncomps(nivel)
c      call redneupsi(npsi, ynew(i), y(i), dpsi, nivel, depth) 
c      end do           
c      end do 
 


      if(modo.eq.1) then
      do i = 1, n
      y(i) = x(i)
      end do
      ncpan = n

c  Compute New y
       do nivel = 1, depth
c      do i = 1, ncomps(nivel)
c      ynew(i) = b(i, nivel)
c      do j = 1, ncpan
c      ynew(i) = ynew(i) + a(i, j, nivel)*y(j)
c      end do
c      end do

      call mavec(a, y, ynew, ncomps(nivel), ncpan, nrowsa, nivel)  

      do i = 1, ncomps(nivel)
      ynew(i) = ynew(i) + b(i, nivel)
      end do



      ncpan = ncomps(nivel)
      do i = 1, ncomps(nivel)
      call redneupsi(npsi, ynew(i), y(i), dpsi, nivel, depth) 
      end do
      end do
           
      do i = 1, ncomps(depth)
      f(i) = y(i)
      end do
      return
      endif
c  This endif corresponds to if(modo.eq.1)
c  From now on we have modo=2, so we compute both F and Jacobian
c 


      mmax = 0
      do i = 1, depth
      mmax = max0(mmax, ncomps(i))
      end do

      nmax = n*ncomps(1)
      do i = 2, depth
      nmax = nmax + ncomps(i-1)*ncomps(i) 
      end do

      do j = 1, nmax
      do i = 1, mmax
      jaca(i, j) = 0.d0
      jacaux(i, j) = 0.d0
      end do
      end do

      do j = 1, n 
      do i = 1, n
      jacx(i, j) = 0.d0
      end do
      jacx(j, j) = 1.d0
      end do




      do nivel = 1, depth 



c      write(*, *) ' Update Jaca at level ', nivel 

c      read(*, *) nada   



      if(nivel.eq.1) then

      do i = 1, n
      y(i) = x(i)
      end do

      ncpan = n
      ncolult = 0
      else
      ncpan = ncomps(nivel-1)
      endif

c  Use chain rule in order to update pre-existing Jacobian
c------------------------------------------------------------------

       
      if(nivel.gt.1) then


c  Columns from j=1 to j = ncolult must be updated

c      write(*, *)' Level ',nivel,' Update cols 1 to ',ncolult,' of jaca'

c      read (*, *) nada

c      write(*, *)' Columns 1 to ', ncolult, ' of jaca before update:'
c      do i = 1, ncpan
c      write(*, *) (jaca(i, j), j = 1, ncolult)
c      end do 

c      read(*, *) nada

c      write(*, *)' Columns 1 to ncomps(nivel-1) of A:'
c      do i = 1, ncomps(nivel)
c      write(*, *) (a(i, j, nivel),j=1,ncpan)
c      end do

c      read (*, *) nada

c  Matrix x Matrix product jacaux = a * jaca
       call maprod(a,jaca,jacaux,ncomps(nivel),ncpan,ncolult,
     *            nrowsa,nivel)  

c      do j = 1, ncolult
c      do i = 1, ncomps(nivel)
c      jacaux(i, j) = 0.d0
c      do k = 1, ncpan
c      jacaux(i, j) = jacaux(i, j) + a(i, k, nivel)*jaca(k, j)
c      end do
c      end do
c      end do 

      do j = 1, ncolult
      do i = 1, ncomps(nivel)
      jaca(i, j) = jacaux(i, j)
      end do
      end do 
c      write(*, *)
c      write(*, *)' Columns 1 to ', ncolult, ' of jaca after update:'
c      do i = 1, ncomps(nivel)
c      write(*, *) (jaca(i, j), j = 1, ncolult)
c      end do 

c      read(*, *) nada
 


      endif 
c  End of "if (nivel .ne. 1)"
c-------------------------------------------------------------------

c-------------------------------------------------------------------
c  Add the columns of jaca that correspond to the new level

c      write(*, *)' Add the new columns of jaca from col', ncolult+1

c      read(*, *) nada

c      write(*, *) ' nivel =', nivel,' ncomps(nivel)=', ncomps(nivel)
c      write(*, *) ' ncpan = ', ncpan

      do k = 1, ncomps(nivel)
      do j = 1, (ncpan+1)*ncomps(nivel)
      jaca(k, ncolult + j)=0.d0
      end do
      end do


      do k = 1, ncomps(nivel)
      do j = 1, ncpan + 1
      if(j.lt.ncpan+1) then
      jaca(k , ncolult + (k-1)*(ncpan+1) + j ) = y(j)
      else
      jaca(k,  ncolult + (k-1)*(ncpan+1) + j) = 1.d0
      endif
      end do
      end do


c      write(*, *) ' First column typed here is number: ', ncolult+1
c      do k = 1, ncomps(nivel)
c      write(*, *) (jaca(k, j), j = ncolult+1, ncolult +
c     *     (ncomps(nivel)-1)*(ncpan+1)+ncpan+1) 
c      end do

c--------------------------------------------------------------------
c      read(*, *) nada 

      k = ncomps(nivel)
c      write(*, *)' Old ncolult:', ncolult
      ncolult =  ncolult + (k-1)*(ncpan+1)+ncpan+1
c      write(*, *)' New ncolult:', ncolult

c      read(*, *) nada

   
c--------------------------------------------------------------------
c  Update Jacobian with respect to x

       call maprod(a, jacx, z, ncomps(nivel), ncpan, n, nrowsa, nivel)  

c      do j = 1, n
c      do i = 1, ncomps(nivel)
c      z(i, j) = 0.d0
c      do k = 1, ncpan
c      z(i, j) = z(i, j) + a(i, k, nivel)*jacx(k, j)
c      end do
c      end do
c      end do

      do j = 1, n
      do i = 1, ncomps(nivel)
      jacx(i, j) = z(i, j)
      end do
      end do

c      write(*, *)' Updated Jacobian wrt x at level ', nivel,' ='
c      do i = 1, ncomps(nivel)
c      write(*, *) (jacx(i, j),j=1,n)
c      end do
c      write(*, *)
c---------------------------------------------------------------------

c      read(*, *) nada

c  Compute New y

      call mavec(a, y, ynew, ncomps(nivel), ncpan, nrowsa, nivel)  
      do i = 1, ncomps(nivel)
      ynew(i) = ynew(i) + b(i, nivel)
      end do

c      do i = 1, ncomps(nivel)
c      ynew(i) = b(i, nivel)
c      do j = 1, ncpan
c      ynew(i) = ynew(i) + a(i, j, nivel)*y(j)
c      end do
c      end do

c  Compute y at nivel and derivative of the corresponding row

      do i = 1, ncomps(nivel)
      call redneupsi(npsi, ynew(i), y(i), dpsi, nivel, depth) 
      do j = 1, n
      jacx(i, j) = jacx(i, j)*dpsi
      end do

      do j = 1, ncolult
      jaca(i, j) = jaca(i, j)*dpsi
      end do

      end do



c      write(*, *)' Updated y at level ', nivel, ' = '
c      write(*, *)(y(i), i=1, ncomps(nivel))

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c      read(*, *) nada

c      write(*, *)

      end do

c This end do corresponds to "do nivel=1,depth"
c---------------------------------------------------------------------

c      write(*, *)' Jacobian jaca is complete. Columns:', ncolult,
c     *  ' Rows:', m

c      do i = 1, m
c      write(*, *)' Row ', i 
c      write(*, *)(jaca(i, j),j=1,ncolult)
c      end do



c      write(*, *)' m = ', m,' ncomps(depth)=', ncomps(depth),' iguales'
      do i = 1, m
      f(i) = y(i)
      end do

c      write(*, *)' F =', (f(i),i=1,m)

c      read(*, *) nada
c      write(*, *)' We return from redneu'
c      read(*, *) nada
      

      return
      end


      subroutine mavec(a, x, b, m, n, mlin, nivel)
      implicit none
      integer n, m, mlin, nivel, i, j 
      double precision a(mlin, mlin, *), x(n), b(m)
c  Given a(m, n) and x(n) computes b(m)
      do i = 1, m
      b(i) = 0.d0
      end do
      do j = 1, n
      do i = 1, m
      b(i) = b(i) + a(i, j, nivel)*x(j)
      end do
      end do
      return
      end      

      subroutine maprod(a, b, c, m, n, p, mlin, nivel)
      implicit none
      integer n, m, p, nlin, mlin, nivel, i, j, k 
      double precision a(mlin, mlin, *), b(mlin, *), c(mlin, *) 
c  Given a(m, n) and b(n, p)  computes c(m, p)
      do k=1,p
      call mavec(a, b(1, k), c(1, k), m, n, mlin, nivel)
      end do
      return
      end 
      
      subroutine redneupsi(npsi, x, psi, dpsi, nivel, depth)
      implicit none
      double precision x, psi, dpsi, z
      integer npsi, nivel, depth

      if(npsi.eq.1) then
      psi = x
      dpsi = 1.d0
      return
      endif
      if(npsi.eq.2) then
c RELU
      if(x.le.0.d0) then
      psi = 0.d0
      dpsi = 0.d0
      else
      psi = x
      dpsi = 1.d0
      endif
      return
      endif

c RELU except for nivel = depth
      if(npsi.eq.3) then
      if(x.le.0.d0) then
      psi = 0.d0
      dpsi = 0.d0
      else
      psi = x
      dpsi = 1.d0
      endif
      if(nivel.eq.depth) then
      psi = x
      dpsi = 1.d0
      endif
      return
      endif


c RELU smoothed
      if(npsi.eq.4) then
      z = dsqrt(x*x + 1.d-3)
      psi = (z + x)/2.d0 
      
      dpsi =  ((x/z) + 1.d0)/2.d0
      return
      endif

c  NN-sin
      if(npsi.eq.5) then
      psi = dsin(x)
      dpsi = dcos(x)
      return
      endif

 
c  NN- 2sin
      if(npsi.eq.6) then
      psi = 2.d0 * dsin(x)
      dpsi = 2.d0* dcos(x)
      return
      endif        

      end 





                       
      subroutine rondo(seed, x)

c   Random between -1 and 1
c
C     This is the random number generator of Schrage:
C
C     L. Schrage, A more portable Fortran random number generator, ACM
C     Transactions on Mathematical Software 5 (1979), 132-138.

      double precision seed, x 

      double precision a,p,b15,b16,xhi,xalo,leftlo,fhi,k
      data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/

      xhi= seed/b16
      xhi= xhi - dmod(xhi,1.d0)
      xalo= (seed-xhi*b16)*a
      leftlo= xalo/b16
      leftlo= leftlo - dmod(leftlo,1.d0)
      fhi= xhi*a + leftlo
      k= fhi/b15
      k= k - dmod(k,1.d0)
      seed= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if (seed.lt.0) seed = seed + p
      x = seed*4.656612875d-10

      x = 2.d0*x - 1.d0

      return

      end 
            

      subroutine randin(seed, menor, mayor, sorteado)
      implicit none
      double precision seed, z
      integer menor, mayor, sorteado
c  Computes a random integer between menor and mayor
c  menor must be less than or equal to mayor
1     call rando(seed, z)
      z = dfloat(menor) + z*dfloat(mayor + 1 -menor)
      sorteado = z
      if(sorteado.lt.menor.or.sorteado.gt.mayor) go to 1
      return
      end 




                       
      subroutine rando(seed, x)

C     This is the random number generator of Schrage:
C
C     L. Schrage, A more portable Fortran random number generator, ACM
C     Transactions on Mathematical Software 5 (1979), 132-138.

      double precision seed, x 

      double precision a,p,b15,b16,xhi,xalo,leftlo,fhi,k
      data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/

      xhi= seed/b16
      xhi= xhi - dmod(xhi,1.d0)
      xalo= (seed-xhi*b16)*a
      leftlo= xalo/b16
      leftlo= leftlo - dmod(leftlo,1.d0)
      fhi= xhi*a + leftlo
      k= fhi/b15
      k= k - dmod(k,1.d0)
      seed= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if (seed.lt.0) seed = seed + p
      x = seed*4.656612875d-10

      return

      end 
             
 	subroutine chole (n, a, posdef, add)
	implicit double precision (a-h,o-z)
	logical posdef
	dimension a(n, n)
	posdef = .true.

c      write(*,*)' Matriz que entra en chole:'
c      do i = 1, n
c      write(*,*)(a(i,j),j=1,n)
c      end do

c      write(*,*)' To continue type 1'
c      read(*,*)nada


	if(a(1,1) .le. 0.d0) then
	posdef = .false.
	return
	endif

	a(1,1) = dsqrt(a(1,1))
	if(n.eq.1)return
	do 1 i=2,n
	do 2 j=1,i-1
	z = 0.d0
	if(j.gt.1)then
	do 3 k=1,j-1
	z = z + a(i,k) * a(j,k)
3     continue
	endif
	a(i,j) = (a(i,j) - z)/a(j,j)
2	continue
	z = 0.d0
	do 4 j=1,i-1
	z = z + a(i,j)**2
4     continue
	if( a(i,i) - z .le. 0.d0) then
	posdef = .false.
	add = z - a(i,i)
	return
	endif

	a(i,i) = dsqrt( a(i,i) - z )
1	continue
	return
	end



 	subroutine sicho (n, a, x, b, aux)
	implicit double precision (a-h,o-z)
	dimension a(n,n),x(n),b(n),aux(n)
	aux(1) = b(1)/a(1,1)
	if(n.gt.1)then
	do 1 i=2,n
	z = 0.d0
	do 2 j=1,i-1
	z = z + a(i,j)*aux(j)
2     continue
	aux(i) = (b(i) - z) / a(i,i)
1	continue
	endif
	x(n) = aux(n)/a(n,n)
	if(n.eq.1)return
	do 3 i=n-1,1,-1
	z = 0.d0
	do 4 j=i+1,n
	z = z + a(j,i)*x(j)
4     continue
	x(i) = (aux(i) - z)/a(i,i)
3	continue
	return
	end


 



      
