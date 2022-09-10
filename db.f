      program dieletricbreakdown
      real*8, dimension(:,:), allocatable :: phi
      integer, dimension(:,:), allocatable :: ground
      real*8, dimension(:), allocatable ::  pot
      integer, dimension(:), allocatable :: xx, yy, prox
      common eta, n, nn
      
C========================================================================
C Il codice simula un fenomeno di rottura dielettrica su griglia quadrata
C Se il parametro i_anim vale 0 allora viene eseguita la simulazione
C salvando solo alla fine, dopo npoints iterazioni, la griglia su file.
C Se invece i_anim assume valore 1, allora viene salvata una griglia
C ogni i_delay iterazioni per un totale di i_time griglie, in modo da poi
C leggerle da un file python ed eseguire un'animazione della curva
C========================================================================
      
      call cpu_time(start)
      call ranstart
	
      open(2, file='datidb.dat',status='unknown') !file dove salvare i risultati
      
      n = 200                     !dimesione griglia
      nn = n**2                   !area della griglia  
      i_anim = 0                  !parametro per l'animalzione
      i_time = 200                !tempo dell'animazione
      i_delay = 10                !iterazioni saltate
      npoints = i_time*i_delay    !tempo totale
      eta = 1                     !parametro del modello
      m = (n+1)/2                 !punto centrale della griglia
      
      !scrivo su file le quantita importanti per il plot
      write(2,*) n
      write(2,*) npoints
      write(2,*) i_delay
      write(2,*) i_anim
      
      !gound è la matrice che contiene la curva di rottura
      !phi è la matrice che continene il potenziale della griglia
      allocate(ground(n,n), phi(n,n))
      
      !inizializzo le matrici
      do i = 1, n
          do j = 1, n
              ground(i, j) = 1
              phi(i, j) = 1.0
          enddo
      enddo
      
      !seme della routtura
      ground(m, m) = 0
      
      !rislouzione equazione di laplace una prima volta
      call solve(3000, ground, phi)
      
      do iter1 = 1, i_time
          do iter2 = 1, i_delay
          
              allocate(xx(nn), yy(nn), prox(nn), pot(nn))
C=================================================
C I seguenti array contengono rispettivamente:
C xx : coordinate x dei primi vicini alla curva
C yy : coordinate y dei primi vicini alla curva
C prox : punti primi vicini alla curva
C pot : distribuzione da campionare per scegliere
C       quale punto aggiungere alla curva
C Per sicurezza sono allocati lunghi n^2, dato che
C non è possibile sapere a priori quanto saranno
C lunghi ad ogni iterazione della simulazione 
C=================================================
                   
              call break(ground, phi, xx, yy, prox, pot, l)
              !considero solo la parte contenete informazioni
              pot = pot(:l)
              prox = prox(:l)
              xx = xx(:l)
              yy = yy(:l)
          
              !distribuzione da campionare
              pot = pot**eta
              norm = sum(pot)
              pot = pot/float(norm)
          
              call choice(l, prox, pot, ic) !campiono la distribuzione
                 
              ground(xx(ic), yy(ic)) = 0    !aggiungo un punto alla curva
          
              call solve(100, ground, phi)  !rislovo con le nuove condizioni al bordo
          
              deallocate(xx, yy, prox, pot)
          enddo
          
          !salvo una grilia ogni i_delay iterazioni
          if (i_anim == 1) then
              do i =1, n
                  do j = 1, n
                      write(2,*) ground(i, j)
                  enddo
              enddo
          endif
                    
      enddo
      
      !salvo solo la griglia finale
      if (i_anim == 0) then
          do i =1, n
              do j = 1, n
                  write(2,*) ground(i, j)
              enddo
          enddo
      endif
      
      call ranfinish
	
      call cpu_time(finish)
      print '("tempo di esecuzione= ", f16.8," secondi.")', finish-start
      
           
      end program
      
      
C======================================================================

      subroutine solve(iters, ground, phi)
      common eta, n, nn
      real*8, dimension(n,n) :: phi
      integer, dimension(n,n) :: ground
      
      !condizioni al bordo
      do i = 2, n-1
          do j = 2, n-1
              if (ground(i, j)==0) phi(i, j)=0
          enddo
      enddo
      
      !risoluzione tramite teorma del valor medio
      do k = 1, iters
          do i = 2, n-1
              do j = 2, n-1
                  force = phi(i, j+1) + phi(i, j-1) +
     &                    phi(i+1, j) + phi(i-1, j)
     
                  if (ground(i, j)==0) then
                      phi(i, j)=0
                  else 
                      phi(i, j) = force/4.0
                  endif
              enddo
          enddo
      enddo

      return
      end
     
C========================================================================

      subroutine break(ground, phi, xx, yy, prox, pot, l)
      common eta, n, nn
      real*8, dimension(n,n) :: phi
      integer, dimension(n,n) :: ground
      real*8, dimension(nn) ::  pot
      integer, dimension(nn) :: xx, yy, prox
      logical :: sn, ds
      
      !inizializzo variabili
      l = 1
      prox = 0
      xx = 0
      yy = 0
      pot = 0
      
      !ciclo sulla griglia per cacoloco dei primi vicini e distribuzione
      do i = 2, n-1
          do j = 2, n-1
              if (phi(i, j) /= 0) then
                  
                  sn = (ground(i,j+1)==0).or.(ground(i,j-1)==0)
                  ds = (ground(i+1,j)==0).or.(ground(i-1,j)==0)
                  
                  if (sn.or.ds) then
                      c = 0
                      if (ground(i+1,j)== 0) c = c+1
                      if (ground(i-1,j)== 0) c = c+1
                      if (ground(i,j+1)== 0) c = c+1
                      if (ground(i,j-1)== 0) c = c+1
                      
          	      prox(l) = l
          	      xx(l) = i
         	      yy(l) = j
                      pot(l) = c**(-eta)*phi(i,j)
                      l = l + 1
                      
                  endif
              endif
          enddo
      enddo
      
      !riaggiusto per ottenere la lunghezza degli array
      l = l - 1
      
      return
      end
      
C========================================================================

      subroutine choice(l, v, p, ic)
      common eta, n, nn
      real*8, dimension(l) :: p, cdf
      integer, dimension(l) :: v
      
      !calcolo della funzione cumulativa
      cdf = 0
      cdf(1) = p(1)
      
      do i = 2, l
          cdf(i) = cdf(i-1) + p(i)
      enddo
      cdf = cdf/(cdf(l)*1.0)
      
      x = ran2() !genero variabile casuale uniforme in [0,1)
      
      !inverisione della cdf, k sarà distribuito secondo p
      k = 1
      do j = 2, l
          if ((x<cdf(j)).and.(x>=cdf(j-1))) k = j
      enddo
      
      ic = v(k)     
      
      return
      end

C============================================================================

c============================================================================
c  RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
c============================================================================
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     &          ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &          ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,
     &          rnmx=1.-eps)
      integer idum2,j,k,iv,iy
      common /dasav/ idum,idum2,iv(ntab),iy
c      save iv,iy,idum2
c      data idum2/123456789/, iv/NTAB*0/, iy/0/

      if(idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
      end

c=============================================================================
      subroutine ranstart
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      read(23,*) idum
      read(23,*,end=117) idum2
      do i=1,32
         read(23,*) iv(i)
      enddo
      read(23,*) iy
      close(23)
      goto 118                          !!takes account of the first start
 117  if(idum.ge.0) idum = -idum -1     !!
      close(23)
 118  continue                          !!

      return
      end

c=============================================================================
      subroutine ranfinish
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      write(23,*) idum
      write(23,*) idum2
      do i=1,32
         write(23,*) iv(i)
      enddo
      write(23,*) iy
      close(23)

      return
      end
c=============================================================================

