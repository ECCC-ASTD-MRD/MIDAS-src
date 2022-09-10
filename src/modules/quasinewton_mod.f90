!---------------------------------------LICENCE BEGIN -----------------------------------
!     Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!     version 3; Last Modified: May 7, 2008.
!     This is free but copyrighted software; you can use/redistribute/modify it under the terms
!     of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!     version 3 or (at your option) any later version that should be found at:
!     http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!     
!     This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!     without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!     See the above mentioned License/Disclaimer for more details.
!     You should have received a copy of the License/Disclaimer along with this software;
!     if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
  ! CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!--------------------------------------LICENCE END --------------------------------------


module quasinewton_mod
  ! MODULE quasinewton_mod (prefix='qna' category='1. High-level functionality')
  !
  ! Purpose: The n1qn3 routine, and its supporting subroutines
  !
      use midasMpi_mod
      implicit none
      save 
      private

      public :: qna_n1qn3

      contains
      
      subroutine dcube(t,f,fp,ta,fa,fpa,tlower,tupper)
!
!:Purpose: [to be completed]
!

!     
!     --- arguments
!     
      real(8) sign,den,anum,t,f,fp,ta,fa,fpa,tlower,tupper
!     
!     --- variables locales
!     
      real(8) z1,b,discri
!     
!     Using f and fp at t and ta, computes new t by cubic formula
!     safeguarded inside [tlower,tupper].
!     
      z1=fp+fpa-3.d0*(fa-f)/(ta-t)
      b=z1+fp
!     
!     first compute the discriminant (without overflow)
!     
      if (dabs(z1).le.1.d0) then
         discri=z1*z1-fp*fpa
      else
         discri=fp/z1
         discri=discri*fpa
         discri=z1-discri
         if (z1.ge.0.d0 .and. discri.ge.0.d0) then
            discri=dsqrt(z1)*dsqrt(discri)
            go to 120
         endif
         if (z1.le.0.d0 .and. discri.le.0.d0) then
            discri=dsqrt(-z1)*dsqrt(-discri)
            go to 120
         endif
         discri=-1.d0
      endif
      if (discri.lt.0.d0) then
         if (fp.lt.0.d0) t=tupper
         if (fp.ge.0.d0) t=tlower
         go to 900
      endif
!     
!     discriminant nonnegative, compute solution (without overflow)
!     
      discri=dsqrt(discri)
 120  if (t-ta.lt.0.d0) discri=-discri
      sign=(t-ta)/dabs(t-ta)
      if (b*sign.gt.0.d+0) then
         t=t+fp*(ta-t)/(b+discri)
      else
         den=z1+b+fpa
         anum=b-discri
         if (dabs((t-ta)*anum).lt.(tupper-tlower)*dabs(den)) then
            t=t+anum*(ta-t)/den
         else
            t=tupper
         endif
      endif
 900  t=dmax1(t,tlower)
      t=dmin1(t,tupper)
!     return
      end subroutine

      subroutine ddd (prosca,dtonb,dtcab,n,sscale,nm,depl,aux,jmin,jmax, &
                      precos,diag,ybar,sbar,izs,rzs,dzs)
      !
      !:Purpose: Calcule le produit H.g ou
      !
      !          * H est une matrice construite par la formule de bfgs inverse
      !            a nm memoires a partir de la matrice diagonale diag
      !            dans un espace hilbertien dont le produit scalaire
      !            est donne par prosca
      !            (cf. J. Nocedal, Math. of Comp. 35/151 (1980) 773-782)
      !          * g est un vecteur de dimension n (en general le gradient)
      !
      !          * la matrice diag apparait donc comme un preconditionneur diagonal
      !
      !          * depl = g (en entree), = H g (en sortie)
      !
      !          * la matrice H est memorisee par les vecteurs des tableaux
      !            ybar, sbar et les pointeurs jmin, jmax
      !
      !          * alpha(nm) est une zone de travail
      !
      !          * izs(1),rzs(1),dzs(1) sont des zones de travail pour prosca
      !

      ! arguments
      logical sscale
      integer n,nm,jmin,jmax,izs(1)
      real rzs(1)
      real(8) depl(n),precos,diag(n),alpha(nm),ybar(n,1), &
          sbar(n,1),aux(n),dzs(1)
      external prosca,dtonb,dtcab
!
!         variables locales
!
      integer jfin,i,j,jp
      real(8) r,ps
!
      jfin=jmax
      if (jfin.lt.jmin) jfin=jmax+nm
!
!         phase de descente
!
      do 100 j=jfin,jmin,-1
          jp=j
          if (jp.gt.nm) jp=jp-nm
          call prosca (n,depl,sbar(1,jp),ps,izs,rzs,dzs)
          r=ps
          alpha(jp)=r
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) 
          do 20 i=1,n
              depl(i)=depl(i)-r*ybar(i,jp)
20        continue
!$OMP END PARALLEL DO
100   continue
!
!         preconditionnement
!
      if (sscale) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) 
          do 150 i=1,n
              depl(i)=depl(i)*precos
  150     continue
!$OMP END PARALLEL DO
      else
          call dtonb (n,depl,aux,izs,rzs,dzs)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) 
          do 151 i=1,n
              aux(i)=aux(i)*diag(i)
  151     continue
!$OMP END PARALLEL DO
          call dtcab (n,aux,depl,izs,rzs,dzs)
      endif
!
!         remontee
!
      do 200 j=jmin,jfin
          jp=j
          if (jp.gt.nm) jp=jp-nm
          call prosca (n,depl,ybar(1,jp),ps,izs,rzs,dzs)
          r=alpha(jp)-ps
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) 
          do 120 i=1,n
              depl(i)=depl(i)+r*sbar(i,jp)
120       continue
!$OMP END PARALLEL DO
200   continue

      end subroutine ddd

      subroutine ddds (prosca,dtonb,dtcab,n,sscale,nm,depl,aux,jmin, &
                       jmax,precos,diag,ybar,sbar,izs,rzs,dzs)
      !
      !
      !:Purpose: This subroutine has the same role as ddd (computation of the
      !          product H.g). It supposes however that the (y,s) pairs are not
      !          stored in core memory, but on a devise chosen by the user.
      !

!
!         arguments
!
      logical sscale
      integer n,nm,jmin,jmax,izs(1)
      real rzs(1)
      real(8) depl(n),precos,diag(n),alpha(nm),ybar(n),sbar(n), &
          aux(n),dzs(1)
      external prosca,dtonb,dtcab
!
!         variables locales
!
      integer jfin,i,j,jp
      real(8) r,ps
!
      jfin=jmax
      if (jfin.lt.jmin) jfin=jmax+nm
!
!         phase de descente
!
      do 100 j=jfin,jmin,-1
          jp=j
          if (jp.gt.nm) jp=jp-nm
          call prosca (n,depl,sbar,ps,izs,rzs,dzs)
          r=ps
          alpha(jp)=r
          do 20 i=1,n
              depl(i)=depl(i)-r*ybar(i)
20        continue
100   continue
!
!         preconditionnement
!
      if (sscale) then
          do 150 i=1,n
              depl(i)=depl(i)*precos
  150     continue
      else
          call dtonb (n,depl,aux,izs,rzs,dzs)
          do 151 i=1,n
              aux(i)=aux(i)*diag(i)
  151     continue
          call dtcab (n,aux,depl,izs,rzs,dzs)
      endif
!
!         remontee
!
      do 200 j=jmin,jfin
          jp=j
          if (jp.gt.nm) jp=jp-nm
          call prosca (n,depl,ybar(1),ps,izs,rzs,dzs)
          r=alpha(jp)-ps
          do 120 i=1,n
              depl(i)=depl(i)+r*sbar(i)
120       continue
200   continue

      end subroutine ddds

      subroutine mupdts (sscale,inmemo,n,m,nrz)

!
!         arguments
!
      logical sscale,inmemo
      integer n,m,nrz
!----
!
!     On entry:
!       sscale: .true. if scalar initial scaling,
!               .false. if diagonal initial scaling
!       n:      number of variables
!
!     This routine has to return:
!       m:      the number of updates to form the approximate Hessien H,
!       inmemo: .true., if the vectors y and s used to form H are stored
!                  in core memory,
!               .false. otherwise (storage of y and s on disk, for
!                  instance).
!
!----
!
      if (n.eq.0) then
         m=0
      else
         if (sscale) then
            m=(nrz-3*n)/(2*n)
         else
            m=(nrz-4*n)/(2*n)
         endif
      endif
      inmemo=.true.
      ! return
      end subroutine mupdts

      subroutine qna_n1qn3 (simul,prosca,dtonb,dtcab,n,x,f,g,dxmin,df1, &
           epsg,impres,io,mode,niter,nsim,iz,dz,ndz, &
           izs,rzs,dzs)
      !     
      !:Purpose: N1QN3, Version 2.0c, June 1995
      !          Jean Charles Gilbert, Claude Lemarechal, INRIA.
      !     
      !          Double precision version of M1QN3.
      !     
      !          N1qn3 has two running modes: the SID (Scalar Initial Scaling) mode
      !          and the DIS (Diagonal Initial Scaling) mode. Both do not require
      !          the same amount of storage, the same subroutines, ...
      !          In the description below, items that differ in the DIS mode with
      !          respect to the SIS mode are given in brakets.
      !     
      !          Use the following subroutines:
      !
      !          * N1QN3A
      !          * DDD, DDDS
      !          * NLIS0 + DCUBE (Dec 88)
      !          * MUPDT
      !     
      !          The following routines are proposed to the user in case the
      !          Euclidean scalar product is used:
      !          DUCLID, DTONBE, DTCABE.
      !     
      !          La sous-routine N1QN3 est une interface entre le programme
      !          appelant et la sous-routine N1QN3A, le minimiseur proprement dit.
      !     
      !          Le module PROSCA est sense realiser le produit scalaire de deux
      !          vecteurs de Rn; le module DTONB est sense realiser le changement
      !          de coordonnees correspondant au changement de bases: base
      !          euclidienne -> base orthonormale (pour le produit scalaire
      !          PROSCA); le module CTBAB fait la transformation inverse: base
      !          orthonormale -> base euclidienne.
      !     
      !          Iz is an integer working zone for N1QN3A, its dimension is 5.
      !          It is formed of 5 scalars that are set by the optimizer:
      !
      !          - the dimension of the problem,
      !          - a identifier of the scaling mode,
      !          - the number of updates,
      !          - two pointers.
      !     
      !          Dz est la zone de travail pour N1QN3A, de dimension ndz.
      !          Elle est subdivisee en
      !
      !          - 3 [ou 4] vecteurs de dimension n: d,gg,[diag,]aux
      !          - m vecteurs de dimension n: ybar
      !          - m vecteurs de dimension n: sbar
      !     
      !          m est alors le plus grand entier tel que
      !
      !          - m*(2*n+1)+3*n .le. ndz [m*(2*n+1)+4*n .le. ndz)]
      !          - soit m := (ndz-3*n) / (2*n+1) [m := (ndz-4*n) / (2*n+1)].
      !          Il faut avoir m >= 1, donc ndz >= 5n+1 [ndz >= 6n+1].
      !     
      !          A chaque iteration la metrique est formee a partir d'un multiple
      !          de l'identite [d'une matrice diagonale] D qui est mise a jour m
      !          fois par la formule de BFGS en utilisant les m couples {y,s} les
      !          plus recents.
      !     

!     
!     arguments
!     
      integer n,impres,io,mode,niter,nsim,iz(5),ndz,izs(1)
      real rzs(1)
      real(8) x(n),f,g(n),dxmin,df1,epsg,dz(ndz),dzs(1)
      external simul,prosca,dtonb,dtcab
!     
!     variables locales
!     
      logical inmemo,sscale
      integer ntravu,id,igg,idiag,iaux,iybar,isbar,m,mmemo,ntotal,ierr
      integer m_max
      real(8) d1,d2,ps
!     
!---- impressions initiales et controle des arguments
!     
      write(io,*) '--------------------------------------------------'
      write(io,*) 'N1QN3: calling modified MPI version of modulopt!!!'
      write(io,*) '--------------------------------------------------'

      call rpn_comm_allreduce(n,ntotal,1,"mpi_integer", &
                              "mpi_max","GRID",ierr)

      if (impres.ge.1) &
           write (io,900) n,dxmin,df1,epsg,niter,nsim,impres
 900  format (/" N1QN3 (Version 2.0c, June 1995): entry point"/ &
           5x,"dimension of the problem (n):",i9/  &
           5x,"absolute precision on x (dxmin):",d9.2/ &
           5x,"expected decrease for f (df1):",d9.2/ &
           5x,"relative precision on g (epsg):",d9.2/ &
           5x,"maximal number of iterations (niter):",i6/ &
           5x,"maximal number of simulations (nsim):",i6/ &
           5x,"printing level (impres):",i4)
      if (ntotal.le.0.or.niter.le.0.or.nsim.le.0.or.dxmin.le.0.d+0.or. &
           epsg.le.0.d+0.or.epsg.gt.1.d+0.or.mode.lt.0.or.  &
           mode.gt.3) then
         mode=2
         write (io,901)
 901     format (/" >>> n1qn3: inconsistent call")
         return
      endif
!     
!---- what method
!     
      if (mod(mode,2).eq.0) then
         if (impres.ge.1) write (io,920)
 920     format (/" n1qn3: Diagonal Initial Scaling mode")
         sscale=.false.
      else
         if (impres.ge.1) write (io,921)
 921     format (/" n1qn3: Scalar Initial Scaling mode")
         sscale=.true.
      endif
!     
      if ( n.gt.0 .and. &
           ((ndz.lt.5*n+1).or.((.not.sscale).and.(ndz.lt.6*n+1))) ) then
         mode=2
         write (io,902)
 902     format (/" >>> n1qn3: not enough memory allocated")
         return
      endif
!     
!---- Compute m
!     
      call mupdts (sscale,inmemo,n,m,ndz)
      call rpn_comm_allreduce(m,m_max,1,"mpi_integer", &
                              "mpi_max","GRID",ierr)
      if (m.ne.m_max) then
         write(io,*) 'replacing value of m, ',m,' with ',m_max
         m=m_max
      endif
!     
!     --- Check the value of m (if (y,s) pairs in core, m will be >= 1)
!     
      if (m.lt.1) then
         mode=2
         write (io,9020)
 9020    format (/" >>> n1qn3: m is set too small in mupdts")
         return
      endif
!     
!     --- mmemo = number of (y,s) pairs in core memory
!     
      mmemo=1
      if (inmemo) mmemo=m
!     
      ntravu=2*(2+mmemo)*n
      if (sscale) ntravu=ntravu-n
      write (io,903) ndz,ntravu,m
 903  format (/5x,"allocated memory (ndz) :",i9/ &
           5x,"used memory :           ",i9/ &
           5x,"number of updates :     ",i9)
      if (ndz.lt.ntravu) then
         mode=2
         write (io,902)
         return
      endif
!     
      if (impres.ge.1) then
         if (inmemo) then
            write (io,907)
         else
            write (io,908)
         endif
      endif
 907  format (5x,"(y,s) pairs are stored in core memory")
 908  format (5x,"(y,s) pairs are stored by the user")
!     
!---- cold start or warm restart ?
!     check iz: iz(1)=n, iz(2)=(0 if DIS, 1 if SIS),
!     iz(3)=m, iz(4)=jmin, iz(5)=jmax
!     
      if (mode/2.eq.0) then
         if (impres.ge.1) write (io,922)
      else
         iaux=0
         if (sscale) iaux=1
!         if (iz(1).ne.n.or.iz(2).ne.iaux.or.iz(3).ne.m.or.iz(4).lt.1 &
!              .or.iz(5).lt.1.or.iz(4).gt.iz(3).or.iz(5).gt.iz(3)) then
         if (iz(2).ne.iaux.or.iz(3).ne.m.or.iz(4).lt.1 &
              .or.iz(5).lt.1.or.iz(4).gt.iz(3).or.iz(5).gt.iz(3)) then
            mode=2
            write(*,*) 'iz=',iz(:)
            write(*,*) 'n,iaux,m=',n,iaux,m
            if (impres.ge.1) write (io,923)
            return
         endif
         if (impres.ge.1) write (io,924)
      endif
 922  format (/" n1qn3: cold start"/,1x)
 923  format (/" >>> n1qn3: inconsistent iz for a warm restart")
 924  format (/" n1qn3: warm restart"/,1x)
      iz(1)=n
      iz(2)=0
      if (sscale) iz(2)=1
      iz(3)=m
!     
!---- split the working zone dz
!     
      idiag=1
      iybar=idiag+n
      if (sscale) iybar=1
      isbar=iybar+n*mmemo
      id=isbar+n*mmemo
      igg=id+n
      iaux=igg+n
!     
!---- call the optimization code
!     
      call n1qn3a (simul,prosca,dtonb,dtcab,n,x,f,g,dxmin,df1,epsg, &
           impres,io,mode,niter,nsim,inmemo,iz(3),iz(4),iz(5), &
           dz(id),dz(igg),dz(idiag),dz(iaux), &
           dz(iybar),dz(isbar),izs,rzs,dzs)
!     
!---- impressions finales
!     
      if (impres.ge.1) write (io,905) mode,niter,nsim,epsg
 905  format (/,1x,79("-")/  &
           /" n1qn3: output mode is ",i2 &
           /5x,"number of iterations: ",i4 &
           /5x,"number of simulations: ",i6 &
           /5x,"realized relative precision on g: ",d9.2)
      call prosca (n,x,x,ps,izs,rzs,dzs)
      d1=dsqrt(ps)
      call prosca (n,g,g,ps,izs,rzs,dzs)
      d2=dsqrt(ps)
      if (impres.ge.1) write (io,906) d1,f,d2
 906  format (5x,"norm of x = ",d15.8 &
           /5x,"f         = ",d15.8 &
           /5x,"norm of g = ",d15.8) 
      !return
      end subroutine qna_n1qn3

      subroutine n1qn3a (simul,prosca,dtonb,dtcab,n,x,f,g,dxmin,df1, &
                         epsg,impres,io,mode,niter,nsim,inmemo,m,jmin, &
                         jmax,d,gg,diag,aux,ybar,sbar,izs,rzs,dzs)
     !
     !:Purpose: Code d'optimisation proprement dit.
     !

!
!         arguments
!
      logical inmemo
      integer n,impres,io,mode,niter,nsim,m,jmin,jmax,izs(1)
      real rzs(1)
      real(8) x(n),f,g(n),dxmin,df1,epsg,d(n),gg(n),diag(n), &
          aux(n),ybar(n,1),sbar(n,1),dzs(1)
      external simul,prosca,dtonb,dtcab
!
!         variables locales
!
      logical sscale,cold,warm,ntotal
      integer i,itmax,moderl,isim,jcour,indic,ierr,impresmax
      real(8) d1,t,tmin,tmin_mpiglobal,tmax,gnorm,eps1,ff, &
           preco,precos,ys,den, &
           dk,dk1,ps,ps2,hp0
!
!         parametres
!
      real(8) rm1,rm2
      parameter (rm1=0.0001d+0,rm2=0.9d+0)
      real(8) pi
      parameter (pi=3.1415927d+0)
      real(8) rmin
!
!---- initialisation
!
      rmin=1.d-20
!
      sscale=.true.
      if (mod(mode,2).eq.0) sscale=.false.
!
      warm=.false.
      if (mode/2.eq.1) warm=.true.
      cold=.not.warm
!
      itmax=niter
      niter=0
      isim=1
      eps1=1.d+0
!
      call rpn_comm_allreduce(impres,impresmax,1,"mpi_integer", &
                              "mpi_max","GRID",ierr)
      call rpn_comm_allreduce(n,ntotal,1,"mpi_integer", &
                              "mpi_max","GRID",ierr)
!
      call prosca (n,g,g,ps,izs,rzs,dzs)
      gnorm=dsqrt(ps)
      if (impres.ge.1) write (io,900) f,gnorm
  900 format (5x,"f         = ",d15.8 &
             /5x,"norm of g = ",d15.8)
      if (gnorm.lt.rmin) then
          mode=2
          if (impres.ge.1) write (io,901)
          goto 1000
      endif
  901 format (/" >>> n1qn3a: initial gradient is too small")
!
!     --- initialisation pour ddd
!
      if (cold) then
          jmin=1
          jmax=0
      endif
      jcour=1
      if (inmemo) jcour=jmax
!
!     --- mise a l'echelle de la premiere direction de descente
!
      if (cold) then
!
!         --- use Fletcher's scaling and initialize diag to 1.
!
          precos=2.d+0*df1/gnorm**2
          do 10 i=1,n
              d(i)=-g(i)*precos
              diag(i)=1.d+0
   10     continue
          if (impres.ge.5) write(io,902) precos
  902     format (/" n1qn3a: descent direction -g: precon = ",d10.3)
      else
!
!         --- use the matrix stored in [diag and] the (y,s) pairs
!
          if (sscale) then
              call prosca (n,ybar(1,jcour),ybar(1,jcour),ps,izs,rzs,dzs)
              precos=1.d+0/ps
          endif
          do 11 i=1,n
              d(i)=-g(i)
  11      continue
          if (inmemo) then
              call ddd (prosca,dtonb,dtcab,n,sscale,m,d,aux,jmin,jmax, &
                       precos,diag,ybar,sbar,izs,rzs,dzs)
          else
              call ddds (prosca,dtonb,dtcab,n,sscale,m,d,aux,jmin,jmax, &
                         precos,diag,ybar,sbar,izs,rzs,dzs)
          endif
      endif
!
      if (impres.eq.3) then
          write(io,903)
          write(io,904)
      endif
      if (impres.eq.4) write(io,903)
  903 format (/,1x,79("-"))
  904 format (1x)
!
!     --- initialisation pour nlis0
!
      tmax=1.d+20
      call prosca (n,d,g,hp0,izs,rzs,dzs)
      if (hp0.ge.0.d+0) then
          mode=7
          if (impres.ge.1) write (io,905) niter,hp0
          goto 1000
      endif
  905 format (/" >>> n1qn3 (iteration ",i2,"): " &
              /5x," the search direction d is not a ", &
               "descent direction: (g,d) = ",d12.5)
!
!     --- compute the angle (-g,d)
!
      if (warm.and.impresmax.ge.5) then
          call prosca (n,g,g,ps,izs,rzs,dzs)
          ps=dsqrt(ps)
          call prosca (n,d,d,ps2,izs,rzs,dzs)
          ps2=dsqrt(ps2)
          ps=hp0/ps/ps2
          ps=dmin1(-ps,1.d+0)
          ps=dacos(ps)
          d1=ps*180.d+0/pi
          if(impres.ge.5) write (io,906) sngl(d1)
      endif
  906 format (/" n1qn3: descent direction d: ", &
              "angle(-g,d) = ",f5.1," degrees")
!
!---- Debut de l'iteration. on cherche x(k+1) de la forme x(k) + t*d,
!     avec t > 0. On connait d.
!
!         Debut de la boucle: etiquette 100,
!         Sortie de la boucle: goto 1000.
!
100   niter=niter+1
      if (impres.lt.0) then
          if (mod(niter,-impres).eq.0) then
              indic=1
              call simul (indic,n,x,f,g,izs,rzs,dzs)
          endif
      endif
      if (impres.ge.5) write(io,903)
      if (impres.ge.4) write(io,904)
      if (impres.ge.3) write (io,910) niter,isim,f,hp0
  910 format (" n1qn3: iter ",i3,", simul ",i3, &
              ", f=",d15.8,", h'(0)=",d12.5)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) 
      do 101 i=1,n
          gg(i)=g(i)
101   continue
!$OMP END PARALLEL DO
      ff=f
!
!     --- recherche lineaire et nouveau point x(k+1)
!
      if (impres.ge.5) write (io,911)
  911 format (/" n1qn3: line search")
!
!         --- calcul de tmin
!
      tmin=0.d+0
      do 200 i=1,n
          tmin=max(tmin,dabs(d(i)))
200   continue
      call rpn_comm_allreduce(tmin,tmin_mpiglobal,1, &
                              "mpi_double_precision", &
                              "mpi_max","GRID",ierr)
      tmin = tmin_mpiglobal

      tmin=dxmin/tmin
      t=1.d+0
      d1=hp0
!
      call nlis0 (n,simul,prosca,x,f,d1,t,tmin,tmax,d,g,rm2,rm1, &
                 impres,io,moderl,isim,nsim,aux,izs,rzs,dzs)
!
!         --- nlis0 renvoie les nouvelles valeurs de x, f et g
!
      if (moderl.ne.0) then
          if (moderl.lt.0) then
!
!             --- calcul impossible
!                 t, g: ou les calculs sont impossibles
!                 x, f: ceux du t_gauche (donc f <= ff)
!
              mode=moderl
          elseif (moderl.eq.1) then
!
!             --- descente bloquee sur tmax
!                 [sortie rare (!!) d'apres le code de nlis0]
!
              mode=3
              if (impres.ge.1) write(io,912) niter
  912         format (/" >>> n1qn3 (iteration ",i3, &
                      "): line search blocked on tmax: ", &
                      "decrease the scaling")
          elseif (moderl.eq.4) then
!
!             --- nsim atteint
!                 x, f: ceux du t_gauche (donc f <= ff)
!
              mode=5
          elseif (moderl.eq.5) then
!
!             --- arret demande par l'utilisateur (indic = 0)
!                 x, f: ceux en sortie du simulateur
!
              mode=0
          elseif (moderl.eq.6) then
!
!             --- arret sur dxmin ou appel incoherent
!                 x, f: ceux du t_gauche (donc f <= ff)
!
              mode=6
          endif
          goto 1000
      endif
!
! NOTE: stopping tests are now done after having updated the matrix, so
! that update information can be stored in case of a later warm restart
!
!     --- mise a jour de la matrice
!
      if (m.gt.0) then
!
!         --- mise a jour des pointeurs
!
          jmax=jmax+1
          if (jmax.gt.m) jmax=jmax-m
          if ((cold.and.niter.gt.m).or.(warm.and.jmin.eq.jmax)) then
              jmin=jmin+1
              if (jmin.gt.m) jmin=jmin-m
          endif
          if (inmemo) jcour=jmax
!
!         --- y, s et (y,s)
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) 
          do 400 i=1,n
              sbar(i,jcour)=t*d(i)
              ybar(i,jcour)=g(i)-gg(i)
400       continue
!$OMP END PARALLEL DO
          if (impresmax.ge.5) then
              call prosca (n,sbar(1,jcour),sbar(1,jcour),ps,izs,rzs,dzs)
              dk1=dsqrt(ps)
              if (impres.ge.5.and.niter.gt.1) write (io,930) dk1/dk
  930         format (/" n1qn3: convergence rate, s(k)/s(k-1) = ", &
                      d12.5)
              dk=dk1
          endif
          call prosca (n,ybar(1,jcour),sbar(1,jcour),ys,izs,rzs,dzs)
          if (ys.le.0.d+0) then
              mode=7
              if (impres.ge.1) write (io,931) niter,ys
  931         format (/" >>> n1qn3 (iteration ",i2, &
                      "): the scalar product (y,s) = ",d12.5 &
                      /27x,"is not positive")
              goto 1000
          endif
!
!         --- ybar et sbar
!
          d1=dsqrt(1.d+0/ys)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) 
          do 410 i=1,n
              sbar(i,jcour)=d1*sbar(i,jcour)
              ybar(i,jcour)=d1*ybar(i,jcour)
  410     continue
!$OMP END PARALLEL DO
!
!         --- compute the scalar or diagonal preconditioner
!
          if (impres.ge.5) write(io,932)
  932     format (/" n1qn3: matrix update:")
!
!             --- Here is the Oren-Spedicato factor, for scalar scaling
!
          if (sscale) then
              call prosca (n,ybar(1,jcour),ybar(1,jcour),ps,izs,rzs,dzs)
              precos=1.d+0/ps
!
              if (impres.ge.5) write (io,933) precos
  933         format (5x,"Oren-Spedicato factor = ",d10.3)
!
!             --- Scale the diagonal to Rayleigh's ellipsoid.
!                 Initially (niter.eq.1) and for a cold start, this is
!                 equivalent to an Oren-Spedicato scaling of the
!                 identity matrix.
!
          else
              call dtonb (n,ybar(1,jcour),aux,izs,rzs,dzs)
              ps=0.d0
              do 420 i=1,n
                  ps=ps+diag(i)*aux(i)*aux(i)
  420         continue
              call mmpi_allreduce_sumreal8scalar(ps,"GRID")
              d1=1.d0/ps
              if (impres.ge.5) then
                  write (io,934) d1
  934             format(5x,"fitting the ellipsoid: factor = ",d10.3)
              endif
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) 
              do 421 i=1,n
                  diag(i)=diag(i)*d1
  421         continue
!$OMP END PARALLEL DO
!
!             --- update the diagonal
!                 (gg is used as an auxiliary vector)
!
              call dtonb (n,sbar(1,jcour),gg,izs,rzs,dzs)
              ps=0.d0
              do 430 i=1,n
                  ps=ps+gg(i)*gg(i)/diag(i)
  430         continue
              call mmpi_allreduce_sumreal8scalar(ps,"GRID")
              den=ps
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) 
              do 431 i=1,n
                  diag(i)=1.d0/ &
                         (1.d0/diag(i)+aux(i)**2-(gg(i)/diag(i))**2/den)
                  if (diag(i).le.0.d0) then
                      if (impres.ge.5) write (io,935) i,diag(i),rmin
                      diag(i)=rmin
                  endif
  431         continue
!$OMP END PARALLEL DO
  935         format (/" >>> n1qn3-WARNING: diagonal element ",i8, &
                       " is negative (",d10.3,"), reset to ",d10.3)
!
              if (impresmax.ge.5) then
                  ps=0.d0
                  do 440 i=1,n
                      ps=ps+diag(i)
  440             continue
                  call mmpi_allreduce_sumreal8scalar(ps,"GRID")
                  ps=ps/ntotal
                  preco=ps
!
                  ps2=0.d0
                  do 441 i=1,n
                      ps2=ps2+(diag(i)-ps)**2
  441             continue
                  call mmpi_allreduce_sumreal8scalar(ps2,"GRID")
                  ps2=dsqrt(ps2/ntotal)
                  if (impres.ge.5) write (io,936) preco,ps2
  936             format (5x,"updated diagonal: average value = ",d10.3, &
                         ", sqrt(variance) = ",d10.3)
              endif
          endif
      endif
!
!     --- tests d'arret
!
      call prosca(n,g,g,ps,izs,rzs,dzs)
      eps1=ps
      eps1=dsqrt(eps1)/gnorm
!
      if (impres.ge.5) write (io,940) eps1
  940 format (/" n1qn3: stopping criterion on g: ",d12.5)
      if (eps1.lt.epsg) then
          mode=1
          goto 1000
      endif
      if (niter.eq.itmax) then
          mode=4
          if (impres.ge.1) write (io,941) niter
  941     format (/" >>> n1qn3 (iteration ",i3, &
                  "): maximal number of iterations")
          goto 1000
      endif
      if (isim.gt.nsim) then
          mode=5
          if (impres.ge.1) write (io,942) niter,isim
  942     format (/" >>> n1qn3 (iteration ",i3,"): ",i6, &
                  " simulations (maximal number reached)")
          goto 1000
      endif
!
!     --- calcul de la nouvelle direction de descente d = - H.g
!
      if (m.eq.0) then
          preco=2.d0*(ff-f)/(eps1*gnorm)**2
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) 
          do 500 i=1,n
              d(i)=-g(i)*preco
  500     continue
!$OMP END PARALLEL DO
      else
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) 
          do 510 i=1,n
              d(i)=-g(i)
  510     continue
!$OMP END PARALLEL DO
          if (inmemo) then
              call ddd (prosca,dtonb,dtcab,n,sscale,m,d,aux,jmin,jmax, &
                        precos,diag,ybar,sbar,izs,rzs,dzs)
          else
              call ddds (prosca,dtonb,dtcab,n,sscale,m,d,aux,jmin,jmax, &
                         precos,diag,ybar,sbar,izs,rzs,dzs)
          endif
      endif
!
!         --- test: la direction d est-elle de descente ?
!             hp0 sera utilise par nlis0
!
      call prosca (n,d,g,hp0,izs,rzs,dzs)
      if (hp0.ge.0.d+0) then
          mode=7
          if (impres.ge.1) write (io,905) niter,hp0
          goto 1000
      endif
      if (impresmax.ge.5) then
          call prosca (n,g,g,ps,izs,rzs,dzs)
          ps=dsqrt(ps)
          call prosca (n,d,d,ps2,izs,rzs,dzs)
          ps2=dsqrt(ps2)
          ps=hp0/ps/ps2
          ps=dmin1(-ps,1.d+0)
          ps=dacos(ps)
          d1=ps
          d1=d1*180.d0/pi
          if (impres.ge.5) write (io,906) sngl(d1)
      endif
!
!---- on poursuit les iterations
!
      goto 100
!
!---- retour
!
 1000 continue
      nsim=isim
      epsg=eps1

      end subroutine n1qn3a

      subroutine nlis0 (n,simul,prosca,xn,fn,fpn,t,tmin,tmax,d,g, &
                        amd,amf,imp,io,logic,nap,napmax,x,izs,rzs,dzs)
      !
      !:Purpose: en sortie logic =
      !
      !          * 0          descente serieuse
      !          * 1          descente bloquee
      !          * 4          nap > napmax
      !          * 5          retour a l'utilisateur
      !          * 6          fonction et gradient pas d'accord
      !          * < 0        contrainte implicite active
      !

!!
!!--- arguments
!!
      external simul,prosca
      integer n,imp,io,logic,nap,napmax,izs(*)
      real rzs(*)
      real(8) xn(n),fn,fpn,t,tmin,tmax,d(n),g(n),amd,amf,x(n),dzs(*)
!
! --- variables locales
!
      logical lfound,lfound2
      integer i,indic,indica,indicd,ierr,ntotal
      real(8) tesf,tesd,tg,fg,fpg,td,ta,fa,fpa,d2,f,fp,ffn,fd, &
       fpd,z,test,barmin,barmul,barmax,barr,gauche,droite,taa,ps
!
 1000 format (/4x,9h nlis0   ,4x,4hfpn=,d10.3,4h d2=,d9.2, &
       7h  tmin=,d9.2,6h tmax=,d9.2)
 1001 format (/4x,6h mlis0,3x,"stop on tmin",8x,  &
         "step",11x,"functions",5x,"derivatives")
 1002 format (4x,6h nlis0,37x,d10.3,2d11.3)
 1003 format (4x,6h nlis0,d14.3,2d11.3)
 1004 format (4x,6h nlis0,37x,d10.3,7h indic=,i3)
 1005 format (4x,6h nlis0,14x,2d18.8,d11.3)
 1006 format (4x,6h nlis0,14x,d18.8,12h      indic=,i3)
 1007 format (/4x,6h mlis0,10x,"tmin forced to tmax")
 1008 format (/4x,6h mlis0,10x,"inconsistent call")
      call rpn_comm_allreduce(n,ntotal,1,"mpi_integer", &
                              "mpi_max","GRID",ierr)
      if (ntotal.gt.0 .and. fpn.lt.0.d0 .and. t.gt.0.d0 &
       .and. tmax.gt.0.d0 .and. amf.gt.0.d0 &
       .and. amd.gt.amf .and. amd.lt.1.d0) go to 5
      logic=6
      go to 999
    5 tesf=amf*fpn
      tesd=amd*fpn
      barmin=0.01d0
      barmul=5.d0
      barmax=0.3d0
      barr=barmin
      td=0.d0
      tg=0.d0
      fg=fn
      fpg=fpn
      ta=0.d0
      fa=fn
      fpa=fpn
      call prosca (n,d,d,ps,izs,rzs,dzs)
      d2=ps
!
!               elimination d'un t initial ridiculement petit
!
      if (t.gt.tmin) go to 20
      t=tmin
      if (t.le.tmax) go to 20
      if (imp.gt.0) write (io,1007)
      tmin=tmax
   20 if (fn+t*fpn.lt.fn+0.9d0*t*fpn) go to 30
      t=2.d0*t
      go to 20
   30 indica=1
      logic=0
      if (t.gt.tmax) then
          t=tmax
          logic=1
      endif
      if (imp.ge.4) write (io,1000) fpn,d2,tmin,tmax
!
!     --- nouveau x
!
      do 50 i=1,n
          x(i)=xn(i)+t*d(i)
   50 continue
!
! --- boucle
!
  100 nap=nap+1
      if(nap.gt.napmax) then
          logic=4
          fn=fg
          do 120 i=1,n
              xn(i)=xn(i)+tg*d(i)
  120     continue
          go to 999
      endif
      indic=4
!
!     --- appel simulateur
!
      call simul(indic,n,x,f,g,izs,rzs,dzs)
      if(indic.eq.0) then
!
!         --- arret demande par l'utilisateur
!
          logic=5
          fn=f
          do 170 i=1,n
              xn(i)=x(i)
  170     continue
          go to 999
      endif
      if(indic.lt.0) then
!
!         --- les calculs n'ont pas pu etre effectues par le simulateur
!
          td=t
          indicd=indic
          logic=0
          if (imp.ge.4) write (io,1004) t,indic
          t=tg+0.1d0*(td-tg)
          go to 905
      endif
!
!     --- les tests elementaires sont faits, on y va
!
      call prosca (n,d,g,ps,izs,rzs,dzs)
      fp=ps
!
!     --- premier test de Wolfe
!
      ffn=f-fn
      if(ffn.gt.t*tesf) then
          td=t
          fd=f
          fpd=fp
          indicd=indic
          logic=0
          if(imp.ge.4) write (io,1002) t,ffn,fp
          go to 500
      endif
!
!     --- test 1 ok, donc deuxieme test de Wolfe
!
      if(imp.ge.4) write (io,1003) t,ffn,fp
      if(fp.gt.tesd) then
          logic=0
          go to 320
      endif
      if (logic.eq.0) go to 350
!
!     --- test 2 ok, donc pas serieux, on sort
!
  320 fn=f
      do 330 i=1,n
          xn(i)=x(i)
  330 continue
      go to 999
!
!
!
  350 tg=t
      fg=f
      fpg=fp
      if(td.ne.0.d0) go to 500
!
!              extrapolation
!
      taa=t
      gauche=(1.d0+barmin)*t
      droite=10.d0*t
      call dcube (t,f,fp,ta,fa,fpa,gauche,droite)
      ta=taa
      if(t.lt.tmax) go to 900
      logic=1
      t=tmax
      go to 900
!
!              interpolation
!
  500 if(indica.le.0) then
          ta=t
          t=0.9d0*tg+0.1d0*td
          go to 900
      endif
      test=barr*(td-tg)
      gauche=tg+test
      droite=td-test
      taa=t
      call dcube (t,f,fp,ta,fa,fpa,gauche,droite)
      ta=taa
      if (t.gt.gauche .and. t.lt.droite) then
          barr=dmax1(barmin,barr/barmul)
!         barr=barmin
        else
          barr=dmin1(barmul*barr,barmax)
      endif
!
! --- fin de boucle
!     - t peut etre bloque sur tmax
!       (venant de l'extrapolation avec logic=1)
!
  900 fa=f
      fpa=fp
  905 indica=indic
!
! --- faut-il continuer ?
!
      if (td.eq.0.d0) go to 950
      if (td-tg.lt.tmin) go to 920
!
!     --- limite de precision machine (arret de secours) ?
!
      lfound=.false.
      do i=1,n
          z=xn(i)+t*d(i)
          if (z.ne.xn(i).and.z.ne.x(i)) then
            lfound=.true.
            exit
          endif
      enddo
      call rpn_comm_allreduce(lfound,lfound2,1,"mpi_logical", &
           "mpi_lor","GRID",ierr)
      if(lfound2) go to 950
!
! --- arret sur dxmin ou de secours
!
  920 logic=6
!
!     si indicd<0, derniers calculs non faits par simul
!
      if (indicd.lt.0) logic=indicd
!
!     si tg=0, xn = xn_depart,
!     sinon on prend xn=x_gauche qui fait decroitre f
!
      if (tg.eq.0.d0) go to 940
      fn=fg
      do 930 i=1,n
  930 xn(i)=xn(i)+tg*d(i)
  940 if (imp.le.0) go to 999
      write (io,1001)
      write (io,1005) tg,fg,fpg
      if (logic.eq.6) write (io,1005) td,fd,fpd
      if (logic.eq.7) write (io,1006) td,indicd
      go to 999
!
!               recopiage de x et boucle
!
  950 do 960 i=1,n
  960 x(i)=xn(i)+t*d(i)
      go to 100
  999 continue
      end subroutine nlis0

      end module quasinewton_mod
