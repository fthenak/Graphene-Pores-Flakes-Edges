      program pore_flake_edge_structure
!     Given: 
!
!     (a) the number of 2-fold coordinated atoms (n2) in the pore,
!         flake or edge boundary, 
!     (b) the sequence of 2-2 paths (seq(n2)) and
!     (c) (for pores only) the width of the supercell (d), 
!
!     the program finds if the given sequence represents a pore,
!     flake or edge. 
!     A. If the sequence represents a pore, the program finds the
!        atomic coordinates of the porous graphene supercell
!        corresponding to that sequence and the lattice vectors of
!        the supercell, and stores them in the file 'structure.xyz'
!        in xyz format. In the comment (second) line, the lattice
!        vector coordinates are stored. Moreover, the program finds
!        the minimum image, the number of vacancies and the number
!        of hexagons in the pore vacuum space.
!     B. If the sequence represents a flake, the program finds the
!        atomic coordinates of the flake atoms corresponding to that
!        sequence and stores them in the file 'structure.xyz' in xyz
!        format. Moreover, the program finds the minimum image, the
!        number of atoms and the number of hexagons of the flake.
!     C. If the sequence represents a (periodic) graphene edge
!        (semi-infinite graphene), the program finds the atomic
!        coordinates of the boundary and the vector representing the
!        period and stores them in the file 'structure.xyz' in xyz
!        format. In the comment (second) line, the coordinates of the
!        vector representing the period of the edge (i.e. latice
!        vector of the one-dimensional ribbon) is stored. Moreover,
!        the program finds the minimum image and the number of the
!        edge atoms in the period.
!        
!
!     Zacharias G. Fthenakis 
!     
!     Istituto Nanoscienze-CNR, Piazza San Silvestro 12, 56127 Pisa, Italy
!     and
!     Theoretical and Physical Chemistry Institute, National Hellenic
!     Research Foundation, GR-11635, Athens, Greece

      implicit none
      integer n2,n21,L0,iswitch,n_vac,n_h,i,hx_min,hx_max
      integer h_vac(1000,2),sum_h,d,system,L
      integer h_flake(1000,2),natom,n_first
      integer, allocatable :: m(:),h(:,:),seq(:),first(:)
      character (LEN=100) FMT

      print *, 'Give the number of 2-fold coordinated atoms'
      print *, 'in the pore boundary'
      read(5,*) n2
      n21=n2+1
      write(FMT,*) n2

      allocate(seq(n2),h(n21,2),m(n21),first(n2))

      print *, 'System: 1 for pore, 2 for flake, 3 for edge'
      read(5,*) system

      print *, 'Give the sequence of lengths of 2-2 paths'
      read(5,*) (seq(i),i=1,n2)

      do i=1,n2
         if (seq(i)<1.or.seq(i)>4) then
         print *, 'Wrong 2-2 path lenght(s) in the pore sequence'
         stop
         endif
      enddo

      if (system==1) then
         print *, 'Give the width of the supercell'
         read(5,*) d

         L0=2*n2+6
      else if (system==2) then
         L0=2*n2-6
      else if (system==3) then
         L0=2*n2
      else
         print *, 'Wrong system'
         stop
      endif

      L=0
      do i=1,n2
         L=L+seq(i)
      enddo
      if (L/=L0) then
         call negative(system)
         stop
      endif

      if (system==3) then
         call multipliers(n2,first,n_first)
      endif

      call test1(n2,seq,iswitch)
      if (iswitch==0) then
         call negative(system)
         stop
      else
         call test2(system,L0,n2,n21,h,m,seq,iswitch)
!          do i=1,n21
!             print *,i,h(i,1),h(i,2)
!          enddo
         if (iswitch==0) then
            call negative(system)
            stop
         else
            call touch(n2,n21,h,m,seq,iswitch)
            if (iswitch==0) then
               call negative(system)
               stop
            else
               write(6,"("//adjustl(FMT)//"I1)") seq(1:n2)
 107           format('Minimum image = ', 100i1)
               if (system==1) then
                  call vacancies(n2,n21,h,m,seq,sum_h,hx_min,hx_max,n_vac,h_vac)
                  call structure(system,n2,n21,d,h,m,seq,hx_min,hx_max,sum_h,h_vac)
                  n_h=(n_vac+n2+2)/2
                  call circle(n2,n21,m,h,seq,sum_h,h_vac)
                  print *, 'Number of vacancies = ', n_vac
                  print *, 'Number of hexagons in the vacuum area = ',n_h,sum_h
               else if (system==2) then
                  call flake(n2,n21,h,m,seq,sum_h,hx_min,hx_max,natom,h_flake)
                  call structure(system,n2,n21,0,h,m,seq,hx_min,hx_max,sum_h,h_flake)
                  n_h=(natom-n2+2)/2
                  print *, 'Number of atoms in the flake = ', natom
                  print *, 'Number of hexagons belonging to the flake = ',n_h,sum_h
               else if (system==3) then
                  call periodicity(n2,seq,first,n_first,iswitch)
                  if (iswitch==0.or.all(seq(:)==2)) then
                     print *, 'The considered as period is a multiple of &
                          & a smaller period. Use the minimum period'
                     stop
                  else
                     call touch_e(n2,n21,h,m,seq,iswitch)
                     if (iswitch==0) then
                        call negative(system)
                        stop
                     else
                        call structure_edge(n2,n21,h,m,seq)
                     endif
                  endif
               endif
               print *, 'Total lenght of the circuit (or path) = ', L0
            endif
         endif
      endif
      end program 

      subroutine definitions(a0,a,b,v)
      implicit none
      real*8 a0,a(2),b(2),v(0:5,2)

      a0=1.42d0
      a(1)=1.5d0*a0
      a(2)=dsqrt(3.d0)/2.d0*a0
      b(1)=0.d0
      b(2)=dsqrt(3.d0)*a0

      v(0,1)=a0
      v(0,2)=0.d0
      v(1,1)=a0/2.d0
      v(1,2)=dsqrt(3.d0)*a0/2.d0
      v(2,1)=-v(1,1)
      v(2,2)=v(1,2)
      v(3,1)=-v(0,1)
      v(3,2)=-v(0,2)
      v(4,1)=-v(1,1)
      v(4,2)=-v(1,2)
      v(5,1)=-v(2,1)
      v(5,2)=-v(2,2)
      return
      end subroutine definitions

      subroutine circle(n2,n21,m,h,seq,sum_h,h_vac)
      implicit none
      integer i,n2,n21,m(n21),h(n21,2),seq(n2),mm,j,sum_h
      integer i1,i2,i3,iswitch,h_vac(1000,2)
      real*8 a0,r2(n2,2),a(2),b(2),v(0:5,2)
      real*8 r21(2),r22(2),r23(2),d21,d22,d23,A12,A13,D,Dx,Dy
      real*8 c(2),radi

      call definitions (a0,a,b,v)

!     r2 = positions of 2-fold coordinated atoms
      write(7,*) n2
      write(7,*) ' '
      do i=1,n2
         r2(i,:)=h(i,1)*a(:)+h(i,2)*b(:)+v(mod(m(i)+1,6),:)
         write(7,*) "C ",r2(i,1),r2(i,2)," 0.0"
      enddo

      print *, 'CIRCUMSCRIBED CIRCLES in the pore'
      print *, 'diameter in a_0 units | center of the circle'
      do i1=1,n2-2
         r21(:)=r2(i1,:)
         do i2=i1+1,n2-1
            r22(:)=r2(i2,:)
            do i3=i2+1,n2
               r23(:)=r2(i3,:)
               call radius(r21,r22,r23,c,radi)
               call in_polygon(c,radi,sum_h,h_vac,iswitch,i1,i2,i3)

               if (iswitch==1) then
                  call in_radius(n2,r2,c,radi,iswitch)
                  if (iswitch==1) then
                     print *, 2.d0*radi/a0,c(1),c(2)
                  endif
               endif
            enddo
         enddo
      enddo
      return
      end subroutine circle

      subroutine in_radius(n2,r2,c,radi,iswitch)
      implicit none
      integer i,n2,iswitch
      real*8 c(2),radi,a0,a(2),b(2),hx,hy
      real*8 r2(n2,2),d2,radi2,tol

      tol=0.0001d0

      radi2=radi*radi
      do i=1,n2
         d2=(r2(i,1)-c(1))*(r2(i,1)-c(1))+(r2(i,2)-c(2))*(r2(i,2)-c(2))
         if (d2+tol<radi2) then
            iswitch=0
            return
         endif
      enddo
      iswitch=1
      return
      end subroutine in_radius


      subroutine in_polygon(c,radi,n_vac,h_vac,iswitch,i1,i2,i3)   
!     Checks if the center c(2) is in the vacumm area of the pore
      implicit none
      integer i,n_vac,iswitch,h_vac(1000,2),i1,i2,i3
      real*8 c(2),radi,a0,a(2),b(2),hx,hy
      real*8 tol,xx1,xx2,yy,v(0:5,2)

      tol=0.0001d0

      call definitions (a0,a,b,v)

      do i=1,n_vac
         hx=h_vac(i,1)*a(1)+h_vac(i,2)*b(1)
         hy=h_vac(i,1)*a(2)+h_vac(i,2)*b(2)
         yy=dabs(c(2)-hy)/dsqrt(3.d0)
         xx1=-(c(1)-hx-a0-tol)
         xx2=c(1)-hx+a0+tol
         if (i1==4.and.i2==5.and.i3==6) then
!            print *, 'in', yy,xx1,xx2,a0/2.d0
         endif
         if (yy<=a0/2.d0+tol.and.yy<=xx1.and.yy<=xx2) then
             iswitch=1
             return
         endif
      enddo
      iswitch=0
      return
      end subroutine in_polygon


      subroutine radius(r21,r22,r23,c,radi)
      implicit none
      real*8 r21(2),r22(2),r23(2),d21,d22,d23,A12,A13,D,Dx,Dy
      real*8 c(2),radi

      d22=r22(1)*r22(1)+r22(2)*r22(2)
      d21=r21(1)*r21(1)+r21(2)*r21(2)
      d23=r23(1)*r23(1)+r23(2)*r23(2)
      A12=(d22-d21)/2.d0
      A13=(d23-d21)/2.d0
      D=(r22(1)-r21(1))*(r23(2)-r21(2))-(r23(1)-r21(1))*(r22(2)-r21(2))
      Dx=A12*(r23(2)-r21(2))-A13*(r22(2)-r21(2))
      Dy=A13*(r22(1)-r21(1))-A12*(r23(1)-r21(1))
      c(1)=Dx/D
      c(2)=Dy/D
      radi=dsqrt((r21(1)-c(1))*(r21(1)-c(1))+(r21(2)-c(2))*(r21(2)-c(2)))
      return
      end subroutine radius
         

      subroutine multipliers(n2,first,n_first)
!     It finds the multipliers of n2, to check if the period of an edge
!     is smaller than the given one.
!
      implicit none
      integer n2,i,first(n2),n_first,sq_n2

      sq_n2=int(sqrt(float(n2)))+1
      n_first=0
      do i=2,sq_n2
         if ((n2/i)*i==n2) then
            n_first=n_first+1
            first(n_first)=i
            print *,'multip', n_first, first(n_first)
            n_first=n_first+1
            first(n_first)=n2/i
            print *,'multip', n_first, first(n_first)
         endif
      enddo
      return
      end subroutine multipliers


      subroutine  periodicity(n2,seq,first,n_first,iswitch)
!     Check if the given period is a multiple of the minimum period
      implicit none
      integer seq(n2),n2,i,j,k,n_first,sq_n2,iswitch,first(n2)
      integer first1,first2,j1,j2,j3,j4


      iswitch=1
      do i=1,n_first,2
         first1=first(i)
         first2=first(i+1)
         j1=1
         j2=first2
         do j=1,first1-1 ! first2 = posible period 
            j3=j1+first2
            j4=j2+first2
            if (all(seq(j1:j2)==seq(j3:j4))) then
               iswitch=0
            else
               iswitch=1
               exit
            endif
           j1=j3
            j2=j4
         enddo
         if (iswitch==0) return
         do j=1,first2-1 ! first1 = posible period 
            j3=j1+first1
            j4=j2+first1
            if (all(seq(j1:j2)==seq(j3:j4))) then
               iswitch=0
            else
               iswitch=1
               exit
            endif
            j1=j3
            j2=j4
         enddo
         if (iswitch==0) return
      enddo
      return
      end subroutine  periodicity


      subroutine negative(system)
      implicit none
      integer system
         if (system==1) then
            print *, 'The sequence does not represent a pore - test 1'
            stop
         else if (system==2) then
            print *, 'The sequence does not represent a flake - test 1'
            stop
         else if (system==3) then
            print *, 'The sequence does not represent an edge - test 1'
            stop
         endif
      return
      end subroutine negative

      subroutine test2(system,L0,n2,n21,h,m,seq,iswitch)
!     a. For pores and flakes it checks  whether or not the minimum image 
!     sequence represents a circuit (i.e. a closed path). 
!     b. It checks wether or not the pore, flake or edge sequence represents
!     a non crossed circuit, according to the first criterion.
!     c. For edges it checks if the doubled edge (two periods) forms a
!     crossed path.

      implicit none
      integer n2,n21,seq(n2),iswitch,n_vac,L0,nn,i,j,k,mm
      integer h(n21,2),m(n21),cos1,cos2,system
      integer, allocatable :: h2(:,:)

      allocate (h2(n21,2))

      iswitch=1
      nn=n2-3

      m(1)=0
      h(1,1)=0
      h(1,2)=0
      i=1
      do i=2,n2
         m(i)=m(i-1)+seq(i)-2 ! m value before the last one in the 2-2 path
         if (m(i)==6) m(i)=0
         if (m(i)==7) m(i)=1
         if (m(i)==-1) m(i)=5
         h(i,1)=h(i-1,1)+cos1(m(i-1))
         h(i,2)=h(i-1,2)+cos2(m(i-1))
      enddo
      m(n21)=m(n2)+seq(1)-2
      
      if (m(n21)==6) m(n21)=0
      if (m(n21)==7) m(n21)=1
      if (m(n21)==-1) m(n21)=5
      h(n2+1,1)=h(n2,1)+cos1(m(n2))
      h(n2+1,2)=h(n2,2)+cos2(m(n2))
      
      if (system==1.or.system==2) then
         if (m(n21)/=m(1).or.h(n21,1)/=h(1,1).or.h(n21,2)/=h(1,2)) then
            iswitch=0
            print *, 'not closed',m(i),m(1),h(i,1),h(1,1),h(i,2),h(1,2)
            return
         endif
      endif

      do i=1,n2-1
         do j=i+1,n2
            if (h(i,1)==h(j,1).and.h(i,2)==h(j,2)) then
               if (seq(i)/=1.or.seq(j)/=1) then
                  iswitch=0 
                  print *, 'same hexagon - 1'
                  return
               else
                  mm=m(i)-m(j)
                  if (mm/=3.and.mm/=-3) then
                     iswitch=0
                     print *, 'same hexagon - 2'
                     return
                  endif
               endif
            endif
         enddo
      enddo

      if (system==3) then
!     Checks if the doubled edge (two periods) forms a crossed path
!         print *, h(n21,:)
         do i=1,n21
            h2(i,:)=h(i,:)+h(n21,:)
!            print *, h2(i,:),h(i,:)
         enddo
         do i=1,n2
            do j=1,n2
               if (h(i,1)==h2(j,1).and.h(i,2)==h2(j,2)) then
                  if (seq(i)/=1.or.seq(j)/=1) then
                     iswitch=0 
!                    print *, 'error 1'
                     return
                  else
                     mm=m(i)-m(j)
                     if (mm/=3.and.mm/=-3) then
                        iswitch=0
!                    print *, 'error 2'
                        return
                     endif
                  endif
               endif
            enddo
         enddo
      endif

      return

      end subroutine test2

      subroutine touch(n2,n21,h,m,seq,iswitch)
!     It checks if a circuit or a path is a crossed circuit or path
!     according to the second criterion.
      implicit none
      integer h(n21,2),m(n21),seq(n2),iswitch
      integer i,j,k,n2,n21,mij1,mij2,d1,d2

      iswitch=1
      do i=1,n2
         do j=i+1,n21
            if (h(j,1)-h(i,1)==1.and.h(j,2)-h(i,2)==0) then
               mij1=0
               d1=m(i)-mij1
               if (d1<0) d1=d1+6
               mij2=3
               d2=m(j)-mij2
               if (d2<0) d2=d2+6
               if (d1<=seq(i)-1.and.d2<=seq(j)-1) then
                  iswitch=0
                  return
               endif
            elseif (h(j,1)-h(i,1)==0.and.h(j,2)-h(i,2)==1) then
               mij1=1
               d1=m(i)-mij1
               if (d1<0) d1=d1+6
               mij2=4
               d2=m(j)-mij2
               if (d2<0) d2=d2+6
               if (d1<=seq(i)-1.and.d2<=seq(j)-1) then
                  iswitch=0
                  return
               endif
            elseif (h(j,1)-h(i,1)==-1.and.h(j,2)-h(i,2)==1) then
               mij1=2
               d1=m(i)-mij1
               if (d1<0) d1=d1+6
               mij2=5
               d2=m(j)-mij2
               if (d2<0) d2=d2+6
               if (d1<=seq(i)-1.and.d2<=seq(j)-1) then
                  iswitch=0
                  return
               endif
            elseif (h(j,1)-h(i,1)==-1.and.h(j,2)-h(i,2)==0) then
               mij1=3
               d1=m(i)-mij1
               if (d1<0) d1=d1+6
               mij2=0
               d2=m(j)-mij2
               if (d2<0) d2=d2+6
               if (d1<=seq(i)-1.and.d2<=seq(j)-1) then
                  iswitch=0
                  return
               endif
            elseif (h(j,1)-h(i,1)==0.and.h(j,2)-h(i,2)==-1) then
               mij1=4
               d1=m(i)-mij1
               if (d1<0) d1=d1+6
               mij2=1
               d2=m(j)-mij2
               if (d2<0) d2=d2+6
               if (d1<=seq(i)-1.and.d2<=seq(j)-1) then
                  iswitch=0
                  return
               endif
            elseif (h(j,1)-h(i,1)==1.and.h(j,2)-h(i,2)==-1) then
               mij1=5
               d1=m(i)-mij1
               if (d1<0) d1=d1+6
               mij2=2
               d2=m(j)-mij2
               if (d2<0) d2=d2+6
               if (d1<=seq(i)-1.and.d2<=seq(j)-1) then
                  iswitch=0
                  return
               endif
            endif
         enddo
      enddo
      return
      end subroutine touch

      subroutine touch_e(n2,n21,h,m,seq,iswitch)
!     It checks if the sequences corresponding to the current and the
!     next period of an edge are crossed with each other
!     according to the second criterion.
      implicit none
      integer h(n21,2),m(n21),seq(n2),iswitch
      integer i,j,k,n2,n21,mij1,mij2,d1,d2
      integer, allocatable :: h2(:,:)

      allocate (h2(n21,2))


      do i=1,n2
         h2(i,:)=h(i,:)+h(n21,:)
      enddo
      h2(n21,:)=h(n2,:)+h(n21,:)

      iswitch=1
      do i=1,n21
         do j=1,n21
            if (h2(j,1)-h(i,1)==1.and.h2(j,2)-h(i,2)==0) then
               mij1=0
               d1=m(i)-mij1
               if (d1<0) d1=d1+6
               mij2=3
               d2=m(j)-mij2
               if (d2<0) d2=d2+6
               if (d1<=seq(i)-1.and.d2<=seq(j)-1) then
                  iswitch=0
                  return
               endif
            elseif (h2(j,1)-h(i,1)==0.and.h2(j,2)-h(i,2)==1) then
               mij1=1
               d1=m(i)-mij1
               if (d1<0) d1=d1+6
               mij2=4
               d2=m(j)-mij2
               if (d2<0) d2=d2+6
               if (d1<=seq(i)-1.and.d2<=seq(j)-1) then
                  iswitch=0
                  return
               endif
            elseif (h2(j,1)-h(i,1)==-1.and.h2(j,2)-h(i,2)==1) then
               mij1=2
               d1=m(i)-mij1
               if (d1<0) d1=d1+6
               mij2=5
               d2=m(j)-mij2
               if (d2<0) d2=d2+6
               if (d1<=seq(i)-1.and.d2<=seq(j)-1) then
                  iswitch=0
                  return
               endif
            elseif (h2(j,1)-h(i,1)==-1.and.h2(j,2)-h(i,2)==0) then
               mij1=3
               d1=m(i)-mij1
               if (d1<0) d1=d1+6
               mij2=0
               d2=m(j)-mij2
               if (d2<0) d2=d2+6
               if (d1<=seq(i)-1.and.d2<=seq(j)-1) then
                  iswitch=0
                  return
               endif
            elseif (h2(j,1)-h(i,1)==0.and.h2(j,2)-h(i,2)==-1) then
               mij1=4
               d1=m(i)-mij1
               if (d1<0) d1=d1+6
               mij2=1
               d2=m(j)-mij2
               if (d2<0) d2=d2+6
               if (d1<=seq(i)-1.and.d2<=seq(j)-1) then
                  iswitch=0
                  return
               endif
            elseif (h2(j,1)-h(i,1)==1.and.h2(j,2)-h(i,2)==-1) then
               mij1=5
               d1=m(i)-mij1
               if (d1<0) d1=d1+6
               mij2=2
               d2=m(j)-mij2
               if (d2<0) d2=d2+6
               if (d1<=seq(i)-1.and.d2<=seq(j)-1) then
                  iswitch=0
                  return
               endif
            endif
         enddo
      enddo
      return
      end subroutine touch_e


      subroutine vacancies(n2,n21,h,m,seq,sum_h,hx_min,hx_max,n_vac,h_vac)
!     It finds the hexagons in the vacuum pore area and stores
!     their coordinates in h_vac. Moreover, it finds the number of
!     vacancies (n_vac).
!     input: h, m, seq - output: h_vac, n_vac

      implicit none
      integer n2,n21,i,j,k,ij,ii,kk,ll,h(n21,2),m(n21),seq(n2)
      integer hx_min,hx_max,sum_h,n_vac,h_vac(1000,2),up,down
      integer mmm
      integer, allocatable :: hh(:),sseq(:),mm(:),signe(:)

      allocate(hh(n2),sseq(n2),mm(n2),signe(n2))

      hx_min=h(1,1)
      hx_max=h(1,1)
      do i=2,n2
         if (h(i,1)>hx_max) hx_max=h(i,1) 
         if (h(i,1)<hx_min) hx_min=h(i,1) 
      enddo
      sum_h=0

      do i=hx_min,hx_max
         k=0
         do j=1,n2
            if (h(j,1)==i) then
               k=k+1
               hh(k)=h(j,2)
               mm(k)=m(j)
               sseq(k)=seq(j)
            endif
         enddo
         if (k>1) then
            call sort(n2,hh,mm,sseq,k)
         endif
! correction start
         do kk=1,k-1
            if (hh(kk)==hh(kk+1)) then
               if (mm(kk)==4) then
                  mmm=mm(kk)
                  mm(kk)=mm(kk+1)
                  mm(kk+1)=mmm
               endif
            endif
         enddo
! correction end
         do kk=1,k
            if (sseq(kk)==4.and.(mm(kk)==1.or.mm(kk)==4)) then
               signe(kk)=0
            elseif ((sseq(kk)==4.and.(mm(kk)==3.or.mm(kk)==2)).or. &
                 &  (sseq(kk)==3.and.(mm(kk)==1.or.mm(kk)==2.or.mm(kk)==3)).or. &           
                 &  (sseq(kk)==2.and.(mm(kk)==1.or.mm(kk)==2)).or. &
                 &  (sseq(kk)==1.and.mm(kk)==1)) then
                 signe(kk)=-1 ! <-----
            elseif ((sseq(kk)==4.and.(mm(kk)==5.or.mm(kk)==0)).or. &
                 &  (sseq(kk)==3.and.(mm(kk)==4.or.mm(kk)==5.or.mm(kk)==0)).or. &           
                 &  (sseq(kk)==2.and.(mm(kk)==4.or.mm(kk)==5)).or. &
                 &  (sseq(kk)==1.and.mm(kk)==4)) then
                 signe(kk)=1 ! ----->
            else
                 signe(kk)=2
            endif
         enddo

         do kk=1,k
            if (signe(kk)==-1) up=hh(kk) 
            if (signe(kk)==1) then
               down=hh(kk)
               do ij=up,down,-1
                   sum_h=sum_h+1
                  h_vac(sum_h,1)=i
                  h_vac(sum_h,2)=ij
               enddo
            endif
            if (signe(kk)==0) then
               sum_h=sum_h+1
               h_vac(sum_h,1)=i
               h_vac(sum_h,2)=hh(kk)
            endif
            if (sum_h>1000) then
               print *,'WARNING: Please increase the dimension of h_vac'
               print *, sum_h
               stop
            endif
         enddo
      enddo
      n_vac=2*sum_h-n2-2
      return
      end subroutine vacancies

      subroutine flake(n2,n21,h,m,seq,sum_h,hx_min,hx_max,natom,h_flake)
!     (a) It finds the hexagons in the flake and stores
!     their coordinates in h_flake.
!     (b) It finds the number of atoms (natom).
!     input: h, m, seq - output: h_flake, natom

      implicit none
      integer n2_max,n21_max
      integer n2,n21,i,j,k,ij,ii,kk,ll,h(n21,2),m(n21)
      integer hx_min,hx_max,sum_h,natom,h_flake(1000,2),up,down
      integer seq(n2)
      integer, allocatable :: hh(:),sseq(:),mm(:),signe(:)

      allocate(hh(n2),sseq(n2),mm(n2),signe(n2))

      hx_min=h(1,1)
      hx_max=h(1,1)
      print *, 1,h(1,1),h(1,2)
      do i=2,n2
         if (h(i,1)>hx_max) hx_max=h(i,1)
         if (h(i,1)<hx_min) hx_min=h(i,1)
         print *, i,h(i,1),h(i,2)
      enddo
      hx_min=hx_min+1
      hx_max=hx_max-1
      sum_h=0
      do i=hx_min,hx_max
         k=0
         do j=1,n2
            if (h(j,1)==i) then
               k=k+1
               hh(k)=h(j,2)
               mm(k)=m(j)
               sseq(k)=seq(j)
            endif
         enddo
         print *, 'col ', i
         if (k>1) then
            call sort(n2,hh,mm,sseq,k)
         endif
         do kk=1,k
            if (sseq(kk)==4.and.(mm(kk)==1.or.mm(kk)==4)) then
               signe(kk)=0
            elseif ((sseq(kk)==4.and.(mm(kk)==3.or.mm(kk)==2)).or. &
                 &  (sseq(kk)==3.and.(mm(kk)==1.or.mm(kk)==2.or.mm(kk)==3)).or.  &
                 &  (sseq(kk)==2.and.(mm(kk)==1.or.mm(kk)==2)).or. &
                 &  (sseq(kk)==1.and.mm(kk)==1)) then
                 signe(kk)=-1 ! <-----
            elseif ((sseq(kk)==4.and.(mm(kk)==5.or.mm(kk)==0)).or. &
                 &  (sseq(kk)==3.and.(mm(kk)==4.or.mm(kk)==5.or.mm(kk)==0)).or.  &
                 &  (sseq(kk)==2.and.(mm(kk)==4.or.mm(kk)==5)).or. &
                 &  (sseq(kk)==1.and.mm(kk)==4)) then
                 signe(kk)=1 ! ----->
            else
                 signe(kk)=2
            endif
            print *, 'hex ', hh(kk),'sign= ',signe(kk)
         enddo
         do kk=1,k
            if (signe(kk)==1) up=hh(kk)
            if (signe(kk)==-1) then
               down=hh(kk)
               do ij=up-1,down+1,-1
                  sum_h=sum_h+1
                  h_flake(sum_h,1)=i
                  h_flake(sum_h,2)=ij
                  print *, sum_h,i,ij
               enddo
            endif
         enddo
      enddo
      natom=2*sum_h+n2-2
         print *, 'natom = ', natom
      return
      end subroutine flake

      subroutine structure_edge(n2,n21,h,m,seq)
!     It finds the atomic coordinates of the boundary atoms in an edge.
      implicit none
      integer m1,mk
      real*8 a(2),b(2),a0,v(0:5,2),rr(2),r(1000,2),dr(2),T(2)
      integer i,j,k1,k2,ii,kk,n2,n21,h(n21,2),m(n21),seq(n2)

      call definitions(a0,a,b,v)

      ii=0
      T=0 ! period
!            print *, 'ii,i,j,m1, seq(k2),k2,mk'
      do k2=1,n2 ! for 2-2 hexagon
         i=h(k2,1)
         j=h(k2,2)
         m1=mod(m(k2)-seq(k2)+6,6)
         do mk=1,seq(k2)
            m1=mod(m1+1,6)
!            print *, ii,i,j,m1, seq(k2),k2,mk
            ii=ii+1
            r(ii,:)=dfloat(i)*a(:)+dfloat(j)*b(:)+v(m1,:)
         enddo
      enddo

      open (unit=1,file='structure.xyz',status='unknown')
      write(1,*) ii
      write(1,*) 'L = (',(h(n21,1)-h(1,1))*a(:)+(h(n21,2)-h(1,2))*b(:),')'
      do i=1,ii
         write(1,111) 'C ',r(i,1),r(i,2), ' 0.0'
      enddo
 111  format(a2,2x,2(f15.6,2x),a4)
      return
      end subroutine structure_edge
      
         

      subroutine structure(system,n2,n21,d,h,m,seq,hx_min,hx_max,sum_h,hh)
!     Structure of the flake or the supercell of the pore
      implicit none
      real*8 a(2),b(2),a0,v(0:5,2),rr(2),r(1000,2),dr(2),tol
      real*8 aa(2),bb(2) ! supercell lattice vectors
      integer hx_min,hx_max,sum_h,n_vac,hh(1000,2)
      integer hy_min,hy_max,hx_min2,hx_max2,hy_min2,hy_max2
      integer i,j,k1,k2,ii,kk,n2,n21,h(n21,2),m(n21),seq(n2),d
      integer iswitch,k,m1,mk,system
!     For pores : hh(i,2) = hexagons in the vacuum area
!     For flakes: hh(i,2) = hexagons of the flake

      tol=1.d-8
      hy_min=h(1,2)
      hy_max=h(1,2)
      do i=2,n2
         if (h(i,2)>hy_max) hy_max=h(i,2) 
         if (h(i,2)<hy_min) hy_min=h(i,2) 
      enddo
      hx_min2=hx_min-d
      hy_min2=hy_min-d
      hx_max2=hx_max+d
      hy_max2=hy_max+d

      call definitions(a0,a,b,v)

      ii=0
      print *, 'x min, max', hx_min2, hx_max2
      print *, 'y min, max', hy_min2, hy_max2

!     For hexagons inside the pore or flake boundary
      do i=hx_min2,hx_max2
         do j=hy_min2,hy_max2
            iswitch=1
            do k1=1,sum_h 
! For pores, it checks if the hexagon (i,j) belongs to the vacuum area.
! For flakes, it checks if the hexagon (i,j) belongs to the flake.

               if (system==1.and.i==hh(k1,1).and.j==hh(k1,2)) then
                     iswitch=0
                     exit
               else if (system==2.and.i==hh(k1,1).and.j==hh(k1,2)) then
                     iswitch=1
                     exit
               else if (system==2) then
                  iswitch=0
               endif
            enddo
            if (iswitch==1) then
               rr(1)=dfloat(i)*a(1)+dfloat(j)*b(1)
               rr(2)=dfloat(i)*a(2)+dfloat(j)*b(2)
               ii=ii+1
               r(ii,1)=rr(1)+v(4,1)
               r(ii,2)=rr(2)+v(4,2)
               ii=ii+1
               r(ii,1)=rr(1)+v(5,1)
               r(ii,2)=rr(2)+v(5,2)
            endif
         enddo
      enddo

!     For 2-2 hexagons
      do k2=1,n2 ! for 2-2 hexagon
         i=h(k2,1)
         j=h(k2,2)
         m1=mod(m(k2)-seq(k2)+6,6)
         do mk=1,seq(k2)+1
            m1=mod(m1+1,6)

            if (m1==4.or.m1==5) then
               rr(1)=dfloat(i)*a(1)+dfloat(j)*b(1)
               rr(2)=dfloat(i)*a(2)+dfloat(j)*b(2)
               ii=ii+1
               r(ii,1)=rr(1)+v(m1,1)
               r(ii,2)=rr(2)+v(m1,2)
            endif 
         enddo
      enddo

      if (system==1) then
         aa(1)=dfloat(hx_max2-hx_min2+1)*a(1)
         aa(2)=dfloat(hx_max2-hx_min2+1)*a(2)
         bb(1)=0.d0
         bb(2)=dfloat(hy_max2-hy_min2+1)*b(2)
      endif
      open (unit=1,file='structure.xyz',status='unknown')
      write(1,*) ii
      if (system==1) then
         write(1,112) 'a=(',aa(1),aa(2),'0.0), b=(',bb(1),bb(2),' 0.0)'
 112     format(a3,2(f15.6,2x),a9,2(f15.6,2x),a5)
      else
         write(1,*) ' '
      endif
      do i=1,ii
         write(1,111) 'C ',r(i,1),r(i,2), ' 0.0'
      enddo
 111  format(a2,2x,2(f15.6,2x),a4)
      do i=1,n21
         write(6,*) h(i,1),h(i,2)
      enddo

      return
      end subroutine structure


      subroutine sort(n2,hh,mm,sseq,k)
!     It sorts out the hh values and permutes the mm and sseq values
!     accordingly, so that the rows of hh(i), mm(i) and sseq(i) values
!     remain the same as the corresponding original rows, but with
!     hh(i) values sorted.

      implicit none
      integer n2,i,k,max_h,max_m,max_seq,j_max,j
      integer hh(n2),mm(n2),sseq(n2)

      do i=1,k-1
         max_h=hh(i)
         max_m=mm(i)
         max_seq=sseq(i)
         j_max=i
         do j=i+1,k
            if (hh(j)>max_h) then
               max_h=hh(j)
               max_m=mm(j)
               max_seq=sseq(j)
               j_max=j
            endif
         enddo
         hh(j_max)=hh(i)
         hh(i)=max_h
         mm(j_max)=mm(i)
         mm(i)=max_m
         sseq(j_max)=sseq(i)
         sseq(i)=max_seq
      enddo
      return
      end subroutine sort
               

      subroutine test1(n2,seq,iswitch)
!     Subroutine test1 finds the minimum image and searches whether
!     or not the section 1111 of the 2-2 path sequence appears in
!     the circuit. If it finds that the 1111 section appears, it
!     returns iswitch=0.

      implicit real*8 (a-h,o-z)
      integer seq(n2)
      integer, allocatable :: seq11(:,:)

      allocate(seq11(2*n2,n2))

      iswitch=1

      nmax=0
      do k=0,n2-1
!        permutations with l_1=l_min
         nmax=nmax+1
         do j=1,n2
            j0=mod(j-1+k,n2)+1
            seq11(nmax,j)=seq(j0)
         enddo
!        inversions with l_1=l_min
         nmax1=nmax+1
         seq11(nmax1,1)=seq11(nmax,1)
         n22=n2+2
         do j=2,n2
            seq11(nmax1,j)=seq11(nmax,n22-j)    
         enddo
         nmax=nmax1
      enddo
!     find the minimum image sequence
      do j=1,n2 ! j digit of the sequence
         l_min=4
         do i=1,nmax
            if (seq11(i,j)<l_min) l_min=seq11(i,j)
         enddo

         kcount=0
         do k=1,nmax ! different images
            if (seq11(k,j)/=l_min) then
               kcount=kcount+1    ! counts positions in the j-th digit
            else
               if (kcount/=0) then
                  seq11(k-kcount,:)=seq11(k,:)
               endif
            endif
         enddo
         nmax=nmax-kcount
         if (nmax==1) then
            exit
         endif
      enddo
      seq(:)=seq11(1,:)

      if (seq(1)==1) then
         if (seq(2)==1.and.seq(3)==1.and.seq(4)==1) iswitch=0
         return
      endif
      return
      end subroutine test1

      integer function cos1(m)
      if (m==0.or.m==3) then
         cos1=0
      elseif (m==1.or.m==2) then
         cos1=-1
      elseif (m==4.or.m==5) then
         cos1=1
      endif
      return
      end

      integer function cos2(m)
      if (m==0.or.m==1) then
         cos2=1
      elseif (m==2.or.m==5) then
         cos2=0
      elseif (m==3.or.m==4) then
         cos2=-1
      endif
      return
      end

