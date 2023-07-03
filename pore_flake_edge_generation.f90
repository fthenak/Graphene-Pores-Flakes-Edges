      program flake_generation
!     The program finds all sequences corresponding to a pore, flake or edge 
!     for a given number n2 of 2-fold coordinated atoms in their boundary.
!     Each row of the output includes (i) the pore, flake or edge number (n2:m),
!     (ii) the minimum image sequence, (iii) the length (L) of the flake circuit
!     and 
!     A. for pores, (iv) the number of vacancies (n_vac) and (v) the number of 
!     hexagons in the vacuum area (n_h) and
!     B. for flakes, (iv) the number of atoms (natom) and (v) hexagons (n_hex) in
!     the flake.

!     Zacharias G. Fthenakis

!     Istituto Nanoscienze-CNR, Piazza San Silvestro 12, 56127 Pisa, Italy
!     and
!     Theoretical and Physical Chemistry Institute, National Hellenic
!     Research Foundation, GR-11635, Athens, Greece

      implicit none
      integer n2_max
      parameter (n2_max=30)
      integer seq(n2_max),count(4,4),LL(n2_max),first(n2_max)
      integer L0,i,L,n2,n21,system,n_first,n_dual
      character (LEN=100) FMT,FMT2

      print *, 'System: 1 for pore, 2 for flake, 3 for edge'
      read(5,*) system

      print *, 'Give the number (n2) of 2-fold coordinated atoms'
      print *, 'at the boundary.' 
      print *, 'For pores, n2>=3. For flakes, n2>=6'

      read(5,*) n2

      if (system==1.and.n2<3) then
         print *, 'Wrong value for n2. n2 must be at least 3 for a pore. Try again'
         stop
      elseif (system==2.and.n2<6) then
         print *, 'Wrong value for n2. n2 must be at least 6 for a flake. Try again'
         stop
      endif

      if (system==1) then
         L0=2*n2+6
      elseif (system==2) then
         L0=2*n2-6
      elseif (system==3) then
         L0=2*n2
      else
         print *, 'Wrong system selection. Try again'
         stop
      endif

      if (system==2) then
         if (n2==6) then
            write(6,*) 'n2:m      sequence N n_hex L'
            write(6,*) '6:1        111111  6 1     6'
            print *, 'Number of flakes having a dual pore = ', 0
            stop
         elseif (n2<6) then
            write(6,*) 'n2 must be >= 6'
            stop
         endif
      endif

      if (system==3) then
         call multipliers(n2,first,n_first)
      endif

      n21=n2+1
      write(FMT,*) n2
      write(6,*) 'n2:m      sequence N n_hex L'
      seq(1)=1
      LL(1)=1
      i=0
      n_dual=0
      if (system==1) then
         call do_loops(system,n2,n21,L0,1,2,LL,seq,L,i,FMT,first,n_first,n_dual)
         seq(1)=2
         LL(1)=2
         call do_loops(system,n2,n21,L0,2,2,LL,seq,L,i,FMT,first,n_first,n_dual)
         if (n2==3) then
            write(6,*) '3:      1  *   444   1   3  12'
            print *, 'Number of pores having a dual flake = ', 1
         elseif (n2==4) then
            write(6,*) '4:      1  *   3434   2   4  14'
            print *, 'Number of pores having a dual flake = ', 1
         elseif (n2==6) then
            write(6,*) '6:      5  *   333333   6   7  18'
            print *, 'Number of pores having a dual flake = ', i-n_dual+1
            print *, 'In the above table, * means that the pore has a dual.'
         else
            print *, 'Number of pores having a dual flake = ', i-n_dual
            print *, 'In the above table, * means that the pore has a dual.'
         endif
      else if (system==2) then
         seq(2)=1
         LL(2)=2
         call do_loops(system,n2,n21,L0,1,3,LL,seq,L,i,FMT,first,n_first,n_dual)
         seq(2)=2
         LL(2)=3
         call do_loops(system,n2,n21,L0,1,3,LL,seq,L,i,FMT,first,n_first,n_dual)
         print *, 'Number of flakes having a dual pore = ', i-n_dual
         print *, 'In the above table, * means that the flake has a dual.'
      elseif (system==3) then
         call do_loops(system,n2,n21,L0,1,2,LL,seq,L,i,FMT,first,n_first,n_dual)
      endif
 
      end program

      recursive subroutine do_loops(system,n2,n21,L0,k,n0,LL,seq,L,i,FMT,first,n_first,n_dual)
!     (a) Search for all posible sequences that may represent a pore,
!     flake or edge. 
!     (b) Checks if (i) L=2*n2-6 for flakes, (ii) L=2*n2+6 for pores and
!     (iii) L=2n2 for edges
!     (c) It calls the subroutines which check if a sequence represents
!     a pore, flake or edge.

      implicit none
      integer n2_max,n21_max,n2,n21,sum_h,n_vac,hx_min,hx_max,L_d
      integer natom
      parameter (n2_max=30,n21_max=31,L_d=36)
      integer k,n,i_start,n0,seq(n2_max),LL(n2_max),L,i
      integer L0,j,iswitch,h(n21_max,2),m(n21_max),system
      integer first(n2_max),n_first,idual,n_dual
      integer seq_d(0:L_d),n3,seq_dd(L_d),ij
      character (LEN=20) FMT

      do j=k,4
         seq(n0)=j
         LL(n0)=LL(n0-1)+seq(n0)
         if (n0==n2) then
            L=LL(n2)
            if (L==L0) then
               call test1(n2,seq,i,iswitch)
               if (iswitch==1) then
                  call test2(system,n2,n21,L0,seq,i,iswitch,h,m)
                  if (iswitch==1) then 
                     call touch(n2,n21,h,m,seq,iswitch)
                     if (iswitch==1) then
                        if (system==1) then
                           i=i+1
                           call vacancies(n2,n21,h,m,seq,sum_h,hx_min,hx_max,n_vac)
                           call dual(n2,n21,h,m,seq,idual)
                           if (idual==1) then ! it does not have a dual
                              n_dual=n_dual+1
             write(6,"(i2,':',i7,2x,"//adjustl(FMT)//"I1,3(2x,i2))") n2,i,seq(1:n2),n_vac,sum_h,L
!             write(6,"(i2,':',i7,2x,' & ',"//adjustl(FMT)//"I1,' & ',3(2x,i2,' & '))") n2,i,seq(1:n2),n_vac,sum_h,L
                           else !it has a dual
!                             call find_dual(n2,L0,seq,seq_d,n3)
!                             seq_dd(1:n3)=seq_d(1:n3)
!                             call find_min_image(n3,seq_dd,iswitch)
             write(6,"(i2,':',i7,2x,' * ',"//adjustl(FMT)//"I1, &
     &       3(2x,i2))") n2,i,seq(1:n2),n_vac,sum_h,L
!             write(6,"(i2,':',i7,2x,' & * & ',"//adjustl(FMT)//"I1, &
!     &      ' & ',3(2x,i2,' & '))") n2,i,seq(1:n2),n_vac,sum_h,L
!             write(6,"(i2,':',i7,2x,' & * & ',"//adjustl(FMT)//"I1, &
!     &      ' & ',20i1,' & ',3(2x,i2,' & '))") n2,i,seq(1:n2),seq_dd(1:n3),n_vac,sum_h,L
!                             print *, 'dual = ',seq_d(1:n3)
!             write(6,"(i2,':',i7,2x,' & * & ',"//adjustl(FMT)//"I1,' & ', &
! & "//adjustl(FMT2)//"I1,' & ',3(2x,i2,' & '))") n2,i,seq(1:n2),seq_d(1:ll_d),n_vac,sum_h,L
                           endif
                        elseif (system==2) then
                           i=i+1
                           call flake(n2,n21,h,m,seq,sum_h,hx_min,hx_max,natom)
                           call dual(n2,n21,h,m,seq,idual)
                           if (idual==1) then ! it does not have a dual
                              n_dual=n_dual+1
             write(6,"(i2,':',i7,2x,"//adjustl(FMT)//"I1,3(2x,i2))") n2,i,seq(1:n2),natom,sum_h,L
!             write(6,"(i2,':',i7,' & ',2x,"//adjustl(FMT)//"I1,' & ',3(2x,i2,' & '))") n2,i,seq(1:n2),natom,sum_h,L
                           else ! it has a dual
             write(6,"(i2,':',i7,' * ',2x,"//adjustl(FMT)//"I1,3(2x,i2))") n2,i,seq(1:n2),natom,sum_h,L
!             write(6,"(i2,':',i7,' & * & ',2x,"//adjustl(FMT)//"I1,' & ',3(2x,i2,' & '))") n2,i,seq(1:n2),natom,sum_h,L
                           endif
                        elseif (system==3) then
                           call periodicity(n2,seq,first,n_first,iswitch)
                           if (iswitch==1) then
                              call touch_e(n2,n21,h,m,seq,iswitch)
                              if (iswitch==1) then 
                                 i=i+1
             write(6,"(i2,':',i10,1x,"//adjustl(FMT)//"I1,1(2x,i2))") n2,i,seq(1:n2),L
!            write(6,"(i2,':',i10,1x,'& ',"//adjustl(FMT)//"I1,' &',1(2x,i2))") n2,i,seq(1:n2),L
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         else
            call do_loops(system,n2,n21,L0,k,n0+1,LL,seq,L,i,FMT,first,n_first,n_dual)
         endif
      enddo
      return
      end

      subroutine test2(system,n2,n21,L0,seq,n,iswitch,h,m)
!     a. Checking whether or not the minimum image sequence represents a
!     circuit (i.e. a closed path). 
!     b. Checking wether or not the pore or flake sequence represents a non crossed
!     circuit, according to the first criterion.

      implicit none
      integer n2_max,n21_max,system
      parameter(n2_max=30,n21_max=n2_max+1)
      integer seq(n2_max),m(n21_max),h(n21_max,2),cos1,cos2
      integer iswitch,n2,n21,L0,n,i,k,j,mm
      integer h2(n21_max,2)

      iswitch=1

      m(1)=0
      h(1,1)=0
      h(1,2)=0
!      write(6,504) 'sequence = ', (seq1(i),i=1,n2)
! 504  format(a11,2x,8i1)
!      print *, '------start--------'
!      write(6,400) 1,m(1),h(1,1),h(1,2)
      do i=2,n2
         m(i)=m(i-1)+seq(i)-2
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
      h(n21,1)=h(n2,1)+cos1(m(n2))
      h(n21,2)=h(n2,2)+cos2(m(n2))
     
      if (system==1.or.system==2) then 
      if (m(n21)/=m(1).or.h(n21,1)/=h(1,1).or.h(n21,2)/=h(1,2)) then
!        print *, 'it does not close'
         iswitch=0
         return
      endif
      endif

      do i=1,n2-1
         do j=i+1,n2
            if (h(i,1)==h(j,1).and.h(i,2)==h(j,2)) then
               if (seq(i).ne.1.or.seq(j).ne.1) then
                  iswitch=0 
!                 print *, 'error 1'
                  return
               else
                  mm=m(i)-m(j)
                  if (mm/=3.and.mm/=-3) then
                     iswitch=0
!                 print *, 'error 2'
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
                  if (seq(i).ne.1.or.seq(j).ne.1) then
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


      subroutine test1(n2,seq,n,iswitch)
!     Checking (i) if the pore, flake or edge sequence candidate 
!                  represents its minimum image, and
!              (ii) if the segment 1111 apears in the sequence  
      implicit none
      integer n2_max
      parameter(n2_max=30)
      integer n2,n,iswitch,l_min,nmax,k,i,j,j0,n22,kcount,nmax1
      integer seq(n2),seq11(2*n2,n2)
      

      iswitch=1

      if (seq(1)==1) l_min=1
      if (seq(1)==2) l_min=2
      nmax=0
      do k=0,n2-1
         if (seq(k+1)==l_min) then
            nmax=nmax+1
!           permutations with l_1=l_min
            do j=1,n2
               j0=mod(j-1+k,n2)+1
               seq11(nmax,j)=seq(j0)
            enddo
!           inversions with l_1=l_min
            nmax1=nmax+1
            seq11(nmax1,1)=seq11(nmax,1)
            n22=n2+2
            do j=2,n2
               seq11(nmax1,j)=seq11(nmax,n22-j)    
            enddo
            nmax=nmax1
         endif
      enddo
!     find minimum sequence
      do j=2,n2 ! j digit of the sequence

         l_min=4
         do i=1,nmax
            if (seq11(i,j)<l_min) l_min=seq11(i,j)
         enddo

         if (seq(j).ne.l_min) then
            iswitch=0
            return
         endif

         kcount=0
         do k=1,nmax ! different images
            if (seq11(k,j)/=l_min) then
               kcount=kcount+1    ! counts positions in the j-th digit
            else
               if (kcount.ne.0) then
                  seq11(k-kcount,:)=seq11(k,:)
               endif
            endif
         enddo
         nmax=nmax-kcount
         if (nmax==1) then
            exit
         endif
      enddo

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


      subroutine vacancies(n2,n21,h,m,seq,sum_h,hx_min,hx_max,n_vac)
!     For pores:
!     (a) It finds the hexagons in the vacuum pore area and stores
!     their coordinates in h_vac.
!     (b) It finds the number of vacancies (n_vac).
!     input: h, m, seq - output: h_vac, n_vac

      implicit none
      integer n2_max,n21_max
      parameter (n2_max=30,n21_max=31)
      integer n2,n21,i,j,k,ij,ijk,ii,kk,ll,h(n21_max,2),m(n21_max)
      integer hx_min,hx_max,sum_h,n_vac,up,down
!     integer h_vac(100,2)
      integer hh(n2_max),sseq(n2_max),mm(n2_max),signe(n2_max)
      integer seq(n2_max),mmm

      hx_min=h(1,1)
      hx_max=h(1,1)
      do i=2,n2
         if (h(i,1)>hx_max) hx_max=h(i,1)
         if (h(i,1)<hx_min) hx_min=h(i,1)
      enddo
      sum_h=0
      ijk=0
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
         enddo
         do kk=1,k
            if (signe(kk)==-1) up=hh(kk)
            if (signe(kk)==1) then
               down=hh(kk)
               sum_h=sum_h+up-hh(kk)+1
!   We don't need h_vac in this code
!              do ij=up,down,-1
!                 ijk=ijk+1
!                 h_vac(ijk,1)=i
!                 h_vac(ijk,2)=ij
!              enddo
            endif
            if (signe(kk)==0) then
               sum_h=sum_h+1
!              ijk=ijk+1
!              h_vac(ijk,1)=i
!              h_vac(ijk,2)=hh(kk)
            endif
         enddo
      enddo
      n_vac=2*sum_h-n2-2
      return
      end subroutine vacancies

      subroutine flake(n2,n21,h,m,seq,sum_h,hx_min,hx_max,natom)
!     (a) It finds the number of hexagons in the flake.
!     (b) It finds the number of atoms (natom).
!     input: h, m, seq - output: h_flake, natom

      implicit none
      integer n2_max,n21_max
      parameter (n2_max=30,n21_max=31)
      integer n2,n21,i,j,k,ij,ijk,ii,kk,ll,h(n21_max,2),m(n21_max)
      integer hx_min,hx_max,sum_h,natom,up,down
!     integer h_flake(1000,2)
      integer hh(n2_max),sseq(n2_max),mm(n2_max),signe(n2_max)
      integer seq(n2_max)

      hx_min=h(1,1)
      hx_max=h(1,1)
      do i=2,n2
         if (h(i,1)>hx_max) hx_max=h(i,1)
         if (h(i,1)<hx_min) hx_min=h(i,1)
      enddo
      hx_min=hx_min+1
      hx_max=hx_max-1
      sum_h=0
      ijk=0
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

         if (k.gt.1) then
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
         enddo
         do kk=1,k
            if (signe(kk)==1) up=hh(kk)
            if (signe(kk)==-1) then
               down=hh(kk)
               sum_h=sum_h+up-down-1
            elseif (signe(kk)==0) then
               down=hh(kk)
               sum_h=sum_h+up-down-1
               up=hh(kk)
            endif
         enddo
      enddo
      natom=2*sum_h+n2-2
      return
      end subroutine flake

      subroutine sort(n2,hh,mm,sseq,k)
!     It sorts out the hh values and permutes the mm and sseq values
!     accordingly, so that the rows of hh(i), mm(i) and sseq(i) values
!     remain the same as the corresponding original rows, but with
!     hh(i) values sorted.

      implicit none
      integer n2_max
      parameter (n2_max=30)
      integer n2,i,k,max_h,max_m,max_seq,j_max,j
      integer hh(n2_max),mm(n2_max),sseq(n2_max)

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

      subroutine multipliers(n2,first,n_first)
!     It finds the multipliers of n2, to check if the period of an edge
!     is smaller than the given one.
!
      implicit none
      integer n2_max,n21_max
      parameter (n2_max=30,n21_max=31)
      integer n2,i,first(n2_max),n_first,sq_n2

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
      integer n2_max,n21_max
      parameter (n2_max=30,n21_max=31)
      integer seq(n2_max),n2,i,j,k,first(n2_max),n_first,sq_n2,iswitch
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


      subroutine dual(n2,n21,h,m,seq,idual)
!     It finds of the pore or flake has a dual
!     If idual=1, the pore or flake does not have a dual.
!     If idual=0, the pore or flake has a dual.
      implicit none
      integer n2_max,n21_max
      parameter(n2_max=30, n21_max=31)
      integer n2,n21,i,j,h(n21_max,2),m(n21_max),seq(n2_max)
      integer idual,mij1,mij2,diff1,diff2

      idual=0

      if (seq(1)==1.and.seq(2)==1.and.seq(3)==1) then
         idual=1
         return
      endif

      do i=1,n2-1
         do j=i+1,n2
            if (h(j,1)==h(i,1)+2.and.h(j,2)==h(i,2)-1) then
               mij1=0
               diff1=m(i)-mij1
               if (diff1<0) diff1=diff1+6
               if (diff1<=seq(i)-2.and.seq(i)>1) then
                  mij2=3
                  diff2=m(j)-mij2
                  if (diff2<0) diff2=diff2+6
                  if (diff2<=seq(j)-2.and.seq(j)>1) then
                     idual=1
                     return
                  endif
               endif
            elseif (h(j,1)==h(i,1)+1.and.h(j,2)==h(i,2)+1) then
               mij1=1
               diff1=m(i)-mij1
               if (diff1<0) diff1=diff1+6
               if (diff1<=seq(i)-2.and.seq(i)>1) then
                  mij2=4
                  diff2=m(j)-mij2
                  if (diff2<0) diff2=diff2+6
                  if (diff2<=seq(j)-2.and.seq(j)>1) then
                     idual=1
                     return
                  endif
               endif
            elseif (h(j,1)==h(i,1)-1.and.h(j,2)==h(i,2)+2) then
               mij1=2
               diff1=m(i)-mij1
               if (diff1<0) diff1=diff1+6
               if (diff1<=seq(i)-2.and.seq(i)>1) then
                  mij2=5
                  diff2=m(j)-mij2
                  if (diff2<0) diff2=diff2+6
                  if (diff2<=seq(j)-2.and.seq(j)>1) then
                     idual=1
                     return
                  endif
               endif
            elseif (h(j,1)==h(i,1)-2.and.h(j,2)==h(i,2)+1) then
               mij1=3
               diff1=m(i)-mij1
               if (diff1<0) diff1=diff1+6
               if (diff1<=seq(i)-2.and.seq(i)>1) then
                  mij2=0
                  diff2=m(j)-mij2
                  if (diff2<0) diff2=diff2+6
                  if (diff2<=seq(j)-2.and.seq(j)>1) then
                     idual=1
                     return
                  endif
               endif
            elseif (h(j,1)==h(i,1)-1.and.h(j,2)==h(i,2)-1) then
               mij1=4
               diff1=m(i)-mij1
               if (diff1<0) diff1=diff1+6
               if (diff1<=seq(i)-2.and.seq(i)>1) then
                  mij2=1
                  diff2=m(j)-mij2
                  if (diff2<0) diff2=diff2+6
                  if (diff2<=seq(j)-2.and.seq(j)>1) then
                     idual=1
                     return
                  endif
               endif
            elseif (h(j,1)==h(i,1)+1.and.h(j,2)==h(i,2)-2) then
               mij1=5
               diff1=m(i)-mij1
               if (diff1<0) diff1=diff1+6
               if (diff1<=seq(i)-2.and.seq(i)>1) then
                  mij2=2
                  diff2=m(j)-mij2
                  if (diff2<0) diff2=diff2+6
                  if (diff2<=seq(j)-2.and.seq(j)>1) then
                     idual=1
                     return
                  endif
               endif
            endif
         enddo
      enddo
      return
      end subroutine dual


      subroutine find_dual(n2,L0,seq,seq_d,n3)
!     Given the sequence of a pore or flake (seq), it finds the sequence
!     of the dual (seq_d).
      implicit none
      integer n2_max,n21_max, L_d
!     L_d = n2_max+6 
      parameter (n2_max=30,n21_max=31,L_d=36)
      integer seq(n2_max),seq_d(0:L_d),i,j,k
      integer n2,sum,L0,d(L_d),n3

      seq_d=0
      n3=L0-n2
      k=0
      do i=1,n2
         k=k+1
         d(k)=3
         if (seq(i)>1) then
            do j=2,seq(i)
               k=k+1
               d(k)=2
            enddo
         endif
      enddo

      k=0
      sum=0
      do i=1,L0
         if (d(i)==2) then
            seq_d(k)=sum
            k=k+1
            sum=1
         elseif (d(i)==3) then
            sum=sum+1
         endif
      enddo
      seq_d(n3)=sum+seq_d(0)
      return
      end subroutine find_dual



      subroutine find_min_image(n2,seq,iswitch)
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
      end subroutine find_min_image

      subroutine touch(n2,n21,h,m,seq,iswitch)
!     It checks if a circuit or a path is a crossed circuit or path
!     according to the second criterion.
      implicit none
      integer n2_max, n21_max
      parameter (n2_max=30,n21_max=31)
      integer h(n21_max,2),m(n21_max),seq(n2_max),iswitch
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
      integer n2_max, n21_max
      parameter (n2_max=30,n21_max=31)
      integer h(n21_max,2),m(n21_max),seq(n2_max),iswitch
      integer h2(n21_max,2)
      integer i,j,k,n2,n21,mij1,mij2,d1,d2

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
        

