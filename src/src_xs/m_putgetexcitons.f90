module m_putgetexcitons
  use modmpi
  use modinput, only: input
  use modbse
  use modscl
  use mod_kpoint, only: vkl
  use modxs, only: vkl0, vqlmt
  use m_getunit
  use m_genfilname

  implicit none

  ! Module status
  logical ::  excitons_allocated = .false.

  ! Read in quantities
  logical :: fcoup_, fesel_
  integer(4) :: nk_max_, nk_bse_
  integer(4) :: hamsize_, nexcstored_, iex1_, iex2_
  integer(4) :: iq_
  integer(4) :: ioref_, iuref_
  integer(4), allocatable :: ik2ikqmtp_(:), ik2ikqmtm_(:), ikqmtm2ikqmtp_(:)
  integer(4), allocatable :: kousize_(:), koulims_(:,:)
  integer(4), allocatable :: smap_(:,:), smap_rel_(:,:)
  real(8), allocatable :: vqlmt_(:)
  integer(4), allocatable :: ngridk_(:), ikmap_(:,:,:)
  real(8), allocatable :: vkl_(:,:), vkl0_(:,:)
  real(8), allocatable :: evals_(:)
  complex(8), allocatable :: rvec_(:, :), avec_(:,:)

  contains

    subroutine clear_excitons

      if(allocated(ik2ikqmtp_)) deallocate(ik2ikqmtp_)
      if(allocated(ik2ikqmtm_)) deallocate(ik2ikqmtm_)
      if(allocated(ikqmtm2ikqmtp_)) deallocate(ikqmtm2ikqmtp_)
      if(allocated(kousize_)) deallocate(kousize_)
      if(allocated(koulims_)) deallocate(koulims_)
      if(allocated(smap_)) deallocate(smap_)
      if(allocated(smap_rel_)) deallocate(smap_rel_)
      if(allocated(ngridk_)) deallocate(ngridk_)
      if(allocated(ikmap_)) deallocate(ikmap_)
      if(allocated(vkl0_)) deallocate(vkl0_)
      if(allocated(vkl_)) deallocate(vkl_)
      if(allocated(vqlmt_)) deallocate(vqlmt_)
      if(allocated(evals_)) deallocate(evals_)
      if(allocated(rvec_)) deallocate(rvec_)
      if(allocated(avec_)) deallocate(avec_)

      excitons_allocated = .false.

    end subroutine clear_excitons

    subroutine put_excitons(evals, rvec, avec, iqmt, a1, a2)
      use mod_kpoint, only: ikmap
      use modinput
      use mod_hdf5
      use xhdf5
      use os_utils

      implicit none

      ! I/O
      real(8), intent(in) :: evals(:)
      complex(8), intent(in) :: rvec(:,:)
      complex(8), intent(in), optional :: avec(:,:)
      integer(4), intent(in), optional :: iqmt, a1, a2

      ! Local
      integer(4) :: stat, unexc, i_exciton
      logical :: fcoup, fesel, fchibarq
      integer(4) :: i1, i2, nexcstored, iq, m, n, ngridk(3)
      real(8), allocatable :: de_sorted(:)

      integer(8) :: mypos

      character(:), allocatable :: group, gname,  char_iqpoint, char_iexciton, tdastring, bsetypestring, scrtypestring
      character(256) :: fname

      type(xhdf5_type) :: h5

      ngridk = input%groundstate%ngridk
      
      m = size(rvec,1)
      n = size(rvec,2)

      ! Input checking
      if(size(evals) /= size(rvec,2)) then 
        write(*,'("Error(put_excitons): Array sizes invalid")')
        write(*,'("  size(evals)=",i8," shape(rvec)=", 2i8)') size(evals), shape(rvec)
        call terminate
      end if
      if(present(avec)) then
        if(any(shape(rvec) /= shape(avec))) then 
          write(*,'("Error(put_excitons): Array sizes invalid")')
          write(*,'("  shape(rvec)=",2i8," shape(avec)=", 2i8)') shape(rvec), shape(avec)
          call terminate
        end if
      end if
      if(present(a1) .and. .not. present(a2) .or. .not. present(a1) .and. present(a2)) then
        write(*,'("Error(put_excitons): Specify upper and lower index.")')
        call terminate
      end if

      if(present(a1) .and. present(a2)) then 
        i1=a1
        i2=a2
        nexcstored = i2-i1+1
      else
        i1=1
        i2=n
        nexcstored = i2-i1+1
      end if

      if(i1 > i2 .or. i1<1 .or. i2<1 .or. nexcstored /= n) then 
        write(*,'("Error(put_excitons): Index range invalid")')
        call terminate
      end if

      if(present(iqmt)) then 
        iq = iqmt
      else
        iq = 1
      end if
      if(iq < 1 .or. iq > size(input%xs%qpointset%qpoint, 2)) then
        write(*,'("Error(put_excitons): Q point index invalid")')
        call terminate
      end if

      ! Was the Coulomb potential truncated?
      fchibarq = input%xs%bse%chibarq

      ! BSE type
      fcoup = input%xs%bse%coupling 
      if(any(input%xs%bse%nstlbse == 0)) then
        fesel = .true.
      else
        fesel = .false.
      end if

      ! Generate file name  
      if(fcoup) then
        tdastring=''
      else
        if(fchibarq) then 
          tdastring='-TDA-BAR'
        else
          tdastring='-TDA'
        end if
      end if
      if(input%xs%bse%bsetype == 'IP') then
        tdastring=''
      end if
      bsetypestring = '-' // trim(input%xs%bse%bsetype) // trim(tdastring)
      scrtypestring = '-' // trim(input%xs%screening%screentype)

      if(input%xs%BSE%brixshdf5) then
        call h5%initialize(fhdf5, mpiglobal)

        allocate(character(len=4) :: char_iqpoint)
        write(char_iqpoint, "(I4.4)") iq

        group = 'eigvec' // bsetypestring // scrtypestring
        call h5%initialize_group('/', group)
        call h5%initialize_group(group, char_iqpoint)
        group = join_paths(group, char_iqpoint)

        call h5%write(group, 'evals', evals, [1], shape(evals))
        de_sorted = de(ensortidx)
        call h5%write(group, "evalsIP", de_sorted, [1], shape(de))

        call h5%initialize_group(group, 'parameters')
        group = join_paths(group, 'parameters')

        call h5%write(group, 'fcoup', fcoup)    ! Was the TDA used?
        call h5%write(group, 'fesel', fesel)
        call h5%write(group, 'nk_max', nk_max)
        call h5%write(group, 'nk_bse', nk_bse)
        call h5%write(group, 'hamsize', hamsize)
        call h5%write(group, 'iq', iq)
        call h5%write(group, 'vqlmt(iq)', vqlmt(:, iq))
        call h5%write(group, 'ngridk', ngridk)
        call h5%write(group, 'ikmap', ikmap(0:, 0:, 0:))
        call h5%write(group, 'vkl0', vkl0)
        call h5%write(group, 'vkl', vkl)
        call h5%write(group, 'ik2ikqmtp', ik2ikqmtp)
        call h5%write(group, 'ik2ikqmtm', ik2ikqmtm)
        call h5%write(group, 'ikqmtm2ikqmtp', ikqmtm2ikqmtp)
        call h5%write(group, 'kousize', kousize)
        call h5%write(group, 'koulims', koulims)
        call h5%write(group, 'smap', smap)
        call h5%write(group, 'smap_rel', smap_rel)
        call h5%write(group, 'nexcstored', nexcstored)
        call h5%write(group, 'i1', i1)
        call h5%write(group, 'i2', i2)
        call h5%write(group, 'ioref', ioref)
        call h5%write(group, 'iuref', iuref)
        call h5%write(group, 'ensortidx', ensortidx)

        group = join_paths('eigvec' // bsetypestring // scrtypestring, char_iqpoint)
        call h5%initialize_group(group, 'rvec')
        group = join_paths(group, 'rvec')

        allocate(character(len=8) :: char_iexciton)
        do i_exciton=1, nexcstored
          write(char_iexciton, '(I8.8)') i_exciton 
          call h5%write(group, char_iexciton, rvec(:, i_exciton))
        end do 

        if (present(avec)) then
          group = join_paths('eigvec' // bsetypestring // scrtypestring, char_iqpoint)
          call h5%initialize_group(group, 'avec')
          group = join_paths(group, 'avec')

          do i_exciton=1, nexcstored
            write(char_iexciton, '(I8.8)') i_exciton 
            call h5%write(group, char_iexciton, avec(:, i_exciton))
          end do 
        end if 

        call h5%finalize()

      else

        ! Set filename to EXCCOEFF_*.OUT
        call genfilname(basename='EXCCOEFF', iqmt=iq, bsetype=trim(bsetypestring),&
        & scrtype=trim(scrtypestring), filnam=fname)

        ! Open stream access file 
        call getunit(unexc)
        open(unexc, file=trim(fname), access='stream',&
          & action='write', form='unformatted', status='replace', iostat=stat)
        if(stat /= 0) then
          write(*,*) stat
          write(*,'("Error(storeexcitons): Error creating file ", a)') trim(fname)
          write(*,*)
          call terminate
        end if
        !   Meta data
        write(unexc)&
          & fcoup,&       ! Was the TDA used?
          & fesel,&       ! Where the participating transitions chosen by energy?
          & nk_max,&      ! Number of non-reduced k-points 
          & nk_bse,&      ! Number of k-points used in the bse hamiltonian
          & hamsize,&     ! Size of the RR block of the BSE hamiltonian and number of considered transitions
          & nexcstored,&  ! Number of saved eigenvectors
          & i1, i2,&      ! Range of saved eigenvectors
          & ioref, iuref,&! Reference absolute state index for occpied and unoccupied index (usually lowest and 1st unoccupied)
          & iq,&          ! Index of momentum transfer vector
          & vqlmt(1:3,iq),& ! Momentum transver vector
          & ngridk,&      ! k-grid spacing
          & ikmap,&       ! k-grid index map 3d -> 1d 
          & vkl0,&        ! Lattice vectors for k=k-qmt/2 grid
          & vkl,&         ! Lattice vectors for k'=k+qmt/2 grid
          & ik2ikqmtm(:),& ! ik -> ik-qmt/2 index map
          & ik2ikqmtp(:),& ! ik -> ik+qmt/2 index map
          & ikqmtm2ikqmtp(:),& ! ik-qmt/2 -> ik+qmt/2 index map
          & kousize,&     ! Number of transitions at each k point
          & koulims,&     ! For each k-point, lower and upper c and v index 
          & smap,&        ! Index map  alpha -> c,v,k (absolute c,v,k indices)
          & smap_rel      ! Index map  alpha -> c,v,k (relative c,v,k indices)
        
        !   EVAL and EVEC

        write(unexc) evals

        write(unexc) rvec

        if(present(avec)) then 
          write(unexc) avec
        end if

        close(unexc)
      end if
    end subroutine put_excitons

    subroutine get_excitons(iqmt, a1, a2, e1, e2)
      integer(4), intent(in), optional :: iqmt, a1, a2
      real(8), intent(in), optional :: e1, e2

      integer(4) :: i1, i2, nexcreq, iq, ivec(3)
      real(8) :: r1, r2, vqlmt(3), erange
      real(8), allocatable :: evalstmp(:)
      real(8), parameter :: epslat = 1.0d-6
      integer(4) :: stat, unexc, cmplxlen
      integer(8) :: mypos, pos1, pos2, lambda, i_exciton
      logical :: fcoup, fesel, fex, useindex, useenergy, fchibarq
      complex(8) :: zdummy
      complex(8), allocatable :: rvec2_(:)

      character(256) :: fname

      character(:), allocatable :: tdastring, bsetypestring, scrtypestring
      character(:), allocatable :: group, char_iqpoint, gname,  pos

      if(present(iqmt)) then 
        iq = iqmt
      else
        iq = 1
      end if
      if(iq < 1) then
        write(*,'("Error(get_excitons): Q point index invalid")')
        call terminate
      end if
      
      ! If neither index nor energy range is specified, get all saved excitions.
      ! If energy range is passed [e1,e2), the corresponding exciton indices
      ! will be computed. If a index interval is passed [i1,i2] those excitons
      ! are read.
      useindex=.false.
      useenergy=.false.
      if((present(a1) .or. present(a2)) .and. (present(e1) .or. present(e2))) then 
        write(*,'("Error(get_excitons): Use either indxed or enegy range.")')
        call terminate
      end if
      if(present(a1) .and. .not. present(a2) .or. .not. present(a1) .and. present(a2)) then
        write(*,'("Error(get_excitons): Specify upper and lower index.")')
        call terminate
      end if
      if(present(e1) .and. .not. present(e2) .or. .not. present(e1) .and. present(e2)) then
        write(*,'("Error(get_excitons): Specify upper and lower energy bounds.")')
        call terminate
      end if
      if(present(a1)) then
        useindex=.true.
        useenergy=.false.
      end if
      if(present(e2)) then 
        useindex=.false.
        useenergy=.true.
      end if

      nexcreq = -1
      if(useindex) then 
        if(present(a1) .and. present(a2)) then 
          i1=a1
          i2=a2
          nexcreq = i2-i1+1
        end if
        if(i1 > i2 .or. i1<1 .or. i2<1) then 
          write(*,'("Error(get_excitons): index range invalid")')
          call terminate
        end if
      end if

      if(useenergy) then 
        r1 = e1
        r2 = e2
        erange= r2-r1
        if(erange < 0.0) then 
          write(*,'("Error(get_excitons): erange negative")')
          call terminate
        end if
      end if

      ! BSE type
      fcoup = input%xs%bse%coupling 
      if(any(input%xs%bse%nstlbse == 0)) then
        fesel = .true.
      else
        fesel = .false.
      end if

      vqlmt = input%xs%qpointset%qpoint(:,iq)
      call r3frac(epslat, vqlmt, ivec)

      ! Was the Coulomb potential truncated?
      fchibarq = input%xs%bse%chibarq
      ! Generate file name  
      if(fcoup) then
        tdastring=''
      else
        if(fchibarq) then 
          tdastring="-TDA-BAR"
        else
          tdastring="-TDA"
        end if
      end if
      bsetypestring = '-' // trim(input%xs%bse%bsetype) // tdastring
      scrtypestring = '-' // trim(input%xs%screening%screentype)

      ! Allocate arrays
      call clear_excitons()
      allocate(vqlmt_(3))
      allocate(ngridk_(3))

      ! Set filename to EXCCOEFF_*.OUT
      call genfilname(basename='EXCCOEFF', iqmt=iq, bsetype=trim(bsetypestring),&
        & scrtype=trim(scrtypestring), filnam=fname)
      ! Check if file exists
      inquire(file=trim(fname), exist=fex)
      if(.not. fex ) then
        write(*,*)
        write(*,'("Error(get_excitons): File does not exist:", a)') trim(fname)
        write(*,*)
        call terminate
      end if

      ! Open stream access file 
      call getunit(unexc)
      open(unexc, file=trim(fname), access='stream',&
        & action='read', form='unformatted', status='old', iostat=stat)
      if(stat /= 0) then
        write(*,*) stat
        write(*,'("Error(get_excitons): Error opening file ", a)') trim(fname)
        write(*,*)
        call terminate
      end if
      
      ! Read Meta data
      read(unexc, pos=1)&
        fcoup_, &         ! Was the TDA used?
        fesel_, &         ! Were the transitions selected by energy?
        nk_max_, &        ! Number of non-reduced k-points 
        nk_bse_, &        ! Number of k-points used in the bse hamiltonian
        hamsize_, &       ! Size of the RR block of the BSE hamiltonian and number of considered transitions
        nexcstored_, &    ! Number of saved eigenvectors
        iex1_, iex2_, &   ! Range of saved eigenvectors
        ioref_, iuref_, & ! Reference absolute state index for occpied and unoccupied index (usually lowest and 1st unoccupied)
        iq_, &            ! Index of momentum transfer vector
        vqlmt_, &         ! Momentum transver vector
        ngridk_           ! k-grid spacing
      inquire(unexc, pos=mypos)
  
      ! Check read parameters against requested ones
      if(fcoup_ .neqv. fcoup) then 
        write(*,*)
        write(*,'("Error(get_excitons): BSE type differs")')
        write(*,'(" Requested: fcoup=", l)') fcoup
        write(*,'(" Stored: fcoup_=", l)') fcoup_
        write(*,*)
        call terminate
      end if

      if(fesel_ .neqv. fesel) then 
        write(*,*)
        write(*,'("Error(get_excitons):  Transition selection schemes differs")')
        write(*,'(" Requested: fensel=", l)') fesel
        write(*,'(" Stored: fensel_=", l)') fesel_
        write(*,*)
        call terminate
      end if

      if(norm2(vqlmt-vqlmt_)>epslat) then 
        write(*,*)
        write(*,'("Error(get_excitons): Differing momentum transfer vector.")')
        write(*,'(" Requested: vqlmt=", 3E10.3)') vqlmt
        write(*,'(" Stored: vqlmt_=", 3E10.3)') vqlmt_
        write(*,*)
        call terminate
      end if

      ! Allocate and read meta data arrays
      allocate(ikmap_(0:ngridk_(1)-1,0:ngridk_(2)-1,0:ngridk_(3)-1))
      allocate(vkl0_(3,nk_max_))
      allocate(vkl_(3,nk_max_))
      allocate(ik2ikqmtm_(nk_max_))
      allocate(ik2ikqmtp_(nk_max_))
      allocate(ikqmtm2ikqmtp_(nk_max_))
      allocate(kousize_(nk_max_))
      allocate(koulims_(4,nk_max_))
      allocate(smap_(3,hamsize_))
      allocate(smap_rel_(3,hamsize_))
      allocate(evalstmp(iex1_:iex2_))

      read(unexc, pos=mypos)&
        & ikmap_,&      ! Non reduced k-grid index map 3d -> 1d 
        & vkl0_,&       ! Lattice vectors for k=k-qmt/2 grid
        & vkl_,&        ! Lattice vectors for k'=k+qmt/2 grid
        & ik2ikqmtm_,&  ! ik -> ik-qmt/2 index map
        & ik2ikqmtp_,&  ! ik -> ik+qmt/2 index map
        & ikqmtm2ikqmtp_,& ! ik-qmt/2 -> ik+qmt/2 index map
        & kousize_,&    ! Number of transitions at each k point
        & koulims_,&    ! For each k-point, lower and upper c and v index 
        & smap_,&       ! Index map  alpha -> c,v,k (absolute c,v,k indices)
        & smap_rel_,&   ! Index map  alpha -> c,v,k (relative c,v,k indices)
        & evalstmp      ! Excitonic enegies
      inquire(unexc, pos=mypos)

      ! Inquire ouput length of a complex number (in units of 4 byte by default)
      inquire(iolength=cmplxlen) zdummy
      cmplxlen=cmplxlen*4

      if(useenergy) then 
        call energy2index(size(evalstmp), size(evalstmp),&
          & evalstmp(iex1_:iex2_), r1, r2, i1, i2)
        i1 = i1+iex1_-1
        i2 = i2+iex1_-1
      else if(.not. useindex) then
        i1=iex1_
        i2=iex2_
      end if

      if(i1 < iex1_ .or. i2 > iex2_) then 
        write(*,*)
        write(*,'("Error(get_excitons): Eigenvector range mismatch.")')
        write(*,'(" Requested: i1=", i6," i2=", i6)') i1, i2
        write(*,'(" Stored: iex1_=", i6," iex2_=", i6)') iex1_, iex2_
        write(*,*)
        call terminate
      end if

      ! Eigenvalues
      allocate(evals_(i1:i2))
      evals_(:) = evalstmp(i1:i2)
      deallocate(evalstmp)

      ! Resonant part of the eigenvectors
      allocate(rvec_(hamsize_, i1:i2))
      pos1=int(i1-iex1_,8)*int(cmplxlen,8)*int(hamsize_,8)+mypos
      read(unexc, pos=pos1) rvec_
      if(fcoup_) then  
        ! Anti-resonant part of the eigenvectors
        allocate(avec_(hamsize_, i1:i2))
        pos1=int(iex2_-iex1_+1,8)*int(cmplxlen,8)*int(hamsize_,8)+mypos
        pos2=int(i1-iex1_,8)*int(cmplxlen,8)*int(hamsize_,8)+pos1
        read(unexc, pos=pos2) avec_
      end if
      close(unexc)

      ! Set stored index range
      iex1_ = i1
      iex2_ = i2

      excitons_allocated = .true.

    end subroutine get_excitons

    subroutine putd_excitons(evals, drvec, davec, iqmt, a1, a2)
      use mod_kpoint, only: ikmap
      use modinput
      use xhdf5, only: xhdf5_type
      use mod_hdf5, only: fhdf5
      use os_utils

      implicit none

      ! I/O
      real(8), intent(in) :: evals(:)
      type(dzmat), intent(in) :: drvec
      type(dzmat), intent(in), optional :: davec
      integer(4), intent(in), optional :: iqmt, a1, a2

      ! Local
      type(dzmat) :: dauxmat
      integer(4) :: stat, unexc
      logical :: fcoup, fesel, fchibarq
      integer(4) :: i1, i2, nexcstored, iq, m, n, m2, n2, i_exciton, ngridk(3)
      logical :: sane, distributed
      real(8) :: de_(hamsize)
      real(8), allocatable :: de_sorted(:)

      character(:), allocatable :: fname, char_iqpoint, char_iexciton
      character(:), allocatable :: tdastring, bsetypestring, scrtypestring
      character(:), allocatable :: group,  gname

      type(xhdf5_type) :: h5

      if(present(davec)) then
        sane = drvec%isdistributed .eqv. davec%isdistributed
      else
        sane = .true.
      end if

      distributed = drvec%isdistributed

      if(.not. distributed .and. sane) then 
        if(present(davec)) then 
          call put_excitons(evals(a1:a2), drvec%za(:,a1:a2), avec=davec%za(:,a1:a2),&
           & iqmt=iqmt, a1=a1, a2=a2)
        else
          call put_excitons(evals(a1:a2), drvec%za(:,a1:a2), iqmt=iqmt, a1=a1, a2=a2)
        end if
        return
      else if(.not. sane) then 
        write(*,*) "putd_excitions: input not sane"
        call terminate
      end if

      ngridk = input%groundstate%ngridk

      m = drvec%nrows
      n = drvec%ncols
  
      ! Input checking
      if(present(a1) .and. .not. present(a2) .or. .not. present(a1) .and. present(a2)) then
        if(bicurrent%isroot) then 
          write(*,'("Error(putd_excitons): Specify upper and lower index.")')
        end if
        call terminate
      end if
      if(present(a1) .and. present(a2)) then 
        i1=a1
        i2=a2
        nexcstored = i2-i1+1
      else
        i1=1
        i2=n
        nexcstored = i2-i1+1
      end if
      if(i1 > i2 .or. i1<1 .or. i2<1 .or. nexcstored > n) then 
        if(bicurrent%isroot) then 
          write(*,'("Error(putd_excitons): Index range invalid")')
          write(*,'("  i1=",i8," i2=",i8," nexcstored=", i8," m=",i8," n=",i8)')&
            & i1,i2,nexcstored,m,n
        end if
        call terminate
      end if
      if(size(evals) /= nexcstored) then 
        if(bicurrent%isroot) then 
          write(*,'("Error(putd_excitons): Array sizes invalid")')
          write(*,'("  size(evals)=",i8," nexcstored=",i8," shape(drvec)=", 2i8)')&
            & size(evals), nexcstored, m, n
        end if
        call terminate
      end if
      if(present(davec)) then
        m2 = davec%nrows
        n2 = davec%ncols
        if(m /= m2 .or. n /= n2) then 
          if(bicurrent%isroot) then 
            write(*,'("Error(putd_excitons): Array sizes invalid")')
            write(*,'("  shape(drvec)=",2i8," shape(davec)=", 2i8)') m, n, m2, n2
          end if
          call terminate
        end if
      end if

      if(present(iqmt)) then 
        iq = iqmt
      else
        iq = 1
      end if
      if(iq < 1 .or. iq > size(input%xs%qpointset%qpoint, 2)) then
        if(mpiglobal%rank==0) then 
          write(*,'("Error(putd_excitons): Q point index invalid")')
        end if
        call terminate
      end if

      ! BSE type
      fcoup = input%xs%bse%coupling 
      if(any(input%xs%bse%nstlbse == 0)) then
        fesel = .true.
      else
        fesel = .false.
      end if

      ! Was the Coulomb potential truncated?
      fchibarq = input%xs%bse%chibarq

      ! Root does the writing to file
      if(bicurrent%isroot) then 

        ! Generate file name  
        if(fcoup) then
          tdastring=''
        else
          if(fchibarq) then 
            tdastring='-TDA-BAR'
          else
            tdastring='-TDA'
          end if
        end if
        bsetypestring = '-' // trim(input%xs%bse%bsetype) // tdastring
        scrtypestring = '-' // trim(input%xs%screening%screentype)

        call new_dzmat(dauxmat, m, 1, bi0d)

        if(input%xs%BSE%brixshdf5) then 
          call h5%initialize(fhdf5, mpiglobal, serial_access=.true.)

          group = 'eigvec' // bsetypestring // scrtypestring
          call h5%initialize_group('/', group)

          allocate(character(len=4) :: char_iqpoint)
          write(char_iqpoint, '(I4.4)') iq ! Generate string out of momentum transfer index
          call h5%initialize_group(group, char_iqpoint)
          group = join_paths(group, char_iqpoint)

          call h5%write(group, 'evals', evals)
          de_sorted = de(ensortidx)
          call h5%write(group, "evalsIP", de_sorted)

          call h5%initialize_group(group, 'parameters')
          group = join_paths(group, 'parameters')

          call h5%write(group, 'fcoup', fcoup)    ! Was the TDA used?
          call h5%write(group, 'fesel', fesel)
          call h5%write(group, 'nk_max', nk_max)
          call h5%write(group, 'nk_bse', nk_bse)
          call h5%write(group, 'hamsize', hamsize)
          call h5%write(group, 'iq', iq)
          call h5%write(group, 'vqlmt(iq)', vqlmt(:, iq))
          call h5%write(group, 'ngridk', ngridk)
          call h5%write(group, 'ikmap', ikmap(0:, 0:, 0:))
          call h5%write(group, 'vkl0', vkl0)
          call h5%write(group, 'vkl', vkl)
          call h5%write(group, 'ik2ikqmtp', ik2ikqmtp)
          call h5%write(group, 'ik2ikqmtm', ik2ikqmtm)
          call h5%write(group, 'ikqmtm2ikqmtp', ikqmtm2ikqmtp)
          call h5%write(group, 'kousize', kousize)
          call h5%write(group, 'koulims', koulims)
          call h5%write(group, 'smap', smap)
          call h5%write(group, 'smap_rel', smap_rel)
          call h5%write(group, 'nexcstored', nexcstored)
          call h5%write(group, 'i1', i1)
          call h5%write(group, 'i2', i2)
          call h5%write(group, 'ioref', ioref)
          call h5%write(group, 'iuref', iuref)
          call h5%write(group, 'ensortidx', ensortidx)

          group = join_paths('eigvec' // bsetypestring // scrtypestring, char_iqpoint)
          call h5%initialize_group(group, 'rvec')
          group = join_paths(group, 'rvec')

          allocate(character(len=8) :: char_iexciton)
          do i_exciton=1, nexcstored
            call dzmat_copy(drvec%context, m, 1, dmata=drvec, dmatb=dauxmat, ra=1, ca=i_exciton+drvec%subj-1)
            write(char_iexciton, '(I8.8)') i_exciton 
            call h5%write(group, char_iexciton, dauxmat%za(:, 1))
          end do 

          if (present(davec)) then
            group = join_paths('eigvec' // bsetypestring // scrtypestring, char_iqpoint)
            call h5%initialize_group(group, 'avec')
            group = join_paths(group, 'avec')

            do i_exciton=1, nexcstored
              call dzmat_copy(davec%context, m, 1, dmata=davec, dmatb=dauxmat, ra=1, ca=i_exciton+davec%subj-1)
              write(char_iexciton, '(I8.8)') i_exciton 
              call h5%write(group, char_iexciton, dauxmat%za(:, 1))
            end do 
          end if 

          call h5%finalize()
        else

          ! Set filename to EXCCOEFF_*.OUT
          call genfilname(basename='EXCCOEFF', iqmt=iq, bsetype=trim(bsetypestring),&
            & scrtype=trim(scrtypestring), filnam=fname)

          ! Open stream access file 
          call getunit(unexc)
          open(unexc, file=trim(fname), access='stream',&
            & action='write', form='unformatted', status='replace', iostat=stat)
          if(stat /= 0) then
            write(*,*) stat
            write(*,'("Error(storeexcitons): Error creating file ", a)') trim(fname)
            write(*,*)
            call terminate
          end if

          ! Write data 
          !   Meta data
          write(unexc)&
            & fcoup,&       ! Was the TDA used?
            & fesel,&       ! Where the participating transitions chosen by energy?
            & nk_max,&      ! Number of non-reduced k-points 
            & nk_bse,&      ! Number of k-points used in the bse hamiltonian
            & hamsize,&     ! Size of the RR block of the BSE hamiltonian and number of considered transitions
            & nexcstored,&  ! Number of saved eigenvectors
            & i1, i2,&      ! Range of saved eigenvectors
            & ioref, iuref,&! Reference absolute state index for occpied and unoccupied index (usually lowest and 1st unoccupied)
            & iq,&          ! Index of momentum transfer vector
            & vqlmt(1:3,iq),& ! Momentum transver vector
            & ngridk,&      ! k-grid spacing
            & ikmap,&       ! k-grid index map 3d -> 1d 
            & vkl0,&        ! Lattice vectors for k=k-qmt/2 grid
            & vkl,&         ! Lattice vectors for k'=k+qmt/2 grid
            & ik2ikqmtm(:),& ! ik -> ik-qmt/2 index map
            & ik2ikqmtp(:),& ! ik -> ik+qmt/2 index map
            & ikqmtm2ikqmtp(:),& ! ik-qmt/2 -> ik+qmt/2 index map
            & kousize,&     ! Number of transitions at each k point
            & koulims,&     ! For each k-point, lower and upper c and v index 
            & smap,&        ! Index map  alpha -> c,v,k (absolute c,v,k indices)
            & smap_rel      ! Index map  alpha -> c,v,k (relative c,v,k indices)
          ! Evals
          write(unexc) evals

          do i_exciton= 1, nexcstored
            ! Copy i_exciton'th column of distributed eigenvector matrix to root 
            call dzmat_copy(drvec%context, m, 1, dmata=drvec, dmatb=dauxmat,&
              & ra=1, ca=i_exciton+drvec%subj-1)
            ! Write eigenvector        
            write(unexc) dauxmat%za(1:m,1)
          end do

          if(present(davec)) then 
            do i_exciton= 1, nexcstored
              call dzmat_copy(davec%context, m, 1, dmata=davec, dmatb=dauxmat,&
                & ra=1, ca=i_exciton+davec%subj-1)
              ! Write eigenvector
              write(unexc) dauxmat%za(1:m,1)
            end do
          end if

          close(unexc)
        end if 

        call del_dzmat(dauxmat)
      ! Send to root
      else

        do i_exciton= 1, nexcstored
          ! Copy i_exciton'th column of distributed eigenvector matrix to root 
          call dzmat_copy(drvec%context, m, 1, dmata=drvec,&
            & ra=1, ca=i_exciton+drvec%subj-1)
        end do

        if(present(davec)) then 
          do i_exciton= 1, nexcstored
            call dzmat_copy(davec%context, m, 1, dmata=davec,&
              & ra=1, ca=i_exciton+davec%subj-1)
          end do
        end if

      end if

    end subroutine putd_excitons
    
end module m_putgetexcitons
