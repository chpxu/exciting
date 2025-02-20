!BOP
! !ROUTINE: scf_cycle
! !INTERFACE:
!
!
subroutine scf_cycle(verbosity)
! !USES:
    use modinput
    use modmain
    use modmpi
    use scl_xml_out_Module
    use TS_vdW_module, only: C6ab, R0_eff_ab
    use sirius_init,   only: sirius_options
    use sirius_api,    only: set_radial_functions_sirius, solve_seceqn_sirius, get_eval_sirius, get_evec_sirius, &
                             put_occ_sirius, generate_density_sirius, get_periodic_function_sirius
    use mod_potential_and_density, only: generate_density_and_magnetization
    use trial_energy_selection, only: select_local_orbital_trial_energies, select_apw_trial_energies
    use lo_recommendation, only: recommend_local_orbital_trial_energies
!

! !DESCRIPTION:
!
! !REVISION HISTORY:
!   Created February 2013 (DIN)
!EOP
!BOC
    Implicit None
    integer, intent(IN) :: verbosity
    Real(8) :: et, fm
    Real(8), Allocatable :: evalfv(:, :)
    Complex (8), Allocatable :: evecfv(:, :, :)
    Complex (8), Allocatable :: evecsv(:, :)
    Logical :: exist
    Integer :: ik, is, ia, idm, id, lmax, nodesmax
    Integer :: n, nwork
    !Integer :: i,j, ias
    Real(8), Allocatable :: v(:),forcesum(:,:)
    Real(8) :: timetot, ts0, ts1, tin1, tin0, ta,tb
    character*(77) :: string, acoord

    Real (8), Allocatable :: rhomtref(:,:,:) ! muffin-tin charge density (reference)
    Real (8), Allocatable :: rhoirref(:)     ! interstitial real-space charge density (reference)

    Type (apw_lo_basis_type) :: mt_basis

    acoord = "lattice"
    if (input%structure%cartesian) acoord = "cartesian"

    If ((verbosity>-1).and.(rank==0)) Then
        write(string,'("Self-consistent loop started")')
        call printbox(60,"+",string)
    End If

! reset density-dependent dispersion coefficients of Tkatchenko-Scheffler method
    If (Allocated(C6ab)) Deallocate(C6ab)
    If (Allocated(R0_eff_ab)) Deallocate(R0_eff_ab)
!
    call MTNullify(mt_hscf)

!_______________________________________________________________
! initialise or read the charge density and potentials from file
    call stopwatch("exciting:init_scf", 1)
!! TIME - Begin of initialisation segment
    Call timesec (ts0)
    If ((task == 1) .or. (task == 3)) Then
        Call readstate
        If ((verbosity>-1).and.(rank==0)) write(60,'(" Potential read in from STATE.OUT")')
    Else If (task == 7) Then
        ! restart from previous HYBRID iteration
        call genmeffig()

       continue
    Else If (task == 200) Then
        Call phveff
        If ((verbosity>-1).and.(rank==0)) write(60,'(" Supercell potential constructed from STATE.OUT")')
    Else
        Call timesec(tin0)
        Call rhoinit
        Call timesec(tin1)
        time_density_init=tin1-tin0
        Call timesec(tin0)
        Call poteff
        Call genveffig
        Call timesec(tin1)
        time_pot_init=tin1-tin0
        If ((verbosity>-1).and.(rank==0)) write(60,'(" Density and potential initialised from atomic data")')
    End If

    Call genmeffig
    If ((verbosity>-1).and.(rank==0)) then
        write (60, *)
        Call flushifc (60)
    end if
    call stopwatch("exciting:init_scf", 0)

    call MTNullify(mt_hscf)

!_____________________________
! convergence vectors

    If (allocated(vcurrentconvergence)) deallocate(vcurrentconvergence)
    Allocate(vcurrentconvergence(input%groundstate%niterconvcheck))
    vcurrentconvergence=0.d0
    If (allocated(vdeltae)) deallocate(vdeltae)
    Allocate(vdeltae(input%groundstate%niterconvcheck))
    vdeltae=0.d0
    If (allocated(vchgdst)) deallocate(vchgdst)
    Allocate(vchgdst(input%groundstate%niterconvcheck))
    vchgdst=0.d0

!_____________________________
! reference density

    If (allocated(rhomtref)) deallocate(rhomtref)
    Allocate(rhomtref(lmmaxvr,nrmtmax,natmtot))
    rhomtref(:,:,:) = rhomt(:,:,:)
    If (allocated(rhoirref)) deallocate(rhoirref)
    Allocate (rhoirref(ngrtot))
    rhoirref(:) = rhoir(:)

    Call timesec (ts1)
    timeinit = timeinit+ts1-ts0
!! TIME - End of initialisation segment

!----------------------------------------------------
!! TIME - Mixer segment
    Call timesec (ts0)
    ! size of mixing vector
    n = lmmaxvr*nrmtmax*natmtot+ngrtot
    If (associated(input%groundstate%spin)) n = n*(1+ndmag)
    If (ldapu .Ne. 0) n = n + 2*lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot
    ! allocate mixing arrays
    Allocate (v(n))
    ! call mixing array allocation functions by setting
    nwork = -1
    ! and call interface
    iscl = 0
    Call packeff (.True., n, v)
    If (rank .Eq. 0) Call mixerifc(input%groundstate%mixernumber, n, v, currentconvergence, nwork)
    Call packeff (.False., n, v)
    Call timesec (ts1)
    timemixer = ts1-ts0+timemixer
!! TIME - End of mixer segment
!----------------------------------------------------

! set last iteration flag
    tlast = .False.
! set stop flag
    tstop = .False.
    engytot = 0.d0
    fm = 0.d0
! delete any existing eigenvector files
    If ((rank .Eq. 0) .And. ((task .Eq. 0) .Or. (task .Eq. 2))) Call delevec()

!! TIME - First IO segment
!----------------------------------------!
! begin the self-consistent loop
!----------------------------------------!
    call stopwatch("exciting:scf", 1)
    Do iscl = 1, input%groundstate%maxscl
        call timesec (ts0)
! exit self-consistent loop if last iteration is complete
        if (tlast) then
            If ((verbosity>-1).and.(rank==0)) Then
                call printline(60," ")
                call printline(60,"+")
                If (input%groundstate%niterconvcheck.ge.2) Then
                   write(string,'("Convergency criteria checked for the last", I2," iterations ")') input%groundstate%niterconvcheck
                Else
                   write(string,'("Convergency criteria checked for the last iteration ")')
                End If
                call printtext(60,"+",string)
                write(string,'("Convergence targets achieved. Performing final SCF iteration")')
                call printtext(60,"+",string)
                call printline(60,"+")
                Call flushifc(60)
            End If
        else
            If (iscl .Ge. input%groundstate%maxscl) Then
                If ((verbosity>-1).and.(rank==0)) Then
                    write(string,'("Reached self-consistent loops maximum : ", I4)') &
                   &  input%groundstate%maxscl
                    call printbox(60,"+",string)
                    call warning('Warning(gndstate): Reached self-consistent loops maximum')
                    Call flushifc(60)
                End If
                tlast = .True.
                goto 10
            End If
            If ((verbosity>-1).and.(rank==0)) Then
                write(string,'("SCF iteration number : ", I4)') iscl
                call printbox(60,"+",string)
                Call flushifc(60)
            End If
        end if

10      continue

        call timesec (ts1)
        timeio = ts1 - ts0 + timeio

!! TIME - End of first IO segment

!! TIME - Muffin-tin segment
        Call timesec (ts0)

        if (task /= 7) then
          ! No updates of core and valence radial functions during hybrids run
          call gencore          ! generate the core wavefunctions and densities
          ! find the first linearization energies 
          if (iscl==1) then 
             call select_local_orbital_trial_energies(nlorb, lorbord, lorbl, lorbn, nspecies, idxas, nrmt, spr, veffmt(1,:,:), lorbe0) 
             call select_apw_trial_energies(maxapword, maxlapw, apwn(:, 0:, :), nspecies, idxas, nrmt, spr, veffmt(1,:,:), apwe0(:, 0:, :))
          endif
          call linengy          ! find the new linearization energies
          if (rank==0) call writelinen
          call genapwfr         ! generate the APW radial functions
          call genlofr(tlast)   ! generate the local-orbital radial functions
          ! compute recommendations for local orbital trial energies
          if ((associated(input%groundstate%lorecommendation)) .and. (tlast)) then
            nodesmax = input%groundstate%lorecommendation%nodesmaxlo
            lmax = input%groundstate%lorecommendation%lmaxlo
            call recommend_local_orbital_trial_energies(nodesmax, lmax, nspecies, spsymb, idxas, nrmt, spr, veffmt(1,:,:))
          endif
          call olprad           ! compute the overlap radial integrals
        end if
        if ( associated(input%groundstate%sirius) ) then
          call set_radial_functions_sirius( input )
        end if
!------------------------------------------------------------
! Effective Hamiltonian Setup: Radial and Angular integrals
!------------------------------------------------------------
        call stopwatch("exciting:rad_int", 1)
        call MTInitAll(mt_hscf)
        call hmlint(mt_hscf)
        call stopwatch("exciting:rad_int", 0)
!________________
! partial charges

        if (input%groundstate%tpartcharges) then
            allocate(chgpart(lmmaxvr,natmtot,nstsv))
            chgpart(:,:,:)=0.d0
        end if

        call timesec (ts1)
        timemt = ts1 - ts0 + timemt

!! TIME - End of muffin-tin segment

!! TIME - Second IO segment

        call timesec (ts0)

!-----------------------------------------------
! Solve Secular Equation for each k-point
!-----------------------------------------------
        call timesec(ta)

        if ( associated(input%groundstate%sirius) .and. sirius_options%use_eigen_states ) then
          call solve_seceqn_sirius()
          call get_eval_sirius()
          call get_evec_sirius()
        else
! start k-point loop
#ifdef MPI
            call barrier()
            If (rank == 0) Call delevec()
            splittfile = .True.
            Do ik = firstofset(rank,nkpt), lastofset(rank,nkpt)
#else
            splittfile = .False.
            Do ik = 1, nkpt
#endif

!____________________________________________
! every thread should allocate its own arrays

                Allocate (evalfv(nstfv, nspnfv))
                Allocate (evecfv(nmatmax, nstfv, nspnfv))
                Allocate (evecsv(nstsv, nstsv))

                if (iscl.le.1) then
                  evecfv=0d0
                elseif (input%groundstate%solver%type.eq.'Davidson') then
                  Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
                endif
!! TIME - seceqn does not belong to IO

                call timesec(ts1)
                timeio = ts1 - ts0 + timeio

!__________________________________________________________
! solve the first- and second-variational secular equations
                call seceqn (ik, evalfv, evecfv, evecsv)

                call timesec(ts0)

!______________________________________
! write the eigenvalues/vectors to file

                Call putevalfv (ik, evalfv)
                Call putevalsv (ik, evalsv(:, ik))
                Call putevecfv (ik, evecfv)
                Call putevecsv (ik, evecsv)

!__________________________
! calculate partial charges
                if (input%groundstate%tpartcharges) call genpchgs(ik,evecfv,evecsv)
                deallocate (evalfv, evecfv, evecsv)

            End Do ! ik

! end k-point loop -------------------------------------------------------------

#ifdef MPI
            call mpi_allgatherv_ifc(nkpt, inplace=.False., rlen=nstsv, rbuf=evalsv)
            if (task==7) call mpi_allgatherv_ifc(nkpt, inplace=.False., rlen=nstfv, rbuf=engyknst)
#endif
        end if

        call timesec(tb)
! Release memory used by the MT Hamiltonian
!        call mt_hscf%release()

        call stopwatch("exciting:occ", 1)
!-----------------------------------------------
! find the occupation numbers and Fermi energy
!-----------------------------------------------
        Call occupy
        If (rank==0) Then
! write out the eigenvalues and occupation numbers
            Call writeeval
! write the Fermi energy to file
            Call writefermi
        End If
!write the occupancies to file
#ifdef MPI
        Do ik = firstofset(rank,nkpt), lastofset(rank,nkpt)
#else
        Do ik = 1, nkpt
#endif
            Call putoccsv (ik, occsv(:, ik))
        End Do
        if ( associated(input%groundstate%sirius) ) then
          call put_occ_sirius()
        end if
        call timesec(ts1)
        timeio = ts1 - ts0 + timeio
        call stopwatch("exciting:occ", 0)

        call stopwatch("exciting:rhomag", 1)
!-----------------------------------------------
! Calculate density and magnetization
!-----------------------------------------------
        if (associated(input%groundstate%sirius) .and. sirius_options%use_density ) then
          call generate_density_sirius()
          call get_periodic_function_sirius(rhoir, ngrid)
          call timesec(ts0)
        else
          call generate_density_and_magnetization()
          call timesec(ts0)
#ifdef MPI
        ! EXX case
          If (input%groundstate%xctypenumber.Lt.0) Call mpiresumeevecfiles()
#endif

          if ((input%groundstate%tpartcharges).and.(rank==0)) then
              ! write out partial charges
              call writepchgs(69,input%groundstate%lmaxvr)
              call flushifc(69)
          end if
          call timesec(ts1)
          timeio = ts1 - ts0 + timeio
!! TIME - End of second IO segment

          call timesec(ts0)
! symmetrise the density
          call symrf(input%groundstate%lradstep, rhomt, rhoir)
! symmetrise the magnetisation
          If (associated(input%groundstate%spin)) Call symrvf(input%groundstate%lradstep, magmt, magir)
! convert the density from a coarse to a fine radial mesh
          call rfmtctof (rhomt)
! convert the magnetisation from a coarse to a fine radial mesh
          Do idm = 1, ndmag
              Call rfmtctof (magmt(:, :, :, idm))
          End Do
        end if

! add the core density to the total density
        Call addrhocr
! calculate the charges
        Call charge
! calculate the moments
        If (associated(input%groundstate%spin)) Call moment
! normalise the density
        Call rhonorm
        call stopwatch("exciting:rhomag", 0)
! LDA+U
        If (ldapu .Ne. 0) Then
! generate the LDA+U density matrix
            Call gendmatlu
! generate the LDA+U potential matrix
            Call genvmatlu
! write the LDA+U matrices to file
            if (rank .eq. 0) Call writeldapu
        End If
! generate charge distance
        call chgdist(rhomtref,rhoirref)
        do id=1, input%groundstate%niterconvcheck-1
           vchgdst(id) = vchgdst(id+1)
        end do
        vchgdst(input%groundstate%niterconvcheck) = chgdst

! store density to reference
        rhoirref(:)=rhoir(:)
        rhomtref(:,:,:)=rhomt(:,:,:)

!-----------------------------------
! Compute the effective potential
!-----------------------------------
        Call poteff
!---------------
! Mixing
!---------------
! pack interstitial and muffin-tin effective potential and field into one array
        Call packeff (.True., n, v)
! mix in the old potential and field with the new
        If (rank .Eq. 0) Then
           Call mixerifc (input%groundstate%mixernumber, n, v, currentconvergence, nwork)
           do id=1, input%groundstate%niterconvcheck-1
              vcurrentconvergence(id) = vcurrentconvergence(id+1)
           end do
           vcurrentconvergence(input%groundstate%niterconvcheck) = currentconvergence

        End If
#ifdef MPI
        Call MPI_bcast (v(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
! unpack potential and field
        Call packeff (.False., n, v)
!---------------
! Fourier transform effective potential to G-space
        Call genveffig
        if (allocated(meffig)) deallocate(meffig)
        if (allocated(m2effig)) deallocate(m2effig)
! add the fixed spin moment effect field
        If (getfixspinnumber() .Ne. 0) Call fsmfield
        Call genmeffig
! reduce the external magnetic fields if required
        If (associated(input%groundstate%spin)) Then
            If (input%groundstate%spin%reducebf .Lt. 1.d0) Then
                input%groundstate%spin%bfieldc(:) = &
               &  input%groundstate%spin%bfieldc(:) * input%groundstate%spin%reducebf
                Do is = 1, nspecies
                    Do ia = 1, natoms (is)
                        input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:) = &
                       &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:) * input%groundstate%spin%reducebf
                    End Do
                End Do
            End If
        End If

!--------------------------------
! compute the energy components
!--------------------------------
        et = engytot
        Call energy
        Call timesec(ts1)
        timepot=ts1-ts0+timepot

!----------------------------------------------
! compute the forces (without IBS corrections)
!----------------------------------------------
        if (input%groundstate%tforce) then
            Call force(.false.)
        end if

!-----------------------------
! Print results
!-----------------------------

!! TIME - Third IO segment
        call timesec(ts0)
        deltae=abs(et-engytot)
        do id=1, input%groundstate%niterconvcheck-1
           vdeltae(id)=vdeltae(id+1)
        end do
        vdeltae(input%groundstate%niterconvcheck)=abs(et-engytot)

        If ((verbosity>-1).and.(rank==0)) Then
! output energy components
            call writeengy(60)
            if (verbosity>0) Write (60,*)
            Write (60, '(" DOS at Fermi energy (states/Ha/cell)",T45 ": ", F18.8)') fermidos
! write DOS at Fermi energy to FERMIDOS.OUT and flush
!            Write (62, '(G18.10)') fermidos
!            Call flushifc (62)
! output charges and moments
            Call writechg (60,input%groundstate%outputlevelnumber)
! write total moment to MOMENT.OUT and flush
            If (associated(input%groundstate%spin)) Then
                Write (63, '(3G18.10)') momtot (1:ndmag)
                Call flushifc (63)
            End If
! output effective fields for fixed spin moment calculations
            If (getfixspinnumber() .Ne. 0) Call writefsm (60)
! output forces to INFO.OUT
!            if (input%groundstate%tforce) call writeforce(60,input%relax%outputlevelnumber)
! write band-gap if the dos at the Fermi energy is smaller than the given threshold
            if (fermidos<1.0d-4) call printbandgap(60)
! check for WRITE file
            Inquire (File='WRITE', Exist=exist)
            If (exist) Then
                Write (60,*)
                Write (60, '(" WRITE file exists - writing STATE.OUT")')
                Call writestate
                Open (50, File='WRITE')
                Close (50, Status='DELETE')
            End If
            Call scl_iter_xmlout ()
            If (associated(input%groundstate%spin)) Call scl_xml_write_moments()
            Call scl_xml_out_write()
        End If
! write STATE.OUT file if required
        If (input%groundstate%nwrite .Ge. 1) Then
            If (Mod(iscl, input%groundstate%nwrite) .Eq. 0) Then
                Call writestate
                if ((verbosity>-1).and.(rank==0)) Then
                    write(60,*)
                    write(60, '(" Wrote STATE.OUT")')
                end if
            End If
        End If
        call timesec(ts1)
        timeio = ts1 - ts0 + timeio
!! TIME - End of third IO segment

! exit self-consistent loop if last iteration is complete
        If (tlast) goto 20

! update convergence criteria
        et = engytot
        if (input%groundstate%tforce) then
            dforcemax=abs(fm-forcemax)
            if (dforcemax .lt. 1.d-10) dforcemax=0.
        end if

!! TIME - Fourth IO segment
        call timesec(ts0)

! output the current total time
        timetot = timeinit+timemat+timefv+timesv+timerho+timepot+timefor+timeio+timemt+timemixer
        if ((verbosity>-1).and.(rank==0)) then
            write(60,*)
            write(60, '(" Wall time (seconds)",T45 ": ", F12.2)') timetot
        end if

! write TOTENERGY.OUT
        if ((verbosity>-1).and.(rank==0)) then
            Write (61, '(G22.12)') engytot
            Call flushifc (61)
        end if

!----------------------
! Convergence tests
!----------------------
        If (iscl .Ge. 2) Then

!...write convergence if only energy
            if ((verbosity>-1).and.(rank==0).and.(input%groundstate%scfconv.eq.'energy')) then
                write(60,*)
                write(60,'(" Absolute change in total energy   (target) : ",G13.6,"  (",G13.6,")")') &
                &     deltae, input%groundstate%epsengy
            end if

!...write convergence if only potential
            if ((verbosity>-1).and.(rank==0).and.(input%groundstate%scfconv.eq.'potential')) then
                if (associated(input%groundstate%OEP)) then
                    write(60,*)
                    write(60, '(" Magnitude of OEP residual",T45 ": ", F18.8)') resoep
                end if
                write(60,*)
                Write(60,'(" RMS change in effective potential (target) : ",G13.6,"  (",G13.6,")")') &
                &     currentconvergence, input%groundstate%epspot
            end if

!...write convergence if only energy
            if ((verbosity>-1).and.(rank==0).and.(input%groundstate%scfconv.eq.'charge')) then
                write(60,*)
                write(60,'(" Charge distance                   (target) : ",G13.6,"  (",G13.6,")")') &
                &     chgdst, input%groundstate%epschg
            end if

!...write convergence if multiple convergence
            if ((verbosity>-1).and.(rank==0).and.(input%groundstate%scfconv.eq.'multiple')) then
                write(60,*)
                Write(60,'(" RMS change in effective potential (target) : ",G13.6,"  (",G13.6,")")') &
               &    currentconvergence, input%groundstate%epspot
                write(60,'(" Absolute change in total energy   (target) : ",G13.6,"  (",G13.6,")")') &
               &    deltae, input%groundstate%epsengy
                write(60,'(" Charge distance                   (target) : ",G13.6,"  (",G13.6,")")') &
               &    chgdst, input%groundstate%epschg
            end if

!...if tforce=true write also convergence of non-IBS forces
            if ((verbosity>-1).and.(rank==0).and.(input%groundstate%tforce)) then
                write(60,'(" Abs. change in max-nonIBS-force   (target) : ",G13.6,"  (",G13.6,")")') &
               &    dforcemax, input%groundstate%epsforcescf
            end if

!...write in RMSDVEFF.OUT and DFSCFMAX.OUT
            if ((verbosity>-2).and.(rank==0)) then
                Write (65, '(G18.10)') currentconvergence
                Call flushifc(65)
                if (input%groundstate%tforce) then
                    Write (67, '(G22.12)') dforcemax
                    Call flushifc(67)
                end if
            end if

!-----------------------
! check for convergence
!-----------------------

            if (input%groundstate%scfconv .eq. 'energy') tlast = all(abs(vdeltae) .lt. input%groundstate%epsengy)

            if (input%groundstate%scfconv .eq. 'potential') tlast = all(abs(vcurrentconvergence) .lt. input%groundstate%epspot)

            if (input%groundstate%scfconv .eq. 'charge') tlast = all(abs(vchgdst) .lt. input%groundstate%epschg)

            if (input%groundstate%scfconv .eq. 'multiple') then
                tlast = all(abs(vcurrentconvergence) .lt. input%groundstate%epspot) .and. &
                &       all(abs(vdeltae) .lt. input%groundstate%epsengy) .and. &
                &       all(abs(vchgdst) .lt. input%groundstate%epschg)
            end if

            if (input%groundstate%tforce) then
                tlast = tlast .and. (dforcemax .lt. input%groundstate%epsforcescf)
                fm = forcemax
            end if

! check for STOP file
            if (rank==0) then
                Inquire (File='STOP', Exist=Exist)
                If (exist) Then
                    write(string,'("STOP file exists - stopping self-consistent loop")')
                    call printbox(60,"+",string)
                    tstop = .True.
                    tlast = .True.
                    Open (50, File='STOP')
                    Close (50, Status='DELETE')
                End If
            end if

        End If ! iscl>2

#ifdef MPI
        Call MPI_bcast (tstop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        Call MPI_bcast (tlast, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
        call timesec(ts1)
        timeio = ts1 - ts0 + timeio
!! TIME - End of fourth IO segment
    End Do ! iscl
! end the self-consistent loop
20  Continue
    call stopwatch("exciting:scf", 0)

    Call timesec(ts0)
    If ((verbosity>-1).and.(rank==0)) Then
        write(string,'("Self-consistent loop stopped")')
        call printbox(60,"+",string)
    end if
    ! write density and potentials to file only if maxscl > 1
    If ((input%groundstate%maxscl > 1)) Then
        Call writestate
        If ((verbosity>-1).and.(rank==0)) Then
            Write (60, '(" STATE.OUT is written")')
        end if
    End If
! delete BROYDEN.OUT
    If (rank==0) then
        Inquire (File='BROYDEN.OUT', Exist=Exist)
        If (exist) Then
            Open (23, File='BROYDEN.OUT')
            Close (23, Status='DELETE')
        End If
    End If
    Call timesec(ts1)
    timeio = ts1 - ts0 + timeio

!------------------
! Compute forces
!------------------
    If (( .Not. tstop) .And. (input%groundstate%tforce)) Then
        Call force(input%groundstate%tfibs)
#ifdef MPI
! For whatever reason each MPI process may produce very slightly different forces.
! At this spot, we equalise them, so that we do not end up with a different geometry
! for every process.
        allocate(forcesum(3,natmtot))
        call MPI_ALLREDUCE(forcetot, forcesum, natmtot*3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        forcetot(1:3,1:natmtot)=forcesum(1:3,1:natmtot)/dble(procs)
        deallocate(forcesum)
#endif
! output forces to INFO.OUT
        if ((verbosity>-1).and.(rank==0)) then
           call printbox(60,"-","Writing atomic positions and forces")
           idm = 0
           write(60,*)
           write(60,'(" Atomic positions (",A,") :")') trim(acoord)
           do is = 1, nspecies
               do ia = 1, natoms (is)
                   idm = idm+1
                   if (input%structure%cartesian) then
                       write(60,'(" atom ",I5,2x,A2,T18,": ",3F14.8)') &
                      &  idm, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
                      &  atposc(:,ia,is)
                   else
                       write(60,'(" atom ",I5,2x,A2,T18,": ",3F14.8)') &
                      &  idm, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
                      &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
                   end if
               end do
           end do
           call writeforce(60,2)
        end if
    End If ! compute forces

    call stopwatch("exciting:post_scf", 1)
    ! set nwork to -2 to tell interface to call the deallocation functions
    If (rank .Eq. 0) Call mixerifc(input%groundstate%mixernumber, n, v, currentconvergence, -2)
    Deallocate(v)
    Call mpiresumeevecfiles()

    if (allocated(rhomtref)) deallocate(rhomtref)
    if (allocated(rhoirref)) deallocate(rhoirref)

    If ((verbosity>-1).and.(rank==0)) Then
! add blank line to TOTENERGY.OUT, FERMIDOS.OUT, MOMENT.OUT and RMSDVEFF.OUT
!      Write (62,*)
      If (associated(input%groundstate%spin)) write (63,*)
! add blank line to DTOTENERGY.OUT, DFORCEMAX.OUT, CHGDIST.OUT and PCHARGE.OUT
!      Write (66,*)
!      If (input%groundstate%tforce) Write (67,*)
!      Write (68,*)
      if (input%groundstate%tpartcharges) write(69,*)
    End If

    If ((verbosity>-2).and.(rank==0)) Then
! write last total energy and add blank line to TOTENERGY.OUT
      Write (61, '(G22.12)') engytot
      Write (61,*)
      Call flushifc(61)
! add blank line to RMSDVEFF.OUT and DFSCFMAX.OUT
      Write (65,*)
      Call flushifc(65)
      if (input%groundstate%tforce) then
          Write (67,*)
          Call flushifc(67)
      end if
    End If

    if (associated(input%groundstate%xsLO).and.(rank==0)) then
      call genxsLOs()
    end if

    call mt_hscf%release()
    call stopwatch("exciting:post_scf", 0)
end subroutine
