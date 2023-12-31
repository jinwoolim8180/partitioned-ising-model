program main
    use OMP_LIB
    use f90getopt
    use hyperparameter
    use montecarlo
    use, intrinsic :: iso_fortran_env
    implicit none

    ! hyperparameters entered as arguments
    ! ------------------------------------
    type(option_s) :: opts(11)

    ! statistical values and spin configuration Monte Carlo
    !--------------------------------------------------
    real, dimension(:), allocatable       :: Ts, Es, Ms, Cs, Xs
    integer, dimension(:, :), allocatable :: spins
    real(real64)                          :: tempE, tempM
    real(real64)                          :: T, beta, E1, M1, E2, M2
    real(real64)                          :: norm, volume

    ! misc settings
    ! -------------
    integer            :: i, j, threadx, thready
    integer            :: xpp, ypp, xnn, ynn
    integer            :: step, num_temps, temp_index
    integer            :: start, end, rate
    real(real64)       :: progress
    character(len=7)   :: size_char, dim_char, eqstep_char, mcstep_char
    character(len=100) :: csv_file

    ! receive argument from the command line
    ! --------------------------------------
    opts(1) = option_s('block_size', .true.,  's')
    opts(2) = option_s('dim',  .true., 'd')
    opts(3) = option_s('init_temp', .true., 'i')
    opts(4) = option_s('final_temp', .true., 'f')
    opts(5) = option_s('temp_step', .true., 't')
    opts(6) = option_s('mcstep', .true., 'm')
    opts(7) = option_s('eqstep', .true., 'e')
    opts(8) = option_s('thread_per_row', .true., 'p')
    opts(9) = option_s('interval', .true., 'v')
    opts(10) = option_s('dir', .true., 'r')
    opts(11) = option_s('help', .false., 'h')

    do 
        select case (getopt("s:d:i:f:t:m:e:p:v:r:h", opts))
            case ('s')
                read (optarg, '(i4)') block_size
            case ('d')
                read (optarg, '(i1)') dim
            case ('i')
                read (optarg, '(f7.5)') init_temp
            case ('f')
                read (optarg, '(f7.5)') final_temp
            case ('t')
                read (optarg, '(f7.6)') temp_step
            case ('m')
                read (optarg, '(i7)') mcstep
            case ('e')
                read (optarg, '(i7)') eqstep
            case ('p')
                read (optarg, '(i3)') thread_per_row
            case ('v')
                read (optarg, '(i4)') interval
            case ('r')
                read (optarg, '(a)') directory
            case ('h')
                call print_help()
                stop
            case default
                exit
        end select
    end do

    ! print number of threads for simulation
    ! --------------------------------------
    num_threads = thread_per_row ** 2
    size = block_size * thread_per_row
    num_temps = int((final_temp - init_temp) / temp_step) + 1
    if (num_threads > get_num_threads()) then
        print '(a)', "Thread setting not available for your CPU"
        stop
    endif

    ! allocate array to store simulation results
    ! ------------------------------------------
    allocate(Ts(num_temps))
    allocate(Es(num_temps))
    allocate(Ms(num_temps))
    allocate(Cs(num_temps))
    allocate(Xs(num_temps))
    Ts = 0
    Es = 0
    Ms = 0
    Cs = 0
    Xs = 0

    ! parallel Monte Carlo simulation
    ! -------------------------------
    temp_index = 0
    progress = 0.
    volume = size ** dim
    norm = mcstep * volume
    progress = 0.
    allocate(spins(size, size))
    101 format(a, "Progress: ", f5.1, '% / ', f5.1, '% (Time: ', f12.6, 's)')
    call system_clock(start, rate)
    do temp_index = 1, num_temps
        ! initialize spin configuration, etc.
        E1 = 0; M1 = 0; E2 = 0; M2 = 0;
        T = init_temp + temp_step * (temp_index - 1)
        beta = 1.0 / T
        spins = 1

        !$OMP PARALLEL PRIVATE(i, j, threadx, thready, xpp, ypp, xnn, ynn, end) &
        !$OMP SHARED(spins, progress) num_threads(num_threads)
        i = int(omp_get_thread_num() / thread_per_row)
        j = omp_get_thread_num() - thread_per_row * i
        threadx = 1 + block_size * i
        thready = 1 + block_size * j
        call get_neighbors_index(i, j, xpp, ypp, xnn, ynn)

        mcmc : block
            integer, dimension(:, :), allocatable :: spin
            integer, dimension(:, :), allocatable :: neighbors
            real, dimension(:, :), allocatable    :: randoms
            real random
            allocate(spin(block_size, block_size))
            allocate(neighbors(2 * dim, block_size))
            allocate(randoms(block_size, block_size))
            call random_number(randoms)
            spin = 2 * int(2 * randoms) - 1
            neighbors = 1

            ! equilibration steps
            do step = 1, eqstep
                call metropolis(spin, neighbors, beta)

                ! show progress bar
                call system_clock(end)
                !$OMP ATOMIC
                progress = progress + 100. / (num_temps * num_threads * (eqstep + mcstep))
                write (*, 101, advance='no') creturn, progress, 100.0, real(end - start) / real(rate)

                spins(threadx:threadx + block_size - 1, thready:thready + block_size - 1) = spin
                call random_number(random)
                if (random < 1.0 / interval) then
                    neighbors(1, :) = spins(block_size * (xnn + 1), thready:thready + block_size - 1)
                    neighbors(2, :) = spins(threadx:threadx + block_size - 1, block_size * (ynn + 1))
                    neighbors(3, :) = spins(1 + block_size * xpp, thready:thready + block_size - 1)
                    neighbors(4, :) = spins(threadx:threadx + block_size - 1, 1 + block_size * ypp)
                end if
            end do

            ! Monte Carlo steps
            do step = 1, mcstep
                call metropolis(spin, neighbors, beta)

                ! show progress bar
                call system_clock(end)
                !$OMP ATOMIC
                progress = progress + 100. / (num_temps * num_threads * (eqstep + mcstep))
                write (*, 101, advance='no') creturn, progress, 100.0, real(end - start) / real(rate)
                spins(threadx:threadx + block_size - 1, thready:thready + block_size - 1) = spin
                call random_number(random)
                if (random < 1.0 / interval) then
                    neighbors(1, :) = spins(block_size * (xnn + 1), thready:thready + block_size - 1)
                    neighbors(2, :) = spins(threadx:threadx + block_size - 1, block_size * (ynn + 1))
                    neighbors(3, :) = spins(1 + block_size * xpp, thready:thready + block_size - 1)
                    neighbors(4, :) = spins(threadx:threadx + block_size - 1, 1 + block_size * ypp)
                end if

                ! calculate physical quantities
                ! $OMP BARRIER
                !$OMP SINGLE
                tempE = calc_energy(spins)
                tempM = calc_magnetization(spins)
                E1 = E1 + tempE / norm
                M1 = M1 + tempM / norm
                E2 = E2 + tempE ** 2 / norm
                M2 = M2 + tempM ** 2 / norm
                !$OMP END SINGLE
            end do
            deallocate(spin)
            deallocate(neighbors)
            deallocate(randoms)
        end block mcmc
        !$OMP END PARALLEL

        ! save results to array
        Ts(temp_index) = T
        Es(temp_index) = E1
        Ms(temp_index) = M1
        Cs(temp_index) = (E2 - volume * E1 ** 2) * beta ** 2
        Xs(temp_index) = (M2 - volume * M1 ** 2) * beta
    end do
    call system_clock(end)
    print *, ""
    print '("Time = ",f12.6," seconds.")', real(end - start) / real(rate)

    ! deallocate allocated memory
    ! ---------------------------
    deallocate(spins)
    print *, ""

    ! save results to the directory
    ! -----------------------------
    write (size_char, '(i4.4)') size
    write (dim_char, '(i1)') dim
    write (eqstep_char, '(i7.7)') eqstep
    write (mcstep_char, '(i7.7)') mcstep
    csv_file = trim(directory) // trim("result_L") // trim(size_char) &
                               // trim("_D") // trim(dim_char) &
                               // trim("_EQ") // trim(eqstep_char) &
                               // trim("_MC") // trim(mcstep_char) // trim(".csv")
    open(unit=1, file=csv_file, status='unknown')
    write (1, '(5(a, x))') 'T', 'E', 'M', 'C', 'X'
    do temp_index = 1, num_temps
        write (1, '(5(g20.10, x))') Ts(temp_index), Es(temp_index), Ms(temp_index), Cs(temp_index), Xs(temp_index)
    end do
    close(1)
    print '(2a, /)', "Results saved to ", csv_file

    deallocate(Ts)
    deallocate(Es)
    deallocate(Ms)
    deallocate(Cs)
    deallocate(Xs)

contains
    subroutine print_help()
        print '(a, /)', 'Command-line options:'
        print '(a)', '  -s, --block_size      size of the partitioned block of a lattice (default: 30)'
        print '(a)', '  -d, --dim             dimension of the lattice                   (default: 3)'
        print '(a)', '  -i, --init_temp       initial temperature of the output          (default: 1.5)'
        print '(a)', '  -f, --final_temp      final temperature of the output            (default: 6.5)'
        print '(a)', '  -t, --temp_step       step size of the temperature               (default: 0.04)'
        print '(a)', '  -m, --mcstep          number of Monte Carlo steps                (default: 1000)'
        print '(a)', '  -e, --eqstep          number of steps for equilibration          (default: 1000)'
        print '(a)', '  -p, --thread_per_row  number of threads per row of a lattice     (default: 1000)'
        print '(a)', '  -r, --dir             directory to save the results              (default: ./results/)'
        print '(a, /)', '  -h, --help            print usage information and exit'
    end subroutine print_help

    integer function get_num_threads()
        integer num_threads
        !$OMP PARALLEL SHARED(num_threads)
        num_threads = omp_get_num_threads()
        !$OMP END PARALLEL
        print '(a, i3)', "Number of threads: ", num_threads
        get_num_threads = num_threads
    end function get_num_threads

    subroutine get_neighbors_index(x, y, xpp, ypp, xnn, ynn)
        integer, intent(in)  :: x, y
        integer, intent(out) :: xpp, ypp, xnn, ynn

        ! obtain neighbor threads
        if (x == 0) then
            xpp = 1
            xnn = thread_per_row - 1
        elseif (x == thread_per_row - 1) then
            xpp = 0; ypp = y
            xnn = thread_per_row - 2
        else
            xpp = x + 1
            xnn = x - 1
        endif

        if (y == 0) then
            ypp = 1
            ynn = thread_per_row - 1
        elseif (y == thread_per_row - 1) then
            ypp = 0
            ynn = thread_per_row - 2
        else
            ypp = y + 1
            ynn = y - 1
        endif
    end subroutine get_neighbors_index

    subroutine get_neighbors(spins, neighbors, threadx, thready, xpp, ypp, xnn, ynn)
        integer, dimension(:, :), allocatable, intent(inout) :: spins
        integer, dimension(:, :), allocatable, intent(inout) :: neighbors
        integer, intent(in)                                  :: threadx, thready
        integer, intent(in)                                  :: xpp, ypp, xnn, ynn
        ! print '(i3, i3, i3, i3)', xpp, ypp, xnn, ynn
        ! print '(i5, i5, i4, i4)', block_size * (xnn + 1), block_size * (ynn + 1), 1 + block_size * xpp, 1 + block_size * ypp
        neighbors(1, :) = spins(block_size * (xnn + 1), thready:thready + block_size - 1)
        neighbors(2, :) = spins(threadx:threadx + block_size - 1, block_size * (ynn + 1))
        neighbors(3, :) = spins(1 + block_size * xpp, thready:thready + block_size - 1)
        neighbors(4, :) = spins(threadx:threadx + block_size - 1, 1 + block_size * ypp)
    end subroutine get_neighbors

end program main