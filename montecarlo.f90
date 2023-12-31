module montecarlo
    use hyperparameter
    use, intrinsic :: iso_fortran_env
    implicit none
    PUBLIC :: init_spin, metropolis, calc_energy, calc_magnetization
    PRIVATE

contains
    ! allocate memory to spin and initialize all spins to 1
    ! -----------------------------------------------------
    subroutine init_spin(spin)
        integer, dimension(:, :), allocatable :: spin
        real, dimension(:, :), allocatable    :: rand
        ! allocate(spin(2 * dim, block_size))
        allocate(rand(2 * dim, block_size))
        call random_number(rand)
        spin = 2 * int(2 * rand) - 1
    end subroutine init_spin

    ! Metropolis-Hastings algorithm for Monte Carlo simulation
    ! --------------------------------------------------------
    subroutine metropolis(spin, neighbors, beta)
        ! argument of the subroutine
        integer, dimension(:, :), allocatable :: spin
        integer, dimension(:, :), allocatable :: neighbors
        real(real64), intent(in)              :: beta
        ! local variables of the subroutine
        integer, dimension(:, :), allocatable :: next_neighbors
        integer i, x, y, s, R
        integer xpp, ypp, xnn, ynn ! neighbor spins
        real x0, y0, rand
        real dH

        ! call init_spin(neighbors)
        ! allocate(next_neighbors(4, block_size))
        ! neighbors(1, :) = spin(block_size, :)
        ! neighbors(2, :) = spin(:, block_size)
        ! neighbors(3, :) = spin(1, :)
        ! neighbors(4, :) = spin(:, 1)
        next_neighbors = neighbors

        ! Metropolis-Hastings algorithm main body
        neighbors = 1
        do i = 1, block_size ** dim
            call random_number(x0)
            call random_number(y0)
            x = int(1 + (block_size - 1) * x0) ! x coordinate of selected spin
            y = int(1 + (block_size - 1) * y0) ! y coordinate of selected spin
            s = spin(x, y)               ! randomly selected spin

            ! obtain neighbor spins
            if (x == 1) then
                xpp = spin(2, y)
                xnn = neighbors(1, y)
            elseif (x == block_size) then
                xpp = neighbors(3, y)
                xnn = spin(block_size - 1, y)
            else
                xpp = spin(x + 1, y)
                xnn = spin(x - 1, y)
            endif

            if (y == 1) then
                ypp = spin(x, 2)
                ynn = neighbors(2, x)
            elseif (y == block_size) then
                ypp = neighbors(4, x)
                ynn = spin(x, block_size - 1)
            else
                ypp = spin(x, y + 1)
                ynn = spin(x, y - 1)
            endif

            ! flip the spin and decide the next state
            R = xpp + ypp + xnn + ynn
            dH = 2 * s * R
            call random_number(rand)
            if (dH < 0.) then
                spin(x, y) = -s
            elseif (rand < exp(-beta * dH)) then
                spin(x, y) = -s
            endif

            ! save neighbors to the matrix
            call random_number(rand)
            if (rand < 1 / interval) then
                neighbors = next_neighbors
                next_neighbors(1, :) = spin(block_size, :)
                next_neighbors(2, :) = spin(:, block_size)
                next_neighbors(3, :) = spin(1, :)
                next_neighbors(4, :) = spin(:, 1)
            endif
        end do
    end subroutine metropolis

    ! calculation of energy of the spin using cshift
    ! ------------------------------
    real(real64) function calc_energy(spin)
        integer, dimension(:, :), allocatable, intent(in) :: spin
        integer, dimension(:, :), allocatable :: xps
        integer, dimension(:, :), allocatable :: yps
        integer, dimension(:, :), allocatable :: xns
        integer, dimension(:, :), allocatable :: yns

        ! shifted configuration of spins for easy computation
        xps = cshift(spin, shift=1, dim=1)
        yps = cshift(spin, shift=1, dim=2)
        xns = cshift(spin, shift=-1, dim=1)
        yns = cshift(spin, shift=-1, dim=2)

        calc_energy = sum(spin * -(xps + yps + xns + yns)) / (2 * dim)
    end function calc_energy

    ! calculation of magnetization of spin
    ! ------------------------------------
    real(real64) function calc_magnetization(spin)
        integer, dimension(:, :), allocatable, intent(in) :: spin
        calc_magnetization = sum(spin)
    end function calc_magnetization

end module