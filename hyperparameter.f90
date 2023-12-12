module hyperparameter
    implicit none
    save
    integer           :: size = 120
    integer           :: dim = 2
    real              :: init_temp = 1.5
    real              :: final_temp = 6.5
    real              :: temp_step = 0.04
    integer           :: mcstep = 1000
    integer           :: eqstep = 1000
    integer           :: num_threads = 36
    integer           :: thread_per_row = 6
    integer           :: block_size = 30
    integer           :: interval = 50
    character(len=20) :: directory = "./results/"
    character*1       :: creturn = achar(13)
end module