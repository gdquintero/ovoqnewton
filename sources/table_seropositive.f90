program table
    implicit none

    integer :: i,inf,sup,n,allocerr
    integer, dimension(29) :: age
    real*8, dimension(29) :: measles,mumps,rubella
    real*8, allocatable :: data(:,:)

    measles(:) = (/&
        0.207d0,0.301d0,0.409d0,0.589d0,0.757d0,0.669d0,0.797d0,0.818d0,0.866d0,0.859d0,&
        0.908d0,0.923d0,0.889d0,0.936d0,0.889d0,0.898d0,0.959d0,0.957d0,0.937d0,0.918d0,&
        0.939d0,0.967d0,0.973d0,0.943d0,0.967d0,0.946d0,0.961d0,0.968d0,0.968d0         &
    /)

    mumps(:) = (/&
        0.115d0,0.147d0,0.389d0,0.516d0,0.669d0,0.768d0,0.786d0,0.798d0,0.878d0,0.861d0,&
        0.844d0,0.881d0,0.895d0,0.882d0,0.869d0,0.895d0,0.911d0,0.920d0,0.915d0,0.950d0,&
        0.909d0,0.873d0,0.880d0,0.915d0,0.906d0,0.933d0,0.917d0,0.898d0,0.839d0         &
    /)
    
    rubella(:) = (/&
        0.126d0,0.171d0,0.184d0,0.286d0,0.400d0,0.503d0,0.524d0,0.634d0,0.742d0,0.664d0,&
        0.735d0,0.815d0,0.768d0,0.842d0,0.760d0,0.869d0,0.844d0,0.852d0,0.907d0,0.935d0,&
        0.921d0,0.896d0,0.890d0,0.949d0,0.899d0,0.955d0,0.937d0,0.933d0,0.917d0         &
    /)

    age(:) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,19,21,23,25,27,29,31,33,35,40,45,55,65/)

    open(UNIT = 100, FILE = "output/latex_sero.txt", ACCESS = "SEQUENTIAL")
    open(UNIT = 110, FILE = "output/latex_mixed_test.txt", ACCESS = "SEQUENTIAL")
    open(UNIT = 200, FILE = "output/fobj_mixed_measles.txt", ACCESS = "SEQUENTIAL")
    open(UNIT = 210, FILE = "output/fobj_mixed_mumps.txt", ACCESS = "SEQUENTIAL")
    open(UNIT = 220, FILE = "output/fobj_mixed_rubella.txt", ACCESS = "SEQUENTIAL")
    open(UNIT = 400, FILE = "output/iterations_mixed_measles.txt", ACCESS = "SEQUENTIAL")
    open(UNIT = 410, FILE = "output/iterations_mixed_mumps.txt", ACCESS = "SEQUENTIAL")
    open(UNIT = 420, FILE = "output/iterations_mixed_rubella.txt", ACCESS = "SEQUENTIAL")
    open(UNIT = 500, FILE = "output/num_mixed_test.txt", ACCESS = "SEQUENTIAL")

    do i = 1,14
        if (i .ne. 14) then
            write(100,10) '$[',age(i),',',age(i+1),')$','&',measles(i),'&',mumps(i),'&',rubella(i),'&',&
                          '$[',age(i+15),',',age(i+16),')$','&',measles(i+15),'&',mumps(i+15),'&',rubella(i+15),'\\'
        else
            write(100,20) '$[',age(i),',',age(i+1),')$','&',measles(i),'&',mumps(i),'&',rubella(i),'&',&
                          '$[',age(29),',+\infty)$','&',measles(29),'&',mumps(29),'&',rubella(29),'\\'
        end if
    end do

    write(100,30) '$[',age(15),',',age(16),')$','&',measles(15),'&',mumps(15),'&',rubella(15),'&','&','&','\\'

    10 format (A2,I2,A1,I2,A2,1X,A1,1X,F5.3,1X,A1,1X,F5.3,1X,A1,1X,F5.3,1X,A1,1X,&
               A2,I2,A1,I2,A2,1X,A1,1X,F5.3,1X,A1,1X,F5.3,1X,A1,1X,F5.3,1X,A2)

    20 format (A2,I2,A1,I2,A2,1X,A1,1X,F5.3,1X,A1,1X,F5.3,1X,A1,1X,F5.3,1X,A1,1X&
               A2,I2,A10,1X,A1,1X,F5.3,1X,A1,1X,F5.3,1X,A1,1X,F5.3,1X,A2)

    30 format (A2,I2,A1,I2,A2,1X,A1,1X,F5.3,1X,A1,1X,F5.3,1X,A1,1X,F5.3,1X,A1,1X,A1,1X,A1,1X,A2)

    read(500,*) inf
    read(500,*) sup

    n = sup - inf + 1

    allocate(data(n,7),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
    end if

    do i = 1, n
        data(i,1) = 29 - i
        read(200,*) data(i,2)
        read(400,*) data(i,3)
        read(210,*) data(i,4)
        read(410,*) data(i,5)
        read(220,*) data(i,6)
        read(420,*) data(i,7)
    enddo

    do i = 1, n
        write(110,40) int(data(i,1)),'&',data(i,2),'&',int(data(i,3)),&
                                     '&',data(i,4),'&',int(data(i,5)),&
                                     '&',data(i,6),'&',int(data(i,7)),'\\'
    enddo

    40 format (I4,1X,A1,1X,ES9.3,1X,A1,1X,I4,1X,A1,1X,ES9.3,1X,A1,1X,I4,1X,A1,1X,ES9.3,1X,A1,1X,I4,1X,A2)

end program
