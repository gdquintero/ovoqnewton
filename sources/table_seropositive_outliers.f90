program table
    implicit none

    integer :: i,inf,sup,n,allocerr,samples
    integer, allocatable :: age(:)
    real*8, allocatable :: measles(:),mumps(:),rubella(:)
    real*8, allocatable :: data(:,:),disease_data(:,:)

    open(UNIT = 1000, FILE = "output/seropositives.txt", ACCESS = "SEQUENTIAL")

    read(1000,*) samples

    allocate(disease_data(samples,5),measles(samples),mumps(samples),rubella(samples),age(samples))

    do i = 1, samples
        read(1000,*) disease_data(i,:)
    enddo

    do i = 1, samples
        age(i) = int(disease_data(i,1))
        measles(i) = disease_data(i,2)
        mumps(i) = disease_data(i,3)
        rubella(i) = disease_data(i,4)
    enddo

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
        data(i,1) = 30 - inf - i
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

    40 format (I5,1X,A1,1X,ES9.3,1X,A1,1X,I5,1X,A1,1X,ES9.3,1X,A1,1X,I5,1X,A1,1X,ES9.3,1X,A1,1X,I4,1X,A2)

end program
