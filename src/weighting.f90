subroutine weighting(weightDIR, weightfile, Q_ktemp, Stemp, Nblk, Nrows, nr, ntheta, numvar, Fr)

    character (len = 160)             :: weightDIR, weightfile, filename
    integer   (kind = 4)              :: Nblk, Nrows, nr, ntheta, numvar, i, j
    complex   (kind = 8)              :: Q_ktemp(Nrows, Nblk), Stemp(Nblk, Nblk), Q_kweight(Nrows, Nblk)
    real      (kind = 8), allocatable :: W(:), Fr
     

!!!!!!!!!!!! Weight matrix calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    write(filename, '(a,a)') trim(weightDIR), trim(weightfile)  

    allocate (W(nr*ntheta))

    open(unit=500, file=filename, status = 'old', form = 'formatted', action = 'read')

    do i = 1, nr*ntheta
        read(500, *) W(i)        
    end do

    close(500)  

    do i = 1, numvar

        if (i .lt. 3) then
            do j = 1, nr*ntheta
                Q_kweight((i-1)*nr*ntheta + j, :) = W(j)*Q_ktemp((i-1)*nr*ntheta + j, :)
            end do
        elseif (i .eq. 4) then
            do j = 1, nr*ntheta
                Q_kweight((i-1)*nr*ntheta + j, :) = W(j)*Q_ktemp((i-1)*nr*ntheta + j, :)/(Fr**2)      !! TAPE - Total available
                !Q_kweight((i-1)*nr*ntheta + j, :) =  Q_ktemp((i-1)*nr*ntheta + j, :)                 !! TAPE - Total available
                                                                                                     !! potential energy
                                                                                                     !! Refer Chongsiripinyo and
                                                                                                     !! Sarkar 2019 
            end do
        endif

    end do

    Stemp = matmul(transpose(conjg(Q_ktemp)),Q_kweight) 
    
    return
end subroutine weighting
