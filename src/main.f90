!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT COMMENTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1. It goes over a set of snapshots to generate 2D SPOD modes
!
! 2. There ghost is removed in the SPOD code itself
!
! 3. Weight matrix and time stamp is provided from an external file
!
! 4. Normalizations are based on the procedure given in 
!    https://github.com/SpectralPOD/spod_matlab
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
  
    use, intrinsic :: iso_c_binding  
    implicit none
    include 'fftw3.f03'  

    !!! Parameters
    real    (kind=8),    parameter ::  pi         = 3.141592653589793238462643383279502884197d0   
    character (len=160), parameter ::  inDIR      = '/mnt/RAID5/sheel/spod_re5e4/fr2/'
    character (len=160), parameter ::  weightDIR  = './'
    character (len=160), parameter ::  weightfile = 'weight_fr2_slice_var_truncated_r_D_10_2d.txt'
    character (len=160), parameter ::  outDIR     = './'
    character (len=160), parameter ::  slice_idx  = '100'
    integer (kind=4)               ::  nr, ntheta, nx, numvar, N, stride, Nfreq, Novlp, idx
    integer (kind=4)               ::  nstart
    real    (kind=8)               ::  dt

    !!! Temporary Variables
    integer (kind=4)               :: i, j, k, iflag, ier, z, info, rst, mode
    real    (kind=8)               :: d, winweight
    real    (kind=8)               :: eps
    character(len=160)             :: fileDIR, filename, basename, folder_name
    integer                        :: namef, start
    real(kind=8),    allocatable   :: Qtemp(:,:), Qtemp1(:,:)
    complex(kind=8), allocatable   :: out1(:), in1(:)
 
    !!! Data Variables
    integer(kind=4)                :: Nrows, Nblk, Nblk_sampled, Nfreq_sampled
    integer(kind=4), allocatable   :: qstart(:), qend(:)                             
    real(kind=8),    allocatable   :: Q(:,:), Q_sub(:,:), Q_mean(:), freq(:), Qblk(:,:,:), P(:,:,:,:)
    real(kind=8),    allocatable   :: Pr(:,:,:), P_trunc(:,:,:), window(:)
    complex(kind=8), allocatable   :: Qout(:,:,:), Q_k(:,:,:), Eigen_V(:,:,:)
    complex(kind=8), allocatable   :: Lambda(:,:,:), Lambda_invsq(:,:,:), Stemp(:,:)
                                
    !!! FFT Variables
    !complex (kind=8), allocatable :: Q_inv(:,:,:), inv1(:,:)
    complex                        :: alpha, alpha1      
    integer (kind=8)               :: plan, plan_forward, plan_backward
   
    !!! LAPACK Variables
    integer                        :: lwork
    real(kind=8),    allocatable   :: rwork(:)
    complex(kind=8), allocatable   :: work(:)
    complex(kind=8), allocatable   :: vl(:,:), vr(:,:), S_Nfreq(:,:), S(:,:,:)                             
    real(kind=8), allocatable      :: L_Nfreq(:)   ! Using zheev
    !real(kind=8), allocatable  ::   L_Nfreq(:)   ! Using zgeev

    !!!!!!!!!!!!!! Reading the input data  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    open(unit=110,file='spod_parameters.input',form="formatted")

    read (110,*) nr
    read (110,*) ntheta
    read (110,*) nx
    read (110,*) numvar
    read (110,*) dt
    read (110,*) N
    read (110,*) nstart
    read (110,*) stride
    read (110,*) Nfreq
    read (110,*) Nblk_sampled
    read (110,*) Nfreq_sampled
    read (110,*) Novlp

    close(110)

    print*,  nr
    print*,  ntheta
    print*,  nx
    print*,  numvar
    print*,  dt
    print*,  N
    print*,  nstart
    print*,  stride
    print*,  Nblk_sampled
    print*,  Nfreq_sampled
    print*,  Nfreq
    print*,  Novlp

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Nrows = nr*ntheta*nx*numvar              ! No. of rows in the snapshot matrix 
   
    !!!!!!!!!!!!! Reading data files in an array !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    allocate (Q(Nrows,N))                    ! Data snpshots arrandged columnwise
    allocate (Q_mean(Nrows))                 ! Mean of data snapshots
    allocate (Q_sub(Nrows,N))                ! Data snapshots after ssubtraction of mean
    allocate (Pr(nr+2,ntheta+2,nx))             ! real part of the data 
    allocate (P_trunc(nr,ntheta,nx))            ! real part of the data 
    allocate (P(nr,ntheta,nx,numvar))     
 
    do rst = 1, N
    
    
        namef = nstart + (rst-1)*stride
        print*, 'Reading the file ', namef
    
        !!!!! Reading the up files
        folder_name = 'up_slice'
        basename = 'up'
        write(filename,'(a,a,a,a,a,a,a,a,i8.8,a,a,a)') trim(inDIR), trim(folder_name), "/", "x_D_", trim(slice_idx), "/", trim(basename), "_", namef, "_", &
                                                 trim(slice_idx), ".5p"
        print*, filename, nr, ntheta, nx
        call io_slice_files(filename, nr, ntheta, nx, Pr)


        do k = 1, nr
             P_trunc(k, 1:ntheta, :) = 0.50d0*(Pr(k+1, 2:ntheta+1, :) + Pr(k, 2:ntheta+1, :))  !!!!!! Centered the u velocity field !!!!!!!!!!!!!!!!!!!!!!!!
        end do

        P(1:nr,1:ntheta,:,1) = P_trunc
                
        !!!!! Reading the vp files
        folder_name = 'vp_slice'
        basename = 'vp'
        write(filename,'(a,a,a,a,a,a,a,a,i8.8,a,a,a)') trim(inDIR), trim(folder_name), "/", "x_D_", trim(slice_idx), "/", trim(basename), "_", namef, "_", &
                                                 trim(slice_idx), ".5p"
        call io_slice_files(filename, nr, ntheta, nx, Pr)

         do k = 1, ntheta
              P_trunc(1:nr, k, :) = 0.50d0*(pr(2:nr+1, k, :) + Pr(2:nr+1, k+1, :))          !!!!!! Centered the v velocity field !!!!!!!!!!!!!!!!!!!!!!!!
         end do
         
        P(1:nr,1:ntheta,:,2) = P_trunc

        !!!!! Reading the wp files
        folder_name = 'wp_slice'
        basename = 'wp'
        write(filename,'(a,a,a,a,a,a,a,a,i8.8,a,a,a)') trim(inDIR), trim(folder_name), "/",  "x_D_", trim(slice_idx), "/", trim(basename), "_", namef, "_", &
                                                 trim(slice_idx), ".5p"
        call io_slice_files(filename, nr, ntheta, nx, Pr)
        P(1:nr,1:ntheta,:,3) = Pr(2:nr+1, 2:ntheta+1, :)                             !!!!!! Centered the w velocity field !!!!!!!!!!!!!!!!!!!!!!!!

        !!!!! Reading the densp files
        folder_name = 'densp_slice'
        basename = 'densp'
        write(filename,'(a,a,a,a,a,a,a,a,i8.8,a,a,a)') trim(inDIR), trim(folder_name), "/", "x_D_", trim(slice_idx), "/", trim(basename), "_", namef, "_", &
                                                 trim(slice_idx), ".5p"
        call io_slice_files(filename, nr, ntheta, nx, Pr)
        P(1:nr,1:ntheta,:,4) = Pr(2:nr+1, 2:ntheta+1, :)                              !!!!!! Centered the dens field !!!!!!!!!!!!!!!!!!!!!!!!


        !!! Arranging all the data in Q - snapshot matrix
        do i=1,numvar
            do k=1,nx
                do j=1,ntheta
                    Q((1 + nr*(j-1) + nr*ntheta*(k-1) + nr*ntheta*nx*(i-1)):(nr*j + nr*ntheta*(k-1) + nr*ntheta*nx*(i-1)),rst) = P(1:nr,j,k,i)
                end do
            end do
        end do

    end do

    deallocate(Pr, P)

    print*, 'Shape of Q ',                     shape(Q)       
    print*, 'minval of real part of Q ',       minval((Q(1:nr*ntheta,:)))     
    print*, 'maxval of real part of Q ',       maxval((Q(1:nr*ntheta,:)))     

    !!!!!!!!!!! Subtracting the mean from snapshot matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call subtract_qmean(Q, Nrows, N, Q_mean, Q_sub)      
    
    deallocate(Q_mean,Q)
    
    print*, 'Rowwise mean subtracted'
    print*, 'minval of real part of  Q_sub ',      minval((Q_sub(:,:)))     
    print*, 'maxval of real part of  Q_sub ',      maxval((Q_sub(:,:)))     
      
    !!!!!!!!!! Dividing data into blocks !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    if (Novlp > Nfreq) then 
        print*, "Error: Overlap too large"
    end if 
 
    !!!!!!!!! Calculating the number of blocks !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    d = ((N-Novlp)/(Nfreq-Novlp))
    Nblk = floor(d)       
    print*, 'Number of blocks ', Nblk
    
    allocate (qstart(Nblk))
    allocate (qend(Nblk))
       
    !!!!!!!! Start and End of Each Block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, Nblk
        qstart(i) = (i-1)*(Nfreq-Novlp)+1
        qend(i)   = qstart(i)+Nfreq-1
    end do

    allocate (Qtemp(Nrows, Nfreq))
    allocate (in1(Nfreq))
    allocate (out1(Nfreq))
    print*, '2. Allocating the most memory consuming array Qblk'
    allocate (Qblk(Nrows, Nfreq, Nblk))
    allocate (window(Nfreq))

    !allocate (Qtemp1(Nfreq,Nrows))
    !allocate (Q_inv(Nfreq,Nrows,Nblk)) 
    !allocate (inv1(Nfreq,Nrows))
    
    !!!!!!!!! Seperating Q into blocks !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, Nblk
        Qblk(:,1:Nfreq,i) = Q_sub(:,qstart(i):qend(i))  
    end do    

    deallocate(Q_sub)

    !!!!!!!! Calling window function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call hamming_window(Nfreq, window)
    print*, 'Constructed the hamming window'    
    print*, 'Minval of window', minval(window(:))
    print*, 'Maxval of window', maxval(window(:))
    
    !!!!!!! Windowing the Q block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    do j = 1, Nblk
        do i = 1, Nfreq
            Qblk(:,i,j) = Qblk(:,i,j)*window(i)
        end do
    end do      
    
    print*, 'Windowed Qblk'    

    !!!!!!! Normalizing the Q block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !winweight = Nfreq/(sum(window))           !!! This is done to preserve the amplitude of windowed signal, doesn't preserve energy
    winweight = 1.587845072201937             !!! This factor for hamming window preserves the energy of windowed signal
    Qblk = (winweight/Nfreq)*Qblk   
    
    print*, 'minval of real Qblk ', minval((Qblk(:,:,:)))     
    print*, 'maxval of real Qblk ', maxval((Qblk(:,:,:)))     

    deallocate(window)  
    
    print*, '1. Allocating the most memory consuming array Qout'
    allocate (Qout(Nrows, Nfreq, Nblk))
 
    !!!!!!!!!!!!!!!! Row wise FFT of individual blocks !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1,Nblk
        print*, 'In Nblk ', i
        Qtemp(1:Nrows,1:Nfreq) = Qblk(1:Nrows,1:Nfreq,i)
        
        do j = 1, Nrows
            in1 = Qtemp(j,:)    
            call dfftw_plan_dft_1d(plan_forward, Nfreq, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE)
            call dfftw_execute(plan_forward, in1, out1)
            call dfftw_destroy_plan(plan_forward)
            Qout(j,:,i) = out1
        end do
        !!!! Backward transform to verify
        !call dfftw_plan_dft_c2r_2d(plan_backward,Nfreq,Nrows,out1,inv1,FFTW_ESTIMATE)
        !call dfftw_execute(plan_backward)
        !call dfftw_destroy_plan(plan_backward) 
        !   Q_inv(:,:,i) = inv1
    end do
   
    
    deallocate (Qtemp, Qblk)       
    deallocate (in1, out1)
    deallocate (qstart, qend) 
    allocate   (Q_k(Nrows,Nblk,Nfreq)) 

    !!!!!!!!!!!! Sorting the FFT matrix at each frequency k !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    do i = 1, Nblk
        do j = 1, Nfreq
            Q_k(:,i,j) = Qout(:,j,i)
        end do
    end do 

    deallocate(Qout)

    print*, 'minval of imag Q_k ', minval(aimag(Q_k(:,:,:)))      
    print*, 'maxval of imag Q_k ', maxval(aimag(Q_k(:,:,:)))      
    print*, 'minval of real Q_k ', minval(real(Q_k(:,:,:)))     
    print*, 'maxval of real Q_k ', maxval(real(Q_k(:,:,:)))     

    !!!!!!!!!!!!!!!! BLOCK FOR DEBUGGING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! filename = 'Q_kreal_1.txt'
    ! open(unit=500, file=filename,  status='replace', & 
    !      form='formatted', access='stream', action='write')

    !  do j = 1,Nblk
    !      do i = 1, Nrows
    !          write(500,*) real(Q_k(i,j,1))
    !      enddo
    !  enddo

    !  close(500)
    
    ! filename = 'Q_kimag_1.txt'
    ! open(unit=500, file=filename,  status='replace', & 
    !      form='formatted', access='stream', action='write')

    !  do j = 1,Nblk
    !      do i = 1,Nrows
    !          write(500,*) aimag(Q_k(i,j,1))
    !      enddo
    !  enddo

    !  close(500)
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!! Calculating the Cross specral Density !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate (S(Nblk,Nblk,Nfreq))
    allocate (Stemp(Nblk, Nblk))
   
    do i = 1, Nfreq
        call weighting(weightDIR, weightfile, Q_k(:,:,i), Stemp, Nblk, Nrows, nr, ntheta, numvar)
        S(:,:,i) = Stemp
    end do

    !!!!!!! Normalizing the Cross Spectral Density !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
    S = S/(Nblk)                             ! For ensuring the weighted cross spectral density that is formed
    
    print*, 'minval of imaginary S ', minval(aimag(S(:,:,:)))      
    print*, 'maxval of imaginary S ', maxval(aimag(S(:,:,:)))      
    print*, 'minval of real S ',      minval(real(S(:,:,:)))     
    print*, 'maxval of real S ',      maxval(real(S(:,:,:)))     
        
    !!!!!!! Calculating the eigenvalues and eigen vectors !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate (Lambda(Nblk,Nblk,Nfreq))
    allocate (Lambda_invsq(Nblk,Nblk,Nfreq))
    allocate (L_Nfreq(Nblk))
    
    Lambda = 0.0d0
    Lambda_invsq = 0.d0

    !! LAPACK Parameters for  eigen values and eigen vectors 
    lwork = 2*Nblk-1              !zheev
    !lwork = 2*Nblk               !zgeev
    !allocate (rwork(2*Nblk))     !zgeev
    allocate (rwork(3*Nblk-2))    !zheev
    allocate (work(lwork))
    allocate (S_Nfreq(Nblk,Nblk))
    allocate (Eigen_V(Nblk,Nblk,Nfreq))
    !allocate (vr(Nblk,Nblk)) ! zheev
    !allocate (vl(Nblk,Nblk)) ! zgeev

    !!!!!!! Performing the Eigenvalue decomposition of CSD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, Nfreq
        
        print*, 'In Nfreq ', i
        
        S_Nfreq(1:Nblk,1:Nblk) = S(1:Nblk,1:Nblk,i) 
                    
        print*, 'minval of imaginary S_Nfreq ', minval(aimag(S_Nfreq(:,:)))      
        print*, 'maxval of imaginary S_Nfreq ', maxval(aimag(S_Nfreq(:,:)))      
        print*, 'minval of real S_Nfreq ',      minval(real(S_Nfreq(:,:)))     
        print*, 'maxval of real S_Nfreq ',      maxval(real(S_Nfreq(:,:)))     
        print*, 'Shape of S_Nfreq ',            shape(S_Nfreq)

        call zheev('V','L', Nblk, S_Nfreq, Nblk, L_Nfreq, work, lwork, rwork, info)    
        !call zgeev('N','V', Nblk, S_Nfreq, Nblk, L_Nfreq, vl, 1, vr, Nblk, work, lwork, rwork, info)
    
        do j = 1, Nblk
            Lambda(j,j,i)       = L_Nfreq(j)
            Lambda_invsq(j,j,i) = 1/sqrt(L_Nfreq(j)) 
        end do

        Eigen_V(:,:,i) = S_Nfreq(:,:)

    end do

    deallocate(S, S_Nfreq)
   
    !!!!!!!!!!!!!!! Checking the eigenvectors !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !filename = 'eigenvector_hf_real.txt'
    !open(unit=500, file=filename,  status='replace', form='formatted', access='stream', action='write')

    !do j = 1,Nblk
    !    do i = 1,Nblk
    !        write(500,*) real(Eigen_V(i,j,hf))
    !    enddo
    !enddo

    !close(500)
    !     
    !filename = 'eigenvector_hf_imag.txt'
    !open(unit=500, file=filename,  status='replace', form='formatted', access='stream', action='write')

    !do j = 1,Nblk
    !    do i = 1,Nblk
    !        write(500,*) aimag(Eigen_V(i,j,hf))
    !    enddo
    !enddo

    !close(500)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!! Writing out eigenvalues, computing eigenmodes and writing them to file !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call io_spod_files(Lambda, Lambda_invsq, Eigen_V, Q_k, outDIR, Nrows, Nblk, Nblk_sampled, Nfreq, Nfreq_sampled, mode)
 
    deallocate (Lambda, Lambda_invsq, L_Nfreq)
    deallocate (Q_k)
    deallocate (rwork)     !zgeev
    deallocate (work)
    deallocate (Eigen_V)
    deallocate (Stemp)

    
end program main
