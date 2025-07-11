module hipFortMap

   ! A BAND-AID MODULE TO MAP HIP C INTERFACES TO FORTRAN
   ! WITH ISO_C_BINDING
   ! ESSENTIALLY THE HIP IS A IMATATION OF THE CUDA INTERFACE BY THE AMD
   ! COMPANY. IT DOES HAVE SOME INCOMPATIBLE PLACES HERE AND THERE -  
   ! PLACES THAT NEED TO BE MANUALLY MODIFIED, PLEASE USE WITH CAUTION
   ! BUT OVERALL YOU ARE PROBABLY FINE JUST REPLACING SOME "cuda***" with 
   ! "hip***"
   ! 
   ! Hao 
   ! 2024.05
   use iso_c_binding
   implicit none
   save
   ! ======================= HIP enumerators ======================== !
   ! for C these are already setup in amdhip/hipblas/hipsparse headers
   ! however, our fortran subroutines know nothing about those
   ! as such they need to be explicitly setup here
   ! memcpy to host or to device
   integer*8         :: hipMemcpyDeviceToHost 
   integer*8         :: hipMemcpyHostToDevice
   integer*8         :: hipMemcpyDeviceToDevice
   parameter (hipMemcpyHostToDevice=1)
   parameter (hipMemcpyDeviceToHost=2)
   parameter (hipMemcpyDeviceToDevice=3)
   integer*4         :: hipSuccess
   parameter (hipSuccess=1)
   ! matrix operation (whether do transpose)
   ! note this is important as the tests confirm the non-transposed
   ! operation is much faster than transposed ones on GPU
   ! need to avoid the transposed operations (that's why I didn't use QMR)
   integer*4         :: HIPSPARSE_OPERATION_NON_TRANSPOSE
   integer*4         :: HIPSPARSE_OPERATION_TRANSPOSE
   integer*4         :: HIPSPARSE_OPERATION_CONJUGATE_TRANSPOSE
   parameter (HIPSPARSE_OPERATION_NON_TRANSPOSE=0)
   parameter (HIPSPARSE_OPERATION_TRANSPOSE=1)
   parameter (HIPSPARSE_OPERATION_CONJUGATE_TRANSPOSE=2)
   ! matrix type (symmetric may save some spaces)
   integer*4         :: HIPSPARSE_MATRIX_TYPE_GENERAL
   integer*4         :: HIPSPARSE_MATRIX_TYPE_SYMMETRIC
   integer*4         :: HIPSPARSE_MATRIX_TYPE_HERMITIAN
   integer*4         :: HIPSPARSE_MATRIX_TYPE_TRIANGULAR
   parameter (HIPSPARSE_MATRIX_TYPE_GENERAL=0)
   parameter (HIPSPARSE_MATRIX_TYPE_SYMMETRIC=1)
   parameter (HIPSPARSE_MATRIX_TYPE_HERMITIAN=2)
   parameter (HIPSPARSE_MATRIX_TYPE_TRIANGULAR=3)
   ! index base - we all use 1 in fortran
   integer*4         :: HIPSPARSE_INDEX_BASE_ZERO, HIPSPARSE_INDEX_BASE_ONE
   parameter (HIPSPARSE_INDEX_BASE_ZERO=0)
   parameter (HIPSPARSE_INDEX_BASE_ONE=1)
   ! fill mode - upper or lower triangular 
   integer*4         :: HIPSPARSE_FILL_MODE_LOWER
   integer*4         :: HIPSPARSE_FILL_MODE_UPPER
   parameter (HIPSPARSE_FILL_MODE_LOWER=0) 
   parameter (HIPSPARSE_FILL_MODE_UPPER=1)
   ! mat index type
   integer*4         :: HIPSPARSE_INDEX_16U, HIPSPARSE_INDEX_32I
   integer*4         :: HIPSPARSE_INDEX_64I
   parameter (HIPSPARSE_INDEX_16U=1)
   parameter (HIPSPARSE_INDEX_32I=2)
   parameter (HIPSPARSE_INDEX_64I=3)
   ! mat value type
   integer*4         :: HIP_R_32F, HIP_C_32F, HIP_R_64F, HIP_C_64F
   parameter (HIP_R_32F=0)
   parameter (HIP_C_32F=4)
   parameter (HIP_R_64F=1)
   parameter (HIP_C_64F=5)
   ! diag type  unit or not
   integer*4         :: HIPSPARSE_DIAG_TYPE_UNIT 
   integer*4         :: HIPSPARSE_DIAG_TYPE_NON_UNIT 
   parameter (HIPSPARSE_DIAG_TYPE_UNIT=1)
   parameter (HIPSPARSE_DIAG_TYPE_NON_UNIT=0)
   ! Algorithm type for SpMV
   integer*4         :: HIPSPARSE_SPMV_ALG_DEFAULT
   integer*4         :: HIPSPARSE_SPMV_COO_ALG1
   integer*4         :: HIPSPARSE_SPMV_CSR_ALG1
   integer*4         :: HIPSPARSE_SPMV_CSR_ALG2
   integer*4         :: HIPSPARSE_SPMV_COO_ALG2
   parameter (HIPSPARSE_SPMV_ALG_DEFAULT=0)
   parameter (HIPSPARSE_SPMV_COO_ALG1=1)
   parameter (HIPSPARSE_SPMV_CSR_ALG1=2)
   parameter (HIPSPARSE_SPMV_CSR_ALG2=3)
   parameter (HIPSPARSE_SPMV_COO_ALG2=4)
   ! Algorithm type for SpSV
   integer*4         :: HIPSPARSE_SPSV_ALG_DEFAULT
   parameter (HIPSPARSE_SPSV_ALG_DEFAULT=0)
   ! policy for hipsparse solve
   integer*4         :: HIPSPARSE_SOLVE_POLICY_NO_LEVEL
   integer*4         :: HIPSPARSE_SOLVE_POLICY_USE_LEVEL
   parameter (HIPSPARSE_SOLVE_POLICY_NO_LEVEL=0)
   parameter (HIPSPARSE_SOLVE_POLICY_USE_LEVEL=1)
   ! matrix attributes for SpMat
   integer*4         :: HIPSPARSE_SPMAT_FILL_MODE
   integer*4         :: HIPSPARSE_SPMAT_DIAG_TYPE
   parameter (HIPSPARSE_SPMAT_FILL_MODE=0)
   parameter (HIPSPARSE_SPMAT_DIAG_TYPE=1)
   ! hipStream flag    
   integer*4         :: hipStreamNonBlocking
   integer*4         :: hipStreamDefault
   parameter (hipStreamNonBlocking=1)
   parameter (hipStreamDefault=0)
   integer, parameter:: sp = kind(0.)
   ! additional C_ONE in FP32 - useful in mixed precision calculations
   complex(kind=sp)  :: C_ONEs
   parameter (C_ONEs=(1.0, 0.0))
   ! ========================= HIP C interfaces ========================= !
   ! I only included some "might-be-useful" interfaces here...
   ! which ranges from basic hip memory manipulation to hipblas and hipsparse
   ! feel free to add more if you need other interfaces, but keep the 
   ! original naming convention.
   ! Hao
   ! 
   ! ==================================================================== !
   ! ======================= BASIC HIP INTERFACES ======================= !
   ! ==================================================================== !
   interface
   ! hipSetDevice
   integer (c_int) function hipSetDevice( device_idx ) &
    &              bind (C, name="hipSetDevice" ) 
     ! set the current working device according to device_idx
     use iso_c_binding
     implicit none
     integer (c_int), value  :: device_idx
   end function hipSetDevice

   ! hipGetErrorString(err));
   type (c_ptr) function hipGetErrorString( ierr ) &
    &       bind (C, name="hipGetErrorString" )
     ! get the number of available devices
     use iso_c_binding
     implicit none
     integer (c_int), value           :: ierr
   end function hipGetErrorString

   ! hipGetDeviceCount
   integer (c_int) function hipGetDeviceCount( dev_count ) &
    &              bind (C, name="hipGetDeviceCount" )
     ! get the number of available devices
     use iso_c_binding
     implicit none
     type(c_ptr), value  :: dev_count
   end function hipGetDeviceCount

   ! hipMemset
   integer (c_int) function hipMemset( devPtr,value, count ) &
    &              bind (C, name="hipMemset" ) 
     ! fill the count bytes of device memory pointed to by devPr with value
     ! doesn't seems quite useful unless you are filling everything with 
     ! zeros
     use iso_c_binding
     implicit none
     ! devPtr is the *device* pointer to memory 
     type (c_ptr), value  :: devPtr
     integer(c_int), value :: value
     integer(c_size_t), value :: count
   end function hipMemset

   ! hipMalloc
   integer (c_int) function hipMalloc ( devPtr, count ) &
    &              bind (C, name="hipMalloc" ) 
     ! allocate count bytes of memory pointed by devPtr 
     use iso_c_binding
     implicit none
     ! devPtr is the *device* pointer to memory 
     type (c_ptr)  :: devPtr
     integer (c_size_t), value :: count
   end function hipMalloc

   ! hipFree
   integer (c_int) function hipFree(devPtr) & 
    &              bind(C, name="hipFree")
     ! free the chunk of memory used by buffer
     use iso_c_binding
     implicit none
     ! devPtr is the *device* pointer to memory that you want to free
     type (c_ptr),value :: devPtr
   end function hipFree

   ! hipMallocHost
   integer (c_int) function hipMallocHost(hostPtr, count ) &
    &              bind (C, name="hipMallocHost" ) 
     ! allocate count bytes of memory pointed by devPtr 
     use iso_c_binding
     implicit none
     ! devPtr is the *host* pointer to memory 
     type (c_ptr)  :: hostPtr
     integer (c_size_t), value :: count
   end function hipMallocHost

   ! hipFreeHost
   integer (c_int) function hipFreeHost(hostPtr) & 
    &              bind(C, name="hipFreeHost")
     ! free the chunk of memory used by buffer
     use iso_c_binding
     implicit none
     ! devPtr is the *host* pointer to memory that you want to free
     type (c_ptr),value :: hostPtr
   end function hipFreeHost

   ! hipMemcpy
   integer (c_int) function hipMemcpy ( dst, src, count, kind ) & 
    &              bind (C, name="hipMemcpy" )
     use iso_c_binding
     implicit none
     ! copy a chunk of memory from *src* to *dst*
     type (c_ptr), value :: dst, src
     ! count is the size of memory (with unit of size_t)
     integer (c_size_t), value :: count, kind 
     ! kind specifies the direction 
     ! hipMemcpyHostToDevice = 1
     ! hipMemcpyDeviceToHost = 2
     ! hipMemcpyDeviceToDevice = 3
   end function hipMemcpy

   ! hipMemcpyAsync
   integer (c_int) function hipMemcpyAsync( dst, src, count, kind ) & 
    &              bind (C, name="hipMemcpy" )
     use iso_c_binding
     implicit none
     ! copy a chunk of memory from *src* to *dst*
     type (c_ptr), value :: dst, src
     ! count is the size of memory (with unit of size_t)
     integer (c_size_t), value :: count, kind 
     ! kind specifies the direction 
     ! hipMemcpyHostToDevice = 1
     ! hipMemcpyDeviceToHost = 2
     ! hipMemcpyDeviceToDevice = 3
   end function hipMemcpyAsync

   ! hipMemGetInfo
   integer (c_int) function hipMemGetInfo(free, total)  &
    &     bind(C, name="hipMemGetInfo")
     use iso_c_binding
     implicit none
     type(c_ptr),value :: free  ! free memory in bytes
     type(c_ptr),value :: total ! total memory in bytes
   end function hipMemGetInfo

   ! hipDeviceSynchronize
   integer(c_int) function hipDeviceSynchronize() &
    &            bind(C,name="hipDeviceSynchronize")
     use iso_c_binding
     implicit none
   end function hipDeviceSynchronize

   ! hipStreamCreate
   integer(c_int) function hipStreamCreate(stream) & 
    &            bind(C,name="hipStreamCreate")
     use iso_c_binding
     implicit none
     type(c_ptr)::stream ! streamid
   end function hipStreamCreate

   ! hipStreamCreateWithFlags
   integer(c_int) function hipStreamCreateWithFlags(stream,flag ) & 
    &            bind(C,name="hipStreamCreate")
     ! this creates cuda stream (with non-blocking flag)
     use iso_c_binding
     implicit none
     type(c_ptr)::stream ! streamid, out
     integer(c_int),value :: flag ! in
     ! can be 
     ! hipStreamDefault     
     ! hipStreamNonBlocking
   end function hipStreamCreateWithFlags

   ! hipStreamDestroy
   integer(c_int) function hipStreamDestroy(stream) & 
    &            bind(C,name="hipStreamDestroy")
     ! this eliminates a cuda stream
     use iso_c_binding
     implicit none
     type(c_ptr), value ::stream ! streamid
   end function hipStreamDestroy

   ! ==================================================================== !
   ! ===================== CUDA EVENT INTERFACES ======================== !
   ! ==================================================================== !
   ! not really useful for computation, needed for performance assessment

   ! hipEventCreate
   integer(c_int) function hipEventCreate(event) &
    &            bind(C,name="hipEventCreate")
    ! this creates a cuda event 
    use iso_c_binding
    implicit none
    type(c_ptr)           :: event ! event ID
   end function hipEventCreate

   ! hipEventCreateWithFlags
   integer(c_int) function hipEventCreateWithFlags(event, flag) &
    &            bind(C,name="hipEventCreateWithFlags")
    ! this creates a cuda event with given flags
    use iso_c_binding
    implicit none
    type(c_ptr)           :: event ! event ID
    integer(c_int), value :: flag  ! flag 
    ! could be hipEventDefault, hipEventBlockingSync
    ! hipEventDisableTiming
   end function hipEventCreateWithFlags

   ! hipEventRecord
   integer(c_int) function hipEventRecord(event, stream) &
    &            bind(C,name="hipEventRecord")
    ! this records a cuda event (in a given stream)
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: event ! event ID
    type(c_ptr), value    :: stream  ! cuda stream ID
   end function hipEventRecord

   ! hipEventElapsedTime
   integer(c_int) function hipEventElapsedTime(time, starts, ends) &
    &            bind(C,name="hipEventElapsedTime")
    ! this calculates the time betweent two recoreded cuda events
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: time  ! (pointer) in milliseconds
    type(c_ptr), value    :: starts! starting event ID
    type(c_ptr), value    :: ends  ! ending event ID
   end function hipEventElapsedTime

   ! hipEventQuery
   integer(c_int) function hipEventQuery(event) &
    &            bind(C,name="hipEventQuery")
    ! this queries the status of the device preceding a cuda event 
    ! that is recorded by hipEventRecord
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: event ! event ID
   end function hipEventQuery

   ! hipEventDestroy
   integer(c_int) function hipEventDestroy(event) &
    &            bind(C,name="hipEventDestroy")
    ! this distroys the specified cuda event 
    use iso_c_binding
    implicit none
    type(c_ptr), value    :: event ! event ID
   end function hipEventDestroy

   ! ==================================================================== !
   ! ====================== HIPBLAS INTERFACES ========================== !
   ! ==================================================================== !
   ! hipblasCreate - need to be called to initialize hipblas
   integer(c_int) function hipblasCreate(handle) &
    &            bind(C,name="hipblasCreate")
     use iso_c_binding
     implicit none
     type(c_ptr):: handle
   end function hipblasCreate

   ! hipblasDestroy - need to be called after hipblas ends
   integer(c_int) function hipblasDestroy(handle) &
    &            bind(C,name="hipblasDestroy")
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
   end function hipblasDestroy

   ! hipblasSetStream
   integer(c_int) function hipblasSetStream(handle,stream) &
    &            bind(C,name="hipblasSetStream")
     ! this sets the stream to be used by the hipblas lib
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     type(c_ptr),value::stream
   end function hipblasSetStream

   ! hipblasGetStream
   integer(c_int) function hipblasGetStream(handle) &
    &            bind(C,name="hipblasGetStream")
     ! this gets the stream to be used by the hipblas lib
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
   end function hipblasGetStream

   ! hipblasGetVector
   integer(c_int) function hipblasGetVector(n, elemSize, x, incx, y, incy) &
    &            bind(C,name="hipblasGetVector")
     ! this copies n elements from a vector x in GPU memory space to a
     ! vector y in the host memory
     ! each element costs elemSize in the memory
     ! note this is column based if the x/y are matrices
     use iso_c_binding
     implicit none
     integer(c_int), value :: n, elemSize
     type(c_ptr),    value :: x, y
     integer(c_int), value :: incx, incy !strides between elements of x/y
   end function hipblasGetVector

   ! hipblasSetVector
   integer(c_int) function hipblasSetVector(n, elemSize, x, incx, y, incy) &
    &            bind(C,name="hipblasSetVector")
     ! this copies n elements from a vector x in host memory space to a
     ! vector y in the GPU memory
     ! each element costs elemSize in the memory
     ! note this is column based if the x/y are matrices
     use iso_c_binding
     implicit none
     integer(c_int), value :: n, elemSize
     type(c_ptr),    value :: x, y
     integer(c_int), value :: incx, incy !strides between elements of x/y
   end function hipblasSetVector

   ! hipblasSetMatrix
   integer(c_int) function hipblasSetMatrix(nrow, ncol, elemSize, &
    &            A, lda, B, ldb) &
    &            bind(C,name="hipblasSetMatrix")
     ! this copies nrow by ncol  elements from a matrix A in host memory 
     ! to a matrix B in the GPU memory
     ! each element costs elemSize in the memory
     use iso_c_binding
     implicit none
     integer(c_int), value :: nrow, ncol, elemSize 
     type(c_ptr),    value :: A, B ! A-> src B-> dest
     ! the lda and ldb are the size of the leading dimension (rows) of A/B 
     integer(c_int)       :: lda, ldb
   end function hipblasSetMatrix

   ! hipblasGetMatrix
   integer(c_int) function hipblasGetMatrix(nrow, ncol, elemSize, &
    &            A, lda, B, ldb) &
    &            bind(C,name="hipblasGetMatrix")
     ! this copies nrow by ncol  elements from a matrix A in GPU memory 
     ! to a matrix B in the host memory
     ! each element costs elemSize in the memory
     use iso_c_binding
     implicit none
     integer(c_int), value :: nrow, ncol, elemSize 
     type(c_ptr),    value :: A, B ! A-> src B-> dest
     ! the lda and ldb are the size of the leading dimension (rows) of A/B 
     integer(c_int)        :: lda, ldb
   end function hipblasGetMatrix

   ! hipblasDaxpy
   integer(c_int) function hipblasDaxpy(handle,n,alpha,x,incx,y,incy) &
    &            bind(C,name="hipblasDaxpy_v2")
     ! compute y = y + a*x with double precision
     ! note that x and y should be located in GPU memory
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x/y
     real(c_double)        ::alpha        ! scaler a
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y in/out
     integer(c_int),value  :: incy ! strides between elements of x
   end function hipblasDaxpy

   ! hipblasZaxpy
   integer(c_int) function hipblasZaxpy(handle,n,alpha,x,incx,y,incy) &
    &             bind(C,name="hipblasZaxpy_v2")
     ! compute y = y + a*x with complex double precision
     ! note that x and y should be located in GPU memory
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x/y
     complex(c_double)::alpha        ! complex scaler a
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y in/out
     integer(c_int),value  :: incy ! strides between elements of x
   end function hipblasZaxpy

   ! hipblasDcopy
   integer(c_int) function hipblasDcopy(handle,n,x,incx,y,incy) &
    &             bind(C,name="hipblasDcopy_v2")
     ! compute y = x with double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x/y
     type(c_ptr), value    :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y out
     integer(c_int),value  :: incy ! strides between elements of x
   end function hipblasDcopy

   ! hipblasZcopy
   integer(c_int) function hipblasZcopy(handle,n,x,incx,y,incy) &
    &             bind(C,name="hipblasZcopy_v2")
     ! compute y = x with complex double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x/y
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y out
     integer(c_int),value  :: incy ! strides between elements of x
   end function hipblasZcopy

   ! hipblasDdot
   integer(c_int) function hipblasDdot(handle,n,x,incx,y,incy,result) &
    &            bind(C,name="hipblasDdot_v2")
     ! compute result = x dot y with double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x/y
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y out
     integer(c_int),value  :: incy ! strides between elements of x
     type(c_ptr),value     :: result  ! output result (can be in host mem)
   end function hipblasDdot

   ! hipblasZdot
   integer(c_int) function hipblasZdot(handle,n,x,incx,y,incy,result) &
    &            bind(C,name="hipblasZdotc_v2")
     ! compute result = x dot y with complex double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x/y
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: y       ! vector y out
     integer(c_int),value  :: incy ! strides between elements of x
     type(c_ptr),value     :: result  ! output result (can be in host mem)
   end function hipblasZdot

   ! hipblasDnrm2
   integer(c_int) function hipblasDnrm2(handle,n,x,incx,norm) &
    &            bind(C,name="hipblasDnrm2_v2")
     ! compute result = norm(x) in double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: norm    ! output result (can be in host mem)
   end function hipblasDnrm2

   ! hipblasZnrm2 there is no such thing like znrm2!
   integer(c_int) function hipblasZnrm2(handle,n,x,incx,norm) &
    &            bind(C,name="hipblasDznrm2_v2")
     ! compute result = norm(x) in complex double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int),value  :: n   ! n --> size of x
     type(c_ptr),value     :: x       ! vector x
     integer(c_int),value  :: incx ! strides between elements of x
     type(c_ptr),value     :: norm    ! output result (can be in host mem)
   end function hipblasZnrm2


   ! hipblasDscal
   integer(c_int) function hipblasDscal(handle,n,alpha,x,incx) &
    &            bind(C,name="hipblasDscal_v2")
     ! this scales the vector x by the scalar alpha 
     ! x = x/alpha in double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int), value :: n   ! n --> size of x
     real(c_double)        :: alpha   ! double scaler alpha
     type(c_ptr),value     :: x       ! vector x (in/out)
     integer(c_int),value  :: incx ! strides between elements of x
   end function hipblasDscal

   ! hipblasZscal
   integer(c_int) function hipblasZscal(handle,n,alpha,x,incx) &
    &            bind(C,name="hipblasZscal_v2")
     ! this scales the vector x by the scalar alpha 
     ! x = x/alpha in double precision
     use iso_c_binding
     implicit none
     type(c_ptr), value    :: handle ! hipblas context
     integer(c_int), value :: n   ! n --> size of x
     complex(c_double):: alpha   ! complex double scaler alpha
     type(c_ptr),value     :: x       ! vector x (in/out)
     integer(c_int),value  :: incx ! strides between elements of x
   end function hipblasZscal

   ! ==================================================================== !
   ! ====================== CUSPARSE INTERFACES ========================= !
   ! ==================================================================== !

   ! hipsparseCreate - need to be called to initialize hipsparse
   integer(c_int) function hipsparseCreate(hipsparseHandle) & 
    &            bind(C,name="hipsparseCreate")
     use iso_c_binding
     implicit none
     type(c_ptr)::hipsparseHandle  ! hipsparse context
   end function hipsparseCreate

   ! hipsparseDestroy - need to be called after the hipsparse operation ends
   integer(c_int) function hipsparseDestroy(hipsparseHandle) &
    &            bind(C,name="hipsparseDestroy")
     use iso_c_binding
     implicit none
     type(c_ptr),value::hipsparseHandle   ! hipsparse context
   end function hipsparseDestroy

   ! hipsparseGetStream
   integer(c_int) function hipsparseGetStream(hipsparseHandle,stream) &
    &            bind(C,name="hipsparseGetStream")
     use iso_c_binding
     implicit none
     type(c_ptr),value :: hipsparseHandle !handle to cuSPARSE context
     type(c_ptr),value :: stream !out, NULL if the stream is not set
   end function hipsparseGetStream

   ! hipsparseSetStream
   integer(c_int) function hipsparseSetStream(hipsparseHandle,stream) &
    &            bind(C,name="hipsparseSetStream")
     ! this sets the stream to be used by the hipsparse lib
     use iso_c_binding
     implicit none
     type(c_ptr),value :: hipsparseHandle !handle to cuSPARSE context
     type(c_ptr),value :: stream !in, stream id 
   end function hipsparseSetStream

   ! hipsparseCreateMatDescr
   integer(c_int) function hipsparseCreateMatDescr(descrA) &
    &            bind(C,name="hipsparseCreateMatDescr")
     ! this creates an empty sparse mat with descrA
     use iso_c_binding
     implicit none
     type(c_ptr):: descrA ! the matrix descr
   end function hipsparseCreateMatDescr

   ! hipsparseDestroyMatDescr
   integer(c_int) function hipsparseDestroyMatDescr(descrA) &
    &            bind(C,name="hipsparseDestroyMatDescr")
     ! this distroys the mat descrA 
     use iso_c_binding
     implicit none
     type(c_ptr), value :: descrA ! the matrix descr
   end function hipsparseDestroyMatDescr

   ! hipsparseDestroySpMat
   integer(c_int) function hipsparseDestroySpMat(SpMatA) &
    &            bind(C,name="hipsparseDestroySpMat")
     ! this distroys the mat descrA and releases the memory
     ! note this is different from the hipsparseDestroyMatDescr...
     use iso_c_binding
     implicit none
     type(c_ptr), value :: SpMatA ! the matrix descr
   end function hipsparseDestroySpMat

   ! hipsparseGetMatType
   integer(c_int) function hipsparseGetMatType(descrA) &
    &            bind(C,name="hipsparseGetMatType")
     ! this gets the MatType from the mat descrA
     ! the type can be:
     ! HIPSPARSE_MATRIX_TYPE_GENERAL
     ! HIPSPARSE_MATRIX_TYPE_SYMMETRIC
     ! HIPSPARSE_MATRIX_TYPE_HERMITIAN
     ! HIPSPARSE_MATRIX_TYPE_TRIANGULAR
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
   end function hipsparseGetMatType

   ! hipsparseSetMatType
   integer(c_int) function hipsparseSetMatType(descrA, type) &
    &            bind(C,name="hipsparseGetMatType")
     ! this sets the MatType for the mat descrA
     ! the type can be:
     ! HIPSPARSE_MATRIX_TYPE_GENERAL
     ! HIPSPARSE_MATRIX_TYPE_SYMMETRIC
     ! HIPSPARSE_MATRIX_TYPE_HERMITIAN
     ! HIPSPARSE_MATRIX_TYPE_TRIANGULAR
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
     integer(c_int), value:: type !in
   end function hipsparseSetMatType

   ! hipsparseGetMatIndexBase
   integer(c_int) function hipsparseGetMatIndexBase(descrA)             &
    &            bind(C,name="hipsparseGetMatIndexBase")
     ! this gets the index base for the mat descrA
     ! could be zero (c convention) or one (fortran convention)
     ! HIPSPARSE_INDEX_BASE_ZERO
     ! HIPSPARSE_INDEX_BASE_ONE
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
   end function hipsparseGetMatIndexBase

   ! hipsparseSetMatIndexBase
   integer(c_int) function hipsparseSetMatIndexBase(descrA,idxbase)      &
    &            bind(C,name="hipsparseSetMatIndexBase")
     ! this sets the index base for the mat descrA
     ! could be zero (c convention) or one (fortran convention)
     ! HIPSPARSE_INDEX_BASE_ZERO
     ! HIPSPARSE_INDEX_BASE_ONE
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
     integer(c_int), value:: idxbase! the matrix index base
   end function hipsparseSetMatIndexBase


   ! hipsparseGetMatFillMode
   integer(c_int) function hipsparseGetMatFillMode(descrA) &
    &            bind(C,name="hipsparseGetMatFillMode")
     ! this gets the fill mode for the mat descrA
     ! could be 
     ! HIPSPARSE_FILL_MODE_LOWER
     ! HIPSPARSE_FILL_MODE_UPPER
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
   end function hipsparseGetMatFillMode

   ! hipsparseSetMatFillMode
   integer(c_int) function hipsparseSetMatFillMode(descrA,fillmode) &
    &            bind(C,name="hipsparseSetMatFillMode")
     ! this sets the fill mode for the mat descrA
     ! could be 
     ! HIPSPARSE_FILL_MODE_LOWER
     ! HIPSPARSE_FILL_MODE_UPPER
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
     integer(c_int),value :: fillmode !in
   end function hipsparseSetMatFillMode

   ! hipsparseGetMatDiagType
   integer(c_int) function hipsparseGetMatDiagType(descrA) &
    &            bind(C,name="hipsparseGetMatDiagType")
     ! this gets the Mat diag type for the mat descrA
     ! could be 
     ! HIPSPARSE_DIAG_TYPE_NON_UNIT
     ! HIPSPARSE_DIAG_TYPE_UNIT
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
   end function hipsparseGetMatDiagType

   ! hipsparseGetMatDiagType
   integer(c_int) function hipsparseSetMatDiagType(descrA,diagtype) &
    &            bind(C,name="hipsparseSetMatDiagType")
     ! this sets the Mat diag type for the mat descrA
     ! could be 
     ! HIPSPARSE_DIAG_TYPE_NON_UNIT
     ! HIPSPARSE_DIAG_TYPE_UNIT
     use iso_c_binding
     implicit none
     type(c_ptr), value:: descrA ! the matrix descr
     integer(c_int),value :: diagtype ! in
   end function hipsparseSetMatDiagType

   ! hipsparseCreateCsr
   integer(c_int) function hipsparseCreateCsr(SpMatA, &
    &           nRow,nCol,nnz,row,col,val,&
    &           rowIndType, ColIndType, idxbase, valueType) &
    &           bind(C,name="hipsparseCreateCsr")
     ! this creates a CSR sparse matrix with given parameters
     use iso_c_binding
     implicit none
     type(c_ptr)         ::SpMatA ! the SP matrix descr (output)
     integer(c_int),value::nRow  ! number of rows 
     integer(c_int),value::nCol  ! number of cols
     integer(c_int),value::nnz   ! number of none zero elements
     type(c_ptr),value::row      ! row indices, size of rows+1 
     type(c_ptr),value::col      ! col indices, size of nnz 
     type(c_ptr),value::val      ! values, size of nnz
     ! could be zero (c convention) or one (fortran convention)
     ! HIPSPARSE_INDEX_BASE_ZERO
     ! HIPSPARSE_INDEX_BASE_ONE
     integer(c_int),value::rowIndType  ! data type of row
     integer(c_int),value::colIndType  ! data type of col
     integer(c_int),value::idxbase     ! base index of col and row 
     integer(c_int),value::valueType   ! datatype of val
   end function hipsparseCreateCsr

   ! hipsparseCsrGet
   integer(c_int) function hipsparseCsrGet(SpMatA, &
    &           nRow,nCol,nnz,row,col,val,&
    &           rowIndType, ColIndType, idxbase, valueType) &
    &           bind(C,name="hipsparseCsrGet")
     ! this gets all the parameters of a CSR sparse matrix 
     use iso_c_binding
     implicit none
     type(c_ptr), value:: SpMatA ! the matrix descr (input)
     integer(c_int),value::nRow  ! number of rows 
     integer(c_int),value::nCol  ! number of cols
     integer(c_int),value::nnz   ! number of none zero elements
     type(c_ptr),value::row      ! row indices, size of rows+1 
     type(c_ptr),value::col      ! col indices, size of nnz 
     type(c_ptr),value::val      ! values, size of nnz
     ! could be zero (c convention) or one (fortran convention)
     ! HIPSPARSE_INDEX_BASE_ZERO
     ! HIPSPARSE_INDEX_BASE_ONE
     integer(c_int),value::rowIndType  ! data type of row
     integer(c_int),value::colIndType  ! data type of col
     integer(c_int),value::idxbase     ! base index of col and row 
     integer(c_int),value::valueType   ! datatype of val
   end function hipsparseCsrGet
   
   ! hipsparseCreateDnVec
   integer(c_int) function hipsparseCreateDnVec(vec, n, &
    &           val, valueType) bind(C,name="hipsparseCreateDnVec")
     ! this creates a Vec datatype used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr)         ::vec    ! the vec descr (out)
     integer(c_int),value::n      ! size of the vec (in)
     type(c_ptr),value   ::val       ! values, size of n (in) on device
     integer(c_int),value::valueType   ! datatype of val (in)
   end function hipsparseCreateDnVec

   ! hipsparseDestroyDnVec
   integer(c_int) function hipsparseDestroyDnVec(vec) &
    &             bind(C,name="hipsparseDestroyDnVec")
     ! this destroys a Vec datatype used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr), value  ::  vec    ! the vec descr (in)
   end function hipsparseDestroyDnVec

   ! hipsparseDnVecGet
   integer(c_int) function hipsparseDnVecGet(vec, n, &
    &           val, valueType) bind(C,name="hipsparseDnVecGet")
     ! this gets all the field of the Vec datatype used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr)   ,value::vec    ! the vec descr (in)
     integer(c_int),value::n      ! size of the vec (out)
     type(c_ptr)         ::val       ! values, size of n (out)
     integer(c_int),value::valueType   ! datatype of val (out)
   end function hipsparseDnVecGet

   ! hipsparseDnVecGetValues
   integer(c_int) function hipsparseDnVecGetValues(vec,  &
    &           val) bind(C,name="hipsparseDnVecGetValues")
     ! this gets the Vec datatype values used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr), value::vec    ! the vec descr (in)
     type(c_ptr)       ::val    ! values, size of n (out)
   end function hipsparseDnVecGetValues

   ! hipsparseDnVecSetValues
   integer(c_int) function hipsparseDnVecSetValues(vec,  &
    &           val) bind(C,name="hipsparseDnVecSetValues")
     ! this sets the Vec datatype values used in SpMV and SpSV
     use iso_c_binding
     implicit none
     type(c_ptr), value :: vec    ! the vec descr (in)
     type(c_ptr), value :: val    ! values, size of n (in)
   end function hipsparseDnVecSetValues


   ! hipsparseSpMatGetSize
   integer(c_int) function hipsparseSpMatGetSize(SpMatA, &
    &           nRow,nCol,nnz)  bind(C,name="hipsparseSpMatGetSize")
     ! this gets the size of a sparse matrix 
     use iso_c_binding
     implicit none
     type(c_ptr), value:: SpMatA ! the matrix descr (input)
     integer(c_int),value::nRow  ! number of rows 
     integer(c_int),value::nCol  ! number of cols
     integer(c_int),value::nnz   ! number of none zero elements
   end function hipsparseSpMatGetSize

   ! hipsparseSpMatGetValues
   integer(c_int) function hipsparseSpMatGetValues(SpMatA, &
    &           val) bind(C,name="hipsparseSpMatGetValues")
     ! this gets the val of a sparse matrix
     use iso_c_binding
     implicit none
     type(c_ptr), value:: SpMatA ! the matrix descr (input)
     type(c_ptr),value::val      ! values, size of nnz (output)
   end function hipsparseSpMatGetValues

   ! hipsparseSpMatSetValues
   integer(c_int) function hipsparseSpMatSetValues(SpMatA, &
    &           val) bind(C,name="hipsparseSpMatSetValues")
     ! this gets the val of a sparse matrix
     use iso_c_binding
     implicit none
     type(c_ptr), value:: SpMatA ! the matrix descr (input)
     type(c_ptr),value::val      ! values, size of nnz (input)
   end function hipsparseSpMatSetValues

   ! hipsparseSpMatSetAttribute
   integer(c_int) function hipsparseSpMatSetAttribute(SpMatA, attribute, &
    &           data,dataSize) bind(C,name="hipsparseSpMatSetAttribute")
     ! this gets the val of a sparse matrix
     use iso_c_binding
     implicit none
     type(c_ptr), value:: SpMatA ! the matrix descr (input)
     integer(c_int), value:: attribute! parameter type
     ! can be
     ! HIPSPARSE_SPMAT_FILL_MODE or HIPSPARSE_SPMAT_DIAG_TYPE
     integer (c_int)        :: data     ! parameter value
     ! can be 
     ! HIPSPARSE_FILL_MODE_LOWER   HIPSPARSE_FILL_MODE_UPPER
     ! HIPSPARSE_DIAG_TYPE_NON_UNIT   HIPSPARSE_DIAG_TYPE_UNIT
     integer (c_int), value :: dataSize ! parameter value size
   end function hipsparseSpMatSetAttribute

   ! hipsparseSpMV_bufferSize
   integer(c_int) function hipsparseSpMV_bufferSize(handle,opA,    &
    &            alpha, matA, vecX, beta, vecY, computeType,      &
    &            alg, bufferSize) bind(C,name="hipsparseSpMV_bufferSize")
     ! this returns the buffersize needed for matrix-vector multiplition
     ! note that we need to manually allocate the pBuffer 
     ! according to the size
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     real(c_double)    ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     real(c_double)    ::beta ! scaler b
     type(c_ptr), value:: vecY ! vector Y (input)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function hipsparseSpMV_bufferSize

   ! hipsparseSpMV_bufferSize_c
   integer(c_int) function hipsparseSpMV_bufferSize_cmplx(handle,opA, &
    &            alpha, matA, vecX, beta, vecY, computeType,      &
    &            alg, bufferSize) bind(C,name="hipsparseSpMV_bufferSize")
     ! this returns the buffersize needed for matrix-vector multiplition
     ! note that we need to manually allocate the pBuffer 
     ! according to the size
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_double)::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     complex(c_double)::beta ! scaler b
     type(c_ptr), value:: vecY ! vector Y (input)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function hipsparseSpMV_bufferSize_cmplx
   
   ! hipsparseSpMV
   integer(c_int) function hipsparseSpMV(handle,opA, &
    &            alpha, matA, vecX, beta, vecY, computeType,      &
    &            alg, pBuffer) bind(C,name="hipsparseSpMV")
     ! this calculates the matrix-vector multiplition
     ! y = alpha*A*x + beta*y (this is called Axpy in most libs)
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by hipsparseSpMV_bufferSize
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     real(c_double)    ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     real(c_double)    ::beta ! scaler b
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value   :: pBuffer
   end function hipsparseSpMV

   ! hipsparseSpMV_cmplx
   integer(c_int) function hipsparseSpMV_cmplx(handle,opA, &
    &            alpha, matA, vecX, beta, vecY, computeType,      &
    &            alg, pBuffer) bind(C,name="hipsparseSpMV")
     ! this calculates the matrix-vector multiplition
     ! y = alpha*A*x + beta*y (this is called Axpy in most libs)
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by hipsparseSpMV_bufferSize
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_double)::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     complex(c_double)::beta ! scaler b
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value   :: pBuffer
   end function hipsparseSpMV_cmplx

   ! hipsparseSpSV_createDescr
   integer(c_int) function hipsparseSpSV_createDescr(spsvDescr  &
    &            ) bind(C,name="hipsparseSpSV_createDescr")
     ! this creates a handle for the triangle solver hipsparseSpSV
     use iso_c_binding
     implicit none
     type(c_ptr):: spsvDescr ! spsv context handle
   end function hipsparseSpSV_createDescr

   ! hipsparseSpSV_destroyDescr
   integer(c_int) function hipsparseSpSV_destroyDescr(spsvDescr &
    &            ) bind(C,name="hipsparseSpSV_destroyDescr")
     ! this destroys a handle for the triangle solver hipsparseSpSV
     ! and releases the memory associated
     use iso_c_binding
     implicit none
     type(c_ptr), value :: spsvDescr ! spsv context handle
   end function hipsparseSpSV_destroyDescr

   ! hipsparseSpSV_buffersize
   integer(c_int) function hipsparseSpSV_bufferSize(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            bufferSize) bind(C,name="hipsparseSpSV_bufferSize")
     ! this computes the buffersize needed by the triangle solver 
     ! hipsparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     real(c_double)    ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value   :: spsvDescr ! spsv context handle
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function hipsparseSpSV_bufferSize

   ! hipsparseSpSV_buffersize_cmplx
   integer(c_int) function hipsparseSpSV_bufferSize_cmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            bufferSize) bind(C,name="hipsparseSpSV_bufferSize")
     ! this computes the buffersize needed by the triangle solver hipsparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size
     use iso_c_binding
     ! complex double version
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_float) ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value   :: spsvDescr ! spsv context handle
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function hipsparseSpSV_bufferSize_cmplx

   ! hipsparseSpSV_buffersize_dcmplx
   integer(c_int) function hipsparseSpSV_bufferSize_dcmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            bufferSize) bind(C,name="hipsparseSpSV_bufferSize")
     ! this computes the buffersize needed by the triangle solver hipsparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size
     use iso_c_binding
     ! complex double version
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_double) ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value   :: spsvDescr ! spsv context handle
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function hipsparseSpSV_bufferSize_dcmplx

   ! hipsparseSpSV_analysis
   integer(c_int) function hipsparseSpSV_analysis(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            pBuffer ) bind(C,name="hipsparseSpSV_analysis")
     ! this does the analysis phase needed by the triangle solver hipsparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by hipsparseSpMV_bufferSize
     ! double version
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     real(c_double)      ::alpha ! scaler a
     type(c_ptr),value   :: matA ! the matrix descr (input)
     type(c_ptr),value   :: vecX ! vector X (input)
     type(c_ptr),value   :: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
     type(c_ptr),value:: pBuffer ! external buffer
   end function hipsparseSpSV_analysis

   ! hipsparseSpSV_analysis_cmplx
   integer(c_int) function hipsparseSpSV_analysis_cmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            pBuffer ) bind(C,name="hipsparseSpSV_analysis")
     ! this does the analysis phase needed by the triangle solver hipsparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by hipsparseSpMV_bufferSize
     ! complex double version
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_float) ::alpha ! scaler a
     type(c_ptr),value   :: matA ! the matrix descr (input)
     type(c_ptr),value   :: vecX ! vector X (input)
     type(c_ptr),value   :: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
     type(c_ptr),value:: pBuffer ! external buffer
   end function hipsparseSpSV_analysis_cmplx

   ! hipsparseSpSV_analysis_dcmplx
   integer(c_int) function hipsparseSpSV_analysis_dcmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, spsvDescr, &
    &            pBuffer ) bind(C,name="hipsparseSpSV_analysis")
     ! this does the analysis phase needed by the triangle solver hipsparseSpSV
     ! note that we need to manually allocate the pBuffer 
     ! according to the size calculated by hipsparseSpSV_bufferSize
     ! complex float version
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_double)   ::alpha ! scaler a
     type(c_ptr),value   :: matA ! the matrix descr (input)
     type(c_ptr),value   :: vecX ! vector X (input)
     type(c_ptr),value   :: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
     type(c_ptr),value:: pBuffer ! external buffer
   end function hipsparseSpSV_analysis_dcmplx

   ! hipsparseSpSV_solve
   integer(c_int) function hipsparseSpSV_solve(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, &
    &            spsvDescr) bind(C,name="hipsparseSpSV_solve")
     ! this solves the triangle system 
     ! y = Lx or y = Ux
     ! double version
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     real(c_double)      ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
   end function hipsparseSpSV_solve

   ! hipsparseSpSV_solve_cmplx
   integer(c_int) function hipsparseSpSV_solve_cmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, &
    &            spsvDescr) bind(C,name="hipsparseSpSV_solve")
     ! this solves the triangle system 
     ! y = Lx or y = Ux
     ! complex float version
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_float) ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
   end function hipsparseSpSV_solve_cmplx

   ! hipsparseSpSV_solve_dcmplx
   integer(c_int) function hipsparseSpSV_solve_dcmplx(handle,opA, &
    &            alpha, matA, vecX, vecY, computeType, alg, &
    &            spsvDescr) bind(C,name="hipsparseSpSV_solve")
     ! this solves the triangle system 
     ! y = Lx or y = Ux 
     ! complex double version
     use iso_c_binding
     implicit none
     type(c_ptr),value :: handle !handle to cuSPARSE context
     integer(c_int),value::opA ! operation
     complex(c_double) ::alpha ! scaler a
     type(c_ptr), value:: matA ! the matrix descr (input)
     type(c_ptr), value:: vecX ! vector X (input)
     type(c_ptr), value:: vecY ! vector Y (output)
     integer(c_int),value:: computeType
     integer(c_int),value:: alg ! algorithm
     type(c_ptr),value:: spsvDescr ! spsv context handle
   end function hipsparseSpSV_solve_dcmplx
   
   
   ! hipsparseCreateCsrilu0Info
   integer(c_int) function hipsparseCreateCsrilu0Info(info) &
    &           bind(C,name="hipsparseCreateCsrilu02Info")
     ! this create the info pointer for the ilu0 context
     use iso_c_binding
     implicit none
     type(c_ptr)::info
   end function hipsparseCreateCsrilu0Info

   ! hipsparseDestroyCsrilu0Info
   integer(c_int) function hipsparseDestroyCsrilu0Info(info) &
    &           bind(C,name="hipsparseDestroyCsrilu02Info")
     ! this destroys the info pointer for the ilu0 context
     ! and frees its memory
     use iso_c_binding
     implicit none
     type(c_ptr), value ::info
   end function hipsparseDestroyCsrilu0Info

   ! hipsparseDcsrilu0_bufferSize
   integer(c_int) function hipsparseDcsrilu0_bufferSize(handle,m,nnz,descrA,&
    &           valA,rowA,colA,info,bufferSize) &
    &           bind(C,name="hipsparseDcsrilu02_bufferSize")
     ! this calculate the buffersize needed for ilu0 operation 
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     integer(c_int),value::m
     integer(c_int),value::nnz
     type(c_ptr), value:: descrA
     type(c_ptr), value::valA
     type(c_ptr), value::rowA
     type(c_ptr), value::colA
     type(c_ptr),value::info
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function hipsparseDcsrilu0_bufferSize

   ! hipsparseZcsrilu0_bufferSize
   integer(c_int) function hipsparseZcsrilu0_bufferSize(handle,m,nnz,descrA,&
    &           valA,rowA,colA,info,bufferSize) &
    &           bind(C,name="hipsparseZcsrilu02_bufferSize")
     ! this calculate the buffersize needed for ilu0 operation 
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     integer(c_int),value::m
     integer(c_int),value::nnz
     type(c_ptr), value:: descrA
     type(c_ptr), value::valA
     type(c_ptr), value::rowA
     type(c_ptr), value::colA
     type(c_ptr),value::info
     integer(c_size_t)   :: bufferSize ! output the buffer size (in bytes)
   end function hipsparseZcsrilu0_bufferSize

   ! hipsparseDcsrilu0_analysis
   integer(c_int) function hipsparseDcsrilu0_analysis(handle,n,nnz,descrA,&
    &           valA,rowA,colA,info,policy,pBuffer) &
    &           bind(C,name="hipsparseDcsrilu02_analysis")
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     integer(c_int),value::n
     integer(c_int),value::nnz
     type(c_ptr), value:: descrA
     type(c_ptr), value::valA
     type(c_ptr), value::rowA
     type(c_ptr), value::colA
     type(c_ptr), value::info
     integer(c_int),value::policy
     type(c_ptr),value   ::pBuffer
   end function hipsparseDcsrilu0_analysis

   ! hipsparseZcsrilu0_analysis
   integer(c_int) function hipsparseZcsrilu0_analysis(handle,n,nnz,descrA,&
    &           valA,rowA,colA,info,policy,pBuffer) &
    &           bind(C,name="hipsparseZcsrilu02_analysis")
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     integer(c_int),value::n
     integer(c_int),value::nnz
     type(c_ptr), value:: descrA
     type(c_ptr), value::valA
     type(c_ptr), value::rowA
     type(c_ptr), value::colA
     type(c_ptr), value::info
     integer(c_int),value::policy
     type(c_ptr),value   ::pBuffer
   end function hipsparseZcsrilu0_analysis

   ! hipsparseDcsrilu0
   integer(c_int) function hipsparseDcsrilu0(handle,m,nnz, &
    &           descrA,valA,rowA,colA,info, policy, pBuffer) &
    &           bind(C,name="hipsparseDcsrilu02")
     ! this performs the incomplete LU factorization A = LU 
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     integer(c_int),value::m
     integer(c_int),value::nnz
     type(c_ptr), value:: descrA
     type(c_ptr), value::valA  ! output
     type(c_ptr), value::rowA
     type(c_ptr), value::colA
     type(c_ptr), value::info
     integer(c_int),value::policy
     ! in - note this is the address (ptr) to the buffer
     type(c_ptr),value   ::pBuffer 
   end function hipsparseDcsrilu0

   ! hipsparseZcsrilu0
   integer(c_int) function hipsparseZcsrilu0(handle,m,nnz, &
    &           descrA,valA,rowA,colA,info, policy, pBuffer) &
    &           bind(C,name="hipsparseZcsrilu02")
     ! this performs the incomplete LU factorization A = LU 
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     integer(c_int),value::m
     integer(c_int),value::nnz
     type(c_ptr), value:: descrA
     type(c_ptr), value::valA  ! output
     type(c_ptr), value::rowA
     type(c_ptr), value::colA
     type(c_ptr), value::info
     integer(c_int),value::policy
     ! in - note this is the address (ptr) to the buffer
     type(c_ptr),value   ::pBuffer 
   end function hipsparseZcsrilu0

   ! hipsparseXcsrilu0_zeroPivot
   integer(c_int) function hipsparseXcsrilu0_zeroPivot(handle,info,loc) &
    &           bind(C,name="hipsparseXcsrilu02_zeroPivot")
     use iso_c_binding
     implicit none
     type(c_ptr),value::handle
     type(c_ptr),value::info
     type(c_ptr),value::loc !output, in hostmem
   end function hipsparseXcsrilu0_zeroPivot

   ! =========================================================================!
   ! ====================== CUSTOM KERNEL INTERFACES =========================!
   ! =========================================================================!

   ! kernelc_s2d 
   subroutine kernelc_s2d( single, double, count ) & 
    &              bind (C, name="kernelc_s2d" )
     use iso_c_binding
     implicit none
     ! force convert a chunk of memory from *single* to *double*
     type (c_ptr), value :: single, double ! pointers 
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value :: count 
   end subroutine kernelc_s2d

   ! kernelc_d2s 
   subroutine kernelc_d2s( double, single, count ) & 
    &              bind (C, name="kernelc_d2s" )
     use iso_c_binding
     implicit none
     ! force convert a chunk of memory from *double* to *single*
     type (c_ptr),value   :: double, single ! pointers
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value :: count 
   end subroutine kernelc_d2s

   ! kernelc_hadar
   subroutine kernelc_hadar( a, b, c, count ) & 
    &              bind (C, name="kernelc_hadar" )
     use iso_c_binding
     implicit none
     ! perform hadamard multiply c = a*b, (real)
     type (c_ptr),value   :: a, b, c ! pointers
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value :: count 
   end subroutine kernelc_hadar

   ! kernelc_hadac
   subroutine kernelc_hadac( a, b, c, count ) & 
    &              bind (C, name="kernelc_hadac" )
     use iso_c_binding
     implicit none
     ! perform hadamard multiply c = a*b, (complex)
     ! force convert two doubles into one complex double
     type (c_ptr),value   :: a, b, c ! pointers
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value :: count 
   end subroutine kernelc_hadac

   ! kernelc_xpbyc
   subroutine kernelc_xpbyc( a, b, c, count ) & 
    &              bind (C, name="kernelc_xpbyc" )
     use iso_c_binding
     implicit none
     ! perform y = x + by, (complex)
     type (c_ptr),value          :: a, c ! pointers
     complex(c_double), value    :: b ! scaler b
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value      :: count 
   end subroutine kernelc_xpbyc

   ! kernelc_update_p
   subroutine kernelc_update_pc( r, v, beta, omega, p, count ) & 
    &              bind (C, name="kernelc_update_pc" )
     use iso_c_binding
     implicit none
     ! perform y = x + by, (complex)
     ! force convert a chunk of memory from *single* to *double*
     type (c_ptr),value          :: p, r, v ! pointers
     complex(c_double), value    :: beta ! scaler beta
     complex(c_double), value    :: omega ! scaler omega
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value      :: count 
   end subroutine kernelc_update_pc

   ! kernelc_update_x
   subroutine kernelc_update_xc( ph, sh, alpha, omega, x, count ) & 
    &              bind (C, name="kernelc_update_xc" )
     use iso_c_binding
     implicit none
     ! perform y = x + by, (complex)
     ! force convert a chunk of memory from *single* to *double*
     type (c_ptr),value          :: ph, sh, x ! pointers
     complex(c_double), value    :: alpha ! scaler alpha
     complex(c_double), value    :: omega ! scaler omega
     ! count is the size of memory (with unit of size_t)
     integer (c_int), value      :: count 
   end subroutine kernelc_update_xc

   ! cf_free
   subroutine cf_free(ptr)  bind (C, name="free" )
     use iso_c_binding
     implicit none
     ! this calls free from the c side to deal with arrays
     ! associated by c_loc
     type (c_ptr),value          :: ptr 
   end subroutine cf_free

  ! cf_hookCtx
   integer(c_int) function cf_hookDev(device_idx) & 
    &              bind (C, name="cf_hookDev" )
     ! bind the current context to a given device
     ! need some extra tweaks to use the "internal" functions instead of 
     ! the NVML lib - this is actually fine for CUDA, but not possible 
     ! for other APIs (AMD ROCM or Intel oneAPI) as they don't have the 
     ! equivalents of NVML
     use iso_c_binding
     implicit none
     ! the current device idx
     integer(c_int),value ::  device_idx
   end function cf_hookDev

  ! cf_resetFlag
   integer(c_int) function cf_resetFlag(device_idx) & 
    &              bind (C, name="cf_resetFlag" )
     ! reset the computation flag of the GPU device
     ! can be "hipDeviceScheduleSpin" or 
     ! "hipDeviceScheduleYield"
     use iso_c_binding
     implicit none
     ! the current device idx
     integer(c_int),value ::  device_idx
   end function cf_resetFlag

   end interface  

end module hipFortMap
