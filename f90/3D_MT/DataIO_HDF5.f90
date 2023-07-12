! *****************************************************************************
!   This is the HDF5 I/O version coded by Spencer Wilbur (USGS) under the direction
!   of Anna Kelbert, Jul. 2023. File intentionally misnamed; it is renamed at the
!   level of the configuration file.
module DataIO
  ! This module contains io routines for reading and writing the data vectors
  ! Version: 3D MT
  use hdf5
  use math_constants
  use file_units
  use utilities
  use dataspace
  use gridcalc
  use transmitters
  use receivers
  use datatypes

  implicit none

  private

!   switch between data formats by leaving uncommented one of the options below
  interface write_dataVectorMTX
	MODULE PROCEDURE write_hdf5_data  
  end interface
  
  interface read_dataVectorMTX 
    !MODULE PROCEDURE read_Z_list
    MODULE PROCEDURE read_hdf5_data
  end interface

  interface deall_dataFileInfo
    MODULE PROCEDURE deall_fileInfo
  end interface

  public     ::  write_dataVectorMTX, read_dataVectorMTX, deall_dataFileInfo

  type :: data_file_block

      ! this block of information constitutes user preferences about the data format;
      ! there is one entry per each transmitter type and data type... (iTxt,iDt)
      ! if there are multiple data blocks of the same transmitter & data types,
      ! the last value is used.
      character(200) :: info_in_file
      character(20)  :: sign_info_in_file
      integer        :: sign_in_file
      character(20)  :: units_in_file
      real           :: origin_in_file(2)
      real           :: geographic_orientation

     ! these lists contain the indices into the data vector for each data type;
     ! they make it possible to sort the data by receiver for output.
     ! no data denoted by zero index; dimensions (nTx) and (nTx,nRx).
     ! these indices are typically allocated as we read the data file
     integer, pointer, dimension(:)   :: tx_index
     integer, pointer, dimension(:)   :: dt_index
     integer, pointer, dimension(:,:) :: rx_index

     ! some transmitter types and data types don't go together
     logical         :: defined

  end type data_file_block

  ! private dictionary of data block info dimension (nTxt,nDt)
  ! where nTxt = number of all possible transmitter types
  !       nDt  = number of all possible data types
  ! number of transmitter types comes from the DICT/txTypes module
  ! and defines the number of conceptually different types of sources
  type (data_file_block), pointer, save, private, dimension(:,:) :: fileInfo

  ! we are converting from an "old format" to a "new format"
  ! the only difference being that in the new format, there is
  ! an additional line in the head that indicates transmitter type.
  ! on output, use the same format as on input. AK 25 May 2018
  logical, save, private  :: old_data_file_format = .true.
  integer(HID_T), private, save   :: file_id, group_id, attr_id, dset_id, dspace_id, atype_id, aspace_id, dtype_id ! file, data set, and dataspace handles
Contains

!**********************************************************************
! Sorts out the data block header

  function DataBlockHeader(txType,dataType) result (header)

    integer, intent(in)         :: txType
    integer, intent(in)         :: dataType
    character(200)              :: header

    select case (dataType)

       case(Full_Impedance,Off_Diagonal_Impedance,Full_Vertical_Components)
          header = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error'

       case(Full_Interstation_TF)
          header = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Ref_Code Ref_Lat '// &
                   'Ref_Lon Ref_X(m) Ref_Y(m) Ref_Z(m) Component Real Imag Error'

       case(Off_Diagonal_Rho_Phase,Phase_Tensor)
          header = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Value Error'

    end select

  end function DataBlockHeader

  ! **************************************************************************
  ! Cleans up and deletes type dictionary at end of program execution
  subroutine init_fileInfo(nTxt,nDt,nTx,nRx)

    integer, intent(in) :: nTxt,nDt
    integer, intent(in), optional :: nTx,nRx
    integer     :: istat,iTxt,iDt

    allocate(fileInfo(nTxt,nDt),STAT=istat)


     do iTxt = 1,nTxt
       do iDt = 1,nDt
         fileInfo(iTxt,iDt)%defined = .false.
         if (present(nTx) .and. present(nRx)) then
           allocate(fileInfo(iTxt,iDt)%tx_index(nTx),STAT=istat)
           allocate(fileInfo(iTxt,iDt)%dt_index(nTx),STAT=istat)
           allocate(fileInfo(iTxt,iDt)%rx_index(nTx,nRx),STAT=istat)
         end if
       end do
     end do

  end subroutine init_fileInfo

  ! **************************************************************************
  ! Cleans up and deletes type dictionary at end of program execution
  subroutine deall_fileInfo()

    integer     :: i,j, istat

    if (associated(fileInfo)) then

     do i = 1,size(fileInfo,1)
       do j = 1,size(fileInfo,2)
          if (associated(fileInfo(i,j)%tx_index)) then
             deallocate(fileInfo(i,j)%tx_index,STAT=istat)
          end if
          if (associated(fileInfo(i,j)%dt_index)) then
             deallocate(fileInfo(i,j)%dt_index,STAT=istat)
          end if
          if (associated(fileInfo(i,j)%rx_index)) then
             deallocate(fileInfo(i,j)%rx_index,STAT=istat)
          end if
       end do
     end do

     deallocate(fileInfo,STAT=istat)

    end if

  end subroutine deall_fileInfo

!**********************************************************************  
 
!********************************************************************** 
   !Subroutine to open the hdf5 for reading 
subroutine open_read_hdf5(cfile)
    character(*), intent(in)                :: cfile
    integer                                 :: hdferr
    logical                                 :: lexist

    inquire(file = cfile, exist = lexist)
    if (lexist) then
        CALL h5open_f(hdferr)
        CALL h5fopen_f(cfile, H5F_ACC_RDONLY_F, file_id, hdferr)
    else 
        write(0,*) 'No HDF5 file to read'
    end if 

end subroutine open_read_hdf5

!********************************************************************** 
   !This subroutine will either create a new hdf5 file based on the name 
   !given in the input or open an already exisiting file to be appended to 
subroutine open_hdf5(cfile)
    character(*), intent(in)                :: cfile
    integer                                 :: hdferr
    CHARACTER(LEN=4), PARAMETER  :: data_group = "Data"
    CHARACTER(LEN=7), PARAMETER  :: data_mt_group = "Data/MT"
    logical                      :: lexist

    inquire(file = cfile, exist = lexist)
    if (lexist) then 
        CALL h5open_f(hdferr)
        CALL h5fopen_f(cfile, H5F_ACC_RDWR_F, file_id, hdferr)
    else 
            CALL h5open_f(hdferr) ! throws error if cannot open 
            CALL h5fcreate_f(cfile, H5F_ACC_TRUNC_F, file_id, hdferr)!create the file using the variable cfile given to the terminal when running the script 
        
            !Create the groups for the dataset
            CALL h5gcreate_f(file_id, data_group, group_id, hdferr) 
            CALL h5gcreate_f(file_id, data_mt_group, group_id, hdferr)
        
    end if
    
end subroutine open_hdf5
!********************************************************************** 
subroutine close_hdf5(cfile)
    character(*), intent(in)                :: cfile
    integer                                 :: hdferr

    ! CALL h5gclose_f(group_id, hdferr)
    CALL h5fclose_f(file_id, hdferr)
    CALL h5close_f(hdferr) 
   

end subroutine close_hdf5
!********************************************************************** 
subroutine write_hdf5_attr(att_name, att_val, d_id)

    character(*), intent(in)                :: att_name
    character(*), intent(in)     :: att_val
    INTEGER(HID_T), intent(in), optional    :: d_id
    integer                                 :: hdferr
   

   
    INTEGER(HSIZE_T)  :: attrlen  ! Length of the attribute string
    INTEGER(HSIZE_T), DIMENSION(1) :: dimsc = 1 ! Scalar or single value string


    attrlen = len(att_val) ! this makes sure that unwanted characters are not stored in each string value 
    
    if (present(d_id)) then
    !Create string from att_name 
        CALL h5screate_simple_f(1, dimsc, aspace_id, hdferr)
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
        CALL h5tset_size_f(atype_id, attrlen, hdferr)

        CALL h5acreate_f(d_id, att_name , atype_id, aspace_id, attr_id, hdferr)
        CALL h5awrite_f(attr_id, atype_id, att_val, dimsc, hdferr)
        CALL h5aclose_f(attr_id, hdferr)
        CALL h5sclose_f(aspace_id, hdferr)

    else
        CALL h5screate_simple_f(1, dimsc, aspace_id, hdferr)
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
        CALL h5tset_size_f(atype_id, attrlen, hdferr)
        CALL h5acreate_f(group_id, att_name , atype_id, aspace_id, attr_id, hdferr)
        CALL h5awrite_f(attr_id, atype_id, att_val, dimsc, hdferr)
        CALL h5aclose_f(attr_id, hdferr)
        CALL h5sclose_f(aspace_id, hdferr)
      
    end if  

end subroutine write_hdf5_attr
!********************************************************************** 
function read_hdf5_attr(att_name, d_id) result( att_data)

    character(*), intent(in)                :: att_name
    INTEGER(HID_T), intent(in),optional       :: d_id
    integer                                 :: hdferr
    INTEGER(8), DIMENSION(1)              :: maxdims !Read buffer dimension
    INTEGER(HSIZE_T)  :: attrlen  ! Length of the attribute string
    INTEGER(HSIZE_T), DIMENSION(1) :: dimsc = 1 ! Scalar or single value string
    CHARACTER(len=100),allocatable, Dimension(1)                 :: att_data(:) !attribute read buffer

    if (present(d_id)) then
        CALL h5aopen_f(d_id, att_name, attr_id, hdferr)
        CALL H5Aget_space_f(attr_id, dspace_id, hdferr)
        CALL H5Sget_simple_extent_dims_f(dspace_id, dimsc, maxdims, hdferr)

        ALLOCATE(att_data(dimsc(1)))

        ! Create the memory datatype.
        CALL H5Tcopy_f(H5T_FORTRAN_S1, atype_id, hdferr)
        CALL H5Tset_size_f(atype_id, attrlen, hdferr)

        ! get the attribute data type 
        CALL H5Aget_type_f(attr_id, atype_id, hdferr)
        CALL h5tget_size_f(atype_id, attrlen, hdferr)

        !Read the attribute
        CALL h5aread_f(attr_id, atype_id, att_data, dimsc, hdferr )
        call h5aclose_f( attr_id, hdferr)
    else

        CALL h5aopen_f(group_id, att_name, attr_id, hdferr)

        ! Get dataspace and allocate memory for read buffer.
        CALL H5Aget_space_f(attr_id, dspace_id, hdferr)
        CALL H5Sget_simple_extent_dims_f(dspace_id, dimsc, maxdims, hdferr)

        ALLOCATE(att_data(dimsc(1)))
       
        ! ! Create the memory datatype.
        CALL H5Tcopy_f(H5T_FORTRAN_S1, atype_id, hdferr)
        CALL H5Tset_size_f(atype_id, attrlen, hdferr)
    
        CALL H5Aget_type_f(attr_id, atype_id, hdferr)
        CALL h5tget_size_f(atype_id, attrlen, hdferr)
    
        !Read the attribute 
        CALL h5aread_f(attr_id, atype_id, att_data, dimsc, hdferr )
        call h5aclose_f( attr_id, hdferr)
    
    end if
 
    
end function read_hdf5_attr


!********************************************************************** 
subroutine write_hdf5_txdict(cfile,order)

    character(*), intent(in)                  :: cfile
    character(*), intent(in), optional        :: order

    ! local
    INTEGER(HSIZE_T)  :: ndat     ! Number of data points in an array
    INTEGER(HSIZE_T)  :: attrlen = 100 ! Length of the attribute string
    INTEGER(HSIZE_T), DIMENSION(1) :: dimsc = 1 ! Scalar or single value string
    INTEGER(HSIZE_T), DIMENSION(1) :: dim1d ! Datasets dimensions for 1D arrays
    INTEGER(HSIZE_T), DIMENSION(2) :: dim2d ! Datasets dimensions for 2D arrays
    integer           :: hdferr, ii, istat

    CHARACTER(LEN=14), parameter :: data_mt_txdict_group = "Data/MT/txdict"
    real(kind=8), allocatable    :: period(:)
    character(len=100)           :: order_attr


    ! first create linear data arrays for writing
    allocate(period(size(txDict)), stat=istat)
    do ii = 1, size(txDict)
        period(ii) = txDict(ii)%period
    end do
    dim1d(1) = size(txDict)
    write(0,*) 'Periods: ',period

    if (present(order)) then
        order_attr = order
    else
        order_attr = 'ascending'
    end if

    ! open cfile and create the txdict group
    CALL open_hdf5(cfile)
    !   !Create the groups for the dataset
    CALL h5gcreate_f(file_id, data_mt_txdict_group, group_id, hdferr)


    ! create attributes
    CALL write_hdf5_attr('order', order_attr)

    ! ! write the linear data arrays
    CALL h5screate_simple_f(1, dim1d, dspace_id, hdferr)
    CALL h5dcreate_f(group_id, 'periods', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, period, dim1d, hdferr)


    ! ! close the data set and data space
    CALL h5dclose_f(dset_id, hdferr)
    CALL h5sclose_f(dspace_id, hdferr)

    ! ! close the HDF5 file
    CALL close_hdf5(cfile) 

    ! ! deallocate temporary linear arrays
    deallocate(period, stat=istat)

end subroutine write_hdf5_txdict
!********************************************************************** 
subroutine write_hdf5_rxdict(cfile,primary_coords)

    character(*), intent(in)                  :: cfile
    character(*), intent(in), optional        :: primary_coords

    ! local
    INTEGER(HSIZE_T)               :: ndat     ! Number of data points in an array
    INTEGER(HSIZE_T)               :: attrlen = 100 ! Length of the attribute string
    INTEGER(HSIZE_T), DIMENSION(1) :: dimsc = 1 ! Scalar or single value string
    INTEGER(HSIZE_T), DIMENSION(1) :: dim1d ! Datasets dimensions for 1D arrays
    INTEGER(HSIZE_T), DIMENSION(2) :: dim2d ! Datasets dimensions for 2D arrays
    integer                        :: hdferr, ii, istat, rz, k
    real(8)                        :: lat, lon
    CHARACTER(5)                   :: codes

    CHARACTER(LEN=18), parameter  :: data_mt_rxdict_group = "Data/MT/rxdict"
    real(8), allocatable          :: rxdict_elv(:),rxdict_lat(:), rxdict_lon(:)

    real(8), allocatable          :: rxdict_x(:),rxdict_y(:), rxdict_z(:)

    real(8), allocatable          :: rxdict_xyz(:,:)
    CHARACTER(len=5), allocatable :: rxdict_codes(:) 
    character(len=100)            :: coord_attr


    
    ! first create linear data arrays for writing
    rz = size(rxDict) 

    allocate(rxdict_codes(rz), rxdict_elv(rz), rxdict_lat(rz), rxdict_lon(rz), STAT = istat)
    allocate(rxdict_x(rz), rxdict_y(rz), rxdict_z(rz), STAT = istat)
    allocate(rxdict_xyz(3,rz), STAT = istat)

    !open the file and create the rxdict group
    CALL open_hdf5(cfile)
    CALL h5gcreate_f(file_id, data_mt_rxdict_group, group_id, hdferr)
    
   
    !Assign dimensions for the xyz dataset beacuse it is 2D and everything else in rxdict is 1D
    dim2d(1) = 3
    dim2d(2) = rz

    do ii= 1,rz

        rxdict_elv = rxDict(ii)%x(3)
        rxdict_codes(ii) = rxDict(ii)%id 
        rxdict_lat(ii) = rxDict(ii)%x(1)
        rxdict_lon(ii) = rxDict(ii)%x(2)
     
        rxdict_xyz(:,ii) = rxDict(ii)%x
      
    end do

    dim1d(1) = size(rxDict)
   
    if (present(primary_coords)) then
        coord_attr = primary_coords
    else
        coord_attr = 'latlon'
    end if


    ! create attributes
    CALL write_hdf5_attr('primarycoords', coord_attr)


    CALL h5screate_simple_f(1, dim1d, dspace_id, hdferr)

    ! write the latitude data array
    CALL h5dcreate_f(group_id, 'lat', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, rxdict_lat, dim1d, hdferr)


    !write longitude array 
    CALL h5dcreate_f(group_id, 'lon', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, rxdict_lon, dim1d, hdferr)


    !write elevation array 
    CALL h5dcreate_f(group_id, 'elv', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, rxdict_elv, dim1d, hdferr)
 
    !write codes array 
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr) !use atype_id here for character strings 
    CALL H5tset_size_f(atype_id, int(5, size_t),hdferr)
    CALL h5dcreate_f(group_id, 'codes', atype_id, dspace_id, dset_id, hdferr)
    CALL h5dwrite_f(dset_id, atype_id,rxdict_codes, dim1d, hdferr)

    !create new data space for xyz data set that is 2D and write data to dataset 
    CALL h5screate_simple_f(2, dim2d, dspace_id, hdferr)
    CALL h5dcreate_f(group_id, 'xyz', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, rxdict_xyz, dim2d, hdferr)

    !Write the Attributes
    CALL write_hdf5_attr('angleunits', 'null', dset_id)
    CALL write_hdf5_attr('falseeasting', 'null', dset_id)
    CALL write_hdf5_attr('falsenorthing', 'null', dset_id)
    CALL write_hdf5_attr('fixedorient', 'null', dset_id)
    CALL write_hdf5_attr('geoid', 'null', dset_id)
    CALL write_hdf5_attr('maplatlimit', 'null', dset_id)
    CALL write_hdf5_attr('maplonlimit', 'null', dset_id)
    CALL write_hdf5_attr('mapprojection', 'null', dset_id)
    CALL write_hdf5_attr('nparalels', 'null', dset_id)
    CALL write_hdf5_attr('origin', '45.276, -119.634', dset_id)
    CALL write_hdf5_attr('scalefactor', 'null', dset_id)
    CALL write_hdf5_attr('userinfo', 'null', dset_id)
    CALL write_hdf5_attr('zone', 'null', dset_id)
    ! close the HDF5 file
    CALL close_hdf5(cfile) 
  


    ! deallocate temporary linear arrays
    deallocate(rxdict_codes,rxdict_lat,rxdict_lon,rxdict_elv, stat=istat)
    deallocate(rxdict_xyz, STAT = istat)

end subroutine write_hdf5_rxdict
!********************************************************************** 
subroutine write_hdf5_typelist(allData, cfile)
    character(*), intent(in)                :: cfile
    type(dataVectorMTX_t), intent(in) :: allData

    ! local
    INTEGER(HSIZE_T)               :: ndat     ! Number of data points in an array
    INTEGER(HSIZE_T)               :: attrlen = 100 ! Length of the attribute string
    INTEGER(HSIZE_T), DIMENSION(1) :: dimsc = 1 ! Scalar or single value string
    INTEGER(HSIZE_T), DIMENSION(1) :: dim1d ! Datasets dimensions for 1D arrays
    INTEGER(HSIZE_T), DIMENSION(2) :: dim2d ! Datasets dimensions for 2D arrays
    INTEGER(HSIZE_T)                :: size_type_t, size_type_z
    integer                        :: hdferr, ii, istat, iDt, ncomp, k
    
    CHARACTER(len=2), allocatable  :: id_t(:) !data array for type components 
    CHARACTER(len=3), allocatable  :: id_z(:) !data array for type components 

    CHARACTER(LEN=18), parameter :: data_typelist = "Data/MT/typelist"
    CHARACTER(LEN=19), parameter :: data_typelist_T = "Data/MT/typelist/T"
    CHARACTER(LEN=19), parameter :: data_typelist_Z = "Data/MT/typelist/Z"
  

    CALL open_hdf5(cfile)
    CALL h5gcreate_f(file_id, data_typelist, group_id, hdferr)
   
        WRITE_DATA_TYPE: do k = 1, alldata%d(1)%ndt !all data, d= one datablock to look foir datatypes, ndt =datatypes 
        iDt = alldata%d(1)%data(k)%datatype
        ! do iDt = 1,size(typeDict) 
            select case (iDt)
                case(Full_Vertical_Components)
                    
                    ncomp = typeDict(iDt)%nComp
                    dim1d = (ncomp/2) ! define the dimensions for the data set using the number of components divided by 2
                        CALL h5gcreate_f(file_id, data_typelist_T, group_id, hdferr)
                            ! Write the Attributes
                        CALL write_hdf5_attr('longname', typeDict(iDt)%name)
                        CALL write_hdf5_attr('units', typeDict(iDt)%units)
                        CALL write_hdf5_attr('complex', '1')
                        CALL write_hdf5_attr('description', 'Vertical Field Transfer Functions(tipper)')
                        CALL write_hdf5_attr('externalurl', 'http://www.iris.edu/dms/products/emtf/tipper.html')
                        CALL write_hdf5_attr('input', 'H')
                        CALL write_hdf5_attr('intention', 'primary data type')
                        CALL write_hdf5_attr('output', 'H')
                        CALL write_hdf5_attr('tag', 'tipper')
                        size_type_t = 2 
                        CALL h5screate_simple_f(1, dim1d, dspace_id, hdferr)
                        CALL h5tcopy_f(H5T_FORTRAN_S1, atype_id, hdferr)
                        CALL H5tset_size_f(atype_id, size_type_t ,hdferr)
                        CALL h5dcreate_f(group_id, "components", atype_id, dspace_id, dset_id, hdferr)


                        allocate(id_t(ncomp/2),STAT=istat)
                        do ii = 1, ncomp/2
                            id_t(ii) = typeDict(iDt)%id(ii) ! Had to allocate the components to a new variable in order to not overwrite the first component using h5dwrite
                           
                        end do
                        CALL h5dwrite_f(dset_id, atype_id, id_t, dim1d, hdferr)
                     
                        deallocate(id_t,STAT=istat)
                case(Full_Impedance)
                    ncomp = typeDict(iDt)%nComp
                    dim1d = (ncomp/2)
                        ! FOR TYPELIST Z
                        CALL h5gcreate_f(file_id, data_typelist_Z, group_id, hdferr) !This creates the Group T under typelists 
                        !Write the Attributes 
                        CALL write_hdf5_attr('longname', typeDict(iDt)%name)
                        CALL write_hdf5_attr('units', typeDict(iDt)%units)
                        CALL write_hdf5_attr('complex', '1')
                        CALL write_hdf5_attr('description', 'MT Impdeance')
                        CALL write_hdf5_attr('externalurl', 'http://www.iris.edu/dms/products/emtf/impedance.ht')
                        CALL write_hdf5_attr('input', 'H')
                        CALL write_hdf5_attr('intention', 'primary data type')
                        CALL write_hdf5_attr('output', 'E')
                        CALL write_hdf5_attr('tag', 'impedance') 
                        
                        size_type_z = 3 
                        CALL h5screate_simple_f(1, dim1d, dspace_id, hdferr)
                        CALL h5tcopy_f(H5T_FORTRAN_S1, atype_id, hdferr)
                        CALL H5tset_size_f(atype_id, size_type_z,hdferr)
                        CALL h5dcreate_f(group_id, "components", atype_id, dspace_id, dset_id, hdferr)


                        allocate(id_z(ncomp/2), STAT=istat)
                        do ii = 1, ncomp/2
                            id_z(ii) = typeDict(iDt)%id(ii) ! Had to allocate the components to a new variable in order to not overwrite the first component using h5dwrite
                            
                        end do 
                        CALL h5dwrite_f(dset_id, atype_id, id_z, dim1d, hdferr)
                       
                        deallocate(id_z,STAT=istat)
            end select 
        end do WRITE_DATA_TYPE
        CALL close_hdf5(cfile)
end subroutine write_hdf5_typelist

!********************************************************************** 
subroutine write_hdf5_datablocks(allData, cfile)
    character(*), intent(in)          :: cfile
    type(dataVectorMTX_t), intent(in) :: allData
      ! local
    INTEGER(HSIZE_T)               :: ndat     ! Number of data points in an array
    INTEGER(HSIZE_T)               :: attrlen = 100 ! Length of the attribute string
    INTEGER(HSIZE_T), DIMENSION(1) :: dimsc = 1 ! Scalar or single value string
    INTEGER(HSIZE_T), DIMENSION(1) :: dim1d ! Datasets dimensions for 1D arrays
    INTEGER(HSIZE_T), DIMENSION(2) :: dim2d ! Datasets dimensions for 2D arrays
    integer                        :: hdferr, ii, istat, iDt, ncomp, iTx, nRx, myint


    CHARACTER(LEN=17), PARAMETER :: data_block_group = 'Data/MT/datablock'
    !DATA BLOCKS NAMES 
    character(len=27)                :: dblk_num, dbTZ, dblk
    character(LEN=3), parameter      :: T = "/T"
    character(len=3), PARAMETER      :: Z = "/Z"



    CALL open_hdf5(cfile)
   
    do iTx= 1, size(alldata%d) !this provides number of data blocks !
        myint = iTx
        write(dblk_num, '(a, a1, I0.2)') data_block_group, '.', myint
        dblk = dblk_num(18:19)
        CALL h5gcreate_f(file_id, dblk_num, group_id, hdferr)

         ! Create attributes for Datablock Headers
        CALL h5screate_simple_f(1, dimsc, aspace_id, hdferr)
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
        CALL h5tset_size_f(atype_id, int(27, size_t), hdferr)

        !Create dataset attribute.
        CALL h5acreate_f(group_id, "Tx", atype_id, aspace_id,attr_id, hdferr)
        CALL h5awrite_f(attr_id, atype_id, dblk, dimsc, hdferr)
    
 


        WRITE_DATA_TYPE: do ii = 1, alldata%d(iTx)%ndt !all data, d=datablocks ndt =datatypes 
        iDt = alldata%d(iTx)%data(ii)%datatype
            select case(iDt) !
                case(Full_Vertical_Components)
                    ncomp = allData%d(iTx)%data(ii)%nComp
                    nRx = allData%d(iTx)%data(ii)%nSite
                    dim2d = (/ncomp,nRx/)
            
                    dim1d = (/nRx/)  !set up dimension size ! THIS ALLOWED FOR SETTING DIMENSION SIZE USING A VARIABLE
                    dbTZ = trim(dblk_num)//trim(T)
                  
                    !this is for T standard dev error 
                    CALL h5gcreate_f(file_id, dbTZ, group_id, hdferr) 
                    ! Write Attributes for Tipper 
                    CALL write_hdf5_attr('column', "component")
                    CALL write_hdf5_attr('row', "Rx")
                    CALL write_hdf5_attr('comment',  "complex values sorted by real/imag pairs ")
          

                    ! Write the standard deviation for the Tipper
                    CALL h5screate_simple_f(2, dim2d, dspace_id, hdferr)
                    CALL h5dcreate_f(group_id, "std", H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
                    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, allData%d(iTx)%data(ii)%error, dim2d, hdferr) ! this is the correct way to do this 
                  

                    ! Write Tipper error value 
                    CALL h5screate_simple_f(2, dim2d, dspace_id, hdferr)
                    CALL h5dcreate_f(group_id, "value", H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
                    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, allData%d(iTx)%data(ii)%value, dim2d, hdferr)
            

                    ! create the dataset for tipper index for codes ! 
                    CALL h5screate_simple_f(1, dim1d, dspace_id, hdferr)
                    CALL h5dcreate_f(group_id, "irx", H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferr)
                    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, allData%d(iTx)%data(ii)%rx, dim1d, hdferr)
                
                  
        
                case(Full_Impedance)
                    ncomp = allData%d(iTx)%data(ii)%nComp
                    nRx = allData%d(iTx)%data(ii)%nSite
                    dim2d = (/ncomp,nRx/)
                    dim1d = (/nRx/)  !set up dimension size ! 
                    dbTZ = trim(dblk_num)//trim(Z)
             
                    CALL h5gcreate_f(file_id, dbTZ, group_id, hdferr)
                  
                    ! Create attributes for Impedance
                    CALL write_hdf5_attr('column', 'component')
                    CALL write_hdf5_attr('row', 'Rx')
                    CALL write_hdf5_attr('comment',  "complex values sorted by real/imag pairs ")

                    !Write the standard deviation for the Impedance
                    CALL h5screate_simple_f(2, dim2d, dspace_id, hdferr)
                    CALL h5dcreate_f(group_id, "std", H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
                    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, allData%d(iTx)%data(ii)%error, dim2d, hdferr)
                   
                    ! Write the error values for Impedance 
                    CALL h5screate_simple_f(2, dim2d, dspace_id, hdferr)
                    CALL h5dcreate_f(group_id, "value", H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
                    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,allData%d(iTx)%data(ii)%value, dim2d, hdferr)
                    

                    CALL h5screate_simple_f(1, dim1d, dspace_id, hdferr)
                    CALL h5dcreate_f(group_id, "irx", H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferr)
                    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, allData%d(iTx)%data(ii)%rx, dim1d, hdferr)
                    
                  
                end select
            end do WRITE_DATA_TYPE
    end do

    CALL close_hdf5(cfile)

end subroutine write_hdf5_datablocks

!**********************************************************************
subroutine read_hdf5_txdict(cfile)
    character(*), intent(in)                  :: cfile

    CHARACTER(LEN=9), PARAMETER  :: data_group = "Data"
    CHARACTER(LEN=9), PARAMETER  :: data_mt_group = "Data/MT"
    CHARACTER(LEN=18), parameter :: data_mt_txdict_group = "Data/MT/txdict"

    ! DATASET DIMENSIONS USED FOR HDF5 
    INTEGER(HSIZE_T), DIMENSION(1) :: dim1d ! Datasets dimensions for 1D arrays
    INTEGER(HSIZE_T), DIMENSION(2) :: dim2d ! Datasets dimensions for T and Z in datablocks 
    INTEGER(HSIZE_T), DIMENSION(1) :: dimsc = 1 ! Scalar or single value string
    ! INTEGER(HSIZE_T)  :: attrlen = 100 

    REAL(KIND=8), DIMENSION(1), allocatable :: rdata(:) !Read data buffer
    integer                        :: i, hdferr, istat, nTx
    INTEGER(HSIZE_T)                :: npoints
   CHARACTER (len=100), Dimension(1) :: att_data !attribute scalar buffer
    
    ! ALLOCATE(rdata(size(txDict)))
  
    write(0,*) 'Reading Transmitter Dictionary'
    CALL open_read_hdf5(cfile) 
    write(0,*) cfile,' is open and ready to read the transmitter dictionary'
 
    CALL h5gopen_f(file_id, data_mt_txdict_group, group_id, hdferr)
    CALL h5dopen_f (group_id, "periods", dset_id, hdferr)

    call h5dget_space_f(dset_id, dspace_id, hdferr)
    call h5sget_simple_extent_npoints_f(dspace_id, npoints, hdferr)
    write(0,*) 'read ',npoints,' periods from file'
    
    !allocate the space for the local variable rdata
    allocate(rdata(npoints))

    !define the number of elements in the array for set_up txdict 
    nTx = npoints

    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE,rdata, dim1d, hdferr)


!    att_data = read_hdf5_attr('order') ! change this from subroutine call everywhere


    write(0,*) 'Transmitter Attribute ', att_data

    CALL close_hdf5(cfile)

    call setup_txDict(nTx,rdata,2) 
    call print_txDict()
    !figure out how to grab nTx from the already written hdf5 dataset

end subroutine read_hdf5_txdict

! !**********************************************************************
!     !READ HDF5 
subroutine read_hdf5_rxdict(cfile)
    character(*), intent(in)                  :: cfile
 
  
    CHARACTER(LEN=9), PARAMETER  :: data_group = "Data"
    CHARACTER(LEN=9), PARAMETER  :: data_mt_group = "Data/MT"
    CHARACTER(LEN=18), parameter :: data_mt_rxdict_group = "Data/MT/rxdict"
  
    ! DATASET DIMENSIONS USED FOR HDF5 
    INTEGER(HSIZE_T), DIMENSION(1) :: dim1d ! Datasets dimensions for 1D arrays
    INTEGER(HSIZE_T), DIMENSION(2) :: dim2d ! Datasets dimensions for T and Z in datablocks 
    INTEGER(HSIZE_T)               :: npoints

    REAL(KIND=8), DIMENSION(1), allocatable :: elv(:), lat(:), lon(:) !Read buffers for siteId
    REAL(KIND=8), DIMENSION(1), allocatable :: siteLocations(:,:) !local variable for storing three comp site locations
    REAL(KIND=8), DIMENSION(2), allocatable :: r2data(:,:) !Read buffer
    CHARACTER(len =5), DIMENSION(1), allocatable :: codes_data(:) !Read buffer
    CHARACTER(len=100), Dimension(1), allocatable :: att_xyz(:)
    integer                        :: ii, hdferr,rz, istat, nSites, attr_num
    INTEGER(SIZE_T)                :: size


    write(0,*) 'Reading Reciever Dictionary'
    CALL open_read_hdf5(cfile) 

    CALL h5gopen_f(file_id, data_mt_rxdict_group, group_id, hdferr)
   
    !Open Elevation dataset 
    CALL h5dopen_f (group_id, "elv", dset_id, hdferr)

    ! get the number of elements, npoints, in reciever dictionary 
    call h5dget_space_f(dset_id, dspace_id, hdferr)
    call h5sget_simple_extent_npoints_f(dspace_id, npoints, hdferr)

    !Use npoints to allocate the size for all read buffer variables
    allocate(elv(npoints), lat(npoints), lon(npoints),codes_data(npoints), STAT=istat)
    allocate(r2data(3, npoints), STAT=istat)

    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, elv, dim1d, hdferr)

    CALL h5dopen_f (group_id, "lat", dset_id, hdferr)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, lat, dim1d, hdferr)

    write(0,*) 'Latitude:', lat

    CALL h5dopen_f (group_id, "lon", dset_id, hdferr)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, lon, dim1d, hdferr)

    write(0,*) 'Longitude:', lon

    CALL h5dopen_f (group_id, "codes", dset_id, hdferr)

    CALL H5Dget_type_f(dset_id, atype_id, hdferr)
    CALL H5Tget_size_f(atype_id, size, hdferr)
    CALL H5Tcopy_f(H5T_FORTRAN_S1, atype_id, hdferr)
    CALL H5Tset_size_f(atype_id, size, hdferr)

    CALL h5dread_f(dset_id, atype_id, codes_data, dim1d, hdferr)
    write(0,*) 'Reading Codes'

    CALL h5dopen_f (group_id, "xyz", dset_id, hdferr)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, r2data, dim2d, hdferr)

    ! Get the number of attributes from dataset 
    ! CALL h5Aget_num_attrs_f(dset_id, attr_num, hdferr)

    ! Read attributes and save to a local variable
    ! allocate(att_xyz(attr_num), STAT = istat)
    ! att_xyz = [read_hdf5_attr('angleunits', dset_id), read_hdf5_attr('falseeasting', dset_id), read_hdf5_attr('falsenorthing', dset_id), & 
    !     read_hdf5_attr('fixedorient', dset_id),read_hdf5_attr('geoid', dset_id), read_hdf5_attr('maplatlimit', dset_id),read_hdf5_attr('maplonlimit', dset_id),& 
    !     read_hdf5_attr('mapprojection', dset_id),read_hdf5_attr('nparalels', dset_id),read_hdf5_attr('origin', dset_id), &
    !     read_hdf5_attr('scalefactor', dset_id),read_hdf5_attr('userinfo', dset_id),read_hdf5_attr('zone', dset_id)]
  

     !allocate the space for the local variable rdata
    allocate(siteLocations(npoints,3),STAT=istat)
    do ii= 1, npoints
        siteLocations(ii,3) = elv(ii)
        siteLocations(ii,2) = lon(ii)
        siteLocations(ii,1) = lat(ii)
    end do
    
    nSites = npoints ! redfine this as int type 4 and not int8
    call setup_rxDict(nSites,siteLocations,codes_data) 
    call print_rxDict()
    deallocate(siteLocations, STAT=istat)
    deallocate(elv, lat, lon,codes_data, STAT=istat)
    deallocate(r2data, STAT=istat)
 
    Call close_hdf5(cfile)
end subroutine read_hdf5_rxdict
! !********************************************************************** 
subroutine read_hdf5_typelist(cfile)
    character(*), intent(in)                  :: cfile
  
    CHARACTER(LEN=18), parameter :: data_typelist = "Data/MT/typelist"
    CHARACTER(LEN=19), parameter :: data_typelist_T = "Data/MT/typelist/T"
    CHARACTER(LEN=19), parameter :: data_typelist_Z = "Data/MT/typelist/Z"
    CHARACTER(len =4), DIMENSION(1), allocatable :: rdata(:) !Read buffer
    
    INTEGER(HSIZE_T), DIMENSION(1) :: dim1d ! Datasets dimensions for 1D arrays
    INTEGER(HSIZE_T)               :: npoints
    integer                        :: hdferr, istat, ncomp, iDt, attr_num, i, tran_num
    INTEGER(SIZE_T)                :: str_len
    character(len=100),allocatable, dimension(1) :: att_T(:), att_Z(:)
    character(len = 30)           :: tran_name
    logical                        :: exists, tran_comp
    call open_read_hdf5(cfile) 
    CALL setup_typeDict()

    ! Read tipper definition, if exists in file
    call h5lexists_f(file_id, data_typelist_T, exists, hdferr)
    if (exists) then
        CALL h5gopen_f(file_id, data_typelist_T, group_id, hdferr)
    
        ! get the number of elements, npoints, in receiver dictionary
        CALL h5dopen_f (group_id, "components", dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        CALL h5sget_simple_extent_npoints_f(dspace_id, npoints, hdferr)
        allocate(rdata(npoints),STAT=istat)

        CALL H5Dget_type_f(dset_id, dtype_id, hdferr)
        CALL H5Tget_size_f(dtype_id, str_len, hdferr)
        CALL H5Tcopy_f(H5T_FORTRAN_S1, dtype_id, hdferr)
        CALL H5Tset_size_f(dtype_id, str_len, hdferr)
        CALL h5dread_f(dset_id, dtype_id, rdata, dim1d, hdferr)
        CALL h5dclose_f(dset_id, hdferr)

        !Begin Reading Attributes for the Group ID
        ! CALL h5gopen_f(file_id, data_typelist_T, group_id, hdferr)
        CALL h5Aget_num_attrs_f(group_id, attr_num, hdferr)

        ! Read attributes and save to a local variable
        allocate(att_T(attr_num), STAT = istat)
        
     
        att_T = [read_hdf5_attr('complex'), read_hdf5_attr('description'), read_hdf5_attr('externalurl'),read_hdf5_attr('input'),read_hdf5_attr('intention'), &
            read_hdf5_attr('longname'),read_hdf5_attr('output'),read_hdf5_attr('tag'),read_hdf5_attr('units')]
    
        tran_name = att_T(6) 
        !convert transmitter complex from attribute string to logical
        read(att_T(1),'(A1, L1)') tran_comp

    
        ! ! Now update the relevant entry in typeDict
        ! typeDict(Full_Vertical_Components)%exists = .true.

        ! ! Here, also update the units, etc with what we just read
        ! typeDict(Full_Vertical_Components)%name = att_T(2)
        ! typeDict(Full_Vertical_Components)%abbrev = 'T'
        ! typeDict(Full_Vertical_Components)%isComplex = tran_comp
        ! typeDict(Full_Vertical_Components)%tfType    = ImpType(trim(tran_name))
        ! typeDict(Full_Vertical_Components)%units     = att_T(9)
        ! typeDict(Full_Vertical_Components)%nComp     = 4

        ! allocate(typeDict(Full_Vertical_Components)%id(2),STAT=istat)
        ! typeDict(Full_Vertical_Components)%id(1)    = rdata(1)
        ! typeDict(Full_Vertical_Components)%id(2)    = rdata(2)
        
        deallocate(rdata,STAT=istat)
        deallocate(att_T, STAT = istat)
    end if

    ! Read impedance definition, if exists in file
    call h5lexists_f(file_id, data_typelist_Z, exists, hdferr)
    if (exists) then
        ! FOR TYPELIST Z
        CALL h5gopen_f(file_id, data_typelist_Z, group_id, hdferr)
    
        ! get the number of elements, npoints, in receiver dictionary
        CALL h5dopen_f (group_id, "components", dset_id, hdferr)
        call h5dget_space_f(dset_id, dspace_id, hdferr)
        CALL h5sget_simple_extent_npoints_f(dspace_id, npoints, hdferr)
        allocate(rdata(npoints),STAT=istat)

        CALL H5Dget_type_f(dset_id, dtype_id, hdferr)
        CALL H5Tget_size_f(dtype_id, str_len, hdferr)
        CALL H5Tcopy_f(H5T_FORTRAN_S1, dtype_id, hdferr)
        CALL H5Tset_size_f(dtype_id, str_len, hdferr)
        CALL h5dread_f(dset_id, dtype_id, rdata, dim1d, hdferr)

        write(0,*) 'Reading Impedance Components'

        ! Now update the relevant entry in typeDict

        CALL h5Aget_num_attrs_f(group_id, attr_num, hdferr)

        !Read Attributes and allocate to local varible
        allocate(att_Z(attr_num), STAT = istat)
        att_Z = [read_hdf5_attr('complex'), read_hdf5_attr('description'), read_hdf5_attr('externalurl'),read_hdf5_attr('input'),read_hdf5_attr('intention'), &
            read_hdf5_attr('longname'),read_hdf5_attr('output'),read_hdf5_attr('tag'),read_hdf5_attr('units')]
       
        !convert transmitter complex from attribute string to logical
        read(att_Z(1),'(A1, L1)') tran_comp
        tran_name = att_Z(6) 
        ! typeDict(Full_Impedance)%exists = .true.

        ! ! Here, also update the units, etc with what we just read
        ! typeDict(Full_Impedance)%name = att_Z(6)
        ! typeDict(Full_Impedance)%abbrev = att_Z(8)
        ! typeDict(Full_Impedance)%isComplex = tran_comp
        ! typeDict(Full_Impedance)%tfType    = ImpType(trim(tran_name))
        ! typeDict(Full_Impedance)%units     = att_Z(9)
        ! typeDict(Full_Impedance)%nComp     = 8
        ! allocate(typeDict(Full_Impedance)%id(4),STAT=istat)
        ! typeDict(Full_Impedance)%id(1)    = rdata(1)
        ! typeDict(Full_Impedance)%id(2)    = rdata(2)
        ! typeDict(Full_Impedance)%id(3)    = rdata(3)
        ! typeDict(Full_Impedance)%id(4)    = rdata(4)
        ! deallocate(rdata,STAT=istat)
        ! deallocate(att_Z, STAT = istat)
    end if
    
    call close_hdf5(cfile) 
end subroutine read_hdf5_typelist
   
! ! !**********************************************************************
subroutine read_hdf5_datablocks(allData, cfile)
    character(*), intent(in)          :: cfile
    type(dataVectorMTX_t), intent(inout) :: allData
   

    ! local
    INTEGER(HSIZE_T)               :: attrlen = 100 ! Length of the attribute string
    INTEGER(HSIZE_T), DIMENSION(1) :: dim1d ! Datasets dimensions for 1D arrays
    
    

    !DATA BLOCKS NAMES  
    character(len=27)                  :: dblk_num, dbTZ, dblk
    character(len=2)                   ::  padded_i
    REAL(8), DIMENSION(2), allocatable :: std(:,:), val(:,:) !Read buffer
    INTEGER(8), DIMENSION(2)           :: dims, maxdims !Read buffer
    INTEGER, DIMENSION(1), allocatable :: idx_data(:) !Read buffer
    INTEGER                            :: txTypeDict
    logical                            :: exists, isComplex
    ! Declare variables
    character(len=16), parameter       :: mt_group_name = "/Data/MT"
    CHARACTER(LEN=18), parameter       :: datablock = '/Data/MT/datablock'
    CHARACTER(LEN=21)                  :: datablock_group
    character(LEN=2), parameter        :: T = "/T"
    character(len=2), PARAMETER        :: Z = "/Z"
    integer                            :: nTx, nDt, nComp, nSite, dataType, hdferr, i, ii, istat, iDt, iTx, nRx, myint


    ! Open the hdf5 file
    call open_read_hdf5(cfile)
    write(0,*) cfile,' is open and ready to read the data blocks'
    
    call h5gopen_f(file_id, '/Data', group_id, hdferr)
   
    call h5lexists_f(file_id, mt_group_name, exists, hdferr)
    if (.not. exists) then
        write(0,*) 'No MT data in the file. Assume that we can only have MT data for now. Exiting...'
        stop
    end if

    ! Open the MT group if it exists
    call h5gopen_f(group_id, mt_group_name, group_id, hdferr)

    ! Create allData for MT data only!
    nTx = size(txDict)

    call create_dataVectorMTX(nTx, allData)

    do iTx = 1, nTx

        ! Name of the data block for this transmitter
        write(padded_i, '(I0.2)') iTx
        datablock_group = datablock//'.'//trim(padded_i)
        write(0,*) datablock_group

        ! Now count the number of datatypes present for this transmitter
        nDt = 0
        call h5lexists_f(file_id, datablock_group//trim(T), exists, hdferr)
        if (exists) then
            nDt = nDt + 1
        end if
        call h5lexists_f(file_id, datablock_group//trim(Z), exists, hdferr)
        if (exists) then
            nDt = nDt + 1
        end if
        write(0,*) 'Data block ',iTx,' has ',nDt,' data types'

        call create_dataVector(nDt, allData%d(iTx)) !this allocates space for two datatypes in datavector "d"
        allData%d(iTx)%tx = iTx
        allData%d(iTx)%txType = MT

        !Read datablock
        do ii = 1, nDt
            ! Open the datablock group
            call h5gopen_f(file_id, datablock_group, group_id, hdferr)
            if (ii == 2) then
                dataType = ImpType('T')
                isComplex = typeDict(dataType)%isComplex
                nComp = typeDict(dataType)%nComp
                ! Check if a T dictionary is present
                call h5lexists_f(file_id, datablock_group//trim(T), exists, hdferr)
                if (exists) then
                    call h5gopen_f(file_id, datablock_group//trim(T), group_id, hdferr)

                    ! Read the data from the T dataset
                    CALL h5dopen_f (group_id, "std", dset_id, hdferr)
                    ! get the dimensions for the datasets an allocate to the std variable
                    call h5dget_space_f(dset_id, dspace_id, hdferr)
                    CALL h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr)
                    nSite = dims(2)
                    ! write(0,*) nSite 
                    allocate(std(dims(1), dims(2)),STAT=istat)
                    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, std, dims, hdferr)
                    write(0,*) dims(1), dims(2)
                    

                    CALL h5dopen_f (group_id, "value", dset_id, hdferr)
                    ! get the dimensions for the datasets an allocate to the value variable
                    call h5dget_space_f(dset_id, dspace_id, hdferr)
                    CALL h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr)
                    allocate(val(dims(1), dims(2)),STAT=istat)
                    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, val, dims, hdferr)
                    

                    CALL h5dopen_f (group_id, "irx", dset_id, hdferr)
                    ! get the dimensions for the datasets an allocate to the value variable
                    call h5dget_space_f(dset_id, dspace_id, hdferr)
                    CALL h5sget_simple_extent_dims_f(dspace_id, dim1d, maxdims, hdferr)
                    allocate(idx_data(dim1d(1)),STAT=istat)
                    CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, idx_data, dim1d, hdferr)

                    

                    call create_dataBlock(nComp, nSite, allData%d(iTx)%data(ii), isComplex, .true.) 
                    allData%d(iTx)%data(ii)%error(:,:) = std
                    allData%d(iTx)%data(ii)%value = val
                    allData%d(iTx)%data(ii)%rx = idx_data
                    allData%d(iTx)%data(ii)%dataType = ImpType('T')
                    write(0,*) 'group ',datablock_group//trim(T),' successfully opened and closed'
                    
                    deallocate(std, STAT=istat)
                    deallocate(val, STAT=istat)
                    deallocate(idx_data, STAT=istat)
                endif
            
        ! Check if a Z dictionary is present
            else 
                dataType = ImpType('Z')
                isComplex = typeDict(dataType)%isComplex
                nComp = typeDict(dataType)%nComp
            
                call h5lexists_f(file_id, datablock_group//trim(Z), exists, hdferr)
                if (exists ) then
                    call h5gopen_f(file_id, datablock_group//trim(Z), group_id, hdferr)

                    ! Read the data from the T dataset
                    CALL h5dopen_f (group_id, "std", dset_id, hdferr)
                    ! get the dimensions for the datasets an allocate to the std variable
                    call h5dget_space_f(dset_id, dspace_id, hdferr)
                    CALL h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr)
                    allocate(std(dims(1), dims(2)),STAT=istat)
                    nSite = dims(2)
                    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, std, dims, hdferr)
                    write(0,*) dims(1), dims(2)
                    CALL h5dopen_f (group_id, "value", dset_id, hdferr)
                    ! get the dimensions for the datasets an allocate to the value variable
                    call h5dget_space_f(dset_id, dspace_id, hdferr)
                    CALL h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr)
                    allocate(val(dims(1), dims(2)),STAT=istat)
                    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, val, dims, hdferr)
                

                    CALL h5dopen_f (group_id, "irx", dset_id, hdferr)
                    ! get the dimensions for the datasets an allocate to the value variable
                    call h5dget_space_f(dset_id, dspace_id, hdferr)
                    CALL h5sget_simple_extent_dims_f(dspace_id, dim1d, maxdims, hdferr)
                    allocate(idx_data(dim1d(1)),STAT=istat)
                    CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, idx_data, dim1d, hdferr)
                  
            
                    call create_dataBlock(nComp, nSite, allData%d(iTx)%data(ii), isComplex, .true.) 
                    allData%d(iTx)%data(ii)%error = std
                    allData%d(iTx)%data(ii)%value = val
                    allData%d(iTx)%data(ii)%rx = idx_data
                    allData%d(iTx)%data(ii)%dataType = ImpType('Z')

                    write(0,*) 'group ',datablock_group//trim(Z),' successfully opened and closed'
                    write(0,*)  allData%d(iTx)%data(1)%value(1, 1)
                
                    deallocate(std, STAT=istat)
                    deallocate(val, STAT=istat)
                    deallocate(idx_data, STAT=istat)
                end if 
            end if 
        end do

        allData%d(iTx)%allocated = .true.

    end do

    allData%allocated = .true.

 CALL close_hdf5(cfile)
   
end subroutine read_hdf5_datablocks

!**********************************************************************
subroutine write_hdf5_data(allData,cfile)
    character(*), intent(in)                  :: cfile
    type(dataVectorMTX_t), intent(in)        :: allData
    
    call write_hdf5_txdict(cfile)
    call write_hdf5_rxdict(cfile)
    call write_hdf5_typelist(allData, cfile)
    call write_hdf5_datablocks(allData, cfile)
    
   
end subroutine

! !**********************************************************************
subroutine read_hdf5_data(allData, cfile)
    character(*), intent(in)                  :: cfile
    type(dataVectorMTX_t), intent(inout)      :: allData

    call read_hdf5_txdict(cfile)
    call read_hdf5_rxdict(cfile)
    call read_hdf5_typelist(cfile)
    call read_hdf5_datablocks(allData, cfile)

end subroutine read_hdf5_data
!**********************************************************************
! reads in the ASCII list data file, sets up all dictionaries
! and the allData structure, including data and error bars.
! logic here is quite complicated, but once written can be used
! to read any kind of data, by adding a new case statement.

subroutine read_Z_list(allData,cfile)

    character(*), intent(in)               :: cfile
    type(dataVectorMTX_t), intent(inout)   :: allData
    ! local variables
    type(dataVectorMTX_t)           :: newData
    integer                         :: nTx,nRx,nDt,ncomp,iRx,iTx,icomp
    integer                         :: countData,countRx
    complex(8), allocatable         :: value(:,:,:) ! (nTx,nRx,ncomp)
    real(8), allocatable            :: error(:,:,:) ! (nTx,nRx,ncomp)
    logical, allocatable            :: exist(:,:,:) ! (nTx,nRx,ncomp)
    integer, allocatable            :: new_TxType(:) ! contains txType indices (nTx)
    integer, allocatable            :: new_Tx(:) ! contains txDict indices (nTx)
    integer, allocatable            :: new_Rx(:) ! contains rxDict indices (nRx)
    character(2)                    :: temp
    character(200)                  :: txTypeName,typeName,typeInfo,typeHeader
    character(50)                   :: siteid,ref_siteid,compid
    integer                         :: nTxt,iTxt,iDt,i,j,k,istat,ios
    character(40)                   :: code,ref_code
    real(8)                         :: x(3),ref_x(3), Period,SI_factor
    real(8)                         :: lat,lon,ref_lat,ref_lon,rx_azimuth
    real(8)                         :: Zreal, Zimag, Zerr
    logical                         :: conjugate, errorBar, isComplex

    ! 2022.09.28, Liu Zhongyin, add azimu variable
    real(kind=prec), allocatable    :: HxAzimuth(:,:), ExAzimuth(:,:), HxAzimuth_ref(:,:) !(nTx,nRx)
    real(kind=prec), allocatable    :: HyAzimuth(:,:), EyAzimuth(:,:), HyAzimuth_ref(:,:) !(nTx,nRx)
    real(kind=prec)                 :: Hxangle, Exangle, Hxangle_ref, Hyangle, Eyangle, Hyangle_ref

    integer                         :: ncount
    character(1000)                 :: tmpline

    !===========================================================================
    !======================================================== New Local Variable
    !===========================================================================
    character(8)                    :: Dipole
    character(40)              	    :: Txid=''
    real(8) 			    :: Moment, Azi, Dip, LatTx, LongTx, Tx(3)
    real(8)                         :: Omega, Amplitude
    type(transmitter_t)             :: aTx

    ! First, set up the data type dictionary, if it's not in existence yet
    call setup_typeDict()

    ! Save the user preferences
    nTxt = 5
    nDt = size(typeDict)
    call init_fileInfo(nTxt,nDt)

    ! Now, read the data file
    open(unit=ioDat,file=cfile,form='formatted',status='old')
      
    ! Read the data blocks for each data type
    READ_DATA_TYPE: do
      
    	read(ioDat,'(a2,a200)',iostat=ios) temp,typeInfo
    	read(ioDat,'(a2,a200)',iostat=ios) temp,typeHeader
    	read(ioDat,'(a2,a100)',iostat=ios) temp,typeName

        ! If transmitter name exists, it precedes the typeName
        if (temp(1:1) == '+') then
            txTypeName = typeName
            read(ioDat,'(a2,a100)',iostat=ios) temp,typeName
            old_data_file_format = .false.
        else
            txTypeName = 'MT'
        end if
        iTxt = tx_type_index(txTypeName)
    	if (ios /= 0) exit
    
    	! Read new data type
    	call compact(typeName)
    	iDt = ImpType(typeName)
    	ncomp = typeDict(iDt)%nComp
    	if (typeDict(iDt)%isComplex) then
        	ncomp = ncomp/2
    	end if

        call compact(typeInfo)
    	fileInfo(iTxt,iDt)%defined = .true.
    	fileInfo(iTxt,iDt)%info_in_file = typeInfo
    	
    	! Sort out the sign convention
    	read(ioDat,'(a2,a20)',iostat=ios) temp,fileInfo(iTxt,iDt)%sign_info_in_file
    	if(index(fileInfo(iTxt,iDt)%sign_info_in_file,'-')>0) then
      		fileInfo(iTxt,iDt)%sign_in_file = - 1
    	else
      		fileInfo(iTxt,iDt)%sign_in_file = 1
    	end if
    	if (fileInfo(iTxt,iDt)%sign_in_file == ISIGN) then
      		conjugate = .false.
    	else
      		conjugate = .true.
    	end if

        read(ioDat,'(a2,a20)',iostat=ios) temp,fileInfo(iTxt,iDt)%units_in_file
        SI_factor = ImpUnits(fileInfo(iTxt,iDt)%units_in_file,typeDict(iDt)%units)

        read(ioDat,*,iostat=ios) temp,fileInfo(iTxt,iDt)%geographic_orientation
        read(ioDat,*,iostat=ios) temp,fileInfo(iTxt,iDt)%origin_in_file(1),fileInfo(iTxt,iDt)%origin_in_file(2)
        read(ioDat,*,iostat=ios) temp,nTx,nRx
        !write(0,'(a6,i5,a18,i8,a24)') 'Found ',nTx,' transmitters and ',nRx,' receivers in data block'


        if (output_level > 3) then
            write(0,*) node_info,'Reading data type: ',trim(typeName)
            write(0,*) node_info,'Sign convention in file: ',trim(fileInfo(iTxt,iDt)%sign_info_in_file)
            write(0,*) node_info,'Units in file: ',trim(fileInfo(iTxt,iDt)%units_in_file)
            write(0,*) node_info,'Number of transmitters: ',nTx
            write(0,*) node_info,'Number of receivers: ',nRx
        end if


        ! Allocate temporary data arrays
        allocate(new_TxType(nTx),new_Tx(nTx),new_Rx(nRx),STAT=istat)
        allocate(value(nTx,nRx,ncomp),error(nTx,nRx,ncomp),exist(nTx,nRx,ncomp),STAT=istat)

        ! 2022.09.28, Liu Zhongyin, add azimu allocation
        allocate(HxAzimuth(nTx,nRx),stat=istat)
        allocate(HyAzimuth(nTx,nRx),stat=istat)
        allocate(ExAzimuth(nTx,nRx),stat=istat)
        allocate(EyAzimuth(nTx,nRx),stat=istat)
        allocate(HxAzimuth_ref(nTx,nRx),stat=istat)
        allocate(HyAzimuth_ref(nTx,nRx),stat=istat)

        new_TxType(:) = 0
        new_Tx(:) = 0
        new_Rx(:) = 0
        value(:,:,:) = dcmplx(0.0d0,0.0d0)
        error(:,:,:) = LARGE_REAL
        exist(:,:,:) = .FALSE.
        countData = 0

        ! 2022.09.28, Liu Zhongyin, add azimu initial
        HxAzimuth(:,:) = R_ZERO
        HyAzimuth(:,:) = R_ZERO
        ExAzimuth(:,:) = R_ZERO
        EyAzimuth(:,:) = R_ZERO
        HxAzimuth_ref(:,:) = R_ZERO
        HyAzimuth_ref(:,:) = R_ZERO

        READ_DATA_LINE: Do

            select case (iDt)
            case(Ex_Field,Ey_Field,Bx_Field,By_Field,Bz_Field)
                if (iTxt == CSEM) then
                    read(ioDat,*,iostat=ios) Dipole, Period, Moment, Azi, Dip, Tx(1), Tx(2), Tx(3), code, x(1), x(2), x(3), compid, Zreal, Zimag, Zerr
                    if (ios /= 0 ) then
                        backspace(ioDat)
                        exit
                    end if
                    ! Find component id for this value
                    icomp = ImpComp(compid,iDt)
                    aTx%Tx_type='CSEM'
                    aTx%nPol=1
                    aTx%Dipole = Dipole
                    aTx%period = Period
                    aTx%omega = 2.0d0*PI/Period
                    aTx%xyzTx = Tx
                    aTx%azimuthTx = Azi
                    aTx%dipTx = Dip
                    atx%moment = Moment
                    !aTx%id = Txid

                else if (iTxt == TIDE) then
                    read(ioDat,*,iostat=ios) Txid, Period, Amplitude, code, lat, lon, x(1), x(2), x(3), compid, Zreal, Zimag, Zerr
                    if (ios /= 0 ) then
                        backspace(ioDat)
                        exit
                    end if
                    ! Find component id for this value
                    icomp = ImpComp(compid,iDt)
                    aTx%Tx_type='TIDE'
                    aTx%nPol=1
                    aTx%omega = 2.0d0*PI/Period
                    aTx%period = Period
                    aTx%amplitude = Amplitude
                    aTx%id = Txid
                    call compact(aTx%id)

                else
                    read(ioDat,*,iostat=ios) Period, code, lat, lon, x(1), x(2), x(3), compid, Zreal, Zimag, Zerr
                    if (ios /= 0) then
                        backspace(ioDat)
                        exit
                    end if
                    ! Find component id for this value
                    icomp = ImpComp(compid,iDt)
                    aTx%Tx_type='MT'
                    aTx%nPol=2
                    aTx%Dipole =''
                    aTx%period = Period
                    aTx%omega = 2.0d0*PI/Period

                end if

                ! Now overwrite aTx%Tx_type with txTypeName... allows for general SFF computation
                aTx%Tx_type = tx_type_name(iTxt)

                iTx = update_txDict(aTx)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
                write(siteid,'(a22,2f9.3)') code
                iRx = update_rxDict(x,siteid)

            case(Exy_Ampli_Phase)
                read(ioDat,*,iostat=ios) Dipole, Period, Moment, Azi, Dip, Tx(1), Tx(2), Tx(3), code, x(1), x(2), x(3), compid, Zreal, Zerr, rx_azimuth
                if (ios /= 0) then
                    backspace(ioDat)
                    exit
                end if
            			
                aTx%Tx_type='CSEM'
                aTx%nPol=1
                aTx%Dipole = Dipole
                aTx%period = Period
                aTx%omega = 2.0d0*PI/Period            		
                aTx%xyzTx = Tx
                aTx%azimuthTx = Azi
                aTx%dipTx = Dip
                atx%moment = Moment
						
                ! Find component id for this value
                icomp = ImpComp(compid,iDt)

                ! Update the transmitter dictionary and the index (sets up if necessary)
                iTx = update_txDict(aTx)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
                write(siteid,'(a20,2f9.3)') code,lat,lon
                iRx = update_rxDict(x,siteid,rx_azimuth)
		
            case(Exy_Field)
                read(ioDat,*,iostat=ios) Dipole, Period, Moment, Azi, Dip, Tx(1), Tx(2), Tx(3), code, x(1), x(2), x(3), compid, Zreal, Zimag, Zerr, rx_azimuth
                if (ios /= 0) then
                    backspace(ioDat)
                    exit
                end if
            			
                aTx%Tx_type='CSEM'
                aTx%nPol=1
                aTx%Dipole = Dipole
                aTx%period = Period
                aTx%omega = 2.0d0*PI/Period            		
                aTx%xyzTx = Tx
                aTx%azimuthTx = Azi
                aTx%dipTx = Dip
                atx%moment = Moment
						
                ! Find component id for this value
                icomp = ImpComp(compid,iDt)

                ! Update the transmitter dictionary and the index (sets up if necessary)
                iTx = update_txDict(aTx)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
                write(siteid,'(a20,2f9.3)') code,lat,lon
                iRx = update_rxDict(x,siteid,rx_azimuth)
		
            case(Full_Impedance,Off_Diagonal_Impedance,Full_Vertical_Components)
                read(ioDat,'(a)',iostat=ios) tmpline

                if ((ios /= 0) .or. (tmpline(1:1)=='#')) then
                    backspace(ioDat)
                    exit
                end if

                ! Liu Zhongyin, 2019.08.27, add new codes for reading data
                backspace(ioDat)
                call strcount(tmpline, ' ', ncount)
                select case (ncount)
                case(11)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr
                    Hxangle = fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = fileInfo(iTxt,iDt)%geographic_orientation
                    Hxangle_ref = 0.0
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(12)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr,Hxangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Hxangle
                    Hxangle_ref = 0.0
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(13)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr,Hxangle,Hyangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Hxangle
                    Hxangle_ref = 0.0
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(14)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr,Hxangle,Hyangle,Exangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Exangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Hxangle_ref = 0.0
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(15)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr,Hxangle,Hyangle,Exangle,Eyangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Exangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Hxangle_ref = 0.0
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Eyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Hyangle_ref = Hxangle_ref + 90.0
                case default
                    exit
                end select
		
                ! Find component id for this value
                icomp = ImpComp(compid,iDt)

                ! Update the transmitter dictionary and the index (sets up if necessary)
                aTx%Tx_type='MT'
                aTx%nPol=2
                aTx%period = Period
                aTx%omega = 2.0d0*PI/Period
                ! Now overwrite aTx%Tx_type with txTypeName... allows for general SFF computation
                aTx%Tx_type = tx_type_name(iTxt)
                iTx = update_txDict(aTx)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
                
                if (FindStr(gridCoords, CARTESIAN)>0) then
                    write(siteid,'(a20,2f9.3)') code,lat,lon
                   
                elseif (FindStr(gridCoords, SPHERICAL)>0) then
                write(siteid,'(a20,2f15.3)') code,x(1),x(2)
                    x(1) = lat
                    x(2) = lon
                end if
                iRx = update_rxDict(x,siteid)

            case(Full_Interstation_TF)
                read(ioDat,'(a)',iostat=ios) tmpline

                if ((ios /= 0) .or. (tmpline(1:1)=='#')) then
                    backspace(ioDat)
                    exit
                end if

                ! Liu Zhongyin, 2019.08.27, add new codes for reading data
                backspace(ioDat)
                call strcount(tmpline, ' ', ncount)
                select case (ncount)
                case(17)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                        ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr
                    Hxangle = fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = 0.0
                    Hxangle_ref = fileInfo(iTxt,iDt)%geographic_orientation
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(18)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                        ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr,Hxangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = 0.0
                    Hxangle_ref = Hxangle
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(19)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                        ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr,Hxangle,Hxangle_ref
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = 0.0
                    Hxangle_ref = Hxangle_ref + fileInfo(iTxt,iDt)%geographic_orientation
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(20)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                        ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr,Hxangle,Hyangle,Hxangle_ref
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = 0.0
                    Hxangle_ref = Hxangle_ref + fileInfo(iTxt,iDt)%geographic_orientation
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(21)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                        ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr,Hxangle,Hyangle,Hxangle_ref,Hyangle_ref
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = 0.0
                    Hxangle_ref = Hxangle_ref + fileInfo(iTxt,iDt)%geographic_orientation
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hyangle_ref + fileInfo(iTxt,iDt)%geographic_orientation
                case default
                    exit
                end select
		
                ! Find component id for this value
                icomp = ImpComp(compid,iDt)

                ! Update the transmitter dictionary and the index (sets up if necessary)
                aTx%Tx_type='MT'
                aTx%nPol=2
                aTx%period = Period
                aTx%omega = 2.0d0*PI/Period
                iTx = update_txDict(aTx)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
		! Note that rx_azimuth is NOT used for MT: instead, we're supporting the
		! possibility that all fields components have different azimuths (stored in allData)
                write(siteid,'(a22,2f9.3)') code,lat,lon
                write(ref_siteid,'(a22,2f9.3)') ref_code,ref_lat,ref_lon
		rx_azimuth = R_ZERO
                iRx = update_rxDict(x,siteid,rx_azimuth,ref_x,ref_siteid)


            case(Off_Diagonal_Rho_Phase,Phase_Tensor)
                read(ioDat,'(a)',iostat=ios) tmpline

                if ((ios /= 0) .or. (tmpline(1:1)=='#')) then
                    backspace(ioDat)
                    exit
                end if

                ! Liu Zhongyin, 2019.08.27, add new codes for reading data
                backspace(ioDat)
                call strcount(tmpline, ' ', ncount)
                select case (ncount)
                case(10)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr
                    Hxangle = fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = fileInfo(iTxt,iDt)%geographic_orientation
                    Hxangle_ref = 0.0
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(11)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr,Hxangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Hxangle
                    Hxangle_ref = 0.0
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(12)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr,Hxangle,Hyangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Hxangle
                    Hxangle_ref = 0.0
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(13)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr,Hxangle,Hyangle,Exangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Exangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Hxangle_ref = 0.0
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(14)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr,Hxangle,Hyangle,Exangle,Eyangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Exangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Hxangle_ref = 0.0
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Eyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Hyangle_ref = Hxangle_ref + 90.0
                case default
                    exit
                end select

                ! Find component id for this value
                icomp = ImpComp(compid,iDt)

                ! For apparent resistivities only, use log10 of the values
                if (index(compid,'RHO')>0) then
                    Zerr  = Zerr/Zreal/dlog(10.0d0) ! Propagation of error
                    Zreal = log10(Zreal)
                end if

            	! For Phase only, using radians but reading degrees [LiuZhongyin 2017.05.27]
            	if (index(compid,'PHS')>0) then
                    if (conjugate) then
                	Zreal = Zreal*D2R
                    else
                        Zreal = -Zreal*D2R
                    endif
                    Zerr  = Zerr*D2R
                end if

                ! Update the transmitter dictionary and the index (sets up if necessary)
                aTx%Tx_type='MT'
                aTx%nPol=2
                aTx%period = Period
                aTx%omega = 2.0d0*PI/Period
                iTx = update_txDict(aTx)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
                write(siteid,'(a22,2f9.3)') code,lat,lon
                iRx = update_rxDict(x,siteid)

            end select

            ! complete transmitter dictionary update
            do i = 1,nTx
                if ((new_Tx(i) == iTx) .or. (new_Tx(i) == 0)) then
                    exit
                end if
            end do
            new_Tx(i) = iTx
            new_TxType(i) = iTxt

            ! complete receiver dictionary update
            do j = 1,nRx
                if ((new_Rx(j) == iRx) .or. (new_Rx(j) == 0)) then
                    exit
                end if
            end do
            new_Rx(j) = iRx

            ! record the value for storage in the data vector
            if (typeDict(iDt)%isComplex) then
                if (conjugate) then
                    value(i,j,icomp) = SI_factor * dcmplx(Zreal,-Zimag)
                else
                    value(i,j,icomp) = SI_factor * dcmplx(Zreal,Zimag)
                end if
            else
                value(i,j,icomp) = SI_factor * Zreal
            end if
            error(i,j,icomp) = SI_factor * Zerr
            exist(i,j,icomp) = .TRUE.

            ! 2022.09.28, Liu Zhongyin, assign angle to azimu
            HxAzimuth(i,j) = Hxangle
            HyAzimuth(i,j) = Hyangle
            ExAzimuth(i,j) = Exangle
            EyAzimuth(i,j) = Eyangle
            HxAzimuth_ref(i,j) = Hxangle_ref
            HyAzimuth_ref(i,j) = Hyangle_ref

            countData = countData + 1

        end do READ_DATA_LINE

	write(0,*) 'Read ',countData,' data values of ',trim(tx_type_name(iTxt)),' type ',trim(typeDict(iDt)%name),' from file'
	call create_dataVectorMTX(nTx,newData)
	newData%allocated = .TRUE.
	errorBar = .TRUE.
        SAVE_DATA: do i = 1,nTx

	       ! Count how many receivers we really have for this transmitter
	       countRx = 0
	       do j = 1,nRx
	        if(count(exist(i,j,:))>0) then
	            countRx = countRx + 1
	        end if
	       end do

	       ! Create a data vector for this transmitter and data type
	       call create_dataVector(1,newData%d(i))
	       newData%d(i)%tx = new_Tx(i)
	       newData%d(i)%txType = new_TxType(i)
	       newData%d(i)%allocated = .TRUE.
	       call create_dataBlock(typeDict(iDt)%nComp,countRx,newData%d(i)%data(1),typeDict(iDt)%isComplex,errorBar)
	       k = 1
	       do j = 1,nRx
	           ! If no data for this receiver, skip it
	           if(count(exist(i,j,:))==0) then
	            cycle
	           end if
	           ! Otherwise, write all components to data vector
	           do icomp = 1,ncomp
	            if(typeDict(iDt)%isComplex) then
	               newData%d(i)%data(1)%value(2*icomp-1,k) = real(value(i,j,icomp))
	               newData%d(i)%data(1)%value(2*icomp  ,k) = imag(value(i,j,icomp))
	               newData%d(i)%data(1)%error(2*icomp-1,k) = error(i,j,icomp)
	               newData%d(i)%data(1)%error(2*icomp  ,k) = error(i,j,icomp)
	               newData%d(i)%data(1)%exist(2*icomp-1,k) = exist(i,j,icomp)
	               newData%d(i)%data(1)%exist(2*icomp  ,k) = exist(i,j,icomp)

	               ! 2022.09.28, Liu Zhongyin, add azimuth
	               newData%d(i)%data(1)%orient(k)%azimuth%Hx = HxAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Hy = HyAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Ex = ExAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Ey = EyAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Hx_ref = HxAzimuth_ref(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Hy_ref = HyAzimuth_ref(i,j)
	            else
	               newData%d(i)%data(1)%value(icomp,k) = real(value(i,j,icomp))
	               newData%d(i)%data(1)%error(icomp,k) = error(i,j,icomp)
	               newData%d(i)%data(1)%exist(icomp,k) = exist(i,j,icomp)

	               ! 2022.09.28, Liu Zhongyin, add azimuth
	               newData%d(i)%data(1)%orient(k)%azimuth%Hx = HxAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Hy = HyAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Ex = ExAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Ey = EyAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Hx_ref = HxAzimuth_ref(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Hy_ref = HyAzimuth_ref(i,j)
	            end if
	           end do
	           newData%d(i)%data(1)%rx(k) = new_Rx(j)
	           k = k+1
	       end do
	       newData%d(i)%data(1)%dataType = iDt
	       newData%d(i)%data(1)%tx = new_Tx(i)
	       newData%d(i)%data(1)%txType = new_TxType(i)
	       newData%d(i)%data(1)%allocated = .TRUE.

        end do SAVE_DATA

	! Merge the new data into the main data vector
	call merge_dataVectorMTX(allData,newData,allData)
    
	! 2022.09.28, Liu Zhongyin, deallocate azimu
	deallocate(HxAzimuth,stat=istat)
	deallocate(HyAzimuth,stat=istat)
	deallocate(ExAzimuth,stat=istat)
	deallocate(EyAzimuth,stat=istat)
	deallocate(HxAzimuth_ref,stat=istat)
	deallocate(HyAzimuth_ref,stat=istat)

	deallocate(value,error,exist,STAT=istat)
	deallocate(new_TxType,new_Tx,new_Rx,STAT=istat)
	call deall_dataVectorMTX(newData)

    end do READ_DATA_TYPE

    close(ioDat)

    ! Finished reading the data: write an empty line to screen
    write(0,*)

    ! Finally, set up the index vectors in the data type dictionary - used for output
    nTxt = 5
    nTx = size(txDict)
    nRx = size(rxDict)
    do iTxt = 1,nTxt
    	do iDt = 1,nDt
		allocate(fileInfo(iTxt,iDt)%tx_index(nTx),STAT=istat)
	        allocate(fileInfo(iTxt,iDt)%dt_index(nTx),STAT=istat)
	        allocate(fileInfo(iTxt,iDt)%rx_index(nTx,nRx),STAT=istat)
	        call index_dataVectorMTX(allData,iTxt,iDt,fileInfo(iTxt,iDt)%tx_index,fileInfo(iTxt,iDt)%dt_index,fileInfo(iTxt,iDt)%rx_index)
	end do
    end do

   end subroutine read_Z_list


end module DataIO
