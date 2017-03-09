module extrema
  use types
  use extrema_mod
  use extrema_types, only : EXT_DATA, CND_CNTRL_TYPE
  use extrema_storage

  implicit none

  real(kind=8), allocatable, dimension(:, :, :) ::eigvec
  real(kind=8), allocatable, dimension(:, :) :: peak_pos, eigval
  real(kind=8), allocatable, dimension(:)    :: peak_values
  integer, allocatable, dimension(:)         :: peak_types, peak_index

  integer :: MOD_NPEAKS, MOD_NBINS(3), MOD_NPROCS = 1
  private :: MOD_NPEAKS, MOD_NBINS, MOD_NPROCS
contains
  !! Set the values
  subroutine set(npeaks, nbins, nthreads)
    integer, intent(in), optional :: npeaks, nbins, nthreads

    if (present(npeaks))   MOD_NPEAKS=npeaks
    if (present(nbins))    MOD_NBINS=nbins
    if (present(nthreads)) MOD_NPROCS=nthreads
  end subroutine set

  !! Call find_extrema from extrema/ files
  !! args:
  !!   - real(dp) field(N1, N2, N3) the field
  !! returns:
  !!   - type(EXT_DATA) ext(NPEAK) the peaks
  !!
  !! Note: this routine can be parallelized (but extrema_compute cannot)
  subroutine extrema_compute_ext(field, ext)
    real(kind=8), intent(in)    :: field(:, :, :)
    type(EXT_DATA), intent(inout) :: ext(:)

    type(CND_CNTRL_TYPE) :: ctrl
    real(kind=4)         :: flattened_field(size(field))

    integer :: ndim, nn(3)

    nn = shape(field)
    ndim = 3

    ctrl%nproc = MOD_NPROCS
    ctrl%justprint = .false.

    ! flatten the field
    flattened_field = reshape(real(field), (/size(field)/))

    print*, 'in extrema_compute_ext', shape(field), shape(ext)
    ! get the extrema
    call find_extrema(flattened_field, nn=nn, ext=ext, nd=ndim, ctrl=ctrl)

  end subroutine extrema_compute_ext

  !! Find the maximum of the field, and store them in the common data of the module.
  !! it is only a f2py compatible wrapper around extrema_compute_ext
  !! args:
  !!   - real(dp) field(N1, N2, N3), the field
  !! returns:
  !!   - integer ndim, npeak
  !! computes:
  !!   - real(dp) peaks(3, Npeak) the locations of the peaks
  !!   - real(dp) eigvect(3, Npeak) the eigvector of the hessian of the peak
  !!   - real(dp) eigvals(Npeak) the eigen value of the hessian of the peak
  !!   - integer peak_type(Npeak) the type of the peak
  subroutine compute(field, ndim, npeak)
    real(kind=8), intent(in) :: field(:, :, :)

    integer, intent(out) :: ndim, npeak

    type(EXT_DATA), dimension(:), allocatable :: ext

    integer :: i, ipeak

    ! Allocate the data container
    allocate(ext(mod_NPEAKS))

    ! get the extrema
    call extrema_compute_ext(field, ext)

    ! store all that in the module vars
    if (allocated(peak_pos))    deallocate(peak_pos)
    if (allocated(peak_values)) deallocate(peak_values)
    if (allocated(eigval))      deallocate(eigval)
    if (allocated(eigvec))      deallocate(eigvec)
    if (allocated(peak_types))  deallocate(peak_types)
    if (allocated(peak_index))  deallocate(peak_index)

    !TODO: remove parallel foobars
    npeak = 0
    do i = 1, MOD_NPEAKS
       if (ext(i)%typ >= 0) then
          npeak = npeak + 1
       end if
    end do

    allocate(peak_pos(3, npeak), peak_values(npeak), eigval(3, npeak), &
         eigvec(3, 3, npeak), peak_types(npeak), peak_index(npeak))

    ipeak = 1
    do i = 1, MOD_NPEAKS
       if (ext(i)%typ >= 0) then
          peak_pos(:, ipeak)= ext(i)%pos(:)
          peak_values(ipeak)= ext(i)%val
          eigval(:, ipeak)  = ext(i)%eigval(:)
          eigvec(:,:,ipeak) = ext(i)%eigvec(:, :)
          peak_types(ipeak) = ext(i)%typ
          peak_index(ipeak) = ext(i)%pix
          ipeak = ipeak + 1
       end if
    end do

    ndim = 3

    deallocate(ext)
  end subroutine compute

  !! Get the peaks
  !! args:
  !!   - integer ndim, npeak
  !! returns:
  !!   - double peakpos(ndim, npeak)    the position of the peak
  !!   - double eigvectors(ndim, npeak) the eigenvector of the hessian
  !!   - double eigvalues(npeak)        the eigenvalue of the hessian
  !!   - integer peaktypes(npeak)       the type of the peak
  !!   - integer index(npeak)           the index of the peak (in contiguous array, you should reshape it)
  subroutine get(ndim, npeak, index, peakpos, peakvals, eigvectors, eigvalues, peaktypes)
    integer, intent(in)                                     :: ndim, npeak
    real(kind=8), intent(out), dimension(ndim, npeak)       :: peakpos, eigvalues
    real(kind=8), intent(out), dimension(ndim, ndim, npeak) :: eigvectors
    real(kind=8), intent(out), dimension(npeak)             :: peakvals
    integer, intent(out), dimension(npeak)                  :: peaktypes, index
    !f2py depend(peaks)                                 :: npeak=shape(npeak, 1)

    peakpos    = peak_pos
    peakvals   = peak_values
    eigvectors = eigvec
    eigvalues  = eigval
    peaktypes  = peak_types
    index      = peak_index

  end subroutine get

end module extrema
