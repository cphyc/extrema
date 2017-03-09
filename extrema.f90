module extrema
  use types
  use extrema_mod
  use extrema_types, only : EXT_DATA, CND_CNTRL_TYPE
  use extrema_storage

  implicit none

  real(kind=8), allocatable, dimension(:, :, :) ::mod_eigvec
  real(kind=8), allocatable, dimension(:, :) :: mod_peaks, mod_eigval
  real(kind=8), allocatable, dimension(:)    :: mod_value
  integer, allocatable, dimension(:)         :: mod_peak_type, mod_index

  integer :: MOD_NPEAKS, MOD_NBINS(3), MOD_NPROCS = 1
  private :: MOD_NPEAKS, MOD_NBINS, MOD_NPROCS
contains
  !! Set the values
  subroutine extrema_set(npeaks, nbins, nthreads)
    integer, intent(in), optional :: npeaks, nbins, nthreads

    if (present(npeaks))   MOD_NPEAKS=npeaks
    if (present(nbins))    MOD_NBINS=NBINS
    if (present(nthreads)) MOD_NPROCS=nthreads
  end subroutine extrema_set

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
  !!   - real(dp) mod_peaks(3, Npeak) the locations of the peaks
  !!   - real(dp) mod_eigvect(3, Npeak) the eigvector of the hessian of the peak
  !!   - real(dp) mod_eigvals(Npeak) the eigen value of the hessian of the peak
  !!   - integer mod_peak_type(Npeak) the type of the peak
  subroutine extrema_compute(field, ndim, npeak)
    real(kind=8), intent(in) :: field(:, :, :)

    integer, intent(out) :: ndim, npeak

    type(EXT_DATA), dimension(:), allocatable :: ext

    integer :: i, ipeak

    ! Allocate the data container
    allocate(ext(MOD_NPEAKS))

    ! get the extrema
    call extrema_compute_ext(field, ext)

    ! store all that in the module vars
    if (allocated(mod_peaks))     deallocate(mod_peaks)
    if (allocated(mod_value))     deallocate(mod_value)
    if (allocated(mod_eigval))    deallocate(mod_eigval)
    if (allocated(mod_eigvec))    deallocate(mod_eigvec)
    if (allocated(mod_peak_type)) deallocate(mod_peak_type)
    if (allocated(mod_index))     deallocate(mod_index)

    !TODO: remove parallel foobars
    npeak = 0
    do i = 1, MOD_NPEAKS
       if (ext(i)%typ >= 0) then
          npeak = npeak + 1
       end if
    end do

    allocate(mod_peaks(3, npeak), mod_value(npeak), mod_eigval(3, npeak), &
         mod_eigvec(3, 3, npeak), mod_peak_type(npeak), mod_index(npeak))

    ipeak = 1
    do i = 1, MOD_NPEAKS
       if (ext(i)%typ >= 0) then
          mod_peaks(:, ipeak)   = ext(i)%pos(:)
          mod_value(ipeak)      = ext(i)%val
          mod_eigval(:, ipeak)  = ext(i)%eigval(:)
          mod_eigvec(:,:,ipeak) = ext(i)%eigvec(:, :)
          mod_peak_type(ipeak)  = ext(i)%typ
          mod_index(ipeak)      = ext(i)%pix
          ipeak = ipeak + 1
       end if
    end do

    ndim = 3

    deallocate(ext)
  end subroutine extrema_compute

  !! Get the peaks at some stepx
  !! args:
  !!   - integer ndim, npeak
  !! returns:
  !!   - double peakpos(ndim, npeak)    the position of the peak
  !!   - double eigvectors(ndim, npeak) the eigenvector of the hessian
  !!   - double eigvalue(npeak)         the eigenvalue of the hessian
  !!   - integer peaktype(npeak)        the type of the peak
  !!   - integer index(npeak)           the index of the peak (in contiguous array, you should reshape it)
  subroutine extrema_get_i(ndim, npeak, istep, index, peakpos, eigvectors, eigvalues, peaktype)
    integer, intent(in) :: ndim, npeak, istep
    real(kind=8), intent(out), dimension(ndim, npeak) :: peakpos
    real(kind=8), intent(out), dimension(ndim, npeak) :: eigvectors
    real(kind=8), intent(out), dimension(npeak)       :: eigvalues
    integer, intent(out), dimension(npeak)            :: peaktype, index

    integer :: i, ipeak
    type(EXT_DATA), dimension(:), allocatable :: ext

    call extrema_pull(ext, istep)

    ipeak = 1
    do i = 1, MOD_NPEAKS
       if (ext(i)%typ >= 0) then
          mod_peaks(:, ipeak)   = ext(i)%pos(:)
          mod_value(ipeak)      = ext(i)%val
          mod_eigval(:, ipeak)  = ext(i)%eigval(:)
          mod_eigvec(:,:,ipeak) = ext(i)%eigvec(:, :)
          mod_peak_type(ipeak)  = ext(i)%typ
          mod_index(ipeak)      = ext(i)%pix
          ipeak = ipeak + 1
       end if
    end do
  end subroutine extrema_get_i


  !! Get the peaks
  !! args:
  !!   - integer ndim, npeak
  !! returns:
  !!   - double peakpos(ndim, npeak)    the position of the peak
  !!   - double eigvectors(ndim, npeak) the eigenvector of the hessian
  !!   - double eigvalues(npeak)        the eigenvalue of the hessian
  !!   - integer peaktypes(npeak)       the type of the peak
  !!   - integer index(npeak)           the index of the peak (in contiguous array, you should reshape it)
  subroutine extrema_get(ndim, npeak, index, peakpos, peakvals, eigvectors, eigvalues, peaktypes)
    integer, intent(in)                                     :: ndim, npeak
    real(kind=8), intent(out), dimension(ndim, npeak)       :: peakpos, eigvalues
    real(kind=8), intent(out), dimension(ndim, ndim, npeak) :: eigvectors
    real(kind=8), intent(out), dimension(npeak)             :: peakvals
    integer, intent(out), dimension(npeak)                  :: peaktypes, index
    !f2py depend(mod_peaks)                                 :: npeak=shape(npeak, 1)

    peakpos    = mod_peaks
    peakvals   = mod_value
    eigvectors = mod_eigvec
    eigvalues  = mod_eigval
    peaktypes  = mod_peak_type
    index      = mod_index

  end subroutine extrema_get

end module extrema
