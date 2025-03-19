program Boltz_plot

    use utilities
    implicit none

    integer, parameter :: num_samples = 1000000
    real(8), parameter :: bin_min = 0, bin_max = 7.0
    real(8), parameter :: bin_width = 0.2
    real*8, parameter :: kbT=50

    integer :: i, bin_index, num_bins, j
    real(8) :: x
    real(8), allocatable :: histogram(:)

    ! Determine number of bins
    num_bins = int((bin_max - bin_min) / bin_width) + 1
    allocate(histogram(num_bins))

    ! Populate bins
    histogram = 0.0
    do i = 1, num_samples
        x = 0.0
        do j=1,3
            x = x + sqrt(kbT)*gaussian_distr()**2
        enddo
        x = sqrt(x)/2.0

        ! bin index of x
        bin_index = int((x - bin_min) / bin_width) + 1
        
        ! Only count if in range
        if (bin_index >= 1 .and. bin_index <= num_bins) then
            histogram(bin_index) = histogram(bin_index) + 1
        end if
    end do

    ! Normalize histogram
    histogram = histogram / num_samples

    ! Print
    print *, "# bin_center Freq"
    do i = 1, num_bins
        print *, bin_min + (i - 1) * bin_width + bin_width / 2.0, histogram(i)
    end do

    deallocate(histogram)

end program Boltz_plot