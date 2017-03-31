module spectra
  implicit none
  integer                                                :: max_points=1000000
  integer                                      , private :: line, lines_total
  integer                                                :: points_in_spectrum
  double precision, dimension(:,:), allocatable, private :: dummyarray
  double precision, dimension(:,:), allocatable          :: spectrum
  
  contains
  ! read the spectra from file and get the number of points
  subroutine get_spectrum_size(filename, points_in_spectrum)
    implicit none
    integer                         , intent(out)          :: points_in_spectrum
    character*80                    , intent(in)           :: filename
    integer                                                :: filestat
    
    allocate(dummyarray(max_points,2))
    open(unit=101, file=filename, status="old", form="formatted")
    
    ! iterate over every line, if it is readable; do so
    ! if it is not readable; throw an error
    ! if there is an EOL (filestat < 0) stop the loop
    do line = 1, max_points
      read(101, *, iostat=filestat) dummyarray(line,:)
      if (filestat > 0) then
        stop "an error has occured reading the spectra"
      else if (filestat < 0) then
        points_in_spectrum = line - 1
        exit
      end if
    end do
    
    deallocate(dummyarray)
    close(unit=101)
  end subroutine get_spectrum_size
  
  ! read a spectrum from file in an array
  subroutine read_spectrum(filename, points_in_spectrum, spectrum)
    implicit none
    character*80                    , intent(in)            :: filename
    integer                         , intent(in)            :: points_in_spectrum
    double precision, dimension(:,:), intent(inout)         :: spectrum
    
    open(unit=101, file=filename, status="old", form="formatted")
    
    do line = 1, points_in_spectrum
      read(unit=101, fmt=*) spectrum(line,:)
    end do
    
  end subroutine read_spectrum
  
  ! take two spectra and find their overlaprange overlap
  ! returns the index from both spectra from which they overlap as well
  ! as the overlaprange range
  ! ind(a,b) holds the first and last index of spectrumhigh in a and similiar for spectrumlow b
  subroutine find_overvap(spectrumhigh, spectrumlow, overlaprange, ind)
    implicit none
    double precision, dimension(:,:), intent(in)           :: spectrumhigh, spectrumlow
    double precision, dimension(2)  , intent(out)          :: overlaprange
    integer         , dimension(2,2), intent(out)          :: ind
    
    ! set values we already know
    ! the end of the overlap in spectrumhigh is its end
    ind(1,2)=size(spectrumhigh, 1)
    ! the start of the overlap in spectrumlow is its start
    ind(2,1)=1
    
    ! go through spectrumhigh line by line
    do line = 1, size(spectrumhigh, 1)
      ! and check if spectrumlow starts to be in the range of spectrumhigh
      if (spectrumhigh(line,1) <= spectrumlow(1,1)) then
        ! then the begin of the overlap in spectrumhigh must be in this line
        ind(1,1)=line
        exit
      end if
    end do
    
    ! got through spectrumlow line by line
    do line = 1, size(spectrumlow, 1)
      ! and check if spectrum2 drops below the lowest value of spectrumhigh
      if (spectrumlow(line,1) <= spectrumhigh(size(spectrumhigh, 1),1)) then
        ! then this must be the last line where spectrums overlap
        ind(2,2)=line
        exit
      end if
    end do
    
    ! and now, find also the overlaprange range
    overlaprange(1)=spectrumlow(1,1)
    overlaprange(2)=spectrumhigh(size(spectrumhigh, 1),1)
  end subroutine find_overvap
end module spectra
