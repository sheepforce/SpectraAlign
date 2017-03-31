! This program takes two spectra with overlapping wavelength/energy/whatever
! ranges from different measurements and shifts the so that the overlap
! is maximal and the relative separation is minimal


program spectraaligner
  use spectra
  use grid_functions
  implicit none
  integer                                                :: line, optcriterion, combinationtype
  character*80                                           :: spectrum1_file, spectrum2_file, configpath
  integer                                                :: points_in_spectrum1, points_in_spectrum2, points_in_spectrumcombined
  double precision, dimension(:,:), allocatable          :: spectrum1, spectrum2, spectrumhigh, spectrumlow, spectrumcombined
  double precision, dimension(2)                         :: overlaprange
  integer         , dimension(2,2)                       :: overlapind
  double precision                                       :: shift
  double precision, dimension(:)  , allocatable          :: scalingtable
  
  
  ! get enviroment variable for loading the config-file
  call getenv("SA_CONFIG", configpath)
  if (configpath == "") then
    ! if no config specified default to robust choises
    optcriterion = 1
    combinationtype = 1
  else
    ! if config is available; read it
    open(unit=102, file=configpath, status="old", form="formatted")
    read(unit=102, fmt=*) optcriterion
    read(unit=102, fmt=*) combinationtype
    close(unit=102)
  end if

  ! get the filenames from the command line
  call get_command_argument(1, spectrum1_file)
  call get_command_argument(2, spectrum2_file)
  
  ! find the number of points in each spectrum
  call get_spectrum_size(spectrum1_file, points_in_spectrum1)
  call get_spectrum_size(spectrum2_file, points_in_spectrum2)
  
  ! allocate array for the two original spectra  
  allocate(spectrum1(points_in_spectrum1,2), spectrum2(points_in_spectrum2,2))
  call read_spectrum(spectrum1_file, points_in_spectrum1, spectrum1)
  call read_spectrum(spectrum2_file, points_in_spectrum2, spectrum2)
  
  ! if spectrum1 starts at the high wavelengths/energies/whatever, than this is spectrumhigh
  if (spectrum1(1,1) >= spectrum2(1,1)) then
    ! allocate memory for the new spectra
    allocate(spectrumhigh(size(spectrum1, 1),2), spectrumlow(size(spectrum2, 1),2))
    ! make the high spectrum to spectrum1
    spectrumhigh = spectrum1
    ! and make the spectrumlow to spectrum2
    spectrumlow = spectrum2
  else
    ! allocate memory for the new spectra
    allocate(spectrumhigh(size(spectrum2, 1),2), spectrumlow(size(spectrum1, 1),2))
    ! make the high spectrum to spectrum2
    spectrumhigh = spectrum2
    ! and make the spectrumlow to spectrum1
    spectrumlow = spectrum1
  end if
  ! spectrum1 and spectrum2 are not needed any longer
  deallocate(spectrum1, spectrum2)
    
  ! now find the overlap in wavelength/energies/whatever and the indizes
  ! in spectrumhigh and spectrumlow where the overlaps start and end
  call find_overvap(spectrumhigh, spectrumlow, overlaprange, overlapind)
  
  ! make a new grid in the range of the overlap of spectrumlow which matches the grid
  ! of spectrum high by interpolating values
  do line = overlapind(2,1), overlapind(2,2)
    ! beginning from the overlap point in spectrumlow start to iterate over spectrumlow
    ! and use the x values from spectrumhigh as new x values in this range fro spectrumlow
    spectrumlow(line,1) = spectrumhigh(overlapind(1,1) + line - 1,1)
    ! now set y values/intensities
    spectrumlow(line,2) = interpolate_points(spectrumlow(line - 1,:), spectrumlow(line,:), spectrumlow(line,1))
  end do
  
  ! now optimize the shift of both spectra
  ! initialize the shift
  shift=0.0
  if (optcriterion == 1) then
    ! use the average shift between all overlaping grid points
    ! add up the shift between all points in range
    do line = 1, (overlapind(1,2) - overlapind(1,1)) + 1
      shift = spectrumhigh(overlapind(1,1) + line - 1,2) - spectrumlow(line,2) + shift
    end do
    ! divide by number of points
    shift = shift / (overlapind(1,2) - overlapind(1,1))
  else if (optcriterion == 2) then
    ! optimize the lower boundary of the spectra
    shift = spectrumhigh(overlapind(1,2),2) - spectrumlow(overlapind(2,2),2)
  else if (optcriterion == 3) then
    ! optimize the upper boundary of the spectra
    shift = spectrumhigh(overlapind(1,1),2) - spectrumlow(overlapind(2,1),2)
  else if (optcriterion == 4) then
    ! optimize the average of lower and upper boundary and ignore in between
    shift = (spectrumhigh(overlapind(1,2),2) - spectrumlow(overlapind(2,2),2) &
           + spectrumhigh(overlapind(1,1),2) - spectrumlow(overlapind(2,1),2)) / 2
  end if
  
  !shift spectrumhigh down by half the calculated shift and spectrumlow up by shift
  spectrumhigh(:,2) = spectrumhigh(:,2) - (shift / 2)
  spectrumlow(:,2) = spectrumlow(:,2) + (shift / 2)
  
  ! combine the shifted spectra
  ! allocate memory for scaling table
  allocate(scalingtable(overlapind(1,2) - overlapind(1,1)))
  ! calculate the scalingtable
  if (combinationtype == 1) then
    ! weigthed linear combination. 1 at the ends where the spectrum continues and 0 where it ends, 0.5 in the middle
    do line = 1, (overlapind(1,2) - overlapind(1,1)) + 1
      scalingtable(line) = dble(line - 1) / dble((overlapind(1,2) - overlapind(1,1)))
    end do
  else if (combinationtype == 2) then
    ! 50:50 mixing of both
    scalingtable = 0.5
  end if
  
  ! allocate memory for the output spectrum (spectrumcombined)
  points_in_spectrumcombined = size(spectrumhigh,1) + (size(spectrumlow, 1) - overlapind(2,2))
  allocate(spectrumcombined(points_in_spectrumcombined,2))
  ! initialize spectrumcombined
  spectrumcombined=0.0
  ! the first part is the shifted spectrumhigh
  do line = 1, overlapind(1,1)
    spectrumcombined(line,:) = spectrumhigh(line,:)
  end do
  ! the middle part is the combined spectrum of both
  do line = overlapind(1,1) + 1, overlapind(1,2)
    spectrumcombined(line,1) = spectrumhigh(line,1)
    spectrumcombined(line,2) = scalingtable(line - overlapind(1,1)) * spectrumlow(line - overlapind(1,1) ,2) &
                             + (1 - scalingtable(line - overlapind(1,1))) * spectrumhigh(line,2)
  end do  
  ! the last part is the shift spectrumlow
  do line = overlapind(1,2) + 1, points_in_spectrumcombined
    spectrumcombined(line,:) = spectrumlow(line - overlapind(1,2) + overlapind(2,2),:)
  end do
  
  deallocate(spectrumhigh, spectrumlow, scalingtable)
  
  do line = 1, points_in_spectrumcombined
    !write(*, *) spectrumcombined(line,:)
  end do
  
  
end program spectraaligner
