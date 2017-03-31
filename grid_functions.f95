module grid_functions
  implicit none
  
  contains
  ! get interpolated value between two point by linear function between this two points
  ! y = m * x + n
  ! m = (y2 - y1) / (x2 - x1)
  double precision function interpolate_points(point1, point2, x)
    implicit none
    double precision, dimension(2)  , intent(in)           :: point1, point2
    double precision                , intent(in)           :: x
    double precision                                       :: m, n
    
    m = (point2(2) - point1(2)) / (point2(1) - point1(2))
    n = point2(2) - m * point2(1)
    
    interpolate_points = m * x + n
  end function interpolate_points
end module grid_functions
