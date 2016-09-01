MODULE randomxtra
  ! a module to containt various array-based random number routines
  ! vector binom deviates
  use random
  implicit none

CONTAINS

  ! 4 array of poisson deviates
  function rporay4( rate, dms )
    ! rates and dimensions
    implicit none
    integer, dimension(4), intent(in) :: dms
    real, dimension(dms(1),dms(2),dms(3),dms(4)), intent(in) :: rate
    integer, dimension(dms(1),dms(2),dms(3),dms(4)) :: rporay4 
    integer :: i1,i2,i3,i4
    do i1=1,dms(1)
       do i2=1,dms(2)
          do i3=1,dms(3)
             do i4=1,dms(4) 
                rporay4(i1,i2,i3,i4) = random_Poisson2(rate(i1,i2,i3,i4))
             enddo
          enddo
       enddo
    enddo
  end function rporay4


  function rbiray4s( nz, prob, dms )
    ! array of numbers n, single probability, dimensions
    implicit none
    integer, dimension(4), intent(in) :: dms
    integer, dimension(dms(1),dms(2),dms(3),dms(4)), intent(in) :: nz
    integer, dimension(dms(1),dms(2),dms(3),dms(4)) :: rbiray4s
    integer :: i1,i2,i3,i4
    real :: prob
    logical :: first = .TRUE.
    do i1=1,dms(1)
       do i2=1,dms(2)
          do i3=1,dms(3)
             do i4=1,dms(4) 
                rbiray4s(i1,i2,i3,i4) = random_binomial2(nz(i1,i2,i3,i4),prob,first)
             enddo
          enddo
       enddo
    enddo
  end function  rbiray4s

  
  function rbiray4( nz, prob, dms )
    ! array of numbers n, single probability, dimensions
    implicit none
    integer, dimension(4), intent(in) :: dms
    integer, dimension(dms(1),dms(2),dms(3),dms(4)), intent(in) :: nz
    integer, dimension(dms(1),dms(2),dms(3),dms(4)) :: rbiray4
    integer :: i1,i2,i3,i4
    real, dimension(dms(1),dms(2),dms(3),dms(4)) :: prob
    logical :: first = .TRUE.
    do i1=1,dms(1)
       do i2=1,dms(2)
          do i3=1,dms(3)
             do i4=1,dms(4) 
                rbiray4(i1,i2,i3,i4) = random_binomial2(nz(i1,i2,i3,i4),prob(i1,i2,i3,i4),first)
             enddo
          enddo
       enddo
    enddo
  end function  rbiray4


END MODULE randomxtra
