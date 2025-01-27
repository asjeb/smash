!%      (MD) Module Differentiated.
!%
!%      Variables
!%      ---------
!%
!%      - sp
!%      - dp
!%      - lchar
!%      - gravity
!%      - manning

module md_constant

    implicit none

    integer, parameter :: sp = 4
    integer, parameter :: dp = 8
    integer, parameter :: lchar = 128
    real(sp) :: gravity = 9.8
    real(sp) :: manning = 0.033

end module md_constant
