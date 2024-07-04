!> Routines to compute the integral of the HSE singularity with 
!> differene methods. The integral around the singularity is 
!> constructed employing the isotropic average method. For more 
!> datails look at Vona <i>et al.</i>, 
!> <i>Adv. Theory Simul.</i>, <b>5</b>, 2100496 (2022). 
module hse_singularity
  Use precision, only: dp 
  implicit none
  private
  public :: hse_singularity_exact_solution, hse_singularity_Taylor_expansion
 
  contains
  !> HSE singularity integral solved exaclty. The exact solution 
  !> of the singularity integral has the following form:
  !> \[
  !>   I_{s}= \frac{16 \pi^2}{V_{\mathbf{k}}}\Bigg(R_{\mathbf{k}}
  !>          -\sqrt{\pi \omega^2}\mathrm{erf}\Big(\sqrt{\frac{1}
  !>          {4\omega^2}}R_{\mathbf{k}}\Big)\Bigg). 
  !> \]
  !> In the equation \(\omega\) is the HSE screening parameter,
  !> \(V_{\mathbf{k}}\) is the spherical volume of radius 
  !> \(R_{\mathbf{k}}\) around the singularity. In particular:
  !> \[ 
  !>    V_{\mathbf{k}}=\frac{\Omega_{\mathrm{BZ}}}{N_{\mathbf{k}}}\;\mathrm{and}
  !>    \; R_{\mathbf{k}}=\Bigg(\frac{3 \Omega_{\mathrm{BZ}}}{4 \pi 
  !>    N_\mathbf{k}}\Bigg)^{1/3},
  !> \]
  !> in which \(\Omega_{\mathrm{BZ}}\) is the Brillouin zone volume and 
  !> \(N_{\mathbf{k}}\) the numer of \(\mathbf{k}\)-points.
  function hse_singularity_exact_solution(omega_hse, volume_unit_cell, number_k_points)
    use constants, only: pi,twopi
    implicit none
    real(dp) :: hse_singularity_exact_solution 
    !> HSE screening parameter
    real(dp) :: omega_hse
    real(dp) :: volume_unit_cell
    integer :: number_k_points
    real(dp) :: omega2
    real(dp) :: BZ_volume
    real(dp) :: vk
    real(dp) :: rk
    real(dp) :: erf_arg
    real(dp) :: parenthesis_arg
    real(dp), parameter :: pi2=pi*pi
 
    omega2 = omega_hse * omega_hse
    BZ_volume = twopi**3 / volume_unit_cell
    vk = BZ_volume / dble(number_k_points)
    rk = ((3._dp * vk) / (4._dp * pi))**(1._dp / 3._dp)
    erf_arg = sqrt(1._dp / (4._dp * omega2)) * rk
    parenthesis_arg = rk - sqrt(omega2 * pi) * erf(erf_arg) 
    hse_singularity_exact_solution = (16._dp * pi2 * parenthesis_arg) / vk
  end function 

  !> HSE singularity integral solved after Taylor expanded its argument. 
  !> By employing the Taylor expansion method the singularity integral
  !> has the following form:
  !> \[
  !>    I_{s}={\frac{\pi}{\omega^{2}}},
  !> \]
  !> in which \(\omega\) is the HSE screening parameter.
  function hse_singularity_Taylor_expansion(omega_hse)
    use constants, only: pi
    implicit none
    real(dp) :: hse_singularity_Taylor_expansion
    !> HSE screening parameter
    real(dp) :: omega_hse
    hse_singularity_Taylor_expansion = pi / (omega_hse**2) 
  end function 
  

end module
