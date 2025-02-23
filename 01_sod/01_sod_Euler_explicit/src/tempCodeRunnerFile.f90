!--------------------------------------------------------------
  ! ヤコビアン対角化行列の計算 R, R_inv, Gam, Gam_abs
  !--------------------------------------------------------------
  subroutine A_pm(ite, R, R_inv, Gam, Gam_abs)
    use globals
    implicit none
    integer, intent(in) :: ite
    real(8), intent(out) :: R(0:2,0:2), R_inv(0:2,0:2), Gam(0:2,0:2), Gam_abs(0:2,0:2)
    real(8) :: h, u, c, b_para, a_para
  
    ! エンタルピー h ： h = (E + p) / ρ  ※ e = Qc(:,2), p = qf(:,2)
    H = (Qc(ite,2) + qf(ite,2)) / Qc(ite,0)
    u = qf(ite,0)
    c = sqrt((gamma - 1.0d0) * (h - 0.5d0*u*u))
    b_para = (gamma - 1.0d0) / (c*c)
    a_para = 0.5d0 * b_para * u*u
  
    ! 左固有ベクトル
    R(0,0) = 1.0d0;    R(0,1) = 1.0d0;     R(0,2) = 1.0d0
    R(1,0) = u - c;    R(1,1) = u;         R(1,2) = u + c
    R(2,0) = H - u*c;  R(2,1) = 0.5d0*u*u;   R(2,2) = H + u*c
  
    ! 右固有ベクトル
    R_inv(0,0) = 0.5d0*(a_para + u/c)
    R_inv(0,1) = 0.5d0*(-b_para*u - 1.0d0/c)
    R_inv(0,2) = 0.5d0*b_para
    R_inv(1,0) = 1.0d0 - a_para
    R_inv(1,1) = b_para*u
    R_inv(1,2) = -b_para
    R_inv(2,0) = 0.5d0*(a_para - u/c)
    R_inv(2,1) = 0.5d0*(-b_para*u + 1.0d0/c)
    R_inv(2,2) = 0.5d0*b_para
  
    ! 対角行列 Gam（固有値）
    Gam = 0.0d0
    Gam(0,0) = u - c
    Gam(1,1) = u
    Gam(2,2) = u + c
  
    ! 対角行列 Gam_abs（固有値の絶対値）
    Gam_abs = 0.0d0
    Gam_abs(0,0) = abs(u - c)
    Gam_abs(1,1) = abs(u)
    Gam_abs(2,2) = abs(u + c)
  end subroutine A_pm