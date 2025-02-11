!----------------------------------------------------
  ! 状態量計算の各関数
  !----------------------------------------------------

  function cal_p3(m, g, p1) result(ans)
    real(kind=8), intent(in) :: m, g, p1
    real(kind=8) :: temp, ans
    temp = 1.0d0 + (2.0d0 * g * (m**2 - 1.0d0))/(g+1.0d0)
    ans = temp * p1
  end function cal_p3

  function cal_p2(p3) result(ans)
    real(kind=8), intent(in) :: p3
    real(kind=8) :: ans
    ans = p3
  end function cal_p2

  function cal_rho3(m, p3, g, a4, c) result(ans)
    real(kind=8), intent(in) :: m, p3, g, a4, c
    real(kind=8) :: temp, a3_local
    temp = 1.0d0 - c * ( m - 1.0d0/m )
    a3_local = a4 * temp
    ans = g * p3 / (a3_local**2)
  end function cal_rho3

  function cal_rho2(m, g, rho1) result(ans)
    real(kind=8), intent(in) :: m, g, rho1
    real(kind=8) :: temp, ans
    temp = (g+1.0d0)*m**2 / ( 2.0d0 + (g-1.0d0)*m**2 )
    ans = rho1 * temp
  end function cal_rho2

  function cal_u2(m, a1, g) result(ans)
    real(kind=8), intent(in) :: m, a1, g
    real(kind=8) :: temp, ans
    temp = 2.0d0 * a1 / (g+1.0d0)
    ans = temp * ( m - 1.0d0/m )
  end function cal_u2

  ! cal_u3は2つの値（u3 と a3）を返すためサブルーチンとする
  subroutine cal_u3(m, a4, c, g, u3, a3_out)
    real(kind=8), intent(in) :: m, a4, c, g
    real(kind=8), intent(out) :: u3, a3_out
    real(kind=8) :: temp, a3_local
    temp = 1.0d0 - c*(m - 1.0d0/m)
    a3_local = a4 * temp
    temp = a4/a3_local - 1.0d0
    u3 = temp * 2.0d0 * a3_local / (g-1.0d0)
    a3_out = a3_local
  end subroutine cal_u3

  function fx(m, g, a, b, c, alpha) result(ans)