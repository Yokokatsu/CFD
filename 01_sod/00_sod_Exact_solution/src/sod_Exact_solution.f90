!==============================================================
! shock tube Exact solution
! 2025/02/11 Yokoyama Katsuyuki
!==============================================================

program shock_tube
  implicit none
  real(kind=8) :: g, rho1, u1, p1, rho4, p4, u4
  real(kind=8) :: t, dt
  integer :: nt
  real(kind=8) :: lenx, dx, x0
  integer :: nx
  real(kind=8) :: a1, a4
  real(kind=8) :: a, b, c, alpha
  real(kind=8) :: m, m_init
  real(kind=8) :: p3, p2, rho2, rho3, u2, u3, a3, us
  integer :: i, tt
  real(kind=8) :: time, temp
  real(kind=8), allocatable :: x(:), rho(:), u(:), p(:)
  character(len=20) :: dir="./../result/"

  ! 初期値
  g    = 1.4d0
  rho1 = 0.125d0
  u1   = 0.0d0
  p1   = 0.1d0

  rho4 = 1.0d0
  p4   = 1.0d0
  u4   = 0.0d0

  nt   = 300
  dt   = 0.001d0
  t    = nt*dt   
  
  nx   = 100
  dx   = 0.01d0
  lenx = nx*dx
  x0   = 0.5d0

  a1 = sqrt(g*p1/rho1)
  a4 = sqrt(g*p4/rho4)

  ! 定数
  a     = p4/p1
  b     = 2.0d0*g/(g+1.0d0)
  c     = (g-1.0d0)/(g+1.0d0) * (a1/a4)
  alpha = 2.0d0*g/(g-1.0d0)

  ! 初期Mach数
  m_init = 5.0d0

  ! Newton法によるMach数の収束
  m = m_init
  do i = 1, 19  !20回としている
     m = m - fx(m, g, a, b, c, alpha)/fx_prime(m, g, a, b, c, alpha)
  end do

  ! 状態量の計算
  p3   = cal_p3(m, g, p1)
  p2   = cal_p2(p3)
  rho2 = cal_rho2(m, g, rho1)
  rho3 = cal_rho3(m, p3, g, a4, c)
  u2   = cal_u2(m, a1, g)
  call cal_u3(m, a4, c, g, u3, a3)
  us   = m * a1

  ! 空間配列の確保と初期化（各セルの中心座標）
  allocate(x(nx), rho(nx), u(nx), p(nx))
  do i = 1, nx
     x(i) = dx/2.0d0 + dx*(i-1)
  end do

  ! 時間発展ループ
  do tt = 0, nt-1
     time = dt * tt
     ! 各セルごとの状態量を計算
     do i = 1, nx
        if ( x0 + us*time < x(i) ) then
           ! 領域1
           rho(i) = rho1
           u(i)   = u1
           p(i)   = p1
        else if ( (x0 + u2*time < x(i)) .and. (x(i) <= x0 + us*time) ) then
           ! 領域2
           rho(i) = rho2
           u(i)   = u2
           p(i)   = p2
        else if ( (x0 + (u3 - a3)*time < x(i)) .and. (x(i) <= x0 + u2*time) ) then
           ! 領域3
           rho(i) = rho3
           u(i)   = u3
           p(i)   = p3
        else if ( (x0 - a4*time < x(i)) .and. (x(i) <= x0 + (u3 - a3)*time) ) then
           ! 膨張波領域 (time=0 のときは選ばれない）
           temp   = a4 * ( 2.0d0/(g+1.0d0) - (g-1.0d0)/(g+1.0d0) * ((x(i)-x0)/(a4*time)) )
           rho(i) = rho4 * ( (temp/a4)**(2.0d0/(g-1.0d0)) )
           u(i)   = 2.0d0*a4/(g+1.0d0) * ( 1.0d0 + ((x(i)-x0)/(a4*time)) )
           p(i)   = p4 * ( (temp/a4)**(2.0d0*g/(g-1.0d0)) )
        else
           ! 領域4
           rho(i) = rho4
           u(i)   = u4
           p(i)   = p4
        end if
     end do

    ! 結果出力
     call write_data(trim(dir), tt, nx, x, rho, u, p, time)

     print *, 'Time step ', tt+1, '/', nt
  end do

  deallocate(x, rho, u, p)

contains

  !----------------------------------------------------
  ! 状態量計算の関数
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
    real(kind=8) :: temp, a3_local,ans
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
    real(kind=8), intent(in) :: m, g, a, b, c, alpha
    real(kind=8) :: ans
    ans = b*(m**2 - 1.0d0) + 1.0d0 - a*( 1.0d0 - c*( m - 1.0d0/m ) )**alpha
  end function fx

  function fx_prime(m, g, a, b, c, alpha) result(ans)
    real(kind=8), intent(in) :: m, g, a, b, c, alpha
    real(kind=8) :: ans, factor
    ! factor = d/dm [1 - c*(m - 1/m)] = -c*(1 + 1/m**2)
    factor = -c * ( 1.0d0 + 1.0d0/(m**2) )
    ans = 2.0d0 * b * m - a * alpha * ( 1.0d0 - c*( m - 1.0d0/m ) )**(alpha-1.0d0) * factor
  end function fx_prime


  !----------------------------------------------------
  !結果出力
  !----------------------------------------------------
  subroutine write_data(dir, tt, n, x, rho, u, p, time)
    character(len=*), intent(in) :: dir
    integer, intent(in) :: tt, n
    real(kind=8), intent(in) :: x(n), rho(n), u(n), p(n), time
    character(len=256) :: filename
    integer :: i, ios

    ! ファイル名
    write(filename, '(A,"/sod_",I4.4,".dat")') trim(dir), tt+1
    
    !ファイル操作
    open(unit=10, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
       print *, 'Error opening file: ', filename
       return
    end if
    ! ヘッダ行
    write(10, '(A, F10.5)') '! Time: ', time
    do i = 1, n
       !結果出力
       write(10, '(F10.5, 1X, F10.5, F10.5, F10.5)') x(i), u(i), rho(i), p(i)
    end do
    close(10)
  end subroutine write_data

end program shock_tube