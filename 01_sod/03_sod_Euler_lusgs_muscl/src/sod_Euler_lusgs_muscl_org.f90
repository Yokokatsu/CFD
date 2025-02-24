! sod_lusgs.f90
module sod_globals
  implicit none
  integer, parameter :: nstep = 300
  integer, parameter :: nx0   = 100
  integer, parameter :: lbound  = 1
  integer, parameter :: nx      = nx0 + 2*lbound

  real(8), parameter :: dt     = 0.002d0
  real(8), parameter :: dx     = 0.01d0
  real(8), parameter :: gamma  = 1.4d0

  real(8), parameter :: k_muscl = 1.0d0/3.0d0
  real(8), parameter :: b_muscl = (3.0d0 - k_muscl) / (1.0d0 - k_muscl)

  real(8), parameter :: norm_ok = 1.0d-6

  character(len=*), parameter :: dir_name      = "sod_lusgs_0.002s_py"
  character(len=*), parameter :: out_name_front = "time"
  character(len=*), parameter :: out_name_back  = "d-3"

  ! グローバル変数（すべて 0-indexed：Python のリストと同様の添字とする）
  real(8), dimension(0:nx-1)         :: x
  real(8), dimension(0:nx-1,0:2)     :: qf, Qc
  real(8), dimension(0:nx-1,0:2)     :: RHS
  real(8), dimension(0:2)           :: bol, bor

  ! LU-SGS 関連
  real(8), dimension(0:nx-1,0:2,0:2) :: Amatrix_plus, Amatrix_minus
  real(8), dimension(0:nx-1)         :: beta_sigma

  ! FVS 用：セル境界での変数
  real(8), dimension(0:nx,0:2)       :: Fplus
  real(8), dimension(0:nx,0:2)       :: qfL, qfR, QcL, QcR

contains
  !----------------------------------------------------------------
  ! setup: 初期状態（基本量，保存量，位置，仮想境界セル）の設定
  !----------------------------------------------------------------
  subroutine setup()
    implicit none
    integer :: i, j
    real(8) :: u(0:nx-1), rho(0:nx-1), p(0:nx-1), e(0:nx-1)

    do i = 0, nx-1
       u(i) = 0.0d0
       if ( real(i) <= 0.5d0*nx ) then
          rho(i) = 1.0d0
          p(i)   = 1.0d0
       else
          rho(i) = 0.125d0
          p(i)   = 0.1d0
       end if
       e(i) = p(i)/(gamma-1.0d0) + 0.5d0*rho(i)*u(i)**2
       x(i) = i*dx - dx/2.0d0
    end do

    do j = 0, 2
       if (j == 0) then
          bol(j) = rho(0)
          bor(j) = rho(nx-1)
       else if (j == 1) then
          bol(j) = u(0)*rho(0)
          bor(j) = u(nx-1)*rho(nx-1)
       else if (j == 2) then
          bol(j) = e(0)
          bor(j) = e(nx-1)
       end if
    end do

    do i = 0, nx-1
       ! 基本量： [u, rho, p]
       qf(i,0) = u(i)
       qf(i,1) = rho(i)
       qf(i,2) = p(i)
       ! 保存量： [rho, rho*u, e]
       Qc(i,0) = rho(i)
       Qc(i,1) = u(i)*rho(i)
       Qc(i,2) = e(i)
    end do
  end subroutine setup

  !----------------------------------------------------------------
  ! cal_Q: 時間発展（内部反復を用いた LU-SGS 法）
  !----------------------------------------------------------------
  subroutine cal_Q()
    implicit none
    call inner_ite(Qc)
  end subroutine cal_Q

  !----------------------------------------------------------------
  ! bound: 境界条件の設定
  !----------------------------------------------------------------
  subroutine bound(lQc)
    implicit none
    real(8), intent(inout) :: lQc(0:nx-1,0:2)
    integer :: j
    do j = 0, 2
       lQc(0,j)    = 2.0d0*bol(j) - lQc(1,j)
       lQc(nx-1,j) = lQc(nx-2,j)
    end do
  end subroutine bound

  !----------------------------------------------------------------
  ! cal_for_lusgs: LU-SGS 用の行列分解に必要な行列・係数の計算
  !----------------------------------------------------------------
  subroutine cal_for_lusgs(lQc, lqf)
    implicit none
    real(8), intent(in) :: lQc(0:nx-1,0:2), lqf(0:nx-1,0:2)
    integer :: i, j
    real(8) :: H, u, c, sigma, beta
    real(8) :: temp(0:2,0:2)
    real(8) :: R(0:2,0:2), R_inv(0:2,0:2), Gam(0:2,0:2), Gam_abs(0:2,0:2)
    real(8) :: A_matrix(0:2,0:2)
    real(8) :: I(0:2,0:2)
    beta = 1.1d0

    I = 0.0d0
    do i = 0, 2
       I(i,i) = 1.0d0
    end do

    do i = 0, nx-1
       H = (lQc(i,2) + lqf(i,2)) / lQc(i,0)
       u = lqf(i,0)
       c = sqrt((gamma-1.0d0)*(H - 0.5d0*u**2))
       sigma = abs(u) + c
       call A_pm(lQc(i,0:2), lqf(i,0:2), R, R_inv, Gam, Gam_abs)
       A_matrix = matmul(R, matmul(Gam, R_inv))
       temp = 0.0d0
       do j = 0, 2
         temp(j,j) = beta*sigma
       end do
       Amatrix_plus(i,0:2,0:2) = 0.5d0*(A_matrix + temp)
       Amatrix_minus(i,0:2,0:2) = 0.5d0*(A_matrix - temp)
       beta_sigma(i) = beta * sigma
    end do
  end subroutine cal_for_lusgs

  !----------------------------------------------------------------
  ! inner_ite: 内部反復（LU-SGS ループ）による Q の更新
  !----------------------------------------------------------------
  subroutine inner_ite(lQc)
    implicit none
    real(8), intent(inout) :: lQc(0:nx-1,0:2)
    integer :: ttt, i, j, ite, con
    real(8) :: s(0:2)
    real(8), dimension(0:nx-1,0:2) :: Qcn, Qcm, delta_Q, delta_Q_temp, lo_R, qfm
    real(8), dimension(0:2) :: sum_b, sum_b_Ax, norm2d
    real(8), dimension(0:nx-1,0:2) :: L, U
    real(8), dimension(0:nx-1) :: D

    Qcn = lQc
    Qcm = lQc

    do ttt = 1, 10
       ! 基本量に変換（保存量→基本量）
       call Qctoqf(Qcm, qfm)
       call cal_RHS(Qcm)
       call cal_for_lusgs(Qcm, qfm)
       delta_Q = 0.0d0
       delta_Q_temp = 0.0d0
       lo_R = RHS

       sum_b = 0.0d0
       do i = lbound, nx-lbound-1
          do j = 0, 2
             sum_b(j) = sum_b(j) + abs(lo_R(i,j))
          end do
       end do

       ite = 0
       con = 0
       do while (con == 0)
          delta_Q_temp = delta_Q

          ! LU-SGS 分解の各項：L, D, U の計算
          do i = 0, nx-1
             L(i,0) = dt*0.5d0*( Amatrix_plus(i,0,0)*delta_Q(i,0) + Amatrix_plus(i,0,1)*delta_Q(i,1) + Amatrix_plus(i,0,2)*delta_Q(i,2) )
             L(i,1) = dt*0.5d0*( Amatrix_plus(i,1,0)*delta_Q(i,0) + Amatrix_plus(i,1,1)*delta_Q(i,1) + Amatrix_plus(i,1,2)*delta_Q(i,2) )
             L(i,2) = dt*0.5d0*( Amatrix_plus(i,2,0)*delta_Q(i,0) + Amatrix_plus(i,2,1)*delta_Q(i,1) + Amatrix_plus(i,2,2)*delta_Q(i,2) )
             D(i) = 2.0d0*dx + dt*beta_sigma(i)
             U(i,0) = dt*0.5d0*( Amatrix_minus(i,0,0)*delta_Q(i,0) + Amatrix_minus(i,0,1)*delta_Q(i,1) + Amatrix_minus(i,0,2)*delta_Q(i,2) )
             U(i,1) = dt*0.5d0*( Amatrix_minus(i,1,0)*delta_Q(i,0) + Amatrix_minus(i,1,1)*delta_Q(i,1) + Amatrix_minus(i,1,2)*delta_Q(i,2) )
             U(i,2) = dt*0.5d0*( Amatrix_minus(i,2,0)*delta_Q(i,0) + Amatrix_minus(i,2,1)*delta_Q(i,1) + Amatrix_minus(i,2,2)*delta_Q(i,2) )
          end do

          do i = lbound, nx-lbound-1
             do j = 0, 2
                RHS(i,j) = -((Qcm(i,j) - Qcn(i,j))*dx) - dt*lo_R(i,j)
             end do
          end do

          do i = lbound, nx-lbound-1
             do j = 0, 2
                if (D(i) /= 0.0d0) then
                   delta_Q(i,j) = (L(i,j) + RHS(i,j)) / D(i)
                else
                   delta_Q(i,j) = 0.0d0
                end if
             end do
          end do

          do i = lbound, nx-lbound-1
             if (i+1 <= nx-1) then
                do j = 0, 2
                   delta_Q(i,j) = delta_Q(i,j) - U(i+1,j)/D(i)
                end do
             end if
          end do

          if (mod(ite+1,100) == 0) then
             sum_b_Ax = 0.0d0
             do i = lbound, nx-lbound-1
                do j = 0, 2
                   sum_b_Ax(j) = sum_b_Ax(j) + abs(delta_Q(i,j) - delta_Q_temp(i,j))
                end do
             end do
             norm2d = 0.0d0
             do j = 0, 2
                if (sum_b(j) /= 0.0d0) then
                   norm2d(j) = sum_b_Ax(j)/sum_b(j)
                else
                   norm2d(j) = 0.0d0
                end if
             end do
             if ( (norm2d(0) < norm_ok) .and. (norm2d(1) < norm_ok) .and. (norm2d(2) < norm_ok) ) then
                con = 1
             end if
          end if
          ite = ite + 1
       end do

       do i = lbound, nx-lbound-1
          do j = 0, 2
             Qcm(i,j) = Qcm(i,j) + delta_Q(i,j)
          end do
       end do

       call bound(Qcm)
       Qcm = Qcm  ! 更新後の Qcm をそのまま用いる
       lQc = Qcm
    end do
  end subroutine inner_ite

  !----------------------------------------------------------------
  ! inv_matrix: 3×3 行列の逆行列をガウスの消去法で求める
  !----------------------------------------------------------------
  subroutine inv_matrix(D)
    implicit none
    real(8), intent(inout) :: D(0:2,0:2)
    integer :: l, m, n
    real(8) :: aa, bb
    do l = 0, 2
       aa = 1.0d0 / D(l,l)
       D(l,l) = 1.0d0
       do m = 0, 2
          D(l,m) = D(l,m) * aa
       end do
       do m = 0, 2
          if (m /= l) then
             bb = D(m,l)
             D(m,l) = 0.0d0
             do n = 0, 2
                D(m,n) = D(m,n) - bb*D(l,n)
             end do
          end if
       end do
    end do
  end subroutine inv_matrix

  !----------------------------------------------------------------
  ! cal_RHS: 境界フラックスから RHS を計算
  !----------------------------------------------------------------
  subroutine cal_RHS(lQc)
    implicit none
    real(8), intent(in) :: lQc(0:nx-1,0:2)
    integer :: i, j
    do i = 0, nx-1
       do j = 0, 2
          RHS(i,j) = 0.0d0
       end do
    end do
    call fvs(lQc)
    do i = 1, nx-2
       do j = 0, 2
          RHS(i,j) = Fplus(i,j) - Fplus(i-1,j)
       end do
    end do
  end subroutine cal_RHS

  !----------------------------------------------------------------
  ! fvs: FVS 法による境界フラックスの計算
  !----------------------------------------------------------------
  subroutine fvs(lQc)
    implicit none
    real(8), intent(in) :: lQc(0:nx-1,0:2)
    integer :: i
    real(8) :: R(0:2,0:2), R_inv(0:2,0:2), Gam(0:2,0:2), Gam_abs(0:2,0:2)
    real(8) :: Ap(0:2,0:2), Am(0:2,0:2)
    real(8) :: tmp1(0:2), tmp2(0:2)
    Fplus = 0.0d0
    call muscl(lQc)
    do i = 0, nx-2
       call A_pm(QcL(i,0:2), qfL(i,0:2), R, R_inv, Gam, Gam_abs)
       Ap = matmul(R, matmul(Gam+Gam_abs, R_inv))
       call A_pm(QcR(i,0:2), qfR(i,0:2), R, R_inv, Gam, Gam_abs)
       Am = matmul(R, matmul(Gam-Gam_abs, R_inv))
       tmp1 = matmul(Ap, QcL(i,0:2))
       tmp2 = matmul(Am, QcR(i,0:2))
       Fplus(i,0:2) = 0.5d0*(tmp1 + tmp2)
    end do
  end subroutine fvs

  !----------------------------------------------------------------
  ! A_pm: ヤコビアンの固有値等の計算
  !----------------------------------------------------------------
  subroutine A_pm(lQc, lqf, R, R_inv, Gam, Gam_abs)
    implicit none
    real(8), intent(in) :: lQc(0:2), lqf(0:2)
    real(8), intent(out) :: R(0:2,0:2), R_inv(0:2,0:2), Gam(0:2,0:2), Gam_abs(0:2,0:2)
    real(8) :: H, u, c, b_para, a_para
    integer :: i, j

    H = (lQc(2) + lqf(2)) / lQc(0)
    u = lqf(0)
    c = sqrt((gamma-1.0d0)*(H - 0.5d0*u**2))
    b_para = (gamma-1.0d0)/(c**2)
    a_para = 0.5d0*b_para*u**2

    R(0,0) = 1.0d0;   R(0,1) = 1.0d0;   R(0,2) = 1.0d0
    R(1,0) = u - c;   R(1,1) = u;     R(1,2) = u + c
    R(2,0) = H - u*c; R(2,1) = 0.5d0*u**2; R(2,2) = H + u*c

    R_inv(0,0) = 0.5d0*(a_para + u/c)
    R_inv(0,1) = 0.5d0*(-b_para*u - 1.0d0/c)
    R_inv(0,2) = 0.5d0*b_para
    R_inv(1,0) = 1.0d0 - a_para
    R_inv(1,1) = b_para*u
    R_inv(1,2) = -b_para
    R_inv(2,0) = 0.5d0*(a_para - u/c)
    R_inv(2,1) = 0.5d0*(-b_para*u + 1.0d0/c)
    R_inv(2,2) = 0.5d0*b_para

    Gam = 0.0d0
    Gam_abs = 0.0d0
    Gam(0,0) = u - c
    Gam(1,1) = u
    Gam(2,2) = u + c

    Gam_abs(0,0) = abs(u - c)
    Gam_abs(1,1) = abs(u)
    Gam_abs(2,2) = abs(u + c)
  end subroutine A_pm

  !----------------------------------------------------------------
  ! muscl: MUSCL 法による再構成（基本量→保存量の境界補間）
  !----------------------------------------------------------------
  subroutine muscl(lQc)
    implicit none
    real(8), intent(in) :: lQc(0:nx-1,0:2)
    integer :: i, j
    real(8), dimension(0:nx-1,0:2) :: lqf
    call Qctoqf(lQc, lqf)
    qfL = 0.0d0
    qfR = 0.0d0

    do i = 1, nx-3
       do j = 0, 2
          real(8) :: dplus_j, dminus_j, dplus_jp, dminus_jp
          dplus_j  = lqf(i+1,j) - lqf(i,j)
          dminus_j = lqf(i,j)   - lqf(i-1,j)
          dplus_jp = lqf(i+2,j) - lqf(i+1,j)
          dminus_jp= lqf(i+1,j) - lqf(i,j)
          qfL(i,j) = lqf(i,j) + 0.25d0*((1.0d0 - k_muscl)*minmod(dminus_j, dplus_j, b_muscl) + &
                                        (1.0d0 + k_muscl)*minmod(dplus_j, dminus_j, b_muscl))
          qfR(i,j) = lqf(i+1,j) - 0.25d0*((1.0d0 - k_muscl)*minmod(dplus_jp, dminus_jp, b_muscl) + &
                                        (1.0d0 + k_muscl)*minmod(dminus_jp, dplus_jp, b_muscl))
       end do
    end do

    do j = 0, 2
       real(8) :: dplus_jp, dminus_j
       dplus_jp = lqf(2,j) - lqf(1,j)
       dminus_j = lqf(1,j) - lqf(0,j)
       qfR(0,j) = lqf(1,j) - 0.25d0*((1.0d0 - k_muscl)*minmod(dplus_jp, dminus_j, b_muscl) + &
                                     (1.0d0 + k_muscl)*minmod(dminus_j, dplus_jp, b_muscl))
       dplus_jp = lqf(nx-1,j) - lqf(nx-2,j)
       dminus_j = lqf(nx-2,j) - lqf(nx-3,j)
       qfL(nx-2,j) = lqf(nx-2,j) + 0.25d0*((1.0d0 - k_muscl)*minmod(dminus_j, dplus_jp, b_muscl) + &
                                          (1.0d0 + k_muscl)*minmod(dplus_jp, dminus_j, b_muscl))
    end do

    call qftoQc(qfL, QcL)
    call qftoQc(qfR, QcR)

    qfL(0,0:2) = lqf(0,0:2)
    QcL(0,0:2) = lQc(0,0:2)
    qfR(nx-2,0:2) = lqf(nx-1,0:2)
    QcR(nx-2,0:2) = lQc(nx-1,0:2)
  end subroutine muscl

  !----------------------------------------------------------------
  ! minmod: スロープリミッタ関数
  !----------------------------------------------------------------
  real(8) function minmod(x, y, b)
    implicit none
    real(8), intent(in) :: x, y, b
    real(8) :: s
    s = sign(1.0d0, x)
    minmod = s * max(0.0d0, min(abs(x), s*y*b))
  end function minmod

  !----------------------------------------------------------------
  ! qftoQc: 基本量（qf）から保存量（Qc）への変換
  !----------------------------------------------------------------
  subroutine qftoQc(qf_in, Qc_out)
    implicit none
    real(8), intent(in) :: qf_in(0:nx-1,0:2)
    real(8), intent(out) :: Qc_out(0:nx-1,0:2)
    integer :: i, j
    do i = 0, nx-1
       do j = 0, 2
          if (j == 0) then
             Qc_out(i,j) = qf_in(i,1)
          else if (j == 1) then
             Qc_out(i,j) = qf_in(i,1) * qf_in(i,0)
          else if (j == 2) then
             Qc_out(i,j) = qf_in(i,2)/(gamma-1.0d0) + 0.5d0*qf_in(i,1)*qf_in(i,0)**2
          end if
       end do
    end do
  end subroutine qftoQc

  !----------------------------------------------------------------
  ! Qctoqf: 保存量（Qc）から基本量（qf）への変換
  !----------------------------------------------------------------
  subroutine Qctoqf(Qc_in, qf_out)
    implicit none
    real(8), intent(in) :: Qc_in(0:nx-1,0:2)
    real(8), intent(out) :: qf_out(0:nx-1,0:2)
    integer :: i, j
    real(8) :: vel
    do i = 0, nx-1
       do j = 0, 2
          if (j == 0) then
             qf_out(i,j) = Qc_in(i,1) / Qc_in(i,0)
          else if (j == 1) then
             qf_out(i,j) = Qc_in(i,0)
          else if (j == 2) then
             vel = Qc_in(i,1)/Qc_in(i,0)
             qf_out(i,j) = (gamma-1.0d0)*(Qc_in(i,2) - 0.5d0*Qc_in(i,0)*vel**2)
          end if
       end do
    end do
  end subroutine Qctoqf

  !----------------------------------------------------------------
  ! output_q: 結果のテキスト出力
  !----------------------------------------------------------------
  subroutine output_q(f_name)
    implicit none
    character(len=*), intent(in) :: f_name
    integer :: i, unit
    character(len=200) :: full_filename
    full_filename = trim(dir_name)//"/"//trim(f_name)
    unit = 10
    open(unit, file=full_filename, status='replace', action='write')
    write(unit,*) "x[m] u[m/s] rho[kg/m3] p[Pa]"
    do i = 0, nx-1
       write(unit, '(F8.4,1x,F8.4,1x,F8.4,1x,F8.4)') x(i), qf(i,0), qf(i,1), qf(i,2)
    end do
    close(unit)
  end subroutine output_q

  !----------------------------------------------------------------
  ! cre_dir: 出力ディレクトリの作成（execute_command_line を利用）
  !----------------------------------------------------------------
  subroutine cre_dir()
    implicit none
    call execute_command_line("mkdir -p " // trim(dir_name))
  end subroutine cre_dir

end module sod_globals

program sod_lusgs
  use sod_globals
  implicit none
  integer :: k
  character(len=30) :: filename
  call cre_dir()
  call setup()

  do k = 0, nstep-1
     print *, "Time step:", k
     call cal_Q()
     call Qctoqf(Qc, qf)
     write(filename, '(A,I0,A)') trim(out_name_front), int(k*dt*1000), trim(out_name_back)
     call output_q(filename)
  end do

end program sod_lusgs
