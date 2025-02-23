program fvs_muscl
  implicit none
  !-------------------------------------
  ! パラメータ・定数
  !-------------------------------------
  integer, parameter :: nstep = 300
  integer, parameter :: nx0   = 100
  integer, parameter :: lvir  = 1
  integer, parameter :: nx    = nx0 + 2*lvir

  real, parameter :: dt     = 0.001
  real, parameter :: dx     = 0.01
  real, parameter :: gamma  = 1.4

  real, parameter :: k_muscl = 1.0/3.0
  real, parameter :: b_muscl = (3.0 - k_muscl) / (1.0 - k_muscl)

  character(len=*), parameter :: dir_name = "./../result"

  !-------------------------------------
  ! グローバル変数（添字は 0 ～ nx-1, 0～2 など）
  !-------------------------------------
  real, dimension(0:nx-1)         :: x
  real, dimension(0:nx-1,0:2)     :: qf      ! 基本量 [u, rho, p]
  real, dimension(0:nx-1,0:2)     :: Qc      ! 保存量 [rho, rho*u, e]
  real, dimension(0:nx-1,0:2)     :: Res

  real, dimension(0:2)           :: bol, bor
  real, dimension(0:nx,0:2)       :: Fplus, qfL, qfR, QcL, QcR

  integer :: k

  !-------------------------------------
  ! メイン処理   !OK
  !-------------------------------------
  call setup()
  call cre_dir()

  do k = 0, nstep-1
     print *, "Time step:", k
     call cal_Q()
     call Qctoqf(Qc, qf)
     call output_q()
  end do

contains

  !--------------------------------------------------------------
  ! 出力用フォルダの作成   !OK 
  !--------------------------------------------------------------
subroutine cre_dir()
  use globals
  implicit none
  integer :: ios
  character(len=100) :: command

  command = "mkdir " // trim(dir_name)
  call execute_command_line(command, wait=.true., exitstat=ios)
  if (ios /= 0) then
    print *, "Warning: Could not create directory ", trim(dir_name)
  end if
end subroutine cre_dir

  !----------------------------------------------------------------
  ! サブルーチン setup   !OK
  !----------------------------------------------------------------
  subroutine setup()
    implicit none
    integer :: i, j
    real :: pos
    real, dimension(0:nx-1) :: u_arr, rho_arr, p_arr, e_arr

    do i = 0, nx-1
       u_arr(i) = 0.0
       if ( real(i) <= 0.5*nx ) then
          rho_arr(i) = 1.0
          p_arr(i)   = 1.0
       else
          rho_arr(i) = 0.125
          p_arr(i)   = 0.1
       end if
       e_arr(i) = p_arr(i)/(gamma-1.0) + 0.5*rho_arr(i)*u_arr(i)**2
       x(i) = i*dx - dx/2.0
    end do

    ! 仮想境界セル（左端、右端）の設定
    bol(0) = rho_arr(0)
    bol(1) = u_arr(0)*rho_arr(0)
    bol(2) = e_arr(0)
    bor(0) = rho_arr(nx-1)
    bor(1) = u_arr(nx-1)*rho_arr(nx-1)
    bor(2) = e_arr(nx-1)

    ! 基本量 qf と保存量 Qc の初期化
    do i = 0, nx-1
       ! qf: [u, rho, p]
       qf(i,0) = u_arr(i)
       qf(i,1) = rho_arr(i)
       qf(i,2) = p_arr(i)
       ! Qc: [rho, rho*u, e]
       Qc(i,0) = rho_arr(i)
       Qc(i,1) = u_arr(i)*rho_arr(i)
       Qc(i,2) = e_arr(i)
    end do

  end subroutine setup

  !----------------------------------------------------------------
  ! サブルーチン cal_Q : 保存量 Qc の時間更新   !OK
  !----------------------------------------------------------------
  subroutine cal_Q()
    implicit none
    integer :: i, j
    real, dimension(0:nx-1,0:2) :: Qc_np

    call cal_Res()
    do i = 1, nx-2
       do j = 0, 2
          Qc(i,j) = Qc(i,j) - (dt/dx)*Res(i,j)
       end do
    end do
    !Qc_np = Qc
    !call cal_Res()
    !do i = 1, nx-2
       !do j = 0, 2
          !Qc(i,j) = 0.5*( Qc_np(i,j) + Qc(i,j) - (dt/dx)*Res(i,j) )
       !end do
    !end do
    call bound()
  end subroutine cal_Q

  !----------------------------------------------------------------
  ! サブルーチン bound : 境界条件の設定   !OK
  !----------------------------------------------------------------
  subroutine bound()
    implicit none
    integer :: j
    do j = 0, 2
       Qc(0,j)    = 2.0*bol(j) - Qc(1,j)
       Qc(nx-1,j) = Qc(nx-2,j)
    end do
  end subroutine bound

  !----------------------------------------------------------------
  ! サブルーチン cal_Res : 残差 Res の計算   !OK
  !----------------------------------------------------------------
  subroutine cal_Res()
    implicit none
    integer :: i, j

    do i = 0, nx-1
       do j = 0, 2
          Res(i,j) = 0.0
       end do
    end do

    call fvs()
    do i = 1, nx-2
       do j = 0, 2
          Res(i,j) = Fplus(i,j) - Fplus(i-1,j)
       end do
    end do
  end subroutine cal_Res

  !----------------------------------------------------------------
  ! サブルーチン fvs : FVS法による境界フラックスの計算
  !----------------------------------------------------------------
  subroutine fvs()
    implicit none
    integer :: i, j
    real, dimension(0:2,0:2) :: R, R_inv, Gam, Gam_abs, Ap, Am

    !call muscl()

    do i = 0, nx-2
       !call A_pm( QcL(i,0:2), qfL(i,0:2), R, R_inv, Gam, Gam_abs )
       call A_pm( i, R, R_inv, Gam, Gam_abs )
       Ap = matmul( matmul(R, (Gam+Gam_abs) ), R_inv )
       !call A_pm( QcR(i,0:2), qfR(i,0:2), R, R_inv, Gam, Gam_abs )
       call A_pm( i, R, R_inv, Gam, Gam_abs )
       Am = matmul( matmul(R, (Gam-Gam_abs) ), R_inv )
       Fplus(i,0:2) = 0.5 * ( matmul(Ap, QcL(i,0:2)) + matmul(Am, QcR(i,0:2)) )
       QcL(i, :) = Qc(i, :)
       QcR(i, :) = Qc(i+1, :)
    end do
  end subroutine fvs

  !--------------------------------------------------------------
  ! ヤコビアン対角化行列の計算 R, R_inv, Gam, Gam_abs
  !--------------------------------------------------------------
  subroutine A_pm(ite, R, R_inv, Gam, Gam_abs)
   !use globals
   implicit none
   integer, intent(in) :: ite
   real(4), intent(out) :: R(0:2,0:2), R_inv(0:2,0:2), Gam(0:2,0:2), Gam_abs(0:2,0:2)
   real(4) :: h, u, c, b_para, a_para
 
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

  !----------------------------------------------------------------
  ! サブルーチン A_pm : ヤコビアンの固有値分解関連量の計算
  !----------------------------------------------------------------
  !subroutine A_pm(lQc, lqf, R, R_inv, Gam, Gam_abs)
    !implicit none
    !real, intent(in)  :: lQc(0:2), lqf(0:2)
    !real, intent(out) :: R(0:2,0:2), R_inv(0:2,0:2), Gam(0:2,0:2), Gam_abs(0:2,0:2)
    !real :: H, u, c, b_para, a_para
    !integer :: i, j

    !H = (lQc(2) + lqf(2)) / lQc(0)
    !u = lqf(0)
    !c = sqrt((gamma-1.0)*(H - 0.5*u**2))
    !b_para = (gamma-1.0)/(c**2)
    !a_para = 0.5*b_para*u**2

    ! R行列の作成
    !R(0,0) = 1.0;   R(0,1) = 1.0;   R(0,2) = 1.0
    !R(1,0) = u-c;   R(1,1) = u;     R(1,2) = u+c
    !R(2,0) = H - u*c; R(2,1) = 0.5*u**2; R(2,2) = H + u*c

    ! R_inv の作成
    !R_inv(0,0) = 0.5*(a_para + u/c)
    !R_inv(0,1) = 0.5*(-b_para*u - 1.0/c)
    !R_inv(0,2) = 0.5*b_para
    !R_inv(1,0) = 1.0 - a_para
    !R_inv(1,1) = b_para*u
    !R_inv(1,2) = -b_para
    !R_inv(2,0) = 0.5*(a_para - u/c)
    !R_inv(2,1) = 0.5*(-b_para*u + 1.0/c)
    !R_inv(2,2) = 0.5*b_para

    ! Gam, Gam_abs（対角行列）の作成
    !do i = 0, 2
       !do j = 0, 2
          !am(i,j) = 0.0
          !Gam_abs(i,j) = 0.0
       !end do
    !end do
    !Gam(0,0) = u - c
    !Gam(1,1) = u
    !Gam(2,2) = u + c

    !Gam_abs(0,0) = abs(u - c)
    !Gam_abs(1,1) = abs(u)
    !Gam_abs(2,2) = abs(u + c)
  !end subroutine A_pm

  !----------------------------------------------------------------
  ! サブルーチン muscl : MUSCL法による再構成
  !----------------------------------------------------------------
  subroutine muscl()
    implicit none
    integer :: i, j
    real :: dplus, dminus, dplus_jp, dminus_jp

    call Qctoqf(Qc, qf)

    ! qfL, qfR を 0 初期化
    do i = 0, nx
       do j = 0, 2
          qfL(i,j) = 0.0
          qfR(i,j) = 0.0
       end do
    end do

    do i = 1, nx-3
       do j = 0, 2
          dplus      = qf(i+1,j) - qf(i,j)
          dminus     = qf(i,j)   - qf(i-1,j)
          dplus_jp   = qf(i+2,j) - qf(i+1,j)
          dminus_jp  = qf(i+1,j) - qf(i,j)
          !qfL(i,j) = qf(i,j) + 0.25 * ( (1.0 - k_muscl)*minmod(dminus, dplus, b_muscl) + &
                                        !(1.0 + k_muscl)*minmod(dplus, dminus, b_muscl) )
          !qfR(i,j) = qf(i+1,j) - 0.25 * ( (1.0 - k_muscl)*minmod(dplus_jp, dminus_jp, b_muscl) + &
                                          !(1.0 + k_muscl)*minmod(dminus_jp, dplus_jp, b_muscl) )
          qfL(i,j) = qf(i, j)
          qfR(i,j) = qf(i+1, j)                      
       end do
    end do

    ! 境界内側用
    do j = 0, 2
       dplus_jp  = qf(2,j) - qf(1,j)
       dminus_jp = qf(1,j) - qf(0,j)
       !qfR(0,j) = qf(1,j) - 0.25 * ( (1.0 - k_muscl)*minmod(dplus_jp, dminus_jp, b_muscl) + &
                                     !(1.0 + k_muscl)*minmod(dminus_jp, dplus_jp, b_muscl) )
       qfR(0,j) = qf(1, j) 
    end do

    do j = 0, 2
       dplus  = qf(nx-1,j) - qf(nx-2,j)
       dminus = qf(nx-2,j) - qf(nx-3,j)
       !qfL(nx-2,j) = qf(nx-2,j) + 0.25 * ( (1.0 - k_muscl)*minmod(dminus, dplus, b_muscl) + &
                                           !(1.0 + k_muscl)*minmod(dplus, dminus, b_muscl) )
       qfL(nx-2,j) = qf(nx-1, j) 
    end do

    ! qfL, qfR から保存量 QcL, QcR へ変換
    call qftoQc(qfL, QcL)
    call qftoQc(qfR, QcR)

    ! 境界外側（風上）の設定
    do j = 0, 2
       qfL(0,j)      = qf(0,j)
       QcL(0,j)      = Qc(0,j)
       qfR(nx-2,j)   = qf(nx-1,j)
       QcR(nx-2,j)   = Qc(nx-1,j)
    end do

  end subroutine muscl

  !----------------------------------------------------------------
  ! 関数 minmod : スロープリミッタ（minmod 型）
  !----------------------------------------------------------------
  real function minmod(x, y, b)
    implicit none
    real, intent(in) :: x, y, b
    if ( x == 0.0 ) then
       minmod = 0.0
    else
      minmod = sign(1.0, x)*max(0.0, min(abs(x), sign(1.0, x)*y*b))
    end if
  end function minmod

  !----------------------------------------------------------------
  ! サブルーチン qftoQc : 基本量 qf から保存量 Qc へ変換
  !----------------------------------------------------------------
  subroutine qftoQc(qf_in, Qc_out)
    implicit none
    real, intent(in)  :: qf_in(0:nx-1,0:2)
    real, intent(out) :: Qc_out(0:nx-1,0:2)
    integer :: i
    do i = 0, nx-1
       Qc_out(i,0) = qf_in(i,1)
       Qc_out(i,1) = qf_in(i,1) * qf_in(i,0)
       Qc_out(i,2) = qf_in(i,2)/(gamma-1.0) + 0.5*qf_in(i,1)*qf_in(i,0)**2
    end do
  end subroutine qftoQc

  !----------------------------------------------------------------
  ! サブルーチン Qctoqf : 保存量 Qc から基本量 qf へ変換
  !----------------------------------------------------------------
  subroutine Qctoqf(Qc_in, qf_out)
    implicit none
    real, intent(in)  :: Qc_in(0:nx-1,0:2)
    real, intent(out) :: qf_out(0:nx-1,0:2)
    integer :: i
    real :: vel
    do i = 0, nx-1
       if ( Qc_in(i,0) /= 0.0 ) then
          qf_out(i,0) = Qc_in(i,1) / Qc_in(i,0)
       else
          qf_out(i,0) = 0.0
       end if
       qf_out(i,1) = Qc_in(i,0)
       if ( Qc_in(i,0) /= 0.0 ) then
          vel = Qc_in(i,1)/Qc_in(i,0)
       else
          vel = 0.0
       end if
       qf_out(i,2) = (gamma-1.0)*(Qc_in(i,2) - 0.5*Qc_in(i,0)*vel**2)
    end do
  end subroutine Qctoqf

  !----------------------------------------------------------------
  ! サブルーチン output_q : 結果のテキスト出力
  !----------------------------------------------------------------
  subroutine output_q()
    implicit none
    integer :: i, unit
    character(len=100) :: filename, fullpath
    
    write(filename, '(A,I4.4)') "sod_", k
    fullpath = trim(dir_name) // "/" // trim(filename) // ".dat"
    unit = 10
    open(unit, file=fullpath, status='replace', action='write')
    write(unit,*) "x[m] u[m/s] rho[kg/m3] p[Pa]"
    do i = 0, nx-1
       write(unit, '(F8.4,1x,F8.4,1x,F8.4,1x,F8.4)') x(i), qf(i,0), qf(i,1), qf(i,2)
    end do
    close(unit)

  end subroutine output_q

end program fvs_muscl