!==============================================================
! shock tube Euler explicit 
! 2025/02/10 Yokoyama Katsuyuki
!==============================================================

module globals
    implicit none
    ! parameter
    integer, parameter :: nstep = 300
    integer, parameter :: interval = 1
    integer, parameter :: nx0   = 100
    integer, parameter :: lpml  = 1
    integer, parameter :: nx    = nx0 + 2*lpml
    real(8), parameter :: dt    = 0.001d0
    real(8), parameter :: dx    = 0.01d0
    real(8), parameter :: gamma = 1.4d0
    real(8), parameter :: k_muscl = 1.0/3.0
    real(8), parameter :: b_muscl = (3.0 - k_muscl) / (1.0 - k_muscl)
    character(len=*), parameter :: dir_name = "./../result"
  
    ! global variable
    real(8), allocatable :: x(:)         ! 座標
    real(8), allocatable :: bol(:)       ! 左側境界
    real(8), allocatable :: bor(:)       ! 右側境界
    real(8), allocatable :: qf(:,:)      ! 基本量 [u, ρ, p]
    real(8), allocatable :: Qc(:,:)      ! 保存量 [ρ, ρu, E]
    real(8), allocatable :: Res(:,:)     ! 対流項
    real(8), allocatable :: Fplus(:,:)   ! 流束

    real(8), allocatable :: Qc_np(:,:)      ! 保存量 [ρ, ρu, E]
    !real(8), allocatable :: QcL(:,:), QcR(:,:), qfL(:,:), qfR(:,:)
    real(8) :: QcL(0:nx-1,0:2), QcR(0:nx-1,0:2), qfL(0:nx-1,0:2), qfR(0:nx-1,0:2)
  end module globals
  !--------------------------------------------------------------
  ! main loop
  !--------------------------------------------------------------
  program main
    use globals
    implicit none
    integer :: k
  
    !call cre_dir()   ! 出力フォルダの作成
    call setup()      ! 初期値の設定
  
    do k = 1, nstep
      write(*,*) "Time step:", k
      call cal_Q()             ! 保存量の更新
      call Qctoqf(Qc, qf)            ! 保存量から基本量へ変換
      call output_q(k)         ! 結果の出力
    end do
  
  end program main

  !--------------------------------------------------------------
  ! 出力用フォルダの作成 
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

  !--------------------------------------------------------------
  ! 初期値の設定
  !--------------------------------------------------------------
  subroutine setup()
    use globals
    implicit none
    integer :: i, j
    real(8), allocatable :: u(:), rho(:), p(:), e(:)
  
    allocate(x(0:nx-1),u(0:nx-1), rho(0:nx-1), p(0:nx-1), e(0:nx-1))
  
    do i = 0, nx-1
      u(i) = 0.0d0
      if ( real(i) <= 0.5d0 * real(nx) ) then
        rho(i) = 1.0d0
        p(i)   = 1.0d0
      else
        rho(i) = 0.125d0
        p(i)   = 0.1d0
      end if
      e(i) = p(i)/(gamma-1.0d0) + 0.5d0 * rho(i)*u(i)**2
      x(i) = i*dx - dx/2.0d0
    end do
  
    allocate(bol(0:2), bor(0:2))
    do j = 0, 2
      select case(j)
      case (0)
        bol(j) = rho(0)
        bor(j) = rho(nx-1)
      case (1)
        bol(j) = u(0)*rho(0)
        bor(j) = u(nx-1)*rho(nx-1)
      case (2)
        bol(j) = e(0)
        bor(j) = e(nx-1)
      end select
    end do
  
    allocate(qf(0:nx-1,0:2), Qc(0:nx-1,0:2))
    do i = 0, nx-1
      ! qf: [u, ρ, p]
      qf(i,0) = u(i)
      qf(i,1) = rho(i)
      qf(i,2) = p(i)
      ! Qc: [ρ, ρu, E]
      Qc(i,0) = rho(i)
      Qc(i,1) = u(i)*rho(i)
      Qc(i,2) = e(i)
    end do
  
    deallocate(u, rho, p, e)
  end subroutine setup

  !--------------------------------------------------------------
  ! 時間進行
  !--------------------------------------------------------------
  subroutine cal_Q()
    use globals
    implicit none
    integer :: i, j
  
    call cal_Res()
  
    do i = 1, nx-2
      do j = 0, 2
        Qc(i,j) = Qc(i,j) - dt/dx * Res(i,j)
      end do
    end do
    Qc_np = Qc
    do i = 1, nx-2
      do j = 0, 2
        Qc(i,j) = 0.5*(Qc_np(i,j)+Qc(i,j) - dt/dx * Res(i,j))
      end do
    end do
  
    call bound()
  end subroutine cal_Q

  !--------------------------------------------------------------
  ! 対流項の評価（風上差分）
  !--------------------------------------------------------------
  subroutine cal_Res()
    use globals
    implicit none
    integer :: i, j
  
    if (.not. allocated(Res)) then
      allocate(Res(0:nx-1,0:2))
    end if
    Res = 0.0d0
  
    call fvs()
  
    do i = 1, nx-2
      do j = 0, 2
        Res(i,j) = Fplus(i,j) - Fplus(i-1,j)
      end do
    end do
  end subroutine cal_Res


  !--------------------------------------------------------------
  ! 境界の流束を評価（FVS）
  !--------------------------------------------------------------
  subroutine fvs()
    use globals
    implicit none
    integer :: i, j
    real(8) :: Ap(0:2,0:2), Am(0:2,0:2)
    real(8) :: R(0:2,0:2), R_inv(0:2,0:2), Gam(0:2,0:2), Gam_abs(0:2,0:2)
    !real(8) :: QcL(0:2), QcR(0:2)
  
    if (.not. allocated(Fplus)) then
      allocate(Fplus(0:nx,0:2))
    end if
    if (.not. allocated(Fplus)) then
      allocate(Fplus(0:nx,0:2))
    end if
    Fplus = 0.0d0

    call muscl()
  
    do i = 0, nx-2
      ! iセルにおけるR, R^-1,Λ,|Λ|
      call A_pm(i, R, R_inv, Gam, Gam_abs)
      Ap = matmul(R, matmul(Gam + Gam_abs, R_inv))
      ! i+1セルにおけるR, R^-1,Λ,|Λ|
      call A_pm(i+1, R, R_inv, Gam, Gam_abs)
      Am = matmul(R, matmul(Gam - Gam_abs, R_inv))
      ! 境界フラックスを計算 A(+)*Q(j) + A(-)*Q(j+1)
      !QcL(i, :) = Qc(i, :)
      !QcR(i, :) = Qc(i+1, :)
      Fplus(i, :) = 0.5d0 * ( matmul(Ap, QcL(i,:)) + matmul(Am, QcR(i,:)) ) 
    end do
  end subroutine fvs

  !----------------------------------------------------------------
  ! サブルーチン muscl : MUSCL法による流束の評価
  !----------------------------------------------------------------
  subroutine muscl()
    use globals
    implicit none
    integer :: i, j
    real :: dplus, dminus, dplus_jp, dminus_jp, minmod

    call Qctoqf(Qc, qf)

    ! n0-MUSCL
    !do i = 0, nx
      !do j = 0, 2
         !QcL(i, :) = Qc(i, :)
         !QcR(i, :) = Qc(i+1, :)
      !end do
    !end do

    ! qfL, qfR を 0 初期化
    do i = 0, nx-1
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
          qfL(i,j) = qf(i,j) + 0.25 * ( (1.0 - k_muscl)*minmod(dminus, dplus, b_muscl) &
                                       +(1.0 + k_muscl)*minmod(dplus, dminus, b_muscl) )
          qfR(i,j) = qf(i+1,j) - 0.25 * ( (1.0 - k_muscl)*minmod(dplus_jp, dminus_jp, b_muscl) &
                                         +(1.0 + k_muscl)*minmod(dminus_jp, dplus_jp, b_muscl) )
       end do
    end do

    ! 境界内側用
    do j = 0, 2
       dplus_jp  = qf(2,j) - qf(1,j)
       dminus_jp = qf(1,j) - qf(0,j)
       qfR(0,j) = qf(1,j) - 0.25 * ( (1.0 - k_muscl)*minmod(dplus_jp, dminus_jp, b_muscl) &
                                    +(1.0 + k_muscl)*minmod(dminus_jp, dplus_jp, b_muscl) )
    end do

    do j = 0, 2
       dplus  = qf(nx-1,j) - qf(nx-2,j)
       dminus = qf(nx-2,j) - qf(nx-3,j)
       qfL(nx-2,j) = qf(nx-2,j) + 0.25 * ( (1.0 - k_muscl)*minmod(dminus, dplus, b_muscl) &
                                          +(1.0 + k_muscl)*minmod(dplus, dminus, b_muscl) )
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

  !--------------------------------------------------------------
  ! 境界条件の設定
  !--------------------------------------------------------------
  subroutine bound()
    use globals
    implicit none
    integer :: j
  
    ! 左端： Qc(0,:) = 2*bol - Qc(1,:)
    do j = 0, 2
      Qc(0,j) = 2.0d0 * bol(j) - Qc(1,j)
    end do
  
    ! 右端： Qc(nx-1,:) = Qc(nx-2,:)
    do j = 0, 2
      Qc(nx-1,j) = Qc(nx-2,j)
    end do
  end subroutine bound


  !----------------------------------------------------------------
  ! 基本量qfから保存量Qcへの変換
  !----------------------------------------------------------------
  subroutine qftoQc(qf_in, Qc_out)
    use globals
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
  ! 保存量Qcから基本量qfへの変換
  !----------------------------------------------------------------
  subroutine Qctoqf(Qc_in, qf_out)
    use globals
    implicit none
    real, intent(in)  :: Qc_in(0:nx-1,0:2)
    real, intent(out) :: qf_out(0:nx-1,0:2)
    integer :: i
    real :: vel
    do i = 1, nx- 1
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

  !--------------------------------------------------------------
  ! 結果の出力
  !--------------------------------------------------------------
  subroutine output_q(k)
    use globals
    implicit none
    integer, intent(in) :: k
    character(len=100) :: filename, fullpath
    integer :: i, ios
    
    if (mod(k,interval)==0)then
      write(filename, '(A,I4.4)') "sod_", k
      fullpath = trim(dir_name) // "/" // trim(filename) // ".dat"
  
      open(unit=10, file=fullpath, status='replace', action='write', iostat=ios)
      if (ios /= 0) then
        print *, "Error opening file: ", fullpath
        return
      end if
  
      write(10,*) "x[m] u[m/s] rho[kg/m3] p[Pa]"
      do i = 0, nx-1
        write(10, '(F10.5,1X,F10.5,1X,F10.5,1X,F10.5)') x(i), qf(i,0), qf(i,1), qf(i,2)
      end do
      close(10)
    end if
  end subroutine output_q