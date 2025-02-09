!==============================================================
! Fortran コード (fvs.f90)
! Python コードの内容を忠実に再現した FVS 法による数値計算プログラム
!==============================================================

module globals
    implicit none
    ! パラメータ定義
    integer, parameter :: nstep = 300
    integer, parameter :: nx0   = 100
    integer, parameter :: lpml  = 1
    integer, parameter :: nx    = nx0 + 2*lpml
    real(8), parameter :: dt    = 0.001d0
    real(8), parameter :: dx    = 0.01d0
    real(8), parameter :: gamma = 1.4d0
    character(len=*), parameter :: dir_name = "./../result" ww
  
    ! グローバル変数（すべて 0-indexed とする）
    real(8), allocatable :: x(:)         ! 位置, サイズ: 0:nx-1
    real(8), allocatable :: bol(:)       ! 左端仮想セル（3成分）: 0:2
    real(8), allocatable :: bor(:)       ! 右端仮想セル（3成分）: 0:2
    real(8), allocatable :: qf(:,:)      ! 基本量 [u, ρ, p], サイズ: (0:nx-1,0:2)
    real(8), allocatable :: Qc(:,:)      ! 保存量 [ρ, ρu, E], サイズ: (0:nx-1,0:2)
    real(8), allocatable :: Res(:,:)     ! 残差（Flux 差分）, サイズ: (0:nx-1,0:2)
    real(8), allocatable :: Fplus(:,:)   ! 境界フラックス, サイズ: (0:nx,0:2)
    
  end module globals
  !--------------------------------------------------------------
  ! メインプログラム
  !--------------------------------------------------------------
  program main
    use globals
    implicit none
    integer :: k
  
    !call cre_dir()   ! 出力フォルダの作成
    call setup()     ! 初期値の設定
  
    do k = 1, nstep
      write(*,*) "Time step:", k
      call cal_Q()             ! 保存量の更新
      call Qctoqf()      ! 保存量から基本量へ変換
      call output_q(k)  ! 結果の出力
    end do
  
  end program main

  !--------------------------------------------------------------
  ! サブルーチン cre_dir
  ! 出力用フォルダの作成 (Unix 系 OS では "mkdir -p" コマンドを利用)
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
  ! サブルーチン setup
  ! 初期値の設定：セルごとに速度 u, 密度 ρ, 圧力 p, エネルギー E, 位置 x の初期化を行う
  ! また、仮想境界セル bol, bor と基本量 qf、保存量 Qc を設定する
  !--------------------------------------------------------------
  subroutine setup()
    use globals
    implicit none
    integer :: i, j
    real(8), allocatable :: u(:), rho(:), p(:), e(:)
  
    allocate(u(0:nx-1), rho(0:nx-1), p(0:nx-1), e(0:nx-1))
    allocate(x(0:nx-1))
  
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
  
    ! 基本量 qf と保存量 Qc の初期設定
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
  ! サブルーチン cal_Q
  ! 保存量 Qc を計算する。cal_Res により残差 Res を計算し、中央差分で更新後、境界条件を課す
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
  
    call bound()
  end subroutine cal_Q

  !--------------------------------------------------------------
  ! サブルーチン cal_Res
  ! 境界フラックス Fplus を計算し、各セル i について Res(i,:) = Fplus(i,:) - Fplus(i-1,:) とする
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
  ! サブルーチン fvs
  ! FVS 法により、各セル境界のフラックス Fplus を計算する
  !--------------------------------------------------------------
  subroutine fvs()
    use globals
    implicit none
    integer :: i, j
    real(8) :: Ap(0:2,0:2), Am(0:2,0:2)
    real(8) :: R(0:2,0:2), R_inv(0:2,0:2), Gam(0:2,0:2), Gam_abs(0:2,0:2)
    real(8) :: vec1(0:2), vec2(0:2)
  
    if (.not. allocated(Fplus)) then
      allocate(Fplus(0:nx,0:2))
    end if
    Fplus = 0.0d0
  
    do i = 0, nx-2
      ! i セル側の固有情報を計算
      call A_pm(i, R, R_inv, Gam, Gam_abs)
      Ap = matmul(R, matmul(Gam + Gam_abs, R_inv))
      ! i+1 セル側の固有情報を計算
      call A_pm(i+1, R, R_inv, Gam, Gam_abs)
      Am = matmul(R, matmul(Gam - Gam_abs, R_inv))
      ! それぞれのセルにおける保存量ベクトル Qc(i,:) と Qc(i+1,:) を用いてフラックスを計算
      vec1 = Qc(i, :)
      vec2 = Qc(i+1, :)
      Fplus(i, :) = 0.5d0 * ( matmul(Ap, vec1) + matmul(Am, vec2) )
    end do
  end subroutine fvs

  !--------------------------------------------------------------
  ! サブルーチン A_pm
  ! 指定セル ite において、ヤコビアンの固有値分解に必要な行列 R, R_inv, Gam, Gam_abs を計算する
  !--------------------------------------------------------------
  subroutine A_pm(ite, R, R_inv, Gam, Gam_abs)
    use globals
    implicit none
    integer, intent(in) :: ite
    real(8), intent(out) :: R(0:2,0:2), R_inv(0:2,0:2), Gam(0:2,0:2), Gam_abs(0:2,0:2)
    real(8) :: H, u, c, b_para, a_para
  
    ! エンタルピー H の計算： H = (E + p) / ρ  ※ E = Qc(:,2), p = qf(:,2)
    H = (Qc(ite,2) + qf(ite,2)) / Qc(ite,0)
    u = qf(ite,0)
    c = sqrt((gamma - 1.0d0) * (H - 0.5d0*u*u))
    b_para = (gamma - 1.0d0) / (c*c)
    a_para = 0.5d0 * b_para * u*u
  
    ! R 行列（各固有ベクトルを列に持つ）
    R(0,0) = 1.0d0;    R(0,1) = 1.0d0;     R(0,2) = 1.0d0
    R(1,0) = u - c;    R(1,1) = u;         R(1,2) = u + c
    R(2,0) = H - u*c;  R(2,1) = 0.5d0*u*u;   R(2,2) = H + u*c
  
    ! 逆行列 R_inv の計算
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
  ! サブルーチン bound
  ! 左端と右端の境界条件を設定する
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

  !--------------------------------------------------------------
  ! サブルーチン qftoQc
  ! 基本量 qf から保存量 Qc への変換 (本コード内では直接は使用していません)
  !--------------------------------------------------------------
  subroutine qftoQc()
    use globals
    implicit none
    !real(8), intent(in) :: qf(0:nx-1,0:2)
    !real(8), intent(out) :: Qc(0:nx-1,0:2)
    integer :: i
  
    do i = 0, nx-1
      Qc(i,0) = qf(i,1)
      Qc(i,1) = qf(i,1) * qf(i,0)
      Qc(i,2) = qf(i,2)/(gamma-1.0d0) + 0.5d0*qf(i,1)*qf(i,0)**2
    end do
  end subroutine qftoQc

  !--------------------------------------------------------------
! サブルーチン Qctoqf
! 保存量 Qc から基本量 qf への変換 (本メインループで用いる)
!--------------------------------------------------------------
subroutine Qctoqf()
    use globals
    implicit none
    !real(8), intent(inout) :: qf(0:nx-1,0:2)
    !real(8), intent(in) :: Qc(0:nx-1,0:2)
    integer :: i
    real(8) :: vel
  
    do i = 0, nx-1
      ! 速度 u = ρu / ρ
      qf(i,0) = Qc(i,1) / Qc(i,0)
      ! 密度 ρ
      qf(i,1) = Qc(i,0)
      ! 圧力 p = (γ-1) [E - 1/2 ρ u^2]
      vel = Qc(i,1) / Qc(i,0)
      qf(i,2) = (gamma - 1.0d0) * (Qc(i,2) - 0.5d0 * Qc(i,0)*vel*vel)
    end do
  end subroutine Qctoqf

  !--------------------------------------------------------------
  ! サブルーチン output_q
  ! 基本量 qf と位置 x をテキスト形式で出力する
  ! ファイル名は "timeXXXXd-3" (XXXX は 4 桁のゼロ埋め整数) とする
  !--------------------------------------------------------------
  subroutine output_q(k)
    use globals
    implicit none
    integer, intent(in) :: k
    character(len=100) :: filename, fullpath
    integer :: i, ios
    

    !if (mod(k,10)==0)then
      ! ファイル名の作成: "timeXXXXd-3"
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
    !end if
  end subroutine output_q