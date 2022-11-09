!-------------------------------------------------------------------------------------------------------------------------------------------
! Tests isotopologue_coordinates_conversion_mod module.
!-------------------------------------------------------------------------------------------------------------------------------------------
program conversion_test
  use isotopologue_coordinates_conversion_mod
  use iso_fortran_env, only: real64
  implicit none
  
  real(real64), parameter :: pi = acos(-1d0)
  real(real64), parameter :: ang_per_bohr = 0.52917721092d0
  real(real64), parameter :: deg_per_rad = 180 / pi
  real(real64), parameter :: oxygen_masses(3) = [15.99491461956d0, 16.99913170d0, 17.9991596129d0]
  real(real64), parameter :: argon_mass = 39.9623831225d0
  real(real64), allocatable :: grid_1(:), grid_2(:), grid_3(:)
  real(real64), allocatable :: iso_coords(:, :), ref_coords(:, :)
  
  grid_1 = generate_real_range(1d0, 5d0, 0.1d0)
  grid_2 = generate_real_range(0d0, 180d0, 2d0) / 180d0 * pi
  grid_3 = generate_real_range(0d0, 359d0, 4d0) / 180d0 * pi
  iso_coords = grids_to_coord_list(grid_1, grid_2, grid_3)
  ref_coords = convert_to_reference_euler(iso_coords)

  print *, ref_coords

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates real grid using specified step.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function generate_real_range(start, end, step) result(grid)
    real(real64), intent(in) :: start, end, step
    real(real64), allocatable :: grid(:)
    real(real64) :: npoints_real
    integer :: npoints, i

    npoints_real = (end - start) / step
    npoints = nint(npoints_real)
    if (npoints < npoints_real .or. npoints - npoints_real > 1d-10) then
      npoints = int(npoints_real)
    end if

    ! +1 because of an extra point at the beginning of interval
    allocate(grid(npoints + 1))
    grid = [(start + i * step, i = 0, npoints)]
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns direct product of the three grids. Each point is a column.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function grids_to_coord_list(grid_1, grid_2, grid_3) result(coords)
    real(real64), intent(in) :: grid_1(:), grid_2(:), grid_3(:)
    real(real64), allocatable :: coords(:, :)
    integer :: total_points, ind_1, ind_2, ind_3, i

    total_points = size(grid_1) * size(grid_2) * size(grid_3)
    allocate(coords(3, total_points))
    i = 1
    do ind_1 = 1, size(grid_1)
      do ind_2 = 1, size(grid_2)
        do ind_3 = 1, size(grid_3)
          coords(:, i) = [grid_1(ind_1), grid_2(ind_2), grid_3(ind_3)]
          i = i + 1
        end do
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Uses euler module to сonvert coordinates from a given isotopomer to reference.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function convert_to_reference_euler(coords_iso) result(coords_ref)
    real(real64), intent(in) :: coords_iso(:, :) ! 3 x N
    real(real64) :: coords_ref(size(coords_iso, 1), size(coords_iso, 2))
    integer :: zx_atom_index_1(2), zx_atom_index_2(2), i
    real(real64), allocatable :: cartesian_input_1(:, :), cartesian_input_2(:, :), cartesian_iso_1(:, :), cartesian_iso_2(:, :)
    real(real64), allocatable :: masses_iso_1(:), masses_iso_2(:), masses_ref_1(:), masses_ref_2(:)
    real(real64) :: coords_iso_full(6), coords_ref_full(6)

    cartesian_input_1 = reshape([0d0, 0d0, 0d0, 1.2717d0, 0d0, 0d0, 1.2717d0 * (1 + cos((180 - 116.84) / deg_per_rad)), 0d0, 1.2717d0 * sin((180 - 116.84) / deg_per_rad)], [3, 3]) / ang_per_bohr
    cartesian_input_2 = reshape([0d0, 0d0, 0d0], [3, 1])
    masses_iso_1 = [oxygen_masses(1), oxygen_masses(3), oxygen_masses(1)]
    masses_iso_2 = [argon_mass]
    masses_ref_1 = [oxygen_masses(1), oxygen_masses(1), oxygen_masses(1)]
    masses_ref_2 = [argon_mass]
    zx_atom_index_1 = [2, 3]
    zx_atom_index_2 = [0, 0]

    cartesian_iso_1 = normalize_cartesian_matrix(cartesian_input_1, masses_iso_1, zx_atom_index_1)
    cartesian_iso_2 = normalize_cartesian_matrix(cartesian_input_2, masses_iso_2, zx_atom_index_2)

    do i = 1, size(coords_iso, 2)
      coords_iso_full = [coords_iso(1, i), 0d0, coords_iso(2, i), coords_iso(3, i), 0d0, 0d0]
      coords_ref_full = convert_coordinates_to_reference(cartesian_iso_1, cartesian_iso_2, coords_iso_full, masses_ref_1, masses_ref_2, zx_atom_index_1, zx_atom_index_2)
      coords_ref(:, i) = [coords_ref_full(1), coords_ref_full(3), coords_ref_full(4)]
    end do
  end function

end program
