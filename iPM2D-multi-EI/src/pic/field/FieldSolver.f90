Module ModuleFieldSolver
    use ModuleControlFlow
    use ModuleFileName
    use ModuleFieldBoundary2DZR
    use ModuleFieldSource
    use ModuleDomain
    use ModuleGeometry
    use ModuleMaterials
    use ModuleExtCircuit

    implicit none

    type FieldSolver
        type(Domain), Pointer :: DM => Null()
        integer(4) :: laplace_num

        real(8), allocatable :: Source(:, :)
        real(8), allocatable :: Solve(:, :)
        real(8), allocatable :: CoeA(:, :), CoeB(:, :), CoeC(:, :), CoeD(:, :), CoeE(:, :)

        real(8), allocatable :: poisson(:, :)
        real(8), allocatable :: laplace(:, :, :)

    contains

        procedure :: Init => InitFieldSolver
        procedure :: Zero => ZeroFieldSolver
        procedure :: Reset => ResetFieldSolver
        procedure :: Destroy => DestroyFieldSolver
        procedure :: SolveOne => SolverPossionOne
        procedure :: Diag => DiagSolveResult
        procedure :: SetBound => SetSolverBoundary
        procedure :: SetSrc => SetFieldSolverSource
        procedure :: SetCoe => SetFieldSolverCoe
        procedure :: SetCoeWN => SetFieldSolverCoeWithoutChi

        procedure :: SolveV => SolverPossionVahedi
        procedure :: SolveL => SolverPossionLaplace
        procedure :: SolveP => SolverPossionZero
        procedure :: SolveE => SolverPossion

    end type FieldSolver

    contains

        subroutine InitFieldSolver(FS, CF, metal_num)
            class(FieldSolver),intent(inout) :: FS
            type(ControlFlow), intent(in) :: CF
            integer(4), intent(in) :: metal_num

            FS%laplace_num = metal_num
            call FS%Reset(CF%DM)

            ! create all petsc solvers
            call init_all_solvers()

            ! init possion solver
            associate (Nz  => CF%DM%GlobalShape(1), Nr         => CF%DM%GlobalShape(2), &
                       dz  => CF%DM%SpaceStep(1),   dr         => CF%DM%SpaceStep(2), &
                       NPz => CF%DM%ImageShape(1),  zShapeList => CF%DM%LocalShapeList(1, :), &
                       NPr => CF%DM%ImageShape(2),  rShapeList => CF%DM%LocalShapeList(2, :))
            
                call init_poisson_solver(Nz, Nr, dz, dr, NPz, zShapeList, NPr, rShapeList)

            end associate 

        end subroutine InitFieldSolver


        subroutine ResetFieldSolver(FS, DM)
            class(FieldSolver),intent(inout) :: FS
            type(Domain), intent(in), target :: DM

            FS%DM => DM
            call FS%Destroy()

            associate(zstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 2))

                allocate(FS%Source(zstart : zend, rstart : rend))
                allocate(FS%Solve(zstart : zend, rstart : rend))
                allocate(FS%CoeA(zstart : zend, rstart : rend))
                allocate(FS%CoeB(zstart : zend, rstart : rend))
                allocate(FS%CoeC(zstart : zend, rstart : rend))
                allocate(FS%CoeD(zstart : zend, rstart : rend))
                allocate(FS%CoeE(zstart : zend, rstart : rend))

                allocate(FS%poisson(zstart : zend, rstart : rend))
                allocate(FS%laplace(zstart : zend, rstart : rend, 1:FS%laplace_num))

            end associate

            call FS%Zero()

        end subroutine ResetFieldSolver


        subroutine ZeroFieldSolver(FS)
            class(FieldSolver),intent(inout) :: FS

            FS%Source = 0.d0
            FS%Solve = 0.d0
            FS%CoeA = 0.d0
            FS%CoeB = 0.d0
            FS%CoeC = 0.d0
            FS%CoeD = 0.d0
            FS%CoeE = 0.d0

            FS%poisson = 0.d0
            FS%laplace = 0.d0

        end subroutine ZeroFieldSolver


        subroutine DestroyFieldSolver(FS)
            class(FieldSolver),intent(inout) :: FS

            if (allocated(FS%Source)) deallocate(FS%Source)
            if (allocated(FS%Solve)) deallocate(FS%Solve)
            if (allocated(FS%CoeA)) deallocate(FS%CoeA)
            if (allocated(FS%CoeB)) deallocate(FS%CoeB)
            if (allocated(FS%CoeC)) deallocate(FS%CoeC)
            if (allocated(FS%CoeD)) deallocate(FS%CoeD)
            if (allocated(FS%CoeE)) deallocate(FS%CoeE)

            if (allocated(FS%poisson)) deallocate(FS%poisson)
            if (allocated(FS%laplace)) deallocate(FS%laplace)

            call destroy_all_solvers()

        end subroutine DestroyFieldSolver


        subroutine SetSolverBoundary(FS, boundary)
            class(FieldSolver), intent(inout) :: FS
            type(FieldBoundary2DZR), intent(in) :: boundary

            associate(Nz     => FS%DM%GlobalShape(1), Nr => FS%DM%GlobalShape(2), &
                      zstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend   => FS%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend   => FS%DM%CornerIndex(BOUNDARY_RIGHT, 2))
                
                FS%Source(zstart:zend, rstart:rend) = merge(boundary%bound(zstart:zend, rstart:rend), FS%Source(zstart:zend, rstart:rend), boundary%bound_mask(zstart:zend, rstart:rend) /= 0)

            end associate

        end subroutine SetSolverBoundary


        subroutine SetFieldSolverSource(FS, FT, MT, Geom)
            class(FieldSolver), intent(inout) :: FS
            type(FieldSource), intent(in) :: FT
            type(Materials), intent(in) :: MT
            type(Geometry), intent(in) :: Geom
            integer(4) :: i, j, point_type_one
            real(8) :: inv_volume

            associate(zstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 2), &
                      dx => FS%DM%SpaceStep(1))

                FS%Source = -1.d0 * dx * dx / Epsilon * (FT%Rho + MT%sigma(zstart:zend, rstart:rend) * Geom%point_inv_volume_weighting(zstart:zend, rstart:rend))

            end associate

        end subroutine SetFieldSolverSource


        subroutine SetFieldSolverCoe(FS, FT, Geom)
            class(FieldSolver), intent(inout) :: FS
            type(FieldSource), intent(in) :: FT
            type(Geometry), intent(in) :: Geom
            integer(4) :: i, j, point_type_one

            associate(Nz => FS%DM%GlobalShape(1), Nr => FS%DM%GlobalShape(2), &
                      zstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 2), &
                      coeA => FS%CoeA, coeB => FS%CoeB, coeC => FS%CoeC, &
                      coeD => FS%CoeD, coeE => FS%CoeE, Chi => FT%Chi, &
                      epsilon_r => Geom%cell_epsilon)

                coeA = 0.d0
                coeB = 0.d0
                coeC = 0.d0
                coeD = 0.d0
                coeE = 0.d0

                do j = rstart, rend                                 ! set coe and source for RZ
                    do i = zstart, zend
                        point_type_one = Geom%point_type(i, j)

                        if (point_type_one > 0) then                    ! 金属介质
                            coeB(i, j) = 1.d0

                        else
                            if (Nr == j .or. 1 == i .or. Nz == i) then  ! 边界 the value is not use           

                            else if (1 == j) then                       ! 中轴
                                ! coeA(i, j) =  1.d0
                                ! coeC(i, j) =  1.d0
                                ! coeD(i, j) =  0.d0
                                ! coeE(i, j) =  4.d0
                                ! coeB(i, j) = -(coeA(i, j) + coeC(i, j) + coeD(i, j) + coeE(i, j))

                                ! coeA(i, j) =  1.d0 + 0.5d0 * (Chi(i, j) + Chi(i-1, j))
                                ! coeC(i, j) =  1.d0 + 0.5d0 * (Chi(i, j) + Chi(i+1, j))
                                ! coeD(i, j) =  0.d0
                                ! coeE(i, j) =  4.d0 * (1.d0 + 0.5d0 * (Chi(i, j) + Chi(i, j+1)))
                                ! coeB(i, j) = -(coeA(i, j) + coeC(i, j) + coeD(i, j) + coeE(i, j))

                                coeA(i, j) =  (1.d0 + 0.5d0 * (Chi(i, j) + Chi(i-1, j))) * 0.5d0 * (epsilon_r(i, j) + epsilon_r(i-1, j))
                                coeC(i, j) =  (1.d0 + 0.5d0 * (Chi(i, j) + Chi(i+1, j))) * 0.5d0 * (epsilon_r(i, j) + epsilon_r(i+1, j))
                                coeE(i, j) =  4.d0 * (1.d0 + 0.5d0 * (Chi(i, j) + Chi(i, j+1))) * 0.5d0 * (epsilon_r(i, j) + epsilon_r(i, j+1))
                                coeB(i, j) = -(coeA(i, j) + coeC(i, j) + coeD(i, j) + coeE(i, j))

                            else
                                ! coeA(i, j) =  1.d0
                                ! coeC(i, j) =  1.d0
                                ! coeD(i, j) =  (j - 1.5d0) / (j - 1.d0)
                                ! coeE(i, j) =  (j - 0.5d0) / (j - 1.d0)
                                ! coeB(i, j) = -(coeA(i, j) + coeC(i, j) + coeD(i, j) + coeE(i, j))

                                ! coeA(i, j) =   1.d0 + 0.5d0 * (Chi(i, j) + Chi(i-1, j))
                                ! coeC(i, j) =   1.d0 + 0.5d0 * (Chi(i, j) + Chi(i+1, j))
                                ! coeD(i, j) =  (1.d0 + 0.5d0 * (Chi(i, j) + Chi(i, j-1))) * (j - 1.5d0) / (j - 1.d0)
                                ! coeE(i, j) =  (1.d0 + 0.5d0 * (Chi(i, j) + Chi(i, j+1))) * (j - 0.5d0) / (j - 1.d0)
                                ! coeB(i, j) = -(coeA(i, j) + coeC(i, j) + coeD(i, j) + coeE(i, j))

                                coeA(i, j) =  (1.d0 + 0.5d0 * (Chi(i, j) + Chi(i-1, j))) * 0.5d0 * (epsilon_r(i-1, j-1) + epsilon_r(i-1, j))
                                coeC(i, j) =  (1.d0 + 0.5d0 * (Chi(i, j) + Chi(i+1, j))) * 0.5d0 * (epsilon_r(i, j-1) + epsilon_r(i, j))
                                coeD(i, j) =  (1.d0 + 0.5d0 * (Chi(i, j) + Chi(i, j-1))) * 0.5d0 * (epsilon_r(i-1, j-1) + epsilon_r(i, j-1)) * (j - 1.5d0) / (j - 1.d0)
                                coeE(i, j) =  (1.d0 + 0.5d0 * (Chi(i, j) + Chi(i, j+1))) * 0.5d0 * (epsilon_r(i-1, j) + epsilon_r(i, j)) * (j - 0.5d0) / (j - 1.d0)
                                coeB(i, j) = -(coeA(i, j) + coeC(i, j) + coeD(i, j) + coeE(i, j))

                            end if
                        end if
                    end do
                end do
            end associate

        end subroutine SetFieldSolverCoe


        subroutine SetFieldSolverCoeWithoutChi(FS, Geom)
            class(FieldSolver), intent(inout) :: FS
            type(Geometry), intent(in) :: Geom
            integer(4) :: i, j, point_type_one

            associate(Nz => FS%DM%GlobalShape(1), Nr => FS%DM%GlobalShape(2), &
                      zstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 2), &
                      coeA => FS%CoeA, coeB => FS%CoeB, coeC => FS%CoeC, &
                      coeD => FS%CoeD, coeE => FS%CoeE, &
                      epsilon_r => Geom%cell_epsilon)

                coeA = 0.d0
                coeB = 0.d0
                coeC = 0.d0
                coeD = 0.d0
                coeE = 0.d0

                do j = rstart, rend                                 ! set coe and source for RZ
                    do i = zstart, zend
                        point_type_one = Geom%point_type(i, j)

                        if (point_type_one > 0) then                    ! 金属介质
                            coeB(i, j) = 1.d0

                        else

                            if (Nr == j .or. 1 == i .or. Nz == i) then  ! 边界 the value is not use

                            else if (1 == j) then                       ! 中轴
                                coeA(i, j) =  0.5d0 * (epsilon_r(i, j) + epsilon_r(i-1, j))
                                coeC(i, j) =  0.5d0 * (epsilon_r(i, j) + epsilon_r(i+1, j))
                                coeE(i, j) =  4.d0 * 0.5d0 * (epsilon_r(i, j) + epsilon_r(i, j+1))
                                coeB(i, j) = -(coeA(i, j) + coeC(i, j) + coeD(i, j) + coeE(i, j))

                            else
                                coeA(i, j) =  0.5d0 * (epsilon_r(i-1, j-1) + epsilon_r(i-1, j))
                                coeC(i, j) =  0.5d0 * (epsilon_r(i, j-1) + epsilon_r(i, j))
                                coeD(i, j) =  0.5d0 * (epsilon_r(i-1, j-1) + epsilon_r(i, j-1)) * (j - 1.5d0) / (j - 1.d0)
                                coeE(i, j) =  0.5d0 * (epsilon_r(i-1, j) + epsilon_r(i, j)) * (j - 0.5d0) / (j - 1.d0)
                                coeB(i, j) = -(coeA(i, j) + coeC(i, j) + coeD(i, j) + coeE(i, j))

                            end if
                        end if
                    end do
                end do
            end associate

        end subroutine SetFieldSolverCoeWithoutChi


        subroutine SolverPossionOne(FS, result)
            class(FieldSolver),intent(inout) :: FS
            real(8), optional, allocatable, intent(inout) :: result(:, :)

            if (present(result)) then
                result = 0.d0
                call poisson_solve(FS%CoeA, FS%CoeB, FS%CoeC, FS%CoeD, FS%CoeE, FS%Source, result)

            else
                FS%Solve = 0.d0
                call poisson_solve(FS%CoeA, FS%CoeB, FS%CoeC, FS%CoeD, FS%CoeE, FS%Source, FS%Solve)

            end if

        end subroutine SolverPossionOne


        subroutine SolverPossionVahedi(FS, CF, FT, Ext, boundary, MT, Geom)
            class(FieldSolver), intent(inout) :: FS
            type(ControlFlow), intent(inout) :: CF
            type(FieldSource), intent(inout) :: FT
            type(ExtCircuits), intent(inout) :: Ext
            type(FieldBoundary2DZR), intent(inout) :: boundary
            type(Materials), intent(inout) :: MT
            type(Geometry), intent(inout) :: Geom
            integer(4) :: i

            if (load_type_plasma == load_type) then
                call FS%SolveL(CF, FT, Ext, boundary, MT, Geom)
                call FS%SolveP(CF, FT, Ext, boundary, MT, Geom)
            end if

            call calq0(FS, FT, MT, Geom, Ext)
            call Ext%update(dble(CF%Timer) * CF%dt, MT)

            FS%Solve = FS%poisson
            if (MT%metal_count > 0) then 
                do i = 1, MT%metal_count
                    if (Ext%circuits_type(i) /= circuit_type_digital_ground) then
                        FS%Solve = FS%Solve + FS%laplace(:, :, i) * MT%metls(i)%voltage
                    end if
                end do

            end if

            call calEd(FS, FT, MT, Geom, Ext)

        end subroutine SolverPossionVahedi


        subroutine SolverPossionLaplace(FS, CF, FT, Ext, boundary, MT, Geom)
            class(FieldSolver), intent(inout) :: FS
            type(ControlFlow), intent(inout) :: CF
            type(FieldSource), intent(inout) :: FT
            type(ExtCircuits), intent(inout) :: Ext
            type(FieldBoundary2DZR), intent(inout) :: boundary
            type(Materials), intent(inout) :: MT
            type(Geometry), intent(inout) :: Geom
            integer(4) :: i

            if (MT%metal_count > 0) then
                FS%laplace = 0.d0
                call FS%SetCoe(FT, Geom)

                do i = 1, MT%metal_count
                    if (Ext%circuits_type(i) /= circuit_type_digital_ground) then
                        call boundary%Zero()
                        call boundary%Update(MT, Geom)
                        boundary%bound = 0.d0
                        boundary%bound = merge(1.d0, boundary%bound, Geom%point_type == i)
                        call boundary%intep(MT, Geom)

                        FS%Source = 0.d0
                        call FS%SetBound(boundary)
                        call FS%SolveOne()
                        FS%laplace(:, :, i) = FS%Solve
                    end if
                end do

            end if

        end subroutine SolverPossionLaplace


        subroutine SolverPossionZero(FS, CF, FT, Ext, boundary, MT, Geom)
            class(FieldSolver), intent(inout) :: FS
            type(ControlFlow), intent(inout) :: CF
            type(FieldSource), intent(inout) :: FT
            type(ExtCircuits), intent(inout) :: Ext
            type(FieldBoundary2DZR), intent(inout) :: boundary
            type(Materials), intent(inout) :: MT
            type(Geometry), intent(inout) :: Geom

            call boundary%Zero()
            call boundary%Update(MT, Geom)
            boundary%bound = 0.d0
            call boundary%intep(MT, Geom)

            call FS%SetCoe(FT, Geom)
            call FS%SetSrc(FT, MT, Geom)
            call FS%SetBound(boundary)
            call FS%SolveOne(FS%poisson)

        end subroutine SolverPossionZero


        subroutine SolverPossion(FS, CF, FT, Ext, boundary, MT, Geom)
            class(FieldSolver), intent(inout) :: FS
            type(ControlFlow), intent(inout) :: CF
            type(FieldSource), intent(inout) :: FT
            type(ExtCircuits), intent(inout) :: Ext
            type(FieldBoundary2DZR), intent(inout) :: boundary
            type(Materials), intent(inout) :: MT
            type(Geometry), intent(inout) :: Geom

            call Ext%update(dble(CF%Timer) * CF%dt, MT)
            call boundary%Zero()
            call boundary%Update(MT, Geom)
            call boundary%intep(MT, Geom)

            call FS%SetCoe(FT, Geom)
            call FS%SetSrc(FT, MT, Geom)
            call FS%SetBound(boundary)
            call FS%SolveOne()

        end subroutine SolverPossion


        subroutine DiagSolveResult(FS, name)
            class(FieldSolver),intent(inout) :: FS
            Character(*), intent(in) :: name
            type(HDF5_PDump) :: hdf5Dump
            type(FileName) 	:: 	IOName
            integer(4) :: zend_on_bound, rend_on_bound

            associate(dz => FS%DM%SpaceStep(1), dr => FS%DM%SpaceStep(2), &
                      Nz => FS%DM%GlobalShape(1), Nr => FS%DM%GlobalShape(2), &
                      Lz => FS%DM%LocalShape(1), Lr => FS%DM%LocalShape(2), &
                      zstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 2), &
                      zRightImageType => FS%DM%NeighborType(BOUNDARY_RIGHT, 1), &
                      rRightImageType => FS%DM%NeighborType(BOUNDARY_RIGHT, 2), &
                      image_shape => FS%DM%ImageShape)

                call IOName%Init(name, PathMode=FILE_PATH_MODE_DIAG, &
                                       ExtensionMode=FILE_EXTENSION_MODE_H5, &
                                       ParallelMode=FILE_PARALLEL_MODE_CLOSE, &
                                       DynamicIndex=FILE_DYNAMIC_MODE_CLOSE)

                call hdf5Dump%init(filename=IOName%FullName%str, mode='write')
                call hdf5Dump%open()

                call hdf5Dump%writeattr('/', 'Nz', Nz)
                call hdf5Dump%writeattr('/', 'Nr', Nr)

                zend_on_bound = zend
                rend_on_bound = rend
                if (image_shape(1) > 1 .and. NEIGHBOR_TYPE_DOMAIN == zRightImageType) zend_on_bound = zend - 1
                if (image_shape(2) > 1 .and. NEIGHBOR_TYPE_DOMAIN == rRightImageType) rend_on_bound = rend - 1
        
                call hdf5Dump%write('/Result', FS%Solve(zstart:zend_on_bound, rstart:rend_on_bound), &
                                    chunkdim=image_shape)
                call hdf5Dump%write('/Resultall', FS%Solve(zstart:zend, rstart:rend), chunkdim=image_shape)
                call hdf5Dump%write('/CoeA', FS%CoeA(zstart:zend_on_bound, rstart:rend_on_bound), &
                                    chunkdim=image_shape)
                call hdf5Dump%write('/CoeB', FS%CoeB(zstart:zend_on_bound, rstart:rend_on_bound), &
                                    chunkdim=image_shape)
                call hdf5Dump%write('/CoeC', FS%CoeC(zstart:zend_on_bound, rstart:rend_on_bound), &
                                    chunkdim=image_shape)
                call hdf5Dump%write('/CoeD', FS%CoeD(zstart:zend_on_bound, rstart:rend_on_bound), &
                                    chunkdim=image_shape)
                call hdf5Dump%write('/CoeE', FS%CoeE(zstart:zend_on_bound, rstart:rend_on_bound), &
                                    chunkdim=image_shape)
                call hdf5Dump%write('/Source', FS%Source(zstart:zend_on_bound, rstart:rend_on_bound), &
                                    chunkdim=image_shape)

                call hdf5Dump%close()

            end associate

        end subroutine DiagSolveResult


        subroutine calq0(FS, FT, MT, Geom, Ext)
            class(FieldSolver), intent(inout) :: FS
            type(FieldSource), intent(in) :: FT
            type(Materials), intent(inout) :: MT
            type(Geometry), intent(in) :: Geom
            type(ExtCircuits), intent(inout) :: Ext
            real(8), allocatable :: pg(:, :), pl(:, :), chig(:, :), chil(:, :)
            integer(4) :: zend_on_bound, rend_on_bound
            integer(4) :: i, k, ierr, iz, ir, id
            real(8) :: vdz, vdr, E1, E2, E3

            associate(Nz => FS%DM%GlobalShape(1), Nr => FS%DM%GlobalShape(2), &
                      zstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 2), &
                      dz => FS%DM%SpaceStep(1), dr => FS%DM%SpaceStep(2), &
                      zRightImageType => FS%DM%NeighborType(BOUNDARY_RIGHT, 1), &
                      rRightImageType => FS%DM%NeighborType(BOUNDARY_RIGHT, 2), &
                      image_shape => FS%DM%ImageShape)

                allocate(pg(Nz, Nr))
                allocate(pl(Nz, Nr))
                allocate(chig(Nz, Nr))
                allocate(chil(Nz, Nr))

                vdz = 1.d0 / dz
                vdr = 1.d0 / dr

                ! com
                zend_on_bound = zend
                rend_on_bound = rend
                if (image_shape(1) > 1 .and. NEIGHBOR_TYPE_DOMAIN == zRightImageType) zend_on_bound = zend - 1
                if (image_shape(2) > 1 .and. NEIGHBOR_TYPE_DOMAIN == rRightImageType) rend_on_bound = rend - 1

                pl = 0.d0
                pl(zstart:zend_on_bound, rstart:rend_on_bound) = FS%poisson(zstart:zend_on_bound, rstart:rend_on_bound)

                chil = 0.d0
                chil(zstart:zend_on_bound, rstart:rend_on_bound) = FT%Chi(zstart:zend_on_bound, rstart:rend_on_bound)

                call MPI_ALLREDUCE(pl, pg, Nz*Nr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
                call MPI_ALLREDUCE(chil, chig, Nz*Nr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)

                ! cal q0
                do k = 1, MT%metal_count
                    if (Ext%circuits_type(k) >= circuit_type_rlc_source) then
                        MT%metls(k)%q0 = 0.d0

                        if (MT%BO(k)%lv > 0) then       ! v
                            do i = 1, MT%BO(k)%lv
                                iz = MT%BO(k)%iv(i, 1)
                                ir = MT%BO(k)%iv(i, 2)
                                id = MT%BO(k)%iv(i, 3)

                                if (id == 1) then
                                    E2 = ((pg(iz, ir+1) + pg(iz+1, ir+1)) * 0.5d0 - (pg(iz, ir+2) + pg(iz+1, ir+2)) * 0.5d0) * vdr
                                    E3 = ((pg(iz, ir+2) + pg(iz+1, ir+2)) * 0.5d0 - (pg(iz, ir+3) + pg(iz+1, ir+3)) * 0.5d0) * vdr
                                    E1 = (E2 + 0.5d0 * (E2 - E3)) * (1 + 0.5d0 * (chig(iz, ir+1) + chig(iz+1, ir+1))) * Geom%cell_epsilon(iz, ir+1) * dble(ir)

                                else if (id == -1) then
                                    E2 = ((pg(iz, ir-1) + pg(iz+1, ir-1)) * 0.5d0 - (pg(iz, ir) + pg(iz+1, ir)) * 0.5d0) * vdr
                                    E3 = ((pg(iz, ir-2) + pg(iz+1, ir-2)) * 0.5d0 - (pg(iz, ir-1) + pg(iz+1, ir-1)) * 0.5d0) * vdr
                                    E1 = (E2 + 0.5d0 * (E2 - E3)) * (1 + 0.5d0 * (chig(iz, ir) + chig(iz+1, ir))) * Geom%cell_epsilon(iz, ir-1) * dble(ir-1)

                                end if

                                MT%metls(k)%q0 = MT%metls(k)%q0 + E1 * dble(id) * 2.d0 * PI * dr * dz
                            end do
                        end if

                        if (MT%BO(k)%lh > 0) then       ! h
                            do i = 1, MT%BO(k)%lh
                                iz = MT%BO(k)%ih(i, 1)
                                ir = MT%BO(k)%ih(i, 2)
                                id = MT%BO(k)%ih(i, 3)

                                if (id == 1) then
                                    E2 = ((pg(iz+1, ir) + pg(iz+1, ir+1)) * 0.5d0 - (pg(iz+2, ir) + pg(iz+2, ir+1)) * 0.5d0) * vdz
                                    E3 = ((pg(iz+2, ir) + pg(iz+2, ir+1)) * 0.5d0 - (pg(iz+3, ir) + pg(iz+3, ir+1)) * 0.5d0) * vdz
                                    E1 = (E2 + 0.5d0 * (E2 - E3)) * (1 + 0.5d0 * (chig(iz+1, ir) + chig(iz+1, ir+1))) * Geom%cell_epsilon(iz+1, ir)

                                else if (id == -1) then
                                    E2 = ((pg(iz-1, ir) + pg(iz-1, ir+1)) * 0.5d0 - (pg(iz, ir) + pg(iz, ir+1)) * 0.5d0) * vdz
                                    E3 = ((pg(iz-2, ir) + pg(iz-2, ir+1)) * 0.5d0 - (pg(iz-1, ir) + pg(iz-1, ir+1)) * 0.5d0) * vdz
                                    E1 = (E2 + 0.5d0 * (E2 - E3)) * (1 + 0.5d0 * (chig(iz, ir) + chig(iz, ir+1))) * Geom%cell_epsilon(iz-1, ir)

                                end if

                                MT%metls(k)%q0 = MT%metls(k)%q0 + E1 * dble(id) * PI * dr * dr * ((ir*ir) - (ir-1)*(ir-1))
                            end do
                        end if

                        MT%metls(k)%q0 = MT%metls(k)%q0 * 8.8542d-12
                        
                    end if
                end do

                deallocate(pg)
                deallocate(pl)
                deallocate(chig)
                deallocate(chil)

            end associate

        end subroutine calq0


        subroutine calEd(FS, FT, MT, Geom, Ext)
            class(FieldSolver), intent(inout) :: FS
            type(FieldSource), intent(in) :: FT
            type(Materials), intent(inout) :: MT
            type(Geometry), intent(in) :: Geom
            type(ExtCircuits), intent(inout) :: Ext
            real(8), allocatable :: pg(:, :), pl(:, :), chig(:, :), chil(:, :)
            integer(4) :: zend_on_bound, rend_on_bound
            integer(4) :: i, k, ierr, iz, ir, id
            real(8) :: vdz, vdr, E1, E2, E3
            real(8) :: Ed(MT%metal_count)

            associate(Nz => FS%DM%GlobalShape(1), Nr => FS%DM%GlobalShape(2), &
                      zstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 1), &
                      zend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 1), &
                      rstart => FS%DM%CornerIndex(BOUNDARY_LEFT, 2), &
                      rend => FS%DM%CornerIndex(BOUNDARY_RIGHT, 2), &
                      dz => FS%DM%SpaceStep(1), dr => FS%DM%SpaceStep(2), &
                      zRightImageType => FS%DM%NeighborType(BOUNDARY_RIGHT, 1), &
                      rRightImageType => FS%DM%NeighborType(BOUNDARY_RIGHT, 2), &
                      image_shape => FS%DM%ImageShape)

                allocate(pg(Nz, Nr))
                allocate(pl(Nz, Nr))
                allocate(chig(Nz, Nr))
                allocate(chil(Nz, Nr))

                vdz = 1.d0 / dz
                vdr = 1.d0 / dr

                ! com
                zend_on_bound = zend
                rend_on_bound = rend
                if (image_shape(1) > 1 .and. NEIGHBOR_TYPE_DOMAIN == zRightImageType) zend_on_bound = zend - 1
                if (image_shape(2) > 1 .and. NEIGHBOR_TYPE_DOMAIN == rRightImageType) rend_on_bound = rend - 1

                pl = 0.d0
                pl(zstart:zend_on_bound, rstart:rend_on_bound) = FS%Solve(zstart:zend_on_bound, rstart:rend_on_bound)

                chil = 0.d0
                chil(zstart:zend_on_bound, rstart:rend_on_bound) = FT%Chi(zstart:zend_on_bound, rstart:rend_on_bound)

                call MPI_ALLREDUCE(pl, pg, Nz*Nr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
                call MPI_ALLREDUCE(chil, chig, Nz*Nr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)

                ! cal q0
                Ed = 0.d0
                do k = 1, MT%metal_count
                    if (Ext%circuits_type(k) >= circuit_type_rlc_source) then
                        Ed(k) = 0.d0

                        if (MT%BO(k)%lv > 0) then       ! v
                            do i = 1, MT%BO(k)%lv
                                iz = MT%BO(k)%iv(i, 1)
                                ir = MT%BO(k)%iv(i, 2)
                                id = MT%BO(k)%iv(i, 3)

                                if (id == 1) then
                                    E2 = ((pg(iz, ir+1) + pg(iz+1, ir+1)) * 0.5d0 - (pg(iz, ir+2) + pg(iz+1, ir+2)) * 0.5d0) * vdr
                                    E3 = ((pg(iz, ir+2) + pg(iz+1, ir+2)) * 0.5d0 - (pg(iz, ir+3) + pg(iz+1, ir+3)) * 0.5d0) * vdr
                                    E1 = (E2 + 0.5d0 * (E2 - E3)) * (1 + 0.5d0 * (chig(iz, ir+1) + chig(iz+1, ir+1))) * Geom%cell_epsilon(iz, ir+1) * dble(ir)

                                else if (id == -1) then
                                    E2 = ((pg(iz, ir-1) + pg(iz+1, ir-1)) * 0.5d0 - (pg(iz, ir) + pg(iz+1, ir)) * 0.5d0) * vdr
                                    E3 = ((pg(iz, ir-2) + pg(iz+1, ir-2)) * 0.5d0 - (pg(iz, ir-1) + pg(iz+1, ir-1)) * 0.5d0) * vdr
                                    E1 = (E2 + 0.5d0 * (E2 - E3)) * (1 + 0.5d0 * (chig(iz, ir) + chig(iz+1, ir))) * Geom%cell_epsilon(iz, ir-1) * dble(ir-1)

                                end if

                                Ed(k) = Ed(k) + E1 * dble(id) * 2.d0 * PI * dr * dz
                            end do
                        end if

                        if (MT%BO(k)%lh > 0) then       ! h
                            do i = 1, MT%BO(k)%lh
                                iz = MT%BO(k)%ih(i, 1)
                                ir = MT%BO(k)%ih(i, 2)
                                id = MT%BO(k)%ih(i, 3)

                                if (id == 1) then
                                    E2 = ((pg(iz+1, ir) + pg(iz+1, ir+1)) * 0.5d0 - (pg(iz+2, ir) + pg(iz+2, ir+1)) * 0.5d0) * vdz
                                    E3 = ((pg(iz+2, ir) + pg(iz+2, ir+1)) * 0.5d0 - (pg(iz+3, ir) + pg(iz+3, ir+1)) * 0.5d0) * vdz
                                    E1 = (E2 + 0.5d0 * (E2 - E3)) * (1 + 0.5d0 * (chig(iz+1, ir) + chig(iz+1, ir+1))) * Geom%cell_epsilon(iz+1, ir)

                                else if (id == -1) then
                                    E2 = ((pg(iz-1, ir) + pg(iz-1, ir+1)) * 0.5d0 - (pg(iz, ir) + pg(iz, ir+1)) * 0.5d0) * vdz
                                    E3 = ((pg(iz-2, ir) + pg(iz-2, ir+1)) * 0.5d0 - (pg(iz-1, ir) + pg(iz-1, ir+1)) * 0.5d0) * vdz
                                    E1 = (E2 + 0.5d0 * (E2 - E3)) * (1 + 0.5d0 * (chig(iz, ir) + chig(iz, ir+1))) * Geom%cell_epsilon(iz-1, ir)

                                end if

                                Ed(k) = Ed(k) + E1 * dble(id) * PI * dr * dr * ((ir*ir) - (ir-1)*(ir-1))
                            end do
                        end if

                        Ed(k) = Ed(k) * 8.8542d-12
                        
                    end if
                end do

                deallocate(pg)
                deallocate(pl)
                deallocate(chig)
                deallocate(chil)

            end associate

            if (0 == FS%DM%MyId) then
                open(10, file="./Ed.txt", position='append')
                    write(10, '(*(es20.12, 1x))') Ed
                close(10)
            end if

        end subroutine calEd

end Module ModuleFieldSolver
