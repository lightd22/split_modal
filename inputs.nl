&inputs
    ! Spatial element parameters
    startRes = 15,      ! Which resolution is run first
    nRuns = 3,          ! How many runs are done
    nScale = 2,         ! Ratio between number of elements in successive runs
    polyDegree = 3,     ! Degree of reconstructing polynomial

    ! Time stepping paramteters
    cflCoeff = 0.9D0    ! Ratio of used CFL number to maximum stable CFL

    ! Outputting parameters
    noutput = 2         ! Number of times to output output, including final time (must be >= 1) (automatically includes ICs)

    ! Testing parmeters
    whichTest = 5       ! 0 = Consistency test
                        ! 1 = Uniform diagonal advection
                        ! 2 = Solid body rotation of cylinder
                        ! 5 = LeVeque deformation of C^3 cosinebell
                        ! 6 = LeVeque deformation of C^5 cosinebell
                        ! 7 = LeVeque deformation of slotted cylinder

    ! Misc parameters
    DEBUG = .FALSE.
    DOTIMETEST = .FALSE.

/

