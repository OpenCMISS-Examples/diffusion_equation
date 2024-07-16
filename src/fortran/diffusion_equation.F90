PROGRAM DiffusionEquation

  USE OpenCMISS

  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------
  
  !Test program parameters
  REAL(OC_RP), PARAMETER :: HEIGHT=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: WIDTH=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: LENGTH=3.0_OC_RP
  
  INTEGER(OC_Intg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: REGION_USER_NUMBER=3
  INTEGER(OC_Intg), PARAMETER :: BASIS_USER_NUMBER=4
  INTEGER(OC_Intg), PARAMETER :: GENERATED_MESH_USER_NUMBER=5
  INTEGER(OC_Intg), PARAMETER :: MESH_USER_NUMBER=6
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSITION_USER_NUMBER=7
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSER_USER_NUMBER=8
  INTEGER(OC_Intg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=9
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=10
  INTEGER(OC_Intg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=11
  INTEGER(OC_Intg), PARAMETER :: MATERIALS_FIELD_USER_NUMBER=12
  INTEGER(OC_Intg), PARAMETER :: ANALYTIC_FIELD_USER_NUMBER=13
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=14
  INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER=15

  !Program types
  
  !Program variables


  !INTEGER(INTG) :: first_global_dof,first_local_dof,first_local_rank,last_global_dof,last_local_dof,last_local_rank,rank_idx
  !INTEGER(INTG) :: EQUATIONS_SET_INDEX !TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOF_MAPPING
  
  !CMISS variables
  TYPE(OC_BasisType) :: basis
  TYPE(OC_BoundaryConditionsType) :: boundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_CoordinateSystemType) :: coordinateSystem
  TYPE(OC_DecompositionType) :: decomposition
  TYPE(OC_DecomposerType) :: decomposer
  TYPE(OC_EquationsType) :: equations
  TYPE(OC_EquationsSetType) :: equationsSet
  TYPE(OC_FieldType) :: geometricField,equationsSetField,dependentField,materialsField,analyticField
  TYPE(OC_FieldsType) :: fields
  TYPE(OC_GeneratedMeshType) :: generatedMesh  
  TYPE(OC_MeshType) :: mesh
  TYPE(OC_ProblemType) :: problem
  TYPE(OC_ControlLoopType) :: controlLoop
  TYPE(OC_RegionType) :: region,worldRegion
  TYPE(OC_SolverType) :: solver,linearSolver
  TYPE(OC_SolverEquationsType) :: solverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  INTEGER(OC_Intg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements
  INTEGER(OC_Intg) :: computationalNodeNumber,numberOfComputationalNodes
  INTEGER(OC_Intg) :: decompositionIndex,equationsSetIndex
  INTEGER(OC_Intg) :: err
  LOGICAL :: exportField

  !STOP

  !Intialise OpenCMISS
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)

  !Get the computational nodes information
  CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  numberOfGlobalXElements=10
  numberOfGlobalYElements=10
  numberOfGlobalZElements=0

  !-----------------------------------------------------------------------------------------------------------
  ! COORDINATE SYSTEM
  !----------------------------------------------------------------------------------------------------------- 
  
  !Start the creation of a new RC coordinate system
  CALL OC_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  IF(numberOfGlobalZElements==0) THEN
    !Set the coordinate system to be 2D
    CALL OC_CoordinateSystem_DimensionSet(coordinateSystem,2,err)
  ELSE
    !Set the coordinate system to be 3D
    CALL OC_CoordinateSystem_DimensionSet(coordinateSystem,3,err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL OC_CoordinateSystem_CreateFinish(coordinateSystem,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------
 
  !Start the creation of the region
  CALL OC_Region_Initialise(region,err)
  CALL OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL OC_Region_LabelSet(region,"diffusion_equation",err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL OC_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL OC_Region_CreateFinish(region,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! BASIS
  !----------------------------------------------------------------------------------------------------------  
  
  !Start the creation of a basis (default is trilinear lagrange)
  CALL OC_Basis_Initialise(basis,err)
  CALL OC_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  IF(numberOfGlobalZElements==0) THEN
    !Set the basis to be a bilinear Lagrange basis
    !CALL OC_Basis_TypeSet(basis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
    CALL OC_Basis_NumberOfXiSet(basis,2,err)
    !CALL OC_Basis_InterpolationXiSet(basis,[3,3],err)
    !CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[4,4],err) 
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL OC_Basis_NumberOfXiSet(basis,3,err)
  ENDIF
  !Finish the creation of the basis
  CALL OC_Basis_CreateFinish(basis,err)

  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_Initialise(generatedMesh,err)
  CALL OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL OC_GeneratedMesh_BasisSet(generatedMesh,basis,err)   
  !Define the mesh on the region
  IF(numberOfGlobalZElements==0) THEN
    CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
  ELSE
    CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
      & numberOfGlobalZElements],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL OC_Mesh_Initialise(mesh,err)
  CALL OC_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! MESH DECOMPOSITION
  !----------------------------------------------------------------------------------------------------------- 
  
  !Create a decomposition
  CALL OC_Decomposition_Initialise(decomposition,err)
  CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  !Finish the decomposition
  CALL OC_Decomposition_CreateFinish(decomposition,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSER
  !-----------------------------------------------------------------------------------------------------------
  
  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------  
  
  !Start to create a default (geometric) field on the region
  CALL OC_Field_Initialise(geometricField,err)
  CALL OC_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  !Set the decomposition to use
  CALL OC_Field_DecompositionSet(geometricField,decomposition,err)
  !Set the domain to be used by the field components.
  CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,1,1,err)
  CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,2,1,err)
  IF(numberOfGlobalZElements/=0) THEN
    CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,3,1,err)
  ENDIF
  !Finish creating the field
  CALL OC_Field_CreateFinish(geometricField,err)
  
  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------  
 
  !Create the equations set and set the equations set to be a generalised diffusion problem

  CALL OC_EquationsSet_Initialise(equationsSet,err)
  CALL OC_Field_Initialise(equationsSetField,err)
  CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & OC_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,OC_EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE], &
    & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
  
  !Finish creating the equations set
  CALL OC_EquationsSet_CreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set dependent field variables
  CALL OC_Field_Initialise(dependentField,err)
  CALL OC_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  !Finish the equations set dependent field variables
  CALL OC_EquationsSet_DependentCreateFinish(equationsSet,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! MATERIAL FIELD
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the equations set material field variables
  CALL OC_Field_Initialise(materialsField,err)
  CALL OC_EquationsSet_MaterialsCreateStart(equationsSet,MATERIALS_FIELD_USER_NUMBER,materialsField,err)
  !Finish the equations set dependent field variables
  CALL OC_EquationsSet_MaterialsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTIC FIELD
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the equations set analytic field variables
  CALL OC_Field_Initialise(analyticField,err)
  IF(numberOfGlobalZElements==0) THEN  
    CALL OC_EquationsSet_AnalyticCreateStart(equationsSet,OC_EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1, &
      & ANALYTIC_FIELD_USER_NUMBER,analyticField,err)
  ELSE
    WRITE(*,'(A)') "Three dimensions is not implemented."
    STOP
  ENDIF
  
  !Finish the equations set analytic field variables
  CALL OC_EquationsSet_AnalyticCreateFinish(equationsSet,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------  

  !Create the equations set equations
  CALL OC_Equations_Initialise(equations,err)
  CALL OC_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  !Set the equations matrices sparsity type
  CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output
  !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT,err)
  !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_TIMING_OUTPUT,err)
  !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_MATRIX_OUTPUT,err)
  !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL OC_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem
  CALL OC_Problem_Initialise(problem,err)
  CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[OC_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & OC_PROBLEM_DIFFUSION_EQUATION_TYPE,OC_PROBLEM_LINEAR_DIFFUSION_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL OC_Problem_CreateFinish(problem,err)

  !Create the problem control
  CALL OC_Problem_ControlLoopCreateStart(problem,err)
  !Get the control loop
  CALL OC_ControlLoop_Initialise(controlLoop,err)
  CALL OC_Problem_ControlLoopGet(problem,OC_CONTROL_LOOP_NODE,controlLoop,err)
  !Set the times
  CALL OC_ControlLoop_TimesSet(controlLoop,0.0_OC_RP,1.001_OC_RP,0.001_OC_RP,err)
  !Finish creating the problem control loop
  CALL OC_Problem_ControlLoopCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER
  !-----------------------------------------------------------------------------------------------------------
  
  !Start the creation of the problem solvers

  CALL OC_Solver_Initialise(solver,err)
  CALL OC_Solver_Initialise(linearSolver,err)
  CALL OC_Problem_SolversCreateStart(problem,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  !CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_NO_OUTPUT,err)
  CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_PROGRESS_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_TIMING_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_SOLVER_OUTPUT,err)
  CALL OC_Solver_DynamicLinearSolverGet(solver,linearSolver,err)
  CALL OC_Solver_LinearIterativeMaximumIterationsSet(linearSolver,1000_OC_Intg,err)
  !Finish the creation of the problem solver
  CALL OC_Problem_SolversCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the problem solver equations
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_SolverEquations_Initialise(solverEquations,err)
  CALL OC_Problem_SolverEquationsCreateStart(problem,err)
  !Get the solve equations
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  CALL OC_Solver_SolverEquationsGet(solver,solverEquations,err)
  !Set the solver equations sparsity
  CALL OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_SPARSE_MATRICES,err)
  !CALL OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_FULL_MATRICES,err)  
  !Add in the equations set
  CALL OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL OC_Problem_SolverEquationsCreateFinish(problem,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the solver equations boundary conditions
  CALL OC_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsAnalytic(solverEquations,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! SOLVE
  !-----------------------------------------------------------------------------------------------------------
  
  !Solve the problem
  CALL OC_Problem_Solve(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! OUTPUT
  !----------------------------------------------------------------------------------------------------------

  !Output Analytic analysis
  Call OC_AnalyticAnalysis_Output(dependentField,"diffusion_equation_analytic",err)

  exportField=.TRUE.
  IF(exportField) THEN
    CALL OC_Fields_Initialise(fields,err)
    CALL OC_Fields_Create(region,fields,err)
    CALL OC_Fields_NodesExport(fields,"diffusion_equation","FORTRAN",err)
    CALL OC_Fields_ElementsExport(fields,"diffusion_equation","FORTRAN",err)
    CALL OC_Fields_Finalise(fields,err)
  ENDIF
  
  !Destroy the context
  CALL OC_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL OC_Finalise(err)
  
  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM DiffusionEquation
