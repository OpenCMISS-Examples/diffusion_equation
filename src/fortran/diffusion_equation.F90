PROGRAM DiffusionEquation

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------
  
  !Test program parameters
  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=3.0_CMISSRP
  
  INTEGER(CMISSIntg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: REGION_USER_NUMBER=3
  INTEGER(CMISSIntg), PARAMETER :: BASIS_USER_NUMBER=4
  INTEGER(CMISSIntg), PARAMETER :: GENERATED_MESH_USER_NUMBER=5
  INTEGER(CMISSIntg), PARAMETER :: MESH_USER_NUMBER=6
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSITION_USER_NUMBER=7
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSER_USER_NUMBER=8
  INTEGER(CMISSIntg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=9
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=10
  INTEGER(CMISSIntg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=11
  INTEGER(CMISSIntg), PARAMETER :: MATERIALS_FIELD_USER_NUMBER=12
  INTEGER(CMISSIntg), PARAMETER :: ANALYTIC_FIELD_USER_NUMBER=13
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=14
  INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER=15

  !Program types
  
  !Program variables


  !INTEGER(INTG) :: first_global_dof,first_local_dof,first_local_rank,last_global_dof,last_local_dof,last_local_rank,rank_idx
  !INTEGER(INTG) :: EQUATIONS_SET_INDEX !TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOF_MAPPING
  
  !CMISS variables
  TYPE(cmfe_BasisType) :: basis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
  TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_DecomposerType) :: decomposer
  TYPE(cmfe_EquationsType) :: equations
  TYPE(cmfe_EquationsSetType) :: equationsSet
  TYPE(cmfe_FieldType) :: geometricField,equationsSetField,dependentField,materialsField,analyticField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh  
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_ControlLoopType) :: controlLoop
  TYPE(cmfe_RegionType) :: region,worldRegion
  TYPE(cmfe_SolverType) :: solver,linearSolver
  TYPE(cmfe_SolverEquationsType) :: solverEquations
  TYPE(cmfe_WorkGroupType) :: worldWorkGroup

  INTEGER(CMISSIntg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements
  INTEGER(CMISSIntg) :: computationalNodeNumber,numberOfComputationalNodes
  INTEGER(CMISSIntg) :: decompositionIndex,equationsSetIndex
  INTEGER(CMISSIntg) :: err
  LOGICAL :: exportField

  STOP

  !Intialise OpenCMISS
  CALL cmfe_Initialise(err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  !Create a context
  CALL cmfe_Context_Initialise(context,err)
  CALL cmfe_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)

  !Get the computational nodes information
  CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
  CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  numberOfGlobalXElements=10
  numberOfGlobalYElements=10
  numberOfGlobalZElements=0

  !-----------------------------------------------------------------------------------------------------------
  ! COORDINATE SYSTEM
  !----------------------------------------------------------------------------------------------------------- 
  
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  IF(numberOfGlobalZElements==0) THEN
    !Set the coordinate system to be 2D
    CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,2,err)
  ELSE
    !Set the coordinate system to be 3D
    CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,3,err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------
 
  !Start the creation of the region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL cmfe_Region_LabelSet(region,"diffusion_equation",err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(region,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! BASIS
  !----------------------------------------------------------------------------------------------------------  
  
  !Start the creation of a basis (default is trilinear lagrange)
  CALL cmfe_Basis_Initialise(basis,err)
  CALL cmfe_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  IF(numberOfGlobalZElements==0) THEN
    !Set the basis to be a bilinear Lagrange basis
    !CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
    CALL cmfe_Basis_NumberOfXiSet(basis,2,err)
    !CALL cmfe_Basis_InterpolationXiSet(basis,[3,3],err)
    !CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[4,4],err) 
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL cmfe_Basis_NumberOfXiSet(basis,3,err)
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(basis,err)

  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,basis,err)   
  !Define the mesh on the region
  IF(numberOfGlobalZElements==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
      & numberOfGlobalZElements],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(mesh,err)
  CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! MESH DECOMPOSITION
  !----------------------------------------------------------------------------------------------------------- 
  
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSER
  !-----------------------------------------------------------------------------------------------------------
  
  CALL cmfe_Decomposer_Initialise(decomposer,err)
  CALL cmfe_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL cmfe_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL cmfe_Decomposer_CreateFinish(decomposer,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------  
  
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(geometricField,err)
  CALL cmfe_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  !Set the decomposition to use
  CALL cmfe_Field_DecompositionSet(geometricField,decomposition,err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
  CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,err)
  IF(numberOfGlobalZElements/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,err)
  ENDIF
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(geometricField,err)
  
  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------  
 
  !Create the equations set and set the equations set to be a generalised diffusion problem

  CALL cmfe_EquationsSet_Initialise(equationsSet,err)
  CALL cmfe_Field_Initialise(equationsSetField,err)
  CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE], &
    & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
  
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(dependentField,err)
  CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! MATERIAL FIELD
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the equations set material field variables
  CALL cmfe_Field_Initialise(materialsField,err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSet,MATERIALS_FIELD_USER_NUMBER,materialsField,err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTIC FIELD
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the equations set analytic field variables
  CALL cmfe_Field_Initialise(analyticField,err)
  IF(numberOfGlobalZElements==0) THEN  
    CALL cmfe_EquationsSet_AnalyticCreateStart(equationsSet,CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1, &
      & ANALYTIC_FIELD_USER_NUMBER,analyticField,err)
  ELSE
    WRITE(*,'(A)') "Three dimensions is not implemented."
    STOP
  ENDIF
  
  !Finish the equations set analytic field variables
  CALL cmfe_EquationsSet_AnalyticCreateFinish(equationsSet,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------  

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_TIMING_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_MATRIX_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem
  CALL cmfe_Problem_Initialise(problem,err)
  CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & CMFE_PROBLEM_DIFFUSION_EQUATION_TYPE,CMFE_PROBLEM_LINEAR_DIFFUSION_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(problem,err)

  !Create the problem control
  CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
  !Get the control loop
  CALL cmfe_ControlLoop_Initialise(controlLoop,err)
  CALL cmfe_Problem_ControlLoopGet(problem,CMFE_CONTROL_LOOP_NODE,controlLoop,err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(controlLoop,0.0_CMISSRP,1.001_CMISSRP,0.001_CMISSRP,err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER
  !-----------------------------------------------------------------------------------------------------------
  
  !Start the creation of the problem solvers

  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Solver_Initialise(linearSolver,err)
  CALL cmfe_Problem_SolversCreateStart(problem,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_NO_OUTPUT,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_TIMING_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_SOLVER_OUTPUT,err)
  CALL cmfe_Solver_DynamicLinearSolverGet(solver,linearSolver,err)
  CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(linearSolver,1000_CMISSIntg,err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_SolverEquations_Initialise(solverEquations,err)
  CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)
  !Get the solve equations
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_SPARSE_MATRICES,err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_FULL_MATRICES,err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the solver equations boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(solverEquations,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! SOLVE
  !-----------------------------------------------------------------------------------------------------------
  
  !Solve the problem
  CALL cmfe_Problem_Solve(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! OUTPUT
  !----------------------------------------------------------------------------------------------------------

  !Output Analytic analysis
  Call cmfe_AnalyticAnalysis_Output(dependentField,"diffusion_equation_analytic",err)

  exportField=.TRUE.
  IF(exportField) THEN
    CALL cmfe_Fields_Initialise(fields,err)
    CALL cmfe_Fields_Create(region,fields,err)
    CALL cmfe_Fields_NodesExport(fields,"diffusion_equation","FORTRAN",err)
    CALL cmfe_Fields_ElementsExport(fields,"diffusion_equation","FORTRAN",err)
    CALL cmfe_Fields_Finalise(fields,err)
  ENDIF
  
  !Destroy the context
  CALL cmfe_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL cmfe_Finalise(err)
  
  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM DiffusionEquation
