SET(ALGORITHM_PROFILE_UID_Namespace "org.mitk")
SET(ALGORITHM_PROFILE_UID_Name "RigidICP.3D.default")
SET(ALGORITHM_PROFILE_UID_Version "1.0.0")

SET(ALGORITHM_PROFILE_Description "Simple 3D rigid (translation and euler angles) registration algorithm using point sets and the iterative closesed points scheme. It tries to minimize the minimal point distance errors (no point paires are assumed). Remark: at least 6 points per point sets are needed; the number must not equal.")
SET(ALGORITHM_PROFILE_Contact "Ralf Floca\; mitk-users@lists.sourceforge.net")
SET(ALGORITHM_PROFILE_Citation "P.J. Besl, N.D. McKay.: A Method for Registration of 3-D Shapes. IEEE Trans. on Pattern Analysis and Machine Intelligence.  Vol. 14 2. pp. 239-256. 1992. doi:10.1109/34.121791.")

SET(ALGORITHM_PROFILE_DataType "Points")
SET(ALGORITHM_PROFILE_ResolutionStyle "Single")
SET(ALGORITHM_PROFILE_DimMoving "3")
SET(ALGORITHM_PROFILE_ModalityMoving "any")
SET(ALGORITHM_PROFILE_DimTarget "3")
SET(ALGORITHM_PROFILE_ModalityTarget "any")
SET(ALGORITHM_PROFILE_Subject "any")
SET(ALGORITHM_PROFILE_Object "any")
SET(ALGORITHM_PROFILE_TransformModel "rigid")
SET(ALGORITHM_PROFILE_TransformDomain "global")
SET(ALGORITHM_PROFILE_ComputationStyle "iterative")
SET(ALGORITHM_PROFILE_Deterministic "1")
SET(ALGORITHM_PROFILE_Keywords "basic" "point sets" "ICP" "rigid")