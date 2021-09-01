/*
programmers can use ASSERT macro to add checks intended to be performed in non release builds.
The purpose of ASSERT is to detect programmer bugs, not bad use of a program: user errors are
supposed to be handled gracefully with a self explanatory error message targetted at the user.
*/
#ifdef  DISABLE_HIB_ASSERT
#define ASSERT(EX)
#define REQUIRE(EX)
#define ENSURE(EX)
#else
#define   ASSERT(EX) if(.not. ( EX )) call fassert(__FILE__, __LINE__)
/*
#define   REQUIRE(EX) if(.not. ( EX )) call FortranAssert('Precondition',#EX, __FILE__, __LINE__)
#define   ENSURE(EX) if(.not. ( EX )) call FortranAssert('Postcondition',#EX, __FILE__, __LINE__)
*/
#endif
