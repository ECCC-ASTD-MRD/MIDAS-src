# Coding Standards

The purpose of this is to summarize coding standards to be used in the
MIDAS fortran code. A shorter list of the "Top 10" things to keep in
mind is provided [here](codingStd_top10.md).

## Examples 

- [random perturbation main program](./examples/program_template.f90)
- [vertical coordinate module](./examples/module_template.f90)

## A long list of coding standard rules (based on the [jules land surface model coding standards](http://jules-lsm.github.io/coding_standards/) with modifications)

### General formatting and layout:

- Use the Fortran 90 free format syntax.
- Indent blocks by 2 characters.
- Use spaces and blank lines where appropriate to format your code to
  improve readability (use genuine spaces rather than tabs, as the tab
  character is not in the Fortran character set).
- Try to confine your line width to 80 characters. This means that
  your code can be viewed easily in any editor on any screen, and can
  be printed easily on standard paper.
- Line up your statements, where appropriate, to improve readability.
- Short and simple Fortran statements are easier to read and
  understand than long and complex ones. Where possible, avoid using
  continuation lines in a statement.
- Avoid putting multiple statements on the same line. It is not good
  for readability.
- Each program unit (module, subroutine, function etc.) should follow
  a structure similar to the examples supplied above. The intended
  behaviour of the unit should be clearly described in the header, and
  the prefix of the module should be defined.
- Each subroutine, function and module should be in a separate file.
  Modules may be used to group related variables, subroutines and
  functions. At this stage, all code except "main" programs should
  be inside modules.
- When naming your variables and program units, always keep in mind
  that Fortran is a case-insensitive language (e.g. `EditOrExit` is the
  same as `EditorExit`.)
- Use only characters in the Fortran character set. In particular,
  accent characters and tabs are not allowed in code, although they
  are usually OK in comments. If your editor inserts tabs
  automatically, you should configure it to switch off the
  functionality when you are editing Fortran source files.
- Although Fortran has no reserved keywords, you should avoid naming
  your program units and variables with names that match an intrinsic
  `FUNCTION` or `SUBROUTINE`. Similarly, you should avoid naming your
  program units and variables with names that match a keyword in a
  Fortran statement.
- To improve readability, write your code using the lower case for all
  Fortran keywords and intrinsic functions/subroutines. The rest of
  the code is written in either lowercase with underscores or,
  preferably using "camelCase" (e.g. `variableNameList` or
  `gsv_getField`). When naming any public entity (variable, subroutine,
  or function) in a module it must begin with the module prefix: e.g.
  `gsv_getField`, where "`gsv`" is the prefix associated with the module
  `gridStateVector_mod`.
- Function or subroutine arguments should be declared separately from,
  and before, local variables, separated by a blank line.

### Commenting:

- Comments should start with a single `!` and be indented with the code.
  A blank line should be left before (but not after) the comment line.
- An important exception to the above rule on comments is for any
  comments that are intended to appear in the automatically generated
  on-line documentation. Refer to the [documentation page](documentationStandards.md) for more details.
- Be generous with comments. State the reason for doing something,
  instead of repeating the Fortran logic in words. However, comments
  are not a replacement for clear code. Better that the code itself is
  easy to understand by the way it is structured and the choice of
  variable and subroutine/function names - it is ok (and encouraged)
  to use longer names that really describe what the variable or
  subroutine/function represents or does (also easier later to
  grep/search for these names as compared with a variable named i, for
  example).

### Provided Functionalities to Use

- Use [`msg_message()`](src/modules/message_mod.f90) instead of naked
  `write(*,*)` to output information: provide the _origin_ of the message (such
  as the caller subroutine, function or program) and a _verbosity level_ from 0
  to 3 that specifies how important is the message:

  * 0 : critical, always printed
  * 1 : default priority; printed in operational context
  * 2 : detailed ouptut, provides extra information
  * 3 : intended for developpers, printed for debugging or specific diagnostcs

### More detailed rules:

- Use the new and clearer syntax for `LOGICAL` comparisons, i.e.:

```fortran
== instead of .EQ.
/= instead of .NE.
&gt; instead of .GT.
&lt; instead of .LT.
&gt;= instead of .GE.
&lt;= instead of .LE.
```

- Positive logic is usually easier to understand. When using an
  `IF-ELSE-END IF` construct you should use positive logic in the IF
  test, provided that the positive and the negative blocks are about
  the same length. It may be more appropriate to use negative logic if
  the negative block is significantly longer than the positive block.
- To improve readability, you should always use the optional space to
  separate the following Fortran keywords:

```fortran
else if
end do
end forall
end function
end if
end interface
end module
end program
end select
end subroutine
end type
end where
select case
```

- If you have a large or complex code block embedding other code
  blocks, you may consider naming some or all of them to improve
  readability.
- Improve readability by always using the full version of the `END`
  statement (i.e. `END SUBROUTINE <name>` or `END FUNCTION <name>`
  instead of just `END`) at the end of each sub-program unit.
- Where possible, consider using `CYCLE`, `EXIT` or a `WHERE`-construct to
  simplify complicated `DO`-loops.
- When writing a `REAL` literal with an integer value, put a `0` after the
  decimal point (i.e. `1.0` as opposed to `1.0`) to improve readability.
  For double precision real literals (`real(8)`) always include "`d0`"
  as in `1.0d0` to ensure the correct precision.
- Where reasonable and sensible to do so, you should try to match the
  names of dummy and actual arguments to a `SUBROUTINE`/`FUNCTION`.
- In an array assignment, it is recommended that you use array
  notations to improve readability, e.g.  
  Avoid this:
  ```fortran
  array1 = 1
  array2 = array1 * scalar
  ```
  Use this instead:
  ```fortran
  array1(:,&nbsp;:) = 1
  array2(:,&nbsp;:) = array1(:,&nbsp;:) * scalar
  ```

- Use `IMPLICIT NONE` in all program units. This forces you to declare
  all your variables explicitly. This helps to reduce bugs in your
  program that will otherwise be difficult to track.
- Design any derived data types carefully and use them to group
  related variables. Appropriate use of derived data types will allow
  you to design modules and procedures with simpler and cleaner
  interfaces.
- Always use a `private` statement at the beginning of a module so
  that all module variables and procedures are by default declared
  private. Any public subroutines, functions or variables are then
  explicity declare public using a `public` statement near the
  beginning of the module to make it clear to a user what is
  accessible from the outside world.
- Where possible, an `ALLOCATE` statement for an `ALLOCATABLE` array (or a
  `POINTER` used as a dynamic array) should be coupled with a `DEALLOCATE`
  within the same scope. If an `ALLOCATABLE` array is a PUBLIC
  MODULE variable, it is highly desirable if its memory allocation and
  deallocation are only performed in procedures within the MODULE in
  which it is declared. You may consider writing specific SUBROUTINES
  within the MODULE to handle these memory managements.
- To avoid memory fragmentation, it is desirable to `DEALLOCATE` in
  reverse order of `ALLOCATE`, as in:

  ```fortran
  ALLOCATE(a(n))
  ALLOCATE(b(n))
  ALLOCATE(c(n))

  ! ... do something ...

  DEALLOCATE(c)
  DEALLOCATE(b)
  DEALLOCATE(a)
  ```

- Inside a function or subroutine, always define a local `POINTER` before using it. 
  DO NOT define a `POINTER` in its declaration by pointing it to the intrinsic
  function `NULL()` since this will invoke an implicit "`save`" attribute which is very
  dangerous! 
  Instead, make sure that your POINTER is defined or
  nullified early on in the program unit. Similarly, NULLIFY a POINTER
  when it is no longer in use, either by using the NULLIFY statement
  or by pointing your `POINTER` to `NULL()`.
  This recommandation does not apply to `POINTER` global to a module
  or program.
  (The reason for this, is that the declaration statement is actually
  executed **only** on the first call, the only one for module or program
  pointer, but for functions or subroutines, subsequent calls will not
  re-declare the pointer to `null()` such that it might already points
  to some value obtained in the previous call.)
- Avoid the `DIMENSION` attribute or statement. Declare the DIMENSION
  with the declared variables. E.g.:

  Avoid this:
  ```fortran
  INTEGER, DIMENSION(10)&nbsp;:: array1
  INTEGER&nbsp;:: array2
  DIMENSION&nbsp;:: array2(20)
  ```
  Instead, do this:
  ```fortran
  INTEGER&nbsp;:: array1(10), array2(20)
  ```

- Avoid `COMMON` blocks and `BLOCK DATA` program units. Instead, use a
  `MODULE` with `PUBLIC` variables.
- Avoid the `EQUIVALENCE` statement. Use a `POINTER` or a derived data
  type, and the `TRANSFER` intrinsic function to convert between types.
- Avoid the `PAUSE` statement, as your program will hang in a batch
  environment. If you need to halt your program for interactive use,
  consider using a `READ*` statement instead.
- Avoid the `ENTRY` statement. Use a `MODULE` or internal `SUBROUTINE`.
- Avoid the `GOTO` statement.
- Avoid numbered statement labels. `DO ... label CONTINUE` constructs
  should be replaced by `DO ... END DO` constructs. Every `DO` loop must
  be terminated with a corresponding `END DO`.
- Never use a `FORMAT` statement - they require the use of labels, and
  obscure the meaning of the I/O statement. The formatting information
  can be placed explicitly within the `READ`, `WRITE` or `PRINT` statement,
  or be assigned to a `CHARACTER` variable in a `PARAMETER` statement in
  the header of the routine for later use in I/O statements. Never
  place output text within the format specifier, i.e. only format
  information may be placed within the `FMT=` part of an I/O statement.
  All variables and literals, including any character literals, must
  be 'arguments' of the I/O routine itself. This improves readability
  by clearly separating what is to be read/written from how to
  read/write it.
- Avoid the `FORALL` statement/construct. Despite what it is supposed to
  do, `FORALL` is often difficult for compilers to optimise (see, for
  example, Implementing the Standards including Fortran 2003 by NAG).
  Stick to the equivalent `DO` construct, `WHERE` statement/construct or
  array assignments unless there are actual performance benefits from
  using `FORALL`.
- A `FUNCTION` should be `PURE`, i.e. it should have no side effects (e.g.
  altering an argument or module variable, or performing I/O). If you
  need to perform a task with side effects, you should use a
  `SUBROUTINE` instead.
- Declare the `INTENT` of all arguments to a subroutine or function.
  This allows checks against unintended access of variables to be done
  at compile time. The above point requiring functions to be pure
  means that all arguments of a `FUNCTION` should be declared as
  `INTENT(IN)`.
- Avoid `RECURSIVE` procedures if possible. `RECURSIVE` procedures are
  usually difficult to understand, and are always difficult to
  optimise in a supercomputer environment.
- Avoid using the specific names of intrinsic procedures. Use the
  generic names of intrinsic procedures where possible.
- Not necessary to use the `ONLY` clause in a `USE <module>` statement,
  since each module should have few public symbols and all should be
  named beginning with the module prefix.
- The use of operator overloading is discouraged, as it can lead to
  confusion. The only acceptable use is to allow the standard
  operators (+, - etc.) to work with derived data types, where this
  makes sense.
- Avoid using archaic Fortran 77 features and features deprecated in
  Fortran 90.
