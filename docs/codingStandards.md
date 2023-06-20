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
  `function` or `subroutine`. Similarly, you should avoid naming your
  program units and variables with names that match a keyword in a
  Fortran statement.
- To improve readability, write your code using the lower case for all
  Fortran keywords and intrinsic functions/subroutines. The rest of
  the code should be written using "camelCase" with no underscores, except
  between a module prefix and a variable name
  (e.g. `variableNameList` or `gsv_getField`). When naming any public entity
  (variable, subroutine, or function) in a module it must begin with
  the module prefix: e.g. `gsv_getField`, where "`gsv`" is the prefix
  associated with the module `gridStateVector_mod`.
- Function or subroutine arguments should be declared separately from,
  and before, local variables, separated by a blank line.

### Programs and modules naming convention
- As for variables, "camelCase" (https://en.wikipedia.org/wiki/Camel_case)
  should be used for all filenames. For modules, an underscore should
  appear only between the module name and the suffix "`mod`. Therefore, 
  a module filename should end with "`_mod.f90`".

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

- Use [`msg()`](src/modules/message_mod.f90) instead of naked
  `write(*,*)` to output information: provide the _origin_ of the message (such
  as the caller subroutine, function or program).
  ```fortran
  call msg('int_tInterp_gsv', 'START', verb_opt=2)
  ! prints a short message on all MPI tiles when verbosity threshold >= 2
  ...
  call msg('int_tInterp_gsv', 'numStepIn='//str(numStepIn)&
       //',numStepOut='//str(numStepOut), mpiAll_opt=.false.)
  ! prints a short message with some numerical values on MPI tile 0 only
  ```
  Optionally, a _verbosity level_ that specifies how important is the message
  can be provided.  The verbosity thresholds are defined as follow:

    * `msg_ALWAYS` : always printed, irrespectively of threshold
    * 0 : critical, should always be printed
    * 1 : default priority, printed in operational context
    * 2 : detailed output, provides extra information
    * 3 : intended for developers, printed for debugging or specific diagnostcs
    * `msg_NEVER` : never printed, irrespectively of threshold

### More detailed rules:

- Use the new and clearer syntax for `logical` comparisons, i.e.:
  ```fortran
  == instead of .eq.
  /= instead of .ne.
  > instead of .gt.
  < instead of .lt.
  >= instead of .ge.
  <= instead of .le.
  ```

- Positive logic is usually easier to understand. When using an
  `if-else-end if` construct you should use positive logic in the IF
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
- Improve readability by always using the full version of the `end`
  statement (i.e. `end subroutine <name>` or `end function <name>`
  instead of just `end`) at the end of each sub-program unit.
- Where possible, consider using `cycle`, `exit` or a `where`-construct to
  simplify complicated `DO`-loops.
- When writing a `real` literal with an integer value, put a `0` after the
  decimal point (i.e. `1.0` as opposed to `1`) to improve readability.
  For double precision real literals (`real(8)`) always include "`d0`"
  as in `1.0d0` to ensure the correct precision.
- Where reasonable and sensible to do so, you should try to match the
  names of dummy and actual arguments to a `subroutine`/`function`.
- In an array assignment, it is recommended that you use array
  notations to improve readability, e.g.  
  Avoid this:
  ```fortran
  array1 = 1
  array2 = array1 * scalar
  ```
  Use this instead:
  ```fortran
  array1(:,:) = 1
  array2(:,:) = array1(:,:) * scalar
  ```

- Use `implicit none` in all program units. This forces you to declare
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
- Where possible, an `allocate` statement for an `allocatable` array (or a
  `pointer` used as a dynamic array) should be coupled with a `deallocate`
  within the same scope. if an `allocatable` array is a public
  `module` variable, it is highly desirable if its memory allocation and
  deallocation are only performed in procedures within the `module` in
  which it is declared. You may consider writing specific `subroutines`
  within the module to handle these memory managements.
- To avoid memory fragmentation, it is desirable to `deallocate` in
  reverse order of `allocate`, as in:
  ```fortran
  allocate(a(n))
  allocate(b(n))
  allocate(c(n))
  ! ... do something ...
  deallocate(c)
  deallocate(b)
  deallocate(a)
  ```

- Inside a function or subroutine, always define a local `pointer` before using it. 
  do not define a `pointer` in its declaration by pointing it to the intrinsic
  function `null()` since this will invoke an implicit "`save`" attribute which is very
  dangerous! 
  Instead, make sure that your `pointer` is defined or
  nullified early on in the program unit. similarly, `nullify` a `pointer`
  when it is no longer in use, either by using the `nullify` statement
  or by pointing your `pointer` to `null()`.
  This recommandation does not apply to `pointer` global to a module
  or program.
  (The reason for this, is that the declaration statement is actually
  executed **only** on the first call, the only one for module or program
  pointer, but for functions or subroutines, subsequent calls will not
  re-declare the pointer to `null()` such that it might already points
  to some value obtained in the previous call.)
- Avoid the `dimension` attribute or statement. declare the `dimension`
  with the declared variables. E.g.:
  Avoid this:
  ```fortran
  integer, dimension(10) :: array1
  integer :: array2
  dimension :: array2(20)
  ```
  Instead, do this:
  ```fortran
  integer :: array1(10), array2(20)
  ```

- Never initialize a local variable on the declaration unless the `save`
  attribute is explicitely present.  
  Avoid this:
  ```fortran
  logical :: trueByDefault = .true.
  ```
  Instead, do this:
  ```fortran
  logical :: trueByDefault
  trueByDefault = .true.
  ```
  If you actually want the local variable to keep it's value after passing out
  of scope, be explicit:
  ```fortran
  logical, save :: thisIsTheFirstCall = .true.
  ```

- Avoid `common` blocks and `block data` program units. instead, use a
  `module` with `public` variables.
- Avoid the `equivalence` statement. Use a `pointer` or a derived data
  type, and the `transfer` intrinsic function to convert between types.
- Avoid the `pause` statement, as your program will hang in a batch
  environment. If you need to halt your program for interactive use,
  consider using a `read*` statement instead.
- Avoid the `entry` statement. Use a `module` or internal `subroutine`.
- Avoid the `goto` statement.
- Avoid numbered statement labels. `do ... label continue` constructs
  should be replaced by `do ... end do` constructs. every `do` loop must
  be terminated with a corresponding `end do`.
- Never use a `format` statement - they require the use of labels, and
  obscure the meaning of the I/O statement. The formatting information
  can be placed explicitly within the `read`, `write` or `print` statement,
  or be assigned to a `character` variable in a `parameter` statement in
  the header of the routine for later use in I/O statements. Never
  place output text within the format specifier, i.e. only format
  information may be placed within the `fmt=` part of an I/O statement.
  All variables and literals, including any character literals, must
  be 'arguments' of the I/O routine itself. This improves readability
  by clearly separating what is to be read/written from how to
  read/write it.
- Avoid the `forall` statement/construct. Despite what it is supposed to
  do, `forall` is often difficult for compilers to optimise (see, for
  example, Implementing the Standards including Fortran 2003 by NAG).
  Stick to the equivalent `do` construct, `where` statement/construct or
  array assignments unless there are actual performance benefits from
  using `forall`.
- A `function` should be `pure`, i.e. it should have no side effects (e.g.
  altering an argument or module variable, or performing I/O). If you
  need to perform a task with side effects, you should use a
  `subroutine` instead.
- Declare the `intent` of all arguments to a subroutine or function.
  This allows checks against unintended access of variables to be done
  at compile time. The above point requiring functions to be pure
  means that all arguments of a `function` should be declared as
  `intent(in)`.
- Avoid `recursive` procedures if possible. `recursive` procedures are
  usually difficult to understand, and are always difficult to
  optimise in a supercomputer environment.
- Avoid using the specific names of intrinsic procedures. Use the
  generic names of intrinsic procedures where possible.
- Not necessary to use the `only` clause in a `use <module>` statement,
  since each module should have few public symbols and all should be
  named beginning with the module prefix.
- The use of operator overloading is discouraged, as it can lead to
  confusion. The only acceptable use is to allow the standard
  operators (+, - etc.) to work with derived data types, where this
  makes sense.
- Avoid using archaic Fortran 77 features and features deprecated in
  Fortran 90.
