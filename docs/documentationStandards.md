# Documentation Standards

The purpose of this page is to summarize how to add to MIDAS code comments that are destined for the Sphinx automatic documentation.

Sphinx reads and interprets header comments (comments directly after declaring a program, module, subroutine, etc.). Because Sphinx will format these comments as does a word processor, don't try to be fancy. Because Sphinx responds to reStructuredText mark-up in the comments, don't use any of the mark-up flags, unless you understand the mark-up language and want to direct Sphinx in the presentation of those comments. 

## 1- Recommended Format 
```fortran
program myprog
  !
  !:Purpose: Describe why the program was created.  Briefly describe what the
  !          program does.

  use myMod_mod
  implicit none

  call mym_mySub()
end program myprog
```
```fortran
module myMod_mod
  ! MODULE myMod_mod (prefix='mym' category='4. Observation operators')
  !
  !:Purpose: Describe why the module was created.  Briefly describe what the
  !          module can be used to do.
  implicit none
  private

  public mym_mySub

  real    :: module_a ! Description of namelist variable(each on a separate line)
  real    :: module_b ! Not a namelist variable.
  integer :: module_c ! This is another namelist variable

  subroutine mym_mySub(a,b,c,d)
    !
    !:Purpose: This is a basic description of the method (subroutine or
    !          function).  Do not use a blank line either before or inside the
    !          purpose block.  The purpose should give a good idea of the
    !          method's behaviour (what it does) and not delve into the
    !          implementation (how it does it). It might be useful to use
    !          bullets.  Bullets are preceeded and followed by a blank comment. 
    !          Blank comments between bullets are optional and will not be
    !          rendered.
    !
    !            - item 1
    !            - item 2
    !            - item 3
    !
    !          Do not give any revision information here. This is handled by GIT.
    !
    !          A literal block might be useful ::
    !                   This is particularly useful for equations
    !               or ASCII art.  New-lines are preserved.
    !               It ends when text returns to the preceding paragraph's indentation
    !          This is not in the literal block.
    !
    !:Arguments:
    !        :a:     Description of A -- If a long description is necessary, it
    !                can be written here.  This is a very, very very, long
    !                description of A.  This is, in fact, an extreeeeemely long
    !                description of A.
    implicit none

    ! Arguments:  (This line is not part of the Sphinx documentation because the
    !             previous line is a blank.  -- Each argument should be on a
    !             separate line.)
    integer, intent(in) :: a ! Description of A (preferred location)
    integer, intent(in) :: b ! Description of B
    real,    intent(out):: c ! Description of C
    real,    intent(out):: d ! Description of D

    ! Locals:
    integer :: e

    namelist /my_namelist/ module_a, module_c

    <your code>
  end subroutine mym_mySub
end module myMod_mod
```


## 2- Things to Avoid within Sphinx-Interpreted Regions

A Sphinx-interpreted region is the comment section that immediately follows a method declaration, or the comment on the same line as an argument declaration. 

* Do not use this syntax: it causes Sphinx to hang:
  ```fortran
  real :: myVar(A+B:C+D)
  ```
  where

    * `A, B, C, D` are arguments
    * both additions are necessary to cause a hang
    * It could be any type (`real`, `integer`, `complex`, other) in the example above

* Do not use any of these character strings. They all have special meaning to Sphinx. Unless, of course, you understand the meaning and want to communicate with Sphinx for special formatting: 
    * `*` (as first non-blank character) 
    * `-` (as first non-blank character) 
    * `+` (as first non-blank character) 
    * `#` (as first non-blank character) 
    * `**` 
    * `++` 
    * `--` 
    * `__` (two underscores) 
    *  `::` 
    * `..` 
    * trailing underscore; e.g. `myVar_` 
* Do not indent your comments; start all comments in the same column.
* Do not split the comments for subroutine variable declaration between multiple lines. Sphinx process only the first line.
* Do not give Sphinx a red herring. Consider the following code: 
  ```fortran
  subroutine hbht_compute(columng, columnhr)
  integer :: columng   ! These are my
  integer :: columnhr  ! real comments
  call hbht_compute_static(columng      & ! RED
                             columnhr     & ! HERRING
                             )
  ```
  Comments like this are allowed but, because of the line continuations, Sphinx does not verify that the comments are associated with the subroutine `hbht_compute`. 
  * If the first word in a line in `:Purpose:` block is subroutine argument, Sphinx mixes up `:Purpose:` block and subroutine argument at this line. Put the word in `<`,`>` to avoid it.

## 3- Things to Do

* Inline comments after argument declarations are also processed by Sphinx, and should be used to supply a human-readable description of each argument; 
  ```fortran
  integer, intent(in) :: a    ! Description of A
  ```
  Note that Sphinx displays arguments one per line; place only one argument on a line that has a comment. In particular, this holds for declarations of variables that are members of a namelist. 

* use `intent` for all subroutine arguments. Do not put this information in the comments 
* Do not put revision or author information in the comments. This information is now stored in git, and is accessed there more efficiently.
* Start the text on a new line in `:Purpose:` block by using a line which only contains `!`.
