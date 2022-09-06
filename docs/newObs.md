# New Observation : A Guideline

These guidelines describe how to add a new observation familly in MIDAS.
It assumes those are encoded in SQLite format.

1. A new observation file line description
  ```
  clvalu(N+1) ='obs{new}'
  ```
  and a familly is associated with another line
  ```
  cfami'(N+1)='NE'
  ```
  are added in in the `subroutine obsf_setupFileNames` (`obsfiles_mod.f90`), 
  where `{new}` is a placeholder for the new observation prefix and `NE` is the 
  two characters identifier of the familly.

2. observation operator subroutines
    * non-linear operator: `subroutine oop_{new}_nl`
    * tangeant linear operator: `subroutine oop_H{new}`
    * adjoint operator: `subroutine oop_HT{new}`
   are added in `obsoperators_mod.f90`

3. add what is needed to include the new observation in cost function evaluation
    * proper calls to `oop_{ppp,sfc,sst,zzz,tovs,chm}_nl` in `inn_computeInnovation` (`innovation_mod.f90`)
    * in `cfn_sumJo` (`costfunction_mod.f90`):
        * `real(8) dljo{new}` variable initialization
        * `dljo{new} = 0.d0`
        * `case('NE'); dljo{new} = dljo{new} + pjo_1`
        * `call mpi_allreduce_sumreal8scalar( dljo{new}, "GRID")`

4. add the new familly in `ofl_familyList` (`obsFamilyList_mod.f90`)
5. the `obs_oer` column of `obspacedata` needs to be defined for the new observation
   by modifying existing subroutines or adding new ones in `obserrors_mod`

6. `obserr` file needs to be modified
    * the new familly need to be added in `bgcheck`
    * a new variable may be needed in the background state
