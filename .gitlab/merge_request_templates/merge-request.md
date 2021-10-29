### Mandatory checklist:

<!--For each point below, choose 'YES' or 'NO' -->

* new functionality?  **YES/NO**
* fix bug in existing functionality? **YES/NO**
* changes to namelist variables? (addition/removal/modification) **YES/NO**
* changes to input and output files? (new file/filename change/removal of file) **YES/NO**
* change to results (contents of output files)? **YES/NO**
* changes respect the coding standards (follow link below to "Coding standards Top 10")? **YES/NO**

    * https://wiki.cmc.ec.gc.ca/wiki/MIDAS/Coding_Standards_Top_10

### Description of changes:

<!--The text here should describe how the change was implemented.-->
<!--Detail here the changes answered as YES in the previous section-->

### Description of the impact on the results

<!--

Describe how the results are affected by the code introduced in this
merge request:

 1. quick explanation of why the results are affected:
   * order of obs changed
   * affecting cost function calculation
 2. which programs are affected (e.g. `midas-var.Abs`)
 3. in which configurations are the impacts seen
   * e.g. all operational NWP systems using `midas-var.Abs` program
 4. how significant are the changes
   * no impact at all (backward compatible)
   * minor, only due to numerical round-off error or
   * major like impacting the meteorological evaluation
      * if so, in which data assimilation experiments were the impact
        of the changes evaluated

-->

Delete this line and replace it with your description of the impact on
the results based on the text above.

### Addition to CHANGELOG:

<!--Some oneliners describing changes for the whole merge-request-->
<!--That information will be added to the 'CHANGELOG.md' file-->
<!--Put any information relevant to the user, especially non-backward compatible changes-->
<!--   * new functionality  -->
<!--   * Namelist variables -->
<!--   * input/output files -->
<!--   * results            -->

**Don't forget to add a description of your changes in the
  'Unreleased' section of `CHANGELOG.md`.**

If the changes have impact on the results, please insert
  * "minor impact on results" or
  * "major impact on results".

### Resolved issues:

<!--Put the list of issues that this merge request resolves-->
Closes #...

<!--(choose one of the following labels)-->
/label ~"Minor Functionality"
/label ~"Major Functionality"
/label ~"Technical change" 
/label ~Optimization
/label ~Documentation
/label ~Bug
