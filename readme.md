## wrapper functions for stargazer

To simplify creating pretty tables for common cases of combining several regressions.

#### Structure of code repository

1. aux_fcts.R  : contains auxiliary functions that manipulate stargazer output
2. example.R   : showcases how to use these functions

#### Explanation
The key auxiliary functions are, in increasing order of complexity:

1. `prep_table.R` just prettifies the stargazer objects, e.g. by replacting variable names with labels. Some of these functionalities can be implemented with stargazer directly, but I preferred my syntax.
2. `regs_by_controls.R` assembles a table summarizing specifications with different control variables.
3. `regs_by_cntrls_by_smpl.R` is based on `regs_by_controls.R` and provides a summary table for different specifications over different subsamples.