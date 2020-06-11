############################# Examples:
source('aux_fcts.R')
lapply(c("stargazer","data.table","stringr"),require,character.only=T)
#******* Example 1:
dt = data.table(y=rnorm(100),x1=rnorm(100),x2=rnorm(100))

reg1 = lm(y ~ x1,data=dt)
reg2 = lm(y ~ x2,data=dt)
regtab = stargazer(reg1,reg2,type='latex')

prep_table(my_stargazer=regtab, # stargazer output, ensure that you used stargazer(..,type='latex')
           my_tabnote = 'This is your table note.', # your table notes.
           scale_notewidth=0.8, # control the width of the table notes
           my_varlabels = c('x1'='Variable 1','x2'='Variable 2'), # named vector providing a dictionary of R-variable names (included in Stargazer table) and the names you want to show up in the latex table
           path = 'Your Path + Your tex-filename.tex')

#******* Example 2:
lapply(c('lmtest','sandwich'),require,character.only=T)
dt = data.table(x2=rnorm(2*100),x3=rnorm(2*100),x4=rnorm(2*100),group=rep(c('A','B'),100))
dt[,x1:=rnorm(2*100)+0.6*x3]
dt[,y:=rnorm(2*100)+0.2*x2+0.9*x3]

regs_by_cntrls_by_smpl(mydepvar='y',             # your independent variable
                       my_key.indep.var='x1',    # the key variable whose coefficient you focus on
                       myindeps=c('x1'),         # independent variables always included
                       myctrls=c('x2','x3','x4'),# control variables, will be added one by one (indivudally or cumulatively)
                       smpl.dim='group',         # the variable indicating the sub-samples (in rows)
                       smpls =c('A','B'),        # the row-labels of sub-samples you want to consider
                       mytype = 'text',
                       mydat=dt)

