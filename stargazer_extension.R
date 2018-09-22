prep_table = function(my_stargazer,my_tabnote='',my_varlabels='',path,scale_notewidth=0.8){
     myn = names(my_varlabels) %>% gsub('_','\\_',.,fixed=T)
     for (i in 1:length(myn)){my_stargazer[grep(myn[i],my_stargazer,fixed=T)] = gsub(myn[i],my_varlabels[i],my_stargazer[grep(myn[i],my_stargazer,fixed=T)],fixed=T) }
     my_stargazer[(grep('Note',my_stargazer)-1)] = ''
     my_stargazer[grep('Note',my_stargazer)] = ''
     my_stargazer[grep('end\\{table\\}',my_stargazer)] = ''
     my_stargazer = c(my_stargazer,c('\\begin{minipage}{',scale_notewidth,'\\textwidth}',my_tabnote,'\\end{minipage}', '\\end{table}'))
     sink(file = path)
     cat(my_stargazer, sep = "\n")        
     sink()
}

### example
lapply(c("stargazer","data.table"),require,character.only=T)
dt = data.table(y=rnorm(100),x1=rnorm(100),x2=rnorm(100))

reg1 = lm(y ~ x1,data=dt)
reg2 = lm(y ~ x2,data=dt)

regtab = stargazer(reg1,reg2,type='latex')

prep_table(my_stargazer=regtab, # stargazer output, ensure that you used stargazer(..,type='latex')
           my_tabnote = 'This is your table note.', # your table notes.
           scale_notewidth=0.8, # control the width of the table notes
           my_varlabels = c('x1'='Variable 1','x2'='Variable 2'), # named vector providing a dictionary of R-variable names (included in Stargazer table) and the names you want to show up in the latex table
           path = 'Path where you want to save your tex-file.')
