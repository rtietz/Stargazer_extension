catchKeyRow = function(strgz,row.label,indep.lab,type='text',stats = c('N','Rsq'),excl.indep.lab=T,text.out='str'){
     # strgz     : stargazer object
     # row.label : label for this row in future main table
     # indep.lab : name of independent variable whose row in stargazer table to catch.
     if (grepl('text',type)){
          ikv = grep(indep.lab,strgz,fixed = T)
          if (length(ikv)==0){stop('indep.lab not found in stargazer object.')}
          if (text.out=='data.table'){
               r1 = strgz[ikv]
               r1s = strsplit(r1,' ') %>% unlist()
               r1s = r1s[!r1s%in%'']
               r2 = strgz[(ikv+1)]
               r2s = strsplit(r2,' ') %>% unlist()
               r2s = r2s[!r2s%in%'']
               ts = rbind(r1s,c(' ',r2s))  
          } else {
               ts = rbind(strgz[ikv],strgz[(ikv+1)])
          }
          ts = cbind(c(paste0(row.label,'  ##|'),' '),ts) %>% as.data.frame()
          
          colnames(ts) = c('Variant','Coefficients')
          options(warn = -1) # just to avoid some meaningless NA-related warnings
          if (sum(grepl('N',stats))>0){
               iob = grep('Observations',strgz,fixed = T) 
               tob = gsub('Observations','',strgz[iob]) %>% str_trim(.) %>% gsub('  ',' ',.) %>% str_split(.,' ') %>% unlist() %>% as.numeric() %>% min(.,na.rm=T)
               ts$Obs = tob
          } 
          if (sum(grepl('Rsq',stats))>0){
               iob = grep('R2',strgz,fixed = T) 
               tob = gsub('R2','',strgz[iob]) %>% str_trim(.) %>% gsub('  ',' ',.) %>% str_split(.,' ') %>% unlist() %>% as.numeric() %>% min(.,na.rm=T)
               ts$R2 = tob
          }
          options(warn = 0)
     }
     if (grepl('latex',type)){
          ikv = grep(indep.lab,strgz,fixed = T)
          if (length(ikv)==0){stop('indep.lab not found in stargazer object.')}
          ts = c(strgz[ikv],strgz[(ikv+1)])
          if (isTRUE(excl.indep.lab)){ ts = ts %>% gsub(indep.lab,' ',.,fixed = T) %>% str_trim() %>% gsub('^&','',.) }
          if (!is.null(row.label)){
               ts[1]=paste0(paste0(row.label,' & '),ts[1])
               ts[2]=paste0(paste0(' ',' & '),ts[2])
          }
          options(warn = -1)
          if (sum(grepl('N',stats))>0){
               iob = grep('Observations',strgz,fixed = T) 
               tob = gsub('Observations','',strgz[iob]) %>% str_trim(.) %>% gsub('(&|\\\\)','',.) %>% str_split(.,' ') %>% unlist() %>% as.numeric() %>% min(.,na.rm=T)
               ts[1] = paste0(gsub('\\\\','',ts[1]),' & ',tob,' \\\\')
               ts[2] = paste0(gsub('\\\\','',ts[2]),' & ',' ',' \\\\')
          } 
          if (sum(grepl('Rsq',stats))>0){
               iob = grep('R$^{2}',strgz,fixed = T) 
               tob = gsub('R$^{2}','',strgz[iob]) %>% str_trim(.) %>% gsub('(&|\\\\)','',.) %>% str_split(.,' ') %>% unlist() %>% as.numeric(.) %>% min(.,na.rm=T)
               ts[1] = paste0(gsub('\\\\','',ts[1]),' & ',tob,' \\\\')
               ts[2] = paste0(gsub('\\\\','',ts[2]),' & ',' ',' \\\\')
          }
          options(warn = 0)
     }
     return(ts)
}
regs_by_controls = function(depvar,indeps,ctrls,regdat,type='text',
                            print.out=T,
                            spec.no.ctrls = T,cum.ctrls = F,
                            col.with.all=c(T,T),method ='lm',fespec='|0|0|0',
                            title='',label='',font.size='normalsize',mySEs=NULL,
                            dep.var.labels.include=F,dep.var.caption='',table.placement='H'){
     # dat            : data to run on
     # main_dep       : name of dependent variable
     # main_indeps    : name(s) of one or more independent variables, always included in regression
     # ctrls          : name(s) of one or more control variables, included one by one or cumulatively.
     # col.with.all   : logical vector, 1st position: incl column with all control variables? 2nd: print it?
     # cum.ctrls      : logical, shall controls be included one by one or cumulative?
     # mySEs          : function specifying the SEs, input: model-object, output: vector with scalar values of standard errors
     regs=list()
     st = ifelse(isTRUE(spec.no.ctrls),0,1)
     for (jj in st:length(ctrls)){
          ccs = ifelse(isFALSE(cum.ctrls),ctrls[jj],paste0(ctrls[1:jj],collapse = '+'))
          regs[[paste0(jj)]] = lm(as.formula(paste0(depvar,'~',paste0(indeps,collapse = '+'),ifelse(jj>0,paste0(' + ',ccs),''))),data=regdat) 
          # if (method=='lm'){
          #      regs[[paste0(jj)]] = lm(as.formula(paste0(depvar,'~',paste0(indeps,collapse = '+'),ifelse(jj>0,paste0(' + ',ccs),''))),data=regdat)    
          # } else if (method=='felm'){
          #      regs[[paste0(jj)]] = felm(as.formula(paste0(depvar,'~',paste0(indeps,collapse = '+'),ifelse(jj>0,paste0(' + ',ccs),''),fespec)),data=regdat)
          # }
     }
     myomit=NULL
     if (col.with.all[1] & isFALSE(cum.ctrls)){
          if (isFALSE(col.with.all[2])){
               setnames(regdat,ctrls,paste0('tmp_',ctrls))
               ctrls=paste0('tmp_',ctrls)
               myomit = '^tmp_' #c(omit,'^tmp_') 
          }
          regs[[paste0(jj+1)]] = lm(as.formula(paste0(depvar,'~',paste0(indeps,collapse = '+'),' + ',paste0(ctrls,collapse = '+'))),data=regdat) 
     }
     if (isFALSE(print.out)){suppressWarnings(sink(''))}
     # define standard errors if not specified:
     if (is.null(mySEs)){ mySEs = function(x) coeftest(x,vcov=NeweyWest(x,prewhite = F,adjust=F,lag = 4))[,2]}
     t = stargazer(regs,type=type,
                   model.numbers = F,dep.var.labels.include = dep.var.labels.include,dep.var.caption = dep.var.caption,omit=myomit,
                   table.placement = table.placement,title = title,label = label,font.size = font.size,
                   #se=lapply(regs,function(x) coeftest(x,vcov=vcovHAC(x))[,2]),
                   se=lapply(regs,mySEs),
                   omit.stat=c("LL","ser","f",'adj.rsq'))
     if (isFALSE(print.out)){sink()}
     return(t)
}
regs_by_cntrls_by_smpl = function(mydepvar,my_key.indep.var,myindeps,myctrls,mydat,mytype='text',
                                  smpl.dim,smpls,mySEs=NULL,
                                  my.spec.no.ctrls = T,spec.all.ctrls=T,cum.ctrls = F,stats = c('N','Rsq'),
                                  my.print.out=F,mytitle='',my_tabnote='',notewidth=0.975,
                                  dep.var.caption='',font.size='normalsize',label='',table.placement='H',path=NULL,
                                  tabnote.fontsize = 'footnotesize',scalebox=NULL,
                                  smpl.dim.label ='Sample',
                                  cntrl.dim.label ='Add-control'){
     # smpl.dim : string, name of variable that is dimension for subsampling
     # smpls    : string vector indicating different sub-samples, each element in smpls is a value of smpl.dim.
     # scalebox : numeric scaling factor if want to have final table wrapped into scalebox.
     if (grepl('text',mytype)){
          tabs=data.table(stringsAsFactors = F)
     } else if (grepl('latex',mytype)){
          tabs = as.character()
     }
     if (is.null(names(smpls))){names(smpls)=smpls}
     if (is.null(names(myctrls))){names(myctrls)=myctrls}
     
     # estimate specifications for different sub-samples, collect model-objects in list.
     for (i in 1:length(smpls)){
          cur.smpl = smpls[i]
          datc = mydat[ get(smpl.dim) %in% cur.smpl ]
          t = regs_by_controls(depvar = mydepvar,indeps=myindeps,ctrls=myctrls,type=mytype,
                               regdat=datc,
                               col.with.all=c(spec.all.ctrls,F),spec.no.ctrls = my.spec.no.ctrls,
                               title=mytitle,label=label,font.size = font.size,mySEs = mySEs,
                               cum.ctrls = cum.ctrls,print.out = my.print.out,
                               dep.var.labels.include=F,dep.var.caption=dep.var.caption,table.placement)
          #rm(datc)
          ts = catchKeyRow(t,names(cur.smpl),indep.lab=my_key.indep.var,type=mytype,stats)
          if (grepl('text',mytype)){
               tabs = rbind(tabs,ts) 
          } 
          if (grepl('latex',mytype)){
               tabs = c(tabs,ts)
          }
     }
     if (grepl('text',mytype)){
          lstr = names(myctrls)
          for (yy in 1:length(lstr)){ if(nchar(lstr[yy])>7) { lstr[yy] = substring(lstr[yy],1,7)}}
          
          if (isTRUE(my.spec.no.ctrls)){lstr = c('None',lstr)}
          if (isTRUE(spec.all.ctrls)){lstr = c(lstr,'All')}
          nspacs = 7*(1/max(c(round(length(lstr)*0.25,0),1)))
          lstr[1] = paste0(paste0(rep(' ',nchar(my_key.indep.var)),collapse = ''),lstr[1])
          lstrdt = c(cntrl.dim.label,paste0(lstr,collapse =paste0(rep(' ',nspacs),collapse = '')),rep('-',length(stats)))%>%t() %>% data.table()
          (tabs= rbind(tabs,lstrdt,use.names=F,fill=F))
          return(tabs)
     }
     if (grepl('latex',mytype)){
          # Row before notes indicating controls:
          lstr = names(myctrls)
          for (yy in 1:length(lstr)){ if(nchar(lstr[yy])>7) { lstr[yy] = substring(lstr[yy],1,7)}}
          if (isTRUE(my.spec.no.ctrls)){lstr = c('None',lstr)}
          if (isTRUE(spec.all.ctrls)&isFALSE(cum.ctrls)){lstr = c(lstr,'All')}
          lstrdt = paste0(c(cntrl.dim.label,lstr %>% gsub('_','.',.,fixed=T),rep(' ',length(stats))),collapse = '&')
          lstrdt = c('\\hline \\\\[-1.8ex]',paste0(lstrdt,'\\\\'),"\\hline ","\\hline \\\\[-1.8ex] ")
          NoCols = length(lstr) + length(stats) + 1
          # Table Header:
          st0  = grep('\\begin{table}',t,fixed = T)
          if (dep.var.caption%in%''){st1  = grep('\\begin{tabular}',t,fixed = T)+2} else (st1  = grep('\\cline',t,fixed = T))
          stt = t[st0:(st1)]
          # insert scalebox:
          if (!is.null(scalebox)){stt = c(stt[1:(grep('\\begin{tabular}',stt,fixed=T)-1)],paste0('\\scalebox{',scalebox,'}{'),stt[grep('\\begin{tabular}',stt,fixed=T):length(stt)])}
          # adjust column-specification for tabular-environment:
          stt[grep('\\begin{tabular}',stt,fixed=T)] = gsub('lc+',paste0(c('l',rep('c',(NoCols-1))),collapse = ''),stt[grep('\\begin{tabular}',stt,fixed=T)])
          if (!dep.var.caption%in%''){ 
               ix = grep('multicolumn',stt)
               stt[ix] = paste0(gsub('\\\\ ','',stt[ix],fixed=T),'&',paste0(stats,collapse=' & '),' \\\\')
               stt[ix] = paste0(smpl.dim.label,stt[ix])
               stt[ix+1] = paste0(stt[ix+1],'\\\\[-1.8ex]') 
          }         
          # Table Ending:
          ed = c('\\end{tabular}',ifelse(!is.null(scalebox),'}',''),c(paste0('\\begin{minipage}{',notewidth,'\\textwidth}\\',tabnote.fontsize),my_tabnote,'\\end{minipage}', '\\end{table}'))
          tex.tab = c(stt,tabs,lstrdt,ed)
          tex.tab = tex.tab[!tex.tab%in%'']
          if (!is.null(path)){
               sink(file = path)
               cat(tex.tab, sep = "\n")        
               sink()  
          }
          return(tex.tab)
     }
}
prep_table = function(my_stargazer,my_tabnote='',my_varlabels='',repl.caption = F,ignore.tabnote=F,path,scale_notewidth=0.95,tbn.font.size='footnotesize'){
     # my_stargazer = stargazer object to be finalized.
     # my_varlabels = named character vector, with variable name and corresponding label to be used in table.
     # repl.caption = default is logical FALSE, may be replaced with string for title. work-around for stargazer bug occurring with argument "keep".
     if (!sum(grepl('\\begin',my_stargazer,fixed = T))>0){stop('Table must be latex table.')}
     myn = names(my_varlabels) %>% gsub('_','\\_',.,fixed=T) #%>% gsub('.','.',.,fixed=T)
     for (i in 1:length(myn)){my_stargazer[grep(myn[i],my_stargazer,fixed=T)] = gsub(myn[i],my_varlabels[i],my_stargazer[grep(myn[i],my_stargazer,fixed=T)],fixed=T) }
     my_stargazer[(grep('Note',my_stargazer)-1)] = ''
     my_stargazer[grep('Note',my_stargazer)] = ''
     if (repl.caption!=F){
          my_stargazer[grep('\\caption',my_stargazer,fixed=T)] = paste0('\\caption{',repl.caption,'}')
     }
     if (ignore.tabnote==F){
          my_stargazer[grep('end\\{table\\}',my_stargazer)] = ''
          my_stargazer = c(my_stargazer,'\\vspace{0.15cm}',c(paste0('\\begin{minipage}{',scale_notewidth,'\\textwidth}\\',tbn.font.size),my_tabnote,'\\end{minipage}', '\\end{table}'))
     }
     my_stargazer = my_stargazer[!my_stargazer%in%'']
     sink(file = path)
     cat(my_stargazer, sep = "\n")        
     sink()
}
