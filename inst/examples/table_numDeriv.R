library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)

theme_set(theme_bw())

data.file <- "inst/examples/numDeriv.Rdata"
table.file <- "vignettes/table_numDeriv.tex"
load(data.file)



W1 <- filter(res, !(bench.expr %in% c("perm"))) %>%
  mutate(ms=bench.time/1000000) %>%
  select(-bench.time) %>%
  spread(bench.expr, ms) %>%
  gather(method, hessian, c(numDeriv, sparse)) %>%
  mutate(M=N*k+k,hessian.f=hessian/f, hessian.df=hessian/df,
         total=ifelse(method=="sparse",setup+hessian,hessian)) %>%
  gather(stat, time, c(setup, f, df, hessian, hessian.f, hessian.df, total)) %>%
  group_by(N, k, method, M, stat)  %>%
  summarize(mean=mean(time), sd=sd(time), lower=mean-1.96*sd, upper=mean+1.96*sd)


P1 <- ggplot(W1, aes(x=M, y=mean, color=method)) %>%
  + geom_line() %>%
  + geom_point() %>%
  + facet_wrap(~stat, scales="free")

##print(P1)


W2 <- filter(W1, stat %in% c("hessian","hessian.df")) %>%
  select(-c(lower,upper)) %>%
  gather(stat2, value, mean:sd) %>%
  dcast(N+k+M~stat+method+stat2,value.var="value") %>%
  arrange(M)

X <- xtable(W2, digits=3, display=c("s",rep("d",3),rep("fg",8)))
print(X, file=table.file, floating=FALSE, include.colnames=FALSE, include.rownames=FALSE)

