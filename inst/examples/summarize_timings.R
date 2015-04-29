library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)

theme_set(theme_bw())

data.file <- "timings3.Rdata"
load(data.file)

tab <- mutate(res, ms=bench.time/1000000) %>%
  dcast(N+k+T+ord+bench.rep+ncolors~bench.expr, value.var="ms")  %>%
  mutate(nvars=N*k+k, hess_f=hess/f, hess_df=hess/df,
         hess_k=hess/k, hess_ncol = hess/ncolors, total=setup+hess,
         total_f=total/f, total_df = total/df, df/nvars)

tab2 <- gather(tab, stat, ms, c(setup:hess,hess_f:total_df)) %>%
  group_by(N, k, T, ord, stat, nvars) %>%
  summarize(mean=mean(ms), sd=sd(ms), lower=mean-1.96*sd, upper=mean+1.96*sd)


## D2 <- subset(data.frame(tab2), stat %in% c("nvars","setup","colors","f","df","hess",
##                                "hess_df","total","total_f","total_df"))

D2 <- data.frame(tab2)

P2 <- ggplot(D2, aes(x=N,y=mean, color=as.factor(k))) %>%
  + geom_line() %>%
  + scale_x_continuous("Number of heterogeneous units") %>%
  + scale_y_continuous("Computation time (milliseconds)") %>%
  + facet_wrap(~stat, scales="free")


P3 <- ggplot(filter(D2,k==4), aes(x=N)) %>%
  + geom_line(aes(y=mean), linetype=1) %>%
  + geom_line(aes(y=upper), linetype=3) %>%
  + geom_line(aes(y=lower), linetype=3) %>%
  + facet_wrap(~stat, scales="free")



print(P2)
