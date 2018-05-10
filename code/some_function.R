#完全缓解 CR: 尿蛋白 < 0.5
#部分缓解 PR: 没有CR 的情况下 尿蛋白改变 > 50% 且 尿蛋白 < 3.5
#复发： 有CR 的情况下， 在后面的额随访情况中 尿蛋白升高到 3.5 g/24h

library(tidyverse)

#====是否有 CR====
to_CR <- function(df){
  df$treat_month[df$urine_protein < 0.3][1]  %>%#小于 0.3g
    if_else(!is.na(.), ., 0)
}

#====是否有PR====
to_PR <- function(df) {
  df_PR <-
    df %>%
    mutate(deta_p = (urine_protein[1] - urine_protein) / urine_protein[1]) %>%
    mutate(C1 = deta_p > 0.5) %>%#改变超过50%
    mutate(C2 = urine_protein < 3.5) %>%
    mutate(c_1_2 = C1 & C2)
  
  df_PR$treat_month[df_PR$c_1_2][1] %>% 
    if_else(!is.na(.), ., 0)
}

#先判断没有CR, 再是否有PR
too_PR <- function(df){
  if_else(to_CR(df)!=0,0,to_PR(df))
}


#==== 随访时间====
fl_time <- function(df){
  idx <- dim(df)[1]
  df$treat_month[idx]
}

#==== 是否有复发====
to_relapse <- function(df){
  mon <- to_CR(df)#完全缓解的月份
  df_R <- 
    df %>% 
    filter(treat_month > mon) %>% 
    mutate(if_R= urine_protein >=3.5)
  
  df_R$treat_month[df_R$if_R][1] %>% 
    if_else(!is.na(.), ., 0)
}

#先判断没有CR, 再是否有复发
too_relapse <- function(df){
  if_else(to_CR(df)==0,0,to_relapse(df))
}

#==== 字符变量 to 因子变量 ====
cha_to_fac <- function(df){
  id <- sapply(df,is.character)
  df[id] <- lapply(df[id],as.factor)
  df
}

#==== 插入中位数或者 平均值 ====
if_impute <- function(x) {
  if(shapiro.test(x)$p.value > 0.05){
    x <- impute(x,fun=mean)
  }
  else{
    x <- impute(x,fun=median)
  }
  return(x)
}

to_impute <- function(df){#🔥 map impute
  id <- sapply(df, is.numeric)
  df[id] <- lapply(df[id], FUN=if_impute)
  df[id] <- lapply(df[id], as.numeric)
  df
}

#===== 获取数值变量 ====
get_num <- function(df){
  id <- sapply(df,is.numeric)
  df[id]
}

#==== 获取非数值变量 =====
get_unum <- function(df){#from Tmisc::unfactor
  id <- sapply(df, is.numeric)
  df[-id]
}

#==== 获取因子变量 =====
get_fac <- function(df){
  id <- sapply(df,is.factor)
  df[id]
}

#==== 两组比较 参数或者非参数检验的P值 ====
get_p <- function(x,fac,paired=FALSE,digits = 3){
  if(shapiro.test(x)$p.value >0.05){
    test<- t.test(x~fac,paired = FALSE)
  }
  else{
    test <- wilcox.test(x~fac,paired = FALSE)
  }
  return(signif(test$p.value,digits = digits))
}


#==== 平均值的统计描述=====
get_stat_mean <- function(x,na.omit=TRUE,digits=3){
  library(stringr)
  library(dplyr)
  if(na.omit)
    x <- x[!is.na(x)]
  m <- signif(mean(x),digits = digits)
  s <- signif(sd(x),digits = digits)
  stat <- str_c(m,"±",s)
  return(stat)
}

#===== 中位数统计描述====
get_stat_median <- function(x,na.omit=TRUE,digits=3){
  library(stringr)
  library(dplyr)
  if(na.omit)
    x <- x[!is.na(x)]
  m <- signif(median(x),digits = digits)
  IQR_1 <- signif(quantile(x)[2],digits = digits)
  IQR_3 <- signif(quantile(x)[4],digits = digits)
  stat <- str_c(m,"(",IQR_1,"-",IQR_3,")")
  return(stat)
}

#====数值变量的统计描述 mean/median====
get_stat_num <- function(x,na.omit=TRUE,digits=3){
  library(stringr)
  library(dplyr)
  if(na.omit)
    x <- x[!is.na(x)]
  if(shapiro.test(x)$p.value >0.05){
    m <- signif(mean(x),digits = digits)
    s <- signif(sd(x),digits = digits)
    stat <- str_c(m,"±",s)
  }
  else{
    m <- signif(median(x),digits = digits)
    IQR_1 <- signif(quantile(x)[2],digits = digits)
    IQR_3 <- signif(quantile(x)[4],digits = digits)
    stat <- str_c(m,"(",IQR_1,"-",IQR_3,")")
  }
  return(stat)
}

# =====是否为正态分布 向量值=====
t_or_w <- function(df){
    as.list(df) %>%
    map(~ shapiro.test(.x)$p.value > 0.05) %>%
    unlist()
}

#====因子变量的统计描述====
get_stat_fac <- function(x,fac=fac){
  a <- table(x) %>%
    prop.table()%>%
    as.matrix() %>%
    apply(FUN=signif,digit=3,MARGIN = 2) %>%
    t() %>%
    t()*100 
  a <- as.data.frame(a)
  a$Freq <- as.data.frame.table(table(x)) %>% .[[2]] 
  a <- a %>%
    mutate(total= str_c(Freq,"(",V1,"%",")")) %>%
    select(total)
  
  b_0 <- as.data.frame.table(xtabs( ~ x + fac)) %>%
    rename(count = Freq)
  b <-
    xtabs(~ x + fac) %>%
    prop.table(., margin = 2) %>%
    as.data.frame.table() %>%
    mutate(Freq = signif(Freq, digit = 3)) %>%
    mutate(Freq = Freq * 100) %>%
    left_join(b_0, by = c("x", "fac")) %>%
    mutate(n_Fre = str_c(count, "(", Freq, "%", ")")) %>%
    select(x, fac, n_Fre) %>%
    spread(key = fac, value = n_Fre)
  
  c <- xtabs( ~ x + fac) %>%
    chisq.test()
  c <- c$p.value
  
  b$Total <- a[[1]]
  b$p_value <- NA
  b$p_value[1] <- c
  b <- select(b, x, Total, everything())
}


# ===== eGFR 计算函数====
eGFR_fun <- function(sex,scr,age) {
  if (sex == "male"){
    X=min(scr/(88.41*0.9),1)^-0.411
    Y=max(scr/(88.41*0.9),1)^-1.029
    Z=0.993^age
    result <- 141*X*Y*Z 
  }
  else{
    X=min(scr/(88.41*0.7),1)^-0.329
    Y=max(scr/(88.41*0.7),1)^-1.029
    Z=0.993^age
    result <- 141*X*Y*Z*1.018
  }
  return(result)
}

#==== reorder factor====
re_ord <- function(fac){
  fac <- factor(fac, levels = c("1","0"))
  return(fac)
}

#====calculate the  month index====
get_idx <- function(df=flow_up, m="12"){
  df %>% 
    filter(df$rep_id!="45") %>% #45 号病人的随访情况太少
    filter(!is.na(urine_protein)) %>% 
    filter(treat_month==m) %>% 
    dplyr::select(rep_id) %>% .[[1]]
}

#====choose time====
choose_time <- function(CR_time,PR_time){
  library(dplyr)
  case_when(
    CR_time > 0 & PR_time == 0 ~ CR_time + PR_time,
    CR_time == 0 & PR_time > 0 ~ CR_time + PR_time,
    CR_time > 0 & PR_time > 0 & CR_time > PR_time  ~ PR_time,
    CR_time > 0 & PR_time > 0 & CR_time < PR_time ~ CR_time,
    CR_time > 0 & PR_time > 0 & CR_time == PR_time ~ PR_time
  )
}
