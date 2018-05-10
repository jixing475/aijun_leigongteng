#raw data clean
#====import libray====
library(tidyverse)
library(readxl)
library(Tmisc)
library(Hmisc)
library(janitor)
library(purrr)
library(clinPK)
source("./code/some_function.R")

#====import data====
rawdata <- read_excel("data/rawdata.xlsx",
                      col_types = c("text", "text", "text",
                                    "numeric", "numeric", "numeric",
                                    "text", "text", "text", "text", "text",
                                    "text", "numeric", "text", "text",
                                    "text", "text", "text", "text", "text",
                                    "blank", "text", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric", "numeric", "numeric",
                                    "numeric", "numeric", "text", "text"))

#====get one time data====
colnames(rawdata)
one_col <- c("姓名","编号","性别","年龄","身高/cm","体重/Kg","吸烟",
              "饮酒", "高血压","糖尿病","乙肝","其他疾病","NS病程/月",
              "既往激素","既往CTX","既往CNI","既往骁悉","爱若华","ACEI/ARB")
one_time <-
  rawdata %>%
  dplyr::select(one_col) %>%
  clean_names() %>%
  remove_empty_rows() %>%
  remove_empty_cols()
glimpse(one_time)
colnames(one_time) <- c("name", "id", "sex", "age", "height", "weight", "smoke",
                        "drink", "hypertension", "diabetes", "Hepatitis_B", "other_disease", "NS_bingcheng",
                        "jiwang_hormone", "jiwang_CTX", "jiwang_CNI", "jiwang_xiaoxi", "airuohua", "ACEI_ARB")
Name <- one_time$name
ID <- one_time$id
## missing value
library(naniar)
gg_miss_var(one_time)

#=====get flow up data====
colnames(rawdata)
flow_col <- c("姓名","编号","激素用量mg/隔日","雷公藤/mg","治疗月 m","血肌酐 umol/L","Cys-C mg/L",
              "UA umol/L","CHOL mmol/L","LDL mmol/L","ALB g/L","GLU mmol/L","尿蛋白定量 g/24h")
flow_up <-
  rawdata %>%
  dplyr::select(flow_col) %>%
  clean_names() %>%
  remove_empty_rows() %>%
  remove_empty_cols()
colnames(flow_up) <- Hmisc::Cs(name,id,hormone,Tripterygium,treat_month,SCreatinine,cys,ua,chol,ldl,alb,glu,urine_protein)

flow_up[1230:1232,1:5]
flow_up <- flow_up[-1231,]#1231是一个补充数据 造成这个病人多了一个随访数据

#判断俩数据集病例数是否相等
dim(one_time)[1]*15==dim(flow_up)[1]


#====缺少初次尿蛋白资料====
flow_up$rep_id <-
  flow_up$id[!is.na(flow_up$id)] %>%
  rep(each=15)

#=================================================
id_no_UrineP <-
  flow_up %>% filter(treat_month == 0) %>%
  filter(is.na(urine_protein)) %>%
  .$id
flow_up  <-
  flow_up %>% filter(!rep_id %in% id_no_UrineP)

#====尿蛋白数据少于两次 ====
id_UrineP_less_than_2 <-
  flow_up %>%
  split(.$rep_id) %>% #切割数据
  map( ~ dplyr::select(.x, treat_month, urine_protein)) %>% #选择变量
  map( ~ filter(.x,!is.na(.x$urine_protein))) %>%
  map( ~ dim(.x)[1] < 2)  %>%
  unlist()
id_UrineP_less_than_2 <- names(id_UrineP_less_than_2)[id_UrineP_less_than_2]
flow_up  <-
  flow_up %>% filter(!rep_id %in% id_UrineP_less_than_2)

#==== create baseline data====
baseline <- inner_join(one_time, flow_up, by="id")
gg_miss_var(baseline)

#====remain data====
#col remain
col_remain <- which(colMeans(is.na(baseline)) < 0.2)
baseline <- baseline[, col_remain]
#remove name.x treat_month Tripterygium
baseline <-
  baseline %>%
  dplyr::select(-name.y,-treat_month, -Tripterygium,-hormone)
#row remain
ID_remain <- which(rowMeans(is.na(baseline)) <0.2)
baseline <- baseline[ID_remain,]
gg_miss_var(baseline)

#=========mice::impute data=============
library(mice)
clean_baseline <- mice(data = baseline,m=5,seed = 123)  %>% complete()
gg_miss_var(clean_baseline)

clean_baseline$jiwang_xiaoxi[which(is.na(clean_baseline$jiwang_xiaoxi))] <- "0"
gg_miss_var(clean_baseline)

#==== mutate eGFR====
clean_baseline$sex <-
  factor(
    clean_baseline$sex,
    levels = c("1", "2"),
    labels = c("male", "female")
  )

clean_baseline$eGFR <- clinPK::calc_egfr(method = "cockcroft_gault",sex = clean_baseline$sex,age = clean_baseline$age,
                                         scr = clean_baseline$SCreatinine,scr_unit = "micromol/L",height = clean_baseline$height,weight = clean_baseline$weight) %>%
  .$value
#==== create group variable: PR CR reverse====
flow_up_list <-
  flow_up %>%
  split(.$rep_id) %>% #切割数据
  map( ~ dplyr::select(.x, treat_month, urine_protein)) %>% #选择变量
  map( ~ filter(.x,!is.na(.x$urine_protein)))

#=========⭐group=======
group <-
  data.frame(
    id = names(flow_up_list),
    CR_time = flow_up_list %>% map( ~ to_CR(.x)) %>% unlist(),
    PR_time = flow_up_list %>% map( ~ too_PR(.x)) %>% unlist(),
    flow_time = flow_up_list %>% map( ~ fl_time(.x)) %>% unlist(),
    relapse_time = flow_up_list %>% map( ~ too_relapse(.x)) %>% unlist()
  )

#====recode column: include rename colname and type of column====
clean_baseline <- clean_baseline %>%
  map_if(is.character, as.factor) %>% as.data.frame()

group <-
  group %>%
  mutate(
    group = case_when(
      CR_time > 0 & PR_time == 0 ~ "CR",
      CR_time == 0 & PR_time > 0 ~ "PR",
      CR_time == 0 & PR_time == 0 ~ "NR"
    )
  )

#NR_id <- group$id[group$group=="NR"] %>% as.character()
group <-
  group %>%
  mutate(time = case_when(
    group == "CR" ~ CR_time,
    group == "PR" ~ PR_time,
    group == "NR" ~ flow_time
  ))

TableOne_data_with_id <-
  inner_join(clean_baseline, group, by="id") %>%
  dplyr::select(-name.x, -rep_id, -CR_time, -PR_time, -flow_time, -relapse_time, -time)
TableOne_data_with_id$group <- as.factor(TableOne_data_with_id$group)
save(TableOne_data_with_id, file = "data/TableOne_data_with_id.Rdata")

TableOne_data <-
  inner_join(clean_baseline, group, by="id") %>%
  dplyr::select(-name.x, -id, -rep_id, -CR_time, -PR_time, -flow_time, -relapse_time, -time)
TableOne_data$group <- as.factor(TableOne_data$group)

#====xray explore====
skimr::skim(clean_baseline)
save.image("data/clean_data.RData")
save(TableOne_data, file = "data/TableOne_data.Rdata")
