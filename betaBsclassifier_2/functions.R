#script con todas las funciones necesarias para el clasificador Bayesiano

pacman::p_load( tidyverse, latticeExtra, RColorBrewer,
               collapsibleTree, stringi, writexl, 
                lobstr, readxl, gtools, pROC, caret)

n.sint <- 10

df <- readxl::read_excel("www/BD aedes_newalgoritmo.xls") %>% 
  mutate_at(vars(1:n.sint), .funs = function(i){ifelse(i >1, 1, i )}) %>% 
  filter(ratioln <= 3.5) %>% 
  filter(leucos <= 20000) %>% 
  filter(plaque <= 500000)

#function to generate the vectors with all permutations depending of the number of symptoms
vectors <- function(n.sint){
  vectors <- gtools::permutations(n=2, r=n.sint, v=c(0,1), repeats.allowed = T)
  vectors <- vectors[nrow(vectors):1, ]
  return(data.frame(vectors))
}

symp_nms <- names(df)[1:n.sint]
dx_nms <- "dxdengue0"
hemo_vars <- names(df)[!names(df) %in% c(symp_nms, dx_nms)]

global_parms <- df %>%
 group_by(dxdengue0) %>% 
  summarise_at(vars(hemo_vars), list("mean", "sd"))


vec <- vectors(n.sint)

comb <- vec %>% 
  tidyr::unite(., comb, sep = "", remove = F)

in_df <- df %>% 
 tidyr::unite(., comb, 1:n.sint, sep = "", remove = F ) %>% 
  mutate(id = as.character( 1:nrow(.)))



arbol_only_symp <- in_df %>%
  dplyr::select(-id) %>% 
  group_by(comb) %>% 
  dplyr::summarise(., dx_pos = sum(dxdengue0 == 1),
                   dx_neg = sum(dxdengue0 == 0),
                   total  = n()) %>% 
  ungroup() %>% 
  full_join(comb, .) %>%
  dplyr::mutate_all(.funs = function(i){ifelse(is.na(i), 0, i)}) %>% 
  dplyr::mutate(score_pos = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
    mean(rbeta(10000, 1+x, 1+y))
  }),
  score_neg = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
    mean(rbeta(10000, 1+y, 1+x))
  }),
  test_res = case_when(
    score_pos > score_neg ~ "1",
    score_neg > score_pos ~ "0",
    TRUE ~ "2"
  ) ) 

 # validacion solo sintomas
global_mtrs <- left_join(in_df, arbol_only_symp %>% dplyr::select(comb, score_pos, score_neg, test_res)) %>%
  dplyr::select(dxdengue0, test_res ) %>% 
  filter(test_res %in% c("0", "1")) %>% 
  dplyr::mutate_all(as.factor) %>%
  confusionMatrix(data = .$test_res, reference = .$dxdengue0, positive = "1")

k <- 5
cv_res <- lapply(1:k, function(k){
  
  test <- in_df %>% 
    sample_n(size = round(0.2*nrow(.)))
  
  train <- in_df %>%
    filter(!id %in% test$id) 
  
  only_symp <- train %>%
    select(-id) %>% 
    group_by(comb) %>% 
    dplyr::summarise(., dx_pos = sum(dxdengue0 == 1),
                     dx_neg = sum(dxdengue0 == 0),
                     total  = n()) %>% 
    ungroup() %>% 
    full_join(vectors, .) %>%
    dplyr::mutate_all(.funs = function(i){ifelse(is.na(i), 0, i)}) %>% 
    dplyr::mutate(score_pos = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
      mean(rbeta(1000, 1+x, 1+y))
    }),
    score_neg = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
      mean(rbeta(1000, 1+y, 1+x))
    }),
    test_res = case_when(
      score_pos > score_neg ~ "1",
      score_neg > score_pos ~ "0",
      TRUE ~ "2"
    ) )
  
  train_mts <- left_join(train, only_symp %>% 
                           select(comb, score_pos, score_neg)) %>% 
    dplyr::mutate(test_res = case_when(
      score_pos > score_neg ~ "1",
      score_neg > score_pos ~ "0",
      TRUE ~ "2"
    ) ) %>% 
    dplyr::select(dxdengue0, test_res ) %>% 
    filter(test_res %in% c("0", "1")) %>% 
    dplyr::mutate_all(as.factor) %>%
    confusionMatrix(data = .$test_res, reference = .$dxdengue0, positive = "1")
  
  
  
  
  test_mts <- left_join(test, only_symp %>% 
                          select(comb, score_pos, score_neg)) %>% 
    dplyr::mutate(test_res = case_when(
      score_pos > score_neg ~ "1",
      score_neg > score_pos ~ "0",
      TRUE ~ "2"
    ) ) %>% 
    dplyr::select(dxdengue0, test_res ) %>% 
    filter(test_res %in% c("0", "1")) %>% 
    dplyr::mutate_all(as.factor) %>%
    confusionMatrix(data = .$test_res, reference = .$dxdengue0, positive = "1")
  

    ret <- bind_rows(  c(train_mts$overall, train_mts$byClass) %>% round(., 3), 
                c(test_mts$overall, test_mts$byClass) %>% round(., 3)) %>% 
      mutate(type = c("train","test"))
    return(ret)
})

### Variables del hemograma

calc_params <- function(x, ...){
  if(length(x[!is.na(x)]) > 1){
    m <- mean(x, ...  )
    desv <- sd(x, ...)
  }else{
    m <- mean(x, ...)
    desv <- m*5
  }
ret <- paste0(m,";", desv)

  return(ret)
}

only_hemo_vars <- purrr::map( .x = list("1", "0") , .f = function(.x ){
  nms <- ifelse(.x == "1", "pos", 'neg')
  
  in_df %>% 
    dplyr::select(comb, dxdengue0:ncol(.)) %>% 
    mutate(dxdengue0 = as.character(dxdengue0), id = as.character(id)) %>% 
    filter(dxdengue0 == .x) %>% 
    group_by(comb) %>% 
    summarise_if(is.numeric, calc_params , na.rm = T) %>% 
    rename_if(is.character, .funs = function(i){
      ifelse(!grepl("comb", i), paste0(i,"_hemo_", nms), i)
    })
  
}) %>% purrr::reduce(., full_join)


 arbol_hemo_vars <- arbol_only_symp %>% 
   full_join(only_hemo_vars) %>% 
   tibble
  

 
 #validacion algoritmo con variables del hemograma
 
 to_val <- in_df %>% 
   left_join(., arbol_hemo_vars) %>%
   mutate_if(.predicate = function(i){ length(table(i))==2 }, as.character) %>% 
   dplyr::select(comb, score_pos, score_neg,  
          !!hemo_vars, 
          !!paste0(hemo_vars,"_hemo_pos"),
          !!paste0(hemo_vars,"_hemo_neg"),
          dxdengue0)
 
 
 scores <- lapply(hemo_vars, function(i){
 x <- to_val %>% 
     dplyr::select(matches(i)) %>%
     apply(., 1, function(k){
       
       q <-  gsub(" ", "", k[1]) %>% as.double()
       
       if(!is.na(q)){
        
         if(!is.na(k[2]) & !grepl(pattern = "NaN", k[2])){
           params_dist_pos <- str_split(k[2], ";", simplify = T)
          
           m_pos <- as.numeric(params_dist_pos[1])
           sd_pos <-  as.numeric(params_dist_pos[2])
           
           #mx_dens <- dnorm(m_pos, m_pos, sd_pos) 
           
           res_pos <- dnorm(q, m_pos, sd_pos, log = F)/dnorm(m_pos, m_pos, sd_pos, log = F)
         
         }else if(is.na(k[2])){
           params <- global_parms %>% 
             filter(dxdengue0 == 1) %>% 
             dplyr::select(matches(i))
           
           res_pos <- dnorm(q, as.numeric(params[1]), as.numeric(params[2]))/dnorm(as.numeric(params[1]), as.numeric(params[1]), as.numeric(params[2]))
           res_pos <- ifelse(res_pos > 1, 1, res_pos)
           
         }else{
           res_pos<- 0.5
         }
         
         
         if(!is.na(k[3]) & !grepl(pattern = "NaN", k[3])){
           params_dist_neg <- str_split(k[3], ";", simplify = T)
           
           m_neg <- as.numeric(params_dist_neg[1])
           sd_neg <-  as.numeric(params_dist_neg[2])
           
           #mx_dens <- dnorm(m_neg, m_neg, sd_neg) 
           
           res_neg <- dnorm(q, m_neg, sd_neg)/dnorm(m_neg, m_neg, sd_neg)
           
          
         }else if(is.na(k[3])){
           params <- global_parms %>% 
             filter(dxdengue0 == 0) %>% 
             dplyr::select(matches(i))
           
           res_neg <- dnorm(q, as.numeric(params[1]), as.numeric(params[2]))/dnorm(as.numeric(params[1]), as.numeric(params[1]), as.numeric(params[2]))
           res_neg <- ifelse(res_neg > 1, 1, res_neg)
           
         }else{
           res_neg<- 0.5
         }
         
         
       }else{
         res_pos <- 0.5
         res_neg<- 0.5
       }
       
       
       return(data.frame( res_pos, res_neg))
     }) %>% bind_rows()
 names(x) <- c(paste0(i, "_score_pos"), paste0(i, "_score_neg"))
 return(x)
 }) %>% 
   bind_cols()%>% 
   tibble
   
# names(scores) <- paste0(hemo_vars, "_score")
 
 to_val <- scores %>% 
   dplyr::select(matches("_pos")) %>%
   mutate(score_hemo_pos = apply(., 1, prod)) %>% 
   bind_cols(to_val, .)
 
 to_val <- scores %>% 
   dplyr::select(matches("_neg")) %>% 
   mutate(score_hemo_neg = apply(., 1, prod, na.rm = T)) %>% 
   bind_cols(to_val, .)
 
 
 to_val <- to_val %>% 
   mutate(score_final_pos  = score_hemo_pos * score_pos,
          score_final_neg  = score_hemo_neg * score_neg,
          test_res         = ifelse(score_final_pos > score_final_neg, "1", "0"),
          test_res         = as.factor(test_res),
          dxdengue0        = as.factor(dxdengue0))  
   
 caret::confusionMatrix(data = to_val$test_res, reference = to_val$dxdengue0, positive = "1")

 k <-5
 
 cv_hemo <- lapply(1:k, function(i){
   
   test <- in_df %>% 
     sample_n(size = round(0.2*nrow(.)))
   
   train <- in_df %>%
     filter(!id %in% test$id) 
   
   only_symp <- train %>%
     dplyr::select(-id) %>% 
     group_by(comb) %>% 
     dplyr::summarise(., dx_pos = sum(dxdengue0 == 1),
                      dx_neg = sum(dxdengue0 == 0),
                      total  = n()) %>% 
     ungroup() %>% 
     full_join(comb, .) %>%
     dplyr::mutate_all(.funs = function(i){ifelse(is.na(i), 0, i)}) %>% 
     dplyr::mutate(score_pos = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
       mean(rbeta(1000, 1+x, 1+y))
     }),
     score_neg = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
       mean(rbeta(1000, 1+y, 1+x))
     }),
     test_res = case_when(
       score_pos > score_neg ~ "1",
       score_neg > score_pos ~ "0",
       TRUE ~ "2"
     ) )
   
   
   only_hemo_vars <- purrr::map( .x = list("1", "0") , .f = function(.x ){
     nms <- ifelse(.x == "1", "pos", 'neg')
     
     train %>% 
       dplyr::select(comb, dxdengue0:ncol(.)) %>% 
       mutate(dxdengue0 = as.character(dxdengue0), id = as.character(id)) %>% 
       filter(dxdengue0 == .x) %>% 
       group_by(comb) %>% 
       summarise_if(is.numeric, calc_params , na.rm = T) %>% 
       rename_if(is.character, .funs = function(i){
         ifelse(!grepl("comb", i), paste0(i,"_hemo_", nms), i)
       })
     
   }) %>% purrr::reduce(., full_join)
   
   
   arbol_hemo_vars <- only_symp %>% 
     full_join(only_hemo_vars) %>% 
     tibble
   
   to_val <- test %>% 
     left_join(., arbol_hemo_vars) %>%
     mutate_if(.predicate = function(i){ length(table(i))==2 }, as.character) %>% 
     dplyr::select(comb, score_pos, score_neg,  
            !!hemo_vars, 
            !!paste0(hemo_vars,"_hemo_pos"),
            !!paste0(hemo_vars,"_hemo_neg"),
            dxdengue0)
   
   
   scores <- lapply(hemo_vars, function(i){
     x <- to_val %>% 
       dplyr::select(matches(i)) %>%
       apply(., 1, function(k){
         
         q <-  gsub(" ", "", k[1]) %>% as.double()
         
         if(!is.na(q)){
           
           if(!is.na(k[2]) & !grepl(pattern = "NaN", k[2])){
             params_dist_pos <- str_split(k[2], ";", simplify = T)
             
             m_pos <- as.numeric(params_dist_pos[1])
             sd_pos <-  as.numeric(params_dist_pos[2])
             
             #mx_dens <- dnorm(m_pos, m_pos, sd_pos) 
             
             res_pos <- dnorm(q, m_pos, sd_pos, log = F)/dnorm(m_pos, m_pos, sd_pos, log = F)
             
           }else if(is.na(k[2])){
             params <- global_parms %>% 
               filter(dxdengue0 == 1) %>% 
               dplyr::select(matches(i))
             
             res_pos <- dnorm(q, as.numeric(params[1]), as.numeric(params[2]))/dnorm(as.numeric(params[1]), as.numeric(params[1]), as.numeric(params[2]))
             res_pos <- ifelse(res_pos > 1, 1, res_pos)
             
           }else{
             res_pos<- 0.5
           }
           
           
           if(!is.na(k[3]) & !grepl(pattern = "NaN", k[3])){
             params_dist_neg <- str_split(k[3], ";", simplify = T)
             
             m_neg <- as.numeric(params_dist_neg[1])
             sd_neg <-  as.numeric(params_dist_neg[2])
             
             #mx_dens <- dnorm(m_neg, m_neg, sd_neg) 
             
             res_neg <- dnorm(q, m_neg, sd_neg)/dnorm(m_neg, m_neg, sd_neg)
             
             
           }else if(is.na(k[3])){
             params <- global_parms %>% 
               filter(dxdengue0 == 0) %>% 
               dplyr::select(matches(i))
             
             res_neg <- dnorm(q, as.numeric(params[1]), as.numeric(params[2]))/dnorm(as.numeric(params[1]), as.numeric(params[1]), as.numeric(params[2]))
             res_neg <- ifelse(res_neg > 1, 1, res_neg)
             
           }else{
             res_neg<- 0.5
           }
           
           
         }else{
           res_pos <- 0.5
           res_neg<- 0.5
         }
         
         
         return(data.frame( res_pos, res_neg))
       }) %>% bind_rows()
     names(x) <- c(paste0(i, "_score_pos"), paste0(i, "_score_neg"))
     return(x)
   }) %>% 
     bind_cols()%>% 
     tibble
   
   
     to_val <- scores %>% 
     dplyr::select(matches("_pos")) %>% 
     mutate(score_hemo_pos = apply(., 1, prod)) %>% 
     bind_cols(to_val, .)
   
   to_val <- scores %>% 
     dplyr::select(matches("_neg")) %>% 
     mutate(score_hemo_neg = apply(., 1, prod)) %>% 
     bind_cols(to_val, .)
   
   
   to_val <- to_val %>% 
     mutate(score_final_pos  = score_hemo_pos * score_pos,
            score_final_neg  = score_hemo_neg * score_neg,
            test_res         = ifelse(score_final_pos > score_final_neg, "1", "0"),
            test_res         = as.factor(test_res),
            dxdengue0        = as.factor(dxdengue0))  
   
  mtrs <- caret::confusionMatrix(data = to_val$test_res, reference = to_val$dxdengue0, positive = "1")
  
  ret <- bind_rows(  c(mtrs$overall, mtrs$byClass) %>% round(., 3)) 
   return(ret)
 })
 
 
 ##### otras coosas por ahi
 
 only_hemo_vars <- in_df %>% 
   select(comb, dxdengue0:ncol(.)) %>% 
   mutate(dxdengue0 = as.character(dxdengue0), id = as.character(id)) %>% 
   #filter(dxdengue0 == .x) %>% 
   group_by(comb) %>% 
   summarise_if(is.numeric, calc_params , na.rm = T) %>% 
   rename_if(is.character, .funs = function(i){
     ifelse(!grepl("comb", i), paste0(i,"_dist"), i)
   })
 

 
 

n <- 40
k <- 30
post_val<- c()
post_pred <- c()
post_pred_dens <- c()
for(i in 1:10000){
post_val[i] <- rbeta(1, 1+k, 1+n-k)
post_pred[i] <- rbinom(1, n, p = post_val[i])
post_pred_dens[i] <- pbinom(post_pred[i], n, p = post_val[i])
}


hist(post_pred)
hist(post_pred_dens)

mean(post_pred)/40
median(post_pred)/40
mean(rbeta(1000, 1+k, 1+n-k))



df <- data.frame(dx = c(rep(1, 30), rep(0, 10)),  )

n < 40
k <-10

Nb = 2000
niter <- 10000

tau <- rep(NA, niter+1)
tau[1] <- 0.75

for(i in 2:niter){
  prop <- tau[i-1]
  x <- rbinom(1, k, prop)
  y <- rbinom(1, n-k, 1-prop)
  tau[i] <- rbeta(1, 1+x, 1+n-x)
    
}

exp(log(0.5)+dbeta(x = 0.75, 30, 10, log = F))

0.5*(10/40)

hist(tau[Nb:niter])

abline(v = mean(tau[Nb:niter]))
abline(v = median(tau[Nb:niter]), col = "red")


#adicionanado una normal

n < 40
k <-30

Nb = 2000
niter <- 10000

tau <- rep(NA, niter+1)
mu1 <- rep(NA, niter+1)
tau[1] <- 0.5
mu1[1] <- 500
for(i in 2:niter){
  prop <- runif(1)
  x <- rbinom(1, k, prop)
  y <- rbinom(1, n-k, 1-prop)
  mu1 <- rnorm(1, mu1[i-1], 10000)
  tau[i] <- dbeta(x = rbeta(1, 1+x, 1+n-x ), 1+x, 1+n-x, log = TRUE ) +
    dnorm(mu1[i], 4200, 3500, log =TRUE)
  
}

hist(tau[Nb:niter])

abline(v = mean(tau[Nb:niter]))
abline(v = median(tau[Nb:niter]), col = "red")


f <- function(theta, data){
  
  x <- theta[1]
  #x^(z)*(1-x)^(n-z) %>% log
  pos <- dbeta(x, 1+data$k, 1+data$n-data$k, log = TRUE)
  return(pos)
}

start <- c(0.5)
data <- list(n = 40, k = 30, prop = k/n)

fn <- LearnBayes::gibbs(logpost = f , start = c(0.5), data , m = 10000, scale = 0.3 )

mean(fn$par, na.rm = T)
hist(fn$par)


logpost = f
start = c(10, 10, 0.02)
p = length(start)
m = 10000
vth = array(0, dim = c(m, p))
f0 = logpost(start, data)
arate = array(0, dim = c(1, p))
th0 = start
for (i in 1:m) {
  for (j in 1:p) {
    th1 = th0
    th1[j] = th0[j] + rnorm(1) * scale[j]
    f1 = logpost(th1, data)
    u = runif(1) < exp(f1 - f0)
    th0[j] = th1[j] * (u == 1) + th0[j] * (u == 0)
    f0 = f1 * (u == 1) + f0 * (u == 0)
    vth[i, j] = th0[j]
    arate[j] = arate[j] + u
  }
}
arate = arate/m
stuff = list(par = vth, accept = arate)

mean(stuff$par, na.rm = T)
hist(stuff$par)

acf(vth)
plot(vth[,1], type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1)
plot(vth[,2], type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1)
plot(vth[,3], type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1)



