setwd(Sys.getenv("PROJECT_DIR", unset=getwd()))
source("ROG_cases10_11_observation.R")
env_num <- function(x,d) {z<-Sys.getenv(x,unset=""); if(nzchar(z)) as.numeric(z) else d}
model <- toupper(Sys.getenv("MODEL",unset="RI")); vb <- if(model=="RS") env_num("VAR_B",0.1) else 0
methods <- strsplit(Sys.getenv("METHODS",unset="ORACLE,OBS,LM,KM,GMM,CPF,MIXREG,CIRG"),",")[[1]]
ans <- run_case10_observation(nloop=as.integer(env_num("NLOOP",50)), p=as.integer(env_num("P",50)),
                              N_test=as.integer(env_num("N_TEST",2500)), K_type=as.integer(env_num("K_TYPE",6)),
                              modes_per_type=as.integer(env_num("MODES_PER_TYPE",2)), Var.a=env_num("VAR_A",2.25),
                              Var.b=vb, Var.e=env_num("VAR_E",9), tau=as.integer(env_num("TAU",51)),
                              sa_max_iter=as.integer(env_num("SA_MAX_ITER",80)), methods=methods)
print(ans$summary)
dir.create("results_observation",showWarnings=FALSE,recursive=TRUE)
saveRDS(ans,file=sprintf("results_observation/case10_%s.rds",model))
