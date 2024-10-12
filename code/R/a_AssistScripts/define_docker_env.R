# Author: Anna Lyubetskaya. Date: 23-04-24
# Capture R environment information

docker_v1 = list(installed.packages(), 
                 BiocManager::repositories(), 
                 sessionInfo(), 
                 "TBio2022-12")

names(docker_v1) = c("installed_packages", "biocm_repositories", "session_info", "starter_env")

saveRDS(docker_v1, "XXXX")
