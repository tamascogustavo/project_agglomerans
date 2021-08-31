library(ggfortify)
library(cluster)

df <- read.csv("/Users/gustavotamasco/project_agglomerans/pantoea_arg.csv")



#Fix variance col = 0

#https://stackoverflow.com/questions/40315227/how-to-solve-prcomp-default-cannot-rescale-a-constant-zero-column-to-unit-var

which(apply(df, 2, var)==0)
df_new<-df[ , which(apply(df, 2, var) != 0)]
row.names(df_new) <- df$Isolates
#new_df <- data.matrix(df, rownames.force = NA)

pca_res <- prcomp(df_new, scale. = TRUE)


summary(pca_res)

autoplot(pca_res, data = df, colours("Isolates"))

autoplot(pam(df_new,4), frame = TRUE, frame.type = "norm", label = TRUE, label.size = 2, loadings = TRUE, loadings.label = TRUE, loadings.label.size  = 4)

autoplot(clara(df_new,3), label = FALSE, label.size = 2, loadings = TRUE, loadings.label = TRUE, loadings.label.size  = 4)
test <-clara(df_new,3)
test$clustering

set.seed(1)
autoplot(kmeans(df_new,4), data = df_new, label = TRUE, label.size = 2, loadings = TRUE, loadings.label = TRUE, loadings.label.size  = 4)
