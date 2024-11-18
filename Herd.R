herd<-read_csv("input/Herd.csv", col_names = FALSE)
model1 <- glm(X2~X1, data=herd, family = binomial(logit))
#summary(model1)
a<-predict(model1, type = 'response')
b<-cbind(herd, a)
ggplot(b,aes(X1, X2)) +
  geom_point() + 
  stat_smooth(method="glm", 
              method.args = list(family="binomial"),
              se=FALSE) +
  theme_bw() +
  theme_classic() +
  xlab("Coverage of vaccine") +
  ylab("Risk reduction for those who are not vaccinated")
ggsave("herd.pdf", width = 6, height = 4, dpi = 300)

#predp <- exp(model1[["coefficients"]][2]*b[,"X1"]+model1[["coefficients"]][1])/
#  (1+exp(model1[["coefficients"]][2]*b[,"X1"]+model1[["coefficients"]][1]))
#c<-cbind(b,predp)