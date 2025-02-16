
#######################################################
#######      Resumenes descriptivos varios ############
#######################################################

asimetria=function(x) {
  m3=mean((x-mean(x))^3)
  skew=m3/(sd(x)^3)
  skew}

#### obtencion del coeficiente de curtosis muestral

kurtosis=function(x) {
  m4=mean((x-mean(x))^4)
  kurt=m4/(sd(x)^4)
  kurt}

###########################################
#### Scatterplot con Histogramas paralelos

scatterhist = function(x, y, xlab="", ylab=""){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  plot(x,y)
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0,
        at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0,
        at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}

##########################################################
#### Función para Gráfico de perfiles para cada variable

makeProfilePlot <- function(mylist,names)
{
  require(RColorBrewer)
  # find out how many variables we want to include
  numvariables <- length(mylist)
  # choose 'numvariables' random colours
  colours <- brewer.pal(numvariables,"Set1")
  # find out the minimum and maximum values of the variables:
  mymin <- 1e+20
  mymax <- 1e-20
  for (i in 1:numvariables)
  {
    vectori <- mylist[[i]]
    mini <- min(vectori)
    maxi <- max(vectori)
    if (mini < mymin) { mymin <- mini }
    if (maxi > mymax) { mymax <- maxi }
  }
  # plot the variables
  for (i in 1:numvariables)
  {
    vectori <- mylist[[i]]
    namei <- names[i]
    colouri <- colours[i]
    if (i == 1) { plot(vectori,col=colouri,type="l",ylim=c(mymin,mymax)) }
    else         { points(vectori, col=colouri,type="l")                                     }
    lastxval <- length(vectori)
    lastyval <- vectori[length(vectori)]
    text((lastxval-10),(lastyval),namei,col="black",cex=0.6)
  }
}


##############################################################
#### Setup of a Correlation Lower Panel in Scatterplot Matrix

myPanel.hist <- function(x, ...){
  usr <- par("usr")
  on.exit(par(usr))
  # Para definir región de graficiación
  par(usr = c(usr[1:2], 0, 1.5) )
  # Para obtener una lista que guarde las marcas de clase y conteos en cada una:
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks;
  nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  # Para dibujar los histogramas
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

###############################################################
#### Setup of a Boxplot Diagonal Panel in Scatterplot Matrix

myPanel.box <- function(x, ...){
  usr <- par("usr", bty = 'n')
  on.exit(par(usr))
  par(usr = c(-1, 1, min(x) - 0.5, max(x) + 0.5))
  b <- boxplot(x, plot = F)
  whisker.i <- b$stats[1,]
  whisker.s <- b$stats[5,]
  hinge.i <- b$stats[2,]
  mediana <- b$stats[3,]
  hinge.s <- b$stats[4,]
  rect(-0.5, hinge.i, 0.5, mediana, col = 'gray')
  segments(0, hinge.i, 0, whisker.i, lty = 2)
  segments(-0.1, whisker.i, 0.1, whisker.i)
  rect(-0.5, mediana, 0.5, hinge.s, col = 'gray')
  segments(0, hinge.s, 0, whisker.s, lty = 2)
  segments(-0.1, whisker.s, 0.1, whisker.s)
}


##############################################################
#### Setup of a Correlation Lower Panel in Scatterplot Matrix

myPanel.cor <- function(x, y, digits = 2, prefix = "", cex.cor){
  usr <- par("usr")
  on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if(missing(cex.cor))
    cex = 0.4/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1 + 1.5*abs(r))
}

# QQ-plot with Shapiro-Wilk normal test
QQnorm <- function(datos){
    lab.plot <- "Normal Q-Q Plot of Datos Crudos"
  shapiro <- shapiro.test(datos)
  shapvalue <- ifelse(shapiro$p.value < 0.001,
      "P value < 0.001", paste("P value = ",
        round(shapiro$p.value, 4), sep = ""))
  shapstat <- paste("W = ", round(shapiro$statistic, 4),
                    sep = "")
  q <- qqnorm(datos, plot.it = FALSE)
  qqnorm(datos, main = lab.plot)
  qqline(datos, lty = 1, col = 2)
  text(min(q$x, na.rm = TRUE), max(q$y,
      na.rm = TRUE)*0.95, pos = 4,
      'Shapiro-Wilk Test', col = "blue", font = 2)
  text(min(q$x, na.rm = TRUE), max(q$y,
      na.rm = TRUE)*0.80, pos = 4, shapstat,
      col = "blue", font = 3)
  text(min(q$x, na.rm = TRUE), max(q$y, na.rm = TRUE)*0.65,
       pos = 4, shapvalue, col = "blue", font = 3)
}

###################################################################
#### QQ-plot with Shapiro-Wilk normal test (datos transformados)

QQnorm_transf <- function(datos){
  lab.plot <- "Normal Q-Q Plot of Datos Transformados"
  shapiro <- shapiro.test(datos)

  shapvalue <- ifelse(shapiro$p.value < 0.001,
          "P value < 0.001", paste("P value = ",
          round(shapiro$p.value, 4), sep = ""))

  shapstat <- paste("W = ", round(shapiro$statistic, 4),
                    sep = "")

  q <- qqnorm(datos, plot.it = FALSE)
  qqnorm(datos, main = lab.plot)
  qqline(datos, lty = 2, col = 2)

  text(min(q$x, na.rm = TRUE),
       max(q$y-0.2, na.rm = TRUE)*0.95, pos = 4,
       'Shapiro-Wilk Test', col = "blue", font = 2)

  text(min(q$x, na.rm = TRUE),
       max(q$y-0.7, na.rm = TRUE)*0.80, pos = 4,
       shapstat, col = "blue", font = 3)

  text(min(q$x, na.rm = TRUE),
       max(q$y-1.5, na.rm = TRUE)*0.65, pos = 4,
       shapvalue, col = "blue", font = 3)
}


##############################################
###### coeficiente de asimetría

asimetria=function(x) {
  m3=mean((x-mean(x))^3)
  skew=m3/(sd(x)^3)
  skew}

#############################################
###### Coeficiente de Kurtosis
kurtosis=function(x) {
  m4=mean((x-mean(x))^4)
  kurt=m4/(sd(x)^4)
  kurt}


######################################################
###### Función para Resumen descriptivo por grupos

resumen_xgrupos <- function(misdatos,grupos)
{
  # se hallan los nombres de las variables
  nombres_misdatos <- c(names(grupos),names(as.data.frame(misdatos)))
  # se halla la media dentro de cada grupo
  grupos <- grupos[,1] # nos aseguramos de que la var grupos no sea una lista
  medias <- aggregate(as.matrix(misdatos) ~ grupos, FUN = mean)
  names(medias) <- nombres_misdatos
  # se hallan las desv-estandar dentro de cada grupos:
  sds <- aggregate(as.matrix(misdatos) ~ grupos, FUN = sd)
  names(sds) <- nombres_misdatos
  # se hallan las varianzas dentro de cada grupos:
  varianzas <- aggregate(as.matrix(misdatos) ~ grupos, FUN = var)
  names(varianzas) <- nombres_misdatos
  # se hallan las medianas dentro de cada grupos:
  medianas <- aggregate(as.matrix(misdatos) ~ grupos, FUN = median)
  names(medianas) <- nombres_misdatos
  # se hallan los tama?os muestrales de cada grupo:
  tamanos_n <- aggregate(as.matrix(misdatos) ~ grupos, FUN = length)
  names(tamanos_n) <- nombres_misdatos
  list(Medias=medias,Desviaciones_Estandar=sds,
       Varianzas=varianzas, Medianas=medianas,
       Tamanos_Muestrales=tamanos_n)
}

##################################################
##################################################
#########           PH-AM               ##########
##################################################
##################################################

##################################################
#### PH Sobre Igualdad de Matrices de Var-Cov ####
##################################################

########################################################
### Función creada para la Prueba M-Box de Matrices 
### de Var-Cov, ie. para Sigam_1=SIgma_2, pob. Normal

prueba_M_Box2=function(x,y,alfa){
  g<-2
  n=nrow(x);m=nrow(y);p=ncol(x)
  s1=var(x);s2=var(y)
  v<-n+m-2
  sp<-( (n-1)*s1+(m-1)*s2 )/v
  M<-v*log( det(sp) )-( (n-1)*log( det(s1) ) + (m-1)*log( det(s2) ) )
  u<-( ( 1/(n-1) ) + ( 1/(m-1) ) - (1/v) )*( (2*p^2 + 3*p - 1)/(6*(p+1)*(g-1)) )
  c<-(1-u)*M
  df=( p*(p+1)*(g-1) )/2  # Grados de liber del num de la chi-cuadrado
  chi_tabla=qchisq(1-alfa,df) # Valor crítico de la chi o Chi-de la tabla
  valor_p=1-pchisq(c,df)  # valor-p de la prueba
  resultados=data.frame(M=M,U=u,C=c,df=df,Chi_Tabla=chi_tabla,Valor_p=valor_p)
  format(resultados, digits = 6)
}


##################################################
#### PH Sobre Igualdad de Vectores de Medias  ####
##################################################

########################################################
#### Función creada para la prueba de igualdad de medias, 
#### ie. para: mu_1-mu_2=mu_0, sigmas iguales, pob. Normal

HT2_sigmas_iguales=function(x,y,mu_0,alfa){
  mux=apply(x,2,mean);muy=apply(y,2,mean)
  sx<-var(x);sy<-var(y)
  n=nrow(x);m=nrow(y);p=ncol(x)
  df1=p;df2<-n+m-p-1 # Grados de libertad del num y denom de la F
  sp<-( (n-1)*sx + (m-1)*sy )/(n+m-2)
  T_2<-( (n*m)/(n+m) )*t(mux-muy-mu_0)%*%solve(sp)%*%(mux-muy-mu_0)
  k<-( (n+m-2)*p )/(n+m-p-1)
  F0<-(1/k)*T_2               # Estad?stica F_0=(1/k)T2
  F_tabla=qf(1-alfa,df1,df2)  # Valor cr?tico de la F o F-de la Tabla
  valor_p=1-pf(F0,df1,df2)    # valor-p de la prueba
  resultados<-data.frame(T2=T_2,k=k,F0=F0,
                   df1=df1,df2=df2,F_Tabla=F_tabla,Valor_p=valor_p)
  cat("El vector mu0 es:", mu_0 )
  format(resultados, digits = 6)
}

##########################################################
#### Función creada para la prueba de igualdad de medias, 
#### ie. para: mu_1-mu_2=mu_0, sigmas iguales, n-grande

HT2_sigmas_iguales_ngrande=function(x,y,mu_0,alfa){
  mux=apply(x,2,mean);muy=apply(y,2,mean)
  sx<-var(x);sy<-var(y)
  n=nrow(x);m=nrow(y);p=ncol(x)
  df=p      # Grados de libertad de la chi-cuadrado
  sp<-( (n-1)*sx + (m-1)*sy )/(n+m-2)
  chi_2<-( (n*m)/(n+m) )*t(mux-muy-mu_0)%*%solve(sp)%*%(mux-muy-mu_0)
  chi_tabla=qchisq(1-alfa,df)  # Valor de la chi_cuadrado, ie. chi_Tabla
  valor_p=1-pchisq(chi_2,df)    # valor-p de la prueba
  resultados<-data.frame(Chi2=chi_2,df=df,
                   Chi_Tabla=chi_tabla,Valor_p=valor_p)
  cat("El vector mu0 es:", mu_0 )
  format(resultados, digits = 6)
}


###########################################################
#### Función para PH de mu_x-mu_y=mu_0, sigmas diferentes 
#### y desconocidas, Poba. Normal
#### Aproximación de: Nel y Van Der Merwe (1986) 
#### para el valro de v

HT2_sigmas_diferentes=function(x,y,mu_0,alfa){
  mux=apply(x,2,mean);muy=apply(y,2,mean)
  sx<-var(x);sy<-var(y)
  n=nrow(x);m=nrow(y);p=ncol(x)
  v1<-(1/n)*sx;v2<-(1/m)*sy
  se<-v1+v2
  v<-( sum(diag(se%*%se)) +
         sum(diag(se))^2 )/( (1/(n-1))*(sum(diag(v1%*%v1)) +
           sum(diag(v1))^2) +
           ( 1/(m-1) )*(sum(diag(v2%*%v2)) +
           sum(diag(v2))^2) )
  v<-ceiling(v)
  df1=p;df2<-v-p+1  # Grados de libertad de la F
  sp<-( (n-1)*sx + (m-1)*sy )/(n+m-2)
  T_2<-t(mux-muy-mu_0)%*%solve(se)%*%(mux-muy-mu_0)
  k<-(v*p)/(v-p+1)
  F0<-(1/k)*T_2
  F_tabla=qf(1-alfa,df1,df2)
  valor_p=1-pf(F0,df1,df2)
  resultados=data.frame(T_2=T_2,v=v,k=k,F0=F0,
                  df1=df1,df2=df2,F_Tabla=F_tabla,Valor_p=valor_p)
  cat("El vector mu0 es:", mu_0 )
  format(resultados, digits = 6)
}


##########################################################
#### Función para PH de mu_x-mu_y=mu_0, sigmas diferentes 
#### y desconocidas, Poba. Normal
#### Aproximación de Krishnamoorty and Yu (2004)
#### texto-Guía con: p+p^2 en el numerador de v

HT2_sigmas_diferentes_texto_guia=function(x,y,mu_0,alfa){
  mux=apply(x,2,mean);muy=apply(y,2,mean)
  sx<-var(x);sy<-var(y)
  n=nrow(x);m=nrow(y);p=ncol(x)
  v1<-(1/n)*sx;v2<-(1/m)*sy
  se<-v1+v2
  numer<-p+(p^2)
  den1<-sum( diag( (v1%*%solve(se))%*%(v1%*%solve(se)) ) )
             + sum( ( diag( v1%*%solve(se) ) )^2 )
  den2<-sum( diag( (v2%*%solve(se))%*%(v2%*%solve(se)) ) )
             + sum( ( diag( v2%*%solve(se) ) )^2 )
  v<-(numer)/( den1/n + den2/m )
v<-ceiling(v)
  df1=p;df2<-v-p+1  # Grados de libertad de la F
  #sp<-( (n-1)*sx + (m-1)*sy )/(n+m-2)
  T_2<-t(mux-muy-mu_0)%*%solve(se)%*%(mux-muy-mu_0)
  k<-(v*p)/(v-p+1)
  F0<-(1/k)*T_2
  F_tabla=qf(1-alfa,df1,df2)
  valor_p=1-pf(F0,df1,df2)
  resultados=data.frame(T_2=T_2,v=v,k=k,F0=F0,
                  df1=df1,df2=df2,F_Tabla=F_tabla,Valor_p=valor_p)
  cat("El vector mu0 es:", mu_0 )
  format(resultados, digits = 6)
}


####################################################
#### PH sobre el vector de medias Poblacional  #####
####################################################

####################################################
#### Función para la PH de mu=mu_0-pob. Normal

HT2_mu0=function(x,mu_0,alfa){
  mu=apply(x,2,mean);s=var(x)
#  mu <- as.vector(mu)
  n=nrow(x);p=ncol(x)
  df1=p;df2=n-p
  T2<-n*t(mu-mu_0)%*%solve(s)%*%(mu-mu_0)
  k<-( (n-1)*p )/(n-p)
  F0<-(1/k)*T2
  F_tabla=qf(1-alfa,df1,df2)
  valor_p=1-pf(F0,df1,df2)
  resultados=data.frame(T2=T2,K=k,F0=F0,df1=df1,df2=df2,
                  F_Tabla=F_tabla,Valor_p=valor_p)
  cat("El vector mu0 es:", mu_0 )
  format(resultados, digits = 6)
}


##################################################
#### Función para la PH de mu=mu_0-n-grande

HT2_mu0_ngrande=function(x,mu_0,alfa){
  mu=apply(x,2,mean);s=var(x)
  n=nrow(x);p=ncol(x)
  df=p
  chi_2<-n*t(mu-mu_0)%*%solve(s)%*%(mu-mu_0)
  chi_tabla=qchisq(1-alfa,df)
  valor_p=1-pchisq(chi_2,df)
   resultados=data.frame(Chi_2=chi_2,df=df,Chi_Tabla=chi_tabla,
                         Valor_p=valor_p)
  format(resultados, digits = 6)
}


####################################################
#######       PH sobre Constrastes    ##############
####################################################

#########################################################
##### Función Creada para la PH de: CU=mu_0-Pob. Normal

HT2_CU=function(x,C,delta_0,alfa){
  mu=as.vector(apply(x,2,mean));s=var(x)
  n=nrow(x);p=ncol(x)
  k<-nrow(C)        ## n?mero de contrastes
  df1=k
  df2=n-k
  T2<-n*t(C%*%mu-delta_0)%*%solve(C%*%s%*%t(C))%*%(C%*%mu-delta_0)
  c<-( (n-1)*k )/(n-k)
  F0<-(1/c)*T2
  F_tabla=qf(1-alfa,df1,df2)
  valor_p=1-pf(F0,df1,df2)
  resultados=data.frame(T2=T2,c=c,F0=F0,df1=df1,df2=df2,
                  F_Tabla=F_tabla,Valor_p=valor_p)
  cat("El vector mu0 es:", delta_0 )
  format(resultados, digits = 6)
}


###########################################################
#### Función Creada para la PH de: CU=mu_0,  n-Grande

HT2_CU_ngrande=function(x,C,delta_0,alfa){
  mu=as.vector(apply(x,2,mean));s=var(x)
  n=nrow(x);p=ncol(x)
  k<-nrow(C)
  df1=k
  chi2<-n*t(C%*%mu-delta_0)%*%solve(C%*%s%*%t(C))%*%(C%*%mu-delta_0)
  chi_tabla=qchisq(1-alfa,df1)
  valor_p=1-pchisq(chi2,df1)
  resultados=data.frame(Chi2=chi2,df1=df1,
                  Chi_Tabla=chi_tabla,Valor_p=valor_p)
  cat("El vector mu0 es:", delta_0 )
  format(resultados, digits = 6)
}


####################################################
#### PH sobre Estructura de Covarianzas    #########
####################################################

########################################################
#### Función para la PH de Razón de Ver. una Matriz de 
#### Var-Cov: ie. Sigma=Sigma_0, n-grande

sigma_sigma0_ngrande=function(x,Sigma_0,alfa){
  x=as.matrix(x)
  Sigma=as.matrix(Sigma_0)
  p=ncol(x);n=nrow(x)
  S=var(x)
  ## Construcción del Estadístico de Prueba
  mesa=S%*%solve(Sigma_0)
  lamda_est= n*sum( diag(mesa) ) - n*log( det(S) ) +
    n*log( det(Sigma_0) ) - n*p
  df=0.5*p*(p+1)    ## grados de libertad de la chi-2
  chi_tabla=qchisq(1-alfa,df)
  valor_p=1-pchisq(lamda_est,df)
  result=data.frame(Landa_est = lamda_est,df=df,
              Chi_Tabla=chi_tabla,Valor_P=valor_p)
  format(result, digits = 6)
}


###########################################################
#### Función para la PH de Razón de Ver. una Matriz de 
#### Var-Cov: ie. Sigma=Sigma_0, n-pequeña

sigma_sigma0_npqna=function(x,Sigma_0,alfa){
  x=as.matrix(x)
  Sigma=as.matrix(Sigma_0)
  p=ncol(x);n=nrow(x)
  S=var(x)
  ## Construcción del Estadístico de Prueba
  mesa=S%*%solve(Sigma_0)
  lamda_est= n*sum( diag(mesa) ) - n*log( det(S) ) +
    n*log( det(Sigma_0) ) - n*p
  c<-1- ( 1/( 6*(n-1) ) )*( 2*p+1-( 2/(p+1) ) )
  lamda_1_est<-c*lamda_est
  df=0.5*p*(p+1)
  chi_tabla=qchisq(1-alfa,df)
  valor_p=1-pchisq(lamda_1_est,df)
  result=data.frame(Lamda1_est=lamda_1_est,c=c, df=df,
              Chi_Tabla=chi_tabla,Valor_P=valor_p)
  format(result, digits = 5)
}


####################################################
#######           MUESTRAS PAREADAS   ##############
####################################################


#############################################################
#### Función creada para la prueba de igualdad de medias, 
#### ie. para: mu_1-mu_2=mu_0, muestras-pareadas, pob. Normal

HT2_pareadas=function(x,y,alfa){
  mux=apply(x,2,mean);muy=apply(y,2,mean)
  sx<-var(x);sy<-var(y)
  n=nrow(x);m=nrow(y);p=ncol(x)
  
  x1 <- x[1:min(n,m),]
  y1 <- y[1:min(n,m),]
  datos_d <- x1-y1
  mu=apply(datos_d,2,mean);s=var(datos_d)
  n=nrow(datos_d);p=ncol(datos_d)
  df1=p;df2=n-p
  T2<-n*t(mu)%*%solve(s)%*%(mu)
  k<-( (n-1)*p )/(n-p)
  F0<-(1/k)*T2
  F_tabla=qf(1-alfa,df1,df2)
  valor_p=1-pf(F0,df1,df2)
  resultados=data.frame(T2=T2,K=k,F0=F0,df1=df1,df2=df2,
                        F_Tabla=F_tabla,Valor_p=valor_p)
  format(resultados, digits = 6)
}

############################################################
#### Función creada para la prueba de igualdad de medias, ie. 
#### para: mu_1-mu_2=mu_0, muestras-pareadas, n-grande


HT2_pareadas_ngrande=function(x,y,alfa){
  mux=apply(x,2,mean);muy=apply(y,2,mean)
  sx<-var(x);sy<-var(y)
  n=nrow(x);m=nrow(y);p=ncol(x)
  
  x1 <- x[1:min(n,m),]
  y1 <- y[1:min(n,m),]
  datos_d <- x1-y1
  mu=apply(datos_d,2,mean);s=var(datos_d)
  
  #  mu <- as.vector(mu)
  n=nrow(datos_d);p=ncol(datos_d)
  df=p
  
  chi_2<-n*t(mu)%*%solve(s)%*%(mu)
  chi_tabla=qchisq(1-alfa,df)
  valor_p=1-pchisq(chi_2,df)
  resultados=data.frame(Chi_2=chi_2,df=df,Chi_Tabla=chi_tabla,
                        Valor_p=valor_p)
  format(resultados, digits = 6)
}





