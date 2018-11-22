############################################################################################
######Olmo, Trabajo Final######IDENTIFICACION DE PERIODOS CALIDOS / FRIOS###################
############################################################################################
rm(list=ls()) ; cat("\014") ; graphics.off()

WORKDIR_UBUNTU <- "/media/matias/Mati/Otros/2C_2018/LPIM/Trabajo_Final";WORKDIR_WINDOWS <- "E:/Otros/2C_2018/LPIM/Trabajo_Final"
setwd(WORKDIR_WINDOWS) #setea el directorio de trabajo
OUTPUTS <- "OUTPUTS/" #carpeta donde se iran guardando los resultados

#Carga las librerias necesarias
library(padr)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(scales)

#Se trabaja con un archivo .txt con los datos de Tx y Tn de la estacion Formosa Aero
#Nombre y clase para las columnas del .txt, luego abro el archivo
columns <- c("Fecha","Tx","Tn") ; class_columns <- c("Date","numeric","numeric")
Formosa_Aero <- read.table("Formosaaero.txt", sep=",",dec=".",col.names=columns,colClasses = class_columns)

######TRATAMIENTO BASICO DE LOS DATOS######
#Se comprueba que haya una fila por cada dia, es decir, que no haya dias salteados, con la liberia "padr"
#que requiere una columna con fechas y controla que no haya saltos en el intervalo de dias (en este caso)
corregidos <- pad(Formosa_Aero) #la funcion pad completa en caso de que existan saltos temporales

#Corrobora que ambos objetos tengan la misma dimension, si no hay saltos en los datos
if(dim(corregidos)[1] == dim(Formosa_Aero)[1]) {
  rm(corregidos,columns,class_columns)} else {
    Formosa_Aero <- corregidos ; rm(corregidos,columns,class_columns)}

#Codifica datos faltantes y conserva los datos de los anios completos 1963-2003 (elimina Nov62 y Dic62 y la ultima fila, unico dia del 2004)
#Esto lo realiza para trabajar con anios completos y no afectar la estadistica de cantidad de eventos por mes
Formosa_Aero[Formosa_Aero==999.9] <- NA
Formosa_Aero <- Formosa_Aero[-which(format(Formosa_Aero[,1],"%Y")==1962),]
Formosa_Aero <- Formosa_Aero[-dim(Formosa_Aero)[1],]

porc_na <- sum(is.na(Formosa_Aero))/dim(Formosa_Aero)[1]*100
print(paste("La base de datos posee",round(porc_na,digits=1),"% de datos faltantes"),quote=F)

#Es posible completar los datos faltantes con el valor medio de los dias contiguos a partir de la siguiente funcion:
completar_NA<-function(datos,n_columna){
  #Para el primer dato
  if (is.na(datos[1,n_columna])){datos[1,n_columna]<-datos[2,n_columna]} 
  #Para el ultimo dato
  if (is.na(datos[length(datos[,1]),n_columna])){datos[length(datos[,1]),n_columna]<-datos[(length(datos[,1])-1),n_columna]}
  #Para el resto de los datos
  for (i in 2:(length(datos[,1])-1)){
    if (is.na(datos[i,n_columna])){datos[i,n_columna]<-(datos[(i+1),n_columna]+datos[(i-1),n_columna])/2 }
  }
  #La funcion devuelve los datos ya completos
  datos[,n_columna]
}

Formosa_Aero$Tx <- completar_NA(Formosa_Aero,2) #Completa Tx
Formosa_Aero$Tn <- completar_NA(Formosa_Aero,3) #Completa Tn

######ESTIMACION DE PERIODOS CALIDOS / FRIOS######
#Selecciona que clase de eventos desea estudiar
periodo <- as.numeric(readline("¿Qué tipo de eventos desea analizar? Complete el numero: CALIDOS (1) o FRIOS (2):"))

repeat{
  if(periodo==1 | periodo==2)break
  periodo <- as.numeric(readline("Error. Intente nuevamente. Complete el número: CALIDOS (1) o FRIOS (2):"))}

if(periodo==1){periodo <- "calidos" ; sg <- 1} ; if(periodo==2){periodo <- "frios" ; sg <- -1}

#Calcula la temperatura media diaria, la media climatologica y el desvio estandar
Formosa_Aero$Tmean <- apply(Formosa_Aero[,2:3],1,mean)
Formosa_Aero$meses <- month(Formosa_Aero[,1])
medias <- aggregate(Tmean~meses,Formosa_Aero,mean,na.rm=T); desvios <- aggregate(Tmean~meses,Formosa_Aero,sd,na.rm=T)

#Repite estos valores para que tengan la misma longitud que los datos de Formosa_Aero y evalua si superaron
#a la media mensual mas/menos un desvio estandar. 
for(i in 1:12){
  Formosa_Aero$Tmsd[Formosa_Aero$meses==i] <- round(medias$Tmean[i] + sg*desvios$Tmean[i],digits = 1)
}

Formosa_Aero$diff <- sg*(Formosa_Aero$Tmean - Formosa_Aero$Tmsd)
Formosa_Aero$diff[is.na(Formosa_Aero$diff)] <- 0 #renombra los NA como ceros para continuar

Formosa_Aero$sg <- NA  #primero crea una nueva columna que se ira llenando con el signo de la diferencia
Formosa_Aero$sg[which(Formosa_Aero$diff >0)] <- 1 ; Formosa_Aero$sg[which(Formosa_Aero$diff <=0)] <- 0

#Genera las fechas de inicio y final de los periodos calidos/frios con la funcion rle, que cuenta cantidad de casos
reclass <- Formosa_Aero$sg #vector aparte 
names(reclass) <- Formosa_Aero$Fecha #los nombres son las fechas
runs <- rle(reclass) #la funcion necesita que los datos sean organizados por valores

endingDates<-names(runs$values[runs$values==1 & runs$lengths >=2])
offset<-runs$lengths[which(names(runs$values) %in% endingDates)]-1 #objeto necesario para originar los inicios
startingDates <- names(reclass)[which(names(reclass) %in% endingDates) - offset]

#Guarda los periodos en un data.frame
longitud_periodos <- runs$lengths[runs$values==1 & runs$lengths >=2] #toma periodos mayores a 2 que superen el umbral de Tmean+-SD
fechas <- data.frame(startingDates,endingDates,longitud_periodos,row.names = NULL)
fechas$Tmmean <- NA ; fechas$Tmabs <- NA
print(paste("Se encontraron",dim(fechas)[1],"eventos",periodo,"para el periodo 1963-2003"),quote=F)

#Calcula la Tx / Tn media de cada periodo y la absoluta, luego las guarda en el mismo data.frame
if(periodo=="calidos"){
for(i in 1:dim(fechas)[1]){
  fechas$Tmmean[i] <- round(mean(Formosa_Aero$Tx[which(Formosa_Aero$Fecha == as.Date(fechas[i,1])):which(Formosa_Aero$Fecha == as.Date(fechas[i,2]))]),digits=1)
  fechas$Tmabs[i] <- max(Formosa_Aero$Tx[which(Formosa_Aero$Fecha == as.Date(fechas[i,1])):which(Formosa_Aero$Fecha == as.Date(fechas[i,2]))])
}
names(fechas)=c("Inicio","Final","Longitud (dias)","Tx media (°C)","Tx abs (°C)")} else{
  
for(i in 1:dim(fechas)[1]){
  fechas$Tmmean[i] <- round(mean(Formosa_Aero$Tn[which(Formosa_Aero$Fecha == as.Date(fechas[i,1])):which(Formosa_Aero$Fecha == as.Date(fechas[i,2]))]),digits=1)
  fechas$Tmabs[i] <- min(Formosa_Aero$Tn[which(Formosa_Aero$Fecha == as.Date(fechas[i,1])):which(Formosa_Aero$Fecha == as.Date(fechas[i,2]))])
}
names(fechas)=c("Inicio","Final","Longitud (dias)","Tn media (°C)","Tn abs (°C)")}

#Genera un archivo .txt con la informacion de la estacion y los periodos
write.table(fechas,file=paste(OUTPUTS,"Periodos",periodo,"en Formosa AERO 1963-2003.txt"),sep="\t",quote=F,row.names = F)

######CANTIDAD DE PERIODOS CALIDOS /FRIOS POR SEMESTRE CALIDO Y FRIO######
#Observacion: se considera que cada evento pertenece al semestre del dia en que se inicio
#Se calcula la cantidad de eventos en todo el periodo correspondientes a cada mes del anio, luego se distingue por semestre calido y frio
#y se calcula la media de cada semestre. Se realiza el calculo tambien para eventos de duracion mayor a 3 dias.
meses_calidos <- medias[order(medias$Tmean),][7:12,1] ; meses_frios <- medias[order(medias$Tmean),][1:6,1]
meses <- c(1:12) ; cantidad <- NA ; cantidad_3dias <- NA

for(i in meses){
  cantidad[i] <- length(which(month(fechas$Inicio) == meses[i]))
  cantidad_3dias[i] <- length(which(month(fechas$Inicio) == meses[i] & fechas$`Longitud (dias)`> 3))
}
print(paste("La media de eventos",periodo, "durante el semestre calido es de",round(mean(cantidad[meses_calidos]))),quote=F)
print(paste("La media de eventos",periodo, "durante el semestre frio es de",round(mean(cantidad[meses_frios]))),quote=F)

######GENERACION DE GRAFICOS DE BARRAS PARA DOS TIPOS DE EVENTOS CALIDOS /FRIOS SEGUN LONGITUD######
datos <- data.frame(meses,cantidad,cantidad_3dias) #datos acomodados para ser leidos por ggplot2
paleta <- rainbow(12,alpha = 0.7) #escala representativa del ciclo anual de la temperatura
my_theme <- theme(plot.title = element_text(hjust = 0.5,size = 12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10)) 

#Periodo +1 dias
p <- ggplot(data=datos,aes(x=datos$meses,y=cantidad)) + geom_bar(stat = "identity",fill=paleta) + theme_light() +
     xlab("Meses") + ylab("Cantidad") + scale_x_continuous(breaks = meses) + 
     ggtitle("Duración mayor a 1 dia") + my_theme +ylim(0,max(cantidad))
#Periodo +3 dias
g <- ggplot(data=datos,aes(x=datos$meses,y=cantidad_3dias)) + geom_bar(stat = "identity",fill=paleta) + theme_light() +
     xlab("Meses") + ylab("Cantidad") + scale_x_continuous(breaks = meses) + 
     ggtitle("Duración mayor a 3 dias")  + my_theme +ylim(0,max(cantidad))
ggsave(filename=paste(OUTPUTS,"barras_periodos",periodo,".jpeg"),
       plot = grid.arrange(p, g, ncol=2, top=paste("FORMOSA AERO: Periodos",periodo ,"por mes 1963-2003")),scale = 2,width=5,height =2.5)
###################################################################################################################
######ESTUDIO DEL CASO DE LONGITUD MAXIMA DE UN EVENTO CALIDO: evolucion de Tx, Tn y diff######
#Identifico el evento cuya longitud fue la maxima encontrada en todo el periodo
ini <- as.Date(fechas$Inicio[which.max(fechas$`Longitud (dias)`)])
fin <- as.Date(fechas$Final[which.max(fechas$`Longitud (dias)`)])
evento <- Formosa_Aero[(which(Formosa_Aero$Fecha==ini)-3):(which(Formosa_Aero$Fecha==fin)+3),]

p <- ggplot(data=evento,aes(x=Fecha)) + geom_col(aes(y=Tx,fill="Tx"),position = "stack") + theme_light() + my_theme
p <- p + geom_col(aes(y=Tn,fill="Tn")) + xlab("Tiempo") + ylab("Temperatura (°C)") + 
     ggtitle("FORMOSA AERO: evolución de la temperatura desde el 27-05-2001 al 18-06-2001") + scale_x_date(breaks=date_breaks("3 days"))
p <- p + geom_point(aes(y=Tmean,colour="Tmean"),size=2.5) + geom_line(aes(y=Tmsd,colour="Umbral"),size=1,linetype="twodash") 
p <- p + scale_fill_manual(values=c("lightskyblue","indianred1"),name=NULL) + scale_colour_manual(values=c("purple","black"),name=NULL)
ggsave(filename=paste(OUTPUTS,"Evento de mayor longitud.jpeg"),plot=p,scale=2,height=2.5,width=5)
###################################################################################################################
######CICLO ANUAL DE LA TEMPERATURA######
#Organizo los datos en un mismo data.frame y ploteo la onda anual y, en sombreado, el desvio estandar
medias$sd <- desvios$Tmean
my_theme <- theme(plot.title = element_text(hjust = 0.5,size = 12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10),legend.position = "bottom") 

p <- ggplot(data=medias,aes(x=meses)) +  geom_ribbon(aes(ymin=Tmean - sd,ymax=Tmean + sd,fill="Desvio estandar"),stat="identity",alpha=0.6) + theme_light()
p <- p + geom_line(aes(y=Tmean,colour="Valor medio"),size=1) + scale_x_continuous(breaks = meses) + xlab("Meses") + ylab("Temperatura (°C)") + my_theme
p <- p + geom_point(aes(y=Tmean),colour="red",size=2.5) + ggtitle("FORMOSA AERO: Onda anual de temperatura (1963-2003)") 
p <- p + scale_colour_manual(values=c("red"),name=NULL) + scale_fill_manual(values="grey77",name=NULL)
ggsave(filename=paste(OUTPUTS,"Onda anual.jpeg"),plot=p,scale=2,height=2.5,width=5)
