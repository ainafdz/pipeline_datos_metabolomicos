library(shiny)
library(data.table)
library(ggplot2)
library(DT)
library(RColorBrewer)
library(DMwR)
library(DiscriMiner)
library("psych")
library(PerformanceAnalytics)
library(leaflet)
library(dplyr)
library(FactoMineR)


load("workspace_3.RData")
dd <- read.csv("caracteristicas_localizaciones.csv", header=TRUE, sep=";", dec=",")

ui <- fluidPage(
  
  tabsetPanel(
    tabPanel("Base de datos original", fluid = TRUE,
               fluidRow(
                 column(4,
                        selectInput("vari1",
                                    "Variedad:",
                                    c("All",
                                      unique(as.character(datos$Variedad))))
                        ),
                   column(4,
                        selectInput("metabolito1", "Metabolito:", 
                                    choices=c("All", colnames(datos)[-1]))
                        )
                 ),
               # Create a new row for the table.
               DT::dataTableOutput("mytable")),
               
    tabPanel("Estudio de valores faltantes", fluid=TRUE,
             sidebarLayout(
               sidebarPanel(
                      selectInput("vari2",
                                  "Variedad:",
                                  c("All",
                                    unique(as.character(datos$Variedad))))
               ),
             mainPanel(
               tabsetPanel(type = "tabs",
                 tabPanel("Valores faltantes reales",
                          textOutput("missingsText"),
                          plotOutput("missingsPlot")),
                 tabPanel( "Valores faltantes por nivel de deteccion",
                           textOutput("zeroText"),
                           plotOutput("zeroPlot"))
               )
             )
             )
  ),
  tabPanel("Analisis descriptivo inicial", fluid=TRUE,
           sidebarLayout(
             sidebarPanel(
               selectInput("metabolito2", "Metabolito:", 
                           choices=c(colnames(datos)[-1]))
             ),
             mainPanel (
              tabsetPanel(type="tabs",
                  tabPanel("Descriptiva por metabolitos",
                               plotOutput("box_met"),
                               DT::dataTableOutput("descriptive2")),
                  tabPanel("Descriptiva separando por variedades",
                               plotOutput("box_var"),
                               DT::dataTableOutput("descriptive_var"))
             )
           )
           )
           
  ),
  tabPanel("Preprocessing", fluid=TRUE,
           sidebarLayout(
             sidebarPanel(
               radioButtons("missings", "Imputacion de missings:",
                            c("No imputar",
                              "KNNImputation")),
               radioButtons("transformacion", "Transformacion:",
                            c("No transformar",
                              "Logaritmica")),
               radioButtons("normscale", "Normalizacion y escalado:",
                            c("No normalizar/escalar",
                              "Escala de Pareto"))
               
             ),
             
             mainPanel(
               plotOutput("box_mult"),
               plotOutput("box_mult_desp")
             )
           )
    ),
  tabPanel("Analisis univariante", fluid=TRUE,
           sidebarLayout(
             sidebarPanel(
               selectInput("metabolito3", "Metabolito:", 
                           choices=c(colnames(datos)[-1]))
             ),
             mainPanel(
             tabsetPanel(
               tabPanel("Descriptiva por metabolitos",
                          plotOutput("box_fin"),
                          DT::dataTableOutput("descriptive3"),
                          DT::dataTableOutput("anova")),
               tabPanel("Descriptiva separando por variedades", 
                        plotOutput("box_var_fin"),
                        DT::dataTableOutput("descriptive3.1")))
               )
           )
    ),
  tabPanel("Analisis multivariante", fluid=TRUE,
           tabsetPanel(
             tabPanel("PCA",
                      sidebarPanel(
                        helpText("NOTA: Se considera la imputacion de missings KNN aunque no este seleccionada, ya que el analisis PCA no se puede realizar si hay valores faltantes.")
                      ),
                      mainPanel(
                      plotOutput("scores_PCA"),
                      plotOutput("load_PCA"))
                      
                      ),
             tabPanel("Heatmap", 
                      sidebarPanel(
                        helpText("NOTA: Si se selecciona el escalado de los datos se debe seleccionar tambien la imputacion de valores faltantes.")
                      ),
                      plotOutput("heatmap")
                      ),
             tabPanel("PLS-DA",
                      sidebarPanel(
                        helpText("NOTA: Se considera la imputacion de missings KNN aunque no este seleccionada, ya que el analisis PLS-DA no se puede realizar si hay valores faltantes.")
                      ),
                      mainPanel(
                      plotOutput("plsda"))
                      ),
             tabPanel("Correlaciones",
                      tabsetPanel(
                        tabPanel("Entre metabolitos",
                                 sidebarPanel(
                                   selectInput("metabolito4", "Metabolito 1:", 
                                               choices=c(colnames(datos)[-1])),
                                   selectInput("metabolito5", "Metabolito 2:", 
                                               choices=c(colnames(datos)[-1]))
                                 ),
                                 mainPanel(
                                   plotOutput("corr_metab")
                                 )
                        ),
                        tabPanel("Entre muestras de vino",
                                 sidebarPanel(
                                   selectInput("vino1", "Muestra 1:", 
                                               choices=c(rownames(datos))),
                                   selectInput("vino2", "Muestra 2:", 
                                               choices=c(rownames(datos)))
                                 ),
                                 mainPanel(
                                   plotOutput("corr_vinos")
                                 )
                        )
                                 
                        )
                      )
                      
                      )
                      
           
  ),
  tabPanel("Integración con las caracteristicas del terreno", fluid=TRUE,
           tabsetPanel(
             tabPanel("Caracteristicas de las localizaciones",
                    mainPanel(
                      leafletOutput("mymap"))),
             tabPanel("Analisis de factores multiples",
                                 sidebarLayout(
                                   sidebarPanel(
                                     selectInput("dim1", "Dimension:", 
                                                 choices=c(1,2,3)),
                                     selectInput("dim2", "Dimension:", 
                                                 choices=c(1,2,3), selected = 2)
                                     
                                   ),
                                   mainPanel(
                                     tabsetPanel(type = "tabs",
                                                 tabPanel("Ejes parciales",
                                                          plotOutput("partialaxes12")),
                                                 tabPanel("Variables",
                                                          plotOutput("variables")),
                                                 tabPanel("Individuos",
                                                          plotOutput("individuals")),
                                                 tabPanel("Puntos parciales",
                                                          plotOutput("puntos"))
                                                 
                                     )
                                   )
                                 
                        )
                      ))
           )
           
    
           
  
  
  
           
)
)


  




  




server <- function(input, output) {
  
  # BASE DE DATOS
  output$mytable = DT::renderDataTable({
    data <- datos
    data2 <- datos
    if (input$vari1 != "All") {
      data <- data[data$Variedad == input$vari1,]
      data2 <- data
    }
    if (input$metabolito1 != "All") {
      data2 <- as.data.frame(data[,input$metabolito1])
      rownames(data2) <- rownames(data)
      names(data2) <- input$metabolito1
    }
    data2
  })
  

  # MISSINGS: GRÁFICO
  output$missingsPlot <- renderPlot({
    if (input$vari2 != "All") {
      datos <- datos[datos$Variedad == input$vari2,]
    }
    buscarNA <- function(x){
      s <- 0
      for(i in 1:length(x)){
        if(is.na(x[i]==TRUE)){
          s <- s  + 1
        }
      }
      return(s*100/length(x))
    }
    
    m <- round(apply(datos, 2, buscarNA),2)
    if(sum(m)<=0){

    }else{
      hist(m, col="steelblue", 
           main="Porcentage de missings por variable",
           ylab="Num. Variables",
           xlab="% Missings")
    }
    
  })
  
  # MISSINGS: Mensaje de num. de missings
  output$missingsText <- renderText({
    if (input$vari2 != "All") {
      datos <- datos[datos$Variedad == input$vari2,]
    }
    buscarNA <- function(x){
      s <- 0
      for(i in 1:length(x)){
        if(is.na(x[i]==TRUE)){
          s <- s  + 1
        }
      }
      return(s*100/length(x))
    }
    total_missings <- round(buscarNA2(datos),2)
    m <- round(apply(datos, 2, buscarNA),2)
    if(sum(m)<=0){
      print("No hay valores missing para esta seleccion.")
    }else{
      paste("En esta seleccion hay un total de", total_missings[1], "missings, que representan el", 
            total_missings[2], "% del total de datos de la seleccion.")
    }
    
  })
  
  
  # VALORES = 0: GRÁFICO
  output$zeroPlot <- renderPlot({
    if (input$vari2 != "All") {
      datos <- datos[datos$Variedad == input$vari2,]
    }
    buscarNA2 <- function(x){
      l <- 0
      for(i in 1:ncol(datos)){
        l <- sum(datos[,i]==0, na.rm=TRUE)
      }
      return(l*100/length(x))
    }
    n <- round(apply(datos, 2, buscarNA2),2)
    if(sum(n)<=0){
      
    }else{
    hist(n, col="steelblue", 
         main="Porcentage de valores cero por variable",
         ylab="Num. Variables",
         xlab="% Valores 0")
    }
  })
  
  # VALORES = 0: TEXTO
  output$zeroText <- renderText({
    if (input$vari2 != "All") {
      datos <- datos[datos$Variedad == input$vari2,]
    }
    buscarNA2 <- function(x){
      l <- 0
      for(i in 1:ncol(datos)){
        l <- sum(datos[,i]==0, na.rm=TRUE)
      }
      return(c(l, l*100/length(x)))
    }
    n <- round(apply(datos, 2, buscarNA2),2)
    total_ceros <- round(buscarNA2(datos), 2)
    if(sum(n)<=0){
      print("No hay valores iguales a cero para esta seleccion.")
    }else{
      paste("En esta seleccion hay un total de", total_ceros[1], "valores cero, que representan el", 
            total_ceros[2], "% del total de datos de la seleccion.")
    }
  
  })
  

  # DESCRIPTIVO POR METABOLITOS
    output$descriptive2 <- DT::renderDataTable({
    descrip2 <- round(estadisticos.metab,2)
    descrip2 <- as.data.frame(descrip2[input$metabolito2,])
    descrip2
  })
    
  # DESCRIPTIVO POR VARIEDADES
    output$descriptive_var <- DT::renderDataTable({
      descrip_vari <- subset(descrip_var, descrip_var$Metabolito==as.character(input$metabolito2))
      descrip_vari <- descrip_vari[,-(1:2)]
      rownames(descrip_vari) <- c("Chardonnay", "Pinot Gris", "Riesling", "Sauvignon Blanc", "Viognier")
      descrip_vari <- round(descrip_vari, 2)
      descrip_vari
      }) 
  
  # BOXPLOTS DE CADA METABOLITO
  output$box_met <- renderPlot({
    boxplot(datos[,input$metabolito2], col=brewer.pal(3,"Blues"), main=as.character(input$metabolito2))
    
  })

  # BOXPLOTS POR VARIEDADES
  output$box_var <- renderPlot({
    boxplot(datos[,input$metabolito2]~datos$Variedad, col=brewer.pal(7,"Blues"), main=as.character(input$metabolito2))
    
  })
  
  # BOXPLOTS MULTIPLES INICIALES
  output$box_mult <- renderPlot({
    boxplot(datos[2:ncol(datos)],datos,horizontal=F,names=FALSE, main="Raw data" )
  })
  

  # BOXPLOTS MULTIPLES DESPUES DEL PREPROCESSING
  output$box_mult_desp <- renderPlot({
    datos_impu <- datos[-1]
    if (input$missings == "KNNImputation") {
      datos_impu <- knnImputation(datos_impu)
    }
    if (input$transformacion == "Logaritmica") {
      datos_impu <- log(datos_impu)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto <- paretoscale(t(datos_impu))
      datos_impu <- t(pareto)
    }
    
    boxplot(datos_impu,horizontal=F,names=FALSE, main="Datos preprocesados" )
    
  })
  
  # MENSAJE: HAY ALGUN DATO FALTANTE
  output$nmissings <- renderText({
    datos_impu2 <- datos
    if (input$missings == "KNNImputation") {
      datos_impu2 <- knnImputation(datos)
    }
    paste("Hay algun valor faltante:", anyNA(datos_impu2))
    
  })
  
  # BOXPLOTS DE CADA METABOLITO DESPUÉS DEL PREPROCESSING
  output$box_fin <- renderPlot({
    datos_impu4.1 <- datos[,-1]
    if (input$missings == "KNNImputation") {
      datos_impu4.1 <- knnImputation(datos_impu4.1)
    }
    if (input$transformacion == "Logaritmica") {
      datos_impu4.1 <- log(datos_impu4.1)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto4.1 <- paretoscale(t(datos_impu4.1))
      datos_impu4.1 <- t(pareto4.1)
    }
    boxplot(datos_impu4.1[,input$metabolito3], col=brewer.pal(3,"Blues"), main=as.character(input$metabolito3))
    
  })
  
  # BOXPLOTS POR VARIEDADES DESPUÉS DEL PREPROCESING
  output$box_var_fin <- renderPlot({
    datos_impu4 <- datos[,-1]
    if (input$missings == "KNNImputation") {
      datos_impu4 <- knnImputation(datos_impu4)
    }
    if (input$transformacion == "Logaritmica") {
      datos_impu4 <- log(datos_impu4)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto4 <- paretoscale(t(datos_impu4))
      datos_impu4 <- as.data.frame(t(pareto4))
    }
    datos_impu4 <- cbind("Variedad" = datos$Variedad, datos_impu4)
    boxplot(datos_impu4[,input$metabolito3]~datos_impu4$Variedad, col=brewer.pal(7,"Blues"), main=as.character(input$metabolito3))
    
  })
  
  # DESCRIPTIVO POR METABOLITOS DESPUÉS DEL PREPROCESING
  output$descriptive3 <- DT::renderDataTable({
  datos_impu3 <- datos[,-1]
  if (input$missings == "KNNImputation") {
    datos_impu3 <- knnImputation(datos_impu3)
  }
  if (input$transformacion == "Logaritmica") {
    datos_impu3 <- log(datos_impu3)
  }
  if (input$normscale == "Escala de Pareto") {
    paretoscale <- function(z) {
      rowmean <- apply(z, 1, mean) # row means
      rowsd <- apply(z, 1, sd)  # row standard deviation
      rowsqrtsd <- sqrt(rowsd) # sqrt of sd
      rv <- sweep(z, 1, rowmean,"-")  # mean center
      rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
      return(rv)
    }
    pareto3 <- paretoscale(t(datos_impu3))
    datos_impu3 <- as.data.frame(t(pareto3))
  }
  sumstats <- function(z) {
     Media <- apply(z, 1, mean, na.rm=TRUE)
     Mediana <- apply(z, 1, median, na.rm=TRUE)
     SD <- apply(z, 1, sd, na.rm=TRUE)
     CV <- apply(z, 1, function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
     result <- data.frame(Media, Mediana, SD, CV)
     return(result)
   }
   
   estadisticos.metab2 <- sumstats(t(datos_impu3))
     descrip3 <- round(estadisticos.metab2,2)
     if (input$metabolito3 != "All") {
       descrip3 <- as.data.frame(descrip3[input$metabolito3,])
     }
     if(descrip3$Media <= 0 ){
       descrip3 <- descrip3[,-4]
     }
    descrip3
  })
  
  # DESCRIPTIVA SEPARANDO POR VARIEDADES DESPUÉS DEL PREPROCESSING
  output$descriptive3.1 <- DT::renderDataTable({
    datos_impu3.1 <- datos[,-1]
    if (input$missings == "KNNImputation") {
      datos_impu3.1 <- knnImputation(datos_impu3.1)
    }
    if (input$transformacion == "Logaritmica") {
      datos_impu3.1 <- log(datos_impu3.1)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto3.1 <- paretoscale(t(datos_impu3.1))
      datos_impu3.1 <- as.data.frame(t(pareto3.1))
    }
    mitj <- data.frame()
    med <- data.frame()
    des.st <- data.frame()
    metab <- rep(colnames(datos)[-1], each=5)
    for(i in 1:ncol(datos_impu3.1)){
      mitj <- rbind(mitj, aggregate(datos_impu3.1[,i], by=list(datos$Variedad),mean, na.rm=TRUE))
      med <- rbind(med, aggregate(datos_impu3.1[,i], by=list(datos$Variedad),median, na.rm=TRUE))
      des.st <- rbind(des.st, aggregate(datos_impu3.1[,i], by=list(datos$Variedad),sd, na.rm=TRUE))
    }
    descrip_var_fin <- cbind(metab, mitj, med, des.st)
    descrip_var_fin <- descrip_var_fin[,-c(4,6)]
    colnames(descrip_var_fin) <- c("Metabolito", "Variedad", "Media", "Mediana", "SD")
    descrip_var_fin$CV <- descrip_var_fin$SD/descrip_var_fin$Media
    descrip_vari_OK <- subset(descrip_var_fin, descrip_var_fin$Metabolito==as.character(input$metabolito3))
    descrip_vari_OK <- descrip_vari_OK[,-(1:2)]
    rownames(descrip_vari_OK) <- c("Chardonnay", "Pinot Gris", "Riesling", "Sauvignon Blanc", "Viognier")
    descrip_vari_OK <- round(descrip_vari_OK, 2)
    descrip_vari_OK
  })   
  # ANOVA SELECCIONANDO METABOLITOS
  output$anova <- DT::renderDataTable({
    datos_impu5 <- datos[,-1]
    if (input$missings == "KNNImputation") {
      datos_impu5<- knnImputation(datos_impu5)
    }
    if (input$transformacion == "Logaritmica") {
      datos_impu5 <- log(datos_impu5)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto5 <- paretoscale(t(datos_impu5))
      datos_impu5 <- as.data.frame(t(pareto5))
    }
    datos_impu5 <- cbind("Variedad" = datos$Variedad, datos_impu5)
    v <- vector()
    for(i in 2:110){
       aux <- summary(aov(datos_impu5[,i] ~ datos_impu5$Variedad))
       v <- c(v, aux[[1]][1,5])
      }
    v.adj <- p.adjust(v)
    df.anova <- matrix(c(v, v.adj), ncol=2)
    colnames(df.anova) <- c("P-value", "Adj. P-value")
    rownames(df.anova) <- colnames(datos_impu5)[-1]
    df.anova <- as.data.frame(df.anova[input$metabolito3,])
    names(df.anova) <- as.character(input$metabolito3)
    df.anova
    
    
  })
  
  # ANOVA TOTAL
  output$anovatot <- DT::renderDataTable({
    datos_impu6 <- datos[,-1]
    if (input$missings == "KNNImputation") {
      datos_impu6<- knnImputation(datos_impu6)
    }
    if (input$transformacion == "Logaritmica") {
      datos_impu6 <- log(datos_impu6)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto6 <- paretoscale(t(datos_impu6))
      datos_impu6 <- as.data.frame(t(pareto6))
    }
    datos_impu6 <- cbind("Variedad" = datos$Variedad, datos_impu6)
    v <- vector()
    for(i in 2:110){
      aux <- summary(aov(datos_impu6[,i] ~ datos_impu6$Variedad))
      v <- c(v, aux[[1]][1,5])
    }
    v.adj <- p.adjust(v)
    df.anova2 <- matrix(c(v, v.adj), ncol=2)
    colnames(df.anova2) <- c("P-value", "Adj. P-value")
    rownames(df.anova2) <- colnames(datos_impu6)[-1]
    df.anova2
    
    
  })
  
  # MENSAJE: los datos estan imputados sin missings (sino no va)
  output$text_impumissings <- renderText({
    "Se considera la imputación de missings KNN aunque no esté seleccionada, ya que el análisis PCA no se puede realizar si hay valores faltantes."
    
  })

  # SCORE PLOT PCA
  output$scores_PCA <- renderPlot({
    datos_impu7 <- datos[,-1]
    datos_impu7 <- knnImputation(datos_impu7)
    if (input$transformacion == "Logaritmica") {
      datos_impu7 <- log(datos_impu7)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto7 <- paretoscale(t(datos_impu7))
      datos_impu7 <- as.data.frame(t(pareto7))
    }
    pc <- prcomp(datos_impu7, retx=TRUE)
    PVE2 <- 100*pc$sdev^2/sum(pc$sdev^2)
    pcaresults <- summary(pc)
    score_data <- as.data.frame(pcaresults$x)
    score_data <- score_data[,1:3]
    score_data <- cbind(score_data, "Variedad"=datos$Variedad)
    score_data <- as.data.frame(score_data)
    
    ggplot(score_data, aes(PC1, PC2)) +
     geom_point(aes(color=Variedad)) + 
     geom_text(aes(label=rownames(score_data), color=Variedad)) +
     stat_ellipse(aes(color=Variedad)) +
     ggtitle("PCA Scores Plot") +
     theme(plot.title=element_text(size=15, vjust=2, face="bold")) +
     geom_hline(yintercept=0, size=0.25) +
     geom_vline(xintercept=0, size=0.25) +
     xlab((paste0("PC1"," ", "(",round(PVE2[1], 1), "%", ")"))) + 
     ylab((paste0("PC2"," ", "(",round(PVE2[2], 1), "%", ")")))


  })
  
  # LOADING PLOT PCA
  output$load_PCA <- renderPlot({
    datos_impu8 <- datos[,-1]
    datos_impu8 <- knnImputation(datos_impu8)
    if (input$transformacion == "Logaritmica") {
      datos_impu8 <- log(datos_impu8)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto8 <- paretoscale(t(datos_impu8))
      datos_impu8 <- as.data.frame(t(pareto8))
    }
    pc1 <- prcomp(datos_impu8, retx=TRUE)
    PVE3 <- 100*pc1$sdev^2/sum(pc1$sdev^2)
    pcaresults1 <- summary(pc1)
    load_data <- as.data.frame(pcaresults1$rotation)
    load_data <- load_data[,1:3]
    load_data <- as.data.frame(load_data)
    
    ggplot(load_data, aes(PC1, PC2)) +
      geom_text(aes(label=rownames(load_data))) +
      ggtitle("PCA Loadings Plot") +
      theme(plot.title=element_text(size=15, vjust=2, face="bold")) +
      geom_hline(yintercept=0, size=0.25) +
      geom_vline(xintercept=0, size=0.25) +
      xlab((paste0("PC1"," ", "(",round(PVE3[1], 1), "%", ")"))) + 
      ylab((paste0("PC2"," ", "(",round(PVE3[2], 1), "%", ")")))
    
    
  })
  
  # HEATMAP
  output$heatmap <- renderPlot({
    datos_impu9 <- datos[,-1]
    if (input$missings == "KNNImputation") {
      datos_impu9 <- knnImputation(datos_impu9)
    }
    if (input$transformacion == "Logaritmica") {
      datos_impu9 <- log(datos_impu9)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto9 <- paretoscale(t(datos_impu9))
      datos_impu9 <- as.data.frame(t(pareto9))
    }
    color <- datos$Variedad
    levels(color) <- c("lightgoldenrod", "lightpink", "darkseagreen", "coral1", "cyan")
    color <- as.character(color)
    
    library(gplots)
    heatmap.2(x = t(datos_impu9), scale = "none", col = bluered(256),
              distfun = function(x){dist(x, method = "euclidean")},
              hclustfun = function(x){hclust(x, method = "average")},
              density.info = "none",
              trace = "none", cexRow = 0.7,
              ColSideColors = color)
  })
  
  # PLS-DA
  output$plsda <- renderPlot({
    datos_impu10 <- datos[,-1]
    datos_impu10 <- knnImputation(datos_impu10)
    
    if (input$transformacion == "Logaritmica") {
      datos_impu10 <- log(datos_impu10)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto10 <- paretoscale(t(datos_impu10))
      datos_impu10 <- as.data.frame(t(pareto10))
    }
    my_pls2 = plsDA(datos_impu10, datos$Variedad, autosel=TRUE)
    plot(my_pls2)
  })
  
  # CORRELACIONES ENTRE METABOLITOS
  output$corr_metab <- renderPlot({
    datos_impu11 <- datos[,-1]
    datos_impu11 <- knnImputation(datos_impu11)
    
    if (input$transformacion == "Logaritmica") {
      datos_impu11 <- log(datos_impu11)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto11 <- paretoscale(t(datos_impu11))
      datos_impu11 <- as.data.frame(t(pareto11))
    }
    df.aux <- data.frame(datos_impu11[,input$metabolito4], datos_impu11[,input$metabolito5])
    names(df.aux) <- c(as.character(input$metabolito4), as.character(input$metabolito5))
    chart.Correlation(df.aux, histogram=F, method = "spearman")
    

  })
  # CORRELACIONES ENTRE MUESTRAS DE VINO
  output$corr_vinos <- renderPlot({
    datos_impu12 <- datos[,-1]
    datos_impu12 <- knnImputation(datos_impu12)
    
    if (input$transformacion == "Logaritmica") {
      datos_impu12 <- log(datos_impu12)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto12 <- paretoscale(t(datos_impu12))
      datos_impu12 <- as.data.frame(t(pareto12))
    }
    datos_impu12_2 <- as.data.frame(t(datos_impu12))
    df.aux2 <- data.frame(datos_impu12_2[,input$vino1], datos_impu12_2[,input$vino2])
    names(df.aux2) <- c(as.character(input$vino1), as.character(input$vino2))
    chart.Correlation(df.aux2, histogram=F, method = "spearman")
    
    
  })
  
  # MAPA
  data <- reactive({
    x <- dd
  })
  output$mymap <- renderLeaflet({
    dd <- data()
    
    m <- leaflet(data = dd) %>%
      addTiles() %>%
      addMarkers(lng = ~longitude,
                 lat = ~latitude,
                 popup = paste("Localizacion:", dd$location, "<br>",
                               "Temperatura: ", dd$Temperatura, "<br>",
                               "Precipitaciones:", dd$Precipitaciones, "<br>",
                               "Altitud:", dd$Altitud, "<br>"))
    m
  })
  
  # PARTIAL AXES
  output$partialaxes12 <- renderPlot({
    carac <- as.data.frame(aux[,c(4:6)])

      datos_impu13 <- datos[,-1]
    if (input$missings == "KNNImputation") {
      datos_impu13 <- knnImputation(datos_impu13)
    }
    if (input$transformacion == "Logaritmica") {
      datos_impu13 <- log(datos_impu13)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto13 <- paretoscale(t(datos_impu13))
      datos_impu13 <- as.data.frame(t(pareto13))
    }
    
    datos_impu13 <- cbind(datos_impu13, datos$Variedad, carac)
    mfa_1 = MFA(datos_impu13, group=c(109,1,3), type=c("s", "n","s"), name.group=c("metabolitos","variedad","caracteristicas"), ncp=3, num.group.sup = 2, graph = FALSE)
    plot(mfa_1, choix="axes", axes=c(as.numeric(input$dim1),as.numeric(input$dim2)))
  })
  
  # VARIABLES
  output$variables <- renderPlot({
    carac <- as.data.frame(datos_int[,c(111:113)])
    
    datos_impu14 <- datos[,-1]
    if (input$missings == "KNNImputation") {
      datos_impu14 <- knnImputation(datos_impu14)
    }
    if (input$transformacion == "Logaritmica") {
      datos_impu14 <- log(datos_impu14)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto14 <- paretoscale(t(datos_impu14))
      datos_impu14 <- as.data.frame(t(pareto14))
    }
    
    datos_impu14 <- cbind(datos_impu14, datos$Variedad, carac)
    mfa_2 = MFA(datos_impu14, group=c(109,1,3), type=c("s", "n","s"), name.group=c("metabolitos","variedad","caracteristicas"), ncp=3, num.group.sup = 2, graph = FALSE)
    plot(mfa_2, habillage=110, cex=0.8, select="contrib 20", axes=c(as.numeric(input$dim1),as.numeric(input$dim2)))
    })
  
  # INDIVIDUALS
  output$individuals <- renderPlot({
    carac <- as.data.frame(datos_int[,c(111:113)])
    
    datos_impu15 <- datos[,-1]
    if (input$missings == "KNNImputation") {
      datos_impu15 <- knnImputation(datos_impu15)
    }
    if (input$transformacion == "Logaritmica") {
      datos_impu15 <- log(datos_impu15)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto15 <- paretoscale(t(datos_impu15))
      datos_impu15 <- as.data.frame(t(pareto15))
    }
    
    datos_impu15 <- cbind(datos_impu15, datos$Variedad, carac)
    mfa_3 = MFA(datos_impu15, group=c(109,1,3), type=c("s", "n","s"), name.group=c("metabolitos","variedad","caracteristicas"), ncp=3, num.group.sup = 2, graph = FALSE)
    plot(mfa_3, choix="var", habillage="group", cex=0.8, select="contrib 10", unselect=1, axes=c(as.numeric(input$dim1),as.numeric(input$dim2)))
    })
  
  # PUNTOS PARCIALES
  output$puntos <- renderPlot({
    carac <- as.data.frame(datos_int[,c(111:113)])
    
    datos_impu16 <- datos[,-1]
    if (input$missings == "KNNImputation") {
      datos_impu16 <- knnImputation(datos_impu16)
    }
    if (input$transformacion == "Logaritmica") {
      datos_impu16 <- log(datos_impu16)
    }
    if (input$normscale == "Escala de Pareto") {
      paretoscale <- function(z) {
        rowmean <- apply(z, 1, mean) # row means
        rowsd <- apply(z, 1, sd)  # row standard deviation
        rowsqrtsd <- sqrt(rowsd) # sqrt of sd
        rv <- sweep(z, 1, rowmean,"-")  # mean center
        rv <- sweep(rv, 1, rowsqrtsd, "/")  # divide by sqrtsd
        return(rv)
      }
      pareto16 <- paretoscale(t(datos_impu16))
      datos_impu16 <- as.data.frame(t(pareto16))
    }
    
    datos_impu16 <- cbind(datos_impu16, datos$Variedad, carac)
    mfa_4 = MFA(datos_impu16, group=c(109,1,3), type=c("s", "n","s"), name.group=c("metabolitos","variedad","caracteristicas"), ncp=3, num.group.sup = 2, graph = FALSE)
    plot(mfa_4, choix="ind", invisible="ind", habillage="group", cex=0.8, partial="all", axes=c(as.numeric(input$dim1),as.numeric(input$dim2)))
    
    
  })
  
}

  



shinyApp(ui = ui, server = server)