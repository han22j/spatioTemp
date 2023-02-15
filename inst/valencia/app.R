#package load###########################################################
#rm(list=ls())
# Package names
library(sp)
library(stlnpp)
library(spatstat)
library(rgdal)
library(plot3D)
library(ggplot2)
library(corrplot)
library(patchwork)

##1. Data load#######################################
####obs_table
obs_table <- read.csv("./data/obs_data_112Valencia.csv", header = T)
load("./data/highway_road_linnet.RData")

xy<-obs_table[,c("crime_lon","crime_lat")]
coordinates(xy)<- c("crime_lon", "crime_lat")
proj4string(xy) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

CRS.new <- CRS("+proj=lcc +lat_1=40 +lat_0=40 +lon_0=0 +k_0=0.9988085293 +x_0=600000 +y_0=600000 +a=6378298.3 +b=6356657.142669561 +pm=madrid +units=m +no_defs")
xyc <- spTransform(xy, CRS.new)

obs_table$proj_lon <- xyc@coords[,1]
obs_table$proj_lat <- xyc@coords[,2]
num_cols <- unlist(lapply(obs_table, is.numeric))         # Identify numeric columns
num_names = names(obs_table[num_cols])

sample_obs = doBy::sampleBy(~crime_type, frac=0.01, data=obs_table) #crime_type,sample


##Intenstity estimation#####################################################################

valencia <-as.stlpp(sample_obs$proj_lon,sample_obs$proj_lat,sample_obs$crime_hour,L=road_linnet)
marks(valencia) = as.factor(sample_obs$crime_type)

lambda <- density(valencia,dimyx=512,diggle=TRUE)
denT = attr(lambda,"tempden")
denS = attr(lambda,"netint")
bwS = attr(lambda,"bw")[1]
bwT = attr(lambda,"bw")[2]

#crime type
a=split(valencia)
lambda1 <- density(a$Agresion,dimyx=512,diggle=TRUE)
lambda2 <- density(a$AlarmasMujer,dimyx=512,diggle=TRUE)
lambda3 <- density(a$Otros,dimyx=512,diggle=TRUE)
lambda4 <- density(a$Sustraccion,dimyx=512,diggle=TRUE)

#year
year_names = levels(factor(sample_obs$year))
marks(valencia) = as.factor(sample_obs$year)
a=split(valencia)

year_lambda_list <-list()
for (i in 1:length(year_names)){
  year_lambda_list[[i]] <- density(a[[i]],dimyx=512,diggle=TRUE)
}
names(year_lambda_list) <- year_names


#week
week_names = c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday")
marks(valencia) = as.factor(sample_obs$week_day_name)
a=split(valencia)

lambda_list <-list()
for (i in week_names){
  lambda_list[[i]] <- density(a[[i]],dimyx=512,diggle=TRUE)
}
names(lambda_list) <- week_names


####connected network########################################

load("./data/connetL.RData")

valenciaL <-as.stlpp(sample_obs$proj_lon,sample_obs$proj_lat,
                     sample_obs$crime_hour,L=connetL)     ##connected
marks(valenciaL) = as.factor(sample_obs$crime_type)

lambda <- density(valencia,dimyx=512,diggle=TRUE)
lambdaL = density(valenciaL,dimyx=512,diggle=TRUE)


##Second order characteristics########################################################

load("./data/kg_1per.RData")

k1 = list_kg_1per[[1]]
Envk = list_kg_1per[[2]]
g = list_kg_1per[[3]]
Envg = list_kg_1per[[4]]

set.seed(100)
time_range = range(valenciaL$data$t)
lam <- npoints(valenciaL)/(volume(valenciaL$domain)*(time_range[2]-time_range[1])) ##correct
sims <- lapply(X=1:99, function(i){
  rpoistlpp(lambda = lam,a=time_range[1],b= time_range[2],L=valenciaL$domain)  ##correct
})

##K function
Envkinvec <- lapply(1:length(sims), function(i){
  as.vector(Envk[[i]]$Kinhom)
})
Envkinmat <- do.call(rbind,Envkinvec)

minEnv <- apply(Envkinmat,2,quantile,prob=0.025)
maxEnv <- apply(Envkinmat,2,quantile,prob=0.975)

minEnv <- matrix(minEnv,byrow = F,nrow = 10)
maxEnv  <- matrix(maxEnv ,byrow = F,nrow = 10)

# par(mfrow=c(1,1))
# persp3D(k1$r,k1$t,k1$Kinhom/1000,theta=30,phi=10,
#         facets=FALSE,col=1,ticktype= "detailed",
#         xlab="r = distance",ylab="t = time",zlab="",
#         main=expression(italic({hat(K)[LI]^{ST}}(r,t))),
#         nticks=6,cex.axis=1.5,cex.lab=2, cex.main=2,
#         zlim= range(c(min(minEnv,k1$Kinhom),max(maxEnv,k1$Kinhom)))/1000,
#         resfac=2)
#
# persp3D(k1$r,k1$t,minEnv/1000,col = "gray80",add=T,facets=T,
#         ticktype= "detailed",
#         nticks=6,cex.axis=1.3,cex.lab=1.7,resfac=2,border = 1)
#
# persp3D(k1$r,k1$t,maxEnv/1000,col = "gray80",add=T,facets=T,
#         ticktype= "detailed",
#         nticks=6,cex.axis=1.3,cex.lab=1.7,resfac=2,border = 1)


##G function
set.seed(102)
sims_g <- lapply(X=1:99, function(i){
  rpoistlpp(lambda = lam,a=time_range[1],b= time_range[2],L=valenciaL$domain) ##correct
})

Envginvec <- lapply(1:length(sims_g), function(i){
  as.vector(Envg[[i]]$ginhom)
})
Envginmat <- do.call(rbind,Envginvec)

minEnvg <- apply(Envginmat,2,quantile,prob=0.025)
maxEnvg <- apply(Envginmat,2,quantile,prob=0.975)

minEnvg <- matrix(minEnvg,byrow = F,nrow = 10)
maxEnvg  <- matrix(maxEnvg ,byrow = F,nrow = 10)



### Global P-value - goodness-of-fit
##k
meansurf <- apply(Envkinmat,2,mean)
sdsurf <- apply(Envkinmat,2,sd)

val <- t((t(Envkinmat)-meansurf)/sdsurf)
valint <- rowSums(val)

valK <- t((t(as.vector(k1$Kinhom))-meansurf)/sdsurf)
valkint <- sum(valK)

pvalueK = (sum(valint>valkint)+1)/(length(valint)+1)
# [1] 0.01

##g
meansurf <- apply(Envginmat,2,mean)
sdsurf <- apply(Envginmat,2,sd)

val <- t((t(Envginmat)-meansurf)/sdsurf)
valint <- rowSums(val)

valg <- t((t(as.vector(g$ginhom))-meansurf)/sdsurf)
valgint <- sum(valg)

pvalueG = (sum(valgint>valgint)+1)/(length(valint)+1)

#ui #######################################################################################
ui <- fluidPage(
  titlePanel("Spatio-temporal dash board"),
  navbarPage("LIST",

             navbarMenu("Statistics",
                        tabPanel("Table",
                                 DT::dataTableOutput("table")
                        ),
                        tabPanel("Summary",
                                 verbatimTextOutput("summary"),
                                 verbatimTextOutput("str")
                        ),
                        tabPanel("Correlation",
                                 sidebarLayout(
                                   sidebarPanel(
                                     selectInput('xcol', 'X Variable', num_names[1:20]),
                                     selectInput('ycol', 'Y Variable', num_names[1:20],
                                                 selected= "crime_hour")
                                   ),
                                   mainPanel(
                                     verbatimTextOutput("cor"),
                                     plotOutput("corrplot")
                                   )
                                 )
                        )
             ),

             tabPanel("Plot",
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons("plotType", "Plot type",
                                       c("Histogram"="h", "Density plot"="d", "BoxPlot"="b", "BarPlot"="bp")
                          )
                        ),
                        mainPanel(
                          plotOutput("plot")
                        )
                      )
             ),

             navbarMenu("Intensity",
                        tabPanel("projection",
                                 plotOutput("projection",width = 800, height = 800)
                        ),
                        tabPanel("spatial",
                                 sidebarLayout(
                                   sliderInput(inputId = "nums", label="spatial bandwith",
                                               min=0, max=800, value=bwS),
                                   mainPanel(
                                     plotOutput("intensityS",width = 800, height = 800),
                                     verbatimTextOutput("statS")
                                   )
                                 )
                        ),
                        tabPanel("temporal",
                                 sidebarLayout(
                                   sliderInput(inputId = "numt", label="temporal bandwith",
                                               min=0, max=20, value=bwT),
                                   mainPanel(
                                     plotOutput("intensityT",width = 800, height = 800),
                                     verbatimTextOutput("statT")
                                   )
                                 )
                        )),

             navbarMenu("Intensity by Group",
                        tabPanel("crime type",
                                 sidebarLayout(
                                   sidebarPanel(
                                     radioButtons("STtype", "Spatial/Temporal",
                                                  c("spatial", "temporal"))
                                   ),
                                   mainPanel(
                                     plotOutput("crimetypePlot",width = 800, height = 800)
                                   )
                                 )
                        ),
                        tabPanel("year",
                                 sidebarLayout(
                                   sidebarPanel(
                                     radioButtons("STtype2", "Spatial/Temporal",
                                                  c("spatial", "temporal"))
                                   ),
                                   mainPanel(
                                     plotOutput("yearPlot",width = 800, height = 800)
                                   )
                                 )
                        ),
                        tabPanel("week",
                                 sidebarLayout(
                                   sidebarPanel(
                                     radioButtons("STtype3", "Spatial/Temporal",
                                                  c("spatial", "temporal"))
                                   ),
                                   mainPanel(
                                     plotOutput("weekPlot")
                                   )
                                 )
                        )
             ),

             navbarMenu("Correlation between pairs of points",
                        tabPanel("k-function",
                                 mainPanel(
                                   plotOutput("plotK",width = 800, height = 800),
                                   verbatimTextOutput("Kfunc")
                                 )
                        ),
                        tabPanel("g-function",
                                 mainPanel(
                                   plotOutput("plotG",width = 800, height = 800),
                                   verbatimTextOutput("Gfunc")
                                 )
                        ))

  )
)

#server###########################################################################

server <- function(input, output, session) {

  ##statistics ######################################
  selectedData <- reactive({
    obs_table[, c(input$xcol, input$ycol)]
  })

  output$table <- DT::renderDataTable({
    DT::datatable(obs_table)
  })

  output$summary <- renderPrint({
    summary(obs_table)
  })

  output$str <- renderPrint({
    str(obs_table)
  })

  output$cor <- renderPrint({
    cor(selectedData(), use="complete.obs")
  })

  output$corrplot <- renderPlot({
    corrplot(cor(obs_table[,num_names[1:20]], use="complete.obs"))
  })

  ##plot ########################################
  output$plot <- renderPlot({

    if(input$plotType == "h"){
      g = ggplot(obs_table, aes(x = crime_hour, fill=crime_type) )
      g = g+ geom_histogram(alpha=0.9)
      g + theme_bw(base_size = 20) + theme(legend.position = "bottom")
    }
    else if(input$plotType == "d"){
      g = ggplot(obs_table, aes(x = crime_hour, fill=crime_type) ) +
        geom_density(alpha = 0.8)
      g + theme_bw(base_size = 20) + theme(legend.position = "bottom")
    }
    else if(input$plotType == "b"){
      p1 = ggplot(obs_table, aes(x=crime_type, y=crime_hour)) + geom_boxplot()
      p2 = ggplot(obs_table, aes(x=crime_type, y=bank_dist)) + geom_boxplot()
      p1/p2
    }
    else if(input$plotType == "bp"){
      p1 = ggplot(obs_table, aes(x=crime_type)) + geom_bar()
      p2 = ggplot(obs_table, aes(x=crime_hour)) + geom_bar()
      p1/p2
    }
  })


  ##intensity #####################################
  #projection
  output$projection <- renderPlot({
    plot(valencia)
  })

  #spatial
  lambdaST <- reactive(density(valencia,lbw=input$nums,tbw=input$numt,dimyx=512,diggle=TRUE))
  output$intensityS <- renderPlot({
    plot(attr(lambdaST(),"netint"),main="Spatial intensity")
  })
  output$statS <- renderPrint({
    print(denS)
  })

  #temporal
  output$intensityT <- renderPlot({
    plot(attr(lambdaST(),"tempden"),main= "Temporal intensity",xlab="hour")
  })
  output$statT <- renderPrint({
    print(denT)
  })

  #intensity by group##################################
  #crime type
  output$crimetypePlot <- renderPlot({
    if (input$STtype=="spatial"){
      par(mfrow=c(2,2))
      plot(attr(lambda1,"netint"), main="'Agresion'")
      plot(attr(lambda2,"netint"), main="'AlarmasMujer'")
      plot(attr(lambda3,"netint"), main="'Otros'")
      plot(attr(lambda4,"netint"), main="'Sustraccion'")
    }
    else {
      par(mfrow=c(2,2))
      plot(attr(lambda1,"tempden"), main="'Agresion'")
      plot(attr(lambda2,"tempden"), main="'AlarmasMujer'")
      plot(attr(lambda3,"tempden"), main="'Otros'")
      plot(attr(lambda4,"tempden"), main="'Sustraccion'")
    }
  })

  #year
  output$yearPlot <- renderPlot({
    if (input$STtype2=="spatial"){
      par(mar=c(3,3,3,3))
      par(mfrow=c(3,3))
      for (i in 2:10){
        plot(attr(year_lambda_list[[i]],"netint"), main=year_names[i],zlim=c(0,1.5e-4))
      }
    }
    else {
      par(mar=c(3,3,3,3))
      par(mfrow=c(3,3))
      for (i in 2:10){
        plot(attr(year_lambda_list[[i]],"tempden"), main=year_names[i],ylim=c(0,0.06))
      }
    }
  })

  #week
  output$weekPlot <- renderPlot({
    if (input$STtype3=="spatial"){
      par(mar=c(3,3,3,3))
      par(mfrow=c(3,3))
      for (i in 1:7){
        plot(attr(lambda_list[[i]],"netint"), main=week_names[i],zlim=c(0,2e-4))
      }
    }
    else {
      par(mar=c(3,3,3,3))
      par(mfrow=c(3,3))
      for (i in 1:7){
        plot(attr(lambda_list[[i]],"tempden"), main=week_names[i],ylim=c(0,0.06))
      }
    }
  })

  ##correlation between ###########################################
  #kfunc
  output$plotK <- renderPlot({
    par(mfrow=c(1,1))
    persp3D(k1$r,k1$t,k1$Kinhom/1000,theta=30,phi=10,
            facets=FALSE,col=1,ticktype= "detailed",
            xlab="r = distance",ylab="t = time",zlab="",
            main=expression(italic({hat(K)[LI]^{ST}}(r,t))),
            nticks=6,cex.axis=1.5,cex.lab=2, cex.main=2,
            zlim= range(c(min(minEnv,k1$Kinhom),max(maxEnv,k1$Kinhom)))/1000,
            resfac=2)

    persp3D(k1$r,k1$t,minEnv/1000,col = "gray80",add=T,facets=T,
            ticktype= "detailed",
            nticks=6,cex.axis=1.3,cex.lab=1.7,resfac=2,border = 1)

    persp3D(k1$r,k1$t,maxEnv/1000,col = "gray80",add=T,facets=T,
            ticktype= "detailed",
            nticks=6,cex.axis=1.3,cex.lab=1.7,resfac=2,border = 1)
  })
  output$Kfunc <- renderPrint({
    pvalueK
  })

  #gfunc
  output$plotG <- renderPlot({
    par(mfrow=c(1,1))
    persp3D(g$r,g$t,g$ginhom,theta=35,phi=10,
            facets=FALSE,col=1,ticktype= "detailed",
            xlab="r = distance",ylab="t = time",zlab="",
            main=expression(italic({hat(g)[LI]^{ST}}(r,t))),
            nticks=6,cex.axis=1.5,cex.lab=2, cex.main=2,
            zlim= range(c(min(minEnvg,g$ginhom),max(maxEnvg,g$ginhom))),
            resfac=2)
    persp3D(g$r,g$t,minEnvg,col = "gray",add=T,facets=T,
            ticktype= "detailed",nticks=6,cex.axis=1.3,cex.lab=1.7,
            resfac=2,border=1)

    persp3D(g$r,g$t,maxEnvg,col = "gray",add=T,facets=T,
            ticktype= "detailed",nticks=6,cex.axis=1.3,cex.lab=1.7,
            resfac=2,border = 1)
  })
  output$Gfunc <- renderPrint({
    pvalueG
  })

}


shinyApp(ui=ui, server=server)
