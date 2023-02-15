library(shiny)
library(rsconnect)
library(ggplot2)
library(gridExtra)
library(stpp)
library(fields)
library(spatstat)
library(mvtnorm)
library(ks)
library(sf)
library(maptools)
library(spatialrisk)
library(gmt)
library(terra)
library(geosphere)
library(tidyverse)


setwd("/Users/seorimy/Dropbox/Shiny/shinyio/param2")

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Covid-19 Hotspots Intensity Map"),
  
  column(4, wellPanel(
    sliderInput(inputId="sval", label="Space Parameter", min=0, max=0.01, value=0.006),
    sliderInput(inputId="tval", label="Time Parameter", min=0, max=0.25, value=0.2),
    numericInput(inputId="rval", label=paste("Time Window:",
                                             "each time window has a length of 4 weeks,",
                                             "the first date is 2020-03-02"),
                 min=1, max=8, value=1)
    )), 
  
  column(6, mainPanel(plotOutput(outputId = "map")))
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  map = read_sf("borde_cali.shp")
  
  timewindow_end <- as.Date(c("2020-03-30","2020-04-27","2020-05-25","2020-06-22","2020-07-20"))
  timewindow_start <- as.Date(c("2020-03-02","2020-03-31","2020-04-28","2020-05-26","2020-06-23"))

  rangeparam <- function(spaceval, timeval, r){
    load("test.RData")
    
    g <- g.list[[r]]
    CC <- CC.list[[r]]
    X <- X.list[[r]]
    q <- probs.g.list[[r]]
    
    Cind <- as.numeric( names(table(g)) )
    par <- as.data.frame(CC[Cind,])
    rownames(par) <- 1:nrow(par)
    par$ind <- 1:nrow(par)
    
    space.range <- spaceval
    time.range <- timeval
    
    X <- data.frame(X)
    
    xx <- st_coordinates(map)[,1]
    yy <- st_coordinates(map)[,2]
    grid.original <- expand.grid(seq(min(xx),max(xx),,200),seq(min(yy),max(yy),,200),seq(0,1,,20))
    colnames(grid.original) <- c('X','Y','T')
    
    pnts_sf <- st_as_sf(grid.original, coords = c('X', 'Y'), crs = st_crs(map))
    map_sf <- st_as_sf(map, coords=c("X","Y"), crs=st_crs(map))
    in.val <- st_intersects(pnts_sf, map_sf)
    inval <- lengths(in.val)>0
    
    S <- grid.original[inval,]
    colnames(S) <- c("x", "y", "t")
    
    n.S=nrow(S)
    intensity_mat4 <- rep(NA,n.S)
    n <- nrow(par)
    density_S <- matrix(NA,n.S,n)
    
    for(i in 1:n){
      density_S[,i] <- dnorm(S[,1],par[i,1],space.range)*dnorm(S[,2],par[i,2],space.range)*dnorm(S[,3],par[i,3],time.range)
    }

    
    intensity_mat4 = apply(density_S, 1, sum, na.rm=T)
    intensity_df <- data.frame(cbind(x=S[,1],y=S[,2], intensity = intensity_mat4))
    intensity_df <- as.data.frame(intensity_df %>% group_by(x, y) %>% summarise(int = sum(intensity)))
    intensity_df$intensity <- intensity_df$int/max(intensity_df$int)
    
    p <-  ggplot(intensity_df)+
      geom_point(aes(x, y, color = intensity, fill = intensity), size = 1) +
      coord_fixed(ratio = 1) +
      scale_color_gradientn(colours = c("blue",  "orange"),
                            #values = c(0, 1),
                            limits = c(0, 1)) +
      scale_fill_gradientn(colours = c("blue",  "orange"),
                           #values = c(0, 1),
                           limits = c(0, 1),
                           guide="none") +
      ggtitle(paste(timewindow_start[r], "~", timewindow_end[r])) + theme_bw() +
      labs(x="longitude", y="latitude") +
      theme(plot.title = element_text(hjust = 0.5, face="bold")) 
    
    print(p)
  }
  
  
  output$map <- renderPlot({
    rangeparam(input$sval, input$tval, input$rval)
  })
}

shinyApp(ui = ui, server = server)
