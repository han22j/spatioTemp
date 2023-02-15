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


# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Covid-19 Incidence Intensity Map"),

  column(4, wellPanel(
    dateInput('date2',
            label = paste('Date input:',
                          'yyyy-mm-dd format, range limited,',
                          'the first incidence is 2020-03-02'),
            value = "2020-03-02",
            min = "2020-03-02", max = "2020-07-20",
            format = "yyyy-mm-dd",
            datesdisabled = c( "2020-04-13", "2020-05-11", "2020-06-08", "2020-07-06"),
            startview = 'year'),
    sliderInput(inputId="sval", label="space parameter", min=0.0001, max=0.01, value=0.006)
  )),

  column(6,
  mainPanel(plotOutput(outputId = "map"))
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$map <- renderPlot({
    map = read_sf("./data/borde_cali.shp")
    X.full <- read.table("./data/Xfull.csv")
    X.full$observe.time <- round(28*(X.full$observe.time))

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
    intensity_mat <- rep(NA,n.S)

    wanttime <- as.numeric(difftime(input$date2, as.Date("2020-03-02"), units = "days"))
    wantX <- X.full[X.full$observe.time==wanttime , ]

    n <- nrow(wantX)
    density_S <- matrix(NA,n.S,n)

    space.range <- input$sval

    for(i in 1:n){
      density_S[,i] <- dnorm(S[,1],wantX[i,1],space.range)*
        dnorm(S[,2],wantX[i,2],space.range)
    }

    intensity_mat <- apply(density_S, 1, sum, na.rm=T)
    intensity_df <- data.frame(cbind(x=S[,1],y=S[,2], intensity = intensity_mat))
    colnames(intensity_df)[3] <- "intensity"
    intensity_df$intensity <-  intensity_df$intensity/max(intensity_df$intensity)


    ggplot(intensity_df)+
      geom_point(aes(x, y, color = intensity), size = 1) +
      #geom_point(data=map, aes(x=x, y=y), size=0.5) +
      coord_fixed(ratio = 1) +
      scale_color_gradientn(colours = c("blue",  "orange"),
                            #values = c(0, 1),
                            limits = c(0, 1)) +
      #geom_point(data=data.frame(X), aes(Latitud,Longitud), size=1, shape=20, col="black") +
      ggtitle(as.Date(wanttime, origin = "2020-03-02", format = "%Y-%m-%d")) + theme_bw() +
      labs(x="longitude", y="latitude") +
      theme(plot.title = element_text(hjust = 0.5, face="bold"))
    })

}

# Run the application
shinyApp(ui = ui, server = server)
