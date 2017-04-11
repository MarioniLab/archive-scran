exploreData <- function(m, pD, fD, reDim, run=TRUE) 
# This function creates an interactive shiny app to explore expression data.
#
# created by Karsten Bach
# written 11 April 2017
# last modified 11 April 2017
{
    # Checks
    if (!identical(rownames(m),rownames(fD))) {
	stop("m and fD need to have identical rownames")
    }
    if (!identical(colnames(m),rownames(pD))) {
	stop("rownames of pD need to correspond to colnames of m")
    }
    # Build data
    covariates <- colnames(pD)
    pD$Dim1 <- reDim[,1]
    pD$Dim2 <- reDim[,2]

    # Set up shiny
    ui <- fluidPage(
	      titlePanel("Explorative Analysis"),
	      sidebarLayout(

		  sidebarPanel(
		      inputPanel(
			  selectInput("colorBy", label = "Color By",
					   choices= covariates, selected=covariates[1]),
			  selectInput("groupBy", label = "Group Boxplot By",
					   choices= covariates, selected=covariates[1]) 
			  ),
		      plotOutput("distPlot1")
		      ),
		  mainPanel(
		      tabsetPanel(tabPanel("Plots",splitLayout(plotOutput("tSNE"),
							     plotOutput("tSNE2")),
			dataTableOutput("table")))
		      )
		  )
	      )

    server <- function(input, output) {
	# Load the gene level data
	output$table <- renderDataTable({
	    out <- fD
	    datatable(out, filter="top", selection=list(mode="single",
							    selected=1))
	    })

	# tSNE plot colored by covariates
	output$tSNE <- renderPlot({
	    tsnPlot <- ggplot(pD, aes_string(x="Dim1", y="Dim2", color=input$colorBy)) +
		geom_point(size=1.5) +
		theme_void()
	    tsnPlot
	    })

	# tSNE plot colored by gene expression
	output$tSNE2 <- renderPlot({
	    s <- input$table_rows_selected
	    if (!is.null(s)) {
		gene <- rownames(fD)[s]
		expr <- log2(t(m)[,gene]+1)
		pD[,gene] <- expr
		pD <- pD[base::order(pD[,gene]),]
		tsnPlot <- ggplot(pD, aes_string(x="Dim1", y="Dim2", color=gene)) +
		    geom_point(size=1.5) +
		    scale_color_viridis() +
		    theme_void()
		tsnPlot
		}
	    })

	# Distribution of gene expression by covariate
	output$distPlot1 <- renderPlot({
	    s <- input$table_rows_selected
	    if (!is.null(s)) {
		gene <- rownames(fD)[s]
		exps <- log2(t(m)[,gene]+1)
		pltDat <- pD
		pltDat$value <- exps
		plt <- ggplot(pltDat, aes_string(x=input$groupBy,y="value")) +
		    geom_boxplot() +
		    geom_point(position="jitter",alpha=0.2,shape=19) +
		    ggtitle(gene) +
		    theme_bw()
		plt
		}
	    })
	}

    app <- shinyApp(ui, server)
       if (run) {
	  return(runApp(app))
       } else {
	  return(app)
       }
}
