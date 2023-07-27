# ####################################
# Global declarations for the analysis
# ####################################


#----------------------------------------
# Datatable default buttons 
#----------------------------------------

# Configure a set of buttons for exporting data from datatables
# (NOTE: an orthogonal named 'export' must be defined in the concerned columns for proper formatting)
exportButtonsListDT = list( 
  list( # Collection containing export buttons with altered 'exportOptions' to use non-formatted/non-ordered data and save full table
    extend = "collection",
    text = "<center>Export<br>Full</center>",
    buttons = list( list( extend  = "excel",
                          exportOptions = list( modifier = list(order    = 'index', # give result in the same order as original table
                                                                page     = 'all',   # do not restrict to currently displayed page
                                                                search   = 'none'), # do not restrict to current filtered rows
                                                # do not define 'selected' here as 'undefined' gives expected behavior (selection ignored, not sure how to set to null or undefined using R)
                                                orthogonal = 'export')), # use 'export' orthogonal source (no formatting, defined in 'columnsDef')
                    list( extend  = "csv",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'none'),
                                                orthogonal = 'export')),
                    list( extend  = "copy",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'none'),
                                                orthogonal = 'export')))),
  list( # Collection containing export buttons with altered 'exportOptions' to use non-formatted/non-ordered data and save filtered table
    extend = "collection",
    text = "<center>Export<br>Filtered</center>",
    buttons = list( list( extend  = "excel",
                          exportOptions = list( modifier = list(order    = 'index',      # give result in the same order as original table
                                                                page     = 'all',        # do not restrict to currently displayed page
                                                                search   = 'applied'),   # restrict to current filtered rows
                                                # do not define 'selected' here as 'undefined' gives expected behavior (selection ignored, not sure how to set to null or undefined using R)
                                                orthogonal = 'export')),
                    list( extend  = "csv",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'applied'),
                                                orthogonal = 'export')),
                    list( extend  = "copy",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'applied'),
                                                orthogonal = 'export')))),
  list( # Collection containing export buttons with altered 'exportOptions' to use non-formatted/non-ordered data and save rows from current page only
    extend = "collection",
    text = "<center>Export<br>Page</center>",
    buttons = list( list( extend  = "excel",
                          exportOptions = list( modifier = list(order    = 'index',      # give result in the same order as original table
                                                                page     = 'current',    # restrict to currently displayed page
                                                                search   = 'none'),      # do not restrict to current filtered rows
                                                # do not define 'selected' here as 'undefined' gives expected behavior (selection ignored, not sure how to set to null or undefined using R)
                                                orthogonal = 'export')),
                    list( extend  = "csv",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'current',
                                                                search   = 'none'),
                                                orthogonal = 'export')),
                    list( extend  = "copy",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'current',
                                                                search   = 'none'),
                                                orthogonal = 'export')))),
  list( # Collection containing export buttons with altered 'exportOptions' to use non-formatted/non-ordered data and save selected rows only (independently of the filtering)
    extend = "collection",
    text = "<center>Export<br>Selected</center>",
    buttons = list( list( extend  = "excel",
                          exportOptions = list( modifier = list(order    = 'index',      # give result in the same order as original table
                                                                page     = 'all',        # do not restrict to currently displayed page
                                                                search   = 'none',       # do not restrict to current filtered rows
                                                                selected = TRUE),        # extension Select
                                                orthogonal = 'export')),
                    list( extend  = "csv",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'none',
                                                                selected = TRUE),
                                                orthogonal = 'export')),
                    list( extend  = "copy",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'none',
                                                                selected = TRUE),
                                                orthogonal = 'export')))),
  # Print buttons
  list( # Collection containing print buttons with altered 'exportOptions' to render visible columns only in displayed rows order, and select desired rows (see comments in export buttons, here it uses default orthogonal that gives formatted data)
    extend = "collection",
    text = "<center>Print<br>(as seen)</center>",
    buttons = list( list( extend  = "print",
                          text = "Full",
                          exportOptions = list( modifier = list(order    = 'current',
                                                                page     = 'all',
                                                                search   = 'none'),
                                                columns = ':visible')),
                    list( extend  = "print",
                          text = "Filtered",
                          exportOptions = list( modifier = list(order    = 'current',
                                                                page     = 'all',
                                                                search   = 'applied'),
                                                columns = ':visible')),
                    list( extend  = "print",
                          text = "Page",
                          exportOptions = list( modifier = list(order    = 'current',
                                                                page     = 'current',
                                                                search   = 'none'),
                                                columns = ':visible')),
                    list( extend  = "print",
                          text = "Selected",
                          exportOptions = list( modifier = list(order    = 'current',
                                                                page     = 'all',
                                                                search   = 'none',
                                                                selected = TRUE),
                                                columns = ':visible')))),
  list( extend = "selectNone",
        text = "<center>Deselect<br>all</center>"));

#----------------------------------------
# Violin/Jitter plot using plotly
#----------------------------------------

# Function that creates a violin plot with jitter points
# To be used for showing #UMIs, #Genes, and %Mito
plotViolinJitter = function(plotData, xAxisFormula, yAxisFormula, colorFormula = NULL, pointSize = 3, pointsColorFormula = "#444444", pointsOpacity = 0.5, xAxisTitle = "", yAxisTitle = "", hoverText = NULL, traceName = "", xTicklabels = c(""), panelWidth = 300, panelHeight = 600, thresholdHigh = NULL, thresholdLow = NULL)
{
  # Generate a jitter that will not be recomputed at each redraw
  constantJitter = jitter(rep(0, nrow(plotData)), amount = 0.15);
  
  # Add it to data so it can be accessed with updated formula (markers)
  plotData = cbind(plotData, jitterCoords = constantJitter + 0.17);
  
  p = plot_ly( plotData,
               x = xAxisFormula,
               y = yAxisFormula,
               name = traceName,
               fillcolor = colorFormula, # 'color' argument does not work with custom values (?)
               height = panelHeight,
               width = panelWidth) %>%
    add_trace( type='violin',
               showlegend=FALSE,
               points = FALSE,
               line = list( width = 1.5, 
                            color = "#44444444"), # impossible to specify individual colors (?)
               meanline = list(visible = TRUE),
               scalemode = 'width',
               width = 0.6,
               side = 'negative',
               hoveron = 'kde',
               hoverinfo = 'skip') %>%
    add_trace( type = "scatter",
               mode = "markers",
               x = update(xAxisFormula, ~ . + jitterCoords), # Update formula adding jitter values
               marker = list(size = pointSize,
                             color = pointsColorFormula,
                             opacity = pointsOpacity),
               text = hoverText,
               hoverinfo = 'text+name') %>%
    layout( yaxis = list( title = yAxisTitle,
                          rangemode = "nonnegative"),
            xaxis = list( title = xAxisTitle,
                          ticktext = xTicklabels,
                          tickvals = ((1:length( xTicklabels)) + 0.01)),
            margin = list( l = 0, r = 0, t = 0, b = 0, pad = 5));
  
  
  # Function creating a horizontal line as a shape to be added to layout (use 'paper' coords ref for x axis)
  hline = function(y, color = "blue") 
  {
    if (is.null( y)) return( NULL);
    list( list( type = "line", 
                x0 = 0, 
                x1 = 1, 
                xref = "paper",
                y0 = y, 
                y1 = y, 
                line = list( color = color)))
  }
  
  # Shapes to be added to layout (for eventual threshold lines)
  shapes = c( hline( thresholdHigh, "red"), hline( thresholdLow, "blue"));
  
  if(length( shapes)) return(p %>% layout( shapes = shapes));
  return(p);
}