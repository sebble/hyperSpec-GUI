library(gWidgets)
library(hyperSpec)
options("guiToolkit"="RGtk2")

### Notes
###  * plot.new() fails for new gui after axis cut in previous and not cleared
 
plotspc.gui <- function(data) {
 
    availDists <- c(Normal="rnorm", Exponential="rexp")
 
 
    updateData <- function(...) {
 
      if(dim(data[[]])[1]<10){
        tmpdata <<- data 
      }else{
        sSize = svalue(sampleSize)
        choose <- sample(1:(dim(data[[]])[1]),size=sSize)
        tmpdata <<- data[choose]
        updatePlot(...)
      }
    }
 
    updatePlot <- function(...) {
 
      #plotspc(paracetamol,wl.range=c(min~1000,1200~1900,2800~max),xoffset=c(100,800))

#        clips <<- sort(clips)  

        # switch to starts (the and ends, these two holding indices!
        ## this should also allow to pre-cut the spectra axis analog to the aggregation. Again this
        ## is something that anyways occurs in the plotspc "preprocessing".
        starts <<- sort (unique (starts))
        ends <<- sort (unique (ends))
        
        range <- mapply (`:`, starts, ends)
        if (length (starts) > 1){
          offset <- tmpdata@wavelength [tail (starts, -1)] -
                    tmpdata@wavelength [head (ends,   -1)]
          offset <- offset - (diff (range (tmpdata@wavelength)) - sum (offset)) / 10
        } else {
          offset <- 0
        }

        range <- mapply (`:`, starts, ends)
        ## for (i in sequence (length (clips) / 2)){ # trap: if length (clips) == 0 => 1 : length (clips) => 1 : 0 !!!
        ##   if(i==1){
        ##     range <- expression (min ~ clips [1])
        ##   }
        ##   if(i==length(clips) / 2){
        ##     range <- c(range,clips[2*i]~max)
        ##   }
        ##   else{
        ##     range <- c(range,clips[2*i-1]~clips[2*i])
        ##   }
        ##   offset <- c(offset,clips[2*i]-clips[2*i-1]+5)
        ## }
        
        if(offset != 0){
          #cutRedraw <<- T
          enabled(cutPlotBtn1) <- F
        }

        plotspc(tmpdata, wl.reverse=svalue(reverseWL), wl.range = range, xoffset = offset, wl.index = TRUE)
        ### If I have ", ...)" above then 'bandwidthAdjust' fails
    }
 
    updateZoom <- function(h, ...) {
 
        print(h)
        ### If I have ", ...)" above then 'bandwidthAdjust' fails
    }
 
    updateCutPlot <- function(h, ...) {

      #if(cutRedraw){
      #  cat("Notice: Cannot cut after redraw, select clear and try again.\n")
      #  return (F)
      #}
      
        pos <- locator(1)
        clips <- pos$x
        abline(v=pos$x,col=2)
        pos <- locator(1) # locator doesn't  work on offset axes
        abline(v=pos$x,col=2)
        clips <- c (clips, pos$x)

        starts <<- c (starts, wl2i (tmpdata, max (clips))) # may need to become data !?
        ends   <<- c (ends,   wl2i (tmpdata, min (clips))) # may need to become data !?

        ## I guess the cutting pre-processing should rather go here so it is not redone every time we
        ## plot: that means the global objects are range and offset rather than starts and ends.  It
        ## has also the advantage that the global variables then correspond better to plotspc's
        ## arguments (which I take to be a sign of good programming)
        
        #updatePlot()
    }
    updateCutPlotClear <- function(h, ...) {
      
      #cutRedraw <<- F
          enabled(cutPlotBtn1) <- T
      
      starts <<- 1
      ends <<- nwl (data)
                                        # clips <<- array()
      updatePlot()
    }
 
    ## gWidgets for Data
    sampleSize <- gradio(c(1,5,10,20,50), handler = updateData)
    reSample <- gbutton("Resample", handler = updateData)
    cutPlotBtn1 <- gbutton("Cut Axis", handler = updateCutPlot)
    cutPlotBtn2 <- gbutton("Redraw", handler = updatePlot)
    cutPlotBtn3 <- gbutton("Clear", handler = updateCutPlotClear)
    ## gWidgets for Plot
    reverseWL <- gcheckbox("Reverse Wavelength", handler = updatePlot)
    bandwidthAdjust <- gslider(from=0,to=2,by=.01,
                               value=1,handler=updatePlot)
 
    ## create GUI window
    ###  +------window-----+
    ###  |+-----BigGrp----+|
    ###  ||+-grp-++-ggfx-+||
    ###  |||  gw ||      |||
    ###  |||  gw ||      |||
    ###  ||+-----++------+||
    ###  |+-----BigGrp----+|
    ###  +------window-----+
    window <- gwindow("plotspc - hyperSpec GUI")
    BigGroup <- ggroup(cont=window)
    group <- ggroup(horizontal=FALSE, container=BigGroup)
                               
    tmp <- gframe("Sample size", container=group)
    add(tmp,sampleSize)
    add(tmp,reSample)
                               
    add(gframe("Some slider", container=group),bandwidthAdjust,
        expand=T)

    add(gframe("Plotting options", container=group),reverseWL)
    tmp <- gframe("Plotting", container=group)
    add(tmp, cutPlotBtn1)
    add(tmp, cutPlotBtn2)
    add(tmp, cutPlotBtn3)
    add(BigGroup, ggraphics(), handler=updateZoom)
 
 
    ## Init data for drawing plot.
    tmpdata <- data; #initialize to data, so everything works from the beginning
    starts <- 1
    ends <- nwl (data)
    #cutRedraw <- F
    updateData(0); ### what is h here?
    updatePlot(0);
}
