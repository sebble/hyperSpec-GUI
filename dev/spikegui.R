library (gWidgetsRGtk2)
lab <- ggraphics("drag me",container=gwindow())
lab
dev.list ()
plot (1:3)

locator ()


lab <- ggraphics("drag me",container=gwindow())
addhandlerclicked(lab, function (h, ...) {
  tag (h$obj, "startx") <- h$x
  tag (h$obj, "starty") <- h$y
  })
plot (1:3)

adddropsource(lab, type = "object")
adddroptarget(lab, type = "object")
adddropmotion (lab, handler = function (h, ...) {
  browser ()
print (h$dropdata)
cat (".............................\n")
}) 


adddropsource(lab)
adddroptarget(gd)
adddropmotion(gf,handler=function(h,...) {
  tag (lab, "start") <- as.numeric (c( h$x, h$y))
  plot (flu)
})

adf = gdf(mtcars, cont = TRUE)
gd = ggraphics(cont = TRUE)
plotData = function() {
  dropvalue = tag(gd,"value")
  theValues = svalue(dropvalue)
  theName = id(dropvalue)
  hist(theValues, xlab=theName, main="")
}

ids = adddroptarget(gd,targetType="object", handler = function(h,...) {
    tag(gd, "value") <- h$dropdata
    plotData()

    if(is.gdataframecolumn(h$dropdata)) {
      view.col = h$dropdata
      id = addhandlerchanged(view.col, handler=function(h,...) plotData())
    }
})
