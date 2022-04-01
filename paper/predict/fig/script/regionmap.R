## see https://soprasteriaanalytics.se/2018/08/15/visualization-of-data-with-maps/

library(raster)
library(rworldmap)

main <- function(save=FALSE) {
    SwedenLevel1 = raster::getData("GADM", country = "Sweden", level = 1)
    ## rworldmap::mapPolys(SwedenLevel1, nameColumnToPlot="NAME_1", addLegend=TRUE)
    ## text(xx$colourVector,xx$NAME_1, cex=0.7)

    library(R.matlab)
    bnumber <- function(x){format(round(as.numeric(x), 1), nsmall=0, big.mark=" ")}

    x <- readMat("../../../../data/Ncounties.mat")
    N <- apply(x$N,2,sum)
                                        # This is the ordering from the population N
    Nnames <- c("Stockholm", "Uppsala", "Södermanland", "Östergötland" ,
                "Jönköping", "Kronoberg", "Kalmar", "Gotland", "Blekinge",
                "Skåne", "Halland", "Västra Götaland", "Värmland",
                "Örebro", "Västmanland", "Dalarna", "Gävleborg",
                "Västernorrland", "Jämtland", "Västerbotten", "Norrbotten")
    plotname <- SwedenLevel1$NAME_1
    ## found problem with Örebro
    elemOrebro <- which(plotname == "Orebro")
    plotname[elemOrebro] <- gsub("O","Ö",plotname[elemOrebro])
    for (i in seq_len(length(plotname))) {
        ## find
        elem <- which(Nnames %in% plotname[i])
        plotname[i] <- paste(plotname[i], bnumber(N[elem]),sep=": ")
    }


    SwedenLevel1$NAME_1 <- plotname


    colors <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA",
                "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744",
                "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77",
                "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455",
                "#DD7788")

    doplot <- function(x){rworldmap::mapPolys(x,
                                              nameColumnToPlot="NAME_1",
                                              addLegend=TRUE,
                                              colourPalette=colors,
                                              borderCol="white",
                                              mapTitle="",xlim=c(2,24))}

    if (save) {
        pdf("../swedenmap.pdf")
        doplot(SwedenLevel1)
        dev.off()
    } else {
        doplot(SwedenLevel1)
    }
}
