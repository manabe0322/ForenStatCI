ForenStatCI <- function(){
  Version <- "ForenStatCI v.1.0.0"

  if(require(tcltk)){
    print("tcltk is loaded correctly")
  }else{
    print("trying to install tcltk...")
    install.packages("tcltk")
    if(require("tcltk")){
      print("tcltk installed and loaded")
    }else{
      stop("could not install tcltk")
    }
  }
  if(require(tcltk2)){
    print("tcltk2 is loaded correctly")
  }else{
    print("trying to install tcltk2...")
    install.packages("tcltk2")
    if(require("tcltk2")){
      print("tcltk2 installed and loaded")
    }else{
      stop("could not install tcltk2")
    }
  }
  if(require(MCMCpack)){
    print("MCMCpack is loaded correctly")
  }else{
    print("trying to install MCMCpack...")
    install.packages("MCMCpack")
    if(require("MCMCpack")){
      print("MCMCpack installed and loaded")
    }else{
      stop("could not install MCMCpack")
    }
  }
  if(require(gtools)){
    print("gtools is loaded correctly")
  }else{
    print("trying to install gtools...")
    install.packages("gtools")
    if(require("gtools")){
      print("gtools installed and loaded")
    }else{
      stop("could not install gtools")
    }
  }

  Sapplywhich <- function(Rightvalue, Leftvec){
    return(which(Leftvec == Rightvalue))
  }

  OpenFile <- function(fp, var, top, filestate, Cursor){
    FileName <- tclvalue(tkgetOpenFile(parent = top, initialdir = tclvalue(fp), 
                         multiple = "true", filetypes = "{{CSV Files} {.csv}}"))
    if(!nchar(FileName)){
      tkmessageBox(message = "No file was selected!", icon = "error", type = "ok")
    }else{
      tmp <- sub("\\}", FileName, replacement = "")
      tmp2 <- sub("\\{", tmp, replacement = "")
      tclvalue(fp) <- tmp2
      foo3 <- strsplit(tmp2, "/")[[1]]
      tclvalue(var) <- strsplit(foo3[length(foo3)], "\\.csv")[[1]][1]
      tclvalue(filestate) <- "normal"
      tclvalue(Cursor) <- "hand2"
    }
  }

  Countdata.make <- function(Countdata.filepath){
    Countdata.File <- read.csv(Countdata.filepath, row.names = 1, header = TRUE)
    Countdata.File <- as.matrix(Countdata.File)
    Locus.names <- colnames(Countdata.File)
    NL <- length(Locus.names)
    Allele.count <- list()
    for(i in 1:NL){
      Allele.count[[i]] <- Countdata.File[!is.na(Countdata.File[, i]), i]
    }
    return(list(Locus.names, Allele.count))
  }

  profileCheck <- function(Countdata.L, Profile, Profile.L){
    if(setequal(Countdata.L, Profile.L)){
      return(Profile[sapply(Countdata.L, Sapplywhich, Leftvec = Profile.L), ])
    }else{
      return(NULL)
    }
  }

  CIrange.calc <- function(CIparcent){
    Quantile1 <- (100 - CIparcent) / 200
    return(c(Quantile1, 1 - Quantile1))
  }

  AlleleCountRemake <- function(Allele.count.oneL, Profiles.oneL){
    AlleleName <- as.numeric(names(Allele.count.oneL))
    Profiles.oneL <- sort(unique(Profiles.oneL))
    NonObsPos <- which(is.element(Profiles.oneL, AlleleName) == FALSE)
    if(length(NonObsPos) > 0){
      AlleleName <- c(AlleleName, Profiles.oneL[NonObsPos])
      Allele.count.oneL <- c(Allele.count.oneL, rep(0, length(NonObsPos)))
      names(Allele.count.oneL) <- AlleleName
      return(Allele.count.oneL)
    }
    return(Allele.count.oneL)
  }

  RandAlFreq <- function(Allele.count.oneL, Nsimu, prior){
    if(prior == "1"){
      Randfreq <- rdirichlet(Nsimu, Allele.count.oneL + 1)
    }else{
      Randfreq <- rdirichlet(Nsimu, Allele.count.oneL + 1 / length(Allele.count.oneL))
    }
    return(Randfreq)
  }

  ExpectAlFreq <- function(Allele.count.oneL, prior){
    if(prior == "1"){
      Expect <- (Allele.count.oneL + 1) / sum(Allele.count.oneL + 1)
    }else{
      Expect <- (Allele.count.oneL + 1 / length(Allele.count.oneL)) / (sum(Allele.count.oneL) + 1)
    }
    return(Expect)
  }

  Freq <- function(){
    Tab1.Freq <- function(){
      tkdestroy(frame1.Freq)
      frame1.Freq <<- tkframe(tab1.Freq)
      countdata.Freq.label <- tklabel(frame1.Freq, text = "Allele count data")
      dataname.Freq.label <- tklabel(frame1.Freq, textvariable = Countdata.Freq.var, width = 30, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      countdata.Freq.butt <- tkbutton(frame1.Freq, text = "    Load    ", cursor = "hand2", command = function() OpenFile(Countdata.Freq.filepath, Countdata.Freq.var, tf.Freq, Countdata.Freq.butt.state, Countdata.Freq.butt.cursor))
      tkgrid(countdata.Freq.label, dataname.Freq.label, countdata.Freq.butt, padx = 20, pady = 20, sticky = "w")
      estimate.butt <- tkbutton(frame1.Freq, text = "    Estimate    ", cursor = "hand2", command = function() Tab2.Freq())
      tkgrid(estimate.butt, padx = 20, pady = 20, sticky = "w")
      tkgrid(frame1.Freq)
    }

    Tab2.Freq <- function(){
      tkdestroy(frame2.Freq)
      frame2.Freq <<- tkframe(tab2.Freq)
      if(tclvalue(Countdata.Freq.filepath) != ""){
        Countdata <- Countdata.make(tclvalue(Countdata.Freq.filepath))

        selectlocus.butt <- tkbutton(frame2.Freq, text = "    Select locus    ", cursor = "hand2", command = function() Tab2.Freq_Locus(Countdata[[1]]))
        selectlocus.label <- tklabel(frame2.Freq, textvariable = Locus.var, width = 30, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
        selectallele.butt <- tkbutton(frame2.Freq, text = "    Select allele    ", cursor = "hand2", command = function() Tab2.Freq_Allele(Countdata))
        selectallele.label <- tklabel(frame2.Freq, textvariable = Allele.var, width = 30, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
        tkgrid(selectlocus.butt, selectlocus.label, padx = 20, pady = 20, sticky = "w")
        tkgrid(selectallele.butt, selectallele.label, padx = 20, pady = 20, sticky = "w")

        graphDL.butt <- tkbutton(frame2.Freq, text = "    Graph    ", cursor = "hand2", command = function() graphDL.make(Countdata))
        exportFreq.butt <- tkbutton(frame2.Freq, text = "    Export expected frequencies    ", cursor = "hand2", command = function() exportFreq.make(Countdata))
        tkgrid(graphDL.butt, exportFreq.butt, padx = 20, pady = 20, sticky = "w")
        tk2notetab.select(tabs.Freq, "Results")
      }
      tkgrid(frame2.Freq)
    }

    exportFreq.make <- function(Countdata){
      NL <- length(Countdata[[2]])
      allAlleles <- sort(unique(as.numeric(unlist(sapply(Countdata[[2]], names)))))
      freqTable <- matrix("", length(allAlleles), NL)
      rownames(freqTable) <- allAlleles
      colnames(freqTable) <- Countdata[[1]]
      for(i in 1:NL){
        Countdata.oneL <- Countdata[[2]][[i]]
        freqTable[which(is.element(allAlleles, names(Countdata.oneL)) == TRUE), i] <- (Countdata.oneL + 1) / sum(Countdata.oneL + 1)
      }
      Save.as <- tkgetSaveFile(filetypes = "{{CSV Files} {.csv}}")
      if(tclvalue(Save.as) != ""){
        if(substr(tclvalue(Save.as), nchar(tclvalue(Save.as)) - 3, nchar(tclvalue(Save.as))) == ".csv"){
          ReportName <- tclvalue(Save.as)
        }else{
          ReportName <- paste(tclvalue(Save.as), ".csv", sep = "")
        }
        write.csv(freqTable, file = ReportName, row.names = TRUE)
      }
    }

    graphDL.make <- function(Countdata){
      if(tclvalue(Locus.var) == ""){
        tkmessageBox(message = "Select locus!", icon = "error", type = "ok")
      }else if(tclvalue(Allele.var) == ""){
        tkmessageBox(message = "Select allele!", icon = "error", type = "ok")
      }else{
        Countdata.oneL <- Countdata[[2]][[which(Countdata[[1]] == tclvalue(Locus.var))]]
        pos.Allele <- which(names(Countdata.oneL) == tclvalue(Allele.var))
        Expect <- (Countdata.oneL + 1) / sum(Countdata.oneL + 1)
        Mode <- Countdata.oneL / (sum(Countdata.oneL + 1) - length(Countdata.oneL))

        Density <- dbeta(seq(0, 1, 0.0001), Countdata.oneL[pos.Allele] + 1, sum(Countdata.oneL[- pos.Allele]) + 1)
        Pos <- which(Density / sum(Density) >= 0.00001)
        Xvalue <- seq(0, 1, 0.0001)[Pos]
        Yvalue <- Density[Pos]
        plot(Xvalue, Yvalue, type = "l", xlab = "Allele probability", ylab = "Density", 
             main = paste(tclvalue(Locus.var), ": allele ", tclvalue(Allele.var), sep = ""))
        abline(v = Mode[pos.Allele], col = 2)
        abline(v = Expect[pos.Allele], col = 3)
        legend(min(Xvalue) + ((max(Xvalue) - min(Xvalue)) * 2) / 3, max(Yvalue), lty = c(1, 1), col = c(3, 2), c("Expected value", "Mode"))
      }
    }

    Tab2.Freq_Locus <- function(Locus.names){
      tf.Tab2FL <- tktoplevel()
      tkwm.title(tf.Tab2FL, "Select locus")
      scr.Locus <- tkscrollbar(tf.Tab2FL, repeatinterval = 4, command = function(...) tkyview(list.Locus, ...))
      list.Locus <- tklistbox(tf.Tab2FL, selectmode = "Single", yscrollcommand = function(...) tkset(scr.Locus, ...), width = 20, background = "white", exportselection = 0)
      tkgrid(list.Locus, scr.Locus)
      tkgrid.configure(scr.Locus, rowspan = 4, sticky = "nsw")
      for(i in 1:length(Locus.names)){
        tkinsert(list.Locus, "end", Locus.names[i])
      }

      onOK_Locus <- function(){
        tclvalue(Locus.var) <- Locus.names[as.numeric(tkcurselection(list.Locus)) + 1]
        tclvalue(Allele.var) <- ""
        Tab2.Freq()
        tkdestroy(tf.Tab2FL)
      }

      Locus.OKbutton <-tk2button(tf.Tab2FL, text = "OK", width = -6, cursor = "hand2", command = function() onOK_Locus())
      tkgrid(Locus.OKbutton)
    }

    Tab2.Freq_Allele <- function(Countdata){
      tf.Tab2FA <- tktoplevel()
      tkwm.title(tf.Tab2FA, "Select allele")
      scr.Allele <- tkscrollbar(tf.Tab2FA, repeatinterval = 4, command = function(...) tkyview(list.Allele, ...))
      list.Allele <- tklistbox(tf.Tab2FA, selectmode = "Single", yscrollcommand = function(...) tkset(scr.Allele, ...), width = 20, background = "white", exportselection = 0)
      tkgrid(list.Allele, scr.Allele)
      tkgrid.configure(scr.Allele, rowspan = 4, sticky = "nsw")
      LID <- which(Countdata[[1]] == tclvalue(Locus.var))
      Allele.names <- names(Countdata[[2]][[LID]])
      for(i in 1:length(Allele.names)){
        tkinsert(list.Allele, "end", Allele.names[i])
      }

      onOK_Allele <- function(){
        tclvalue(Allele.var) <- Allele.names[as.numeric(tkcurselection(list.Allele)) + 1]
        Tab2.Freq()
        tkdestroy(tf.Tab2FA)
      }

      Allele.OKbutton <-tk2button(tf.Tab2FA, text = "OK", width = -6, cursor = "hand2", command = function() onOK_Allele())
      tkgrid(Allele.OKbutton)
    }

    Countdata.Freq.filepath <- tclVar("")
    Countdata.Freq.var <- tclVar("")
    Countdata.Freq.butt.state <- tclVar("disabled")
    Countdata.Freq.butt.cursor <- tclVar("arrow")
    Locus.var <- tclVar("")
    Allele.var <- tclVar("")

    tf.Freq <- tktoplevel()
    tkwm.title(tf.Freq, "Estimate expected allele frequencies")
    tabs.Freq <- tk2notebook(tf.Freq, tabs = c("Files", "Results"))
    tkpack(tabs.Freq, fill = "both", expand = 1)
    tab1.Freq <- tk2notetab(tabs.Freq, "Files")
    tab2.Freq <- tk2notetab(tabs.Freq, "Results")
    frame1.Freq <- tkframe(tab1.Freq)
    Tab1.Freq()
    frame2.Freq <- tkframe(tab2.Freq)
    Tab2.Freq()
  }

  RMP <- function(){
    RMP.simu.oneL <- function(Allele.count.oneL, Profile.oneL, theta, Nsimu, prior){
      RMPvec <- rep(0, Nsimu)
      Allele.count.oneL <- AlleleCountRemake(Allele.count.oneL, Profile.oneL)
      pos.Allele <- which(is.element(as.numeric(names(Allele.count.oneL)), Profile.oneL) == TRUE)
      Randfreq <- RandAlFreq(Allele.count.oneL, Nsimu, prior)
      Randfreq <- Randfreq[, pos.Allele, drop = FALSE]
      if(ncol(Randfreq) == 2){
        RMPvec <- 2 * (theta + (1 - theta) * Randfreq[, 1]) * (theta + (1 - theta) * Randfreq[, 2]) / ((1 + theta) * (1 + 2 * theta))
      }else{
        RMPvec <- (3 * theta + (1 - theta) * Randfreq[, 1]) * (2 * theta + (1 - theta) * Randfreq[, 1]) / ((1 + theta) * (1 + 2 * theta))
      }
      return(RMPvec)
    }

    RMP.simu.total <- function(Allele.count, Profile, theta, Nsimu, prior){
      NL <- length(Allele.count)
      RMPmat <- matrix(0, Nsimu, NL)
      for(i in 1:NL){
        RMPmat[, i] <- RMP.simu.oneL(Allele.count[[i]], Profile[i, ], theta, Nsimu, prior)
      }
      return(RMPmat)
    }

    RMP.expected.calc <- function(Allele.count, Profile, theta, prior){
      NL <- length(Allele.count)
      RMPvec <- rep(0, NL)
      for(i in 1:NL){
        Allele.count.oneL <- Allele.count[[i]]
        Profile.oneL <- Profile[i, ]
        Allele.count.oneL <- AlleleCountRemake(Allele.count.oneL, Profile.oneL)
        Expect <- ExpectAlFreq(Allele.count.oneL, prior)
        pos.Allele <- which(is.element(as.numeric(names(Allele.count.oneL)), Profile.oneL) == TRUE)
        if(length(unique(Profile.oneL)) == 2){
          alFreqs <- Expect[pos.Allele]
          RMPvec[i] <- 2 * (theta + (1 - theta) * alFreqs[1]) * (theta + (1 - theta) * alFreqs[2]) / ((1 + theta) * (1 + 2 * theta))
        }else if(length(unique(Profile.oneL)) == 1){
          RMPvec[i] <- (3 * theta + (1 - theta) * Expect[pos.Allele]) * (2 * theta + (1 - theta) * Expect[pos.Allele]) / ((1 + theta) * (1 + 2 * theta))
        }
      }
      return(RMPvec)
    }

    graphRMP.make <- function(RMPsimu, RMPexpected, RMP.CI, CIrange.RMP){
      RMP.Density <- density(log10(RMPsimu))
      plot(RMP.Density, xlab = expression(paste(Log[10]," (Random Match Probability)")), las = 1, main = "")
      abline(v = log10(RMPexpected), col = 3)
      par(xpd = TRUE)
      text(log10(RMPexpected), max(RMP.Density[[2]]) * 1.075, "Expected", col = 3)
      par(xpd = FALSE)
      abline(v = log10(RMP.CI[1]), col = 4)
      abline(v = log10(RMP.CI[2]), col = 4)
      segments(log10(RMP.CI[1]), max(RMP.Density[[2]]) * 0.0625, log10(RMP.CI[2]), max(RMP.Density[[2]]) * 0.0625, lwd = 4, col = 4)
      segments(log10(RMP.CI[1]), max(RMP.Density[[2]]) * 0.0425, log10(RMP.CI[1]), max(RMP.Density[[2]]) * 0.0825, lwd = 4, col = 4)
      segments(log10(RMP.CI[2]), max(RMP.Density[[2]]) * 0.0425, log10(RMP.CI[2]), max(RMP.Density[[2]]) * 0.0825, lwd = 4, col = 4)
      text(log10(RMP.CI[1]) + (log10(RMP.CI[2]) - log10(RMP.CI[1])) / 2, max(RMP.Density[[2]]) * 0.1, paste(round(100 * (CIrange.RMP[2] - CIrange.RMP[1]), 0), "% CI", sep = ""), col = 4)
      par(xpd = TRUE)
      text(log10(RMP.CI[1]), max(RMP.Density[[2]]) * 1.075, paste(CIrange.RMP[1] * 100, "% quantile", sep = ""), col = 4)
      text(log10(RMP.CI[2]), max(RMP.Density[[2]]) * 1.075, paste(CIrange.RMP[2] * 100, "% quantile", sep = ""), col = 4)
      par(xpd = FALSE)
    }

    reportRMP.make <- function(RMP.CI, CIrange.RMP, RMPsimu.quantile, RMPexpected){
      Save.as <- tkgetSaveFile(filetypes = "{{CSV Files} {.csv}}")
      if(tclvalue(Save.as) != ""){
        if(substr(tclvalue(Save.as), nchar(tclvalue(Save.as)) - 3, nchar(tclvalue(Save.as))) == ".csv"){
          ReportName <- tclvalue(Save.as)
        }else{
          ReportName <- paste(tclvalue(Save.as), ".csv", sep = "")
        }
        write.table("======== Software version ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)
        write.table(Version, file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Imput files ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Allele count data : ", tclvalue(Countdata.RMP.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Single source profile : ", tclvalue(Profile.RMP.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Parameters ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Theta value : ", tclvalue(theta.RMP.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Number of simulations : ", tclvalue(Nsimu.RMP.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Credible interval (%) : ", tclvalue(CI.RMP.var), " (two way)", sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Prior : ", tclvalue(Prior.RMP), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("======== ", tclvalue(CI.RMP.var), "% credible interval of RMP ========", sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c(paste(CIrange.RMP[1] * 100, "% quantile", sep = ""), RMP.CI[1]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c(paste(CIrange.RMP[2] * 100, "% quantile", sep = ""), RMP.CI[2]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Quantiles of RMP ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("Minimum", RMPsimu.quantile[1]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("1%", RMPsimu.quantile[2]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("5%", RMPsimu.quantile[3]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("25%", RMPsimu.quantile[4]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("Median", RMPsimu.quantile[5]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("75%", RMPsimu.quantile[6]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("95%", RMPsimu.quantile[7]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("99%", RMPsimu.quantile[8]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("Maximum", RMPsimu.quantile[9]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Expected RMP ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(RMPexpected, file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      }
    }

    alldataRMP.make <- function(RMPsimu){
      Save.as <- tkgetSaveFile(filetypes = "{{CSV Files} {.csv}}")
      if(tclvalue(Save.as) != ""){
        if(substr(tclvalue(Save.as), nchar(tclvalue(Save.as)) - 3, nchar(tclvalue(Save.as))) == ".csv"){
          ReportName <- tclvalue(Save.as)
        }else{
          ReportName <- paste(tclvalue(Save.as), ".csv", sep = "")
        }
        write.table(as.matrix(RMPsimu, ncol = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)
      }
    }

    Tab1.RMP <- function(){
      tkdestroy(frame1.RMP)
      frame1.RMP <<- tkframe(tab1.RMP)
      countdata.RMP.label <- tklabel(frame1.RMP, text = "Allele count data")
      dataname.RMP.label <- tklabel(frame1.RMP, textvariable = Countdata.RMP.var, width = 30, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      countdata.RMP.butt <- tkbutton(frame1.RMP, text = "    Load    ", cursor = "hand2", command = function() OpenFile(Countdata.RMP.filepath, Countdata.RMP.var, tf.RMP, Countdata.RMP.butt.state, Countdata.RMP.butt.cursor))
      tkgrid(countdata.RMP.label, dataname.RMP.label, countdata.RMP.butt, padx = 20, pady = 20, sticky = "w")

      profile.RMP.label <- tklabel(frame1.RMP, text = "Single source profile")
      profilename.RMP.label <- tklabel(frame1.RMP, textvariable = Profile.RMP.var, width = 30, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      profile.RMP.butt <- tkbutton(frame1.RMP, text = "    Load    ", cursor = "hand2", command = function() OpenFile(Profile.RMP.filepath, Profile.RMP.var, tf.RMP, Profile.RMP.butt.state, Profile.RMP.butt.cursor))
      tkgrid(profile.RMP.label, profilename.RMP.label, profile.RMP.butt, padx = 20, pady = 20, sticky = "w")

      theta.RMP.label <- tklabel(frame1.RMP, text = "Theta")
      theta.RMP.entry <- tkentry(frame1.RMP, textvariable = theta.RMP.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      tkgrid(theta.RMP.label, theta.RMP.entry, padx = 20, pady = 20, sticky = "w")

      Nsimu.RMP.label <- tklabel(frame1.RMP, text = "Number of simulations")
      Nsimu.RMP.entry <- tkentry(frame1.RMP, textvariable = Nsimu.RMP.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      tkgrid(Nsimu.RMP.label, Nsimu.RMP.entry, padx = 20, pady = 20, sticky = "w")

      CI.RMP.label <- tklabel(frame1.RMP, text = "Credible interval (%)")
      CI.RMP.entry <- tkentry(frame1.RMP, textvariable = CI.RMP.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      tkgrid(CI.RMP.label, CI.RMP.entry, padx = 20, pady = 20, sticky = "w")

      frame1_2.RMP <- tkframe(frame1.RMP)
      prior.RMP.rb1 <- tkradiobutton(frame1_2.RMP)
      prior.RMP.rb2 <- tkradiobutton(frame1_2.RMP)
      tkconfigure(prior.RMP.rb1, variable = Prior.RMP, value = "1")
      tkconfigure(prior.RMP.rb2, variable = Prior.RMP, value = "1 / K")
      tkgrid(prior.RMP.rb1, tklabel(frame1_2.RMP, text = "1"), prior.RMP.rb2, tklabel(frame1_2.RMP, text = "1 / K"))
      tkgrid(tklabel(frame1.RMP, text = "Prior"), frame1_2.RMP, padx = 20, pady = 20, sticky = "w")

      calculate.RMP.butt <- tkbutton(frame1.RMP, text = "    Calculate    ", cursor = "hand2", command = function() Tab2.RMP())
      tkgrid(calculate.RMP.butt, padx = 20, pady = 20, sticky = "w")
      tkgrid(frame1.RMP)
    }

    Tab2.RMP <- function(){
      tkdestroy(frame2.RMP)
      frame2.RMP <<- tkframe(tab2.RMP)
      if((tclvalue(Countdata.RMP.filepath) != "") && (tclvalue(Profile.RMP.filepath) != "")){
        Countdata <- Countdata.make(tclvalue(Countdata.RMP.filepath))
        Profile.info <- read.csv(tclvalue(Profile.RMP.filepath), header = TRUE)
        Profile.info <- as.matrix(Profile.info)
        Profile.L <- Profile.info[, 2]
        Profile <- matrix(as.numeric(Profile.info[, 3:4]), ncol = 2)
        Profile <- profileCheck(Countdata[[1]], Profile, Profile.L)
        if(length(Profile) == 0){
          tkmessageBox(message = "Locus names are invalid!", icon = "error", type = "ok")
        }else{
          RMPmat <- RMP.simu.total(Countdata[[2]], Profile, as.numeric(tclvalue(theta.RMP.var)), as.numeric(tclvalue(Nsimu.RMP.var)), tclvalue(Prior.RMP))
          RMPsimu <- apply(RMPmat, 1, prod)
          RMPsimu.quantile <- quantile(RMPsimu, probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1))
          RMPexpected.vec <- RMP.expected.calc(Countdata[[2]], Profile, as.numeric(tclvalue(theta.RMP.var)), tclvalue(Prior.RMP))
          CIrange.RMP <- CIrange.calc(as.numeric(tclvalue(CI.RMP.var)))
          RMP.CI <- quantile(RMPsimu, probs = CIrange.RMP)

          tkgrid(tklabel(frame2.RMP, text = paste(as.numeric(tclvalue(CI.RMP.var)), "% credible interval", sep = ""), font = "Helvetica 10 bold"), sticky = "w")
          tkgrid(tklabel(frame2.RMP, text = paste(signif(RMP.CI[1], digits = 3), "     -     ", signif(RMP.CI[2], digits = 3), sep = "")), sticky = "w")
          tkgrid(tklabel(frame2.RMP, text = ""), sticky = "w")

          tkgrid(tklabel(frame2.RMP, text = "Quantile", font = "Helvetica 10 bold"), sticky = "w")
          RMPmin.label <- tklabel(frame2.RMP, text = "Minimum")
          RMP1.label <- tklabel(frame2.RMP, text = "1%")
          RMP5.label <- tklabel(frame2.RMP, text = "5%")
          RMP25.label <- tklabel(frame2.RMP, text = "25%")
          RMPmedian.label <- tklabel(frame2.RMP, text = "Median")
          RMP75.label <- tklabel(frame2.RMP, text = "75%")
          RMP95.label <- tklabel(frame2.RMP, text = "95%")
          RMP99.label <- tklabel(frame2.RMP, text = "99%")
          RMPmax.label <- tklabel(frame2.RMP, text = "Max")
          RMPmin.val.label <- tklabel(frame2.RMP, text = paste(signif(RMPsimu.quantile[1], digits = 3)))
          RMP1.val.label <- tklabel(frame2.RMP, text = paste(signif(RMPsimu.quantile[2], digits = 3)))
          RMP5.val.label <- tklabel(frame2.RMP, text = paste(signif(RMPsimu.quantile[3], digits = 3)))
          RMP25.val.label <- tklabel(frame2.RMP, text = paste(signif(RMPsimu.quantile[4], digits = 3)))
          RMPmedian.val.label <- tklabel(frame2.RMP, text = paste(signif(RMPsimu.quantile[5], digits = 3)))
          RMP75.val.label <- tklabel(frame2.RMP, text = paste(signif(RMPsimu.quantile[6], digits = 3)))
          RMP95.val.label <- tklabel(frame2.RMP, text = paste(signif(RMPsimu.quantile[7], digits = 3)))
          RMP99.val.label <- tklabel(frame2.RMP, text = paste(signif(RMPsimu.quantile[8], digits = 3)))
          RMPmax.val.label <- tklabel(frame2.RMP, text = paste(signif(RMPsimu.quantile[9], digits = 3)))
          RMPexpected.val.label <- tklabel(frame2.RMP, text = paste(signif(prod(RMPexpected.vec), digits = 3)))
          tkgrid(RMPmin.label, RMPmin.val.label, sticky = "w")
          tkgrid(RMP1.label, RMP1.val.label, sticky = "w")
          tkgrid(RMP5.label, RMP5.val.label, sticky = "w")
          tkgrid(RMP25.label, RMP25.val.label, sticky = "w")
          tkgrid(RMPmedian.label, RMPmedian.val.label, sticky = "w")
          tkgrid(RMP75.label, RMP75.val.label, sticky = "w")
          tkgrid(RMP95.label, RMP95.val.label, sticky = "w")
          tkgrid(RMP99.label, RMP99.val.label, sticky = "w")
          tkgrid(RMPmax.label, RMPmax.val.label, sticky = "w")
          tkgrid(tklabel(frame2.RMP, text = ""), tklabel(frame2.RMP, text = ""), sticky = "w")
          tkgrid(tklabel(frame2.RMP, text = "Expected value", font = "Helvetica 10 bold"), sticky = "w")
          tkgrid(RMPexpected.val.label, sticky = "w")

          graph.RMP.butt <- tkbutton(frame2.RMP, text = "    Graph    ", cursor = "hand2", command = function() graphRMP.make(RMPsimu, prod(RMPexpected.vec), RMP.CI, CIrange.RMP))
          report.RMP.butt <- tkbutton(frame2.RMP, text = "    Report    ", cursor = "hand2", command = function() reportRMP.make(RMP.CI, CIrange.RMP, RMPsimu.quantile, prod(RMPexpected.vec)))
          alldata.RMP.butt <- tkbutton(frame2.RMP, text = "    Save all simulated RMP    ", cursor = "hand2", command = function() alldataRMP.make(RMPsimu))
          tkgrid(graph.RMP.butt, report.RMP.butt, alldata.RMP.butt, padx = 20, pady = 20, sticky = "w")
          tk2notetab.select(tabs.RMP, "Results")
        }
      }
      tkgrid(frame2.RMP)
    }

    Countdata.RMP.filepath <- tclVar("")
    Countdata.RMP.var <- tclVar("")
    Countdata.RMP.butt.state <- tclVar("disabled")
    Countdata.RMP.butt.cursor <- tclVar("arrow")
    Profile.RMP.filepath <- tclVar("")
    Profile.RMP.var <- tclVar("")
    Profile.RMP.butt.state <- tclVar("disabled")
    Profile.RMP.butt.cursor <- tclVar("arrow")
    theta.RMP.var <- tclVar("0")
    Nsimu.RMP.var <- tclVar("10000")
    CI.RMP.var <- tclVar("95")
    Prior.RMP <- tclVar("1")

    tf.RMP <- tktoplevel()
    tkwm.title(tf.RMP, "Random match probability")
    tabs.RMP <- tk2notebook(tf.RMP, tabs = c("Files", "Results"))
    tkpack(tabs.RMP, fill = "both", expand = 1)
    tab1.RMP <- tk2notetab(tabs.RMP, "Files")
    tab2.RMP <- tk2notetab(tabs.RMP, "Results")
    frame1.RMP <- tkframe(tab1.RMP)
    Tab1.RMP()
    frame2.RMP <- tkframe(tab2.RMP)
    Tab2.RMP()
  }

  CPI <- function(){
    CPI.simu.oneL <- function(Allele.count.oneL, Mixture.oneL, theta, Nsimu, prior){
      Allele.count.oneL <- AlleleCountRemake(Allele.count.oneL, Mixture.oneL)
      pos.Allele <- which(is.element(as.numeric(names(Allele.count.oneL)), Mixture.oneL) == TRUE)
      Randfreq <- RandAlFreq(Allele.count.oneL, Nsimu, prior)
      Randfreq <- Randfreq[, pos.Allele, drop = FALSE]
      if(length(Mixture.oneL) == 1){
        CPIcomp <- Randfreq[, 1] * theta + Randfreq[, 1]^2 * (1 - theta)
        return(CPIcomp)
      }else{
        pos.Hetero <- combn(length(Mixture.oneL), 2)
        CPIcomp <- matrix(0, Nsimu, ncol(pos.Hetero) + length(Mixture.oneL))
        for(i in 1:ncol(pos.Hetero)){
          posA1 <- pos.Hetero[1, i]
          posA2 <- pos.Hetero[2, i]
          CPIcomp[, i] <- 2 * Randfreq[, posA1] * Randfreq[, posA2] * (1 - theta)
        }
        for(i in 1:length(Mixture.oneL)){
          CPIcomp[, ncol(pos.Hetero) + i] <- Randfreq[, i] * theta + Randfreq[, i]^2 * (1 - theta)
        }
        return(apply(CPIcomp, 1, sum))
      }
    }

    CPI.simu.total <- function(Allele.count, Mixture, theta, Nsimu, prior){
      NL <- length(Allele.count)
      CPImat <- matrix(0, Nsimu, NL)
      Error <- FALSE
      for(i in 1:NL){
        Mixture.oneL <- Mixture[i, !is.na(Mixture[i, ])]
        if(length(Mixture.oneL) == 0){
          Error <- TRUE
          break
        }else{
          CPImat[, i] <- CPI.simu.oneL(Allele.count[[i]], Mixture.oneL, theta, Nsimu, prior)
        }
      }
      if(Error == TRUE){
        return("Error")
      }else{
        return(CPImat)
      }
    }

    CPI.expected.calc <- function(Allele.count, Mixture, theta, prior){
      NL <- length(Allele.count)
      CPIvec <- rep(0, NL)
      Error <- FALSE
      for(i in 1:NL){
        Allele.count.oneL <- Allele.count[[i]]
        Mixture.oneL <- Mixture[i, !is.na(Mixture[i, ])]
        Allele.count.oneL <- AlleleCountRemake(Allele.count.oneL, Mixture.oneL)
        Expect <- ExpectAlFreq(Allele.count.oneL, prior)
        if(length(Mixture.oneL) == 0){
          Error <- TRUE
          break
        }else if(length(Mixture.oneL) == 1){
          pos.Allele <- which(is.element(as.numeric(names(Allele.count.oneL)), Mixture.oneL) == TRUE)
          CPIvec[i] <- Expect[pos.Allele] * theta + Expect[pos.Allele]^2 * (1 - theta)
        }else{
          pos.Allele <- which(is.element(as.numeric(names(Allele.count.oneL)), Mixture.oneL) == TRUE)
          pos.Hetero <- combn(length(Mixture.oneL), 2)
          CPIcomp <- rep(0, ncol(pos.Hetero) + length(Mixture.oneL))
          for(j in 1:ncol(pos.Hetero)){
            A1 <- Mixture.oneL[pos.Hetero[1, j]]
            A2 <- Mixture.oneL[pos.Hetero[2, j]]
            posA <- which(is.element(as.numeric(names(Allele.count.oneL)), c(A1, A2)) == TRUE)
            CPIcomp[j] <- 2 * Expect[posA[1]] * Expect[posA[2]] * (1 - theta)
          }
          for(j in 1:length(Mixture.oneL)){
            A1 <- Mixture.oneL[j]
            posA <- which(is.element(as.numeric(names(Allele.count.oneL)), A1) == TRUE)
            CPIcomp[ncol(pos.Hetero) + j] <- Expect[posA] * theta + Expect[posA]^2 * (1 - theta)
          }
          CPIvec[i] <- sum(CPIcomp)
        }
      }
      if(Error == TRUE){
        return("Error")
      }else{
        return(CPIvec)
      }
    }

    graphCPI.make <- function(CPIsimu, CPIexpected, CPI.CI, CIrange.CPI){
      CPI.Density <- density(log10(CPIsimu))
      plot(CPI.Density, xlab = expression(paste(Log[10]," (Combined Probability of Inclusion)")), las = 1, main = "")
      abline(v = log10(CPIexpected), col = 3)
      par(xpd = TRUE)
      text(log10(CPIexpected), max(CPI.Density[[2]]) * 1.075, "Expected", col = 3)
      par(xpd = FALSE)
      abline(v = log10(CPI.CI[1]), col = 4)
      abline(v = log10(CPI.CI[2]), col = 4)
      segments(log10(CPI.CI[1]), max(CPI.Density[[2]]) * 0.0625, log10(CPI.CI[2]), max(CPI.Density[[2]]) * 0.0625, lwd = 4, col = 4)
      segments(log10(CPI.CI[1]), max(CPI.Density[[2]]) * 0.0425, log10(CPI.CI[1]), max(CPI.Density[[2]]) * 0.0825, lwd = 4, col = 4)
      segments(log10(CPI.CI[2]), max(CPI.Density[[2]]) * 0.0425, log10(CPI.CI[2]), max(CPI.Density[[2]]) * 0.0825, lwd = 4, col = 4)
      text(log10(CPI.CI[1]) + (log10(CPI.CI[2]) - log10(CPI.CI[1])) / 2, max(CPI.Density[[2]]) * 0.1, paste(round(100 * (CIrange.CPI[2] - CIrange.CPI[1]), 0), "% CI", sep = ""), col = 4)
      par(xpd = TRUE)
      text(log10(CPI.CI[1]), max(CPI.Density[[2]]) * 1.075, paste(CIrange.CPI[1] * 100, "% quantile", sep = ""), col = 4)
      text(log10(CPI.CI[2]), max(CPI.Density[[2]]) * 1.075, paste(CIrange.CPI[2] * 100, "% quantile", sep = ""), col = 4)
      par(xpd = FALSE)
    }

    reportCPI.make <- function(CPI.CI, CIrange.CPI, CPIsimu.quantile, CPIexpected){
      Save.as <- tkgetSaveFile(filetypes = "{{CSV Files} {.csv}}")
      if(tclvalue(Save.as) != ""){
        if(substr(tclvalue(Save.as), nchar(tclvalue(Save.as)) - 3, nchar(tclvalue(Save.as))) == ".csv"){
          ReportName <- tclvalue(Save.as)
        }else{
          ReportName <- paste(tclvalue(Save.as), ".csv", sep = "")
        }
        write.table("======== Software version ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)
        write.table(Version, file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Imput files ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Allele count data : ", tclvalue(Countdata.CPI.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Mixture profile : ", tclvalue(mixture.CPI.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Parameters ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Theta value : ", tclvalue(theta.CPI.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Number of simulations : ", tclvalue(Nsimu.CPI.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Credible interval (%) : ", tclvalue(CI.CPI.var), " (two way)", sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Prior : ", tclvalue(Prior.CPI), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("======== ", tclvalue(CI.CPI.var), "% credible interval of CPI ========", sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c(paste(CIrange.CPI[1] * 100, "% quantile", sep = ""), CPI.CI[1]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c(paste(CIrange.CPI[2] * 100, "% quantile", sep = ""), CPI.CI[2]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Quantiles of CPI ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("Minimum", CPIsimu.quantile[1]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("1%", CPIsimu.quantile[2]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("5%", CPIsimu.quantile[3]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("25%", CPIsimu.quantile[4]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("Median", CPIsimu.quantile[5]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("75%", CPIsimu.quantile[6]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("95%", CPIsimu.quantile[7]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("99%", CPIsimu.quantile[8]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("Maximum", CPIsimu.quantile[9]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Expected CPI ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(CPIexpected, file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      }
    }

    alldataCPI.make <- function(CPIsimu){
      Save.as <- tkgetSaveFile(filetypes = "{{CSV Files} {.csv}}")
      if(tclvalue(Save.as) != ""){
        if(substr(tclvalue(Save.as), nchar(tclvalue(Save.as)) - 3, nchar(tclvalue(Save.as))) == ".csv"){
          ReportName <- tclvalue(Save.as)
        }else{
          ReportName <- paste(tclvalue(Save.as), ".csv", sep = "")
        }
        write.table(as.matrix(CPIsimu, ncol = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)
      }
    }

    Tab1.CPI <- function(){
      tkdestroy(frame1.CPI)
      frame1.CPI <<- tkframe(tab1.CPI)

      countdata.CPI.label <- tklabel(frame1.CPI, text = "Allele count data")
      dataname.CPI.label <- tklabel(frame1.CPI, textvariable = Countdata.CPI.var, width = 30, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      countdata.CPI.butt <- tkbutton(frame1.CPI, text = "    Load    ", cursor = "hand2", command = function() OpenFile(Countdata.CPI.filepath, Countdata.CPI.var, tf.CPI, Countdata.CPI.butt.state, Countdata.CPI.butt.cursor))
      tkgrid(countdata.CPI.label, dataname.CPI.label, countdata.CPI.butt, padx = 20, pady = 20, sticky = "w")

      mixture.CPI.label <- tklabel(frame1.CPI, text = "Mixture profile")
      mixtureName.CPI.label <- tklabel(frame1.CPI, textvariable = mixture.CPI.var, width = 30, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      mixture.CPI.butt <- tkbutton(frame1.CPI, text = "    Load    ", cursor = "hand2", command = function() OpenFile(mixture.CPI.filepath, mixture.CPI.var, tf.CPI, mixture.CPI.butt.state, mixture.CPI.butt.cursor))
      tkgrid(mixture.CPI.label, mixtureName.CPI.label, mixture.CPI.butt, padx = 20, pady = 20, sticky = "w")

      theta.CPI.label <- tklabel(frame1.CPI, text = "Theta")
      theta.CPI.entry <- tkentry(frame1.CPI, textvariable = theta.CPI.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      tkgrid(theta.CPI.label, theta.CPI.entry, padx = 20, pady = 20, sticky = "w")

      Nsimu.CPI.label <- tklabel(frame1.CPI, text = "Number of simulations")
      Nsimu.CPI.entry <- tkentry(frame1.CPI, textvariable = Nsimu.CPI.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      tkgrid(Nsimu.CPI.label, Nsimu.CPI.entry, padx = 20, pady = 20, sticky = "w")

      CI.CPI.label <- tklabel(frame1.CPI, text = "Credible interval (%)")
      CI.CPI.entry <- tkentry(frame1.CPI, textvariable = CI.CPI.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      tkgrid(CI.CPI.label, CI.CPI.entry, padx = 20, pady = 20, sticky = "w")

      frame1_2.CPI <- tkframe(frame1.CPI)
      prior.CPI.rb1 <- tkradiobutton(frame1_2.CPI)
      prior.CPI.rb2 <- tkradiobutton(frame1_2.CPI)
      tkconfigure(prior.CPI.rb1, variable = Prior.CPI, value = "1")
      tkconfigure(prior.CPI.rb2, variable = Prior.CPI, value = "1 / K")
      tkgrid(prior.CPI.rb1, tklabel(frame1_2.CPI, text = "1"), prior.CPI.rb2, tklabel(frame1_2.CPI, text = "1 / K"))
      tkgrid(tklabel(frame1.CPI, text = "Prior"), frame1_2.CPI, padx = 20, pady = 20, sticky = "w")

      calculate.CPI.butt <- tkbutton(frame1.CPI, text = "    Calculate    ", cursor = "hand2", command = function() Tab2.CPI())
      tkgrid(calculate.CPI.butt, padx = 20, pady = 20, sticky = "w")
      tkgrid(frame1.CPI)
    }

    Tab2.CPI <- function(){
      tkdestroy(frame2.CPI)
      frame2.CPI <<- tkframe(tab2.CPI)
      if((tclvalue(Countdata.CPI.filepath) != "") && (tclvalue(mixture.CPI.filepath) != "")){
        Countdata <- Countdata.make(tclvalue(Countdata.CPI.filepath))
        Mixture.info <- read.csv(tclvalue(mixture.CPI.filepath), header = TRUE)
        Mixture.info <- as.matrix(Mixture.info)
        Mixture.L <- Mixture.info[, 2]
        Mixture <- matrix(as.numeric(Mixture.info[, 3:ncol(Mixture.info)]), nrow = length(Mixture.L))
        Mixture <- profileCheck(Countdata[[1]], Mixture, Mixture.L)
        if(length(Mixture) == 0){
          tkmessageBox(message = "Locus names are invalid!", icon = "error", type = "ok")
        }else{
          CPImat <- CPI.simu.total(Countdata[[2]], Mixture, as.numeric(tclvalue(theta.CPI.var)), as.numeric(tclvalue(Nsimu.CPI.var)), tclvalue(Prior.CPI))
          if(CPImat[1] == "Error"){
            tkmessageBox(message = "CPI cannot be calculated because of allelic drop-out!", icon = "error", type = "ok")
          }else{
            CPIsimu <- apply(CPImat, 1, prod)
            CPIsimu.quantile <- quantile(CPIsimu, probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1))
            CPIexpected.vec <- CPI.expected.calc(Countdata[[2]], Mixture, as.numeric(tclvalue(theta.CPI.var)), tclvalue(Prior.CPI))

            CIrange.CPI <- CIrange.calc(as.numeric(tclvalue(CI.CPI.var)))
            CPI.CI <- quantile(CPIsimu, probs = CIrange.CPI)

            tkgrid(tklabel(frame2.CPI, text = paste(as.numeric(tclvalue(CI.CPI.var)), "% credible interval", sep = ""), font = "Helvetica 10 bold"), sticky = "w")
            tkgrid(tklabel(frame2.CPI, text = paste(signif(CPI.CI[1], digits = 3), "     -     ", signif(CPI.CI[2], digits = 3), sep = "")), sticky = "w")
            tkgrid(tklabel(frame2.CPI, text = ""), sticky = "w")

            tkgrid(tklabel(frame2.CPI, text = "Quantile", font = "Helvetica 10 bold"), sticky = "w")
            CPImin.label <- tklabel(frame2.CPI, text = "Minimum")
            CPI1.label <- tklabel(frame2.CPI, text = "1%")
            CPI5.label <- tklabel(frame2.CPI, text = "5%")
            CPI25.label <- tklabel(frame2.CPI, text = "25%")
            CPImedian.label <- tklabel(frame2.CPI, text = "Median")
            CPI75.label <- tklabel(frame2.CPI, text = "75%")
            CPI95.label <- tklabel(frame2.CPI, text = "95%")
            CPI99.label <- tklabel(frame2.CPI, text = "99%")
            CPImax.label <- tklabel(frame2.CPI, text = "Max")
            CPImin.val.label <- tklabel(frame2.CPI, text = paste(signif(CPIsimu.quantile[1], digits = 3)))
            CPI1.val.label <- tklabel(frame2.CPI, text = paste(signif(CPIsimu.quantile[2], digits = 3)))
            CPI5.val.label <- tklabel(frame2.CPI, text = paste(signif(CPIsimu.quantile[3], digits = 3)))
            CPI25.val.label <- tklabel(frame2.CPI, text = paste(signif(CPIsimu.quantile[4], digits = 3)))
            CPImedian.val.label <- tklabel(frame2.CPI, text = paste(signif(CPIsimu.quantile[5], digits = 3)))
            CPI75.val.label <- tklabel(frame2.CPI, text = paste(signif(CPIsimu.quantile[6], digits = 3)))
            CPI95.val.label <- tklabel(frame2.CPI, text = paste(signif(CPIsimu.quantile[7], digits = 3)))
            CPI99.val.label <- tklabel(frame2.CPI, text = paste(signif(CPIsimu.quantile[8], digits = 3)))
            CPImax.val.label <- tklabel(frame2.CPI, text = paste(signif(CPIsimu.quantile[9], digits = 3)))
            CPIexpected.val.label <- tklabel(frame2.CPI, text = paste(signif(prod(CPIexpected.vec), digits = 3)))
            tkgrid(CPImin.label, CPImin.val.label, sticky = "w")
            tkgrid(CPI1.label, CPI1.val.label, sticky = "w")
            tkgrid(CPI5.label, CPI5.val.label, sticky = "w")
            tkgrid(CPI25.label, CPI25.val.label, sticky = "w")
            tkgrid(CPImedian.label, CPImedian.val.label, sticky = "w")
            tkgrid(CPI75.label, CPI75.val.label, sticky = "w")
            tkgrid(CPI95.label, CPI95.val.label, sticky = "w")
            tkgrid(CPI99.label, CPI99.val.label, sticky = "w")
            tkgrid(CPImax.label, CPImax.val.label, sticky = "w")
            tkgrid(tklabel(frame2.CPI, text = ""), tklabel(frame2.CPI, text = ""), sticky = "w")
            tkgrid(tklabel(frame2.CPI, text = "Expected value", font = "Helvetica 10 bold"), sticky = "w")
            tkgrid(CPIexpected.val.label, sticky = "w")

            graph.CPI.butt <- tkbutton(frame2.CPI, text = "    Graph    ", cursor = "hand2", command = function() graphCPI.make(CPIsimu, prod(CPIexpected.vec), CPI.CI, CIrange.CPI))
            report.CPI.butt <- tkbutton(frame2.CPI, text = "    Report    ", cursor = "hand2", command = function() reportCPI.make(CPI.CI, CIrange.CPI, CPIsimu.quantile, prod(CPIexpected.vec)))
            alldata.CPI.butt <- tkbutton(frame2.CPI, text = "    Save all simulated CPI    ", cursor = "hand2", command = function() alldataCPI.make(CPIsimu))
            tkgrid(graph.CPI.butt, report.CPI.butt, alldata.CPI.butt, padx = 20, pady = 20, sticky = "w")
            tk2notetab.select(tabs.CPI, "Results")
          }
        }
      }
      tkgrid(frame2.CPI)
    }

    Countdata.CPI.filepath <- tclVar("")
    Countdata.CPI.var <- tclVar("")
    Countdata.CPI.butt.state <- tclVar("disabled")
    Countdata.CPI.butt.cursor <- tclVar("arrow")
    mixture.CPI.filepath <- tclVar("")
    mixture.CPI.var <- tclVar("")
    mixture.CPI.butt.state <- tclVar("disabled")
    mixture.CPI.butt.cursor <- tclVar("arrow")
    theta.CPI.var <- tclVar("0")
    Nsimu.CPI.var <- tclVar("10000")
    CI.CPI.var <- tclVar("95")
    Prior.CPI <- tclVar("1")

    tf.CPI <- tktoplevel()
    tkwm.title(tf.CPI, "Combined probability of inclusion")
    tabs.CPI <- tk2notebook(tf.CPI, tabs = c("Files", "Results"))
    tkpack(tabs.CPI, fill = "both", expand = 1)
    tab1.CPI <- tk2notetab(tabs.CPI, "Files")
    tab2.CPI <- tk2notetab(tabs.CPI, "Results")
    frame1.CPI <- tkframe(tab1.CPI)
    Tab1.CPI()
    frame2.CPI <- tkframe(tab2.CPI)
    Tab2.CPI()
  }

  mixLR <- function(){
    Tab1.mixLR <- function(){
      tkdestroy(frame1.mixLR)
      frame1.mixLR <<- tkframe(tab1.mixLR)

      countdata.mixLR.label <- tklabel(frame1.mixLR, text = "Allele count data")
      dataname.mixLR.label <- tklabel(frame1.mixLR, textvariable = Countdata.mixLR.var, width = 30, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      countdata.mixLR.butt <- tkbutton(frame1.mixLR, text = "    Load    ", cursor = "hand2", command = function() OpenFile(Countdata.mixLR.filepath, Countdata.mixLR.var, tf.mixLR, Countdata.mixLR.butt.state, Countdata.mixLR.butt.cursor))
      tkgrid(countdata.mixLR.label, dataname.mixLR.label, countdata.mixLR.butt, padx = 20, pady = 20, sticky = "w")

      mixture.mixLR.label <- tklabel(frame1.mixLR, text = "Mixture profile")
      mixtureName.mixLR.label <- tklabel(frame1.mixLR, textvariable = mixture.mixLR.var, width = 30, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      mixture.mixLR.butt <- tkbutton(frame1.mixLR, text = "    Load    ", cursor = "hand2", command = function() OpenFile(mixture.mixLR.filepath, mixture.mixLR.var, tf.mixLR, mixture.mixLR.butt.state, mixture.mixLR.butt.cursor))
      tkgrid(mixture.mixLR.label, mixtureName.mixLR.label, mixture.mixLR.butt, padx = 20, pady = 20, sticky = "w")

      profile.mixLR.label <- tklabel(frame1.mixLR, text = "Reference profiles")
      profilename.mixLR.label <- tklabel(frame1.mixLR, textvariable = Profile.mixLR.var, width = 30, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      profile.mixLR.butt <- tkbutton(frame1.mixLR, text = "    Load    ", cursor = "hand2", command = function() OpenFile(Profile.mixLR.filepath, Profile.mixLR.var, tf.mixLR, Profile.mixLR.butt.state, Profile.mixLR.butt.cursor))
      tkgrid(profile.mixLR.label, profilename.mixLR.label, profile.mixLR.butt, padx = 20, pady = 20, sticky = "w")

      theta.mixLR.label <- tklabel(frame1.mixLR, text = "Theta")
      theta.mixLR.entry <- tkentry(frame1.mixLR, textvariable = theta.mixLR.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      tkgrid(theta.mixLR.label, theta.mixLR.entry, padx = 20, pady = 20, sticky = "w")

      Nsimu.mixLR.label <- tklabel(frame1.mixLR, text = "Number of simulations")
      Nsimu.mixLR.entry <- tkentry(frame1.mixLR, textvariable = Nsimu.mixLR.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      tkgrid(Nsimu.mixLR.label, Nsimu.mixLR.entry, padx = 20, pady = 20, sticky = "w")

      CI.mixLR.label <- tklabel(frame1.mixLR, text = "Credible interval (%)")
      CI.mixLR.entry <- tkentry(frame1.mixLR, textvariable = CI.mixLR.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      tkgrid(CI.mixLR.label, CI.mixLR.entry, padx = 20, pady = 20, sticky = "w")

      frame1_2.mixLR <- tkframe(frame1.mixLR)
      prior.mixLR.rb1 <- tkradiobutton(frame1_2.mixLR)
      prior.mixLR.rb2 <- tkradiobutton(frame1_2.mixLR)
      tkconfigure(prior.mixLR.rb1, variable = Prior.mixLR, value = "1")
      tkconfigure(prior.mixLR.rb2, variable = Prior.mixLR, value = "1 / K")
      tkgrid(prior.mixLR.rb1, tklabel(frame1_2.mixLR, text = "1"), prior.mixLR.rb2, tklabel(frame1_2.mixLR, text = "1 / K"))
      tkgrid(tklabel(frame1.mixLR, text = "Prior"), frame1_2.mixLR, padx = 20, pady = 20, sticky = "w")

      next.mixLR.butt <- tkbutton(frame1.mixLR, text = "    Next    ", cursor = "hand2", command = function() Tab2.mixLR())
      tkgrid(next.mixLR.butt, padx = 20, pady = 20, sticky = "w")
      tkgrid(frame1.mixLR)
    }

    Tab2.mixLR <- function(){
      tkdestroy(frame2.mixLR)
      frame2.mixLR <<- tkframe(tab2.mixLR)
      if((tclvalue(Countdata.mixLR.filepath) != "") && (tclvalue(mixture.mixLR.filepath) != "") && (tclvalue(Profile.mixLR.filepath) != "")){
        Countdata <- Countdata.make(tclvalue(Countdata.mixLR.filepath))
        Mixture.info <- read.csv(tclvalue(mixture.mixLR.filepath), header = TRUE)
        Mixture.info <- as.matrix(Mixture.info)
        Mixture.L <- Mixture.info[, 2]
        Mixture <- matrix(as.numeric(Mixture.info[, 3:ncol(Mixture.info)]), nrow = length(Mixture.L))
        Mixture <- profileCheck(Countdata[[1]], Mixture, Mixture.L)
        Error <- FALSE
        if(length(Mixture) == 0){
          Error <- TRUE
        }else{
          Profile.info <- read.csv(tclvalue(Profile.mixLR.filepath), header = TRUE)
          Profile.info <- as.matrix(Profile.info)
          knownP <- unique(Profile.info[, 1])
          Profile <- matrix(0, nrow(Profile.info), 2)
          for(i in 1:length(knownP)){
            knownPonePos <- which(Profile.info[, 1] == knownP[i])
            Profile.L <- Profile.info[knownPonePos, 2]
            Profile.one <- matrix(as.numeric(Profile.info[knownPonePos, 3:4]), nrow = length(Profile.L))
            Profile.one <- profileCheck(Countdata[[1]], Profile.one, Profile.L)
            if(length(Profile) == 0){
              Error <- TRUE
              break
            }else{
              Profile[(nrow(Profile.one) * (i - 1) + 1):(nrow(Profile.one) * i), ] <- Profile.one
            }
          }
        }
        if(Error){
          tkmessageBox(message = "Locus names are invalid!", icon = "error", type = "ok")
        }else{
          tclvalue(fileOK.mixLR.var) <- "1"
          frame2_Hp <- tkframe(frame2.mixLR, relief = "groove", borderwidth = 2)
          frame2_Hd <- tkframe(frame2.mixLR, relief = "groove", borderwidth = 2)
          tkgrid(tklabel(frame2_Hp, text = "Prosecutor hypothesis", font = "Helvetica 10 bold"), padx = 20, pady = 20, sticky = "w")
          tkgrid(tklabel(frame2_Hd, text = "Defense hypothesis", font = "Helvetica 10 bold"), padx = 20, pady = 20, sticky = "w")
          Hp.var.list <- Hd.var.list <- list()
          for(i in 1:length(knownP)){
            Hp.var.list[[i]] <- tclVar("0")
            Hd.var.list[[i]] <- tclVar("0")
            tkgrid(tkcheckbutton(frame2_Hp, text = paste(knownP[i]), variable = Hp.var.list[[i]]), padx = 20, pady = 20, sticky = "w")
            tkgrid(tkcheckbutton(frame2_Hd, text = paste(knownP[i]), variable = Hd.var.list[[i]]), padx = 20, pady = 20, sticky = "w")
          }
          Hp.unk.label <- tklabel(frame2_Hp, text = "Unknown contributor(s)")
          Hp.unk.entry <- tkentry(frame2_Hp, textvariable = NunkHp.mixLR.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
          Hd.unk.label <- tklabel(frame2_Hd, text = "Unknown contributor(s)")
          Hd.unk.entry <- tkentry(frame2_Hd, textvariable = NunkHd.mixLR.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
          tkgrid(Hp.unk.label, Hp.unk.entry, padx = 20, pady = 20, sticky = "w")
          tkgrid(Hd.unk.label, Hd.unk.entry, padx = 20, pady = 20, sticky = "w")
          tkgrid(frame2_Hp, frame2_Hd, padx = 20, pady = 20, sticky = "w")
          calc.mixLR.butt <- tkbutton(frame2.mixLR, text = "    Calculate    ", cursor = "hand2", command = function() Tab3.mixLR(Countdata, Mixture, Profile, Hp.var.list, Hd.var.list, knownP))
          tkgrid(calc.mixLR.butt, padx = 20, pady = 20, sticky = "w")
          tk2notetab.select(tabs.mixLR, "Hypotheses")
        }
      }
      tkgrid(frame2.mixLR)
    }

    Genotype.comb.make <- function(MixoneL, HNC){
      NumA <- length(MixoneL)
      if(NumA == 1){
        Gcomb <- matrix(1, 1, 2)
      }else{
        Gcomb <- rbind(cbind(1:NumA, 1:NumA), combinations(n = NumA, r = 2))
      }
      Ncomb <- nrow(Gcomb)
      CGcomb <- permutations(n = Ncomb, r = HNC, repeats.allowed = TRUE)
      Possible.combID <- numeric(0)
      for(i in 1:nrow(CGcomb)){
        GenoID <- as.vector(t(Gcomb[CGcomb[i, ], ]))
        if(all(is.element(1:NumA, GenoID))){
          Possible.combID <- c(Possible.combID, i)
        }
      }
      if(length(Possible.combID) == 0){
        return(NULL)
      }else{
        return(matrix(as.vector(t(Gcomb[as.vector(t(CGcomb[Possible.combID, ])), ])), ncol = 2 * HNC, byrow = TRUE))
      }
    }

    mixLike.simuE.oneL <- function(Allele.count.oneL, Mixture.oneL, Profile.oneL, Genotype.comb, KCpos, theta, Nsimu, prior){
      mixLRvec <- rep(0, Nsimu)
      Allele.count.oneL <- AlleleCountRemake(Allele.count.oneL, Mixture.oneL)
      pos.Allele <- which(is.element(as.numeric(names(Allele.count.oneL)), Mixture.oneL) == TRUE)
      Randfreq <- RandAlFreq(Allele.count.oneL, Nsimu, prior)
      Randfreq <- Randfreq[, pos.Allele, drop = FALSE]
      Expect <- ExpectAlFreq(Allele.count.oneL, prior)
      Expect <- Expect[pos.Allele]

      extract.count <- rep(0, length(Mixture.oneL))
      if(length(Profile.oneL) != 0){
        ProfileAl.count <- table(as.vector(Profile.oneL))
        for(i in 1:length(ProfileAl.count)){
          extract.pos <- which(Mixture.oneL == as.numeric(names(ProfileAl.count[i])))
          extract.count[extract.pos] <- extract.count[extract.pos] + ProfileAl.count[i]
        }
      }

      Error <- FALSE
      if(length(KCpos) != 0){
        combpos <- 1:nrow(Genotype.comb)
        for(i in 1:length(KCpos)){
          Profile.KC1 <- Profile.oneL[i, ]
          Profile.KC1.pos <- which(is.element(Mixture.oneL, Profile.KC1) == TRUE)
          if(length(Profile.KC1.pos) == 1){
            Profile.KC1.pos <- c(Profile.KC1.pos, Profile.KC1.pos)
          }
          combpos2 <- intersect(which(Genotype.comb[, 2 * (KCpos[i] - 1) + 1] == Profile.KC1.pos[1]), which(Genotype.comb[, 2 * KCpos[i]] == Profile.KC1.pos[2]))
          combpos <- intersect(combpos, combpos2)
          if(length(combpos) == 0){
            Error <- TRUE
            break
          }
        }
        UnkGenoComb <- Genotype.comb[combpos, -c(2 * KCpos - 1, 2 * KCpos), drop = FALSE]
      }else{
        UnkGenoComb <- Genotype.comb
      }
      if(Error){
        return(NULL)
      }else{
        likeTable <- matrix(0, Nsimu, nrow(UnkGenoComb))
        likeTableExpect <- rep(0, nrow(UnkGenoComb))
        if(length(UnkGenoComb) != 0){
          for(i in 1:nrow(UnkGenoComb)){
            UnkGenoComb.one <- UnkGenoComb[i, ]
            m <- sum(extract.count)
            EC <- extract.count
            P <- 1
            likeOneComb <- matrix(0, Nsimu, length(UnkGenoComb.one) / 2)
            likeOneCombExpect <- rep(0, length(UnkGenoComb.one) / 2)
            for(j in 1:(length(UnkGenoComb.one) / 2)){
              A1 <- UnkGenoComb.one[2 * j - 1]
              P1 <- (EC[A1] * theta + (1 - theta) * Randfreq[, A1]) / ((1 - theta) + m * theta)
              P1_Expect <- (EC[A1] * theta + (1 - theta) * Expect[A1]) / ((1 - theta) + m * theta)
              m <- m + 1
              EC[A1] <- EC[A1] + 1
              A2 <- UnkGenoComb.one[2 * j]
              P2 <- (EC[A2] * theta + (1 - theta) * Randfreq[, A2]) / ((1 - theta) + m * theta)
              P2_Expect <- (EC[A2] * theta + (1 - theta) * Expect[A2]) / ((1 - theta) + m * theta)
              m <- m + 1
              EC[A2] <- EC[A2] + 1
              if(A1 == A2){
                likeOneComb[, j] <- P1 * P2
                likeOneCombExpect[j] <- P1_Expect * P2_Expect
              }else{
                likeOneComb[, j] <- 2 * P1 * P2
                likeOneCombExpect[j] <- 2 * P1_Expect * P2_Expect
              }
            }
            likeTable[, i] <- apply(likeOneComb, 1, prod)
            likeTableExpect[i] <- prod(likeOneCombExpect)
          }
          return(list(apply(likeTable, 1, sum), sum(likeTableExpect)))
        }else{
          return(list(rep(1, Nsimu), 1))
        }
      }
    }

    mixLR.simuE.total <- function(Allele.count, Mixture, Profile, KCpos_Hp, KCpos_Hd, HNC_Hp, HNC_Hd, theta, Nsimu, prior){
      NL <- length(Allele.count)
      mixLikeSimuHp <- mixLikeSimuHd <- matrix(0, Nsimu, NL)
      mixLikeExpectHp <- mixLikeExpectHd <- rep(0, NL)
      Error <- FALSE
      if(HNC_Hp == HNC_Hd){
        for(i in 1:NL){
          Mixture.oneL <- Mixture[i, !is.na(Mixture[i, ])]
          Profile.oneL <- Profile[NL * (sort(unique(KCpos_Hp, KCpos_Hd)) - 1) + i, , drop = FALSE]
          Genotype.comb.Hp <- Genotype.comb.Hd <- Genotype.comb.make(Mixture.oneL, HNC_Hp)
          if(length(Genotype.comb.Hp) != 0){
            likedata.Hp <- mixLike.simuE.oneL(Allele.count[[i]], Mixture.oneL, Profile.oneL, Genotype.comb.Hp, KCpos_Hp, theta, Nsimu, prior)
            if(length(likedata.Hp) == 0){
              Error <- TRUE
              break
            }
            mixLikeSimuHp[, i] <- likedata.Hp[[1]]
            mixLikeExpectHp[i] <- likedata.Hp[[2]]
            likedata.Hd <- mixLike.simuE.oneL(Allele.count[[i]], Mixture.oneL, Profile.oneL, Genotype.comb.Hd, KCpos_Hd, theta, Nsimu, prior)
            if(length(likedata.Hd) == 0){
              Error <- TRUE
              break
            }
            mixLikeSimuHd[, i] <- likedata.Hd[[1]]
            mixLikeExpectHd[i] <- likedata.Hd[[2]]
          }else{
            Error <- TRUE
            break
          }
        }
      }else{
        for(i in 1:NL){
          Mixture.oneL <- Mixture[i, ]
          Profile.oneL <- Profile[NL * (sort(unique(KCpos_Hp, KCpos_Hd)) - 1) + i, , drop = FALSE]
          Genotype.comb.Hp <- Genotype.comb.make(Mixture.oneL, HNC_Hp)
          Genotype.comb.Hd <- Genotype.comb.make(Mixture.oneL, HNC_Hd)
          if(length(Genotype.comb.Hp) != 0){
            likedata.Hp <- mixLike.simu.oneL(Allele.count[[i]], Mixture.oneL, Profile.oneL, Genotype.comb.Hp, KCpos_Hp, theta, Nsimu, prior)
            if(length(likedata.Hp) == 0){
              Error <- TRUE
              break
            }
            mixLikeSimuHp[, i] <- likedata.Hp[[1]]
            mixLikeExpectHp[i] <- likedata.Hp[[2]]
            likedata.Hd <- mixLike.simu.oneL(Allele.count[[i]], Mixture.oneL, Profile.oneL, Genotype.comb.Hd, KCpos_Hd, theta, Nsimu, prior)
            if(length(likedata.Hd) == 0){
              Error <- TRUE
              break
            }
            mixLikeSimuHd[, i] <- likedata.Hd[[1]]
            mixLikeExpectHd[i] <- likedata.Hd[[2]]
          }else{
            Error <- TRUE
            break
          }
        }
      }
      if(Error){
        return(NULL)
      }else{
        return(list(mixLikeSimuHp, mixLikeExpectHp, mixLikeSimuHd, mixLikeExpectHd))
      }
    }

    graphmixLR.make <- function(mixLRsimu, mixLRexpected, mixLR.CI, CIrange.mixLR){
      mixLR.Density <- density(log10(mixLRsimu))
      plot(mixLR.Density, xlab = expression(paste(Log[10]," (LR)")), las = 1, main = "")
      abline(v = log10(mixLRexpected), col = 3)
      par(xpd = TRUE)
      text(log10(mixLRexpected), max(mixLR.Density[[2]]) * 1.075, "Expected", col = 3)
      par(xpd = FALSE)
      abline(v = log10(mixLR.CI[1]), col = 4)
      abline(v = log10(mixLR.CI[2]), col = 4)
      segments(log10(mixLR.CI[1]), max(mixLR.Density[[2]]) * 0.0625, log10(mixLR.CI[2]), max(mixLR.Density[[2]]) * 0.0625, lwd = 4, col = 4)
      segments(log10(mixLR.CI[1]), max(mixLR.Density[[2]]) * 0.0425, log10(mixLR.CI[1]), max(mixLR.Density[[2]]) * 0.0825, lwd = 4, col = 4)
      segments(log10(mixLR.CI[2]), max(mixLR.Density[[2]]) * 0.0425, log10(mixLR.CI[2]), max(mixLR.Density[[2]]) * 0.0825, lwd = 4, col = 4)
      text(log10(mixLR.CI[1]) + (log10(mixLR.CI[2]) - log10(mixLR.CI[1])) / 2, max(mixLR.Density[[2]]) * 0.1, paste(round(100 * (CIrange.mixLR[2] - CIrange.mixLR[1]), 0), "% CI", sep = ""), col = 4)
      par(xpd = TRUE)
      text(log10(mixLR.CI[1]), max(mixLR.Density[[2]]) * 1.075, paste(CIrange.mixLR[1] * 100, "% quantile", sep = ""), col = 4)
      text(log10(mixLR.CI[2]), max(mixLR.Density[[2]]) * 1.075, paste(CIrange.mixLR[2] * 100, "% quantile", sep = ""), col = 4)
      par(xpd = FALSE)
    }

    reportmixLR.make <- function(mixLR.CI, CIrange.mixLR, mixLRsimu.quantile, mixLRexpected, KCname_Hp, KCname_Hd){
      Save.as <- tkgetSaveFile(filetypes = "{{CSV Files} {.csv}}")
      if(tclvalue(Save.as) != ""){
        if(substr(tclvalue(Save.as), nchar(tclvalue(Save.as)) - 3, nchar(tclvalue(Save.as))) == ".csv"){
          ReportName <- tclvalue(Save.as)
        }else{
          ReportName <- paste(tclvalue(Save.as), ".csv", sep = "")
        }
        write.table("======== Software version ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)
        write.table(Version, file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Imput files ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Allele count data : ", tclvalue(Countdata.mixLR.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Mixture profile : ", tclvalue(mixture.mixLR.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Reference profile : ", tclvalue(Profile.mixLR.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Parameters ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Theta value : ", tclvalue(theta.mixLR.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Number of simulations : ", tclvalue(Nsimu.mixLR.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Credible interval (%) : ", tclvalue(CI.mixLR.var), " (two way)", sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Prior : ", tclvalue(Prior.mixLR), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Hypotheses ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Prosecutor hypothesis (Hp) : ", KCname_Hp, " + ", tclvalue(NunkHp.mixLR.var), " unknown contributor(s)", sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Defense hypothesis (Hd) : ", KCname_Hd, " + ", tclvalue(NunkHd.mixLR.var), " unknown contributor(s)", sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("======== ", tclvalue(CI.mixLR.var), "% credible interval of mixLR ========", sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c(paste(CIrange.mixLR[1] * 100, "% quantile", sep = ""), mixLR.CI[1]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c(paste(CIrange.mixLR[2] * 100, "% quantile", sep = ""), mixLR.CI[2]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Quantiles of mixLR ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("Minimum", mixLRsimu.quantile[1]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("1%", mixLRsimu.quantile[2]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("5%", mixLRsimu.quantile[3]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("25%", mixLRsimu.quantile[4]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("Median", mixLRsimu.quantile[5]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("75%", mixLRsimu.quantile[6]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("95%", mixLRsimu.quantile[7]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("99%", mixLRsimu.quantile[8]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("Maximum", mixLRsimu.quantile[9]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Expected mixLR ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(mixLRexpected, file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      }
    }

    alldatamixLR.make <- function(mixLRsimu){
      Save.as <- tkgetSaveFile(filetypes = "{{CSV Files} {.csv}}")
      if(tclvalue(Save.as) != ""){
        if(substr(tclvalue(Save.as), nchar(tclvalue(Save.as)) - 3, nchar(tclvalue(Save.as))) == ".csv"){
          ReportName <- tclvalue(Save.as)
        }else{
          ReportName <- paste(tclvalue(Save.as), ".csv", sep = "")
        }
        write.table(as.matrix(mixLRsimu, ncol = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)
      }
    }

    Tab3.mixLR <- function(Countdata = NULL, Mixture = NULL, Profile = NULL, Hp.var.list = NULL, Hd.var.list = NULL, knownP = NULL){
      tkdestroy(frame3.mixLR)
      frame3.mixLR <<- tkframe(tab3.mixLR)
      if(as.numeric(tclvalue(fileOK.mixLR.var)) == 1){
        KCpos_Hp <- which(as.numeric(sapply(Hp.var.list, tclvalue)) == 1)
        KCname_Hp <- paste(knownP[KCpos_Hp], collapse = " + ")
        NKC_Hp <- length(KCpos_Hp)
        NUC_Hp <- as.numeric(tclvalue(NunkHp.mixLR.var))
        HNC_Hp <- NKC_Hp + NUC_Hp
        KCpos_Hd <- which(as.numeric(sapply(Hd.var.list, tclvalue)) == 1)
        KCname_Hd <- paste(knownP[KCpos_Hd], collapse = " + ")
        NKC_Hd <- length(KCpos_Hd)
        NUC_Hd <- as.numeric(tclvalue(NunkHd.mixLR.var))
        HNC_Hd <- NKC_Hd + NUC_Hd
        likelihoods.data <- mixLR.simuE.total(Countdata[[2]], Mixture, Profile, KCpos_Hp, KCpos_Hd, HNC_Hp, HNC_Hd, as.numeric(tclvalue(theta.mixLR.var)), as.numeric(tclvalue(Nsimu.mixLR.var)), tclvalue(Prior.mixLR))
        if(length(likelihoods.data) == 0){
          tkmessageBox(message = "Hypotheses setting was wrong!", icon = "error", type = "ok")
        }else{
          mixLikeSimuHp <<- likelihoods.data[[1]]
          mixLikeExpectHp <<- likelihoods.data[[2]]
          mixLikeSimuHd <<- likelihoods.data[[3]]
          mixLikeExpectHd <<- likelihoods.data[[4]]
          mixLRsimu.perL <- mixLikeSimuHp / mixLikeSimuHd
          mixLRsimu <- apply(mixLRsimu.perL, 1, prod)
          mixLRsimu.quantile <- quantile(mixLRsimu, probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1))
          mixLRexpected <- mixLikeExpectHp / mixLikeExpectHd
          CIrange.mixLR <- CIrange.calc(as.numeric(tclvalue(CI.mixLR.var)))
          mixLR.CI <- quantile(mixLRsimu, probs = CIrange.mixLR)

          tkgrid(tklabel(frame3.mixLR, text = paste(as.numeric(tclvalue(CI.mixLR.var)), "% credible interval", sep = ""), font = "Helvetica 10 bold"), sticky = "w")
          tkgrid(tklabel(frame3.mixLR, text = paste(signif(mixLR.CI[1], digits = 3), "     -     ", signif(mixLR.CI[2], digits = 3), sep = "")), sticky = "w")
          tkgrid(tklabel(frame3.mixLR, text = ""), sticky = "w")

          tkgrid(tklabel(frame3.mixLR, text = "Quantile", font = "Helvetica 10 bold"), sticky = "w")
          mixLRmin.label <- tklabel(frame3.mixLR, text = "Minimum")
          mixLR1.label <- tklabel(frame3.mixLR, text = "1%")
          mixLR5.label <- tklabel(frame3.mixLR, text = "5%")
          mixLR25.label <- tklabel(frame3.mixLR, text = "25%")
          mixLRmedian.label <- tklabel(frame3.mixLR, text = "Median")
          mixLR75.label <- tklabel(frame3.mixLR, text = "75%")
          mixLR95.label <- tklabel(frame3.mixLR, text = "95%")
          mixLR99.label <- tklabel(frame3.mixLR, text = "99%")
          mixLRmax.label <- tklabel(frame3.mixLR, text = "Max")
          mixLRmin.val.label <- tklabel(frame3.mixLR, text = paste(signif(mixLRsimu.quantile[1], digits = 3)))
          mixLR1.val.label <- tklabel(frame3.mixLR, text = paste(signif(mixLRsimu.quantile[2], digits = 3)))
          mixLR5.val.label <- tklabel(frame3.mixLR, text = paste(signif(mixLRsimu.quantile[3], digits = 3)))
          mixLR25.val.label <- tklabel(frame3.mixLR, text = paste(signif(mixLRsimu.quantile[4], digits = 3)))
          mixLRmedian.val.label <- tklabel(frame3.mixLR, text = paste(signif(mixLRsimu.quantile[5], digits = 3)))
          mixLR75.val.label <- tklabel(frame3.mixLR, text = paste(signif(mixLRsimu.quantile[6], digits = 3)))
          mixLR95.val.label <- tklabel(frame3.mixLR, text = paste(signif(mixLRsimu.quantile[7], digits = 3)))
          mixLR99.val.label <- tklabel(frame3.mixLR, text = paste(signif(mixLRsimu.quantile[8], digits = 3)))
          mixLRmax.val.label <- tklabel(frame3.mixLR, text = paste(signif(mixLRsimu.quantile[9], digits = 3)))
          mixLRexpected.val.label <- tklabel(frame3.mixLR, text = paste(signif(prod(mixLRexpected), digits = 3)))
          tkgrid(mixLRmin.label, mixLRmin.val.label, sticky = "w")
          tkgrid(mixLR1.label, mixLR1.val.label, sticky = "w")
          tkgrid(mixLR5.label, mixLR5.val.label, sticky = "w")
          tkgrid(mixLR25.label, mixLR25.val.label, sticky = "w")
          tkgrid(mixLRmedian.label, mixLRmedian.val.label, sticky = "w")
          tkgrid(mixLR75.label, mixLR75.val.label, sticky = "w")
          tkgrid(mixLR95.label, mixLR95.val.label, sticky = "w")
          tkgrid(mixLR99.label, mixLR99.val.label, sticky = "w")
          tkgrid(mixLRmax.label, mixLRmax.val.label, sticky = "w")
          tkgrid(tklabel(frame3.mixLR, text = ""), tklabel(frame3.mixLR, text = ""), sticky = "w")
          tkgrid(tklabel(frame3.mixLR, text = "Expected value", font = "Helvetica 10 bold"), sticky = "w")
          tkgrid(mixLRexpected.val.label, sticky = "w")

          graph.mixLR.butt <- tkbutton(frame3.mixLR, text = "    Graph    ", cursor = "hand2", command = function() graphmixLR.make(mixLRsimu, prod(mixLRexpected), mixLR.CI, CIrange.mixLR))
          report.mixLR.butt <- tkbutton(frame3.mixLR, text = "    Report    ", cursor = "hand2", command = function() reportmixLR.make(mixLR.CI, CIrange.mixLR, mixLRsimu.quantile, prod(mixLRexpected), KCname_Hp, KCname_Hd))
          alldata.mixLR.butt <- tkbutton(frame3.mixLR, text = "    Save all simulated LR    ", cursor = "hand2", command = function() alldatamixLR.make(mixLRsimu))
          tkgrid(graph.mixLR.butt, report.mixLR.butt, alldata.mixLR.butt, padx = 20, pady = 20, sticky = "w")

          tk2notetab.select(tabs.mixLR, "Results")
        }
      }
      tkgrid(frame3.mixLR)
    }

    Countdata.mixLR.filepath <- tclVar("")
    Countdata.mixLR.var <- tclVar("")
    Countdata.mixLR.butt.state <- tclVar("disabled")
    Countdata.mixLR.butt.cursor <- tclVar("arrow")
    mixture.mixLR.filepath <- tclVar("")
    mixture.mixLR.var <- tclVar("")
    mixture.mixLR.butt.state <- tclVar("disabled")
    mixture.mixLR.butt.cursor <- tclVar("arrow")
    Profile.mixLR.filepath <- tclVar("")
    Profile.mixLR.var <- tclVar("")
    Profile.mixLR.butt.state <- tclVar("disabled")
    Profile.mixLR.butt.cursor <- tclVar("arrow")
    theta.mixLR.var <- tclVar("0")
    Nsimu.mixLR.var <- tclVar("10000")
    CI.mixLR.var <- tclVar("95")
    Prior.mixLR <- tclVar("1")
    NunkHp.mixLR.var <- tclVar("0")
    NunkHd.mixLR.var <- tclVar("1")
    fileOK.mixLR.var <- tclVar("0")

    tf.mixLR <- tktoplevel()
    tkwm.title(tf.mixLR, "Likelihood ratio for mixture interpretation")
    tabs.mixLR <- tk2notebook(tf.mixLR, tabs = c("Files", "Hypotheses", "Results"))
    tkpack(tabs.mixLR, fill = "both", expand = 1)
    tab1.mixLR <- tk2notetab(tabs.mixLR, "Files")
    tab2.mixLR <- tk2notetab(tabs.mixLR, "Hypotheses")
    tab3.mixLR <- tk2notetab(tabs.mixLR, "Results")
    frame1.mixLR <- tkframe(tab1.mixLR)
    Tab1.mixLR()
    frame2.mixLR <- tkframe(tab2.mixLR)
    Tab2.mixLR()
    frame3.mixLR <- tkframe(tab3.mixLR)
    Tab3.mixLR()
  }

  Kinship <- function(){
    Wenk.simu <- function(Query, Reference, Allele.count, Nsimu, theta, k2, k1, k0, prior){
      NL <- length(Allele.count)
      KI <- matrix(0, Nsimu, NL)
      for(i in 1:NL){
        qgt <- sort(Query[i, ])
        rgt <- sort(Reference[i, ])
        Allele.count.oneL <- Allele.count[[i]]
        Allele.count.oneL <- AlleleCountRemake(Allele.count.oneL, c(qgt, rgt))
        qgtpos <- which(is.element(names(Allele.count.oneL), qgt) == TRUE)
        Randfreq <- RandAlFreq(Allele.count.oneL, Nsimu, prior)
        if(length(qgtpos) == 1){
          Randfreq <- cbind(Randfreq[, qgtpos], Randfreq[, qgtpos])
        }else{
          Randfreq <- Randfreq[, qgtpos, drop = FALSE]
        }
        if(!any(is.element(qgt, rgt))){
          KI[, i] <- k0
        }else if(setequal(qgt, rgt)){
          if(qgt[1] == qgt[2]){
            A3 <- (2 * theta + (1 - theta) * Randfreq[, 1]) / (1 + theta)
            A4 <- (3 * theta + (1 - theta) * Randfreq[, 1]) / (1 + 2 * theta)
            H1 <- k2 + 2 * k1 * A4 + k0 * A3 * A4
            H2 <- A3 * A4
            KI[, i] <- H1 / H2
          }else{
            A4 <- (theta + (1 - theta) * Randfreq[, 1]) / (1 + 2 * theta)
            B4 <- (theta + (1 - theta) * Randfreq[, 2]) / (1 + 2 * theta)
            AB34 <- 2 * (theta + (1 - theta) * Randfreq[, 1]) / (1 + theta) * (theta + (1 - theta) * Randfreq[, 2]) / (1 + 2 * theta)
            H1 <- k2 + k1 * A4 + k1 * B4 + k0 * AB34
            H2 <- AB34
            KI[, i] <- H1 / H2
          }
        }else if(length(unique(c(qgt, rgt))) == 3){
          if(is.element(qgt[1], rgt)){
            H1 <- k1 + 2 * k0 * (theta + (1 - theta) * Randfreq[, 1]) / (1 + theta)
            H2 <- 2 * (theta + (1 - theta) * Randfreq[, 1]) / (1 + theta)
            KI[, i] <- H1 / H2
          }else{
            H1 <- k1 + 2 * k0 * (theta + (1 - theta) * Randfreq[, 2]) / (1 + theta)
            H2 <- 2 * (theta + (1 - theta) * Randfreq[, 2]) / (1 + theta)
            KI[, i] <- H1 / H2
          }
        }else{
          if(is.element(qgt[1], rgt)){
            H1 <- k1 + k0 * (2 * theta + (1 - theta) * Randfreq[, 1]) / (1 + theta)
            H2 <- (2 * theta + (1 - theta) * Randfreq[, 1]) / (1 + theta)
            KI[, i] <- H1 / H2
          }else{
            H1 <- k1 + k0 * (2 * theta + (1 - theta) * Randfreq[, 2]) / (1 + theta)
            H2 <- (2 * theta + (1 - theta) * Randfreq[, 2]) / (1 + theta)
            KI[, i] <- H1 / H2
          }
        }
      }
      return(KI)
    }

    Wenk.expected <- function(Query, Reference, Allele.count, theta, k2, k1, k0, prior){
      NL <- length(Allele.count)
      KI <- rep(0, NL)
      for(i in 1:NL){
        qgt <- sort(Query[i, ])
        rgt <- sort(Reference[i, ])
        Allele.count.oneL <- Allele.count[[i]]
        Allele.count.oneL <- AlleleCountRemake(Allele.count.oneL, c(qgt, rgt))
        qgtpos <- which(is.element(names(Allele.count.oneL), qgt) == TRUE)
        Expect <- ExpectAlFreq(Allele.count.oneL, prior)
        if(length(qgtpos) == 1){
          Expect <- cbind(Expect[qgtpos], Expect[qgtpos])
        }else{
          Expect <- Expect[qgtpos]
        }
        if(!any(is.element(qgt, rgt))){
          KI[i] <- k0
        }else if(setequal(qgt, rgt)){
          if(qgt[1] == qgt[2]){
            A3 <- (2 * theta + (1 - theta) * Expect[1]) / (1 + theta)
            A4 <- (3 * theta + (1 - theta) * Expect[1]) / (1 + 2 * theta)
            H1 <- k2 + 2 * k1 * A4 + k0 * A3 * A4
            H2 <- A3 * A4
            KI[i] <- H1 / H2
          }else{
            A4 <- (theta + (1 - theta) * Expect[1]) / (1 + 2 * theta)
            B4 <- (theta + (1 - theta) * Expect[2]) / (1 + 2 * theta)
            AB34 <- 2 * (theta + (1 - theta) * Expect[1]) / (1 + theta) * (theta + (1 - theta) * Expect[2]) / (1 + 2 * theta)
            H1 <- k2 + k1 * A4 + k1 * B4 + k0 * AB34
            H2 <- AB34
            KI[i] <- H1 / H2
          }
        }else if(length(unique(c(qgt, rgt))) == 3){
          if(is.element(qgt[1], rgt)){
            H1 <- k1 + 2 * k0 * (theta + (1 - theta) * Expect[1]) / (1 + theta)
            H2 <- 2 * (theta + (1 - theta) * Expect[1]) / (1 + theta)
            KI[i] <- H1 / H2
          }else{
            H1 <- k1 + 2 * k0 * (theta + (1 - theta) * Expect[2]) / (1 + theta)
            H2 <- 2 * (theta + (1 - theta) * Expect[2]) / (1 + theta)
            KI[i] <- H1 / H2
          }
        }else{
          if(is.element(qgt[1], rgt)){
            H1 <- k1 + k0 * (2 * theta + (1 - theta) * Expect[1]) / (1 + theta)
            H2 <- (2 * theta + (1 - theta) * Expect[1]) / (1 + theta)
            KI[i] <- H1 / H2
          }else{
            H1 <- k1 + k0 * (2 * theta + (1 - theta) * Expect[2]) / (1 + theta)
            H2 <- (2 * theta + (1 - theta) * Expect[2]) / (1 + theta)
            KI[i] <- H1 / H2
          }
        }
      }
      return(KI)
    }

    graphLRkin.make <- function(LRkin.simu, LRkin.expected, LRkin.CI, CIrange.LRkin){
      LRkin.Density <- density(log10(LRkin.simu))
      plot(LRkin.Density, xlab = expression(paste(Log[10]," (LR)")), las = 1, main = "")
      abline(v = log10(LRkin.expected), col = 3)
      par(xpd = TRUE)
      text(log10(LRkin.expected), max(LRkin.Density[[2]]) * 1.075, "Expected", col = 3)
      par(xpd = FALSE)
      abline(v = log10(LRkin.CI[1]), col = 4)
      abline(v = log10(LRkin.CI[2]), col = 4)
      segments(log10(LRkin.CI[1]), max(LRkin.Density[[2]]) * 0.0625, log10(LRkin.CI[2]), max(LRkin.Density[[2]]) * 0.0625, lwd = 4, col = 4)
      segments(log10(LRkin.CI[1]), max(LRkin.Density[[2]]) * 0.0425, log10(LRkin.CI[1]), max(LRkin.Density[[2]]) * 0.0825, lwd = 4, col = 4)
      segments(log10(LRkin.CI[2]), max(LRkin.Density[[2]]) * 0.0425, log10(LRkin.CI[2]), max(LRkin.Density[[2]]) * 0.0825, lwd = 4, col = 4)
      text(log10(LRkin.CI[1]) + (log10(LRkin.CI[2]) - log10(LRkin.CI[1]))/2, max(LRkin.Density[[2]]) * 0.1, paste(round(100 * (CIrange.LRkin[2] - CIrange.LRkin[1]), 0), "% CI", sep = ""), col = 4)
      par(xpd = TRUE)
      text(log10(LRkin.CI[1]), max(LRkin.Density[[2]]) * 1.075, paste(CIrange.LRkin[1] * 100, "% quantile", sep = ""), col = 4)
      text(log10(LRkin.CI[2]), max(LRkin.Density[[2]]) * 1.075, paste(CIrange.LRkin[2] * 100, "% quantile", sep = ""), col = 4)
      par(xpd = FALSE)
    }

    reportLRkin.make <- function(LRkin.CI, CIrange.LRkin, LRkin.simu.quantile, LRkin.expected){
      Save.as <- tkgetSaveFile(filetypes = "{{CSV Files} {.csv}}")
      if(tclvalue(Save.as) != ""){
        if(substr(tclvalue(Save.as), nchar(tclvalue(Save.as)) - 3, nchar(tclvalue(Save.as))) == ".csv"){
          ReportName <- tclvalue(Save.as)
        }else{
          ReportName <- paste(tclvalue(Save.as), ".csv", sep = "")
        }
        write.table("======== Software version ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)
        write.table(Version, file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Imput files ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Allele count data : ", tclvalue(Countdata.Kinship.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Query profile : ", tclvalue(Query.Kinship.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Reference profile : ", tclvalue(Reference.Kinship.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Parameters ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Theta value : ", tclvalue(theta.Kinship.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Number of simulations : ", tclvalue(Nsimu.Kinship.var), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Credible interval (%) : ", tclvalue(CI.Kinship.var), " (two way)", sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("Prior : ", tclvalue(Prior.Kinship), sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(paste("======== ", tclvalue(CI.Kinship.var), "% credible interval of Kinship ========", sep = ""), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c(paste(CIrange.LRkin[1] * 100, "% quantile", sep = ""), LRkin.CI[1]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c(paste(CIrange.LRkin[2] * 100, "% quantile", sep = ""), LRkin.CI[2]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Quantiles of LR (kinship) ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("Minimum", LRkin.simu.quantile[1]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("1%", LRkin.simu.quantile[2]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("5%", LRkin.simu.quantile[3]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("25%", LRkin.simu.quantile[4]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("Median", LRkin.simu.quantile[5]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("75%", LRkin.simu.quantile[6]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("95%", LRkin.simu.quantile[7]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("99%", LRkin.simu.quantile[8]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(matrix(c("Maximum", LRkin.simu.quantile[9]), nrow = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table("======== Expected LR (kinship) ========", file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        write.table(LRkin.expected, file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      }
    }

    alldataLRkin.make <- function(LRkin.simu){
      Save.as <- tkgetSaveFile(filetypes = "{{CSV Files} {.csv}}")
      if(tclvalue(Save.as) != ""){
        if(substr(tclvalue(Save.as), nchar(tclvalue(Save.as)) - 3, nchar(tclvalue(Save.as))) == ".csv"){
          ReportName <- tclvalue(Save.as)
        }else{
          ReportName <- paste(tclvalue(Save.as), ".csv", sep = "")
        }
        write.table(as.matrix(LRkin.simu, ncol = 1), file = ReportName, sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)
      }
    }

    Tab1.Kinship <- function(){
      tkdestroy(frame1.Kinship)
      frame1.Kinship <<- tkframe(tab1.Kinship)

      countdata.Kinship.label <- tklabel(frame1.Kinship, text = "Allele count data")
      dataname.Kinship.label <- tklabel(frame1.Kinship, textvariable = Countdata.Kinship.var, width = 30, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      countdata.Kinship.butt <- tkbutton(frame1.Kinship, text = "    Load    ", cursor = "hand2", command = function() OpenFile(Countdata.Kinship.filepath, Countdata.Kinship.var, tf.Kinship, Countdata.Kinship.butt.state, Countdata.Kinship.butt.cursor))
      tkgrid(countdata.Kinship.label, dataname.Kinship.label, countdata.Kinship.butt, padx = 20, pady = 20, sticky = "w")

      query.Kinship.label <- tklabel(frame1.Kinship, text = "Query profile")
      queryName.Kinship.label <- tklabel(frame1.Kinship, textvariable = Query.Kinship.var, width = 30, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      query.Kinship.butt <- tkbutton(frame1.Kinship, text = "    Load    ", cursor = "hand2", command = function() OpenFile(Query.Kinship.filepath, Query.Kinship.var, tf.Kinship, Query.Kinship.butt.state, Query.Kinship.butt.cursor))
      tkgrid(query.Kinship.label, queryName.Kinship.label, query.Kinship.butt, padx = 20, pady = 20, sticky = "w")

      reference.Kinship.label <- tklabel(frame1.Kinship, text = "Reference profile")
      referenceName.Kinship.label <- tklabel(frame1.Kinship, textvariable = Reference.Kinship.var, width = 30, highlightthickness = 1, relief = "groove", justify = "center", background = "white")
      reference.Kinship.butt <- tkbutton(frame1.Kinship, text = "    Load    ", cursor = "hand2", command = function() OpenFile(Reference.Kinship.filepath, Reference.Kinship.var, tf.Kinship, Reference.Kinship.butt.state, Reference.Kinship.butt.cursor))
      tkgrid(reference.Kinship.label, referenceName.Kinship.label, reference.Kinship.butt, padx = 20, pady = 20, sticky = "w")

      frame1_2.Kinship <- tkframe(frame1.Kinship)
      k2.entry <- tkentry(frame1_2.Kinship, textvariable = k2.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      k1.entry <- tkentry(frame1_2.Kinship, textvariable = k1.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      k0.entry <- tkentry(frame1_2.Kinship, textvariable = k0.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      tkgrid(k2.entry, k1.entry, k0.entry, padx = 20, pady = 20, sticky = "w")
      tkgrid(tklabel(frame1.Kinship, text = "Kinship coefficient (k2, k1, k0)"), frame1_2.Kinship, padx = 20, pady = 20, sticky = "w")

      theta.Kinship.label <- tklabel(frame1.Kinship, text = "Theta")
      theta.Kinship.entry <- tkentry(frame1.Kinship, textvariable = theta.Kinship.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      tkgrid(theta.Kinship.label, theta.Kinship.entry, padx = 20, pady = 20, sticky = "w")

      Nsimu.Kinship.label <- tklabel(frame1.Kinship, text = "Number of simulations")
      Nsimu.Kinship.entry <- tkentry(frame1.Kinship, textvariable = Nsimu.Kinship.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      tkgrid(Nsimu.Kinship.label, Nsimu.Kinship.entry, padx = 20, pady = 20, sticky = "w")

      CI.Kinship.label <- tklabel(frame1.Kinship, text = "Credible interval (%)")
      CI.Kinship.entry <- tkentry(frame1.Kinship, textvariable = CI.Kinship.var, width = 10, highlightthickness = 1, relief = "solid", justify = "center", background = "white")
      tkgrid(CI.Kinship.label, CI.Kinship.entry, padx = 20, pady = 20, sticky = "w")

      frame1_2.Kinship <- tkframe(frame1.Kinship)
      prior.Kinship.rb1 <- tkradiobutton(frame1_2.Kinship)
      prior.Kinship.rb2 <- tkradiobutton(frame1_2.Kinship)
      tkconfigure(prior.Kinship.rb1, variable = Prior.Kinship, value = "1")
      tkconfigure(prior.Kinship.rb2, variable = Prior.Kinship, value = "1 / K")
      tkgrid(prior.Kinship.rb1, tklabel(frame1_2.Kinship, text = "1"), prior.Kinship.rb2, tklabel(frame1_2.Kinship, text = "1 / K"))
      tkgrid(tklabel(frame1.Kinship, text = "Prior"), frame1_2.Kinship, padx = 20, pady = 20, sticky = "w")

      calculate.Kinship.butt <- tkbutton(frame1.Kinship, text = "    Calculate    ", cursor = "hand2", command = function() Tab2.Kinship())
      tkgrid(calculate.Kinship.butt, padx = 20, pady = 20, sticky = "w")
      tkgrid(frame1.Kinship)
    }

    Tab2.Kinship <- function(){
      tkdestroy(frame2.Kinship)
      frame2.Kinship <<- tkframe(tab2.Kinship)
      if(all(tclvalue(Countdata.Kinship.filepath) != "", tclvalue(Query.Kinship.filepath) != "", tclvalue(Reference.Kinship.filepath) != "")){
        Countdata <- Countdata.make(tclvalue(Countdata.Kinship.filepath))
        Query.info <- read.csv(tclvalue(Query.Kinship.filepath), header = TRUE)
        Query.info <- as.matrix(Query.info)
        Query.L <- Query.info[, 2]
        Query <- matrix(as.numeric(Query.info[, 3:4]), ncol = 2)
        Query <- profileCheck(Countdata[[1]], Query, Query.L)
        Reference.info <- read.csv(tclvalue(Reference.Kinship.filepath), header = TRUE)
        Reference.info <- as.matrix(Reference.info)
        Reference.L <- Reference.info[, 2]
        Reference <- matrix(as.numeric(Reference.info[, 3:4]), ncol = 2)
        Reference <- profileCheck(Countdata[[1]], Reference, Reference.L)
        if((length(Query) == 0) || (length(Reference) == 0)){
          tkmessageBox(message = "Locus names are invalid!", icon = "error", type = "ok")
        }else{
          LR.Kinship.mat <- Wenk.simu(Query, Reference, Countdata[[2]], as.numeric(tclvalue(Nsimu.Kinship.var)), as.numeric(tclvalue(theta.Kinship.var)), as.numeric(tclvalue(k2.var)), as.numeric(tclvalue(k1.var)), as.numeric(tclvalue(k0.var)), tclvalue(Prior.Kinship))
          LRkin.simu <- apply(LR.Kinship.mat, 1, prod)
          LRkin.simu.quantile <- quantile(LRkin.simu, probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1))
          LRkin.expected.vec <- Wenk.expected(Query, Reference, Countdata[[2]], as.numeric(tclvalue(theta.Kinship.var)), as.numeric(tclvalue(k2.var)), as.numeric(tclvalue(k1.var)), as.numeric(tclvalue(k0.var)), tclvalue(Prior.Kinship))
          CIrange.LRkin <- CIrange.calc(as.numeric(tclvalue(CI.Kinship.var)))
          LRkin.CI <- quantile(LRkin.simu, probs = CIrange.LRkin)

          tkgrid(tklabel(frame2.Kinship, text = paste(as.numeric(tclvalue(CI.Kinship.var)), "% credible interval", sep = ""), font = "Helvetica 10 bold"), sticky = "w")
          tkgrid(tklabel(frame2.Kinship, text = paste(signif(LRkin.CI[1], digits = 3), "     -     ", signif(LRkin.CI[2], digits = 3), sep = "")), sticky = "w")
          tkgrid(tklabel(frame2.Kinship, text = ""), sticky = "w")

          tkgrid(tklabel(frame2.Kinship, text = "Quantile", font = "Helvetica 10 bold"), sticky = "w")
          LRkin.min.label <- tklabel(frame2.Kinship, text = "Minimum")
          LRkin.1.label <- tklabel(frame2.Kinship, text = "1%")
          LRkin.5.label <- tklabel(frame2.Kinship, text = "5%")
          LRkin.25.label <- tklabel(frame2.Kinship, text = "25%")
          LRkin.median.label <- tklabel(frame2.Kinship, text = "Median")
          LRkin.75.label <- tklabel(frame2.Kinship, text = "75%")
          LRkin.95.label <- tklabel(frame2.Kinship, text = "95%")
          LRkin.99.label <- tklabel(frame2.Kinship, text = "99%")
          LRkin.max.label <- tklabel(frame2.Kinship, text = "Max")
          LRkin.min.val.label <- tklabel(frame2.Kinship, text = paste(signif(LRkin.simu.quantile[1], digits = 3)))
          LRkin.1.val.label <- tklabel(frame2.Kinship, text = paste(signif(LRkin.simu.quantile[2], digits = 3)))
          LRkin.5.val.label <- tklabel(frame2.Kinship, text = paste(signif(LRkin.simu.quantile[3], digits = 3)))
          LRkin.25.val.label <- tklabel(frame2.Kinship, text = paste(signif(LRkin.simu.quantile[4], digits = 3)))
          LRkin.median.val.label <- tklabel(frame2.Kinship, text = paste(signif(LRkin.simu.quantile[5], digits = 3)))
          LRkin.75.val.label <- tklabel(frame2.Kinship, text = paste(signif(LRkin.simu.quantile[6], digits = 3)))
          LRkin.95.val.label <- tklabel(frame2.Kinship, text = paste(signif(LRkin.simu.quantile[7], digits = 3)))
          LRkin.99.val.label <- tklabel(frame2.Kinship, text = paste(signif(LRkin.simu.quantile[8], digits = 3)))
          LRkin.max.val.label <- tklabel(frame2.Kinship, text = paste(signif(LRkin.simu.quantile[9], digits = 3)))
          LRkin.expected.val.label <- tklabel(frame2.Kinship, text = paste(signif(prod(LRkin.expected.vec), digits = 3)))
          tkgrid(LRkin.min.label, LRkin.min.val.label, sticky = "w")
          tkgrid(LRkin.1.label, LRkin.1.val.label, sticky = "w")
          tkgrid(LRkin.5.label, LRkin.5.val.label, sticky = "w")
          tkgrid(LRkin.25.label, LRkin.25.val.label, sticky = "w")
          tkgrid(LRkin.median.label, LRkin.median.val.label, sticky = "w")
          tkgrid(LRkin.75.label, LRkin.75.val.label, sticky = "w")
          tkgrid(LRkin.95.label, LRkin.95.val.label, sticky = "w")
          tkgrid(LRkin.99.label, LRkin.99.val.label, sticky = "w")
          tkgrid(LRkin.max.label, LRkin.max.val.label, sticky = "w")
          tkgrid(tklabel(frame2.Kinship, text = ""), tklabel(frame2.Kinship, text = ""), sticky = "w")
          tkgrid(tklabel(frame2.Kinship, text = "Expected value", font = "Helvetica 10 bold"), sticky = "w")
          tkgrid(LRkin.expected.val.label, sticky = "w")

          graph.LRkin.butt <- tkbutton(frame2.Kinship, text = "    Graph    ", cursor = "hand2", command = function() graphLRkin.make(LRkin.simu, prod(LRkin.expected.vec), LRkin.CI, CIrange.LRkin))
          report.LRkin.butt <- tkbutton(frame2.Kinship, text = "    Report    ", cursor = "hand2", command = function() reportLRkin.make(LRkin.CI, CIrange.LRkin, LRkin.simu.quantile, prod(LRkin.expected.vec)))
          alldata.LRkin.butt <- tkbutton(frame2.Kinship, text = "    Save all simulated LR    ", cursor = "hand2", command = function() alldataLRkin.make(LRkin.simu))
          tkgrid(graph.LRkin.butt, report.LRkin.butt, alldata.LRkin.butt, padx = 20, pady = 20, sticky = "w")
          tk2notetab.select(tabs.Kinship, "Results")
        }
      }
      tkgrid(frame2.Kinship)
    }

    Countdata.Kinship.filepath <- tclVar("")
    Countdata.Kinship.var <- tclVar("")
    Countdata.Kinship.butt.state <- tclVar("disabled")
    Countdata.Kinship.butt.cursor <- tclVar("arrow")
    Query.Kinship.filepath <- tclVar("")
    Query.Kinship.var <- tclVar("")
    Query.Kinship.butt.state <- tclVar("disabled")
    Query.Kinship.butt.cursor <- tclVar("arrow")
    Reference.Kinship.filepath <- tclVar("")
    Reference.Kinship.var <- tclVar("")
    Reference.Kinship.butt.state <- tclVar("disabled")
    Reference.Kinship.butt.cursor <- tclVar("arrow")
    k2.var <- tclVar("0")
    k1.var <- tclVar("0.5")
    k0.var <- tclVar("0")
    theta.Kinship.var <- tclVar("0")
    Nsimu.Kinship.var <- tclVar("10000")
    CI.Kinship.var <- tclVar("95")
    Prior.Kinship <- tclVar("1")

    tf.Kinship <- tktoplevel()
    tkwm.title(tf.Kinship, "Likelihood ratio for kinship analysis")
    tabs.Kinship <- tk2notebook(tf.Kinship, tabs = c("Files", "Results"))
    tkpack(tabs.Kinship, fill = "both", expand = 1)
    tab1.Kinship <- tk2notetab(tabs.Kinship, "Files")
    tab2.Kinship <- tk2notetab(tabs.Kinship, "Results")
    frame1.Kinship <- tkframe(tab1.Kinship)
    Tab1.Kinship()
    frame2.Kinship <- tkframe(tab2.Kinship)
    Tab2.Kinship()
  }

  tf <- tktoplevel()
  tkwm.title(tf, Version)
  frame.top <- tkframe(tf)
  Freq.butt <- tkbutton(frame.top, text = "    Estimate expected allele frequencies    ", cursor = "hand2", command = function() Freq())
  RMP.butt <- tkbutton(frame.top, text = "    Random match probability    ", cursor = "hand2", command = function() RMP())
  CPI.butt <- tkbutton(frame.top, text = "    Combined probability of inclusion    ", cursor = "hand2", command = function() CPI())
  mixLR.butt <- tkbutton(frame.top, text = "    Likelihood ratio for mixture interpretation    ", cursor = "hand2", command = function() mixLR())
  Kinship.butt <- tkbutton(frame.top, text = "    Likelihood ratio for kinship analysis    ", cursor = "hand2", command = function() Kinship())
  tkgrid(Freq.butt, padx = 20, pady = 20, sticky = "w")
  tkgrid(RMP.butt, padx = 20, pady = 20, sticky = "w")
  tkgrid(CPI.butt, padx = 20, pady = 20, sticky = "w")
  tkgrid(mixLR.butt, padx = 20, pady = 20, sticky = "w")
  tkgrid(Kinship.butt, padx = 20, pady = 20, sticky = "w")
  tkgrid(frame.top)
}



