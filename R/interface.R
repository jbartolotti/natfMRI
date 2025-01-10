

#' @export
natfMRI_templates = function(writedir = getwd()){
  file.copy(
    system.file('extdata','demoEvents.tsv', package = 'natfMRI'),
    file.path(writedir,'demoEvents.tsv'))
  file.copy(
    system.file('extdata','demoTimecourses.tsv', package = 'natfMRI'),
    file.path(writedir,'demoTimecourses.tsv'))
  file.copy(
    system.file('extdata','config.tsv', package = 'natfMRI'),
    file.path(writedir,'config.tsv'))
  file.copy(
    system.file('extdata','README.txt', package = 'natfMRI'),
    file.path(writedir,'README.txt'))


}


LOAD.file <- function(input, file_query = ''){
  # If it's not already a dataframe, load it
  if(!is.data.frame(input)){
    loadfromfile <- 'no'

    # if NA is provided (default behavior), then ask the user to select the file
    if(is.na(input)){
      if (!is.list(file_query) && is.character(file_query)){
        message(sprintf('Select File: %s', file_query))
        }
      name <- UTILS.getFilename('gui')
    } else {
      name <- input
    }

    # Now, if "name" is a string, try to load the filename
    if(is.character(name)){
      spl <- strsplit(name,'[.]')[[1]]
      extension <- tolower(spl[length(spl)])
      if(extension == 'tsv'){
        mydelim <- '\t'
      } else if(extension == 'csv'){
        mydelim <- ','
      } else{
        warning(sprintf('Could not determine delimiter from extension "%s" in file %s. Assuming tab-delimited.',extension, name))
        mydelim <- '\t'
      }
      imported_data <- read.delim(name, sep = mydelim)
    }
  } else {
    # It is a dataframe already, just return it
    imported_data <- input
  }

  return(imported_data)
}

#' @export
natfMRI <- function(events = NA, timecourses = NA, config = NA, savedir = getwd()){

  dat <- list()
  dat$ev <- LOAD.file(events, file_query = 'Events')
  dat$tc <- LOAD.file(events, file_query = 'Timecourses')

  return(dat)

}


doTheThing <- function(folder, config){
  mydata <- loadFolder(folder)

  # Configs: windowsize
  config <- list(
    windowsize = 20,
    windowtype = 'centered',
    windowoffset = 5
  )
  if(config$windowtype == 'centered'){
    config$window_adjust <- 0 - floor(config$windowsize/2) + config$windowoffset
  }

  # data frame with timecourses for each participant for each region.

  # data frame with events labeled by:
  # participant, group, eventnumber, eventtype(onset/offset), eventTR, condition(recall/forget)

  alldat <- demodata()

  # for a given event, extract the timecourse with the given windowsize for all participants in that group+condition (e.g., young+remember).
  # within-subject: correlate seed timecourse with all target timecourses, then z transform, then mean those together.

  procdat <- PROCESS.analyze(alldat,config)

 # end values:
  # for each combination of event, on/off, group, condition: zcor

  procdat$group = factor(procdat$group, levels = c('young','old'))
  procdat$type = factor(procdat$type, levels = c('onset','offset'))

  ggplot(procdat, aes(x = condition, y = fz, color = pid)) +
    theme_bw() +
    geom_quasirandom() +
    facet_wrap(group~type)

    ggplot(procdat, aes(x = condition, y = fz, color = condition)) +
    theme_bw() +
    stat_summary(aes(group = pid), fun='mean', color = 'black', alpha = .25) +
    stat_summary(aes(group = pid), fun='mean', geom = 'line', color = 'black', alpha = .25) +
    stat_summary(data = subset(procdat, condition == 'forget'), fun.data='mean_cl_normal', geom='errorbar', width = .1, size = 1, position = position_nudge(x = -.1)) +
    stat_summary(data = subset(procdat, condition == 'forget'), fun='mean', geom='point', size = 4, position = position_nudge(x = -.1)) +
      stat_summary(data = subset(procdat, condition == 'recall'), fun.data='mean_cl_normal', geom='errorbar', width = .1, size = 1, position = position_nudge(x = .1)) +
      stat_summary(data = subset(procdat, condition == 'recall'), fun='mean', geom='point', size = 4, position = position_nudge(x = .1)) +
      facet_wrap(group~type) +
      scale_color_manual(values = c('red','blue'), guide = 'none') +
      labs(y = 'HC-mPFC FC (z)', x = '')
    ggsave('hc-mpfc.png', width = 8, height = 6)


}


PROCESS.analyze <- function(dat, config){


  ev <- dat$ev
  tc <- dat$tc

  dd <- data.frame()

  events <- unique(ev$evnum)
  participants <- unique(ev$pid)
  etypes <- unique(ev$type)
  regions <- unique(tc$region)
  for(pp in participants){
    gg <- unique(ev$group[ev$pid == pp])
    for(ee in events){
      cc <- unique(ev$condition[ev$pid == pp & ev$evnum == ee])
      for(tt in etypes){
        ev_onset <- ev$time[ev$pid == pp & ev$evnum == ee & ev$type == tt]
        onset <- ev_onset + config$window_adjust
        duration <- config$windowsize

        tcsub <- subset(tc, pid == pp & tr >= onset & tr < (onset+duration))

        rcor <- cor(tcsub$value[tcsub$region == regions[1]], tcsub$value[tcsub$region == regions[2]])

        fz <- UTIL.fisherz(rcor)

        line <- c(pp, gg, ee, cc, tt, fz)
        dd <- rbind(dd,line)

      }
    }
  }
  colnames(dd) <- c('pid','group','evnum','condition','type','fz')

  dd$pid <- factor(dd$pid)
  dd$group <- factor(dd$group)
  dd$evnum <- as.numeric(dd$evnum)
  dd$condition <- factor(dd$condition)
  dd$type <- factor(dd$type)
  dd$fz <- as.numeric(dd$fz)

  return(dd)
}

UTIL.fisherz <- function(r){
  z <- (log(1+r) - log(1-r)) / 2
  return(z)
  }


demodata <- function(){
  timepoints <- 1260
  regions <- c('hippo','mpfc')
  participants <- 40
  groups <- c('young','old')
  events <- 20
  evdur <- 10
  conditions <- c('recall','forget')
  spacing <- floor((timepoints/(events+1)))
  evtimes_on <- (1:events)*spacing
  evtimes_off <- (1:events)*spacing+evdur



  evdat = data.frame()
  for(pp in 1:participants){
      gg <- groups[(pp %% length(groups))+1]
      for(ee in 1:events){
        cc <- conditions[(ee %% length(conditions))+1]

        evtime_on <- evtimes_on[ee]
        line <- c(pp,gg,ee,cc,'onset',evtime_on)
        evdat <- rbind(evdat,line)

        evtime_off <- evtimes_off[ee]
        line <- c(pp,gg,ee,cc,'offset',evtime_off)
        evdat <- rbind(evdat,line)

      }
  }
  colnames(evdat) <- c('pid','group','evnum','condition','type','time')
  evdat$evnum <- as.numeric(evdat$evnum)
  evdat$time <- as.numeric(evdat$time)

  tdat <- data.frame(
    pid = rep(1:participants,each=length(regions)*timepoints),
    region = rep(regions, participants, each = timepoints),
    tr = rep(1:timepoints,participants*length(regions)),
    value = runif(participants*length(regions)*timepoints)
    )
  alldat<-list(ev = evdat, tc = tdat)

  return(alldat)
}


loadFolder <- function(folder){
  mydata <- list()
  for(file in dir(folder)){
    mydata$file <- loadFile(file)
    }
  return(mydata)
}


