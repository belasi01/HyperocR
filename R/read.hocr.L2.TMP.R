#' Read HOCR files processed at L2 (*.dat) using Prosoft and converted into ASCII
#'
#' @param filen is the file name (*.dat)
#' @author Simon Bélanger
#'
#'@export

read.hocr.L2.TMP <- function(filen, RADIOMETERS=NA, APPLY.DARKS=TRUE){

  print(paste("Reading:", filen))

  if (is.na(RADIOMETERS)) {
    print("RADIOMETERS parameter must be specified!!")
    print("Accepted values are: ")
    print("Es for surface downwelling irradiance")
    print("EdZ for in-water downwelling irradiance")
    print("LuZ for in-water upwelling radiance")
    print("Lu0 for in-water upwelling radiance at surface")
    print("Lt for total surface radiance in air")
    print("... they should be separated by an underscore")
    print("Example: Es_Lu0_LuZ in the same order in which they are strored in the *.dat file")
    print("Abort reading")
    return(0)
  }
  instruments <- unlist(strsplit(RADIOMETERS, "_"))
  ninst <- length(instruments)
  print(paste("Number of instruments in the file:",
              ninst))
  print(instruments)
  if (ninst > 1) {
    print("WARNING: instrument list")
    print("provided in RADIOMETERS parameter")
    print("must be in the same order as they")
    print("appear in the L2 *.dat file...")
  }

  nDarks <- rep(0,ninst)
  nRecs  <- rep(0,ninst)
  list.Darks <- list()
  list.Darks.Time <- list()
  list.Meas <- list()
  list.Meas.Time <- list()
  mean.Darks <- matrix(NA, ncol=ninst, nrow = 137)
  waves <- matrix(NA, ncol=ninst, nrow = 137)


  # Count the number of header lines
  id = file(filen, "r")
  line = strsplit(readLines(con=id, n =1), " ") # Reads the first header line
  nrec = 1
  while (line != "character(0)"){
    line = unlist(strsplit(readLines(con=id, n =1), " "))
    nrec <- nrec+1
    if (line[1] == "CAL_FILE_NAMES") {
      calfiles <- unlist(strsplit(line[3],","))
      ncal = length(calfiles)
      print(paste("Number of calibration files is: ", ncal))
      print(calfiles)
      # Extract the first three charcters of the cal name
      Cal.Ids =  toupper(str_sub(calfiles,1,3))
      Block.Type = str_sub(Cal.Ids, 3,3)
      id.Dark = which(Block.Type =="D")
      id.Meas = which(Block.Type !="D")
    }
  }
  nHeaderLines = nrec - 1
  nrec = 0
  print(paste("Number of header lines:", nHeaderLines))

  # Count the number of dark lines of the first block
 # for (i in 1:ninst) {
 #   line = strsplit(readLines(con=id, n =1), " ")
 #  while (line != "character(0)"){
 #      line = strsplit(readLines(con=id, n =1), " ")
 #      nrec <- nrec+1
 #    }
 #    nDarks[i] = nrec
 #    nrec = 0
#    print(paste("Number of darks for", instruments[i],
 #               ":",nDarks[i]))
  #}

  #for (i in 1:ninst) {
   # line = strsplit(readLines(con=id, n =1), " ")
    #while (line != "character(0)"){
    #  line = strsplit(readLines(con=id, n =1), " ")
    #  nrec <- nrec+1
    #}
    #nRecs[i] <- nrec
  #  nrec = 0
  #  print(paste("Number of measurements for", instruments[i],
  #              ":",nRecs[i]))
  #}

  # Count the number of dark lines of each block
  nRecs = rep(0,ncal)
  mean.Darks <- matrix(NA, ncol=ncal, nrow = 137)
  waves <- matrix(NA, ncol=ncal, nrow = 137)

  for (i in 1:ncal) {
    line = strsplit(readLines(con=id, n =1), " ")
    while (line != "character(0)"){
          line = strsplit(readLines(con=id, n =1), " ")
          nrec <- nrec+1
    }
    nRecs[i] = nrec
    nrec = 0

    print(paste("Number of records for", calfiles[i],
                   ":",nRecs[i]))
  }

  close(id)

  # Reads and store the header
  id = file(filen, "r")
  Header = rep("NA", nHeaderLines)
  for (i in 1:(nHeaderLines)) {
    Header[i] = readLines(con=id, n =1)
  }
  close(id)

  # Adjust the actual numbers of lines for each parameters
  #nDarks <- nDarks - 2
  nRecs  <- nRecs  - 2

  # Reads darks and extract the wavelenght and time
  for (i in 1:ninst) {
    if (i==1) {
      skip = nHeaderLines + 2
      df.Darks =  read.table(filen, skip = skip, nrows = nDarks[i], header=T)
    } else
      {
      skip = skip + nDarks[i-1] + 3
      df.Darks =  read.table(filen, skip = skip , nrows = nDarks[i], header=T)
    }
    XLambda = names(df.Darks)[3:139]
    waves[,i] = as.numeric(str_sub(XLambda, 2,7))

    # extract darks
    Darks = as.matrix(df.Darks[,3:139])

    # extract time
    DateDay = df.Darks[,dim(df.Darks)[2]-2]
    Year = as.numeric(str_sub(DateDay, 1,4))
    DOY = as.numeric(str_sub(DateDay, 5,7))
    Hour = df.Darks[,dim(df.Darks)[2]-1]
    Darks.Time = as.POSIXct(paste(Year,DOY,Hour),
                               format="%Y %j %H:%M:%S",tz="GMT")

    list.Darks[[i]] <- Darks
    list.Darks.Time[[i]] <- Darks.Time

  }

  # Read the measurements
#  for (i in 1:ninst) {
  for (i in 1:ncal) {
    if (i==1) {
#      skip = skip + nDarks[ninst] + 3
      skip = skip  + 2
      df =  read.table(filen, skip = skip, nrows = nRecs[i], header=T)
    } else
    {
      skip = skip + nRecs[i-1] + 3
      df =  read.table(filen, skip = skip , nrows = (nRecs[i]), header=T)
    }

    XLambda = names(df)[3:139]
    waves[,i] = as.numeric(str_sub(XLambda, 2,7))

    Meas = as.matrix(df[,3:139])
    DateDay = df[,dim(df)[2]-2]
    Year = as.numeric(str_sub(DateDay, 1,4))
    DOY = as.numeric(str_sub(DateDay, 5,7))
    Hour = df[,dim(df)[2]-1]
    Meas.Time = as.POSIXct(paste(Year,DOY,Hour),
                         format="%Y %j %H:%M:%S",tz="GMT")

    list.Meas[[i]] <- Meas
    list.Meas.Time[[i]] <- Meas.Time
  }

  # compute mean darks
  for (i in 1:ninst) {
    mean.Darks[,i] = apply(list.Darks[[i]],2,mean, na.rm=T)
    if (APPLY.DARKS) {
      mean.Darks.m = matrix(mean.Darks[,i], ncol=137, nrow=nRecs[i], byrow = T)
      list.Meas[[i]] = list.Meas[[i]] - mean.Darks.m
    }
  }



  hocr <- list(instrument=instruments,
               Meas.Time=list.Meas.Time,
               Meas=list.Meas,
               Darks.Time=list.Darks.Time,
               Darks=list.Darks,
               waves=waves)

  return(hocr)

}

read.hocr.L2.Lu <- function(filen){


  nLuLines = nrec
  print(nLuLines)
  close(id)

  id = file(filen, "r")
  Header = rep("NA", nHeaderLines-1)
  for (i in 1:(nHeaderLines-1)) {
    Header[i] = readLines(con=id, n =1)
  }
  close(id)

########

  if (nLuLines < nDarkLines) {
    print("Darks writen after Lu data")
    x = nLuLines
    nLuLines = nDarkLines
    nDarkLines = x

    # Reads the Lu and extract the
    df.Lu =  read.table(filen, skip = nHeaderLines, nrows = (nLuLines-3), header=T)
    if (names(df.Lu)[1] != "LU") {
      print("WARNING: Not an Lu file")
      return(0)
    }
    XLambda = names(df.Lu)[3:139]
    Lu.wl = as.numeric(str_sub(XLambda, 2,7))

    Lu = as.matrix(df.Lu[,3:139])

    if (ncol(df.Lu) == 148){
      DateDay = df.Lu[,146]
      Year = as.numeric(str_sub(DateDay, 1,4))
      DOY = as.numeric(str_sub(DateDay, 5,7))
      Hour = df.Lu[,147]

      Lu.Time = as.POSIXct(paste(Year,DOY,Hour),
                              format="%Y %j %H:%M:%S",tz="GMT")
    } else {
      Lu.Timer = df.Lu[,144]
    }


    # Read the Darks data

    df.Lu.Darks = read.table(filen, skip = nHeaderLines+nLuLines, nrows = (nDarkLines-1), header=T)

    Darks = as.matrix(df.Lu.Darks[,3:139])

    if (ncol(df.Lu.Darks) == 148){
      DateDay = df.Lu.Darks[,146]
      Year = as.numeric(str_sub(DateDay, 1,4))
      DOY = as.numeric(str_sub(DateDay, 5,7))
      Hour = df.Lu.Darks[,147]
      Darks.Time = as.POSIXct(paste(Year,DOY,Hour),
                           format="%Y %j %H:%M:%S",tz="GMT")

      Lu = list(Header=Header, Lu.wl=Lu.wl, Lu=Lu, Darks=Darks, Lu.Time =Lu.Time, Darks.Time = Darks.Time)

    } else {
      print('WARNING : No Date and Time TAGS in the file')
      Darks.Timer = df.Lu.Darks[,144]
      Lu = list(Header=Header, Lu.wl=Lu.wl, Lu=Lu, Darks=Darks, Lu.Timer = Lu.Timer, Darks.Timer=Darks.Timer)

    }

    return(Lu)
  }  else {

  # Reads the darks and extract the
  df.Lu.Darks =  read.table(filen, skip = nHeaderLines, nrows = (nDarkLines-3), header=T)
  if (names(df.Lu.Darks)[1] != "LU") {
    print("WARNING: Not an Lu file")
    return(0)
  }
  XLambda = names(df.Lu.Darks)[3:139]
  Lu.wl = as.numeric(str_sub(XLambda, 2,7))

  Darks = as.matrix(df.Lu.Darks[,3:139])

  if (ncol(df.Lu.Darks) == 148){
    DateDay = df.Lu.Darks[,146]
    Year = as.numeric(str_sub(DateDay, 1,4))
    DOY = as.numeric(str_sub(DateDay, 5,7))
    Hour = df.Lu.Darks[,147]

    Darks.Time = as.POSIXct(paste(Year,DOY,Hour),
                            format="%Y %j %H:%M:%S",tz="GMT")
  } else {
    Darks.Timer = df.Lu.Darks[,143]
  }


  # Read the Lu data

  df.Lu = read.table(filen, skip = nHeaderLines+nDarkLines, nrows = (nLuLines-1), header=T)

  Lu = as.matrix(df.Lu[,3:139])

  if (ncol(df.Lu) == 148){
    DateDay = df.Lu[,146]
    Year = as.numeric(str_sub(DateDay, 1,4))
    DOY = as.numeric(str_sub(DateDay, 5,7))
    Hour = df.Lu[,147]
    Lu.Time = as.POSIXct(paste(Year,DOY,Hour),
                         format="%Y %j %H:%M:%S",tz="GMT")

    Lu = list(Header=Header, Lu.wl=Lu.wl, Lu=Lu, Darks=Darks, Lu.Time =Lu.Time, Darks.Time = Darks.Time)

  } else {
    print('WARNING : No Date and Time TAGS in the file')
    Lu.Timer = df.Lu[,144]
    Lu = list(Header=Header, Lu.wl=Lu.wl, Lu=Lu, Darks=Darks, Lu.Timer = Lu.Timer, Darks.Timer=Darks.Timer)

  }

  return(Lu)
  }

}
