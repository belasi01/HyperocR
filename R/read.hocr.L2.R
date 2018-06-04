# Read Ed files from HOCR processed at L2 using Prosoft

# The possol function is used to compute the sun geometry
#source('~/Copy/R/Cops/R/possol.R', echo=TRUE)

read.hocr.L2.Ed <- function(fn){
  print(paste("Reading:", fn))
  id = file(fn, "r")
  line = strsplit(readLines(con=id, n =1), " ") # Reads the first header line
  nrec = 1

  while (line != "character(0)"){
    line = strsplit(readLines(con=id, n =1), " ")
    nrec <- nrec+1

  }

  nHeaderLines = nrec + 1
  nrec = 0
  print(nHeaderLines)

  line = strsplit(readLines(con=id, n =1), " ")

  while (line != "character(0)"){
    line = strsplit(readLines(con=id, n =1), " ")
    nrec <- nrec+1

  }

  nDarkLines = nrec + 1
  nrec = 0
  print(nDarkLines)

  line = strsplit(readLines(con=id, n =1), " ")

  while (line != "character(0)"){
    line = strsplit(readLines(con=id, n =1), " ")
    nrec <- nrec+1

  }
  nEdLines = nrec
  print(nEdLines)
  close(id)

  id = file(fn, "r")
  Header = rep("NA", nHeaderLines-1)
  for (i in 1:(nHeaderLines-1)) {
    Header[i] = readLines(con=id, n =1)
  }
  close(id)

  # Reads the darks and extract the
  df.Ed.Darks =  read.table(fn, skip = nHeaderLines, nrows = (nDarkLines-3), header=T)
  if (names(df.Ed.Darks)[1] != "ES") {
    print("WARNING: Not an Ed file")
    return(0)
  }
  XLambda = names(df.Ed.Darks)[3:139]
  Ed.wl = as.numeric(str_sub(XLambda, 2,7))

  Darks = as.matrix(df.Ed.Darks[,3:139])

  if (ncol(df.Ed.Darks) == 148){
    DateDay = df.Ed.Darks[,146]
    Year = as.numeric(str_sub(DateDay, 1,4))
    DOY = as.numeric(str_sub(DateDay, 5,7))
    Hour = df.Ed.Darks[,147]

    Darks.Time = as.POSIXct(paste(Year,DOY,Hour),
                               format="%Y %j %H:%M:%S",tz="GMT")
  } else {
    Darks.Timer = df.Ed.Darks[,143]
  }


  # Read the Ed data

  df.Ed = read.table(fn, skip = nHeaderLines+nDarkLines, nrows = (nEdLines-1), header=T)

  Ed = as.matrix(df.Ed[,3:139])

  if (ncol(df.Ed) == 148){
    DateDay = df.Ed[,146]
    Year = as.numeric(str_sub(DateDay, 1,4))
    DOY = as.numeric(str_sub(DateDay, 5,7))
    Hour = df.Ed[,147]
    Ed.Time = as.POSIXct(paste(Year,DOY,Hour),
                         format="%Y %j %H:%M:%S",tz="GMT")

    Ed = list(Header=Header, Ed.wl=Ed.wl, Ed=Ed, Darks=Darks, Ed.Time =Ed.Time, Darks.Time = Darks.Time)

  } else {
    print('WARNING : No Date and Time TAGS in the file')
    Ed.Timer = df.Ed[,143]
    Ed = list(Header=Header, Ed.wl=Ed.wl, Ed=Ed, Darks=Darks, Ed.Timer = Ed.Timer, Darks.Timer = Darks.Timer)

  }

  return(Ed)

}

read.hocr.L2.Lu <- function(fn){
  print(paste("Reading:", fn))
  id = file(fn, "r")
  line = strsplit(readLines(con=id, n =1), " ") # Reads the first header line
  nrec = 1

  while (line != "character(0)"){
    line = strsplit(readLines(con=id, n =1), " ")
    nrec <- nrec+1

  }

  nHeaderLines = nrec + 1
  nrec = 0
  print(nHeaderLines)

  line = strsplit(readLines(con=id, n =1), " ")

  while (line != "character(0)"){
    line = strsplit(readLines(con=id, n =1), " ")
    nrec <- nrec+1

  }

  nDarkLines = nrec + 1
  nrec = 0
  print(nDarkLines)

  line = strsplit(readLines(con=id, n =1), " ")

  while (line != "character(0)"){
    line = strsplit(readLines(con=id, n =1), " ")
    nrec <- nrec+1

  }
  nLuLines = nrec
  print(nLuLines)
  close(id)

  id = file(fn, "r")
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
    df.Lu =  read.table(fn, skip = nHeaderLines, nrows = (nLuLines-3), header=T)
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

    df.Lu.Darks = read.table(fn, skip = nHeaderLines+nLuLines, nrows = (nDarkLines-1), header=T)

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
  df.Lu.Darks =  read.table(fn, skip = nHeaderLines, nrows = (nDarkLines-3), header=T)
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

  df.Lu = read.table(fn, skip = nHeaderLines+nDarkLines, nrows = (nLuLines-1), header=T)

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


