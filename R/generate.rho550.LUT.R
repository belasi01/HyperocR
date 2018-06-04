generate.rho550.LUT <- function(filen){

  rho550_LUT <- array(NA, dim=c(10,13,10,10))

  nhead = 21
  block = 0
  for (wind in 1:10) {
    for (thetaS in 1:10) {
      # skip a given number of lines
      nskip=nhead+(block*119)

      # read the block of data
      tmp = scan(filen, skip=nskip, nlines=118)

      # extract the rho values from the block
      # rho is now a vector of 118 elements
      rho = tmp[(1:118)*3]

      # store the value in the 4D array
      # for thetaV = 0, repeat the value 13 times for each Dphi
      rho550_LUT[1,,wind, thetaS] = rep(rho[1],13)

      # fill the rest of the table
      # rearrage the rho vector into a 2D array
      rho550_LUT[2:10,,wind, thetaS] = matrix(rho[2:118], nrow=9, ncol=13, byrow = TRUE)


      block = block+1
    }
  }

  dimnames(rho550_LUT)<- list(thetaV=c("0","10","20","30","40","50","60","70","80","87.5"),
                              phiV=c("0","15","30","45","60","75","90","105","120","135","150","165","180"),
                              windspeed=c("0","2","4","5","6","8","10","12","14","15"),
                              thetaS=c("0","10","20","30","40","50","60","70","80","87.5"))

  return(rho550_LUT)

}
