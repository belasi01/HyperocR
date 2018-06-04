




compute.PP.from.HyperSAS <- function (RRS, DAILY.PAR, Pb.max=2, Ek=20, MODEL="BELANGER") {




  Rrs443=mean(RRS[[1]]$Rrs[13:14])
  Rrs490=RRS[[1]]$Rrs[23]
  Rrs510=RRS[[1]]$Rrs[27]
  Rrs555=RRS[[1]]$Rrs[36]

  CHL = OC4v6(Rrs443, Rrs490, Rrs510, Rrs555)



}
