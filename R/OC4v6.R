# TODO: Add comment
# 
# Author: simonbelanger
###############################################################################


OC4v6_fcn = function(X)
{
	X = log10(X)
	a = c(0.3272, -2.994, 2.7218, -1.2259, -0.5683)
	Y = 10^(a[1] + a[2]*X + a[3]*X^2 + a[4]*X^3 + a[5]*X^4)
	
	return(Y)
}
#Function Implementing the OC4v4 version 4 (operational SeaWIFS algorithm) 
OC4v6 = function(Rrs443, Rrs489, Rrs510, Rrs555)
{
	X = pmax(Rrs443, Rrs489, Rrs510)/Rrs555
	OC4v6 = OC4v6_fcn(X)
	
	return(OC4v6)
}
