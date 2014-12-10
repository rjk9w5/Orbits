function nu = eccen2true(E,eccen)
	sqrt((1+eccen)/(1-eccen))*tan(E/2)
	nu = 2*atan(sqrt((1+eccen)/(1-eccen))*tan(E/2));
end